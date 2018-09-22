MODULE tomosubs_v2

! Copyright Mats Lindroos and Steve Hancock, CERN, Switzerland, 1998, v1
! Copyright Mats Lindroos and Steve Hancock, CERN, Switzerland, 2000, v2

  USE nrstuff
  USE long_v2
  IMPLICIT NONE

!-----------------------------------------
! maps: array in three dimensions. The first refer to the projection. The two last to 
!       the i and j coordinates of the physical square area in which the picture 
!       will be reconstructed. The map contains an integer
!       which is is the index of the arrays, maps and mapsweight which in turn
!       holds the actual data.
! mapsi: array in number of active points in maps and depth, mapsi holds the i in which 
!        mapsweight number of the orginally tracked Npt**2 number of points ended up.
! mapsweight: array in active points in maps and depth, see mapsi
! mapsix: array in number of extended mapsi vectors and depth, mapsix holds the overflow
!         from mapsi (last element dimension depth in mapsi holds the first index for
!         mapsix
! mapsweightx: array in number of extended mapsweight vectors and depth, holds the
!              overflow from mapsweight
! reverseweight: array in profilecount and profilelength, holds the sum of tracked
!                points at a certain i divided by the total number of tracked points 
!                launched (and still within valid limits)
! fmlistlength: initial depth of maps
! xlength: the number of extended vectors in maps
! xunit: the number of vectors additionally allocated every time a new allocation
!        takes place
! xindex: A global in this module holding the index of the next free vector
!-----------------------------------------

  INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: maps
  INTEGER(I2B), DIMENSION(:,:), ALLOCATABLE :: mapsi,mapsweight
  INTEGER(I2B), DIMENSION(:,:), ALLOCATABLE :: mapsix,mapsweightx
  REAL(SP), DIMENSION(:,:), ALLOCATABLE :: reverseweight
  INTEGER fmlistlength,xlength,xunit,xindex

CONTAINS

  SUBROUTINE extend_maps(II,JJ,PP,icountin,ix,weightx)
!-----------------------------------------
!   Routine to save overflow maps data.
!   xindex: is a module global holding the next free vector
!   xlength: is the presently allocated length of mapsx and mapsweightx
!   mapsix: is the extended data array mapsi
!   mapsweight: is the extended data array mapsweight
!-----------------------------------------
    INTEGER, INTENT(IN) :: icountin,II,JJ,PP,ix,weightx
    INTEGER icount,xfmlist
    icount=icountin
    IF (icount.EQ.fmlistlength) THEN
      xindex=xindex+1
      mapsi(maps(PP,II,JJ),fmlistlength)=-1*xindex
    ENDIF
    IF (xindex.GT.xlength) THEN
      write(*,*) 'Extending maps!'
      CALL reallocate
    ENDIF
    IF (icount.LT.Npt*Npt+2) THEN
      xfmlist=ABS(mapsi(maps(PP,II,JJ),fmlistlength))
      mapsix(xfmlist,icount-fmlistlength+1)=ix
      mapsweightx(xfmlist,icount-fmlistlength+1)=weightx
    ELSE
      STOP 'Ran out of indices in extend_maps!'
    ENDIF    
  END SUBROUTINE extend_maps

  SUBROUTINE reallocate()
!-----------------------------------------
!   Reallocation routine
!   cmapsi,cmapsweight: local copies during reallocation
!-----------------------------------------
    INTEGER(I2B), DIMENSION(SIZE(mapsix,1),SIZE(mapsix,2)) :: &
             cmapsweight,cmapsi
    cmapsi=mapsix
    cmapsweight=mapsweightx
    DEALLOCATE(mapsix,mapsweightx)
    xlength=xlength+xunit
    ALLOCATE(mapsix(xlength,Npt*Npt-fmlistlength+1),&
             mapsweightx(xlength,Npt*Npt-fmlistlength+1))
    mapsix=-1
    mapsweightx=0
    mapsix(1:xlength-xunit,:)=cmapsi
    mapsweightx(1:xlength-xunit,:)=cmapsweight
  END SUBROUTINE reallocate

  FUNCTION discrepancy(diffprofiles)
!-----------------------------------------
!   Calculate the discrepancy between projections and profiles
!-----------------------------------------
    USE nrstuff
    IMPLICIT NONE
    REAL(SP) discrepancy
    REAL(SP), DIMENSION(:,:) :: diffprofiles
    discrepancy=SQRT(SUM(diffprofiles**2)/(profilecount*profilelength)) 
  END FUNCTION discrepancy


  FUNCTION project(picture)
!-----------------------------------------
!   Project phasespace onto profilecount profiles
!-----------------------------------------
    USE nrstuff
    IMPLICIT NONE
    REAL(SP), DIMENSION(profilecount,profilelength) :: project
    REAL(SP), DIMENSION(:,:) :: picture
    INTEGER fl,pp,ii,jj,xfmlist
    REAL(SP) numpts
    project=0
    DO pp=1,profilecount
      DO ii=imin(reconstruct_p),imax(reconstruct_p)
        DO jj=jmin(reconstruct_p,ii),jmax(reconstruct_p,ii)
          numpts=REAL(SUM(mapsweight(maps(pp,ii,jj),:)),SP)
          IF (mapsi(maps(pp,ii,jj),fmlistlength).LT.-1) THEN
            xfmlist=ABS(mapsi(maps(pp,ii,jj),fmlistlength))
            numpts=numpts+REAL(SUM(mapsweightx(xfmlist,:)),SP)
          ENDIF
          DO fl=1,Npt*Npt
            IF (fl.LT.fmlistlength) THEN
              IF (mapsi(maps(pp,ii,jj),fl).GT.0) THEN
                project(pp,mapsi(maps(pp,ii,jj),fl))= &
                project(pp,mapsi(maps(pp,ii,jj),fl))+ &
                REAL(mapsweight(maps(pp,ii,jj),fl),SP)/numpts*picture(ii,jj)
              ELSE
                 EXIT
              END IF
            ELSE
              IF (mapsi(maps(pp,ii,jj),fmlistlength).LT.-1) THEN
                IF (mapsix(xfmlist,fl-fmlistlength+1).GT.0) THEN
                  project(pp,mapsix(xfmlist,fl-fmlistlength+1))= &
                  project(pp,mapsix(xfmlist,fl-fmlistlength+1))+ &
                  REAL(mapsweightx(xfmlist,fl-fmlistlength+1),SP)&
                  /numpts*picture(ii,jj)
                ELSE
                  EXIT
                END IF
              ELSE
                EXIT
              END IF          
            END IF
          END DO !fl
        END DO !jj
      END DO !ii
    END DO !pp   
  END FUNCTION project

  FUNCTION backproject(profiles)
!-----------------------------------------
!   Backproject profilecount profiles onto profilelength**2 phasespace
!-----------------------------------------
    USE nrstuff
    IMPLICIT NONE
    REAL(SP), DIMENSION(profilelength,profilelength) :: backproject
    REAL(SP), DIMENSION(:,:) :: profiles
    INTEGER i,j,p,fl,xfmlist
    REAL(SP) numpts
    backproject=0
    DO p=1,profilecount
      DO i=imin(reconstruct_p),imax(reconstruct_p)
        DO j=jmin(reconstruct_p,i),jmax(reconstruct_p,i)
          numpts=REAL(SUM(mapsweight(maps(p,i,j),:)),SP)
          IF (mapsi(maps(p,i,j),fmlistlength).LT.-1) THEN
            xfmlist=ABS(mapsi(maps(p,i,j),fmlistlength))
            numpts=numpts+REAL(SUM(mapsweightx(xfmlist,:)),SP)
          END IF
          DO fl=1,Npt*Npt
            IF (fl.LT.fmlistlength) THEN
              IF (mapsi(maps(p,i,j),fl).GT.0) THEN
                IF (reverseweight(p,mapsi(maps(p,i,j),fl)).LE.0_SP) THEN
                  WRITE(*,*) 'p,i,j,fl,mapsi,maps,mapsweight: ',p,i,j,fl,&
                              mapsi(maps(p,i,j),fl),maps(p,i,j),&
                              mapsweight(maps(p,i,j),fl)
                  STOP 'Would have divided by zero in backproject.'
                END IF
                backproject(i,j)=backproject(i,j)&
                +REAL(mapsweight(maps(p,i,j),fl),SP)/numpts&
                /reverseweight(p,mapsi(maps(p,i,j),fl))&
                *profiles(p,mapsi(maps(p,i,j),fl))
              ELSE
                EXIT
              END IF
            ELSE
              IF (mapsi(maps(p,i,j),fmlistlength).LT.-1) THEN
                IF (mapsix(xfmlist,fl-fmlistlength+1).GT.0) THEN
                  IF (reverseweight(p,mapsix(xfmlist,fl-fmlistlength+1))&
                      .LE.0_SP) THEN
                    STOP 'Would have divided by zero in backproject.'
                  END IF
                  backproject(i,j)=backproject(i,j)&
                  +REAL(mapsweightx(xfmlist,fl-fmlistlength+1),SP)/numpts&
                  /reverseweight(p,mapsix(xfmlist,fl-fmlistlength+1))&
                  *profiles(p,mapsix(xfmlist,fl-fmlistlength+1))
                ELSE
                  EXIT
                END IF
              ELSE
                EXIT
              END IF          
            END IF
          END DO !fl
        END DO !j
      END DO !i
    END DO !p
  END FUNCTION backproject

END MODULE tomosubs_v2



