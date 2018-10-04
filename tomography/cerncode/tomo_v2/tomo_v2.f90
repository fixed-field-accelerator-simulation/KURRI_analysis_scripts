PROGRAM tomo_v2

! Copyright Mats Lindroos and Steve Hancock, CERN, Switzerland, 1998, v1
! Copyright Mats Lindroos and Steve Hancock, CERN, Switzerland, 2000, v2

! Numerical recipies routines
  USE nrstuff
! Longitudinal phase space 
  USE long_v2
! Generic subroutines
  USE tomosubs_v2

  IMPLICIT NONE

  TYPE pt  
     REAL(SP)   x
     REAL(SP)   y
  END TYPE pt

  TYPE(pt), DIMENSION(:,:), ALLOCATABLE :: points

! The data format for indata files

  REAL(SP), DIMENSION(:), ALLOCATABLE :: indata
  REAL(SP) renpt

! Variables for the MAPS

  INTEGER  i, j, barray, k, k2, ii, jj, pp, i1, j1
  REAL(SP), DIMENSION(:,:), ALLOCATABLE :: phasespace,profiles,&
                                           diffprofiles,&
                                           copyprofiles,sprofiles,dsprofiles
  REAL(SP), DIMENSION(:), ALLOCATABLE :: indarr,binsum

  INTEGER   i01,i02,fl

! Variables for smoothing

  INTEGER   gs_order,gs_nl,gs_nr,gs_np
  REAL(SP), DIMENSION(:,:), ALLOCATABLE :: resp
  INTEGER, DIMENSION(1) :: gs_shape
  REAL(SP), DIMENSION(1) :: gs_pad 
  REAL(SP), DIMENSION(:), ALLOCATABLE :: cp_of_p


! Variables for longtrack
  REAL(SP), DIMENSION(:), ALLOCATABLE :: xp,yp
  REAL(SP) cxp,cyp

! Sort for maps
  LOGICAL, DIMENSION(:), ALLOCATABLE :: xlog ! DIM will be (Npt*Npt)
  INTEGER, DIMENSION(:), ALLOCATABLE :: numapsarr
  INTEGER, DIMENSION(:), ALLOCATABLE :: xvec,xnumb ! DIM will be (Npt*Npt)
  INTEGER icount,xet,ll,actmaps,numaps,lengthxyp
  INTEGER lowlim,uplim,endprofile,twice,xfmlist,iioffset,nowturn
  INTEGER film,direction
  REAL(SP) numpts,is_out

! For rebinning
  INTEGER newplength

! For filmoutput
! lmyext gives 10**lmyext-1 as maximum number of pictures in a film  
  INTEGER, PARAMETER :: lmyext=3
  INTEGER i10
  CHARACTER(LEN=lmyext) myext 

! Variables for discrepancy calculation
  REAL(SP), DIMENSION(:), ALLOCATABLE :: darray
 
  CALL get_parameters

  ALLOCATE(indarr(profilelength))
  DO i=1,profilelength
     indarr(i)=REAL(i,SP)
  END DO

! Get data from file

  ALLOCATE(indata(alldata),profiles(profilecount,profilelength))
  CALL get_indata(indata)

  CAll subtract_baseline(indata)

! Get the profiles
  DO i=1,profilecount
    profiles(i,:) = indata((skipcount+i-1)*framelength+preskiplength+1:&
                           (skipcount+i)*framelength-postskiplength)
  END DO
  DEALLOCATE(indata)

  IF (rebin.gt.1) THEN
    ALLOCATE(copyprofiles(profilecount,profilelength),binsum(profilecount))
    copyprofiles=profiles
    DEALLOCATE(profiles)
    IF (MOD(profilelength,rebin).EQ.0) THEN
      newplength=profilelength/rebin
    ELSE
      newplength=profilelength/rebin + 1
    END IF
    ALLOCATE(profiles(profilecount,newplength))
    DO i=1,newplength-1
      binsum=0
      DO j=1,rebin
        binsum=binsum+copyprofiles(:,(i-1)*rebin+j)
      END DO
      profiles(:,i)=binsum
    END DO
    binsum=0
    DO j=(newplength-1)*rebin+1,profilelength
      binsum=binsum+copyprofiles(:,j)
    END DO
    profiles(:,newplength)=binsum*&
      REAL(rebin,SP)/REAL(profilelength-(newplength-1)*rebin,SP)
    profilelength=newplength
    DEALLOCATE(binsum,copyprofiles)
  END IF

! Make sure there are no negative values
  DO i=1,profilecount
    WHERE (profiles(i,:) < 0)
      profiles(i,:)=0.0_SP
    END WHERE
  END DO

! Calculate contents of reference profile
  CALL perimage(profiles(beam_ref_frame,:))

  DO i=1,profilecount
! Normalize to unit contents
    profiles(i,:)=profiles(i,:)/SUM(profiles(i,:))
  END DO

  CALL find_xat0(profiles(beam_ref_frame,:))

  yat0=REAL(profilelength,SP)/2.0_SP

  CALL out_profiles(profiles,'profiles.data')

  IF (self_field_flag.EQ.1) THEN 
!####################################################
!Calculate and store derivatives of smoothed profiles
!####################################################
    gs_order=4
    gs_nl=3
    gs_nr=3
    gs_np=gs_nl+gs_nr+1
    gs_shape(1)=2**CEILING(&
                LOG(REAL(profilelength+MAX(gs_nl,gs_nr),sp))/LOG(2.0_sp))
    ALLOCATE(resp(gs_np,profilecount),&
           dsprofiles(profilecount,gs_shape(1)),&
           cp_of_p(gs_shape(1)))
    gs_pad(1)=0.0_SP
!   Derivative of smoothed profiles 
    resp=SPREAD(savgol(gs_nl,gs_nr,1,gs_order),2,profilecount)
    DO k=1,profilecount
      cp_of_p=RESHAPE(profiles(k,:),gs_shape,gs_pad)
      dsprofiles(k,:)=convlv(cp_of_p,resp(:,k),1)
    END DO
    CALL calculate_self(dsprofiles)
  END IF
!####################################################
!Build maps
!####################################################  
  ALLOCATE(jmax(profilecount,profilelength),jmin(profilecount,profilelength))
  ALLOCATE(allbinmax(profilecount),allbinmin(profilecount))
  ALLOCATE(imin(profilecount),imax(profilecount))
  CALL ijlimits

  CALL out_plotinfo;

! Each bin to be populated by npt*npt bins for tracking
  ALLOCATE(points(Npt,Npt))
  points%x=SPREAD((2.0_SP*indarr(1:Npt)-1.0_SP)/(2.0_SP*Npt),2,Npt)
  points(1,:)%y=(2.0_SP*indarr(1:Npt)-1.0_SP)/(2.0_SP*Npt)
  points%y=SPREAD(points(1,:)%y,1,Npt)

! Allocate space for maps
! Calculate how many we need.

  ALLOCATE(numapsarr(profilecount))
  numapsarr=0
  DO pp=filmstart,filmstop,filmstep
    numapsarr(pp)=SUM(jmax(pp,imin(pp):imax(pp))-&
                  jmin(pp,imin(pp):imax(pp))+1)
  END DO
  numaps=profilecount*MAXVAL(numapsarr)
  DEALLOCATE(numapsarr)
! Calculate depth of maps
  fmlistlength=MIN(MAX(4,FLOOR(0.1_SP*REAL(profilelength,SP))),&
                   MIN(Npt*Npt,profilelength))
! Allocate space for maps and length of extended maps
  xlength=CEILING(0.1_SP*REAL(numaps,SP))
  xunit=xlength
  ALLOCATE(maps(profilecount,profilelength,profilelength))
  ALLOCATE(mapsi(numaps,fmlistlength))
  ALLOCATE(mapsweight(numaps,fmlistlength))
  ALLOCATE(mapsix(xlength,Npt*Npt-fmlistlength+1),&
           mapsweightx(xlength,Npt*Npt-fmlistlength+1))
  ALLOCATE(xp(numaps*Npt*Npt/profilecount),yp(numaps*Npt*Npt/profilecount))
  ALLOCATE(reverseweight(profilecount,profilelength))
  ALLOCATE(xlog(Npt*Npt),xnumb(Npt*Npt),xvec(Npt*Npt))  

! Make a sequence of reconstructed phase space pictures

  DO film=filmstart,filmstop,filmstep
    write(*,'(a,I3.3,a)') 'Image',film,':'

!   Make extension possible

    xindex=1

!   Initiate to negative values
    mapsi=-1
    mapsweight=0
    mapsix=-1
    mapsweightx=0
    maps=0
    reconstruct_p=film
!   Do first map
    actmaps=0
    DO ii=imin(reconstruct_p),imax(reconstruct_p)
      DO jj=jmin(reconstruct_p,ii),jmax(reconstruct_p,ii)
        actmaps=actmaps+1
        maps(reconstruct_p,ii,jj)=actmaps
        mapsi(maps(reconstruct_p,ii,jj),1)=ii
        mapsweight(maps(reconstruct_p,ii,jj),1)=Npt**2
      END DO
    END DO

!   Set initial conditions for points to be tracked    
    direction=1
    endprofile=profilecount
    DO twice=1,2
      nowturn=(reconstruct_p-1)*dturns
      k=0
      DO ii=imin(reconstruct_p),imax(reconstruct_p)
        DO JJ=jmin(reconstruct_p,ii),jmax(reconstruct_p,ii)
          DO i=1,Npt
            DO j=1,Npt
              k=k+1  
              xp(k)=REAL(ii-1,SP)+points(i,j)%x
              yp(k)=REAL(jj-1,SP)+points(i,j)%y
            END DO
          END DO
        END DO
      END DO

      DO PP=reconstruct_p+direction,endprofile,direction
        is_out=0
        IF (self_field_flag.EQ.1) THEN
          call longtrack_self(direction,dturns,yp(1:k),xp(1:k),nowturn)
        ELSE
          call longtrack(direction,dturns,yp(1:k),xp(1:k),nowturn)
        END IF          
!       calculate weight factors  
        iioffset=0      
        IILOOP:DO ii=imin(reconstruct_p),imax(reconstruct_p)
          JJLOOP:DO jj=jmin(reconstruct_p,ii),jmax(reconstruct_p,ii)
            actmaps=actmaps+1
            maps(PP,ii,jj)=actmaps
            lowlim=(jj-jmin(reconstruct_p,ii))*Npt*Npt+1+iioffset
            uplim=(jj-jmin(reconstruct_p,ii)+1)*Npt*Npt+iioffset
            xvec=CEILING(xp(lowlim:uplim))
            xlog=.TRUE.
            icount=0
            renpt=REAL(Npt*Npt,SP)
            DO ll=1,Npt*Npt
              IF (xlog(ll)) THEN
                xet=xvec(ll)
                xlog(ll)=.FALSE.
                xnumb=0
                xnumb(ll)=1
                WHERE(xlog.AND.(xvec.EQ.xet))
                  xlog=.FALSE.
                  xnumb=1
                END WHERE  
                IF ((xet.LT.1).OR.(xet.GT.profilelength)) THEN
                  is_out=is_out+1
                ELSE
                  icount=icount+1
                  IF (icount.LT.fmlistlength) THEN
                    mapsi(maps(PP,II,JJ),icount)=xet
                    mapsweight(maps(PP,II,JJ),icount)=SUM(xnumb)
                  ELSE
                    CALL extend_maps(II,JJ,PP,icount,xet,SUM(xnumb))
                  END IF
                END IF
              END IF
            END DO !ll
          END DO JJLOOP       
          iioffset=uplim
        END DO IILOOP
!        WRITE(*,'(a,I3,a,I3,a,F7.3,a)')&
!          ' Tracking from time slice ',reconstruct_p,' to ',PP,', ',&
!          100_SP*is_out/REAL(k,SP),'% went outside the image width.'
      END DO 
      direction=-1
      endprofile=1
    END DO

!   Just store the total weight factor (summed) for reversemaps to save memory
    reverseweight=0_SP
    DO PP=1,profilecount
      DO II=imin(reconstruct_p),imax(reconstruct_p)
        DO JJ=jmin(reconstruct_p,II),jmax(reconstruct_p,II)
          numpts=REAL(SUM(mapsweight(maps(PP,II,JJ),:)),SP)
          IF (mapsi(maps(PP,II,JJ),fmlistlength).LT.-1) THEN
            xfmlist=ABS(mapsi(maps(PP,II,JJ),fmlistlength))
            numpts=numpts+REAL(SUM(mapsweightx(xfmlist,:)),SP)
          ENDIF
          DO fl=1,Npt*Npt
            IF (fl.LT.fmlistlength) THEN
              IF (mapsi(maps(PP,II,JJ),fl).GT.0) THEN
                reverseweight(PP,mapsi(maps(PP,II,JJ),fl))= &
                reverseweight(PP,mapsi(maps(PP,II,JJ),fl))&
                +REAL(mapsweight(maps(PP,II,JJ),fl),SP)/numpts
              ELSE
                EXIT
              ENDIF
            ELSE
              IF (mapsi(maps(PP,II,JJ),fmlistlength).LT.-1) THEN
                IF (mapsix(xfmlist,fl-fmlistlength+1).GT.0) THEN
                  reverseweight(PP,mapsix(xfmlist,fl-fmlistlength+1))= &
                  reverseweight(PP,mapsix(xfmlist,fl-fmlistlength+1))&
                  +REAL(mapsweightx(xfmlist,fl-fmlistlength+1),SP)/numpts
                ELSE
                  EXIT
                END IF
              ELSE  
                EXIT
              END IF
            END IF
          END DO !fl
        END DO !JJ
      END DO !II
    END DO !PP

    reverseweight=reverseweight * REAL(profilecount,SP)

!####################################################
!...and finally we get around to do tomography
!####################################################

    ALLOCATE(phasespace(profilelength,profilelength))
    ALLOCATE(darray(niter+1))
    phasespace=backproject(profiles)

    DO i10=lmyext,1,-1
      myext(lmyext-i10+1:lmyext-i10+1)=&
        CHAR(48+INT(MODULO(reconstruct_p,10**(i10))/10**(i10-1)))
    END DO

    ALLOCATE(diffprofiles(profilecount,profilelength))
!   Improve picture by iterations and write "norm" for each step to array
    write(*,'(a)') ' Iterating...'
    DO ii=1,niter
!     Project and find difference from last projection
      diffprofiles=profiles-project(phasespace)
      darray(ii)=discrepancy(diffprofiles)
!     Improve picture
      phasespace=phasespace+backproject(diffprofiles)
!     Supress zeroes and normalize phasespace
      WHERE (phasespace<0.0_SP)
        phasespace=0_SP
      END WHERE
      IF (SUM(phasespace).LE.0_SP) THEN
        STOP 'Phasespace reduced to zeroes!'
      ELSE
        phasespace=phasespace/SUM(phasespace)
      END IF
    END DO

!   Calculate discrepancy for the last projection and write to file
    diffprofiles=profiles-project(phasespace)
    darray(niter+1)=discrepancy(diffprofiles)  
    CALL out_array(darray,'d'//myext//'.data')
!   Write final picture to file
    CALL out_picture(phasespace,&
         'image'//myext//'.data')

    DEALLOCATE(diffprofiles,darray,phasespace)


  END DO !film

  DEALLOCATE(yp,xp)
  DEALLOCATE(points)

END PROGRAM tomo_v2










