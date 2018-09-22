MODULE long_v2

! Copyright Mats Lindroos and Steve Hancock, CERN, Switzerland, 1998, v1
! Copyright Mats Lindroos and Steve Hancock, CERN, Switzerland, 2000, v2

  USE nrstuff

  TYPE jlim    !Type declaration for upper limits of j (buildmaps)
     INTEGER u
     INTEGER l
  END TYPE jlim

!-----------------------------------------
! There are two reference systems: i) the reconstructed phase space 
! coordinate system which has its origin at some minimum energy and minimum 
! phase, with x-bins dtbin wide and y-bins dEbin wide; and ii) the physical 
! coordinates expressed in energy and rf phase wrt the synchronous particle.
!
! xat0: synchronous phase in bins in beam_ref_frame 
! yat0: synchronous energy (0 in relative terms) in reconstructed phase 
!       space coordinate system 
! rebin: rebinning factor
! dturns: number of machine turns between each measurement
! preskiplength: subtract this number of bins from the beginning of
!                the 'raw' input traces
! postskiplength: subtract this number of bins from the end of
!                 the 'raw' input traces
! skipcount: skip this number of traces from the beginning of the 'raw'
!            input file
! framecount: number of traces in the 'raw' input file
! framelength: length of each trace in the 'raw' input file
! alldata: total number of datapoints in the 'raw' input file
! allbinmin,allbinmax: imin and imax as calculated from jmax and jmin
! Npt: Npt**2 is the number of test particles tracked from each pixel of 
!      reconstructed phase space
! niter: number of iterations in reconstruction process
! reconstruct_p: profile at which phase space is reconstructed
! machine_ref_frame: frame to which machine parameters are referenced (B0,VRF1,VRF2)
! beam_ref_frame: frame to which beam parameters are referenced (baseline, phi0)
! filmstep: step between consecutive reconstructions for the profiles from
!           filmstart to filmstop
!-----------------------------------------

  REAL(SP)     xat0,yat0,xorigin
  INTEGER      preskiplength,postskiplength,skipcount
  INTEGER      framecount,framelength
  INTEGER      profilecount,profilelength
  INTEGER      alldata,iminin,imaxin
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: jmax,jmin
  INTEGER, DIMENSION(:), ALLOCATABLE :: imin,imax,allbinmin,allbinmax
  INTEGER      Npt,niter,reconstruct_p,machine_ref_frame,beam_ref_frame
  INTEGER      filmstart,filmstop,filmstep,rebin,dturns,self_field_flag

  INTERFACE trajectoryheight
     MODULE PROCEDURE trajectoryheight_var,trajectoryheight_arr
  END INTERFACE

  INTERFACE rfvolt
     MODULE PROCEDURE rfvolt_r,rfvolt_arr
  END INTERFACE

  INTERFACE realft
     MODULE PROCEDURE realft_sp,realft_dp
  END INTERFACE

  INTERFACE four1
     MODULE PROCEDURE four1_sp,four1_dp
  END INTERFACE

  INTERFACE fourrow
     MODULE PROCEDURE fourrow_sp,fourrow_dp
  END INTERFACE

!-----------------------------------------
! B0: B-field at machine_ref_frame
! Bdot: time derivative of B-field (considered constant)
! Rnom: mean orbit radius
! rhonom: bending radius
! gammatnom: transition gamma
! Erest: rest energy of accelerated particle
! q: charge state of accelerated particle
! c: speed of light in vacuum
! eunit: elementary electric charge
! h: principal harmonic number
! hratio: ratio of harmonics between the two RF systems
! phi12: phase difference between the two RF systems (considered constant)
! self_field_flag: flag to include self-fields in the tracking
! zpu: effective pick-up sensitivity
! gcoupling: space charge coupling coefficient
! zwallovern: magnitude of Zwall/n
! eperimage: total charge in beam_ref_profile
! vself: space charge voltage
! VRF1: peak voltage of first RF system at machine_ref_frame
! VRF2: peak voltage of second RF system at machine_ref_frame
! VRF1dot,VRF2dot: time derivatives of the RF voltages (considered constant)
! dtbin: pixel width
! dEbin: pixel height
! dEmax: maximum energy of reconstructed phase space
! E0: total energy of synchronous particle at each turn
! beta0: Lorenz beta factor (v/c) at each turn
! eta0: phase slip factor at each turn
! omegarev0: revolution frequency at each turn
! tatturn: time at each turn relative to machine_ref_frame
! phi0: synchronous phase angle at each turn
! fitxat0: synchronous phase in bins from foot tangent fit
! proffile: file with 'raw' input frames
! odir: directory into which the ouput is directed
! full_pp_flag: if set, all pixels in reconstructed phase space will be tracked
! wraplength: maximum number of bins to cover an integer number of rf periods
!-----------------------------------------

  REAL(SP), PRIVATE ::  h, Bdot, Rnom, rhonom, gammatnom, Erest, q, eperimage
  REAL(SP), PRIVATE ::  VRF1,VRF2,eunit,fitxat0,bunchphaselength
  REAL(SP), PRIVATE ::  VRF1dot,VRF2dot,zwallovern,gcoupling
  REAL(SP), PRIVATE ::  dtbin,dEbin,energy,B0,c,hratio,phi12,dEmax,zpu
  REAL(SP), PRIVATE ::  tangentfootl,tangentfootu,phiwrap
  REAL(SP), DIMENSION(:), ALLOCATABLE, PRIVATE :: beta0,eta0,&
                          omegarev0,tatturn,phi0,E0,c1,c2,c3
  REAL(SP), DIMENSION(:,:), ALLOCATABLE :: vself
  INTEGER, PRIVATE :: full_pp_flag,wraplength
  CHARACTER*60, PRIVATE :: proffile,odir

CONTAINS

  SUBROUTINE find_xat0(profile)
!   Subroutine to find the synchronous phase from the beam reference profile. 

    IMPLICIT NONE
    REAL(SP), DIMENSION(:), INTENT(IN) :: profile
    INTEGER   i,i01,tangentbinl,tangentbinu,turnnow
    REAL(SP)  al,bl,au,bu,chil,chiu,ql,qu,sigal,sigbl,sigau,sigbu
    REAL(SP)  thold,bunchduration,phil,dradbin
    REAL(SP), DIMENSION(profilelength) :: indarr
    INTEGER, DIMENSION(1) :: maxbin
!
    IF (xat0.LT.0) THEN
      DO i=1,profilelength
        indarr(i)=REAL(i,SP)-0.5_SP
      END DO
      turnnow=(beam_ref_frame-1)*dturns
      thold=0.15_SP*MAXVAL(profile)
      maxbin=MAXLOC(profile)
    
      DO i=maxbin(1),1,-1
        IF (profile(i).LT.thold) THEN
          tangentbinl=i+1
          EXIT
        END IF
      END DO
      DO i=maxbin(1),profilelength,1
        IF (profile(i).LT.thold) THEN
          tangentbinu=i-1
          EXIT
        END IF
      END DO
      CALL fit(indarr(tangentbinl-2:tangentbinl+1),&
               profile(tangentbinl-2:tangentbinl+1),al,bl,sigal,sigbl,chil,ql)
      CALL fit(indarr(tangentbinu-1:tangentbinu+2),&
               profile(tangentbinu-1:tangentbinu+2),au,bu,sigau,sigbu,chiu,qu)
      tangentfootl=-1.0_SP*al/bl
      tangentfootu=-1.0_SP*au/bu
      bunchduration=(tangentfootu-tangentfootl)*dtbin
      bunchphaselength=h*omegarev0(turnnow)*bunchduration
      phil=rtnewt(phi0(turnnow)-bunchphaselength,phi0(turnnow),&
                phi0(turnnow)-bunchphaselength/2_SP,0.0001_SP,phaselow,turnnow)
      fitxat0=tangentfootl+(phi0(turnnow)-phil)/(h*omegarev0(turnnow)*dtbin)
      xat0=fitxat0
    ELSE
      tangentfootl=0.0_SP
      tangentfootu=0.0_SP
      fitxat0=0.0_SP
    END IF
!   Calculate the absolute difference in bins between phase=0 and origin of
!   the reconstructed phase space coordinate system.
    i01=(beam_ref_frame-1)*dturns
    xorigin=phi0(i01)/(h*omegarev0(i01)*dtbin)-xat0
!   Calculate the number of bins in the first integer number of rf periods
!   larger than the image width.
    IF (Bdot .GT. 0.0_sp) THEN
      dradbin=h*omegarev0((profilecount-1)*dturns)*dtbin
    ELSE
      dradbin=h*omegarev0(0)*dtbin
    END IF
    phiwrap=REAL(CEILING(profilelength*dradbin/twopi),sp)*twopi
    wraplength=CEILING(phiwrap/dradbin)
  END SUBROUTINE find_xat0

  SUBROUTINE phaselow(phase,hamil,dhamil,rfturn)
!-----------------------------------------------
!  This function is needed by the Newton-Raphson
!  root finder to calculate fitxat0.
!-----------------------------------------------
    USE nrstuff
    IMPLICIT NONE
    REAL(SP), INTENT(IN) :: phase
    REAL(SP), INTENT(OUT) :: hamil,dhamil
    INTEGER, INTENT(IN) :: rfturn
    hamil=VRF2*(COS(hratio*(phase+bunchphaselength-phi12))-&
                COS(hratio*(phase-phi12)))/hratio+&
          VRF1*(COS(phase+bunchphaselength)-COS(phase))+&
          bunchphaselength*rfvolt(phi0(rfturn),rfturn)
    dhamil=-1.0_SP*VRF2*(SIN(hratio*(phase+bunchphaselength-phi12))-&
           SIN(hratio*(phase-phi12)))-&
           VRF1*(SIN(phase+bunchphaselength)-SIN(phase))
  END SUBROUTINE phaselow

  SUBROUTINE perimage(profile)
    IMPLICIT NONE
    REAL(SP), DIMENSION(:), INTENT(IN) :: profile
!
!   Calculate the total charge in profile.
    eperimage=SUM(profile)*dtbin/(rebin*eunit*zpu)
  END SUBROUTINE perimage

  SUBROUTINE calculate_self(dsprofiles)
    IMPLICIT NONE
    REAL(SP), DIMENSION(:,:), INTENT(IN) :: dsprofiles
    INTEGER p
    ALLOCATE(vself(profilecount-1,wraplength))
!
!   Calculate self-field voltage.
    vself=0.0_sp
    DO p=1,profilecount-1
      vself(p,1:profilelength)=0.5_sp*eperimage*(&
                               c3(p)*dsprofiles(p,1:profilelength)+&
	                       c3(p+1)*dsprofiles(p+1,1:profilelength))
    END DO
    CALL out_profiles(vself(:,1:profilelength),'vself.data')
  END SUBROUTINE calculate_self    

  SUBROUTINE subtract_baseline(indata)
!   Find the baseline from the first 5% percent of the beam reference profile.
    IMPLICIT NONE
    INTEGER barray
    REAL(SP), DIMENSION(:), INTENT(INOUT) :: indata
    REAL(SP) baseline
!
!   Skip some of the frame data to the start of the reference profile.
    barray=(skipcount+beam_ref_frame-1)*framelength + 1 + preskiplength
!
    baseline=SUM(indata(barray:&
                 FLOOR(barray+0.05*REAL(profilelength,SP))))&
             /REAL(FLOOR(0.05*REAL(profilelength,SP)+1.0_SP),SP)
!   Offset by baseline.
    indata=indata-baseline
  END SUBROUTINE subtract_baseline

  SUBROUTINE get_indata(indata)
    USE nrstuff
    IMPLICIT NONE
    REAL(SP), DIMENSION(:), INTENT(INOUT) :: indata
    INTEGER i
    OPEN(1,file=proffile,status='old',access='sequential',&
         form='formatted')
    DO i=1,alldata
      read(1,*) indata(i)  
    END DO
    CLOSE(1)
  END SUBROUTINE get_indata 

  SUBROUTINE get_parameters()
!-----------------------------------------
!  Read input data and initialise parameters.
!-----------------------------------------
    USE nrstuff
    IMPLICIT NONE
    CHARACTER*60 comment
    INTEGER allturn,i,iminskip,imaxskip,turnhere
    OPEN(9,file="input_v2.dat",status='old',&
           access='sequential',form='formatted')
    DO i=1,11
      READ(9,'(a60)') comment
    END DO
!   
    READ(9,'(a60)') comment
    READ(9,'(a60)') proffile
!   
    READ(9,'(a60)') comment
    READ(9,'(a60)') odir
!   
    READ(9,'(a60)') comment
    READ(9,*) framecount    
!   
    READ(9,'(a60)') comment
    READ(9,*) skipcount
!
    READ(9,'(a60)') comment
    READ(9,*) framelength
!
    READ(9,'(a60)') comment
    READ(9,*) dtbin
!
    READ(9,'(a60)') comment
    READ(9,*) dturns
!
    READ(9,'(a60)') comment
    READ(9,*) preskiplength
!
    READ(9,'(a60)') comment
    READ(9,*) postskiplength
!
    READ(9,'(a60)') comment
    READ(9,'(a60)') comment
    READ(9,*) iminskip
!
    READ(9,'(a60)') comment
    READ(9,'(a60)') comment
    READ(9,*) imaxskip
!
    READ(9,'(a60)') comment
    READ(9,*) rebin
!
    READ(9,'(a60)') comment
    READ(9,'(a60)') comment
    READ(9,*) xat0
!
    READ(9,'(a60)') comment
    READ(9,*) dEmax
!
    READ(9,'(a60)') comment
    READ(9,*) filmstart
!
    READ(9,'(a60)') comment
    READ(9,*) filmstop
!
    READ(9,'(a60)') comment
    READ(9,*) filmstep
!
    READ(9,'(a60)') comment
    READ(9,*) niter
!
    READ(9,'(a60)') comment
    READ(9,*) Npt
!
    READ(9,'(a60)') comment
    READ(9,*) full_pp_flag
!
    READ(9,'(a60)') comment
    READ(9,*) beam_ref_frame
!
    READ(9,'(a60)') comment
    READ(9,*) machine_ref_frame
!
    READ(9,'(a60)') comment
    READ(9,'(a60)') comment
!
    READ(9,'(a60)') comment
    READ(9,*) VRF1
!
    READ(9,'(a60)') comment
    READ(9,*) VRF1dot
!
    READ(9,'(a60)') comment
    READ(9,*) VRF2
!
    READ(9,'(a60)') comment
    READ(9,*) VRF2dot
!
    READ(9,'(a60)') comment
    READ(9,*) h
!
    READ(9,'(a60)') comment
    READ(9,*) hratio
!
    READ(9,'(a60)') comment
    READ(9,*) phi12
!
    READ(9,'(a60)') comment
    READ(9,*) B0
!
    READ(9,'(a60)') comment
    READ(9,*) Bdot
!
    READ(9,'(a60)') comment
    READ(9,*) Rnom
!
    READ(9,'(a60)') comment
    READ(9,*) rhonom
!
    READ(9,'(a60)') comment
    READ(9,*) gammatnom
!
    READ(9,'(a60)') comment
    READ(9,*) Erest
!
    READ(9,'(a60)') comment
    READ(9,*) q
!
    READ(9,'(a60)') comment
    READ(9,'(a60)') comment
    READ(9,'(a60)') comment
    READ(9,*) self_field_flag
!
    READ(9,'(a60)') comment
    READ(9,*) gcoupling
!
    READ(9,'(a60)') comment
    READ(9,*) zwallovern
!
    READ(9,'(a60)') comment
    READ(9,*) zpu
!
    CLOSE(9)
!
    allturn=(framecount-skipcount-1)*dturns
    ALLOCATE(tatturn(0:allturn),omegarev0(0:allturn),phi0(0:allturn),&
             c1(0:allturn),c2(0:allturn),&
             beta0(0:allturn),eta0(0:allturn),E0(0:allturn))
!   Speed of light in vacuum.
    c=2.99792458e8_SP
!   Elementary electric charge.
    eunit=1.60217733e-19_sp
!   Lorenz beta (relativistic factor).
    CALL pbttrack(phi0,beta0,tatturn,E0,Erest,allturn)
!   Phase slip factor.
    eta0=(1.0_SP-beta0**2)-gammatnom**(-2)
    c1=twopi*h*eta0/(E0*beta0**2)
!   Revolution frequency.
    omegarev0=beta0*c/Rnom
!   Time binning.
    dtbin=REAL(rebin,SP)*dtbin
    xat0=xat0/REAL(rebin,SP)
    profilecount=framecount-skipcount
!   Here, profilelength is still in frame bins.
    profilelength=framelength-preskiplength-postskiplength
    iminin=iminskip/rebin + 1
    IF (MOD(profilelength-imaxskip,rebin).EQ.0) THEN
      imaxin=(profilelength-imaxskip)/rebin
    ELSE
      imaxin=(profilelength-imaxskip)/rebin + 1
    END IF
    alldata=framecount*framelength
!   Self-field coefficient.
    ALLOCATE(c3(profilecount))
    DO i=1,profilecount
      turnhere=(i-1)*dturns
      c3(i)=(eunit/omegarev0(turnhere))*&
         ((1_sp/beta0(turnhere)-beta0(turnhere))*&
         gcoupling*pi*2.e-7_sp*c-zwallovern)/dtbin**2
    END DO
  END SUBROUTINE get_parameters

  SUBROUTINE pbttrack(phi0,beta0,tatturn,E0,Erest,allturn)
!-----------------------------------------
!   Subroutine which calculates synchronous phase, Lorenz beta, energy and time
!   for the synchronous particle as a function of turn. Values are calculated
!   immediately after the 'single' cavity of the ring.
!
!   Erest: rest energy of accelerated particle
!   E0: total energy at the end of the turn
!   beta0: Lorenz beta at the end of the turn
!   phi0: synchronous phase at the end of the turn
!   tatturn: time relative to machine_ref_frame at the end of the turn
!   allturn: total number of turns
!   i,i0: integer variables for loops
!-----------------------------------------
    IMPLICIT NONE
    REAL(SP), INTENT(IN) :: Erest
    REAL(SP) phistart,phil,phiu
    INTEGER allturn,i,i0
    REAL(SP), DIMENSION(0:allturn), INTENT(OUT) :: phi0,beta0,tatturn,E0
    i0=(machine_ref_frame-1)*dturns
    tatturn(i0)=0
    E0(i0)=BtoEnom(B0)
    beta0(i0)=SQRT(1.0_SP-(Erest/E0(i0))**2)
    IF (q*(gammatnom-E0(i0)/Erest)>0) THEN
      phil=-1.0_SP*pi
      phiu=pi
    ELSE
      phil=0.0_SP
      phiu=twopi
    END IF
    phistart=rtnewt(phil,phiu,(phil+phiu)/2.0_SP,0.001_SP,rfvrf1only,i0)
!   Synchronous phase of a particle on the nominal orbit.
    phi0(i0)=rtnewt(phil,phiu,phistart,0.001_SP,rfvanddv,i0)
    DO i=i0+1,allturn
      tatturn(i)=tatturn(i-1)+twopi*Rnom/(beta0(i-1)*c)
      phi0(i)=rtnewt(phil,phiu,phi0(i-1),0.001_SP,rfvanddv,i)
      E0(i)=E0(i-1)+q*rfvolt(phi0(i),i)
      beta0(i)=SQRT(1.0_SP-(Erest/E0(i))**2)
      c2(i)=E0(i)-E0(i-1)
    END DO
    DO i=i0-1,0,-1
      E0(i)=E0(i+1)-q*rfvolt(phi0(i+1),i+1)
      beta0(i)=SQRT(1.0_SP-(Erest/E0(i))**2)
      c2(i)=E0(i+1)-E0(i)
      tatturn(i)=tatturn(i+1)-twopi*Rnom/(beta0(i)*c)
      phi0(i)=rtnewt(phil,phiu,phi0(i+1),0.001_SP,rfvanddv,i)
    END DO
  END SUBROUTINE pbttrack

  SUBROUTINE out_plotinfo
!-----------------------------------------
!  Writes data needed for plots on a file.
!-----------------------------------------
    USE nrstuff
    IMPLICIT NONE
    INTEGER p
    OPEN(8,file=TRIM(odir)//'plotinfo.data',status='new',access='sequential',&
         form='formatted')
    write(8,'(a)') 'Number of profiles used in each reconstruction,'
    write(8,'(a,I3)') ' profilecount = ',profilecount
    write(8,'(a)') 'Width (in pixels) of each image = length (in bins) of each profile,'
    write(8,'(a,I3)') ' profilelength = ',profilelength
    write(8,'(a)') 'Width (in s) of each pixel = width of each profile bin,'
    write(8,'(a,1P,E11.4)') ' dtbin =',dtbin
    write(8,'(a)') 'Height (in eV) of each pixel,'
    write(8,'(a,1P,E11.4)') ' dEbin =',dEbin
    write(8,'(a)') 'Number of elementary charges in each image,'
    write(8,'(a,1P,E10.3)') ' eperimage =',eperimage
    write(8,'(a)') 'Position (in pixels) of the reference synchronous point:'
    write(8,'(a,F8.3)') ' xat0 =',xat0
    write(8,'(a,F8.3)') ' yat0 =',yat0
    write(8,'(a)') 'Foot tangent fit results (in bins):'
    write(8,'(a,F8.3)') ' tangentfootl = ',tangentfootl
    write(8,'(a,F8.3)') ' tangentfootu = ',tangentfootu
    write(8,'(a,F8.3)') ' fit xat0 =',fitxat0
    write(8,'(a)') 'Synchronous phase (in radians):'
    DO p=filmstart,filmstop,filmstep
      write(8,'(a,I3,a,F7.4)') ' phi0(',p,') =',phi0((p-1)*dturns)
    END DO
    write(8,'(a)') 'Horizontal range (in pixels) of the region in phase space of map elements:'
    DO p=filmstart,filmstop,filmstep
      write(8,'(a,i3,a,i3,a,i3,a,i3)') ' imin(',p,') = ',imin(p),&
          ' and imax(',p,') = ',imax(p)
    END DO
    CLOSE(8)
  END SUBROUTINE out_plotinfo

  SUBROUTINE out_profiles(profiles,fname)
!-----------------------------------------
!  Writes profiles to a file.
!-----------------------------------------
    USE nrstuff
    IMPLICIT NONE
    INTEGER P,I
    CHARACTER*(*) :: fname
    REAL(SP), DIMENSION(:,:), INTENT(IN) :: profiles
    OPEN(8,file=TRIM(odir)//fname,status='new',access='sequential',&
           form='formatted')
    DO P=1,SIZE(profiles)/profilelength
      DO I=1,profilelength
        write(8,'(E14.7)') profiles(P,I)
      END DO
    END DO
    CLOSE(8)
  END SUBROUTINE out_profiles

  SUBROUTINE out_picture(phasespace,fname)
!-----------------------------------------
!  Writes phasepace to a file.
!-----------------------------------------
   USE nrstuff
   IMPLICIT NONE
   REAL(SP), DIMENSION(:,:), INTENT(IN)  :: phasespace
   INTEGER i,j
   CHARACTER*(*) :: fname
   OPEN(8,file=TRIM(odir)//fname,status='new',access='sequential',&
        form='formatted')
   DO i=1,profilelength
     DO j=1,profilelength
       write(8,*) phasespace(i,j)
     END DO
   END DO
   CLOSE(8)
  END SUBROUTINE out_picture

  SUBROUTINE out_array(onearray,fname)
!-----------------------------------------
!  Writes any array to a file.
!-----------------------------------------
   USE nrstuff
   IMPLICIT NONE
   REAL(SP), DIMENSION(:), INTENT(IN)  :: onearray
   INTEGER i
   CHARACTER*(*) :: fname
   CLOSE(8)
   OPEN(8,file=TRIM(odir)//fname,status='new',access='sequential',&
        form='formatted')
   DO i=1,SIZE(onearray)
     write(8,*) i-1,onearray(i)
   END DO
   CLOSE(8)
  END SUBROUTINE out_array

  SUBROUTINE ijlimits
!-----------------------------------------
!  This procedure calculates and sets the limits in i (phase) and j (energy).
!  The i-j coordinate system is the one used locally (reconstructed phase space).
!
!  LOCAL VARIABLES:
!  turnnow: machine turn (=0 at profile=1)
!  jmaxul: array holding calculated max values (left and right of bin)
!  phases: phase at the edge of each bin along the i-axis
!  energiesl,energiesu: maximum energy starting from lowest phase and 
!                       uppermost phase respectively
!  indarray: local array holding real numbers from 1 to profilelength
!-----------------------------------------
    USE nrstuff
    IMPLICIT NONE
    INTEGER i,p,turnnow
    TYPE(jlim), DIMENSION(profilelength+1) :: jmaxul
    REAL(SP),DIMENSION(profilelength+1)::phases,energiesl,energiesu
    REAL(SP), DIMENSION(profilelength) :: indarr
    turnnow=(beam_ref_frame-1)*dturns
    DO i=1,profilelength
      indarr(i)=REAL(i,SP)
    END DO
    phases(1)=xorigin*dtbin*h*omegarev0(turnnow)
    phases(2:profilelength+1)=(xorigin+indarr)*dtbin*h*omegarev0(turnnow)
!   Energy binning.
    IF (dEmax.LT.0.0_SP) THEN
      IF (VRFt(VRF2,VRF2dot,turnnow).NE.0.0_SP) THEN
        energiesl=trajectoryheight(phases,phases(1),0.0_SP,turnnow)
        energiesu=trajectoryheight(phases,phases(profilelength+1),&
                  0.0_SP,turnnow)
        dEbin=MIN(MAXVAL(energiesl),MAXVAL(energiesu))/(profilelength-yat0)
      ELSE
        dEbin=beta0(turnnow)*SQRT(E0(turnnow)*q*VRFt(VRF1,VRF1dot,turnnow)*&
         COS(phi0(turnnow))/(twopi*h*eta0(turnnow)))*dtbin*h*omegarev0(turnnow)
      END IF
    ELSE
      dEbin=dEmax/(profilelength-yat0)
    ENDIF
!   Extrema of active pixels in energy (j-index or y-index) direction.
    energy=0.0_SP
    IF (full_pp_flag.EQ.1) THEN
      jmax=profilelength
      jmin=CEILING(2.0_SP*yat0-jmax+0.5_SP)
      allbinmin=1
      allbinmax=profilelength
    ELSE 
      DO p=filmstart,filmstop,filmstep
        turnnow=(p-1)*dturns
        phases(1)=xorigin*dtbin*h*omegarev0(turnnow)
        phases(2:profilelength+1)=(xorigin+indarr)*dtbin*h*omegarev0(turnnow)
        DO i=1,profilelength+1
          jmaxul(i)%l=FLOOR(yat0+trajectoryheight(phases(i),phases(1),&
                        energy,turnnow)/dEbin)
          jmaxul(i)%u=FLOOR(yat0+trajectoryheight(phases(i),&
                        phases(profilelength+1),energy,turnnow)/dEbin)
        END DO
        DO i=1,profilelength
          jmax(p,i)=MIN(jmaxul(i)%u,jmaxul(i+1)%u,jmaxul(i)%l,jmaxul(i+1)%l,&
                profilelength)
        END DO
        jmin(p,:)=MAX(1,CEILING(2.0_SP*yat0-jmax(p,:)+0.5_SP))
        DO i=1,profilelength
          IF ((jmax(p,i)-jmin(p,i)).GE.0) THEN
            allbinmin(p)=i
            EXIT 
          END IF
        END DO
        DO i=profilelength,1,-1
          IF ((jmax(p,i)-jmin(p,i)).GE.0) THEN
            allbinmax(p)=i
            EXIT
          END IF
        END DO
      END DO
    END IF
    OPEN(8,file=TRIM(odir)//'jmax.data',status='new',&
         access='sequential',form='formatted')
    DO p=filmstart,filmstop,filmstep
      imin(p)=allbinmin(p)
      IF (iminin.GT.allbinmin(p) .OR. full_pp_flag.EQ.1) THEN
        imin(p)=iminin
        jmax(p,1:iminin-1)=FLOOR(yat0)
        jmin(p,:)=CEILING(2.0_SP*yat0-jmax(p,:)+0.5_SP)
      END IF
      imax(p)=allbinmax(p)
      IF (imaxin.LT.allbinmax(p) .OR. full_pp_flag.EQ.1) THEN
        imax(p)=imaxin
        jmax(p,imaxin+1:profilelength)=FLOOR(yat0)
        jmin(p,:)=CEILING(2.0_SP*yat0-jmax(p,:)+0.5_SP)
      END IF
      DO i=1,profilelength
        write(8,*) i,jmax(p,i)
      END DO
    END DO 
    CLOSE(8)
    DEALLOCATE(E0,eta0,beta0)
  END SUBROUTINE ijlimits

  FUNCTION trajectoryheight_var(phi,phiknown,deltaEknown,turnnow)
!-----------------------------------------
!  This is a trajectory height calculator
!  given a phase and energy.
!-----------------------------------------
    USE nrstuff
    IMPLICIT NONE
    COMPLEX complexheight
    REAL(SP) trajectoryheight_var,phi,phiknown,deltaEknown,aa,bb,cc,dd
    INTEGER turnnow
    aa=deltaEknown**2
    bb=2.0_SP*q/c1(turnnow)
    cc=VRFt(VRF1,VRF1dot,turnnow)*(COS(phi)-COS(phiknown))+&
       VRFt(VRF2,VRF2dot,turnnow)*(COS(hratio*(phi-phi12))-&
       COS(hratio*(phiknown-phi12)))/hratio+(phi-phiknown)*&
       rfvolt(phi0(turnnow),turnnow)
    dd=aa+bb*cc
    complexheight=SQRT(CMPLX(dd))
    trajectoryheight_var=REAL(complexheight,SP)
  END FUNCTION

  FUNCTION trajectoryheight_arr(phi,phiknown,deltaEknown,turnnow)
!-----------------------------------------
!  This is a trajectory height calculator
!  given a phase and energy.
!-----------------------------------------
    USE nrstuff
    IMPLICIT NONE
    REAL(SP) phiknown,deltaEknown,aa,bb
    INTEGER turnnow
    REAL(SP), DIMENSION(:) :: phi
    REAL(SP),DIMENSION(SIZE(phi))::trajectoryheight_arr,cc,dd
    COMPLEX, DIMENSION(SIZE(phi))::complexheight
    aa=deltaEknown**2
    bb=2.0_SP*q/c1(turnnow)
    cc=VRFt(VRF1,VRF1dot,turnnow)*(COS(phi)-COS(phiknown))+&
       VRFt(VRF2,VRF2dot,turnnow)*(COS(hratio*(phi-phi12))-&
       COS(hratio*(phiknown-phi12)))/hratio+(phi-phiknown)*&
       rfvolt(phi0(turnnow),turnnow)
    dd=aa+bb*cc
    complexheight=SQRT(CMPLX(dd))
    trajectoryheight_arr=REAL(complexheight,SP)
  END FUNCTION

  FUNCTION BtoEnom(B)
!-----------------------------------------
!  Calculates the energy for a particle in
!  a circular machine at dipole field B.
!-----------------------------------------
    USE nrstuff
    IMPLICIT NONE
    REAL(SP) BtoEnom,B
    BtoEnom=SQRT((q*B*rhonom*c)**2 + Erest**2)
  END FUNCTION

SUBROUTINE rfvanddv(phi,rfv,drfv,rfturn)
!-----------------------------------------
!  The rfvolt function and its derivative.
!  This function is needed by the Newton-Raphson
!  root finder to calculate phi0.
!-----------------------------------------
  USE nrstuff
  IMPLICIT NONE
  REAL(SP), INTENT(IN) :: phi
  REAL(SP), INTENT(OUT) :: rfv,drfv
  INTEGER, INTENT(IN) :: rfturn
  REAL(SP) V1,V2
  V1=VRFt(VRF1,VRF1dot,rfturn)
  V2=VRFt(VRF2,VRF2dot,rfturn)
  rfv=V1*SIN(phi)+V2*SIN(hratio*(phi-phi12))-twopi*Rnom*rhonom*Bdot*SIGN(1.0_SP,q)
  drfv=V1*COS(phi)+hratio*V2*COS(hratio*(phi-phi12))
END SUBROUTINE rfvanddv

SUBROUTINE rfvrf1only(phi,rfv,drfv,rfturn)
!-----------------------------------------
!  The rfvolt function and its derivative.
!  This function is needed by the Newton-Raphson
!  root finder to calculate phi0.
!-----------------------------------------
  USE nrstuff
  IMPLICIT NONE
  REAL(SP), INTENT(IN) :: phi
  REAL(SP), INTENT(OUT) :: rfv,drfv
  INTEGER, INTENT(IN) :: rfturn
  REAL(SP) V1
  V1=VRFt(VRF1,VRF1dot,rfturn)
  rfv=V1*SIN(phi)-twopi*Rnom*rhonom*Bdot*SIGN(1.0_SP,q)
  drfv=V1*COS(phi)
END SUBROUTINE rfvrf1only

SUBROUTINE longtrack(direction,nrep,yp,xp,turnnow)
!---------------------------------------------------------------------------
! h: principal harmonic number
! eta0: phase slip factor
! E0: energy of synchronous particle 
! beta0: relativistic beta of synchronous particle
! phi0: synchronous phase
! q: charge state of particles
! dphi: phase difference between considered particle and synchronous one
! denergy: energy difference between considered particle and synchronous one
! nrep: pass cavity nrep times before returning data
! direction: to inverse the time advance (rotation in the bucket), 1 or -1
! xp and yp: time and energy in pixels
! dtbin and dEbin: GLOBAL time and energy pixel size in s and MeV
! omegarev0: revolution frequency
! VRF1,VRF2,VRF1dot,VRF2dot: GLOBAL RF voltages and derivatives of volts
! turnnow: present turn
!---------------------------------------------------------------------------
  USE nrstuff
  IMPLICIT NONE
  REAL(SP), DIMENSION(:), INTENT(INOUT) :: xp,yp
  REAL(SP), DIMENSION(SIZE(xp)) :: dphi,denergy
!HPF$ distribute dphi(block)
!HPF$ align with dphi :: denergy
integer :: mm
  INTEGER :: i,nrep,direction,turnnow
  dphi=(xp+xorigin)*h*omegarev0(turnnow)*dtbin-phi0(turnnow)
  denergy=(yp-yat0)*dEbin
  IF (direction.GT.0) THEN
    DO i=1,nrep
      forall(mm=1:size(xp)) dphi(mm)=dphi(mm)-c1(turnnow)*denergy(mm)
      turnnow=turnnow+1
      forall(mm=1:size(xp)) denergy(mm)=denergy(mm)+q*(&
	(VRF1+VRF1dot*tatturn(turnnow))*SIN(dphi(mm)+phi0(turnnow))+&
        (VRF2+VRF2dot*tatturn(turnnow))*&
	SIN(hratio*(dphi(mm)+phi0(turnnow)-phi12)))-c2(turnnow)
    END DO
  ELSE
    DO i=1,nrep
      forall(mm=1:size(xp)) denergy(mm)=denergy(mm)-q*(&
	(VRF1+VRF1dot*tatturn(turnnow))*SIN(dphi(mm)+phi0(turnnow))+&
        (VRF2+VRF2dot*tatturn(turnnow))*&
	SIN(hratio*(dphi(mm)+phi0(turnnow)-phi12)))+c2(turnnow)
      turnnow=turnnow-1
      forall(mm=1:size(xp)) dphi(mm)=dphi(mm)+c1(turnnow)*denergy(mm)
    END DO  
  END IF
  xp=(dphi+phi0(turnnow))/(h*omegarev0(turnnow)*dtbin)-xorigin
  yp=denergy/dEbin+yat0
END SUBROUTINE longtrack

SUBROUTINE longtrack_self(direction,nrep,yp,xp,turnnow)
!-------------------------------------------------------------------------
! h: principal harmonic number
! eta0: phase slip factor
! E0: energy of synchronous particle 
! beta0: relativistic beta of synchronous particle
! phi0: synchronous phase
! q: charge state of particles
! dphi: phase difference between considered particle and synchronous one
! denergy: energy difference between considered particle and synchronous one
! nrep: pass cavity nrep times before returning data
! direction: to inverse the time advance (rotation in the bucket), 1 or -1
! xp and yp: time and energy in pixels
! dtbin and dEbin: GLOBAL time and energy pixel size in s and MeV
! omegarev0: revolution frequency
! VRF1,VRF2,VRF1dot,VRF2dot: GLOBAL RF voltages and derivatives of volts
! turnnow: present turn
!---------------------------------------------------------------------------
  USE nrstuff
  IMPLICIT NONE
  REAL(SP), DIMENSION(:), INTENT(INOUT) :: xp,yp
  REAL(SP), DIMENSION(SIZE(xp)) :: dphi,denergy,selfvolt
!HPF$ distribute dphi(block)
!HPF$ align with dphi :: denergy,selfvolt,xp
integer :: mm
  INTEGER :: i,j,p,nrep,direction,turnnow,xx
  dphi=(xp+xorigin)*h*omegarev0(turnnow)*dtbin-phi0(turnnow)
  denergy=(yp-yat0)*dEbin
  IF (direction.GT.0) THEN
    p=turnnow/dturns+1
    DO i=1,nrep
      forall(mm=1:size(xp)) dphi(mm)=dphi(mm)-c1(turnnow)*denergy(mm)
      turnnow=turnnow+1
      forall(mm=1:size(xp)) xp(mm)=dphi(mm)+phi0(turnnow)-&
                                   xorigin*h*omegarev0(turnnow)*dtbin
      forall(mm=1:size(xp)) xp(mm)=(xp(mm)-&
        phiwrap*FLOOR(xp(mm)/phiwrap))/(h*omegarev0(turnnow)*dtbin)
      forall(mm=1:size(xp)) selfvolt(mm)=vself(p,FLOOR(xp(mm))+1)
      forall(mm=1:size(xp)) denergy(mm)=denergy(mm)+q*((&
        (VRF1+VRF1dot*tatturn(turnnow))*SIN(dphi(mm)+phi0(turnnow))+&
	(VRF2+VRF2dot*tatturn(turnnow))*&
        SIN(hratio*(dphi(mm)+phi0(turnnow)-phi12)))+selfvolt(mm))-c2(turnnow)
    END DO
  ELSE
    p=turnnow/dturns
    DO i=1,nrep
      forall(mm=1:size(xp)) selfvolt(mm)=vself(p,FLOOR(xp(mm))+1)
      forall(mm=1:size(xp)) denergy(mm)=denergy(mm)-q*((&
        (VRF1+VRF1dot*tatturn(turnnow))*SIN(dphi(mm)+phi0(turnnow))+&
	(VRF2+VRF2dot*tatturn(turnnow))*&
        SIN(hratio*(dphi(mm)+phi0(turnnow)-phi12)))+selfvolt(mm))+c2(turnnow)
      turnnow=turnnow-1
      forall(mm=1:size(xp)) dphi(mm)=dphi(mm)+c1(turnnow)*denergy(mm)
      forall(mm=1:size(xp)) xp(mm)=dphi(mm)+phi0(turnnow)-&
                                   xorigin*h*omegarev0(turnnow)*dtbin
      forall(mm=1:size(xp)) xp(mm)=(xp(mm)-&
        phiwrap*FLOOR(xp(mm)/phiwrap))/(h*omegarev0(turnnow)*dtbin)
    END DO
  END IF
  yp=denergy/dEbin+yat0
END SUBROUTINE longtrack_self

FUNCTION rfvolt_r(phi,rfturn)
!-----------------------------------------
!  For a new RF function change rfvolt here,
!  but pay attention for optmization purpose
!  rfvolt is implemeted inline as well.
!-----------------------------------------
  USE nrstuff
  IMPLICIT NONE
  REAL(SP) phi,rfvolt_r
  INTEGER rfturn
  rfvolt_r=VRFt(VRF1,VRF1dot,rfturn)*SIN(phi)+&
           VRFt(VRF2,VRF2dot,rfturn)*SIN(hratio*(phi-phi12))
END FUNCTION

FUNCTION rfvolt_arr(phi,rfturn)
! For a new RF function change rfvolt here.
  USE nrstuff
  IMPLICIT NONE
  INTEGER rfturn
  REAL(SP), DIMENSION(:) :: phi
  REAL(SP), DIMENSION(SIZE(phi)) :: rfvolt_arr
  rfvolt_arr=VRFt(VRF1,VRF1dot,rfturn)*SIN(phi)+&
             VRFt(VRF2,VRF2dot,rfturn)*SIN(hratio*(phi-phi12))
END FUNCTION

FUNCTION rtnewt(x1,x2,xstart,xacc,funcd,rfturn)
!--------------------------------------------
!  This is a Newton-Raphson root finder.
!  Modified from Numerical Recipies to start at
!  xstart and calculate rfvolt at rfturn.
!--------------------------------------------
        USE nrstuff
        IMPLICIT NONE
        REAL(SP) x1,x2,xstart,xacc
        REAL(SP) :: rtnewt
        INTEGER rfturn
        INTERFACE
                SUBROUTINE funcd(x,fval,fderiv,rfturn)
                USE nrstuff
                IMPLICIT NONE
                REAL(SP), INTENT(IN) :: x
                REAL(SP), INTENT(OUT) :: fval,fderiv
                INTEGER, INTENT(IN) :: rfturn
                END SUBROUTINE funcd
        END INTERFACE
        INTEGER(I4B), PARAMETER :: MAXIT=20
        INTEGER(I4B) :: j
        REAL(SP) :: df,dx,f
        rtnewt=xstart
        do j=1,MAXIT
                call funcd(rtnewt,f,df,rfturn)
                dx=f/df
                rtnewt=rtnewt-dx
                if ((x1-rtnewt)*(rtnewt-x2) < 0.0)&
                        call nrerror('rtnewt: values jumped out of brackets')
                if (abs(dx) < xacc) RETURN
        end do
        call nrerror('rtnewt exceeded maximum iterations')
END FUNCTION rtnewt

FUNCTION VRFt(VRF,VRFdot,rfturn)
!--------------------------------------------
!  This function calculates the RF peak voltage at turn rfturn
!  assuming a linear voltage function VRFt=VRFdot*time+VRF.
!  time=0 at machine_ref_frame.
!
!  rfturn: turn for which the RF voltage should be calculated
!  VRF: Volts
!  VRFdot: Volts/s
!--------------------------------------------
  REAL(SP) VRFt,VRF,VRFdot
  INTEGER rfturn
  VRFt=VRF+VRFdot*tatturn(rfturn)
END FUNCTION VRFt

  SUBROUTINE fit(x,y,a,b,siga,sigb,chi2,q,sig)
    USE nrstuff
    IMPLICIT NONE
    REAL(SP), DIMENSION(:), INTENT(IN) :: x,y
    REAL(SP), INTENT(OUT) :: a,b,siga,sigb,chi2,q
    REAL(SP), DIMENSION(:), OPTIONAL, INTENT(IN) :: sig
    INTEGER(I4B) :: ndum
    REAL(SP) :: sigdat,ss,sx,sxoss,sy,st2
    REAL(SP), DIMENSION(size(x)), TARGET :: t
    REAL(SP), DIMENSION(:), POINTER :: wt
    if (present(sig)) then
      ndum=assert_eq(size(x),size(y),size(sig),'fit')
      wt=>t
      wt(:)=1.0_sp/(sig(:)**2)
      ss=sum(wt(:))
      sx=dot_product(wt,x)
      sy=dot_product(wt,y)
    else
      ndum=assert_eq(size(x),size(y),'fit')
      ss=real(size(x),sp)
      sx=sum(x)
      sy=sum(y)
    end if
    sxoss=sx/ss
    t(:)=x(:)-sxoss
    if (present(sig)) then
      t(:)=t(:)/sig(:)
      b=dot_product(t/sig,y)
    else
      b=dot_product(t,y)
    end if
    st2=dot_product(t,t)
    b=b/st2
    a=(sy-sx*b)/ss
    siga=sqrt((1.0_sp+sx*sx/(ss*st2))/ss)
    sigb=sqrt(1.0_sp/st2)
    t(:)=y(:)-a-b*x(:)
    if (present(sig)) then
      t(:)=t(:)/sig(:)
      chi2=dot_product(t,t)
      q=gammq(0.5_sp*(size(x)-2),0.5_sp*chi2)
    else
      chi2=dot_product(t,t)
      q=1.0
      sigdat=sqrt(chi2/(size(x)-2))
      siga=siga*sigdat
      sigb=sigb*sigdat
    end if
  END SUBROUTINE fit

  FUNCTION savgol(nl,nr,ld,m)
    IMPLICIT NONE
    INTEGER(I4B), INTENT(IN) :: nl,nr,ld,m
    REAL(SP), DIMENSION(nl+nr+1) :: savgol
    INTEGER(I4B) :: imj,ipj,mm,np
    INTEGER(I4B), DIMENSION(m+1) :: indx
    REAL(SP) :: d,sm
    REAL(SP), DIMENSION(m+1) :: b
    REAL(SP), DIMENSION(m+1,m+1) :: a
    INTEGER(I4B) :: irng(nl+nr+1)
    call assert(nl >= 0, nr >= 0, ld <= m, nl+nr >= m, 'Bad savgol args!')
    do ipj=0,2*m
      sm=sum(arth(1.0_sp,1.0_sp,nr)**ipj)+&
         sum(arth(-1.0_sp,-1.0_sp,nl)**ipj)
      if (ipj == 0) sm=sm+1.0_sp
      mm=min(ipj,2*m-ipj)
      do imj=-mm,mm,2
        a(1+(ipj+imj)/2,1+(ipj-imj)/2)=sm
      end do
    end do
    call ludcmp(a(:,:),indx(:),d)
    b(:)=0.0
    b(ld+1)=1.0
    call lubksb(a(:,:),indx(:),b(:))
    savgol(:)=0.0
    np=nl+nr+1
    irng(:)=arth(-nl,1,np)
    savgol(mod(np-irng(:),np)+1)=poly(real(irng(:),sp),b(:))
  END FUNCTION savgol

  SUBROUTINE lubksb(a,indx,b)
    IMPLICIT NONE
    REAL(SP), DIMENSION(:,:), INTENT(IN) :: a
    INTEGER(I4B), DIMENSION(:), INTENT(IN) :: indx
    REAL(SP), DIMENSION(:), INTENT(INOUT) :: b
    INTEGER(I4B) :: i,n,ii,ll
    REAL(SP) :: summ
    n=assert_eq(size(a,1),size(a,2),size(indx),'lubksb')
    ii=0
    do i=1,n
      ll=indx(i)
      summ=b(ll)
      b(ll)=b(i)
      if (ii /= 0) then
        summ=summ-dot_product(a(i,ii:i-1),b(ii:i-1))
      else if (summ /= 0.0) then
        ii=i
      end if
      b(i)=summ
    end do
    do i=n,1,-1
      b(i) = (b(i)-dot_product(a(i,i+1:n),b(i+1:n)))/a(i,i)
    end do
  END SUBROUTINE lubksb

  SUBROUTINE ludcmp(a,indx,d)
    IMPLICIT NONE
    REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
    INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: indx
    REAL(SP), INTENT(OUT) :: d
    REAL(SP), DIMENSION(size(a,1)) :: vv
    REAL(SP), PARAMETER :: TINY=1.0e-20_sp
    INTEGER(I4B) :: j,n,imax
    n=assert_eq(size(a,1),size(a,2),size(indx),'ludcmp')
    d=1.0
    vv=maxval(abs(a),dim=2)
    if (any(vv == 0.0)) call nrerror('Singular matrix in ludcmp!')
    vv=1.0_sp/vv
    do j=1,n
        imax=(j-1)+imaxloc(vv(j:n)*abs(a(j:n,j)))
        if (j /= imax) then
            call swap(a(imax,:),a(j,:))
            d=-d
            vv(imax)=vv(j)
        end if
        indx(j)=imax
        if (a(j,j) == 0.0) a(j,j)=TINY
        a(j+1:n,j)=a(j+1:n,j)/a(j,j)
        a(j+1:n,j+1:n)=a(j+1:n,j+1:n)-outerprod(a(j+1:n,j),a(j,j+1:n))
    end do
  END SUBROUTINE ludcmp

  FUNCTION convlv(data,respns,isign)
    IMPLICIT NONE
    REAL(SP), DIMENSION(:), INTENT(INOUT) :: data
    REAL(SP), DIMENSION(:), INTENT(IN) :: respns
    INTEGER(I4B), INTENT(IN) :: isign
    REAL(SP), DIMENSION(size(data)) :: convlv
    INTEGER(I4B) :: no2,n,m
    COMPLEX(SPC), DIMENSION(size(data)/2) :: tmpd,tmpr
    n=size(data)
    m=size(respns)
    call assert(iand(n,n-1)==0, 'n must be a power of 2 in convlv!')
    call assert(mod(m,2)==1, 'm must be odd in convlv!')
    convlv(1:m)=respns(:)
    convlv(n-(m-3)/2:n)=convlv((m+3)/2:m)
    convlv((m+3)/2:n-(m-1)/2)=0.0
    no2=n/2
    call realft(data,1,tmpd)
    call realft(convlv,1,tmpr)
    if (isign == 1) then
        tmpr(1)=cmplx(real(tmpd(1))*real(tmpr(1))/no2, &
            aimag(tmpd(1))*aimag(tmpr(1))/no2, kind=spc)
        tmpr(2:)=tmpd(2:)*tmpr(2:)/no2
    else if (isign == -1) then
        if (any(abs(tmpr(2:)) == 0.0) .or. real(tmpr(1)) == 0.0 &
            .or. aimag(tmpr(1)) == 0.0) call nrerror &
            ('deconvolving at response zero in convlv')
        tmpr(1)=cmplx(real(tmpd(1))/real(tmpr(1))/no2, &
            aimag(tmpd(1))/aimag(tmpr(1))/no2, kind=spc)
        tmpr(2:)=tmpd(2:)/tmpr(2:)/no2
    else
        call nrerror('No meaning for isign in convlv!')
    end if
    call realft(convlv,-1,tmpr)
  END FUNCTION convlv

  SUBROUTINE realft_sp(data,isign,zdata)
    IMPLICIT NONE
    REAL(SP), DIMENSION(:), INTENT(INOUT) :: data
    INTEGER(I4B), INTENT(IN) :: isign
    COMPLEX(SPC), DIMENSION(:), OPTIONAL, TARGET :: zdata
    INTEGER(I4B) :: n,ndum,nh,nq
    COMPLEX(SPC), DIMENSION(size(data)/4) :: w
    COMPLEX(SPC), DIMENSION(size(data)/4-1) :: h1,h2
    COMPLEX(SPC), DIMENSION(:), POINTER :: cdata
    COMPLEX(SPC) :: z
    REAL(SP) :: c1=0.5_sp,c2
    n=size(data)
    call assert(iand(n,n-1)==0, 'n must be a power of 2 in realft_sp!')
    nh=n/2
    nq=n/4
    if (present(zdata)) then
        ndum=assert_eq(n/2,size(zdata),'realft_sp')
        cdata=>zdata
        if (isign == 1) cdata=cmplx(data(1:n-1:2),data(2:n:2),kind=spc)
    else
        allocate(cdata(n/2))
        cdata=cmplx(data(1:n-1:2),data(2:n:2),kind=spc)
    end if
    if (isign == 1) then
        c2=-0.5_sp
        call four1(cdata,+1)
    else
        c2=0.5_sp
    end if
    w=zroots_unity(sign(n,isign),n/4)
    w=cmplx(-aimag(w),real(w),kind=spc)
    h1=c1*(cdata(2:nq)+conjg(cdata(nh:nq+2:-1)))
    h2=c2*(cdata(2:nq)-conjg(cdata(nh:nq+2:-1)))
    cdata(2:nq)=h1+w(2:nq)*h2
    cdata(nh:nq+2:-1)=conjg(h1-w(2:nq)*h2)
    z=cdata(1)
    if (isign == 1) then
        cdata(1)=cmplx(real(z)+aimag(z),real(z)-aimag(z),kind=spc)
    else
        cdata(1)=cmplx(c1*(real(z)+aimag(z)),c1*(real(z)-aimag(z)),kind=spc)
        call four1(cdata,-1)
    end if
    if (present(zdata)) then
        if (isign /= 1) then
            data(1:n-1:2)=real(cdata)
            data(2:n:2)=aimag(cdata)
        end if
    else
        data(1:n-1:2)=real(cdata)
        data(2:n:2)=aimag(cdata)
        deallocate(cdata)
    end if
   END SUBROUTINE realft_sp


   SUBROUTINE realft_dp(data,isign,zdata)
    IMPLICIT NONE
    REAL(DP), DIMENSION(:), INTENT(INOUT) :: data
    INTEGER(I4B), INTENT(IN) :: isign
    COMPLEX(DPC), DIMENSION(:), OPTIONAL, TARGET :: zdata
    INTEGER(I4B) :: n,ndum,nh,nq
    COMPLEX(DPC), DIMENSION(size(data)/4) :: w
    COMPLEX(DPC), DIMENSION(size(data)/4-1) :: h1,h2
    COMPLEX(DPC), DIMENSION(:), POINTER :: cdata
    COMPLEX(DPC) :: z
    REAL(DP) :: c1=0.5_dp,c2
    n=size(data)
    call assert(iand(n,n-1)==0, 'n must be a power of 2 in realft_dp!')
    nh=n/2
    nq=n/4
    if (present(zdata)) then
        ndum=assert_eq(n/2,size(zdata),'realft_dp')
        cdata=>zdata
        if (isign == 1) cdata=cmplx(data(1:n-1:2),data(2:n:2),kind=spc)
    else
        allocate(cdata(n/2))
        cdata=cmplx(data(1:n-1:2),data(2:n:2),kind=spc)
    end if
    if (isign == 1) then
        c2=-0.5_dp
        call four1(cdata,+1)
    else
        c2=0.5_dp
    end if
    w=zroots_unity(sign(n,isign),n/4)
    w=cmplx(-aimag(w),real(w),kind=dpc)
    h1=c1*(cdata(2:nq)+conjg(cdata(nh:nq+2:-1)))
    h2=c2*(cdata(2:nq)-conjg(cdata(nh:nq+2:-1)))
    cdata(2:nq)=h1+w(2:nq)*h2
    cdata(nh:nq+2:-1)=conjg(h1-w(2:nq)*h2)
    z=cdata(1)
    if (isign == 1) then
        cdata(1)=cmplx(real(z)+aimag(z),real(z)-aimag(z),kind=dpc)
    else
        cdata(1)=cmplx(c1*(real(z)+aimag(z)),c1*(real(z)-aimag(z)),kind=dpc)
        call four1(cdata,-1)
    end if
    if (present(zdata)) then
        if (isign /= 1) then
            data(1:n-1:2)=real(cdata)
            data(2:n:2)=aimag(cdata)
        end if
    else
        data(1:n-1:2)=real(cdata)
        data(2:n:2)=aimag(cdata)
        deallocate(cdata)
    end if
  END SUBROUTINE realft_dp

  SUBROUTINE four1_sp(data,isign)
    IMPLICIT NONE
    COMPLEX(SPC), DIMENSION(:), INTENT(INOUT) :: data
    INTEGER(I4B), INTENT(IN) :: isign
    COMPLEX(SPC), DIMENSION(:,:), ALLOCATABLE :: dat,temp
    COMPLEX(DPC), DIMENSION(:), ALLOCATABLE :: w,wp
    REAL(DP), DIMENSION(:), ALLOCATABLE :: theta
    INTEGER(I4B) :: n,m1,m2,j
    n=size(data)
    call assert(iand(n,n-1)==0, 'n must be a power of 2 in four1_sp!')
    m1=2**ceiling(0.5_sp*log(real(n,sp))/log(2.0_sp))
    m2=n/m1
    allocate(dat(m1,m2),theta(m1),w(m1),wp(m1),temp(m2,m1))
    dat=reshape(data,shape(dat))
    call fourrow(dat,isign)
    theta=arth(0,isign,m1)*TWOPI_D/n
    wp=cmplx(-2.0_dp*sin(0.5_dp*theta)**2,sin(theta),kind=dpc)
    w=cmplx(1.0_dp,0.0_dp,kind=dpc)
    do j=2,m2
        w=w*wp+w
        dat(:,j)=dat(:,j)*w
    end do
    temp=transpose(dat)
    call fourrow(temp,isign)
    data=reshape(temp,shape(data))
    deallocate(dat,w,wp,theta,temp)
  END SUBROUTINE four1_sp

  SUBROUTINE four1_dp(data,isign)
    IMPLICIT NONE
    COMPLEX(DPC), DIMENSION(:), INTENT(INOUT) :: data
    INTEGER(I4B), INTENT(IN) :: isign
    COMPLEX(DPC), DIMENSION(:,:), ALLOCATABLE :: dat,temp
    COMPLEX(DPC), DIMENSION(:), ALLOCATABLE :: w,wp
    REAL(DP), DIMENSION(:), ALLOCATABLE :: theta
    INTEGER(I4B) :: n,m1,m2,j
    n=size(data)
    call assert(iand(n,n-1)==0, 'n must be a power of 2 in four1_dp!')
    m1=2**ceiling(0.5_sp*log(real(n,sp))/log(2.0_sp))
    m2=n/m1
    allocate(dat(m1,m2),theta(m1),w(m1),wp(m1),temp(m2,m1))
    dat=reshape(data,shape(dat))
    call fourrow(dat,isign)
    theta=arth(0,isign,m1)*TWOPI_D/n
    wp=cmplx(-2.0_dp*sin(0.5_dp*theta)**2,sin(theta),kind=dpc)
    w=cmplx(1.0_dp,0.0_dp,kind=dpc)
    do j=2,m2
        w=w*wp+w
        dat(:,j)=dat(:,j)*w
    end do
    temp=transpose(dat)
    call fourrow(temp,isign)
    data=reshape(temp,shape(data))
    deallocate(dat,w,wp,theta,temp)
  END SUBROUTINE four1_dp

  SUBROUTINE fourrow_sp(data,isign)
    IMPLICIT NONE
    COMPLEX(SPC), DIMENSION(:,:), INTENT(INOUT) :: data
    INTEGER(I4B), INTENT(IN) :: isign
    INTEGER(I4B) :: n,i,istep,j,m,mmax,n2
    REAL(DP) :: theta
    COMPLEX(SPC), DIMENSION(size(data,1)) :: temp
    COMPLEX(DPC) :: w,wp
    COMPLEX(SPC) :: ws
    n=size(data,2)
    call assert(iand(n,n-1)==0, 'n must be a power of 2 in fourrow_sp!')
    n2=n/2
    j=n2
    do i=1,n-2
        if (j > i) call swap(data(:,j+1),data(:,i+1))
        m=n2
        do
            if (m < 2 .or. j < m) exit
            j=j-m
            m=m/2
        end do
        j=j+m
    end do
    mmax=1
    do
        if (n <= mmax) exit
        istep=2*mmax
        theta=PI_D/(isign*mmax)
        wp=cmplx(-2.0_dp*sin(0.5_dp*theta)**2,sin(theta),kind=dpc)
        w=cmplx(1.0_dp,0.0_dp,kind=dpc)
        do m=1,mmax
            ws=w
            do i=m,n,istep
                j=i+mmax
                temp=ws*data(:,j)
                data(:,j)=data(:,i)-temp
                data(:,i)=data(:,i)+temp
            end do
            w=w*wp+w
        end do
        mmax=istep
    end do
  END SUBROUTINE fourrow_sp

  SUBROUTINE fourrow_dp(data,isign)
    IMPLICIT NONE
    COMPLEX(DPC), DIMENSION(:,:), INTENT(INOUT) :: data
    INTEGER(I4B), INTENT(IN) :: isign
    INTEGER(I4B) :: n,i,istep,j,m,mmax,n2
    REAL(DP) :: theta
    COMPLEX(DPC), DIMENSION(size(data,1)) :: temp
    COMPLEX(DPC) :: w,wp
    COMPLEX(DPC) :: ws
    n=size(data,2)
    call assert(iand(n,n-1)==0, 'n must be a power of 2 in fourrow_dp!')
    n2=n/2
    j=n2
    do i=1,n-2
        if (j > i) call swap(data(:,j+1),data(:,i+1))
        m=n2
        do
            if (m < 2 .or. j < m) exit
            j=j-m
            m=m/2
        end do
        j=j+m
    end do
    mmax=1
    do
        if (n <= mmax) exit
        istep=2*mmax
        theta=PI_D/(isign*mmax)
        wp=cmplx(-2.0_dp*sin(0.5_dp*theta)**2,sin(theta),kind=dpc)
        w=cmplx(1.0_dp,0.0_dp,kind=dpc)
        do m=1,mmax
            ws=w
            do i=m,n,istep
                j=i+mmax
                temp=ws*data(:,j)
                data(:,j)=data(:,i)-temp
                data(:,i)=data(:,i)+temp
            end do
            w=w*wp+w
        end do
        mmax=istep
    end do
  END SUBROUTINE fourrow_dp

END MODULE long_v2





