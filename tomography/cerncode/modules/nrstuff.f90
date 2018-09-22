MODULE nrstuff
! From Numerical recepies
! From nrtype.f90
  INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
  INTEGER, PARAMETER :: I2B = SELECTED_INT_KIND(4)
  INTEGER, PARAMETER :: I1B = SELECTED_INT_KIND(2)
  INTEGER, PARAMETER :: SP = KIND(1.0)
  INTEGER, PARAMETER :: DP = KIND(1.0D0)
  INTEGER, PARAMETER :: SPC = KIND((1.0,1.0))
  INTEGER, PARAMETER :: DPC = KIND((1.0D0,1.0D0))
  INTEGER, PARAMETER :: LGT = KIND(.true.)
  REAL(SP), PARAMETER :: PI=3.141592653589793238462643383279502884197_sp
  REAL(SP), PARAMETER :: PIO2=1.57079632679489661923132169163975144209858_sp
  REAL(SP), PARAMETER :: TWOPI=6.283185307179586476925286766559005768394_sp
  REAL(SP), PARAMETER :: SQRT2=1.41421356237309504880168872420969807856967_sp
  REAL(SP), PARAMETER :: EULER=0.5772156649015328606065120900824024310422_sp
  REAL(DP), PARAMETER :: PI_D=3.141592653589793238462643383279502884197_dp
  REAL(DP), PARAMETER :: PIO2_D=1.57079632679489661923132169163975144209858_dp
  REAL(DP), PARAMETER :: TWOPI_D=6.283185307179586476925286766559005768394_dp

! From nrutil.f90
  INTEGER(I4B), PARAMETER :: NPAR_ARTH=16,NPAR2_ARTH=8
  INTEGER(I4B), PARAMETER :: NPAR_GEOP=4,NPAR2_GEOP=2
  INTEGER(I4B), PARAMETER :: NPAR_CUMSUM=16
  INTEGER(I4B), PARAMETER :: NPAR_CUMPROD=8
  INTEGER(I4B), PARAMETER :: NPAR_POLY=8
  INTEGER(I4B), PARAMETER :: NPAR_POLYTERM=8

  INTERFACE assert
     MODULE PROCEDURE assert1,assert2,assert3,assert4,assert_v
  END INTERFACE
  INTERFACE assert_eq
     MODULE PROCEDURE assert_eq2,assert_eq3,assert_eq4,assert_eqn
  END INTERFACE
  INTERFACE arth
     MODULE PROCEDURE arth_r, arth_d, arth_i
  END INTERFACE
  INTERFACE diagmult
     MODULE PROCEDURE diagmult_rv,diagmult_r
  END INTERFACE
  INTERFACE gammq
     MODULE PROCEDURE gammq_s,gammq_v
  END INTERFACE
  INTERFACE gammln
     MODULE PROCEDURE gammln_s,gammln_v
  END INTERFACE
  INTERFACE gser
     MODULE PROCEDURE gser_s,gser_v
  END INTERFACE
  INTERFACE gcf
     MODULE PROCEDURE gcf_s,gcf_v
  END INTERFACE
  INTERFACE imaxloc
     MODULE PROCEDURE imaxloc_r,imaxloc_i
  END INTERFACE
  INTERFACE swap
     MODULE PROCEDURE swap_i,swap_r,swap_d,swap_rv,swap_c, &
        swap_cv,swap_cm,swap_z,swap_zv,swap_zm, &
	masked_swap_rs,masked_swap_rv,masked_swap_rm
  END INTERFACE
  INTERFACE outerprod
     MODULE PROCEDURE outerprod_r,outerprod_d
  END INTERFACE
  INTERFACE poly
     MODULE PROCEDURE poly_rr,poly_rrv,poly_dd,poly_ddv,&
                      poly_rc,poly_cc,poly_msk_rrv,poly_msk_ddv
  END INTERFACE
  INTERFACE vabs
     MODULE PROCEDURE vabs_r,vabs_d
  END INTERFACE
CONTAINS
  !BL
  FUNCTION assert_eq2(n1,n2,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, INTENT(IN) :: n1,n2
    INTEGER :: assert_eq2
    if (n1 == n2) then
       assert_eq2=n1
    else
       write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       STOP 'program terminated by assert_eq2'
    end if
  END FUNCTION assert_eq2
  !BL
  FUNCTION assert_eq3(n1,n2,n3,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, INTENT(IN) :: n1,n2,n3
    INTEGER :: assert_eq3
    if (n1 == n2 .and. n2 == n3) then
       assert_eq3=n1
    else
       write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       STOP 'program terminated by assert_eq3'
    end if
  END FUNCTION assert_eq3
  !BL
  FUNCTION assert_eq4(n1,n2,n3,n4,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, INTENT(IN) :: n1,n2,n3,n4
    INTEGER :: assert_eq4
    if (n1 == n2 .and. n2 == n3 .and. n3 == n4) then
       assert_eq4=n1
    else
       write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       STOP 'program terminated by assert_eq4'
    end if
  END FUNCTION assert_eq4
  !BL
  FUNCTION assert_eqn(nn,string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    INTEGER, DIMENSION(:), INTENT(IN) :: nn
    INTEGER :: assert_eqn
    if (all(nn(2:) == nn(1))) then
       assert_eqn=nn(1)
    else
       write (*,*) 'nrerror: an assert_eq failed with this tag:', &
            string
       STOP 'program terminated by assert_eqn'
    end if
  END FUNCTION assert_eqn

!BL
  FUNCTION arth_r(first,increment,n)
        REAL(SP), INTENT(IN) :: first,increment
        INTEGER(I4B), INTENT(IN) :: n
        REAL(SP), DIMENSION(n) :: arth_r
        INTEGER(I4B) :: k,k2
        REAL(SP) :: temp
        if (n > 0) arth_r(1)=first
        if (n <= NPAR_ARTH) then
                do k=2,n
                        arth_r(k)=arth_r(k-1)+increment
                end do
        else
                do k=2,NPAR2_ARTH
                        arth_r(k)=arth_r(k-1)+increment
                end do
                temp=increment*NPAR2_ARTH
                k=NPAR2_ARTH
                do
                        if (k >= n) exit
                        k2=k+k
                        arth_r(k+1:min(k2,n))=temp+arth_r(1:min(k,n-k))
                        temp=temp+temp
                        k=k2
                end do
        end if
  END FUNCTION arth_r
!BL
  FUNCTION arth_d(first,increment,n)
        REAL(DP), INTENT(IN) :: first,increment
        INTEGER(I4B), INTENT(IN) :: n
        REAL(DP), DIMENSION(n) :: arth_d
        INTEGER(I4B) :: k,k2
        REAL(DP) :: temp
        if (n > 0) arth_d(1)=first
        if (n <= NPAR_ARTH) then
                do k=2,n
                        arth_d(k)=arth_d(k-1)+increment
                end do
        else
                do k=2,NPAR2_ARTH
                        arth_d(k)=arth_d(k-1)+increment
                end do
                temp=increment*NPAR2_ARTH
                k=NPAR2_ARTH
                do
                        if (k >= n) exit
                        k2=k+k
                        arth_d(k+1:min(k2,n))=temp+arth_d(1:min(k,n-k))
                        temp=temp+temp
                        k=k2
                end do
        end if
  END FUNCTION arth_d
!BL
  FUNCTION arth_i(first,increment,n)
        INTEGER(I4B), INTENT(IN) :: first,increment,n
        INTEGER(I4B), DIMENSION(n) :: arth_i
        INTEGER(I4B) :: k,k2,temp
        if (n > 0) arth_i(1)=first
        if (n <= NPAR_ARTH) then
                do k=2,n
                        arth_i(k)=arth_i(k-1)+increment
                end do
        else
                do k=2,NPAR2_ARTH
                        arth_i(k)=arth_i(k-1)+increment
                end do
                temp=increment*NPAR2_ARTH
                k=NPAR2_ARTH
                do
                        if (k >= n) exit
                        k2=k+k
                        arth_i(k+1:min(k2,n))=temp+arth_i(1:min(k,n-k))
                        temp=temp+temp
                        k=k2
                end do
        end if
  END FUNCTION arth_i
!BL
  SUBROUTINE assert1(n1,string)
        CHARACTER(LEN=*), INTENT(IN) :: string
        LOGICAL, INTENT(IN) :: n1
        if (.not. n1) then
                write (*,*) 'nrerror: an assertion failed with this tag:', &
                        string
                STOP 'program terminated by assert1'
        end if
  END SUBROUTINE assert1
!BL
  SUBROUTINE assert2(n1,n2,string)
        CHARACTER(LEN=*), INTENT(IN) :: string
        LOGICAL, INTENT(IN) :: n1,n2
        if (.not. (n1 .and. n2)) then
                write (*,*) 'nrerror: an assertion failed with this tag:', &
                        string
                STOP 'program terminated by assert2'
        end if
  END SUBROUTINE assert2
!BL
  SUBROUTINE assert3(n1,n2,n3,string)
        CHARACTER(LEN=*), INTENT(IN) :: string
        LOGICAL, INTENT(IN) :: n1,n2,n3
        if (.not. (n1 .and. n2 .and. n3)) then
                write (*,*) 'nrerror: an assertion failed with this tag:', &
                        string
                STOP 'program terminated by assert3'
        end if
  END SUBROUTINE assert3
!BL
  SUBROUTINE assert4(n1,n2,n3,n4,string)
        CHARACTER(LEN=*), INTENT(IN) :: string
        LOGICAL, INTENT(IN) :: n1,n2,n3,n4
        if (.not. (n1 .and. n2 .and. n3 .and. n4)) then
                write (*,*) 'nrerror: an assertion failed with this tag:', &
                        string
                STOP 'program terminated by assert4'
        end if
  END SUBROUTINE assert4
!BL
  SUBROUTINE assert_v(n,string)
        CHARACTER(LEN=*), INTENT(IN) :: string
        LOGICAL, DIMENSION(:), INTENT(IN) :: n
        if (.not. all(n)) then
                write (*,*) 'nrerror: an assertion failed with this tag:', &
                        string
                STOP 'program terminated by assert_v'
        end if
  END SUBROUTINE assert_v
!BL
  SUBROUTINE diagmult_rv(mat,diag)
    REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: mat
    REAL(SP), DIMENSION(:), INTENT(IN) :: diag
    INTEGER(I4B) :: j,n
    n = assert_eq2(size(diag),min(size(mat,1),size(mat,2)),'diagmult_rv')
    do j=1,n
       mat(j,j)=mat(j,j)*diag(j)
    end do
  END SUBROUTINE diagmult_rv
  !BL
  SUBROUTINE diagmult_r(mat,diag)
    REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: mat
    REAL(SP), INTENT(IN) :: diag
    INTEGER(I4B) :: j,n
    n = min(size(mat,1),size(mat,2))
    do j=1,n
       mat(j,j)=mat(j,j)*diag
    end do
  END SUBROUTINE diagmult_r
!BL
  FUNCTION gammln_s(xx)
        IMPLICIT NONE
        REAL(SP), INTENT(IN) :: xx
        REAL(SP) :: gammln_s
        REAL(DP) :: tmp,x
        REAL(DP) :: stp = 2.5066282746310005_dp
        REAL(DP), DIMENSION(6) :: coef = (/76.18009172947146_dp,&
                -86.50532032941677_dp,24.01409824083091_dp,&
                -1.231739572450155_dp,0.1208650973866179e-2_dp,&
                -0.5395239384953e-5_dp/)
        call assert(xx > 0.0, 'gammln_s arg')
        x=xx
        tmp=x+5.5_dp
        tmp=(x+0.5_dp)*log(tmp)-tmp
        gammln_s=tmp+log(stp*(1.000000000190015_dp+&
                sum(coef(:)/arth(x+1.0_dp,1.0_dp,size(coef))))/x)
   END FUNCTION gammln_s
!BL
   FUNCTION gammln_v(xx)
        IMPLICIT NONE
        INTEGER(I4B) :: i
        REAL(SP), DIMENSION(:), INTENT(IN) :: xx
        REAL(SP), DIMENSION(size(xx)) :: gammln_v
        REAL(DP), DIMENSION(size(xx)) :: ser,tmp,x,y
        REAL(DP) :: stp = 2.5066282746310005_dp
        REAL(DP), DIMENSION(6) :: coef = (/76.18009172947146_dp,&
                -86.50532032941677_dp,24.01409824083091_dp,&
                -1.231739572450155_dp,0.1208650973866179e-2_dp,&
                -0.5395239384953e-5_dp/)
        if (size(xx) == 0) RETURN
        call assert(all(xx > 0.0), 'gammln_v arg')
        x=xx
        tmp=x+5.5_dp
        tmp=(x+0.5_dp)*log(tmp)-tmp
        ser=1.000000000190015_dp
        y=x
        do i=1,size(coef)
                y=y+1.0_dp
                ser=ser+coef(i)/y
        end do
        gammln_v=tmp+log(stp*ser/x)
  END FUNCTION gammln_v
!BL
  FUNCTION gammq_s(a,x)
        IMPLICIT NONE
        REAL(SP), INTENT(IN) :: a,x
        REAL(SP) :: gammq_s
        call assert( x >= 0.0,  a > 0.0, 'gammq_s args')
        if (x<a+1.0_sp) then
                gammq_s=1.0_sp-gser(a,x)
        else
                gammq_s=gcf(a,x)
        end if
   END FUNCTION gammq_s
!BL
   FUNCTION gammq_v(a,x)
        IMPLICIT NONE
        REAL(SP), DIMENSION(:), INTENT(IN) :: a,x
        REAL(SP), DIMENSION(size(a)) :: gammq_v
        LOGICAL(LGT), DIMENSION(size(x)) :: mask
        INTEGER(I4B) :: ndum
        ndum=assert_eq(size(a),size(x),'gammq_v')
        call assert( all(x >= 0.0),  all(a > 0.0), 'gammq_v args')
        mask = (x<a+1.0_sp)
        gammq_v=merge(1.0_sp-gser(a,merge(x,0.0_sp,mask)), &
                gcf(a,merge(x,0.0_sp,.not. mask)),mask)
  END FUNCTION gammq_v
!BL
  FUNCTION gcf_s(a,x,gln)
        IMPLICIT NONE
        REAL(SP), INTENT(IN) :: a,x
        REAL(SP), OPTIONAL, INTENT(OUT) :: gln
        REAL(SP) :: gcf_s
        INTEGER(I4B), PARAMETER :: ITMAX=100
        REAL(SP), PARAMETER :: EPS=epsilon(x),FPMIN=tiny(x)/EPS
        INTEGER(I4B) :: i
        REAL(SP) :: an,b,c,d,del,h
        if (x == 0.0) then
                gcf_s=1.0
                RETURN
        end if
        b=x+1.0_sp-a
        c=1.0_sp/FPMIN
        d=1.0_sp/b
        h=d
        do i=1,ITMAX
                an=-i*(i-a)
                b=b+2.0_sp
                d=an*d+b
                if (abs(d) < FPMIN) d=FPMIN
                c=b+an/c
                if (abs(c) < FPMIN) c=FPMIN
                d=1.0_sp/d
                del=d*c
                h=h*del
                if (abs(del-1.0_sp) <= EPS) exit
        end do
        if (i > ITMAX) call nrerror('a too large, ITMAX too small in gcf_s')
        if (present(gln)) then
                gln=gammln(a)
                gcf_s=exp(-x+a*log(x)-gln)*h
        else
                gcf_s=exp(-x+a*log(x)-gammln(a))*h
        end if
  END FUNCTION gcf_s
!BL
  FUNCTION gcf_v(a,x,gln)
        IMPLICIT NONE
        REAL(SP), DIMENSION(:), INTENT(IN) :: a,x
        REAL(SP), DIMENSION(:), OPTIONAL, INTENT(OUT) :: gln
        REAL(SP), DIMENSION(size(a)) :: gcf_v
        INTEGER(I4B), PARAMETER :: ITMAX=100
        REAL(SP), PARAMETER :: EPS=epsilon(x),FPMIN=tiny(x)/EPS
        INTEGER(I4B) :: i
        REAL(SP), DIMENSION(size(a)) :: an,b,c,d,del,h
        LOGICAL(LGT), DIMENSION(size(a)) :: converged,zero
        i=assert_eq(size(a),size(x),'gcf_v')
        zero=(x == 0.0)
        where (zero)
                gcf_v=1.0
        elsewhere
                b=x+1.0_sp-a
                c=1.0_sp/FPMIN
                d=1.0_sp/b
                h=d
        end where
        converged=zero
        do i=1,ITMAX
                where (.not. converged)
                        an=-i*(i-a)
                        b=b+2.0_sp
                        d=an*d+b
                        d=merge(FPMIN,d, abs(d)<FPMIN )
                        c=b+an/c
                        c=merge(FPMIN,c, abs(c)<FPMIN )
                        d=1.0_sp/d
                        del=d*c
                        h=h*del
                        converged = (abs(del-1.0_sp)<=EPS)
                end where
                if (all(converged)) exit
        end do
        if (i > ITMAX) call nrerror('a too large, ITMAX too small in gcf_v')
        if (present(gln)) then
                if (size(gln) < size(a)) call &
                        nrerror('gser: Not enough space for gln')
                gln=gammln(a)
                where (.not. zero) gcf_v=exp(-x+a*log(x)-gln)*h
        else
                where (.not. zero) gcf_v=exp(-x+a*log(x)-gammln(a))*h
        end if
  END FUNCTION gcf_v
!BL
  FUNCTION gser_s(a,x,gln)
        IMPLICIT NONE
        REAL(SP), INTENT(IN) :: a,x
        REAL(SP), OPTIONAL, INTENT(OUT) :: gln
        REAL(SP) :: gser_s
        INTEGER(I4B), PARAMETER :: ITMAX=100
        REAL(SP), PARAMETER :: EPS=epsilon(x)
        INTEGER(I4B) :: n
        REAL(SP) :: ap,del,summ
        if (x == 0.0) then
                gser_s=0.0
                RETURN
        end if
        ap=a
        summ=1.0_sp/a
        del=summ
        do n=1,ITMAX
                ap=ap+1.0_sp
                del=del*x/ap
                summ=summ+del
                if (abs(del) < abs(summ)*EPS) exit
        end do
        if (n > ITMAX) call nrerror('a too large, ITMAX too small in gser_s')
        if (present(gln)) then
                gln=gammln(a)
                gser_s=summ*exp(-x+a*log(x)-gln)
        else
                gser_s=summ*exp(-x+a*log(x)-gammln(a))
        end if
  END FUNCTION gser_s
!BL
  FUNCTION gser_v(a,x,gln)
        IMPLICIT NONE
        REAL(SP), DIMENSION(:), INTENT(IN) :: a,x
        REAL(SP), DIMENSION(:), OPTIONAL, INTENT(OUT) :: gln
        REAL(SP), DIMENSION(size(a)) :: gser_v
        INTEGER(I4B), PARAMETER :: ITMAX=100
        REAL(SP), PARAMETER :: EPS=epsilon(x)
        INTEGER(I4B) :: n
        REAL(SP), DIMENSION(size(a)) :: ap,del,summ
        LOGICAL(LGT), DIMENSION(size(a)) :: converged,zero
        n=assert_eq(size(a),size(x),'gser_v')
        zero=(x == 0.0)
        where (zero) gser_v=0.0
        ap=a
        summ=1.0_sp/a
        del=summ
        converged=zero
        do n=1,ITMAX
                where (.not. converged)
                        ap=ap+1.0_sp
                        del=del*x/ap
                        summ=summ+del
                        converged = (abs(del) < abs(summ)*EPS)
                end where
                if (all(converged)) exit
        end do
        if (n > ITMAX) call nrerror('a too large, ITMAX too small in gser_v')
        if (present(gln)) then
                if (size(gln) < size(a)) call &
                        nrerror('gser: Not enough space for gln')
                gln=gammln(a)
                where (.not. zero) gser_v=summ*exp(-x+a*log(x)-gln)
        else
                where (.not. zero) gser_v=summ*exp(-x+a*log(x)-gammln(a))
        end if
  END FUNCTION gser_v
!BL
  SUBROUTINE swap_i(a,b)
    INTEGER(I4B), INTENT(INOUT) :: a,b
    INTEGER(I4B) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_i
!BL
  SUBROUTINE swap_r(a,b)
    REAL(SP), INTENT(INOUT) :: a,b
    REAL(SP) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_r
!BL
  SUBROUTINE swap_d(a,b)
    REAL(DP), DIMENSION(:), INTENT(INOUT) :: a,b
    REAL(DP), DIMENSION(SIZE(a)) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_d
!BL
  SUBROUTINE swap_rv(a,b)
    REAL(SP), DIMENSION(:), INTENT(INOUT) :: a,b
    REAL(SP), DIMENSION(SIZE(a)) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_rv
!BL
  SUBROUTINE swap_c(a,b)
    COMPLEX(SPC), INTENT(INOUT) :: a,b
    COMPLEX(SPC) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_c
!BL
  SUBROUTINE swap_cv(a,b)
    COMPLEX(SPC), DIMENSION(:), INTENT(INOUT) :: a,b
    COMPLEX(SPC), DIMENSION(SIZE(a)) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_cv
!BL
  SUBROUTINE swap_cm(a,b)
    COMPLEX(SPC), DIMENSION(:,:), INTENT(INOUT) :: a,b
    COMPLEX(SPC), DIMENSION(size(a,1),size(a,2)) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_cm
!BL
  SUBROUTINE swap_z(a,b)
    COMPLEX(DPC), INTENT(INOUT) :: a,b
    COMPLEX(DPC) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_z
!BL
  SUBROUTINE swap_zv(a,b)
    COMPLEX(DPC), DIMENSION(:), INTENT(INOUT) :: a,b
    COMPLEX(DPC), DIMENSION(SIZE(a)) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_zv
!BL
  SUBROUTINE swap_zm(a,b)
    COMPLEX(DPC), DIMENSION(:,:), INTENT(INOUT) :: a,b
    COMPLEX(DPC), DIMENSION(size(a,1),size(a,2)) :: dum
    dum=a
    a=b
    b=dum
  END SUBROUTINE swap_zm
!BL
  SUBROUTINE masked_swap_rs(a,b,mask)
    REAL(SP), INTENT(INOUT) :: a,b
    LOGICAL(LGT), INTENT(IN) :: mask
    REAL(SP) :: swp
    if (mask) then
       swp=a
       a=b
       b=swp
    end if
  END SUBROUTINE masked_swap_rs
!BL
  SUBROUTINE masked_swap_rv(a,b,mask)
    REAL(SP), DIMENSION(:), INTENT(INOUT) :: a,b
    LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: mask
    REAL(SP), DIMENSION(size(a)) :: swp
    where (mask)
       swp=a
       a=b
       b=swp
    end where
  END SUBROUTINE masked_swap_rv
!BL
  SUBROUTINE masked_swap_rm(a,b,mask)
    REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a,b
    LOGICAL(LGT), DIMENSION(:,:), INTENT(IN) :: mask
    REAL(SP), DIMENSION(size(a,1),size(a,2)) :: swp
    where (mask)
       swp=a
       a=b
       b=swp
    end where
  END SUBROUTINE masked_swap_rm
!BL
  FUNCTION outerprod_r(a,b)
    REAL(SP), DIMENSION(:), INTENT(IN) :: a,b
    REAL(SP), DIMENSION(size(a),size(b)) :: outerprod_r
    outerprod_r = spread(a,dim=2,ncopies=size(b)) * &
         spread(b,dim=1,ncopies=size(a))
  END FUNCTION outerprod_r
!BL
  FUNCTION outerprod_d(a,b)
    REAL(DP), DIMENSION(:), INTENT(IN) :: a,b
    REAL(DP), DIMENSION(size(a),size(b)) :: outerprod_d
    outerprod_d = spread(a,dim=2,ncopies=size(b)) * &
         spread(b,dim=1,ncopies=size(a))
  END FUNCTION outerprod_d
!BL
  FUNCTION vabs_r(v)
    REAL(SP), DIMENSION(:), INTENT(IN) :: v
    REAL(SP) :: vabs_r
    vabs_r=sqrt(dot_product(v,v))
  END FUNCTION vabs_r
!BL
  FUNCTION vabs_d(v)
    REAL(DP), DIMENSION(:), INTENT(IN) :: v
    REAL(DP) :: vabs_d
    vabs_d=sqrt(dot_product(v,v))
  END FUNCTION vabs_d
!BL
  FUNCTION outerand(a,b)
    LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: a,b
    LOGICAL(LGT), DIMENSION(size(a),size(b)) :: outerand
    outerand = spread(a,dim=2,ncopies=size(b)) .and. &
         spread(b,dim=1,ncopies=size(a))
  END FUNCTION outerand
!BL
  SUBROUTINE nrerror(string)
    CHARACTER(LEN=*), INTENT(IN) :: string
    write (*,*) 'nrerror: ',string
    STOP 'program terminated by nrerror'
  END SUBROUTINE nrerror
!BL
  FUNCTION zroots_unity(n,nn)
    INTEGER(I4B), INTENT(IN) :: n,nn
    COMPLEX(SPC), DIMENSION(nn) :: zroots_unity
    INTEGER(I4B) :: k
    REAL(SP) :: theta
    zroots_unity(1)=1.0
    theta=TWOPI/n
    k=1
    do
      if (k >= nn) exit
      zroots_unity(k+1)=cmplx(cos(k*theta),sin(k*theta),SPC)
      zroots_unity(k+2:min(2*k,nn))=zroots_unity(k+1)*&
                   zroots_unity(2:min(k,nn-k))
      k=2*k
    end do
  END FUNCTION zroots_unity
!BL
  FUNCTION imaxloc_r(arr)
    REAL(SP), DIMENSION(:), INTENT(IN) :: arr
    INTEGER(I4B) :: imaxloc_r
    INTEGER(I4B), DIMENSION(1) :: imax
    imax=maxloc(arr(:))
    imaxloc_r=imax(1)
  END FUNCTION imaxloc_r
!BL
  FUNCTION imaxloc_i(iarr)
    INTEGER(I4B), DIMENSION(:), INTENT(IN) :: iarr
    INTEGER(I4B), DIMENSION(1) :: imax
    INTEGER(I4B) :: imaxloc_i
    imax=maxloc(iarr(:))
    imaxloc_i=imax(1)
  END FUNCTION imaxloc_i
!BL
  FUNCTION poly_rr(x,coeffs)
        REAL(SP), INTENT(IN) :: x
        REAL(SP), DIMENSION(:), INTENT(IN) :: coeffs
        REAL(SP) :: poly_rr
        REAL(SP) :: pow
        REAL(SP), DIMENSION(:), ALLOCATABLE :: vec
        INTEGER(I4B) :: i,n,nn
        n=size(coeffs)
        if (n <= 0) then
                poly_rr=0.0_sp
        else if (n < NPAR_POLY) then
                poly_rr=coeffs(n)
                do i=n-1,1,-1
                        poly_rr=x*poly_rr+coeffs(i)
                end do
        else
                allocate(vec(n+1))
                pow=x
                vec(1:n)=coeffs
                do
                        vec(n+1)=0.0_sp
                        nn=ishft(n+1,-1)
                        vec(1:nn)=vec(1:n:2)+pow*vec(2:n+1:2)
                        if (nn == 1) exit
                        pow=pow*pow
                        n=nn
                end do
                poly_rr=vec(1)
                deallocate(vec)
        end if
  END FUNCTION poly_rr
!BL
  FUNCTION poly_dd(x,coeffs)
        REAL(DP), INTENT(IN) :: x
        REAL(DP), DIMENSION(:), INTENT(IN) :: coeffs
        REAL(DP) :: poly_dd
        REAL(DP) :: pow
        REAL(DP), DIMENSION(:), ALLOCATABLE :: vec
        INTEGER(I4B) :: i,n,nn
        n=size(coeffs)
        if (n <= 0) then
                poly_dd=0.0_dp
        else if (n < NPAR_POLY) then
                poly_dd=coeffs(n)
                do i=n-1,1,-1
                        poly_dd=x*poly_dd+coeffs(i)
                end do
        else
                allocate(vec(n+1))
                pow=x
                vec(1:n)=coeffs
                do
                        vec(n+1)=0.0_dp
                        nn=ishft(n+1,-1)
                        vec(1:nn)=vec(1:n:2)+pow*vec(2:n+1:2)
                        if (nn == 1) exit
                        pow=pow*pow
                        n=nn
                end do
                poly_dd=vec(1)
                deallocate(vec)
        end if
  END FUNCTION poly_dd
!BL
  FUNCTION poly_rc(x,coeffs)
        COMPLEX(SPC), INTENT(IN) :: x
        REAL(SP), DIMENSION(:), INTENT(IN) :: coeffs
        COMPLEX(SPC) :: poly_rc
        COMPLEX(SPC) :: pow
        COMPLEX(SPC), DIMENSION(:), ALLOCATABLE :: vec
        INTEGER(I4B) :: i,n,nn
        n=size(coeffs)
        if (n <= 0) then
                poly_rc=0.0_sp
        else if (n < NPAR_POLY) then
                poly_rc=coeffs(n)
                do i=n-1,1,-1
                        poly_rc=x*poly_rc+coeffs(i)
                end do
        else
                allocate(vec(n+1))
                pow=x
                vec(1:n)=coeffs
                do
                        vec(n+1)=0.0_sp
                        nn=ishft(n+1,-1)
                        vec(1:nn)=vec(1:n:2)+pow*vec(2:n+1:2)
                        if (nn == 1) exit
                        pow=pow*pow
                        n=nn
                end do
                poly_rc=vec(1)
                deallocate(vec)
        end if
  END FUNCTION poly_rc
!BL
  FUNCTION poly_cc(x,coeffs)
        COMPLEX(SPC), INTENT(IN) :: x
        COMPLEX(SPC), DIMENSION(:), INTENT(IN) :: coeffs
        COMPLEX(SPC) :: poly_cc
        COMPLEX(SPC) :: pow
        COMPLEX(SPC), DIMENSION(:), ALLOCATABLE :: vec
        INTEGER(I4B) :: i,n,nn
        n=size(coeffs)
        if (n <= 0) then
                poly_cc=0.0_sp
        else if (n < NPAR_POLY) then
                poly_cc=coeffs(n)
                do i=n-1,1,-1
                        poly_cc=x*poly_cc+coeffs(i)
                end do
        else
                allocate(vec(n+1))
                pow=x
                vec(1:n)=coeffs
                do
                        vec(n+1)=0.0_sp
                        nn=ishft(n+1,-1)
                        vec(1:nn)=vec(1:n:2)+pow*vec(2:n+1:2)
                        if (nn == 1) exit
                        pow=pow*pow
                        n=nn
                end do
                poly_cc=vec(1)
                deallocate(vec)
        end if
  END FUNCTION poly_cc
!BL
  FUNCTION poly_rrv(x,coeffs)
        REAL(SP), DIMENSION(:), INTENT(IN) :: coeffs,x
        REAL(SP), DIMENSION(size(x)) :: poly_rrv
        INTEGER(I4B) :: i,n,m
        m=size(coeffs)
        n=size(x)
        if (m <= 0) then
                poly_rrv=0.0_sp
        else if (m < n .or. m < NPAR_POLY) then
                poly_rrv=coeffs(m)
                do i=m-1,1,-1
                        poly_rrv=x*poly_rrv+coeffs(i)
                end do
        else
                do i=1,n
                        poly_rrv(i)=poly_rr(x(i),coeffs)
                end do
        end if
  END FUNCTION poly_rrv
!BL
  FUNCTION poly_ddv(x,coeffs)
        REAL(DP), DIMENSION(:), INTENT(IN) :: coeffs,x
        REAL(DP), DIMENSION(size(x)) :: poly_ddv
        INTEGER(I4B) :: i,n,m
        m=size(coeffs)
        n=size(x)
        if (m <= 0) then
                poly_ddv=0.0_dp
        else if (m < n .or. m < NPAR_POLY) then
                poly_ddv=coeffs(m)
                do i=m-1,1,-1
                        poly_ddv=x*poly_ddv+coeffs(i)
                end do
        else
                do i=1,n
                        poly_ddv(i)=poly_dd(x(i),coeffs)
                end do
        end if
  END FUNCTION poly_ddv
!BL
  FUNCTION poly_msk_rrv(x,coeffs,mask)
        REAL(SP), DIMENSION(:), INTENT(IN) :: coeffs,x
        LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: mask
        REAL(SP), DIMENSION(size(x)) :: poly_msk_rrv
        poly_msk_rrv=unpack(poly_rrv(pack(x,mask),coeffs),mask,0.0_sp)
  END FUNCTION poly_msk_rrv
!BL
  FUNCTION poly_msk_ddv(x,coeffs,mask)
        REAL(DP), DIMENSION(:), INTENT(IN) :: coeffs,x
        LOGICAL(LGT), DIMENSION(:), INTENT(IN) :: mask
        REAL(DP), DIMENSION(size(x)) :: poly_msk_ddv
        poly_msk_ddv=unpack(poly_ddv(pack(x,mask),coeffs),mask,0.0_dp)
  END FUNCTION poly_msk_ddv
!BL
END MODULE nrstuff
