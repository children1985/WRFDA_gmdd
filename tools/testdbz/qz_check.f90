program linear_check

implicit none

integer,parameter             :: nx=449,ny=449,nz=42
real,dimension(nx,ny,nz)      :: p,tc,qv,qi,qc,qr,qs,qg,nr,ns,ng
real,dimension(nx,ny,nz)      :: dqr,dqs,dqg,dnr,dns,dng,dtc,dqv
real,dimension(nx,ny,nz)      :: dbz_org,dbz_prt,zmm_org,zmm_prt
real,dimension(nx,ny,nz)      :: dbz_lin,zmm_lin

real,dimension(nx,ny,nz)      :: tc_org,qv_org,qr_org,qs_org,qg_org
real,dimension(nx,ny,nz)      :: tc_prt,qv_prt,qr_prt,qs_prt,qg_prt

real,dimension(nx,ny,nz)      :: dbz_grd,zmm_grd

real                          :: rn0_r,rn0_s,rn0_g
real                          :: rhos,rhog

integer                       :: tlopt0=2,tlopt
integer                       :: gropt0=2,gropt
integer                       :: keepconstopt=234  ! 1 rho is constant during minimization
                                                 ! 2 F (fraction) is constant during minimization
                                                 ! 3 SD of the canting angle for hail is constant during minimization
												 ! 4 rhom is constatnt
						                         ! 12 denotes the combination of 1 & 2, so does the 13,and 23,234
integer                       :: i,j,k,ii,jj
real,dimension(nx,ny,nz)      :: temra01,temra02,temra03

integer                       :: ivar,seed
real                          :: stdt,stdqv,stdqr,stdqs,stdqg
real                          :: ratio ! sum[h(x+dx)-H(dx)]/sum[h(x)]
real                          :: ratio2 ! sum(h(x+dx)-h(x))/sum(H(dx))
real                          :: mag_std
real                          :: min_mag=1e-5,mag_int
integer                       :: is,ie,js,je,ks,ke
integer                       :: print_lvl=5
real                          :: mag_thres=0.00
integer                       :: print_init
integer                       :: radii=100
integer                       :: seed_init=200

real                          :: pmean,tmean,qvmean

integer                       :: count_neg_right,count_neg_wrong
integer                       :: count_pos_right,count_pos_wrong
real                          :: temr01,temr02,temr03,temr04,temr05

real                          :: grd_chk_lt,grd_chk_rt  !! left <y,y> right <x,(A*)y> 

character(len=256)            :: diagfilename
character(len=128)            :: wrtfmt
integer                       :: icount
character(len=1)              :: temc1_01

integer,parameter             :: nqbin=301
real                          :: binint=0.02e-3  ! kg/kg
real,dimension(nqbin)         :: qrtest,qstest,qgtest,qrz,qsz,qgz
real,dimension(nqbin,nqbin)   :: qrqsz,qrqgz,qtem
real                          :: ibin

!----------------------------------------------------
!    begin executable codes
!---------------------------------------------------

rn0_r=8e6
rn0_s=3e6
rn0_g=4e5
rhos=100.0
rhog=400.0
print_init=1

qrtest=0
qstest=0
qgtest=0
qrqsz=0
qrqgz=0
qrz=0
qsz=0
qgz=0

do i=1,nqbin
  qrtest(i)=(i-1)*binint
enddo
qstest=qrtest
qgtest=qrtest

diagfilename="dbzdata_qz_check"
open(1234,file=trim(diagfilename)//'.dat',form='binary')

icount=0
print*,"test begin"


call readdata()

p=p*100.0
tc=tc+273.15
qr=qr
qs=qs
qg=qg

nr=0
ns=0
ng=0


is=250  !nx/2-radii
ie=300  !nx/2+radii
js=200  !ny/2-radii
je=250  !ny/2+radii
ks=15
ke=15

if(print_init==1) print('(A,E12.5)'),"qr mean",sum(qr(is:ie,js:je,ks:ke))/((ie-is+1)*(je-js+1)*(ke-ks+1))
if(print_init==1) print('(A,E12.5)'),"qs mean",sum(qs(is:ie,js:je,ks:ke))/((ie-is+1)*(je-js+1)*(ke-ks+1))
if(print_init==1) print('(A,E12.5)'),"qg mean",sum(qg(is:ie,js:je,ks:ke))/((ie-is+1)*(je-js+1)*(ke-ks+1))

print_init=0

if(print_lvl>=10) print*,"calc original dbz"
tc_org=tc
qv_org=qv
qr_org=qr
qs_org=qs
qg_org=qg

tlopt=0
gropt=0
do k=ks,ke

  pmean=sum(p(is:ie,js:je,ks:ke))/((ie-is+1)*(je-js+1)*(ke-ks+1))
  tmean=sum(tc(is:ie,js:je,ks:ke))/((ie-is+1)*(je-js+1)*(ke-ks+1))
  qvmean=sum(qv(is:ie,js:je,ks:ke))/((ie-is+1)*(je-js+1)*(ke-ks+1))
  print*,"mean p,t,qv",pmean,tmean,qvmean
  print*,"max qrsg",maxval(qrtest,1)
  i=is
  j=js

  ! test qr z
  do ii=1,nqbin
!      call oudbzcalc_lin(qv(i,j,k),qr(i,j,k),qs(i,j,k),qg(i,j,k),nr(i,j,k),           &
!                         ns(i,j,k),ng(i,j,k),tc(i,j,k),p(i,j,k),dbz_org(i,j,k),       &
!                         1,1,1,0,0,0,rn0_r,rn0_s,rn0_g,                               &
!                         rhos,rhog,dtc(i,j,k),dqv(i,j,k),dqr(i,j,k),dqs(i,j,k),dqg(i,j,k),                  &
!                         dnr(i,j,k),dns(i,j,k),dng(i,j,k),zmm_org(i,j,k),             &
!                         tlopt,gropt,i,j,k,zmm_org(i,j,k),keepconstopt)

      call oudbzcalc_lin(qvmean,qrtest(ii),qstest(ii)*0.0,qgtest(ii)*0.0,nr(i,j,k),                  &
                         ns(i,j,k),ng(i,j,k),tmean,pmean,dbz_org(i,j,k),                          &
                         1,1,1,0,0,0,rn0_r,rn0_s,rn0_g,                                           &
                         rhos,rhog,dtc(i,j,k),dqv(i,j,k),dqr(i,j,k),dqs(i,j,k),dqg(i,j,k),        &
                         dnr(i,j,k),dns(i,j,k),dng(i,j,k),zmm_org(i,j,k),                         &
                         tlopt,gropt,i,j,k,zmm_org(i,j,k),keepconstopt,0)						 
      qrz(ii)=dbz_org(i,j,k)    
  enddo
  
  ! test qs z
  do ii=1,nqbin

      call oudbzcalc_lin(qvmean,qrtest(ii)*0.0,qstest(ii),qgtest(ii)*0.0,nr(i,j,k),                  &
                         ns(i,j,k),ng(i,j,k),tmean,pmean,dbz_org(i,j,k),                          &
                         1,1,1,0,0,0,rn0_r,rn0_s,rn0_g,                                           &
                         rhos,rhog,dtc(i,j,k),dqv(i,j,k),dqr(i,j,k),dqs(i,j,k),dqg(i,j,k),        &
                         dnr(i,j,k),dns(i,j,k),dng(i,j,k),zmm_org(i,j,k),                         &
                         tlopt,gropt,i,j,k,zmm_org(i,j,k),keepconstopt,0)						 
      qsz(ii)=dbz_org(i,j,k)
  enddo  
  
  ! test qg z
  do ii=1,nqbin

      call oudbzcalc_lin(qvmean,qrtest(ii)*0.0,qstest(ii)*0.0,qgtest(ii),nr(i,j,k),                  &
                         ns(i,j,k),ng(i,j,k),tmean,pmean,dbz_org(i,j,k),                          &
                         1,1,1,0,0,0,rn0_r,rn0_s,rn0_g,                                           &
                         rhos,rhog,dtc(i,j,k),dqv(i,j,k),dqr(i,j,k),dqs(i,j,k),dqg(i,j,k),        &
                         dnr(i,j,k),dns(i,j,k),dng(i,j,k),zmm_org(i,j,k),                         &
                         tlopt,gropt,i,j,k,zmm_org(i,j,k),keepconstopt,0)						 
      qgz(ii)=dbz_org(i,j,k)
  enddo  
  

  ! test qr qs z
  do ii=1,nqbin
    do jj=1,nqbin
      call oudbzcalc_lin(qvmean,qrtest(ii),qstest(jj),qgtest(ii)*0.0,nr(i,j,k),                      &
                         ns(i,j,k),ng(i,j,k),tmean,pmean,dbz_org(i,j,k),                          &
                         1,1,1,0,0,0,rn0_r,rn0_s,rn0_g,                                           &
                         rhos,rhog,dtc(i,j,k),dqv(i,j,k),dqr(i,j,k),dqs(i,j,k),dqg(i,j,k),        &
                         dnr(i,j,k),dns(i,j,k),dng(i,j,k),zmm_org(i,j,k),                         &
                         tlopt,gropt,i,j,k,zmm_org(i,j,k),keepconstopt,1)						 
      qrqsz(ii,jj)=dbz_org(i,j,k)
	enddo
  enddo  

  ! test qr qg z
  do ii=1,nqbin
    do jj=1,nqbin
      call oudbzcalc_lin(qvmean,qrtest(ii),qstest(ii)*0.0,qgtest(jj),nr(i,j,k),                      &
                         ns(i,j,k),ng(i,j,k),tmean,pmean,dbz_org(i,j,k),                          &
                         1,1,1,0,0,0,rn0_r,rn0_s,rn0_g,                                           &
                         rhos,rhog,dtc(i,j,k),dqv(i,j,k),dqr(i,j,k),dqs(i,j,k),dqg(i,j,k),        &
                         dnr(i,j,k),dns(i,j,k),dng(i,j,k),zmm_org(i,j,k),                         &
                         tlopt,gropt,i,j,k,zmm_org(i,j,k),keepconstopt,2)						 
      qrqgz(ii,jj)=dbz_org(i,j,k)
	enddo
  enddo  

 
  
enddo
if(print_lvl>=5) print*,"done calc original dbz"
if(print_lvl>=5) print*,"qr z max min",maxval(qrz,1),minval(qrz,1)
if(print_lvl>=5) print*,"qs z max min",maxval(qsz,1),minval(qsz,1)
if(print_lvl>=5) print*,"qg z max min",maxval(qgz,1),minval(qgz,1)
if(print_lvl>=5) print*,"qr qs z max min",maxval(maxval(qrqsz,1),1),minval(minval(qrqsz,1),1) 
if(print_lvl>=5) print*,"qr qg z max min",maxval(maxval(qrqgz,1),1),minval(minval(qrqgz,1),1)  
if(print_lvl>=5) print*,"------------------------------------------------"

do i=1,nqbin
  qtem(:,i)=qrz(:)
enddo
write(1234) qtem
do i=1,nqbin
  qtem(:,i)=qsz(:)
enddo
write(1234) qtem
do i=1,nqbin
  qtem(:,i)=qgz(:)
enddo
write(1234) qtem
write(1234) qrqsz
write(1234) qrqgz

close(1234)

open(1001,file=trim(diagfilename)//'.ctl')
write(1001,'(100a)')  'DSET    ./'//trim(diagfilename)//'.dat'
write(1001,'(100a)')  "TITLE   tangent linear check  "
write(1001,'(100a)')  "OPTIONS little_endian  "
write(1001,'(100a)')  "UNDEF   -999.0  "
write(1001,'(a10,i4,100a)')  "XDEF      ",nqbin,"  LINEAR    1 1  "
write(1001,'(a10,i4,100a)')  "YDEF      ",nqbin,"  LINEAR    1 1  "
write(1001,'(a10,i4,100a)')  "ZDEF      ",1,"  LINEAR    1 1  "
write(1001,'(a10,i10,100a)') "TDEF      ",1,"  LINEAR  07:00Z07jul2013        01MN  "
write(1001,'(100a)')  "VARS  5  "
write(1001,'(a,i4,100a)')  "qr          ",1," 99          qr    "
write(1001,'(a,i4,100a)')  "qs          ",1," 99          qs    "
write(1001,'(a,i4,100a)')  "qg          ",1," 99          qg    "
write(1001,'(a,i4,100a)')  "qrqs        ",1," 99          qrqs  "
write(1001,'(a,i4,100a)')  "qrqg        ",1," 99          qrqg  "
write(1001,'(100a)')  "ENDVARS  "
close(1001)

      
contains

subroutine readdata()

implicit none

character(len=256)  :: filename='./test_data.dat'
integer             :: nvar=14
integer             :: ivar

! 5, 6, 8, 9, 10, 11, 12, 13
! t, p, qv,qc,qr, qi, qs, qg

open(1001,file=trim(filename),form='binary')
if(print_init==1) print*,"filename:",trim(filename)
do ivar=1,nvar
  read(1001) temra01
  if(ivar==5) then
    tc=temra01
    if(print_init==1) print('(A,2E12.5)'),"tc max,min:",maxval(maxval(maxval(temra01(3:nx-2,3:ny-2,3:20),1),1),1), &
                         minval(minval(minval(temra01(3:nx-2,3:ny-2,3:20),1),1),1)
  endif
  if(ivar==6) then
    p=temra01
    if(print_init==1) print('(A,2E12.5)'),"p  max,min:",maxval(maxval(maxval(temra01(3:nx-2,3:ny-2,3:20),1),1),1), &
                         minval(minval(minval(temra01(3:nx-2,3:ny-2,3:20),1),1),1)
  endif
  if(ivar==8) then
    qv=temra01
    if(print_init==1) print('(A,2E12.5)'),"qv max,min:",maxval(maxval(maxval(temra01(3:nx-2,3:ny-2,3:20),1),1),1), &
                         minval(minval(minval(temra01(3:nx-2,3:ny-2,3:20),1),1),1)
  endif
  if(ivar==9) then
    qc=temra01
    if(print_init==1) print('(A,2E12.5)'),"qc max,min:",maxval(maxval(maxval(temra01(3:nx-2,3:ny-2,3:20),1),1),1), &
                         minval(minval(minval(temra01(3:nx-2,3:ny-2,3:20),1),1),1)
  endif
  if(ivar==10) then
    qr=temra01
    if(print_init==1) print('(A,2E12.5)'),"qr max,min:",maxval(maxval(maxval(temra01(3:nx-2,3:ny-2,3:20),1),1),1), &
                         minval(minval(minval(temra01(3:nx-2,3:ny-2,3:20),1),1),1)
  endif
  if(ivar==11) then
    qi=temra01
    if(print_init==1) print('(A,2E12.5)'),"qi max,min:",maxval(maxval(maxval(temra01(3:nx-2,3:ny-2,3:20),1),1),1), &
                         minval(minval(minval(temra01(3:nx-2,3:ny-2,3:20),1),1),1)
  endif
  if(ivar==12) then
    qs=temra01
    if(print_init==1) print('(A,2E12.5)'),"qs max,min:",maxval(maxval(maxval(temra01(3:nx-2,3:ny-2,3:20),1),1),1), &
                         minval(minval(minval(temra01(3:nx-2,3:ny-2,3:20),1),1),1)
  endif
  if(ivar==13) then
    qg=temra01
    if(print_init==1) print('(A,2E12.5)'),"qg max,min:",maxval(maxval(maxval(temra01(3:nx-2,3:ny-2,3:20),1),1),1), &
                         minval(minval(minval(temra01(3:nx-2,3:ny-2,3:20),1),1),1)
    exit
  endif
enddo

close(1001)
end subroutine readdata


FUNCTION GASDEV(IDUM)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Generates a random number from normal distribution by feeding
!  a negative integer iseed.
!
!  Added by M.Tong
!  Reference: Numerical Recipes
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    IDUM      an arbitrary negative integer as a seed for a
!              sequence of random numbers
!
!  OUTPUT:
!
!    GASDEV    A random number from Gaussian distribution with mean of 0
!              and standard deviation of 1.
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE        ! Force explicit declarations

  INTEGER :: IDUM        ! The seed for random number generation
  REAL :: GASDEV         ! The function to generate random number.
!
!-----------------------------------------------------------------------
!
!  Miscellaneous local variables:
!
!-----------------------------------------------------------------------
!
   INTEGER,SAVE::ISET
   REAL,SAVE::GSET
   REAL :: V1, V2, R
   REAL :: RAN1
   REAL :: FAC
   DATA ISET/0/
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
  IF (ISET.EQ.0) THEN
  1 V1=2.*RAN1(IDUM)-1.
    V2=2.*RAN1(IDUM)-1.
    R=V1**2+V2**2
    IF(R.GE.1.)GO TO 1
    FAC=SQRT(-2.*LOG(R)/R)
    GSET=V1*FAC
    GASDEV=V2*FAC
    ISET=1
  ELSE
    GASDEV=GSET
    ISET=0
  ENDIF

  RETURN
END FUNCTION GASDEV

end program linear_check

!
!##################################################################
!##################################################################
!######                                                      ######
!######                  FUNCTION RAN1                       ######
!######                                                      ######
!######                     Developed by                     ######
!######     Center for Analysis and Prediction of Storms     ######
!######                University of Oklahoma                ######
!######                                                      ######
!##################################################################
!##################################################################
!

  FUNCTION RAN1(IDUM)
!
!-----------------------------------------------------------------------
!
!  PURPOSE:
!
!  Generates a random number between 0 and 1 by feeding
!  a negative integer iseed.
!
!  Added by M.Tong
!  Reference: "Seminumerical Algorithms" by Donald Knuth
!
!-----------------------------------------------------------------------
!
!  INPUT :
!
!    IDUM      an arbitrary negative integer as a seed for a
!              sequence of random numbers
!
!  OUTPUT:
!
!    RAN1      A random number between 0 and 1.
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
!  Variable Declarations:
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE        ! Force explicit declarations
  INTEGER :: IDUM        ! The seed for random number generation
  REAL :: RAN1           ! The function to generate random number.
!
!-----------------------------------------------------------------------
!
!  Miscellaneous local variables:
!
!-----------------------------------------------------------------------
!
  REAL,SAVE :: R(97)
  INTEGER :: IX1,IX2,IX3,J,IFF
  INTEGER :: M1,M2,M3,IA1,IA2,IA3,IC1,IC2,IC3
  REAL :: RM1,RM2
  SAVE IX1,IX2,IX3

  PARAMETER (M1=259200,IA1=7141,IC1=54773,RM1=3.8580247E-6)
  PARAMETER (M2=134456,IA2=8121,IC2=28411,RM2=7.4373773E-6)
  PARAMETER (M3=243000,IA3=4561,IC3=51349)
  DATA IFF /0/
!
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!  Beginning of executable code...
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!
!----------------------------------------------------------------------
!
!  Initialize the sequence of random numbers between 0 and 1,
!  using iseed.
!
!----------------------------------------------------------------------
!
  IF (IDUM.LT.0.OR.IFF.EQ.0) THEN
    IFF=1
    IX1=MOD(IC1-IDUM,M1)
    IX1=MOD(IA1*IX1+IC1,M1)
    IX2=MOD(IX1,M2)
    IX1=MOD(IA1*IX1+IC1,M1)
    IX3=MOD(IX1,M3)
    DO J=1,97
      IX1=MOD(IA1*IX1+IC1,M1)
      IX2=MOD(IA2*IX2+IC2,M2)
      R(J)=(FLOAT(IX1)+FLOAT(IX2)*RM2)*RM1
    ENDDO
    IDUM=1
  ENDIF
  IX1=MOD(IA1*IX1+IC1,M1)
  IX2=MOD(IA2*IX2+IC2,M2)
  IX3=MOD(IA3*IX3+IC3,M3)
  J=1+(97*IX3)/M3
  IF(J.GT.97.OR.J.LT.1)THEN
    WRITE(*,*)'J is greater than 97 or less than 1','IDUM=',IDUM
    STOP
  ENDIF
  RAN1=R(J)
  R(J)=(FLOAT(IX1)+FLOAT(IX2)*RM2)*RM1

  RETURN
  END FUNCTION RAN1
