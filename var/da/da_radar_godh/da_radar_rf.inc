subroutine da_radar_rf (ra,rs,rg,tc1,ro1,rze)

   !-----------------------------------------------------------------------
   ! Purpose: calculate background/model radar equivalent reflectivity factor
   !-----------------------------------------------------------------------

   implicit none

   real, intent(in)    :: ra, rs, rg, tc1, ro1
   real, intent(inout) :: rze

   real                :: zrr, zds, zws, zgr
   real                :: zerr, zews, zeds, zegr

   if (trace_use) call da_trace_entry("da_radar_rf")

   ! ro1 use wrfda value, do not need to scale by 0.001

   zrr = 3.63*1.00e+9  ! rainwater
   zws = 4.26*1.00e+11 ! wet snow
   zds = 9.80*1.00e+8  ! dry snow
   zgr = 4.33*1.00e+10 ! graupel

   zerr = zrr*(ro1*max(0.0,ra))**1.75
   zews = zws*(ro1*max(0.0,rs))**1.75
   zeds = zds*(ro1*max(0.0,rs))**1.75
   zegr = zgr*(ro1*max(0.0,rg))**1.75

   if (tc1.lt.0.0) then
      rze = zerr + zeds + zegr
   else
      rze = zerr + zews + zegr
   end if

   if (rze.lt.1.0e-20) rze=1.0e-20

   if (trace_use) call da_trace_exit("da_radar_rf")

end subroutine da_radar_rf


!===================================================================
!      Following Jung et al., 2008
!      only for test
!===================================================================


subroutine RadZIceVar(qvp0,qra0,qsn0,qgr0,qnr0,qns0,qng0,tmk0,prs,dbz,                 &
                      miy,mjx,mkzh,in0r,in0s,in0g,rn0_r,rn0_s,rn0_g,                   &
                      rhos,rhog,dtmk,dqvp,dqra,dqsn,dqgr,dqnr,dqns,dqng,zmm,tlopt,     &
                      gropt,zmm_ref)

  implicit none
  
  integer,intent(in) :: miy,mjx,mkzh,in0r,in0s,in0g
  integer,intent(in) :: tlopt,gropt  ! tlopt=0 nonlinear,=1 tangent linear; gropt=0 tangent linear; >=1 adjoint

  integer,parameter  :: dualpol_opt=0 ! completed but not yet been tested

  real ::  rn0_r,rn0_s,rn0_g
  real ::  qra0(miy,mjx,mkzh),qsn0(miy,mjx,mkzh),qgr0(miy,mjx,mkzh), &
           qnr0(miy,mjx,mkzh),qns0(miy,mjx,mkzh),qng0(miy,mjx,mkzh), &
		   tmk0(miy,mjx,mkzh),qvp0(miy,mjx,mkzh)

  real ::  qra(miy,mjx,mkzh),qvp(miy,mjx,mkzh),                      &
           qsn(miy,mjx,mkzh),qgr(miy,mjx,mkzh),tmk(miy,mjx,mkzh),    &
           prs(miy,mjx,mkzh),dbz(miy,mjx,mkzh),qnr(miy,mjx,mkzh),    &
           qns(miy,mjx,mkzh),qng(miy,mjx,mkzh),ref(miy,mjx,mkzh),    &
           zrs(miy,mjx,mkzh),zss(miy,mjx,mkzh),zhs(miy,mjx,mkzh),    &
           zmm(miy,mjx,mkzh),zmm_ref(miy,mjx,mkzh),                  &
          ! above: background, below: increments
           dqra(miy,mjx,mkzh),dqsn(miy,mjx,mkzh),dqgr(miy,mjx,mkzh), &
           dqnr(miy,mjx,mkzh),dqns(miy,mjx,mkzh),dqng(miy,mjx,mkzh), &
		   dtmk(miy,mjx,mkzh),dqvp(miy,mjx,mkzh)

  real :: rgas=287.04
  real :: z_e

  character meth*1
  integer :: i,j,k,ii,jj,kk
      
!      include 'comconst'

!     ou operator constant
  real,parameter :: rdrwave = 107.0  ! unit mm S band
  real,parameter :: lambda = rdrwave
  real,parameter :: Kw2 = 0.93
  real           :: pi = 3.1415926
  real,parameter :: mm3todBZ=1.0E+9
!  real :: gamma  ! function
  real :: rhor = 1000   ! kg m^-3
  real :: rhos != 100.   ! kg m^-3
  real :: rhog != 100.   ! kg m^-3
  real :: rhoair

  real :: rhows,rhowg  !  wet snow, wet graupel

!     temporal mixing ratio
  real :: prain  ! pure rain mixing ratio
  real :: dsnow  ! dry snow mixing ratio
  real :: dgr    ! dry graupel mixing ration 
  real :: wsnow  ! wet snow mixing ratio
  real :: wgr    ! wet graupel mixing ratio
!  real :: upper_f ! function of Fraction fmax*min(qice/qr,qr/qice)**0.3

!     parameters for rain            
  real :: alpha_ra = 4.28e-4
  real :: alpha_rb = 4.28e-4
  real :: beta_ra = 3.04
  real :: beta_rb = 2.77
  real :: alphar = 0
  real :: dr = 3
  real :: cr = 3.1415926/6
  real :: zrh,zrv
     
  real :: para10r  ! =(1.-(2.*betarx+1.0)/4.)
  real :: para14r  ! mm3todBZ*(4*lambda**4*alpharx**2/(pi**4*Kw2))*(pi*rhor/rhoa)**(p2/4.)**(n0r)**para10r*gamma(-(2.*betarx+1.0)*1d0)
  real :: ronv,gonv,sonv

!     parameters for snow graupel or hail

  integer,parameter :: npara_alpharxa=7           
  real :: para_alpha_rxa(npara_alpharxa,2) ! second dimension: 1 for snow 
!                                            2 for hail/graupel
  real,save :: para_alpha_rxb(npara_alpharxa,2)
  real :: phimean=0,sigma,ice_abc(3),ice_bac(3)
  real :: fw
  real :: pxabk_all(3)
  real :: pxkh,pxkv
  real :: pxkh_tlr,pxkv_tlr,pxkh_tlx,pxkv_tlx
  real :: zsh,zsv,zgh,zgv,zdsh,zdsv,zdgh,zdgv

  real :: para1sg   ! para1sg=mm3todBZ*gamma(7)*lambda**4/pi**4/Kw2
  real :: para2sg   ! para2sg=(pi*rhox/rhoa)**1.75

  real :: cs=3.1415926/6
  real :: cice=440. 
  real :: ds=3.
  real :: dg=3. 
  real :: alphas = 0
  real :: alphag = 0

  real :: alpha_rdsa=0.194*10.**(-4)
  real :: alpha_rdsb=0.191*10.**(-4)
  real :: alpha_rdha=0.191*10.**(-3)
  real :: alpha_rdhb=0.165*10.**(-3)

  real :: zh_e,zv_e
  
  REAL,PARAMETER :: sigmahd = 1.0472
  REAL,PARAMETER :: Ahd = 0.4308
  REAL,PARAMETER :: Bhd = 0.3192
  REAL,PARAMETER :: Chd = 0.1250
  REAL,PARAMETER :: Dhd = 0.3750
  REAL,PARAMETER :: Ckhd = 0.1116
  ! for dry snow
  REAL,PARAMETER :: sigmas = 0.3491
  REAL,PARAMETER :: Asd = 0.8140
  REAL,PARAMETER :: Bsd = 0.0303
  REAL,PARAMETER :: Csd = 0.0778
  REAL,PARAMETER :: Dsd = 0.4221
  REAL,PARAMETER :: Cksd = 0.7837

  real           :: ice_abc_d(3),ice_bac_d(3)
  real :: pdfrrs,pdfrrg,pdfsrs,pdfgrg   
  ! partial derivative of F (r&s,r&g) with respect to (r,s,g) 
  real :: pdfrhot,pdfrhoq
  ! partial derivative of rho for t and qv
  
  real :: prain_coef,dsnow_coef,dgr_coef
  real :: qthres
!  real :: fw_tl_thres
!  real :: qinv_thres

  real :: temr01,temr02,temr03
  integer :: temi01,temi02,temi03
  
  integer :: print_chk=0  

  data para_alpha_rxa(:,1)/0.194e-4,7.094e-4,2.135e-4,-5.225e-4,0,0,0/
  data para_alpha_rxb(:,1)/0.191e-4,6.916e-4,-2.841e-4,-1.160e-4,0,0,0/
  data para_alpha_rxa(:,2)/0.191e-3,2.39e-3,-12.57e-3,38.71e-3,-65.53e-3,56.16e-3,  &
                          -18.98e-3/
  data para_alpha_rxb(:,2)/0.165e-3,1.72e-3,-9.92e-3,32.15e-3,-56.0e-3,48.83e-3,   &
                          -16.69e-3/
  qthres=rf_qthres


!  print*,"RadZIceVar test"
!
  qra=qra0
  qsn=qsn0
  qgr=qgr0
  tmk=tmk0
  qvp=qvp0
  do k=1,mkzh
    do j=1,mjx
      do i=1,miy
        if(qra0(i,j,k)<qthres) qra(i,j,k)=qthres
        if(qsn0(i,j,k)<qthres) qsn(i,j,k)=qthres
        if(qgr0(i,j,k)<qthres) qgr(i,j,k)=qthres
      enddo
    enddo
  enddo

  if(tlopt>=1.and.gropt>=1) then
    dqra=0
    dqgr=0
    dqsn=0
    dqvp=0
    dtmk=0
    dqnr=0
    dqns=0
    dqng=0
  endif

  do k=1,mkzh
    do j=1,mjx
      do i=1,miy
!
         zrh=0
		 zsh=0
		 zgh=0
		 zdsh=0
		 zdgh=0
		 zrv=0
		 zsv=0
		 zgv=0
		 zdsv=0
		 zdgv=0

         if(qra(i,j,k)<qthres.and.qsn(i,j,k)<qthres.and.qgr(i,j,k)<qthres) then
           zrs(i,j,k)=0
           zss(i,j,k)=0
           zhs(i,j,k)=0
           ref(i,j,k)=0
           zmm(i,j,k)=0
           dbz(i,j,k)=0
           cycle
         endif

         rhoair=prs(i,j,k)*1./    &
           (rgas*virtual(tmk(i,j,k),qvp(i,j,k)))      ! air density
         if(tlopt>=1) then
		   pdfrhot=-prs(i,j,k)*1./    &
            (rgas*virtual(tmk(i,j,k),qvp(i,j,k)))*(1./tmk(i,j,k))
		   pdfrhoq=-prs(i,j,k)*1./    &
            (rgas*virtual(tmk(i,j,k),qvp(i,j,k)))*(0.61/(1+0.61*qvp(i,j,k)))		   	
		 endif
        pdfrhot=0
        pdfrhoq=0		 
        print_chk=0

!
!      Calculate variable intercept parameters if wanted
!
         if (in0s.eq.1) then  
            sonv = rn0_s  ! for diagnosing n0 from number concentration, not yet completed
         else
            sonv = rn0_s
         endif
!
         if (in0g.eq.1) then  
            gonv = rn0_g  ! for diagnosing n0 from number concentration, not yet completed
         else
            gonv = rn0_g
         endif
!
         if (in0r.eq.1) then  
            ronv = rn0_r  ! for diagnosing n0 from number concentration, not yet completed
         else
            ronv = rn0_r
         endif

         ronv=ronv*1.e-12   ! convert unit
         gonv=gonv*1.e-12   ! convert unit
         sonv=sonv*1.e-12   ! convert unit
 
         prain_coef=1-upper_f(qra(i,j,k),qgr(i,j,k),qthres,2)-upper_f(qra(i,j,k),qsn(i,j,k),qthres,1)
         if(prain_coef<0.0) then
           prain_coef=0
         endif     
         prain=prain_coef*qra(i,j,k)

		 dsnow_coef=1-upper_f(qra(i,j,k),qsn(i,j,k),qthres,1)
         if(dsnow_coef<0.0) then
           dsnow_coef=0
         endif
         dsnow=dsnow_coef*qsn(i,j,k)
         wsnow=(1-dsnow_coef)*(qsn(i,j,k)+qra(i,j,k))

		 dgr_coef=1-upper_f(qra(i,j,k),qgr(i,j,k),qthres,2)
         if(dgr_coef<0.0) then
           dgr_coef=0.0
         endif     
         dgr=dgr_coef*qgr(i,j,k)
         wgr=(1-dgr_coef)*(qgr(i,j,k)+qra(i,j,k))

         if(tlopt>=1) then
		   pdfrrs=0
		   pdfrrg=0
		   pdfsrs=0
		   pdfgrg=0 
         endif	
		 
         if(tlopt>=1.and.gropt==2) then
            if(zmm(i,j,k)>1.0) then
              zmm_ref(i,j,k)=10./(zmm(i,j,k)*log(10.0))*zmm_ref(i,j,k)
            else
              zmm_ref(i,j,k)=0.0
            endif
         endif
		 
!    ==================FOR RAIN=============================
         !print*,"----------------for rain----------------------"
         call parameter_zrx(para10r,para14r,1.0,rhor,cr,dr,alphar,     &
                            beta_ra,alpha_ra,mm3todBZ,lambda,Kw2,pi,   &
                            ronv)

         if(in0r.eq.1) then
           ! for two moment microphysics scheme, not yet completed
         else
           if(tlopt==0) zrh=para14r*(rhoair*prain)**(1.-para10r)
           if(tlopt>=1.and.gropt==0) zrh=para14r*(1.-para10r)*(rhoair*prain)**(-para10r)*rhoair*prain_coef*dqra(i,j,k)                &
		                                +para14r*(1.-para10r)*(rhoair*prain)**(-para10r)*prain*pdfrhot*dtmk(i,j,k)                    &
										+para14r*(1.-para10r)*(rhoair*prain)**(-para10r)*prain*pdfrhoq*dqvp(i,j,k)                    &
		                                -para14r*(1.-para10r)*(rhoair*prain)**(-para10r)*(pdfrrs+pdfrrg)*qra(i,j,k)*dqra(i,j,k)       &
										-para14r*(1.-para10r)*(rhoair*prain)**(-para10r)*pdfsrs*qra(i,j,k)*dqsn(i,j,k)                &
										-para14r*(1.-para10r)*(rhoair*prain)**(-para10r)*pdfgrg*qra(i,j,k)*dqgr(i,j,k)
           if(tlopt>=1.and.gropt>=1) then
             ! here, zmm_ref serves as the HX (or yb) for gradient calculation
             dqra(i,j,k)=dqra(i,j,k)+para14r*(1.-para10r)*(rhoair*prain)**(-para10r)*rhoair*prain_coef*zmm_ref(i,j,k) &
			                        -para14r*(1.-para10r)*(rhoair*prain)**(-para10r)*(pdfrrs+pdfrrg)*qra(i,j,k)*zmm_ref(i,j,k)
             dqsn(i,j,k)=dqsn(i,j,k)-para14r*(1.-para10r)*(rhoair*prain)**(-para10r)*pdfsrs*qra(i,j,k)*zmm_ref(i,j,k)
             dqgr(i,j,k)=dqgr(i,j,k)-para14r*(1.-para10r)*(rhoair*prain)**(-para10r)*pdfgrg*qra(i,j,k)*zmm_ref(i,j,k)		
             dtmk(i,j,k)=dtmk(i,j,k)+para14r*(1.-para10r)*(rhoair*prain)**(-para10r)*prain*pdfrhot*zmm_ref(i,j,k)
             dqvp(i,j,k)=dqvp(i,j,k)+para14r*(1.-para10r)*(rhoair*prain)**(-para10r)*prain*pdfrhoq*zmm_ref(i,j,k)			 
           endif
         endif
         if(print_chk>=2) then
           print*,"zrh,10.*log10(zrh)"
           print*,zrh,max(0.,10.*log10(zrh)),dqra(i,j,k),prain_coef
         endif

         if(dualpol_opt==1) then
           call parameter_zrx(para10r,para14r,1.0,rhor,cr,dr,alphar,     &
                              beta_rb,alpha_rb,mm3todBZ,lambda,Kw2,pi,   &
                              ronv)


           if(in0r.eq.1) then
             ! for two moment microphysics scheme, not yet completed
           else
             if(tlopt==0) zrv=para14r*(rhoair*prain)**(1.-para10r)
             if(tlopt>=1.and.gropt==0) zrv=para14r*(1.-para10r)*(rhoair*prain)**(-para10r)*rhoair*prain_coef*dqra(i,j,k)          &
                                          +para14r*(1.-para10r)*(rhoair*prain)**(-para10r)*prain*pdfrhot*dtmk(i,j,k)              &
                                          +para14r*(1.-para10r)*(rhoair*prain)**(-para10r)*prain*pdfrhoq*dqvp(i,j,k)              &		   
                                          -para14r*(1.-para10r)*(rhoair*prain)**(-para10r)*(pdfrrs+pdfrrg)*qra(i,j,k)*dqra(i,j,k) &
                                          -para14r*(1.-para10r)*(rhoair*prain)**(-para10r)*pdfsrs*qra(i,j,k)*dqsn(i,j,k)          &
                                          -para14r*(1.-para10r)*(rhoair*prain)**(-para10r)*pdfgrg*qra(i,j,k)*dqgr(i,j,k)
             if(tlopt>=1.and.gropt>=1) then
             ! here, zmm_ref serves as the HX (or yb) for gradient calculation (not available yet)
               dqra(i,j,k)=dqra(i,j,k)+para14r*(1.-para10r)*(rhoair*prain)**(-para10r)*rhoair*prain_coef*zmm_ref(i,j,k)   &
                                      -para14r*(1.-para10r)*(rhoair*prain)**(-para10r)*(pdfrrs+pdfrrg)*qra(i,j,k)*zmm_ref(i,j,k)
               dqsn(i,j,k)=dqsn(i,j,k)-para14r*(1.-para10r)*(rhoair*prain)**(-para10r)*pdfsrs*qra(i,j,k)*zmm_ref(i,j,k)
               dqgr(i,j,k)=dqgr(i,j,k)-para14r*(1.-para10r)*(rhoair*prain)**(-para10r)*pdfgrg*qra(i,j,k)*zmm_ref(i,j,k)	
			   dtmk(i,j,k)=dtmk(i,j,k)+para14r*(1.-para10r)*(rhoair*prain)**(-para10r)*prain*pdfrhot*zmm_ref(i,j,k)
               dqvp(i,j,k)=dqvp(i,j,k)+para14r*(1.-para10r)*(rhoair*prain)**(-para10r)*prain*pdfrhoq*zmm_ref(i,j,k)		
             endif
           endif
         endif
         zrs(i,j,k)=zrh 
         if(print_chk>=2) then
           print*,"zrv,10.*log10(zrv)"  
           print*,zrv,max(0.,10.*log10(zrv) ),dqra(i,j,k),prain_coef
         endif
!    ==================FOR snow=============================
       if(rf_noice==0) then
        ! print*,"----------------for snow----------------------"
         fw=waterfraction(qra(i,j,k),qsn(i,j,k))
         rhows=rhos*(1.-fw**2)+rhor*fw**2

         call parameter_zxx(para1sg,para2sg,1.0,rhows,                   &
                            cs,ds,alphas,mm3todBZ,lambda,Kw2,pi,         &
                            sonv)

         sigma=sigma_in_abc(qsn(i,j,k),fw,1)  ! for snow
        ! print*,"sigma in abc for snow"
        ! print*,sigma

         call calc_ice_abc(phimean,sigma,ice_abc)
        ! print*,"A,B,C:"
        ! print*,ice_abc

         ice_bac(1)=ice_abc(2)
         ice_bac(2)=ice_abc(1)
         ice_bac(3)=ice_abc(3)
         pxkh=0
         pxkv=0
         pxkh_tlr=0
         pxkv_tlr=0
         pxkh_tlx=0
         pxkv_tlx=0
         do kk=0,2*npara_alpharxa-1
           pxabk_all(1)=pxabk(para_alpha_rxa(:,1),para_alpha_rxa(:,1), &
                             kk,npara_alpharxa)
           pxabk_all(2)=pxabk(para_alpha_rxb(:,1),para_alpha_rxb(:,1), &
                             kk,npara_alpharxa)
           pxabk_all(3)=pxabk(para_alpha_rxa(:,1),para_alpha_rxb(:,1), &
                             kk,npara_alpharxa)
           pxkh=pxkh+pkx(ice_abc,pxabk_all)*fw**kk
           pxkv=pxkv+pkx(ice_bac,pxabk_all)*fw**kk
           if(tlopt==2.and.kk>=1) then
             pxkh_tlr=pxkh_tlr+pkx(ice_abc,pxabk_all)*kk*fw**kk*(1./qra(i,j,k)-1./(qra(i,j,k)+qsn(i,j,k)))
             pxkv_tlr=pxkv_tlr+pkx(ice_bac,pxabk_all)*kk*fw**kk*(1./qra(i,j,k)-1./(qra(i,j,k)+qsn(i,j,k)))
             pxkh_tlx=pxkh_tlx+pkx(ice_abc,pxabk_all)*kk*fw**kk*(-1./(qra(i,j,k)+qsn(i,j,k)))
             pxkv_tlx=pxkv_tlx+pkx(ice_bac,pxabk_all)*kk*fw**kk*(-1./(qra(i,j,k)+qsn(i,j,k)))
           endif
         enddo
        ! print*,"pxkh,pxkv"
        ! print*,pxkh,pxkv


         if(in0s.eq.1) then
           ! for two moment microphysics scheme, not yet completed
         else
           if(tlopt==0) zsh= para1sg*para2sg*sonv**(-0.75)*(rhoair*wsnow)**(1.75)*pxkh
           if(tlopt>=1.and.gropt==0) zsh= 1.75*para1sg*para2sg*sonv**(-0.75)*(rhoair*wsnow)**(0.75)*pxkh*rhoair*(1-dsnow_coef)*dqsn(i,j,k)   &
		                    +1.75*para1sg*para2sg*sonv**(-0.75)*(rhoair*wsnow)**(0.75)*pxkh*rhoair*(1-dsnow_coef)*dqra(i,j,k)                &
							+1.75*para1sg*para2sg*sonv**(-0.75)*(rhoair*wsnow)**(0.75)*pxkh*wsnow*pdfrhot*dtmk(i,j,k)                        &
							+1.75*para1sg*para2sg*sonv**(-0.75)*(rhoair*wsnow)**(0.75)*pxkh*wsnow*pdfrhoq*dqvp(i,j,k)                        &
		                    +1.75*para1sg*para2sg*sonv**(-0.75)*(rhoair*wsnow)**(0.75)*pxkh*pdfrrs*qsn(i,j,k)*dqra(i,j,k)                    &
							+1.75*para1sg*para2sg*sonv**(-0.75)*(rhoair*wsnow)**(0.75)*pxkh*pdfsrs*qsn(i,j,k)*dqsn(i,j,k)                    &
                            +para1sg*para2sg*sonv**(-0.75)*(rhoair*wsnow)**(1.75)*pxkh_tlr*dqra(i,j,k)                                       &
                            +para1sg*para2sg*sonv**(-0.75)*(rhoair*wsnow)**(1.75)*pxkh_tlx*dqsn(i,j,k)     
           if(tlopt>=1.and.gropt>=1) then
             ! here, zmm_ref serves as the HX (or yb) for gradient calculation
             dqsn(i,j,k)=dqsn(i,j,k)+1.75*para1sg*para2sg*sonv**(-0.75)*(rhoair*wsnow)**(0.75)*pxkh*rhoair*(1-dsnow_coef)*zmm_ref(i,j,k) &
			                        +1.75*para1sg*para2sg*sonv**(-0.75)*(rhoair*wsnow)**(0.75)*pxkh*pdfsrs*qsn(i,j,k)*zmm_ref(i,j,k)
             dqra(i,j,k)=dqra(i,j,k)+para1sg*para2sg*sonv**(-0.75)*(rhoair*wsnow)**(1.75)*pxkh_tlr*zmm_ref(i,j,k) &
			                        +1.75*para1sg*para2sg*sonv**(-0.75)*(rhoair*wsnow)**(0.75)*pxkh*pdfrrs*qsn(i,j,k)*zmm_ref(i,j,k) &
									+1.75*para1sg*para2sg*sonv**(-0.75)*(rhoair*wsnow)**(0.75)*pxkh*rhoair*(1-dsnow_coef)*zmm_ref(i,j,k)
             dqsn(i,j,k)=dqsn(i,j,k)+para1sg*para2sg*sonv**(-0.75)*(rhoair*wsnow)**(1.75)*pxkh_tlx*zmm_ref(i,j,k)
			 dtmk(i,j,k)=dtmk(i,j,k)+1.75*para1sg*para2sg*sonv**(-0.75)*(rhoair*wsnow)**(0.75)*pxkh*wsnow*pdfrhot*zmm_ref(i,j,k)
			 dqvp(i,j,k)=dqvp(i,j,k)+1.75*para1sg*para2sg*sonv**(-0.75)*(rhoair*wsnow)**(0.75)*pxkh*wsnow*pdfrhoq*zmm_ref(i,j,k)
           endif
         endif
         if(print_chk>=2) then
           print*,"zsh,10.*log10(zsh)"
           print*,zsh,max(0.,10.*log10(zsh)),dqsn(i,j,k),dsnow_coef
         endif

         if(dualpol_opt==1) then
           if(in0s.eq.1) then
             ! for two moment microphysics scheme, not yet completed
           else
             if(tlopt==0) zsv= para1sg*para2sg*sonv**(-0.75)*(rhoair*wsnow)**(1.75)*pxkv
             if(tlopt>=1.and.gropt==0) zsv= 1.75*para1sg*para2sg*sonv**(-0.75)*(rhoair*wsnow)**(0.75)*pxkv*rhoair*(1-dsnow_coef)*dqsn(i,j,k)   &
                                          +1.75*para1sg*para2sg*sonv**(-0.75)*(rhoair*wsnow)**(0.75)*pxkv*rhoair*(1-dsnow_coef)*dqra(i,j,k)    &
                                          +1.75*para1sg*para2sg*sonv**(-0.75)*(rhoair*wsnow)**(0.75)*pxkh*wsnow*pdfrhot*dtmk(i,j,k)            &
                                          +1.75*para1sg*para2sg*sonv**(-0.75)*(rhoair*wsnow)**(0.75)*pxkh*wsnow*pdfrhoq*dqvp(i,j,k)            &							
                                          +1.75*para1sg*para2sg*sonv**(-0.75)*(rhoair*wsnow)**(0.75)*pxkv*pdfrrs*qsn(i,j,k)*dqra(i,j,k)        &
                                          +1.75*para1sg*para2sg*sonv**(-0.75)*(rhoair*wsnow)**(0.75)*pxkv*pdfsrs*qsn(i,j,k)*dqsn(i,j,k)        &		   
                                          +para1sg*para2sg*sonv**(-0.75)*(rhoair*wsnow)**(1.75)*pxkv_tlr*dqra(i,j,k)                           &
                                          +para1sg*para2sg*sonv**(-0.75)*(rhoair*wsnow)**(1.75)*pxkv_tlx*dqsn(i,j,k)
             if(tlopt>=1.and.gropt>=1) then
               ! here, zmm_ref serves as the HX (or yb) for gradient calculation(not available yet)
               dqsn(i,j,k)=dqsn(i,j,k)+1.75*para1sg*para2sg*sonv**(-0.75)*(rhoair*wsnow)**(0.75)*pxkv*rhoair*(1-dsnow_coef)*zmm_ref(i,j,k) &
                                      +1.75*para1sg*para2sg*sonv**(-0.75)*(rhoair*wsnow)**(0.75)*pxkv*pdfsrs*qsn(i,j,k)*zmm_ref(i,j,k)
               dqra(i,j,k)=dqra(i,j,k)+para1sg*para2sg*sonv**(-0.75)*(rhoair*wsnow)**(1.75)*pxkv_tlr*zmm_ref(i,j,k)                        &
                                      +1.75*para1sg*para2sg*sonv**(-0.75)*(rhoair*wsnow)**(0.75)*pxkv*pdfrrs*qsn(i,j,k)*zmm_ref(i,j,k)     &
                                      +1.75*para1sg*para2sg*sonv**(-0.75)*(rhoair*wsnow)**(0.75)*pxkv*rhoair*(1-dsnow_coef)*zmm_ref(i,j,k)
               dqsn(i,j,k)=dqsn(i,j,k)+para1sg*para2sg*sonv**(-0.75)*(rhoair*wsnow)**(1.75)*pxkv_tlx*zmm_ref(i,j,k)
			   dtmk(i,j,k)=dtmk(i,j,k)+1.75*para1sg*para2sg*sonv**(-0.75)*(rhoair*wsnow)**(0.75)*pxkh*wsnow*pdfrhot*zmm_ref(i,j,k)
			   dqvp(i,j,k)=dqvp(i,j,k)+1.75*para1sg*para2sg*sonv**(-0.75)*(rhoair*wsnow)**(0.75)*pxkh*wsnow*pdfrhoq*zmm_ref(i,j,k)			 
             endif
           endif
         endif
         if(print_chk>=2) then
           print*,"zsv,10.*log10(zsv)"
           print*,zsv,max(0.,10.*log10(zsv)),dqsn(i,j,k),dsnow_coef
         endif
       !  print*,"----------------for dry snow----------------------"
         call parameter_zxx(para1sg,para2sg,1.0,rhos,                    &
                            cs,ds,alphas,mm3todBZ,lambda,Kw2,pi,         &
                            sonv)

         ice_abc_d(1)=Asd
         ice_abc_d(2)=Bsd
         ice_abc_d(3)=Csd
         ice_bac_d(1)=Bsd
         ice_bac_d(2)=Asd						
						
         pxabk_all(1)=alpha_rdsa**2
         pxabk_all(2)=alpha_rdsb**2
         pxabk_all(3)=alpha_rdsb*alpha_rdsa
         pxkh=pkx(ice_abc_d,pxabk_all)
         pxkv=pkx(ice_bac_d,pxabk_all)
        ! print*,"pxkh,pxkv"
        ! print*,pxkh,pxkv


         if(in0s.eq.1) then
           ! for two moment microphysics scheme, not yet completed        
         else
           if(tlopt==0) zdsh=para1sg*para2sg*sonv**(-0.75)*(rhoair*dsnow)**(1.75)*pxkh
           if(tlopt>=1.and.gropt==0) zdsh=1.75*para1sg*para2sg*sonv**(-0.75)*(rhoair*dsnow)**(0.75)*pxkh*rhoair*dsnow_coef*dqsn(i,j,k)          &
		                                 +1.75*para1sg*para2sg*sonv**(-0.75)*(rhoair*dsnow)**(0.75)*pxkh*dsnow*pdfrhot*dtmk(i,j,k)             & 
										 +1.75*para1sg*para2sg*sonv**(-0.75)*(rhoair*dsnow)**(0.75)*pxkh*dsnow*pdfrhoq*dqvp(i,j,k)             & 
		                                 -1.75*para1sg*para2sg*sonv**(-0.75)*(rhoair*dsnow)**(0.75)*pxkh*pdfrrs*qsn(i,j,k)*dqra(i,j,k)   &
										 -1.75*para1sg*para2sg*sonv**(-0.75)*(rhoair*dsnow)**(0.75)*pxkh*pdfsrs*qsn(i,j,k)*dqsn(i,j,k)
           if(tlopt>=1.and.gropt>=1) then
             ! here, zmm_ref serves as the HX (or yb) for gradient calculation
             dqsn(i,j,k)=dqsn(i,j,k)+1.75*para1sg*para2sg*sonv**(-0.75)*(rhoair*dsnow)**(0.75)*pxkh*rhoair*dsnow_coef*zmm_ref(i,j,k) &
			                        -1.75*para1sg*para2sg*sonv**(-0.75)*(rhoair*dsnow)**(0.75)*pxkh*pdfsrs*qsn(i,j,k)*zmm_ref(i,j,k)
             dqra(i,j,k)=dqra(i,j,k)-1.75*para1sg*para2sg*sonv**(-0.75)*(rhoair*dsnow)**(0.75)*pxkh*pdfrrs*qsn(i,j,k)*zmm_ref(i,j,k)
             dtmk(i,j,k)=dtmk(i,j,k)+1.75*para1sg*para2sg*sonv**(-0.75)*(rhoair*dsnow)**(0.75)*pxkh*dsnow*pdfrhot*zmm_ref(i,j,k)
             dqvp(i,j,k)=dqvp(i,j,k)+1.75*para1sg*para2sg*sonv**(-0.75)*(rhoair*dsnow)**(0.75)*pxkh*dsnow*pdfrhoq*zmm_ref(i,j,k)			 
           endif 
         endif
         if(print_chk>=2) then
           print*,"zdsh,10.*log10(zdsh)"
           print*,zdsh,max(0.,10.*log10(zdsh)),dqsn(i,j,k),dsnow_coef
         endif

         if(dualpol_opt==1) then
           if(in0s.eq.1) then
             ! for two moment microphysics scheme, not yet completed
           else
             if(tlopt==0) zdsv=para1sg*para2sg*sonv**(-0.75)*(rhoair*dsnow)**(1.75)*pxkv
             if(tlopt>=1.and.gropt==0) zdsv=1.75*para1sg*para2sg*sonv**(-0.75)*(rhoair*dsnow)**(0.75)*pxkv*rhoair*dsnow_coef*dqsn(i,j,k)         &
		                                   +1.75*para1sg*para2sg*sonv**(-0.75)*(rhoair*dsnow)**(0.75)*pxkh*dsnow*pdfrhot*dtmk(i,j,k)             & 
										   +1.75*para1sg*para2sg*sonv**(-0.75)*(rhoair*dsnow)**(0.75)*pxkh*dsnow*pdfrhoq*dqvp(i,j,k)             & 		   
		                                   -1.75*para1sg*para2sg*sonv**(-0.75)*(rhoair*dsnow)**(0.75)*pxkv*pdfrrs*qsn(i,j,k)*dqra(i,j,k)         &
										   -1.75*para1sg*para2sg*sonv**(-0.75)*(rhoair*dsnow)**(0.75)*pxkv*pdfsrs*qsn(i,j,k)*dqsn(i,j,k) 
             if(tlopt>=1.and.gropt>=1) then
               ! here, zmm_ref serves as the HX (or yb) for gradient calculation(not available yet)
               dqsn(i,j,k)=dqsn(i,j,k)+1.75*para1sg*para2sg*sonv**(-0.75)*(rhoair*dsnow)**(0.75)*pxkv*rhoair*dsnow_coef*zmm_ref(i,j,k) &
			                          -1.75*para1sg*para2sg*sonv**(-0.75)*(rhoair*dsnow)**(0.75)*pxkv*pdfsrs*qsn(i,j,k)*zmm_ref(i,j,k)
			   dqra(i,j,k)=dqra(i,j,k)-1.75*para1sg*para2sg*sonv**(-0.75)*(rhoair*dsnow)**(0.75)*pxkv*pdfrrs*qsn(i,j,k)*zmm_ref(i,j,k)
               dtmk(i,j,k)=dtmk(i,j,k)+1.75*para1sg*para2sg*sonv**(-0.75)*(rhoair*dsnow)**(0.75)*pxkh*dsnow*pdfrhot*zmm_ref(i,j,k)
               dqvp(i,j,k)=dqvp(i,j,k)+1.75*para1sg*para2sg*sonv**(-0.75)*(rhoair*dsnow)**(0.75)*pxkh*dsnow*pdfrhoq*zmm_ref(i,j,k)	
             endif
           endif
         endif
         zss(i,j,k)=zsh+zdsh
         if(print_chk>=2) then
           print*,"zdsv,10.*log10(zdsv)"
           print*,zdsv,max(0.,10.*log10(zdsv)),dqsn(i,j,k),dsnow_coef
         endif

!    ==================FOR graupel==========================
        ! print*,"----------------for graupel----------------------"

         fw=waterfraction(qra(i,j,k),qgr(i,j,k))
         rhowg=rhog*(1.-fw**2)+rhor*fw**2

         call parameter_zxx(para1sg,para2sg,1.0,rhowg,                    &
                            cice,dg,alphag,mm3todBZ,lambda,Kw2,pi,        &
                            gonv)

         sigma=sigma_in_abc(qgr(i,j,k),fw,2)  ! for snow
       !  print*,"sigma in abc for snow"
       !  print*,sigma

         call calc_ice_abc(phimean,sigma,ice_abc)
       !  print*,"A,B,C:"
       !  print*,ice_abc

         ice_bac(1)=ice_abc(2)
         ice_bac(2)=ice_abc(1)
         ice_bac(3)=ice_abc(3)
         pxkh=0
         pxkv=0
         pxkh_tlr=0
         pxkv_tlr=0
         pxkh_tlx=0
         pxkv_tlx=0
         do kk=0,2*npara_alpharxa-1
           pxabk_all(1)=pxabk(para_alpha_rxa(:,2),para_alpha_rxa(:,2),  &
                             kk,npara_alpharxa)
           pxabk_all(2)=pxabk(para_alpha_rxb(:,2),para_alpha_rxb(:,2),  &
                             kk,npara_alpharxa)
           pxabk_all(3)=pxabk(para_alpha_rxa(:,2),para_alpha_rxb(:,2),  &
                             kk,npara_alpharxa)
           pxkh=pxkh+pkx(ice_abc,pxabk_all)*fw**kk
           pxkv=pxkv+pkx(ice_bac,pxabk_all)*fw**kk
           if(tlopt==2.and.kk>=1) then
             pxkh_tlr=pxkh_tlr+pkx(ice_abc,pxabk_all)*kk*fw**kk*(1./qra(i,j,k)-1./(qra(i,j,k)+qgr(i,j,k)))
             pxkv_tlr=pxkv_tlr+pkx(ice_bac,pxabk_all)*kk*fw**kk*(1./qra(i,j,k)-1./(qra(i,j,k)+qgr(i,j,k)))
             pxkh_tlx=pxkh_tlx+pkx(ice_abc,pxabk_all)*kk*fw**kk*(-1./(qra(i,j,k)+qgr(i,j,k)))
             pxkv_tlx=pxkv_tlx+pkx(ice_bac,pxabk_all)*kk*fw**kk*(-1./(qra(i,j,k)+qgr(i,j,k)))
             if(print_chk>=2.and.kk==2*npara_alpharxa-1) print*,pxkh_tlr,pxkv_tlr,pxkh_tlx,pxkv_tlx,pxkh,pxkv
          endif
         enddo
       !  print*,"pxkh,pxkv"
       !  print*,pxkh,pxkv

         if(in0g.eq.1) then
           ! for two moment microphysics scheme, not yet completed
         else
           if(tlopt==0) zgh= para1sg*para2sg*gonv**(-0.75)*(rhoair*wgr)**(1.75)*pxkh
           if(tlopt>=1.and.gropt==0) zgh= 1.75*para1sg*para2sg*gonv**(-0.75)*(rhoair*wgr)**(0.75)*pxkh*rhoair*(1-dgr_coef)*dqgr(i,j,k)   &
		                    +1.75*para1sg*para2sg*gonv**(-0.75)*(rhoair*wgr)**(0.75)*pxkh*rhoair*(1-dgr_coef)*dqra(i,j,k)                &
							+1.75*para1sg*para2sg*gonv**(-0.75)*(rhoair*wgr)**(0.75)*pxkh*wgr*pdfrhot*dtmk(i,j,k)                        &
							+1.75*para1sg*para2sg*gonv**(-0.75)*(rhoair*wgr)**(0.75)*pxkh*wgr*pdfrhoq*dqvp(i,j,k)                        &
		                    +1.75*para1sg*para2sg*gonv**(-0.75)*(rhoair*wgr)**(0.75)*pxkh*pdfrrg*qgr(i,j,k)*dqra(i,j,k)                  & 
							+1.75*para1sg*para2sg*gonv**(-0.75)*(rhoair*wgr)**(0.75)*pxkh*pdfgrg*qgr(i,j,k)*dqgr(i,j,k)                  & 
                            +para1sg*para2sg*gonv**(-0.75)*(rhoair*wgr)**(1.75)*pxkh_tlr*dqra(i,j,k)                                     &
                            +para1sg*para2sg*gonv**(-0.75)*(rhoair*wgr)**(1.75)*pxkh_tlx*dqgr(i,j,k)
           if(tlopt>=1.and.gropt>=1) then
             ! here, zmm_ref serves as the HX (or yb) for gradient calculation
             dqgr(i,j,k)=dqgr(i,j,k)+1.75*para1sg*para2sg*gonv**(-0.75)*(rhoair*wgr)**(0.75)*pxkh*rhoair*(1-dgr_coef)*zmm_ref(i,j,k) &
			                        +1.75*para1sg*para2sg*gonv**(-0.75)*(rhoair*wgr)**(0.75)*pxkh*pdfgrg*qgr(i,j,k)*zmm_ref(i,j,k)
             dqra(i,j,k)=dqra(i,j,k)+para1sg*para2sg*gonv**(-0.75)*(rhoair*wgr)**(1.75)*pxkh_tlr*zmm_ref(i,j,k)   &
			                        +1.75*para1sg*para2sg*gonv**(-0.75)*(rhoair*wgr)**(0.75)*pxkh*pdfrrg*qgr(i,j,k)*zmm_ref(i,j,k) &
									+1.75*para1sg*para2sg*gonv**(-0.75)*(rhoair*wgr)**(0.75)*pxkh*rhoair*(1-dgr_coef)*zmm_ref(i,j,k)   
             dqgr(i,j,k)=dqgr(i,j,k)+para1sg*para2sg*gonv**(-0.75)*(rhoair*wgr)**(1.75)*pxkh_tlx*zmm_ref(i,j,k)
			 dtmk(i,j,k)=dtmk(i,j,k)+1.75*para1sg*para2sg*gonv**(-0.75)*(rhoair*wgr)**(0.75)*pxkh*wgr*pdfrhot*zmm_ref(i,j,k)
			 dqvp(i,j,k)=dqvp(i,j,k)+1.75*para1sg*para2sg*gonv**(-0.75)*(rhoair*wgr)**(0.75)*pxkh*wgr*pdfrhoq*zmm_ref(i,j,k)
           endif
         endif
         if(print_chk>=2) then
           print*,"zgh,10.*log10(zgh)"
           print*,zgh,max(0.,10.*log10(zgh)),dqgr(i,j,k),dgr_coef
         endif

         if(dualpol_opt==1) then
           if(in0g.eq.1) then
             ! for two moment microphysics scheme, not yet completed  
           else
             if(tlopt==0) zgv= para1sg*para2sg*sonv**(-0.75)*(rhoair*wgr)**(1.75)*pxkv
             if(tlopt>=1.and.gropt==0) zgv= 1.75*para1sg*para2sg*gonv**(-0.75)*(rhoair*wgr)**(0.75)*pxkv*rhoair*(1-dgr_coef)*dqgr(i,j,k)       &
		                      +1.75*para1sg*para2sg*gonv**(-0.75)*(rhoair*wgr)**(0.75)*pxkv*rhoair*(1-dgr_coef)*dqra(i,j,k)                    &
							  +1.75*para1sg*para2sg*gonv**(-0.75)*(rhoair*wgr)**(0.75)*pxkh*wgr*pdfrhot*dtmk(i,j,k)                            &
							  +1.75*para1sg*para2sg*gonv**(-0.75)*(rhoair*wgr)**(0.75)*pxkh*wgr*pdfrhoq*dqvp(i,j,k)                            &
		                      +1.75*para1sg*para2sg*gonv**(-0.75)*(rhoair*wgr)**(0.75)*pxkv*pdfrrg*qgr(i,j,k)*dqra(i,j,k)                      & 
							  +1.75*para1sg*para2sg*gonv**(-0.75)*(rhoair*wgr)**(0.75)*pxkv*pdfgrg*qgr(i,j,k)*dqgr(i,j,k)                      & 
                              +para1sg*para2sg*gonv**(-0.75)*(rhoair*wgr)**(1.75)*pxkv_tlr*dqra(i,j,k)                                         &
                              +para1sg*para2sg*gonv**(-0.75)*(rhoair*wgr)**(1.75)*pxkv_tlx*dqgr(i,j,k)
             if(tlopt>=1.and.gropt>=1) then
               ! here, zmm_ref serves as the HX (or yb) for gradient calculation(not available yet)
               dqgr(i,j,k)=dqgr(i,j,k)+1.75*para1sg*para2sg*gonv**(-0.75)*(rhoair*wgr)**(0.75)*pxkv*rhoair*(1-dgr_coef)*zmm_ref(i,j,k) &
			                          +1.75*para1sg*para2sg*gonv**(-0.75)*(rhoair*wgr)**(0.75)*pxkv*pdfgrg*qgr(i,j,k)*zmm_ref(i,j,k)
               dqra(i,j,k)=dqra(i,j,k)+para1sg*para2sg*gonv**(-0.75)*(rhoair*wgr)**(1.75)*pxkv_tlr*zmm_ref(i,j,k)  &
			                          +1.75*para1sg*para2sg*gonv**(-0.75)*(rhoair*wgr)**(0.75)*pxkv*rhoair*(1-dgr_coef)*zmm_ref(i,j,k) &
			                          +1.75*para1sg*para2sg*gonv**(-0.75)*(rhoair*wgr)**(0.75)*pxkv*pdfrrg*qgr(i,j,k)*zmm_ref(i,j,k)
               dqgr(i,j,k)=dqgr(i,j,k)+para1sg*para2sg*gonv**(-0.75)*(rhoair*wgr)**(1.75)*pxkv_tlx*zmm_ref(i,j,k)
			   dtmk(i,j,k)=dtmk(i,j,k)+1.75*para1sg*para2sg*gonv**(-0.75)*(rhoair*wgr)**(0.75)*pxkh*wgr*pdfrhot*zmm_ref(i,j,k)
			   dqvp(i,j,k)=dqvp(i,j,k)+1.75*para1sg*para2sg*gonv**(-0.75)*(rhoair*wgr)**(0.75)*pxkh*wgr*pdfrhoq*zmm_ref(i,j,k)			 
             endif 
           endif
         endif
         if(print_chk>=2) then
           print*,"zgv,10.*log10(zgv)"
           print*,zgv,max(0.,10.*log10(zgv)),dqgr(i,j,k),dgr_coef
         endif
!         print*,"----------------for dry graupel----------------------"
         call parameter_zxx(para1sg,para2sg,1.0,rhog,                     &
                            cice,dg,alphag,mm3todBZ,lambda,Kw2,pi,        &
                            gonv)

         ice_abc_d(1)=Ahd
         ice_abc_d(2)=Bhd
         ice_abc_d(3)=Chd
         ice_bac_d(1)=Bhd
         ice_bac_d(2)=Ahd						
						
         pxabk_all(1)=alpha_rdha**2
         pxabk_all(2)=alpha_rdhb**2
         pxabk_all(3)=alpha_rdhb*alpha_rdha
         pxkh=pkx(ice_abc_d,pxabk_all)
         pxkv=pkx(ice_bac_d,pxabk_all)
       !  print*,"pxkh,pxkv"
       !  print*,pxkh,pxkv


         if(in0s.eq.1) then
           ! for two moment microphysics scheme, not yet completed
         else
           if(tlopt==0) zdgh=para1sg*para2sg*gonv**(-0.75)*(rhoair*dgr)**(1.75)*pxkh
           if(tlopt>=1.and.gropt==0) zdgh=1.75*para1sg*para2sg*gonv**(-0.75)*(rhoair*dgr)**(0.75)*pxkh*rhoair*dgr_coef*dqgr(i,j,k)           &
		                                 +1.75*para1sg*para2sg*gonv**(-0.75)*(rhoair*dgr)**(0.75)*pxkh*dgr*pdfrhot*dtmk(i,j,k)               &
										 +1.75*para1sg*para2sg*gonv**(-0.75)*(rhoair*dgr)**(0.75)*pxkh*dgr*pdfrhoq*dqvp(i,j,k)               &
                                         -1.75*para1sg*para2sg*gonv**(-0.75)*(rhoair*dgr)**(0.75)*pxkh*pdfrrg*qgr(i,j,k)*dqra(i,j,k)         &
                                         -1.75*para1sg*para2sg*gonv**(-0.75)*(rhoair*dgr)**(0.75)*pxkh*pdfgrg*qgr(i,j,k)*dqgr(i,j,k)
           if(tlopt>=1.and.gropt>=1) then
             ! here, zmm_ref serves as the HX (or yb) for gradient calculation
             dqgr(i,j,k)=dqgr(i,j,k)+1.75*para1sg*para2sg*gonv**(-0.75)*(rhoair*dgr)**(0.75)*pxkh*rhoair*dgr_coef*zmm_ref(i,j,k) &
			                        -1.75*para1sg*para2sg*gonv**(-0.75)*(rhoair*dgr)**(0.75)*pxkh*pdfgrg*qgr(i,j,k)*zmm_ref(i,j,k)
			 dqra(i,j,k)=dqra(i,j,k)-1.75*para1sg*para2sg*gonv**(-0.75)*(rhoair*dgr)**(0.75)*pxkh*pdfrrg*qgr(i,j,k)*zmm_ref(i,j,k)	
             dtmk(i,j,k)=dtmk(i,j,k)+1.75*para1sg*para2sg*gonv**(-0.75)*(rhoair*dgr)**(0.75)*pxkh*dgr*pdfrhot*zmm_ref(i,j,k)
             dqvp(i,j,k)=dqvp(i,j,k)+1.75*para1sg*para2sg*gonv**(-0.75)*(rhoair*dgr)**(0.75)*pxkh*dgr*pdfrhoq*zmm_ref(i,j,k)             			 
           endif           
         endif
         if(print_chk>=2) then
           print*,"zdgh,10.*log10(zdsh)"
           print*,zdgh,max(0.,10.*log10(zdgh)),dqgr(i,j,k),dgr_coef
         endif

         if(dualpol_opt==1) then
           if(in0s.eq.1) then
             ! for two moment microphysics scheme, not yet completed
           else
             if(tlopt==0) zdgv=para1sg*para2sg*gonv**(-0.75)*(rhoair*dgr)**(1.75)*pxkv
             if(tlopt>=1.and.gropt==0) zdgv=1.75*para1sg*para2sg*gonv**(-0.75)*(rhoair*dgr)**(0.75)*pxkv*rhoair*dgr_coef*dqgr(i,j,k)           &
		                                   +1.75*para1sg*para2sg*gonv**(-0.75)*(rhoair*dgr)**(0.75)*pxkh*dgr*pdfrhot*dtmk(i,j,k)               &
										   +1.75*para1sg*para2sg*gonv**(-0.75)*(rhoair*dgr)**(0.75)*pxkh*dgr*pdfrhoq*dqvp(i,j,k)               &		   
                                           -1.75*para1sg*para2sg*gonv**(-0.75)*(rhoair*dgr)**(0.75)*pxkv*pdfrrg*qgr(i,j,k)*dqra(i,j,k)         &
                                           -1.75*para1sg*para2sg*gonv**(-0.75)*(rhoair*dgr)**(0.75)*pxkv*pdfgrg*qgr(i,j,k)*dqgr(i,j,k)  
             if(tlopt>=1.and.gropt>=1) then
               ! here, zmm_ref serves as the HX (or yb) for gradient calculation(not available yet)
               dqgr(i,j,k)=dqgr(i,j,k)+1.75*para1sg*para2sg*gonv**(-0.75)*(rhoair*dgr)**(0.75)*pxkv*rhoair*dgr_coef*zmm_ref(i,j,k) &
			                          -1.75*para1sg*para2sg*gonv**(-0.75)*(rhoair*dgr)**(0.75)*pxkv*pdfgrg*qgr(i,j,k)*zmm_ref(i,j,k)
			   dqra(i,j,k)=dqra(i,j,k)-1.75*para1sg*para2sg*gonv**(-0.75)*(rhoair*dgr)**(0.75)*pxkv*pdfrrg*qgr(i,j,k)*zmm_ref(i,j,k)
               dtmk(i,j,k)=dtmk(i,j,k)+1.75*para1sg*para2sg*gonv**(-0.75)*(rhoair*dgr)**(0.75)*pxkh*dgr*pdfrhot*zmm_ref(i,j,k)
               dqvp(i,j,k)=dqvp(i,j,k)+1.75*para1sg*para2sg*gonv**(-0.75)*(rhoair*dgr)**(0.75)*pxkh*dgr*pdfrhoq*zmm_ref(i,j,k)    			 
             endif
           endif
         endif
         if(print_chk>=2) then
           print*,"zdgv,10.*log10(zdgv)"
           print*,zdgv,max(0.,10.*log10(zdgv)),dqgr(i,j,k),dgr_coef
         endif
         zhs(i,j,k)=zgh+zdgh
!    ==================done=================================
       endif  !(if rf_noice == 0)

         z_e =zrh+zsh+zgh+zdsh+zdgh
         zv_e=zrv+zsv+zgv+zdsv+zdgv
!
!      Convert to dBZ
!
         if(tlopt==0) then
           dbz(i,j,k) = max(0.,10. * log10(z_e))
         else
           if(zmm(i,j,k)>1.0) then
             dbz(i,j,k) = 10./(zmm(i,j,k)*log(10.0))*z_e !ref(i,j,k)
           else
             dbz(i,j,k) =0.0
           endif
		   !print*,z_e,zmm(i,j,k)
         endif

         if(abs(dbz(i,j,k))>5.0.and.print_chk>=1) then
           print*,"------------------total dbz---------------------"
           print*,z_e,dbz(i,j,k),zmm(i,j,k),zmm_ref(i,j,k)
         endif
!        save z_e mm^6 mm^-3
         zmm(i,j,k)=z_e
!        stop
!
      enddo
    enddo
  enddo
!
  contains

  subroutine parameter_zrx(p10,p14,rhoa,rhor,cr,dr,alphar,betarx,    &
                           alpharx,mm3todBZ,lambda,Kw2,pi,n0r)

  real :: p1,p2,p3,p4,p5,p7,p8,p9,p10,p11,p12,p13,p14
  real :: rhoa,rhor,cr,dr,alphar,betarx,alpharx
  real :: mm3todBZ,lambda,Kw2,pi,n0r
      
  p1=mm3todBZ*(4*lambda**4*alpharx**2/(pi**4*Kw2))
  p2=-(2.*betarx+1.0)
  p7=(pi*rhor/rhoa)**(p2/4.)
  p10=(1.+p2/4.)
  p14=p1*p7*(n0r)**p10*gamma(-p2*1d0)

  end subroutine parameter_zrx

  function lower_lambdax(rhox,dx,alphax,cx,qnx,qx)

  real :: rhox,dx,alphax,cx,qnx,qx
  real :: gamma
  real :: lower_lambdax
      
  lower_lambdax=( (gamma((1.0+dx+alphax))*1d0*cx*qnx)   &
                 /(gamma((1.0+alphax)*1d0)*rhox*qx)    &
                  )**(1.0/dx)

  end function lower_lambdax
      
  function n0x(qnx,lower_lambda,alphax)
      
  real :: qnx,lower_lambda,alphax
  real :: gamma
  real :: n0x
      
  n0x=qnx*lower_lambda**(1.0+alphax)/gamma(1.0+alphax)
      
  end function n0x

  function upper_f(qr,qice,qthres,flg)

  real :: qr, qice,qthres
  real :: upper_f
  real :: fmax=0.5
  integer :: flg

  if(flg==1) fmax=0.5 ! for snow
  if(flg==2) fmax=0.3 ! for hail/graupel
      
  upper_f=0
  if(qr>=qthres.and.qice>=qthres) then
    upper_f=fmax*min(qice/qr,qr/qice)**0.3
  endif

  end function upper_f

  function waterfraction(qr,qice)

  real :: qr,qice
  real :: waterfraction
 
  if(qr<1.e-8.and.qice<1.e-8) then
    waterfraction=1.e-8
  else
    waterfraction=qr/(qr+qice) 
  endif

  end function waterfraction


  subroutine calc_ice_abc(phimean,sigma,ice_abc)

  real :: phimean,sigma
  real :: ice_abc(3)  ! 1 for A, 2 for B, 3 for C
  real,parameter :: pi=3.1415926
  real :: a2a=pi/180.

  ice_abc(1)=1./8.*(3.+4.*cos(2.*phimean*pi)*exp(-2.*sigma**2)  &
              +cos(4.*phimean)*exp(-8.*sigma**2))
  ice_abc(2)=1./8.*(3.-4.*cos(2.*phimean*pi)*exp(-2.*sigma**2) &
              +cos(4.*phimean)*exp(-8.*sigma**2))
  ice_abc(3)=1./8.*(1.                                          &
              -cos(4.*phimean)*exp(-8.*sigma**2)) 

  end subroutine calc_ice_abc

  function sigma_in_abc(qg,fw,gsflag)

  real :: qg,fw
  integer :: gsflag
  real :: sigma_in_abc
  real,parameter :: pi=3.1415926
  real :: a2a=pi/180.

  if(gsflag==1) then ! for snow
    sigma_in_abc=20.*a2a
  elseif(gsflag==2) then ! for graupel
    if(qg>0.2*1e-3) then 
      sigma_in_abc=60.*a2a*(1-0.8*fw)
    else
      sigma_in_abc=60.*a2a*(1-4.*qg*1e3*fw)
    endif
  endif

  end function sigma_in_abc

  function alpha_rxx(fw,alpha_rxab,n)

  integer :: n
  real :: fw
  real :: alpha_rxab(n)
  real :: alpha_rxx
  integer :: i,j,k
      
  alpha_rxx=0
  do i=0,n-1
    alpha_rxx=1.e-3*alpha_rxx+alpha_rxab(i+1)*fw**i
  enddo

  end function alpha_rxx

  function pxabk(alpha_rxa,alpha_rxb,k,n)

  integer :: n
  integer :: k
  real :: alpha_rxa(n),alpha_rxb(n)
  real :: pxabk
  integer :: i,j

  pxabk=0
  do j=0,n-1
    do i=0,n-1
      if(i+j==k) then
        pxabk=pxabk+alpha_rxa(i+1)*alpha_rxb(j+1)
      endif
    enddo
  enddo

  end function pxabk

  function pkx(ice_abc,pxabk_all)

  real :: ice_abc(3),pxabk_all(3)
  real :: pkx
  integer :: i

  pkx=0
  do i=1,3
    if(i<3) then
      pkx=pkx+ice_abc(i)*pxabk_all(i)
    else
      pkx=pkx+2.*ice_abc(i)*pxabk_all(i)
    endif
  enddo

  end function pkx

  subroutine parameter_zxx(p1,p2,rhoa,rhox,cx,dx,alphax,      &
                           mm3todBZ,lambda,Kw2,pi,n0x)  

  real :: p1,p2
  real :: rhoa,rhox,cx,dx,alphax
  real :: mm3todBZ,lambda,Kw2,pi,n0x
  real :: gamma

  p1=mm3todBZ*4.*gamma(7.)*lambda**4/(pi**4*Kw2)
  p2=(pi*rhox/rhoa)**(-1.75)

  end subroutine parameter_zxx

  function sum_pxkfw(pxk_all,fw,n)

  integer :: n
  real :: pxk_all(2*n-1)
  real :: fw
  integer :: k
  real :: sum_pxkfw

  sum_pxkfw=0
  do k=0,2*(n-1)
     sum_pxkfw=sum_pxkfw+pxk_all(k+1)*fw**k
  enddo 

  end function sum_pxkfw

  function virtual(temp,ratmix)
!
!   This function returns virtual temperature in K, given temperature
!      in K and mixing ratio in kg/kg.
!
  real,parameter :: eps=0.622
  real :: virtual
  real :: temp
  real :: ratmix
!
  !virtual=temp*(eps+ratmix)/(eps*(1.+ratmix))
  virtual=temp*(1+0.61*ratmix)
  return
  end function virtual

  FUNCTION gamma(xx)

!  Modified from "Numerical Recipes"

  IMPLICIT NONE

! PASSING PARAMETERS:
  DOUBLE PRECISION, INTENT(IN) :: xx

! LOCAL PARAMETERS:
  DOUBLE PRECISION  :: gamma
  INTEGER  :: j
  DOUBLE PRECISION  :: ser,stp,tmp,x,y,cof(6)


  SAVE cof,stp
  DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,               &
       24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,  &
       -.5395239384953d-5,2.5066282746310005d0/
  x=xx
  y=x
  tmp=x+5.5d0
  tmp=(x+0.5d0)*log(tmp)-tmp
  ser=1.000000000190015d0
! do j=1,6   !original
  do j=1,4
!!do j=1,3   !gives result to within ~ 3 %
     y=y+1.d0
     ser=ser+cof(j)/y
  enddo
  gamma=tmp+log(stp*ser/x)
  gamma= exp(gamma)

  END FUNCTION gamma

end subroutine RadZIceVar
