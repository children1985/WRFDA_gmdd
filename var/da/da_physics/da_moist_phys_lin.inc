subroutine da_moist_phys_lin(grid)

   !---------------------------------------------------------------------------
   !  Purpose: Partition of the hydrometeors via the moist explicit scheme.
   !           A warm rain process is used in this subroutine.
   !           This is the tangent linear code of the scheme.
   !
   !  Method: The warm rain process is according to Hsie and Anthes (1984)
   !          and Dudhia (1989)
   !
   !  Assumptions: 1) Model level stored top down.
   !---------------------------------------------------------------------------

   implicit none

   type(domain), intent(inout)               :: grid

   real :: T_OLD(ims:ime,jms:jme,kms:kme),T_NEW(ims:ime,jms:jme,kms:kme)
   real :: Q_OLD(ims:ime,jms:jme,kms:kme),Q_NEW(ims:ime,jms:jme,kms:kme)
   real :: QCW_OLD(ims:ime,jms:jme,kms:kme),QCW_NEW(ims:ime,jms:jme,kms:kme)
   real :: QRN_OLD(ims:ime,jms:jme,kms:kme),QRN_NEW(ims:ime,jms:jme,kms:kme)

   real    :: EES(kms:kme)
   real    :: QVSS(kms:kme)
   real    :: EES9(kms:kme)
   real    :: QVSS9(kms:kme)
   real    :: DT(its:ite,jts:jte,kms:kme)
   real    :: QVT(kms:kme)
   real    :: QCT(kms:kme)
   real    :: QRT(kms:kme)
   real    :: TTT(kms:kme)
   real    :: QVT9(kms:kme)
   real    :: QCT9(kms:kme)
   real    :: QRT9(kms:kme)
   real    :: TTT9(kms:kme)
   real    :: SCR2(kms:kme)
   real    :: SCR3(kms:kme)
   real    :: SCR4(kms:kme)
   real    :: SCR5(kms:kme)
   real    :: SCR6(kms:kme)
   real    :: DUM31(kms:kme)
   real    :: PRA(kms:kme)
   real    :: PRC(kms:kme)
   real    :: PRD(kms:kme)
   real    :: PRE(kms:kme)
   real    :: SCR31(kms:kme)
   real    :: SCR42(kms:kme)
   real    :: SCR71(kms:kme)
   real    :: DUM112(kms:kme)
   real    :: DUM113(kms:kme)
   real    :: DUM211(kms:kme)
   real    :: DUM411(kms:kme)
   real    :: SCR29(kms:kme)
   real    :: SCR39(kms:kme)
   real    :: SCR49(kms:kme)
   real    :: SCR59(kms:kme)
   real    :: SCR69(kms:kme)
   real    :: DUM319(kms:kme)
   real    :: PRA9(kms:kme)
   real    :: PRC9(kms:kme)
   real    :: PRD9(kms:kme)
   real    :: PRE9(kms:kme)
   real    :: SCR319(kms:kme)
   real    :: SCR429(kms:kme)
   real    :: SCR719(kms:kme)
   real    :: DUM1129(kms:kme)
   real    :: DUM4119(kms:kme)
   real    :: TMP(kms:kme)

   integer, parameter :: qcth = 0.5e-3
   integer, parameter :: qrth = 1.0e-6
   integer, parameter :: alpha = 0.001
   integer, parameter :: beta  = 0.0486
   integer, parameter :: gamma = 0.002

   integer :: i, j, k

   real    :: qrth1

   if (trace_use) call da_trace_entry("da_moist_phys_lin")

   qrth1 = (QRTH*1.0e3)**0.875

   T_OLD  (its:ite,jts:jte,kts:kte) = grid%xa%t  (its:ite,jts:jte,kts:kte)
   Q_OLD  (its:ite,jts:jte,kts:kte) = grid%xa%q  (its:ite,jts:jte,kts:kte)
   QCW_OLD(its:ite,jts:jte,kts:kte) = grid%xa%qcw(its:ite,jts:jte,kts:kte)
   QRN_OLD(its:ite,jts:jte,kts:kte) = grid%xa%qrn(its:ite,jts:jte,kts:kte)

   !  Preparation

   grid%xa%q(its:ite,jts:jte,kts:kte) =grid%xa%qt(its:ite,jts:jte,kts:kte) - grid%xa%qcw(its:ite,jts:jte,kts:kte) - grid%xa %qrn(its:ite,jts:jte,kts:kte)
   DT(its:ite,jts:jte,kts:kte) = grid%xb%delt(its:ite,jts:jte,kts:kte)

   do j=jts,jte
      do i=its,ite
         do K=kts,kte

            if (dt(i,j,k) <= 0.0) cycle

            if ( grid%xb%t(I,J,K) > TO )then
               EES(K)=SVP1*EXP(SVP2*(grid%xb%t(I,J,K)-SVPT0)/(grid%xb%t(I,J,K)-SVP3))
               EES9(K)=EES(K)*SVP2*(SVPT0-SVP3)/((grid%xb%t(I,J,K)-SVP3) * (grid%xb%t(I,J,K)-SVP3))*grid%xa%t(I,J,K)
            else
               EES(K)=.611*EXP(22.514-6.15E3/grid%xb%t(I,J,K))
               EES9(K)=EES(K)*6.15E3/(grid%xb%t(I,J,K)*grid%xb%t(I,J,K))*grid%xa%t(I,J,K)
            end if

            TMP(K)=622.0/((grid%xb%p(I,J,K)-EES(K))**2)
            QVSS9(K)=TMP(K)*grid%xb%p(I,J,K)*EES9(K) - TMP(K)*EES(K)*grid%xa%p(I,J,K)
            QVSS(K)=622.0*EES(K)/(grid%xb%p(I,J,K)-EES(K))

            SCR49(K)=grid%xa%q(I,J,K)/QVSS(K)-grid%xb%q(I,J,K)/QVSS(K)**2*QVSS9(K)
            SCR4(K)=grid%xb%q(I,J,K)/QVSS(K)

            if (grid%xb%qcw(I,J,K) > 0.0) then
               SCR29(K)=grid%xa%qcw(I,J,K)
               SCR2(K)=grid%xb%qcw(I,J,K)
            else
               SCR29(K)=0.0
               SCR2(K)=0.0
            end if
            if (grid%xb%qrn(I,J,K) > 1.0e-25) then
               SCR39(K)=grid%xa%qrn(I,J,K)
               SCR3(K)=grid%xb%qrn(I,J,K)
            else
               SCR39(K)=0.0
               SCR3(K)=1.0E-25
            end if
            SCR59(K)=grid%xa%q(I,J,K)/SCR4(K)-grid%xb%q(I,J,K)/SCR4(K)**2*SCR49(K)
            SCR5(K)=grid%xb%q(I,J,K)/SCR4(K)

            SCR69(K)=grid%xa%p(I,J,K)/(gas_constant*grid%xb%t(I,J,K))-grid%xb%p(I,J,K)/  &
                     (gas_constant*grid%xb%t(I,J,K)**2)*grid%xa%t(I,J,K)
            SCR6(K)=grid%xb%p(I,J,K)/(gas_constant*grid%xb%t(I,J,K))

            DUM319(K)=-XLV1*grid%xa%t(I,J,K) 
            DUM31(K)=3.1484E6-XLV1*grid%xb%t(I,J,K)
 
            ! Auto conversion

            if (scr2(k) >= qcth) then
               prc9(k) = alpha * scr29(k)
               prc(k) = alpha * (scr2(k) - qcth)
            else
               prc9(k) = 0.0
               prc(k) = 0.0
            end if

            ! Accretion

            if (SCR2(k) > 0.0 .and. SCR3(k) > QRTH ) then
               PRA9(K) = gamma * 0.875 * SCR2(k) * (SCR3(K)*1.0e3)**(-0.125) * 1.0e3 * SCR39(K)  &
                            + gamma * SCR29(k) * (SCR3(K)*1.0e3)**0.875
               PRA(k) = gamma * SCR2(k) * (SCR3(k)*1.0e3)**0.875
            else if (SCR2(k) > 0.0 .and. SCR3(k) <= QRTH ) then
               PRA9(K) = gamma * SCR29(k) * qrth1
               PRA(k) = gamma * SCR2(k) * qrth1
            else
               PRA9(K) = 0.0
               PRA(k) = 0.0
            end if

         end do

         call da_evapo_lin(DT(i,j,:),SCR3,SCR5,grid%xb%q(I,J,:),PRE,SCR6,  &
                           SCR39,SCR59,grid%xa%q(I,J,:),PRE9,SCR69, &
                           kts,kte,kms,kme)

         do K=kts, kte

            if (dt(i,j,k) <= 0.0) cycle

            !  Readjust

            DUM112(K)=(PRC(k)+PRA(k))*dt(i,j,k)
            if (DUM112(K) > SCR2(k)) then
               DUM1129(K)=(PRC9(k)+PRA9(k))*dt(i,j,k)
               PRC9(K)=SCR29(K)*PRC(K)/DUM112(K)  &
                      +PRC9(K)*SCR2(K)/DUM112(K)  &
                      -SCR2(K)*PRC(K)/(DUM112(K)*DUM112(K))*DUM1129(K)
               PRC(k)=SCR2(K)*PRC(K)/DUM112(K)
               PRA9(K)=SCR29(K)*PRA(K)/DUM112(K)  &
                      +PRA9(K)*SCR2(K)/DUM112(K)  &
                      -SCR2(K)*PRA(K)/(DUM112(K)*DUM112(K))*DUM1129(K)
               PRA(k)=SCR2(K)*PRA(K)/DUM112(K)
            end if
            QVT9(K)=-PRE9(K)
            QVT(K)=-PRE(K)
            QCT9(K)=-PRC9(K)-PRA9(K)
            QCT(K)=-PRC(K)-PRA(K)
            QRT9(K)=PRC9(K)+PRA9(K)+PRE9(K)
            QRT(K)=PRC(K)+PRA(K)+PRE(K)
            if (grid%xb%t(I,J,K).GT.TO)then
               DUM4119(K)=DUM319(K)
               DUM411(K)=DUM31(K)
            else
               DUM4119(K)=0.0
               DUM411(K)=XLS
            end if
            PRD9(K)=cp*0.887*grid%xa%q(I,J,K)
            PRD(K)=cp*(1.0+0.887*grid%xb%q(I,J,K))
            TTT9(K)=-DUM4119(K)*QVT(K)/PRD(K)  &
                   -QVT9(K)*DUM411(K)/PRD(K)  &
                   +DUM411(K)*QVT(K)/(PRD(K)*PRD(K))*PRD9(K)
            TTT(K)=-DUM411(K)*QVT(K)/PRD(K)

            DUM113(K)=grid%xb%q(I,J,K)+dt(i,j,k)*QVT(K)
            if (DUM113(K) > 1.0e-12 ) then
               SCR429(K)=grid%xa%q(I,J,K)+dt(i,j,k)*QVT9(K)
               SCR42(K)=DUM113(K)
            else
               SCR429(K)=0.0
               SCR42(K)=1.0e-12
            end if
            DUM211(K)=grid%xb%qcw(I,J,K)+QCT(K)*dt(i,j,k)
            if (DUM211(K) > 0.0) then
               SCR319(K)=grid%xa%qcw(I,J,K)+QCT9(K)*dt(i,j,k)
               SCR31(K)=DUM211(K)
            else
               SCR319(K)=0.0
               SCR31(K)=0.0
            end if
            SCR719(K)=grid%xa%t(I,J,K)+TTT9(K)*dt(i,j,k)
            SCR71(K)=grid%xb%t(I,J,K)+TTT(K)*dt(i,j,k)
         end do

         call da_condens_lin(DT(i,j,:),SCR31,SCR42,SCR71,DUM31,PRD,         &
                             QVT,QCT,QRT,TTT,                        &
                             grid%xb%p(I,J,:),grid%xb%t(I,J,:),grid%xb%q(I,J,:),    &
                             grid%xb%qcw(I,J,:),grid%xb%qrn(I,J,:),            &
                             SCR319,SCR429,SCR719,DUM319,PRD9,       &
                             QVT9,QCT9,QRT9,TTT9,                    &
                             grid%xa%p(I,J,:),grid%xa%t(I,J,:),grid%xa%q(I,J,:),    &
                             grid%xa%qcw(I,J,:),grid%xa%qrn(I,J,:),kts,kte)
      end do
   end do

   T_NEW  (its:ite,jts:jte,kds:kde) = grid%xa%t   (its:ite,jts:jte,kds:kde) - T_OLD  (its:ite,jts:jte,kds:kde)
   Q_NEW  (its:ite,jts:jte,kds:kde) = grid%xa%q   (its:ite,jts:jte,kds:kde) - Q_OLD  (its:ite,jts:jte,kds:kde)
   QCW_NEW(its:ite,jts:jte,kds:kde) = grid%xa%qcw (its:ite,jts:jte,kds:kde) - QCW_OLD(its:ite,jts:jte,kds:kde)
   QRN_NEW(its:ite,jts:jte,kds:kde) = grid%xa%qrn (its:ite,jts:jte,kds:kde) - QRN_OLD(its:ite,jts:jte,kds:kde)

   call da_filter(grid, t_new)
   call da_filter(grid, q_new)
   call da_filter(grid, qcw_new)
   call da_filter(grid, qrn_new)

   grid%xa%t   (its:ite,jts:jte,kds:kde) = T_NEW  (its:ite,jts:jte,kds:kde) + T_OLD  (its:ite,jts:jte,kds:kde)
   grid%xa%q   (its:ite,jts:jte,kds:kde) = Q_NEW  (its:ite,jts:jte,kds:kde) + Q_OLD  (its:ite,jts:jte,kds:kde)
   grid%xa%qcw (its:ite,jts:jte,kds:kde) = QCW_NEW(its:ite,jts:jte,kds:kde) + QCW_OLD(its:ite,jts:jte,kds:kde)
   grid%xa%qrn (its:ite,jts:jte,kds:kde) = QRN_NEW(its:ite,jts:jte,kds:kde) + QRN_OLD(its:ite,jts:jte,kds:kde)

#ifdef DM_PARALLEL
#include "HALO_XA_CLOUD.inc"
#endif

   if (trace_use) call da_trace_exit("da_moist_phys_lin")

end subroutine da_moist_phys_lin


subroutine thermoist_trans(grid,gropt)

  implicit none

  type(domain), intent(inout)      :: grid

  integer                          :: nvar,nz
  real                             :: bakvar(kts:kte,7)
  real                             :: prtvar(kts:kte,7)
  real                             :: grdvar(kts:kte,7)
  real                             :: kalman(kts:kte,kts:kte,3)
  real                             :: kvar(kts:kte,3)
  integer                          :: i,j,k,k1,k2,ivar
  integer                          :: gropt
  
  character(len=80)                :: filepath
  logical                          :: alive
  
  nvar=7
  nz=kte-kts+1
  filepath='be_dvar_ascii.dat'
  !print*,kts,kte,nz
  inquire(file=trim(filepath),exist=alive) 
  if(.not.alive) return
 
  !open(1006,file=trim(filepath),form='binary')
  !read(1006) kalman
  !close(1006)

  open(1006,file=trim(filepath))
  read(1006,*) kalman
  close(1006)
  
  if(gropt==0) then

    do j=jts,jte ! 1,ny
      do i=its,ite  ! 1,nx
	  
        prtvar(:,1)=grid%xa%t  (i,j,kts:kte)	
		prtvar(:,2)=grid%xa%q  (i,j,kts:kte)	
		prtvar(:,3)=grid%xa%qcw(i,j,kts:kte)	
		prtvar(:,4)=grid%xa%qci(i,j,kts:kte)	
		prtvar(:,5)=grid%xa%qrn(i,j,kts:kte)	
		prtvar(:,6)=grid%xa%qsn(i,j,kts:kte)	
		prtvar(:,7)=grid%xa%qgr(i,j,kts:kte)	
		
        bakvar(:,1)=grid%xb%t  (i,j,kts:kte)	
		bakvar(:,2)=grid%xb%q  (i,j,kts:kte)	
		bakvar(:,3)=grid%xb%qcw(i,j,kts:kte)	
		bakvar(:,4)=grid%xb%qci(i,j,kts:kte)	
		bakvar(:,5)=grid%xb%qrn(i,j,kts:kte)	
		bakvar(:,6)=grid%xb%qsn(i,j,kts:kte)	
		bakvar(:,7)=grid%xb%qgr(i,j,kts:kte)		

	    kvar=0  
	    call transform(bakvar,kalman,gropt,nz,nvar,prtvar,kvar,grdvar)	

        grid%xa%qrn(i,j,kts:kte)=kvar(:,1)
		grid%xa%qsn(i,j,kts:kte)=kvar(:,2)
		grid%xa%qgr(i,j,kts:kte)=kvar(:,3)
		
	  enddo
    enddo 
	
  else
	
    do j=jts,jte ! 1,ny
      do i=its,ite  ! 1,nx
	  
        prtvar=0	
		prtvar(:,5)=grid%xa%qrn(i,j,kts:kte)	
		prtvar(:,6)=grid%xa%qsn(i,j,kts:kte)	
		prtvar(:,7)=grid%xa%qgr(i,j,kts:kte)

		bakvar(:,1)=grid%xb%t  (i,j,kts:kte)	
		bakvar(:,2)=grid%xb%q  (i,j,kts:kte)	
		bakvar(:,3)=grid%xb%qcw(i,j,kts:kte)	
		bakvar(:,4)=grid%xb%qci(i,j,kts:kte)	
		bakvar(:,5)=grid%xb%qrn(i,j,kts:kte)	
		bakvar(:,6)=grid%xb%qsn(i,j,kts:kte)	
		bakvar(:,7)=grid%xb%qgr(i,j,kts:kte)
		grdvar=0	
		
	    call transform(bakvar,kalman,gropt,nz,nvar,prtvar,kvar,grdvar)

		grid%xa%t  (i,j,kts:kte)=grdvar(:,1)+grid%xa%t  (i,j,kts:kte)	
		grid%xa%q  (i,j,kts:kte)=grdvar(:,2)+grid%xa%q  (i,j,kts:kte)
		grid%xa%qcw(i,j,kts:kte)=grdvar(:,3)+grid%xa%qcw(i,j,kts:kte)
		grid%xa%qci(i,j,kts:kte)=grdvar(:,4)+grid%xa%qci(i,j,kts:kte)
        grid%xa%qrn(i,j,kts:kte)=grdvar(:,5)
		grid%xa%qsn(i,j,kts:kte)=grdvar(:,6)
		grid%xa%qgr(i,j,kts:kte)=grdvar(:,7)
		
	  enddo
    enddo     
  endif

contains
  
subroutine transform(bakvar,kalman,gropt,nz,nvar,prtvar,kvar,grdvar)

implicit none

integer        :: nz,nvar
integer        :: gropt

real           :: bakvar(nz,nvar)
real           :: prtvar(nz,nvar)
real           :: grdvar(nz,nvar)
real           :: kvar  (nz,3)
real           :: paravp(nz,nvar)
real           :: kalman(nz,nz,3)

real           :: ratio_vp
integer        :: i,j,k,k1,k2

  ratio_vp=0.5

  do k=1,nz
    paravp(k,1)=2.e-1/(bakvar(k,1)+300.0)
    paravp(k,2)=0.61
    paravp(k,3:7)=1.0	
  enddo

  call transform_phy(prtvar,kvar,grdvar,kalman,paravp,nvar,3,nz,gropt,ratio_vp)	



end subroutine transform

subroutine transform_phy(prtvar,kvar,grdvar,kalman,para,nvar,nvar2,nz,gropt,ratio_vp)

implicit none

integer                     :: gropt
integer                     :: nvar,nz
integer                     :: nvar2
real                        :: ratio_vp

real,dimension(nz,nvar)     :: prtvar,grdvar
real,dimension(nz,nvar2)    :: kvar
real,dimension(nz,nvar)     :: para
real,dimension(nz,nz,nvar2) :: kalman

integer                     :: i,j,k,k1,k2,ivar

real                        :: vp_val(nz),vpt_val(nz,nvar)
real                        :: kvp_val(nz,nvar2),kvpt_val(nz)
real                        :: i0_val(nz,nvar2),i0t_val(nz,nvar)
real                        :: temra1d01(nvar)

  if(gropt==0) then
      kvpt_val=0
	  vp_val=0
	  do k=1,nz	        
	    call diagnostic_var(para(k,:),prtvar(k,:),nvar,vp_val(k),vpt_val(k,:),gropt)
      enddo
	  
	  kvp_val=0
	  call kvp(kalman,nz,1,nvar2,vp_val,ratio_vp,gropt,kvp_val,kvpt_val)
	  
	  call trans_i0(prtvar(:,5:7),i0_val,nz,1,nvar2,ratio_vp,gropt)
  
	  do ivar=1,nvar2
	  	kvar(:,ivar)=kvp_val(:,ivar)+i0_val(:,ivar)
	  enddo
  else
  
	  kvp_val=prtvar(:,5:7)
	  kvpt_val=0
	  call kvp(kalman,nz,1,nvar2,vp_val,ratio_vp,gropt,kvp_val,kvpt_val)  

	  do k=1,nz	  
        temra1d01=0      
	    call diagnostic_var(para(k,:),prtvar(k,:),nvar,kvpt_val(k),temra1d01,gropt)	
        vpt_val(k,:)=temra1d01		
      enddo
 
	  i0t_val=0
	  call trans_i0(i0t_val,prtvar,nz,1,nvar,ratio_vp,gropt)

	  do ivar=1,nvar
	  	grdvar(:,ivar)=vpt_val(:,ivar)+i0t_val(:,ivar)
	  enddo  
 
  endif

end subroutine transform_phy


subroutine kvp(kalman,nz,isv1,iev1,vp,ratio_vp,gropt,newv,grdv)

implicit none

integer                       :: nz,isv1,iev1
integer                       :: gropt
real                          :: kalman(nz,nz,isv1:iev1)
real                          :: ratio_vp
real,dimension(nz)            :: vp,grdv
real,dimension(nz,isv1:iev1)  :: newv

integer                   :: k,k1,k2,ivar


  do k=1,nz
    do ivar=isv1,iev1
      do k1=1,nz
        if(gropt==0) newv(k,ivar)=newv(k,ivar)+ratio_vp*kalman(k,k1,ivar)*vp(k1)
        if(gropt>=1) grdv(k1)=grdv(k1)+ratio_vp*kalman(k,k1,ivar)*newv(k,ivar)	    
      enddo
    enddo
  enddo


end subroutine kvp

subroutine trans_i0(dvar,newv,nz,isv1,iev1,ratio_vp,gropt)

implicit none

integer                        :: nz,isv1,iev1
real,dimension(nz,isv1:iev1)   :: dvar,newv
real                           :: ratio_vp
integer                        :: gropt

integer                        :: k,ivar

do ivar=isv1,iev1
  do k=1,nz
    if(gropt==0) newv(k,ivar)=(1-ratio_vp)*dvar(k,ivar)
    if(gropt>=1) dvar(k,ivar)=(1-ratio_vp)*newv(k,ivar)
  enddo
enddo


end subroutine trans_i0


subroutine diagnostic_var(para,var,nvar,diagvar,gvar,gropt)


implicit none
integer                   :: nvar 
real,dimension(nvar)      :: para,var,gvar
real                      :: diagvar
integer                   :: gropt

integer                   :: i

do i=1,nvar
  if(gropt==0) diagvar=diagvar+para(i)*var(i)
  if(gropt>=1) gvar(i)=gvar(i)+para(i)*diagvar
enddo


end subroutine diagnostic_var
  
end subroutine thermoist_trans

