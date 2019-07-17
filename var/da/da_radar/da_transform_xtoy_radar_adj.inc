subroutine da_transform_xtoy_radar_adj(grid, iv, jo_grad_y, jo_grad_x)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !    Updated for Analysis on Arakawa-C grid
   !    Author: Syed RH Rizvi,  MMM/ESSL/NCAR,  Date: 10/22/2008
   !-----------------------------------------------------------------------

   !------------------------------------------------------------------------
   ! This subroutine is the adjoint of Doppler radar observation operators.
   !------------------------------------------------------------------------

   implicit none

   type (domain),  intent(in)    :: grid
   type (iv_type), intent(in)    :: iv          ! obs. inc vector (o-b).
   type (y_type) , intent(inout) :: jo_grad_y   ! grad_y(jo)
   type (x_type) , intent(inout) :: jo_grad_x   ! grad_x(jo)

   integer :: k  ! Index dimension.

   integer :: n

   real, allocatable :: model_p(:,:)
   real, allocatable :: model_u(:,:)
   real, allocatable :: model_v(:,:)
   real, allocatable :: model_w(:,:)
   real, allocatable :: model_qrn(:,:)
   real, allocatable :: model_qcw(:,:)
   real, allocatable :: model_qrnb(:,:)
   real, allocatable :: model_ps(:)

   !------------------------
   !  for jung et al 2008
   !------------------------   
   real, allocatable :: model_qsnb(:,:)
   real, allocatable :: model_qgrb(:,:)      
   real, allocatable :: model_qnrb(:,:)  
   real, allocatable :: model_qnsb(:,:) 
   real, allocatable :: model_qngb(:,:) 
   real, allocatable :: model_qnr(:,:)  
   real, allocatable :: model_qns(:,:) 
   real, allocatable :: model_qng(:,:) 
   !------------------------ 

   real, allocatable :: model_qsn(:,:)
   real, allocatable :: model_qgr(:,:)
   real, allocatable :: model_qv(:,:)
   real, allocatable :: model_qvb(:,:)
   real, allocatable :: model_t(:,:)
   real, allocatable :: model_tb(:,:)

   real    :: xr,yr,zr

   real    :: alog10,qrn1,qsn1,qgr1

   !------------------------
   !  for jung et al 2008
   !------------------------
   real    :: qvp(1,1,1),qra(1,1,1),qsn(1,1,1),qgr(1,1,1) ! mixing ratio
                                                          ! to match the type of argument in subroutine
   real    :: dqra(1,1,1),dqsn(1,1,1),dqgr(1,1,1),dtmk(1,1,1),dqvp(1,1,1)  
   real    :: dqnr(1,1,1),dqns(1,1,1),dqng(1,1,1) 
   real    :: zmm(1,1,1),zmm_ref(1,1,1)                                                       
   real    :: qnr(1,1,1),qns(1,1,1),qng(1,1,1)            ! number concentration of rain snow and graupel
   real    :: tmk(1,1,1),prs(1,1,1)                       ! temperature and pressure
   real    :: dbz(1,1,1)                                  ! reflectivity in dBZ
   real    :: rn0_r,rn0_s,rn0_g                           ! intercept parameter of rain snow and graupel
   real    :: rhos,rhog                                   ! density of snow and graupel
   !------------------------   

   if (trace_use) call da_trace_entry("da_transform_xtoy_radar_adj")

   alog10= alog(10.0)

   !------------------------
   !  for graupel  !GODDARD HAIL
   !------------------------
   qnr=0
   qns=0
   qng=0
   rn0_r=8e6
   rn0_s=3e6  !1.6e7  !3e6
   rn0_g=4e5  !2.0e5  !4e5
   rhos=100.0
   rhog=400.0 !917.0   !400.0
   !---hail-------------
   !rn0_r=8e6
   !rn0_s=1.67e7  !3e6
   !rn0_g=2.0e5   !4e5
   !rhos=100.0
   !rhog=917.0    !400.0
   !------------------------  

   allocate (model_p(iv%info(radar)%max_lev,iv%info(radar)%n1:iv%info(radar)%n2))
   allocate (model_u(iv%info(radar)%max_lev,iv%info(radar)%n1:iv%info(radar)%n2))
   allocate (model_v(iv%info(radar)%max_lev,iv%info(radar)%n1:iv%info(radar)%n2))
   allocate (model_w(iv%info(radar)%max_lev,iv%info(radar)%n1:iv%info(radar)%n2))
   allocate (model_qcw(iv%info(radar)%max_lev,iv%info(radar)%n1:iv%info(radar)%n2))
   allocate (model_qrn(iv%info(radar)%max_lev,iv%info(radar)%n1:iv%info(radar)%n2))
   allocate (model_qrnb(iv%info(radar)%max_lev,iv%info(radar)%n1:iv%info(radar)%n2))
   allocate (model_ps(iv%info(radar)%n1:iv%info(radar)%n2))

   !------------------------
   !  for jung et al 2008
   !------------------------   
   allocate (model_qsnb(iv%info(radar)%max_lev,iv%info(radar)%n1:iv%info(radar)%n2))
   allocate (model_qgrb(iv%info(radar)%max_lev,iv%info(radar)%n1:iv%info(radar)%n2))
   allocate (model_qnrb(iv%info(radar)%max_lev,iv%info(radar)%n1:iv%info(radar)%n2))
   allocate (model_qnsb(iv%info(radar)%max_lev,iv%info(radar)%n1:iv%info(radar)%n2))
   allocate (model_qngb(iv%info(radar)%max_lev,iv%info(radar)%n1:iv%info(radar)%n2))
   allocate (model_qnr(iv%info(radar)%max_lev,iv%info(radar)%n1:iv%info(radar)%n2))  
   allocate (model_qns(iv%info(radar)%max_lev,iv%info(radar)%n1:iv%info(radar)%n2))  
   allocate (model_qng(iv%info(radar)%max_lev,iv%info(radar)%n1:iv%info(radar)%n2))        
   !------------------------   
   
   allocate (model_qsn(iv%info(radar)%max_lev,iv%info(radar)%n1:iv%info(radar)%n2))
   allocate (model_qgr(iv%info(radar)%max_lev,iv%info(radar)%n1:iv%info(radar)%n2))
   allocate (model_qv(iv%info(radar)%max_lev,iv%info(radar)%n1:iv%info(radar)%n2))
   allocate (model_qvb(iv%info(radar)%max_lev,iv%info(radar)%n1:iv%info(radar)%n2))
   allocate (model_t(iv%info(radar)%max_lev,iv%info(radar)%n1:iv%info(radar)%n2))
   allocate (model_tb(iv%info(radar)%max_lev,iv%info(radar)%n1:iv%info(radar)%n2))

   !basic states
   call da_interp_lin_3d (grid%xb%t,   iv%info(radar), model_tb)
   call da_interp_lin_3d (grid%xb%q,   iv%info(radar), model_qvb)

   ! Needed
   model_u = 0.0
   model_v = 0.0
   model_w = 0.0
   model_qrn = 0.0
   model_qsn = 0.0
   model_qgr = 0.0
   model_qv  = 0.0
   model_qcw = 0.0
   model_t   = 0.0



   ! W_HALF is vertical velocity at half-sigma levels.

   model_ps(iv%info(radar)%n1:iv%info(radar)%n2) = iv%radar(iv%info(radar)%n1:iv%info(radar)%n2)%model_ps 

   do n=iv%info(radar)%n1,iv%info(radar)%n2

      ! [1.7] Calculate rv and rf at OBS location

      xr = grid%xb%ds * (iv%info(radar)%x(1,n) - iv%radar(n)%stn_loc%x)
      yr = grid%xb%ds * (iv%info(radar)%y(1,n) - iv%radar(n)%stn_loc%y)

      model_qrnb(1:iv%info(radar)%levels(n),n) = iv%radar(n)%model_qrn(1:iv%info(radar)%levels(n))
      model_p   (1:iv%info(radar)%levels(n),n) = iv%radar(n)%model_p(1:iv%info(radar)%levels(n))

      !------------------------
      !  for jung et al 2008
      !------------------------           
      model_qsnb(1:iv%info(radar)%levels(n),n) = iv%radar(n)%model_qsn(1:iv%info(radar)%levels(n))
      model_qgrb(1:iv%info(radar)%levels(n),n) = iv%radar(n)%model_qgr(1:iv%info(radar)%levels(n))
      model_qnrb(1:iv%info(radar)%levels(n),n) = 0
      model_qnsb(1:iv%info(radar)%levels(n),n) = 0
      model_qngb(1:iv%info(radar)%levels(n),n) = 0
      !------------------------    

      do k = 1,iv%info(radar)%levels(n)
         if (iv % radar(n) % height_qc(k) /= below_model_surface .and.  &
              iv % radar(n) % height_qc(k) /= above_model_lid) then

            if ((use_radar_rf.or.use_radar_rhv).and.radar_rf_opt>=2) then
               if (iv % radar(n) % rf(k) % qc >= obs_qc_pointer) then
                 select case (radar_rf_opt)

				 case(2)
                   !print*,"calculate dbz at",iv%info(radar)%x(1,n),iv%info(radar)%y(1,n),iv%radar(n)%height(k)
                   qvp(1,1,1)=model_qvb(k,n)
                   qra(1,1,1)=model_qrnb(k,n)
                   qsn(1,1,1)=model_qsnb(k,n)
                   qgr(1,1,1)=model_qgrb(k,n)
                   qnr(1,1,1)=0
                   qns(1,1,1)=0
                   qng(1,1,1)=0
                   tmk(1,1,1)=model_tb(k,n) 
                   prs(1,1,1)=model_p(k,n)
                   dqra(1,1,1)=model_qrn(k,n)
                   dqsn(1,1,1)=model_qsn(k,n)
                   dqgr(1,1,1)=model_qgr(k,n)
                   dtmk(1,1,1)=model_t(k,n)
                   dqvp(1,1,1)=model_qv(k,n)					   
                   dqnr=0
                   dqns=0
                   dqng=0
                   zmm=iv % radar(n) % zmm(k) % inv
                   zmm_ref=radar_rf_rscl*jo_grad_y%radar(n)%rf(k)

                   call RadZIceVar(qvp,qra,qsn,qgr,qnr,qns,qng,tmk,prs,dbz,        &
                                   1,1,1,0,0,0,rn0_r,rn0_s,rn0_g,                  &
                                   rhos,rhog,dtmk,dqvp,dqra,dqsn,dqgr,dqnr,dqns,dqng,zmm,2,  &
                                   2,zmm_ref)                                  
                   model_qrn(k,n) = dqra(1,1,1)
                   model_qsn(k,n) = dqsn(1,1,1)
                   model_qgr(k,n) = dqgr(1,1,1)
				   model_t(k,n)   = dtmk(1,1,1)
				   model_qv(k,n)  = dqvp(1,1,1)
                  !print*,"tangent linear dbz value:",dbz(1,1,1) 
                 case(3) 
				   model_qrn(k,n) = model_qrn(k,n) + radar_rf_rscl*leh2/(model_qrnb(k,n)*alog10) * jo_grad_y%radar(n)%rf(k) 			  
                 end select				 
               end if
            end if			  
			  
			  
            if (use_radar_rf .and. .not.use_radar_rhv.and.radar_rf_opt==1) then
               if (iv % radar(n) % rf(k) % qc >= obs_qc_pointer) then                
                 model_qrn(k,n) = model_qrn(k,n) + leh2/(model_qrnb(k,n)*alog10) * jo_grad_y%radar(n)%rf(k)                                
               end if
            end if

            if (.not.use_radar_rf .and. use_radar_rhv.and.radar_rf_opt==1) then
               if (iv % radar(n) % rrn(k) % qc >= obs_qc_pointer) then
                  model_qrn(k,n) = model_qrn(k,n) + jo_grad_y%radar(n)%rrn(k)
               end if
               if (iv % radar(n) % rsn(k) % qc >= obs_qc_pointer) then
                  model_qsn(k,n) = model_qsn(k,n) + jo_grad_y%radar(n)%rsn(k)
               end if
               if (iv % radar(n) % rgr(k) % qc >= obs_qc_pointer) then
                  model_qgr(k,n) = model_qgr(k,n) + jo_grad_y%radar(n)%rgr(k)
               end if
            end if

            if (use_radar_rqv.and.radar_rf_opt==1) then
               if (iv % radar(n) % rqv(k) % qc >= obs_qc_pointer) then
               !TL  y%radar(n)%rqv(k) = y%radar(n)%rqv(k) +( es_beta*es_gamma/(model_tb(k,n)+es_gamma)**2.0 )*model_qvb(k,n)*model_t(k,n)
                  model_qv(k,n) = model_qv(k,n) + jo_grad_y%radar(n)%rqv(k)
                  model_t(k,n)  = model_t(k,n)  + (es_beta*es_gamma/(model_tb(k,n)+es_gamma)**2.0)*model_qvb(k,n)*jo_grad_y%radar(n)%rqv(k)
               end if
            end if


            if (use_radar_rv) then
               if (iv % radar(n) % rv(k) % qc >= obs_qc_pointer) then
                  zr=iv%radar(n)%height(k) - iv%radar(n)%stn_loc%elv
                  jo_grad_y%radar(n)%rv(k)=jo_grad_y%radar(n)%rv(k)*radar_rv_rscl
                  call da_radial_velocity_adj(jo_grad_y%radar(n)%rv(k), &
                     model_p(k,n), model_u(k,n), model_v(k,n), model_w(k,n),  &
                     model_qrn(k,n), model_ps(n), xr, yr, zr, model_qrnb(k,n))
                  
               end if
            end if
         end if
      end do
      jo_grad_y%radar(n)%rv(:) = 0.0
      jo_grad_y%radar(n)%rf(:) = 0.0
      jo_grad_y%radar(n)%rrn(:)= 0.0
      jo_grad_y%radar(n)%rsn(:)= 0.0
      jo_grad_y%radar(n)%rgr(:)= 0.0
      jo_grad_y%radar(n)%rqv(:)= 0.0
   end do ! n

   ! [1.6] Interpolate horizontally from crs points:

   call da_interp_lin_3d_adj (jo_grad_x % wh,  iv%info(radar), model_w)
   if ( cloud_cv_options >= 1 ) then
      call da_interp_lin_3d_adj (jo_grad_x % qrn, iv%info(radar), model_qrn)
      call da_interp_lin_3d_adj (jo_grad_x % qcw, iv%info(radar), model_qcw)
      if ( cloud_cv_options >= 2 ) then
         call da_interp_lin_3d_adj (jo_grad_x % qsn, iv%info(radar), model_qsn)
         call da_interp_lin_3d_adj (jo_grad_x % qgr, iv%info(radar), model_qgr)
      end if
   end if
   call da_interp_lin_3d_adj (jo_grad_x % q,   iv%info(radar), model_qv)
   call da_interp_lin_3d_adj (jo_grad_x % t,   iv%info(radar), model_t)
#ifdef A2C
   call da_interp_lin_3d_adj (jo_grad_x % v,   iv%info(radar), model_v,'u')
   call da_interp_lin_3d_adj (jo_grad_x % u,   iv%info(radar), model_u,'v')
#else
   call da_interp_lin_3d_adj (jo_grad_x % v,   iv%info(radar), model_v)
   call da_interp_lin_3d_adj (jo_grad_x % u,   iv%info(radar), model_u)
#endif

   deallocate (model_p)
   deallocate (model_u)
   deallocate (model_v)
   deallocate (model_w)
   deallocate (model_qrn)
   deallocate (model_qrnb)
   deallocate (model_ps)
   deallocate (model_qv)
   deallocate (model_qvb)
   deallocate (model_t)
   deallocate (model_tb)
   deallocate (model_qsn)
   deallocate (model_qgr)
   deallocate (model_qcw)

   deallocate (model_qsnb)
   deallocate (model_qgrb)
   deallocate (model_qnr)
   deallocate (model_qns)
   deallocate (model_qng)
   deallocate (model_qnrb)
   deallocate (model_qnsb)
   deallocate (model_qngb)   
   print*,"pass da_transform_xtoy_radar"
   
   if (trace_use) call da_trace_exit("da_transform_xtoy_radar_adj")

end subroutine da_transform_xtoy_radar_adj


