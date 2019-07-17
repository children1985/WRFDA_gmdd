module modeltype

  integer,parameter  :: nvar=7
  integer,parameter  :: ntime=2
  character(len=16)  :: varname(nvar)
  integer            :: dim(nvar,3)
  type var3dtype
    real,allocatable :: val(:,:,:)
  end type var3dtype
  type(var3dtype)   :: var3d(nvar)
  type(var3dtype)   :: var4d(nvar,ntime)

  data varname/'T','QVAPOR','QCLOUD','QICE','QRAIN','QSNOW','QGRAUP'/

end module modeltype

program thermoist_cov

use modeltype
implicit none

include 'netcdf.inc'

character(len=256) :: filepath,be_dvar
real,allocatable   :: vartem(:,:,:)
real,allocatable   :: cov(:,:,:),cov_dvar(:,:),stdx(:,:),kalman(:,:,:)
real,allocatable   :: dvar(:,:,:),gvar(:,:,:)

real               :: mean2d
logical            :: alive
integer            :: i,j,k,k1,k2,k3,ivar,itime
integer            :: scanopt
integer            :: vertloc=1
real               :: ratio_vp=0.5

real               :: dvarstd
integer            :: temi01,temi02,temi03,temi04,temi05
real               :: temr01,temr02,temr03,temr04,temr05 
real,allocatable   :: temra1d01(:),temra1d02(:),temra1d03(:),temra1d04(:)

real,allocatable   :: temra2d01(:,:),temra2d02(:,:),temra2d03(:,:),temra2d04(:,:)
real,allocatable   :: temra3d01(:,:,:),temra3d02(:,:,:),temra3d03(:,:,:)
real,allocatable   :: temra4d01(:,:,:,:)
real,allocatable   :: para(:)

real,allocatable,dimension(:,:)    :: prtvar,grdvar
real,allocatable,dimension(:,:)    :: kvar
real,allocatable,dimension(:,:)    :: paravp


integer            :: itest,jtest 
integer            :: gropt

real               :: chkadj_rhs,chkadj_lhs
real               :: para1,para2,para3

type(var3dtype)    :: prt3d(nvar)
type(var3dtype)    :: vp3d(nvar)
type(var3dtype)    :: new3d(nvar)

!!! read model data
  filepath="/glade/u/home/wangs/WRFV3_mwsm6/run_mwsm6/wrfout_d01_0001-01-01_01:50:00"

  allocate(vartem(1,1,1))

  do itime=1,ntime
    open(1,file='filelist.dat')
    if(itime==1) j=20
	if(itime==2) j=21
    do i=1,j
	  read(1,'(A)') filepath
	enddo
	!read(1,*)  filepath
	print*,trim(filepath)
	close(1)
    do i=1,nvar
      print*,trim(varname(i))
      scanopt=1
      call  readwrf(filepath,varname(i),dim(i,1),dim(i,2),dim(i,3),vartem,scanopt)
      if(itime==1)allocate(var3d(i)%val(dim(i,1),dim(i,2),dim(i,3)))
      allocate(var4d(i,itime)%val(dim(i,1),dim(i,2),dim(i,3)))
      scanopt=0
      call  readwrf(filepath,varname(i),dim(i,1),dim(i,2),dim(i,3),var3d(i)%val,scanopt)
	  var4d(i,itime)%val=var3d(i)%val
      if(itime==1)allocate(prt3d(i)%val(dim(i,1),dim(i,2),dim(i,3)))
      if(itime==1)allocate(vp3d(i)%val(dim(i,1),dim(i,2),dim(i,3)))
      if(itime==1)allocate(new3d(i)%val(dim(i,1),dim(i,2),dim(i,3)))
    enddo 
  enddo
  itime=2
  do i=1,nvar
    var3d(i)%val=var4d(i,itime)%val-var4d(i,itime-1)%val
  enddo
!  var3d(i)%val=var4d(i,1)%val
  allocate(temra1d01(nvar))
!  do i=1,2
!    do k=1,dim(1,3)
!	  temra1d01(i)=mean2d(var4d(i,1)%val(:,:,k),dim(1,1),dim(1,2))
!      var3d(i)%val(:,:,k)=var4d(i,1)%val(:,:,k)-temra1d01(i)
!	enddo
!  enddo  
  

  itest=maxloc(maxval(var4d(7,1)%val(:,:,15),2),1)
  jtest=maxloc(maxval(var4d(7,1)%val(:,:,15),1),1)

  print*,"test point:",itest,jtest

  allocate(dvar(dim(1,1),dim(1,2),dim(1,3)))
  allocate(gvar(dim(1,1),dim(1,2),dim(1,3)))
  allocate(cov(dim(1,3),dim(1,3),nvar))
  allocate(cov_dvar(dim(1,3),dim(1,3)))
  allocate(stdx(dim(1,3),nvar))

  allocate(kalman(dim(1,3),dim(1,3),nvar))

! /'T','QVAPOR','QCLOUD','QICE','QRAIN','QSNOW','QGRAUP'/

  allocate(temra1d02(dim(1,3)))
  allocate(temra2d01(dim(1,1),dim(1,2)))
  allocate(para(nvar))
  gvar=0
  gropt=0
  para1=2.e-1
  
  do k=1,dim(1,3)

    do j=1,dim(1,2)
	  do i=1,dim(1,1)
	    para(1)=para1/(var4d(1,1)%val(i,j,k)+300.0)  !var4d(2,1)%val(i,j,k)*sum(var4d(3:5,1)%val(i,j,k))   !1./(var4d(1,1)%val(i,j,k)+300.0)
		para(2)=0.61  !(var4d(1,1)%val(i,j,k)+300.0)*sum(var4d(3:5,1)%val(i,j,k))   !0.61
		para(3:7)=1  !(var4d(1,1)%val(i,j,k)+300.0)*var4d(2,1)%val(i,j,k) !-1.0
		dvar(i,j,k)=0
	    call diagnostic_var(para,var3d(:)%val(i,j,k),nvar,dvar(i,j,k),gvar(i,j,k),gropt)
      enddo
    enddo
  enddo


  cov=0
  cov_dvar=0
  stdx=0
  dvarstd=0

  do k1=1,dim(1,3)
    !print*,"processing covariance for level",k1
    temr01=mean2d(dvar(:,:,k1),dim(1,1),dim(1,2))
    do k2=1,dim(1,3)
      if(abs(k1-k2)<vertloc) then
        temr02=mean2d(dvar(:,:,k2),dim(1,1),dim(1,2))
        do ivar=1,nvar
          temra1d01(ivar)=mean2d(var3d(ivar)%val(:,:,k2),dim(1,1),dim(1,2))
        enddo
        do j=1,dim(1,2)
          do i=1,dim(1,1)
            cov_dvar(k1,k2)=cov_dvar(k1,k2)+(dvar(i,j,k1)-temr01)*(dvar(i,j,k2)-temr02)
            do ivar=1,nvar
              cov(k1,k2,ivar)=cov(k1,k2,ivar)+(dvar(i,j,k1)-temr01)*(var3d(ivar)%val(i,j,k2)-temra1d01(ivar))
            enddo
          enddo
        enddo         
      endif
    enddo
    do ivar=1,nvar
      temra1d01(ivar)=mean2d(var3d(ivar)%val(:,:,k1),dim(1,1),dim(1,2))
    enddo
    do j=1,dim(1,2)
      do i=1,dim(1,1)
        do ivar=1,nvar
           stdx(k1,ivar)=stdx(k1,ivar)+(var3d(ivar)%val(i,j,k1)-temra1d01(ivar))**2
        enddo
      enddo
    enddo
!    print('(I4,7F20.8)'),k1,stdx(k1,:)
  enddo

  do k1=1,dim(1,3)
    do k2=1,dim(1,3)
      if(k1==k2) dvarstd=dvarstd+cov_dvar(k1,k2) 
    enddo
  enddo
  dvarstd=(dvarstd/dim(1,3))**0.5
  print*,"dvar std",dvarstd

  print*,"--------correlation on diag line--------"
  print('(10X,7A15)'),varname
  do k1=dim(1,3),1,-1
    do ivar=1,nvar
      if(cov_dvar(k1,k1)==0.0.or.stdx(k1,ivar)==0.0) then
        temra1d01(ivar)=0.0
      else
        temra1d01(ivar)=cov(k1,k1,ivar)/(stdx(k1,ivar)**0.5*cov_dvar(k1,k1)**0.5)
      endif
    enddo
    print('(I4,7F15.8)'),k1,temra1d01   !cov(k1,k1,:)/(stdx(k1,:)**0.5*cov_dvar(k1,k1)**0.5)
  enddo 

  kalman=0
  do k1=1,dim(1,3)
    do k2=1,dim(1,3)
      do ivar=1,nvar
        !kalman(k1,k2,ivar)=exp(-(1.0*(k1-k2)/vertloc)**2)*cov(k1,k1,ivar)/(cov_dvar(k1,k1)+(0.1*dvarstd)**2)
		kalman(k1,k2,ivar)=exp(-(1.0*(k1-k2)/vertloc)**2)*cov(k1,k2,ivar)/(cov_dvar(k1,k2)+(0.1*dvarstd)**2)
      enddo
    enddo
  enddo

  print*,"--------kalman on diag line--------"
  print('(10X,7A15)'),varname
  do k1=dim(1,3),1,-1
     print('(I4,7F15.8)'),k1,kalman(k1,k1,:)
  enddo 

  open(1,file='diag.log')
  do ivar=1,nvar
    write(1,*) "----------------kalman for ",trim(varname(ivar)),"-------------------------"
    do k1=dim(1,3),1,-1
       write(1,'(I4,52F15.8)') k1,kalman(k1,:,ivar)
    enddo
  enddo


  be_dvar='./be_dvar.dat'
  open(1,file=trim(be_dvar),form='binary')
  write(1) kalman(:,:,5:7)
  close(1)
  
  be_dvar='./be_dvar_ascii.dat'
  open(1,file=trim(be_dvar))
  write(1,*) kalman(:,:,5:7)
  close(1)  
  
  print*,size(kalman(:,:,5:7),1),size(kalman(:,:,5:7),2),size(kalman(:,:,5:7),3)
 ! stop("test")

  !========================================================================================
  !========================================================================================
  !            tangent linear and adjoint test
  !========================================================================================
  !========================================================================================
  do ivar=1,nvar
    prt3d(ivar)%val(:,:,:)=0.0
  enddo
  do ivar=5,nvar
    prt3d(ivar)%val(:,:,:)=1e-1*var4d(ivar,1)%val(:,:,:) !var3d(ivar)%val(:,:,:)
    new3d(ivar)%val(:,:,:)=0
  enddo


  allocate(prtvar(dim(1,3),nvar))  
  allocate(grdvar(dim(1,3),nvar))
  allocate(kvar(dim(1,3),3))
  allocate(paravp(dim(1,3),nvar)) 
  

 ! stop("test diagnostic_var") 
  gropt=0
 
  do j=1,dim(1,2) !jtest,jtest !  1,dim(1,2)
    do i=1,dim(1,1) !itest,itest!  1,dim(1,1)


	  do k=1,dim(1,3)
	    paravp(k,1)=para1/(var4d(1,1)%val(i,j,k)+300.0)
		paravp(k,2)=0.61
		paravp(k,3:7)=1.0	
        prtvar(k,:)=prt3d(:)%val(i,j,k)		
	  enddo

      call transform_phy(prtvar,kvar,grdvar,kalman(:,:,5:7),paravp,nvar,3,dim(1,3),gropt,ratio_vp)
	
      k=1
	  do ivar=5,7
	  	new3d(ivar)%val(i,j,:)=kvar(:,k)
		k=k+1
	  enddo

    enddo
  enddo
  !stop("test")
  i=itest
  j=jtest
  !print*,"-------------dvar------------"
  !print*,vp3d(1)%val(i,j,:)
  !print*,"---------------------------"
  
  print*,"---------qr new,old------------"
  do k=dim(1,3),1,-1
     print('(I4,2F12.8)'),k,new3d(5)%val(i,j,k)*1e3,prt3d(5)%val(i,j,k)*1e3
  enddo
  print*,"---------qs new,old------------"
  do k=dim(1,3),1,-1
     print('(I4,2F12.8)'),k,new3d(6)%val(i,j,k)*1e3,prt3d(6)%val(i,j,k)*1e3
  enddo
  print*,"---------qg new,old------------"
  do k=dim(1,3),1,-1
     print('(I4,2F12.8)'),k,new3d(7)%val(i,j,k)*1e3,prt3d(7)%val(i,j,k)*1e3
  enddo
  !print*,"---------theta new,old------------"
  !do k=dim(1,3),1,-1
  !   print('(I4,2F12.8)'),k,new3d(1)%val(i,j,k),prt3d(1)%val(i,j,k)
  !enddo
  !print*,"---------qv new,old------------"
  !do k=dim(1,3),1,-1
  !   print('(I4,2F12.8)'),k,new3d(2)%val(i,j,k)*1e3,prt3d(2)%val(i,j,k)*1e3
  !enddo
  
  chkadj_rhs=0    !<x,A*y>
  chkadj_lhs=0    !<y,y>
  
  print*,dim(1,:)
  
  do j=1,dim(1,2) !jtest,jtest !  1,dim(1,2)
    do i=1,dim(1,1) !itest,itest!  1,dim(1,1)
	  do k=1,dim(1,3)
	    do ivar=1,nvar
          chkadj_lhs=chkadj_lhs+new3d(ivar)%val(i,j,k)**2
		enddo
	  enddo
	enddo
  enddo
  print*,"check adjoint lhs <y,y>:",chkadj_lhs
  !stop
  gropt=1

  do j=1,dim(1,2) !  1,dim(1,2),jtest,jtest
    do i=1,dim(1,1)!  1,dim(1,1),itest,itest

	
	  do k=1,dim(1,3)
	    paravp(k,1)=para1/(var4d(1,1)%val(i,j,k)+300.0)
		paravp(k,2)=0.61
		paravp(k,3:7)=1.0	
        prtvar(k,:)=new3d(:)%val(i,j,k)	
	  enddo

      call transform_phy(prtvar,kvar,grdvar,kalman(:,:,5:7),paravp,nvar,3,dim(1,3),gropt,ratio_vp)	
	
	  do ivar=1,nvar
	  	new3d(ivar)%val(i,j,:)=grdvar(:,ivar)  !grdvar(:,ivar) !temra2d03(:,ivar)+temra2d04(:,ivar)
	  enddo

    enddo
  enddo
!stop
  chkadj_rhs=0
  do j=1,dim(1,2) !  1,dim(1,2),jtest,jtest
    do i=1,dim(1,1)!  1,dim(1,1),itest,itest
	  do k=1,dim(1,3)
	    do ivar=1,nvar
          chkadj_rhs=chkadj_rhs+prt3d(ivar)%val(i,j,k)*new3d(ivar)%val(i,j,k)
		enddo
	  enddo
	enddo
  enddo
  print*,"check adjoint rhs <x,A*y>:",chkadj_rhs

end program thermoist_cov


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

	  !print*,"kvpt_val",kvpt_val
	  
	  do k=1,nz	  
        temra1d01=0      
	    call diagnostic_var(para(k,:),prtvar(k,:),nvar,kvpt_val(k),temra1d01,gropt)	
        vpt_val(k,:)=temra1d01		
      enddo

      !print*,"vpt_val",vpt_val
	  
	  i0t_val=0
	  call trans_i0(i0t_val,prtvar,nz,1,nvar,ratio_vp,gropt)
      !print*,i0t_val
	  !stop
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


real function mean2d(var,nx,ny)

implicit none

integer :: nx,ny
real    :: var(nx,ny)
integer :: i,j

mean2d=0.0
do j=1,ny
  do i=1,nx
    mean2d=mean2d+var(i,j)
  enddo
enddo
mean2d=mean2d/(nx*ny)

end function mean2d

subroutine readwrf(filepath,varname,nx,ny,nz,var,scanopt)

implicit none

include 'netcdf.inc'

integer            :: scanopt
integer            :: nx,ny,nz
character(len=256) :: filepath
logical            :: alive
integer            :: ncid4,varid4
character(len=16)  :: varname
integer            :: ivtype,nDims,dimids(10),nAtts, dims(3) 
integer            :: istart(4),iend(4)
real               :: var(nx,ny,nz)

real,allocatable   :: vartem(:,:,:)
integer            :: i,j,k
integer            :: status
integer            :: nxmin,nymin,nzmin

   inquire(file=trim(filepath),exist=alive) 
   if(.not.alive) stop
   status=nf_open(trim(filepath),nf_write,ncid4) 
   status=nf_inq_varid(ncid4,trim(varname),varid4) 
   dims =1 !first give all dimension value 1
   status=nf_inq_var(ncid4,varid4,trim(varname),ivtype,nDims,dimids,nAtts)
   do i=1,nDims
     status=nf_inq_dimlen(ncid4,dimids(i),dims(i)) !selectively give dimension value
   enddo
   istart        = 1 !first give all dimension value 1
   iend          = 1 !first give all dimension value 1
   do i = 1,nDims 
     iend(i)= dims(i) !selectively give dimension value
   enddo
   if(scanopt==1) then
     print*,iend
     nx=iend(1)
     ny=iend(2)
     nz=iend(3) 
     status=nf_close(ncid4)
     return
   endif
   if(nx.ne.iend(1).or.ny.ne.iend(2).or.nz.ne.iend(3)) then
     print*,"dimension not match"
     print*,"nx,ny,nz",nx,ny,nz
     print*,"iend(1),iend(2),iend(3)",iend(1),iend(2),iend(3)
   endif

   if(allocated(vartem)) deallocate(vartem)
   allocate(vartem(iend(1),iend(2),iend(3)))
   
   status=nf_get_vara_real(ncid4,varid4,istart,iend,vartem)

   status=nf_close(ncid4)
   print*,"max,min:",maxval(maxval(maxval(vartem,1),1),1),minval(minval(minval(vartem,1),1),1)


   var=vartem
   deallocate(vartem)

end subroutine readwrf
