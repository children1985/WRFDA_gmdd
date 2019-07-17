program thermoist_trans


implicit none

integer,parameter         :: nx=79,ny=79,nz=51,nvar=7
character(len=256)        :: filepath,be_dvar
character(len=16)         :: varname(nvar)

real                      :: bakvar(nx,ny,nz,nvar)
real                      :: prtvar(nz,nvar)
real                      :: grdvar(nx,ny,nz,nvar)
real                      :: kalman(nz,nz,3)
real                      :: kvar(nx,ny,nz,3)
integer                   :: i,j,k,k1,k2,ivar
integer                   :: gropt

real                      :: chkadj_rhs,chkadj_lhs


  data varname/'T','QVAPOR','QCLOUD','QICE','QRAIN','QSNOW','QGRAUP'/

  open(1,file='filelist.dat')
  do i=1,40
    read(1,'(A)') filepath
  enddo
  print*,trim(filepath)
  close(1)
  do i=1,nvar
    print*,trim(varname(i))
    call  readwrf(filepath,varname(i),nx,ny,nz,bakvar(:,:,:,i),0)
  enddo 

!  open(1,file='be_dvar.dat',form='binary')
!  read(1) kalman
!  close(1)
  
  open(1,file='be_dvar_ascii.dat')
  read(1,*) kalman
  close(1)
  
  gropt=0
  grdvar=0	
  kvar=0
  do j=1,ny  !1,ny
    do i=1,nx  !1,nx
      prtvar=0	
	  do ivar=5,7
	    prtvar(:,ivar)=1e-1*bakvar(i,j,:,ivar)
	  enddo
	  kvar(i,j,:,:)=0  
	  call transform(bakvar(i,j,:,:),kalman,gropt,nz,nvar,prtvar,kvar(i,j,:,:),grdvar(i,j,:,:))	 
	enddo
  enddo
  chkadj_lhs=0
  do j=1,ny  !1,ny
    do i=1,nx  !1,nx
	  do k=1,nz
	    do ivar=1,3
          chkadj_lhs=chkadj_lhs+kvar(i,j,k,ivar)**2
		  !if(kvar(i,j,k,ivar).ne.1e-1*bakvar(i,j,k,ivar+4)) print*,"i,j,k,ivar,val",i,j,k,ivar,kvar(i,j,k,ivar),1e-1*bakvar(i,j,k,ivar+4)
		enddo
	  enddo
	enddo
  enddo
  print*,"check adjoint lhs <y,y>:",chkadj_lhs
  print*,"check ",sum((1e-1*bakvar(:,:,:,5:7))**2)
  gropt=1
  grdvar=0	
  do j=1,ny  !1,ny
    do i=1,nx  !1,nx
      prtvar=0	
	  prtvar(:,5:7)=kvar(i,j,:,:)
	  call transform(bakvar(i,j,:,:),kalman,gropt,nz,nvar,prtvar,kvar(i,j,:,:),grdvar(i,j,:,:))	  
	enddo
  enddo

  chkadj_rhs=0
  do j=1,ny  !1,ny
    do i=1,nx  !1,nx
	  prtvar=0	
	  do ivar=5,7
	    prtvar(:,ivar)=1e-1*bakvar(i,j,:,ivar)
	  enddo
	  do k=1,nz
	    do ivar=1,nvar
          chkadj_rhs=chkadj_rhs+prtvar(k,ivar)*grdvar(i,j,k,ivar)
		  !if(prtvar(k,ivar).ne.grdvar(i,j,k,ivar)) print*,"i,j,k,ivar,val",i,j,k,ivar,prtvar(k,ivar),grdvar(i,j,k,ivar)
		enddo
	  enddo
	enddo
  enddo
  print*,"check adjoint rhs <x,A*y>:",chkadj_rhs

end program thermoist_trans


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

  ratio_vp=0.0

  do k=1,nz
    paravp(k,1)=1./(bakvar(k,1)+300.0)
    paravp(k,2)=0.61
    paravp(k,3:7)=-1.0	
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
	  vpt_val=0
	  call kvp(kalman,nz,1,nvar2,vp_val,ratio_vp,gropt,kvp_val,kvpt_val)  

	  do k=1,nz	  
        temra1d01=0      
	    call diagnostic_var(para(k,:),prtvar(k,:),nvar,kvpt_val(k),temra1d01,gropt)	
        vpt_val(k,:)=temra1d01		
      enddo
 
	  i0t_val=0
	  call trans_i0(i0t_val(:,5:7),prtvar(:,5:7),nz,1,nvar2,ratio_vp,gropt)

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
    if(gropt==0) newv(k,ivar)=(1.0-ratio_vp)*dvar(k,ivar)
    if(gropt>=1) dvar(k,ivar)=dvar(k,ivar)+(1.0-ratio_vp)*newv(k,ivar)
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