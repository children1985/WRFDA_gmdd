
module modeltype

  integer,parameter  :: nvar=2
  character(len=32)  :: varname(nvar)
  integer            :: dim(nvar,3)
  type var3dtype
    real,allocatable :: val(:,:,:)
  end type var3dtype
  type(var3dtype)   :: var3d(nvar)

  data varname/'XLONG','XLAT'/

end module modeltype

module st4data

  integer,parameter  :: st4nvar=3
  character(len=32)  :: st4name(st4nvar)
  integer            :: st4dim(st4nvar,3)
  type st4type
    real,allocatable :: val(:,:,:)
  end type st4type
  type(st4type)   :: st4info(st4nvar)

  data st4name/'A_PCP_GDS5_SFC_acc1h','g5_lat_0','g5_lon_1'/

end module st4data

program st4towrf

  use st4data
  use modeltype
  implicit none
  
  include 'netcdf.inc'

  character(len=256) :: filepath,st4path
  real,allocatable   :: vartem(:,:,:)
  real,allocatable   :: st4m(:,:,:)

  integer            :: nx,ny,nz
  logical            :: alive,alive1,alive2
  integer            :: i,j,k
  integer            :: scanopt

  integer            :: temi01,temi02,temi03,temi04,temi05
  real               :: temr01,temr02,temr03,temr04,temr05 


!--------------------------------------------------------!
!                                                        ! 
!                     Read wrf data                      !
!                                                        !
!--------------------------------------------------------!

  inquire(file='./config.dat',exist=alive)
  inquire(file='./wrfgrid.nc',exist=alive1)
  inquire(file='./cwbobs.dat',exist=alive2)  
  if(alive2) then
    !filepath="wrfgrid.nc"
    st4path="cwbobs.dat"
  endif
  if(alive1) then
    filepath="wrfgrid.nc"
    print('("read in file:",A)'),trim(filepath)
  endif    
  if(alive) then
    open(1,file='config.dat')
    read(1,*) st4path
    close(1)
  endif
  if(.not.alive1.and..not.alive2.and..not.alive) then
    print*,"need config.dat or fg and cwbobs,stop"
    stop
  endif

  allocate(vartem(1,1,1))

  do i=1,nvar
    print*,trim(varname(i))
    scanopt=1
    call  readnc(filepath,varname(i),dim(i,1),dim(i,2),dim(i,3),vartem,scanopt)
    allocate(var3d(i)%val(dim(i,1),dim(i,2),dim(i,3)))
    scanopt=0
    call  readnc(filepath,varname(i),dim(i,1),dim(i,2),dim(i,3),var3d(i)%val,scanopt)
  enddo 
  allocate(st4m(dim(1,1),dim(1,2),dim(1,3)))
  print*,dim(1,1),dim(1,2),dim(1,3)
!--------------------------------------------------------!
!                                                        ! 
!                     Read obs data                      !
!                                                        !
!--------------------------------------------------------!

  !do i=1,st4nvar
  !  print*,trim(st4name(i))
  !  scanopt=1
  !  call  readnc(st4path,st4name(i),st4dim(i,1),st4dim(i,2),st4dim(i,3),vartem,scanopt)
  !  allocate(st4info(i)%val(st4dim(i,1),st4dim(i,2),st4dim(i,3)))
  !  scanopt=0
  !  call  readnc(st4path,st4name(i),st4dim(i,1),st4dim(i,2),st4dim(i,3),st4info(i)%val,scanopt)
  !enddo 
  do i=1,st4nvar
    st4dim(i,1)=441
    st4dim(i,2)=561
    st4dim(i,3)=1
    allocate(st4info(i)%val(st4dim(i,1),st4dim(i,2),st4dim(i,3)))
  enddo 

  nx=441
  ny=561
  do j=1,ny
    do i=1,nx
       st4info(2)%val(i,j,1)=20.0+(j-1)*0.0125
       st4info(3)%val(i,j,1)=118.0+(i-1)*0.0125
    enddo
  enddo

  print*,"max lat,min lat",maxval(maxval(st4info(2)%val(:,:,1),1),1), &
                           minval(minval(st4info(2)%val(:,:,1),1),1)

  print*,"max lon,min lon",maxval(maxval(st4info(3)%val(:,:,1),1),1), &
                           minval(minval(st4info(3)%val(:,:,1),1),1)

  open(1001,file=trim(st4path))
  do j=1,ny
    read(1001,*) st4info(1)%val(:,j,1)
  enddo
  print*,"max cwbrain,min cwb",maxval(maxval(st4info(1)%val(:,:,1),1),1), &
                               minval(minval(st4info(1)%val(:,:,1),1),1)
!  stop("check read cwb rainobs")
!--------------------------------------------------------!
!                                                        ! 
!               interpolate obs to wrf grid              !
!                                                        !
!--------------------------------------------------------!
  call intpl_obs2wrf(st4m,dim(1,1),dim(1,2),dim(1,3))
  print*,"max st4m,min st4m",maxval(maxval(st4m(:,:,1),1),1), &
                             minval(minval(st4m(:,:,1),1),1)
!--------------------------------------------------------!
!                                                        ! 
!               write out interpolated obs               !
!                                                        !
!--------------------------------------------------------!
  
  print*,"write out interpolated obs"
  open(100,file='st4obs.bin',form='unformatted')
  write(100) st4m
  close(100)

  print*,"write out model lat lon"
  open(100,file='model_latlon.bin',form='unformatted')
  write(100) var3d(2)%val
  write(100) var3d(1)%val
  close(100)
  
  print*,"Finished"
   
end program st4towrf


subroutine intpl_obs2wrf(var,nx,ny,nz) 

  use st4data
  use modeltype
  implicit none  

  integer            :: nx,ny,nz
  real               :: var(nx,ny,nz)

  real               :: weia,weib
  real               :: zp
  real               :: dh
  integer            :: i,j,ii,jj
  real               :: thres_dh=10000.0
  integer            :: collocate_point,ic,jc

  real               :: mlatmax,mlatmin,mlonmax,mlonmin
  real               :: olatmax,olatmin,olonmax,olonmin
  integer            :: ios,ioe,jos,joe
  real               :: dismin
  real               :: value_close  

  print*,"interpolate obs to wrf grid"
!'XLONG_M','XLAT_M'
!'A_PCP_GDS5_SFC_acc1h','g5_lat_0','g5_lon_1'
  mlatmax=maxval(var3d(2)%val(:,ny,1))
  mlatmin=minval(var3d(2)%val(:,1,1))
  mlonmax=maxval(var3d(1)%val(nx,:,1))
  mlonmin=minval(var3d(1)%val(1,:,1))  
  print*,"model mlatmax,mlatmin,mlonmax,mlonmin:",mlatmax,mlatmin,mlonmax,mlonmin

  do i=st4dim(1,2),1,-1
    olatmax=minval(st4info(2)%val(:,i,1))
    !print*,i,olatmax,mlatmax
    if(mlatmax>olatmax) then
      joe=i
      exit
    endif
  enddo
  do i=1,st4dim(1,2)
    olatmin=maxval(st4info(2)%val(:,i,1))
    if(mlatmin<olatmin) then
      jos=i
      exit
    endif
  enddo
  do i=st4dim(1,1),1,-1
    olonmax=minval(st4info(3)%val(i,:,1))
    if(mlonmax>olonmax) then
      ioe=i
      exit
    endif
  enddo
  do i=1,st4dim(1,1)
    olonmin=maxval(st4info(3)%val(i,:,1))
    if(mlonmin<olonmin) then
      ios=i
      exit
    endif
  enddo
  
  print*,"ios,ioe,jos,joe:",ios,ioe,jos,joe
  print*,"olatmax,olatmin,olonmax,olonmin:",olatmax,olatmin,olonmax,olonmin

  var=0
  do j=1,ny
    do i=1,nx
      weia=0.0
      weib=0.0
      collocate_point=0
      ic=1
      jc=1
      dismin=9999999.0
      value_close=-999.0
      do jj=jos,joe  !1,st4dim(1,2)
        do ii=ios,ioe  !1,st4dim(1,1)
          if(abs(var3d(2)%val(i,j,1)-st4info(2)%val(ii,jj,1))>0.05.or.  &
             abs(var3d(1)%val(i,j,1)-st4info(3)%val(ii,jj,1))>0.05.or.  &
             st4info(1)%val(ii,jj,1)<0.0) then
            cycle
          endif
          zp=0.0
          dh=dstn(var3d(2)%val(i,j,1),var3d(1)%val(i,j,1),zp         &
                 ,st4info(2)%val(ii,jj,1),st4info(3)%val(ii,jj,1),zp     & 
                 )            
          if(dh>thres_dh) then
            cycle
          else
            if(dh<dismin) then
              dismin=dh
              value_close=st4info(1)%val(ii,jj,1)
            endif
          endif
          
          !if(dh>1.0) then
          !  weia=weia+st4info(1)%val(ii,jj,1)*(1./dh)**2
          !  weib=weib+(1.0/dh)**2
          !else
           ! collocate_point=1
           ! ic=ii
           ! jc=jj
           ! exit
          ! weia=weia+st4info(1)%val(ii,jj,1)*(1./1.)**2
          ! weib=weib+(1.0/1.0)**2
          !endif
        enddo
        !if(collocate_point==1) exit
      enddo
      var(i,j,1)=value_close
      !if(weib>0.0) then
      !  var(i,j,1)=weia/weib         
      !else if(collocate_point==1) then
      !  var(i,j,1)=st4info(1)%val(ic,jc,1)
      !endif
      !if(var(i,j,1)>10.0)print*,i,j,var(i,j,1) 
    enddo
  enddo
 
  contains

  real function dstn(lat1,lon1,z1,lat2,lon2,z2) 

    implicit none

    real      :: lat1,lon1,lat2,lon2,z1,z2
    real      :: earthr=6371e3
    real      :: dishorz,disvert
    real      :: c
    real      :: pi

    pi=2*asin(1.0)

    !C = sin(MLatA)*sin(MLatB)*cos(MLonA-MLonB) + cos(MLatA)*cos(MLatB)
    !Distance = R*Arccos(C)*Pi/180

    c=sin(lat1)*sin(lat2)*cos(lon1-lon2)+cos(lat1)*cos(lat2)    
    dishorz=earthr*acos(c)*pi/180.0

    disvert=abs(z1-z2)

    dstn=(dishorz**2+disvert**2)**0.5

  end function dstn

end subroutine intpl_obs2wrf


subroutine readnc(filepath,varname,nx,ny,nz,var,scanopt)

implicit none

include 'netcdf.inc'

integer            :: scanopt
integer            :: nx,ny,nz
character(len=256) :: filepath
logical            :: alive
integer            :: ncid4,varid4
character(len=32)  :: varname
integer            :: ivtype,nDims,dimids(10),nAtts, dims(3) 
integer            :: istart(4),iend(4)
real               :: var(nx,ny,nz)

character*16       :: varname_tem
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

end subroutine readnc





