
module modeltype

  integer,parameter  :: nvar=1
  character(len=32)  :: varname(nvar)
  integer            :: dim(nvar,3)
  type var3dtype
    real,allocatable :: val(:,:,:)
  end type var3dtype
  type(var3dtype)   :: var3d(nvar)

  data varname/'RAINNC'/

end module modeltype


program st4vswrf_fss

  use modeltype
  implicit none
  
  include 'netcdf.inc'

  character(len=256) :: filepath(2),st4path
  real,allocatable   :: vartem(:,:,:)
  real,allocatable   :: st4m(:,:,:),fcst(:,:,:)


  logical            :: alive1,alive2,alive3
  integer            :: i,j,k
  integer            :: scanopt

  integer            :: temi01,temi02,temi03,temi04,temi05
  real               :: temr01,temr02,temr03,temr04,temr05 
  integer,parameter  :: maxthres=10
  integer            :: nthres
  real               :: prthres(maxthres)  
  real               :: dx=3e3,dy=3e3,radii=50e3
  real               :: fss(1)  
  character*12       :: timlab ! ccyymmddhhnn
  integer            :: unbiasopt=0
  namelist/options/  nthres,prthres,dx,dy,radii,timlab,unbiasopt
!--------------------------------------------------------!
!                                                        ! 
!                     Read wrf data                      !
!                                                        !
!--------------------------------------------------------!

  !data prthres/1.0,5.0,10.,20.,50./
  prthres=9999.0
  
  inquire(file='./wrfhist1.nc',exist=alive1)
  inquire(file='./wrfhist2.nc',exist=alive2)
  inquire(file='./st4obs.bin',exist=alive3)  
  if(alive1.and.alive2.and.alive3) then
    filepath(1)="wrfhist1.nc"
	filepath(2)="wrfhist2.nc"
	st4path="st4obs.nc"
  else  
    print*,"need data files,stop"
    stop
  endif
  !print('("read in file:",A)'),trim(filepath)

  open(100,file="namelist.input")
  read(100,nml=options)
  !print*,nthres,prthres,dx,dy,radii,timlab,unbiasopt
  close(100)
  
  allocate(vartem(1,1,1))

  do i=1,nvar
    !print*,trim(varname(i))
    scanopt=1
    call  readnc(filepath(2),varname(i),dim(i,1),dim(i,2),dim(i,3),vartem,scanopt)
    allocate(var3d(i)%val(dim(i,1),dim(i,2),dim(i,3)))
    scanopt=0
    call  readnc(filepath(2),varname(i),dim(i,1),dim(i,2),dim(i,3),var3d(i)%val,scanopt)
	allocate(fcst(dim(1,1),dim(1,2),dim(1,3)))
	fcst=var3d(i)%val
	call  readnc(filepath(1),varname(i),dim(i,1),dim(i,2),dim(i,3),var3d(i)%val,scanopt)
	fcst=fcst-var3d(i)%val
  enddo 
  allocate(st4m(dim(1,1),dim(1,2),dim(1,3)))

!--------------------------------------------------------!
!                                                        ! 
!                     Read obs data                      !
!                                                        !
!--------------------------------------------------------!

  !print*,"read interpolated obs"
  open(100,file='st4obs.bin',form='unformatted')
  read(100) st4m
  close(100)
  !print*,"ST4 max,min:",maxval(maxval(st4m(:,:,1),1),1),minval(minval(st4m(:,:,1),1),1)

!--------------------------------------------------------!
!                                                        ! 
!                     Compute FSS                        !
!                                                        !
!--------------------------------------------------------!  
  !timlab="201806200400"
  do i=1,nthres
    call ensfss(dim(1,1),dim(1,2),1,fcst(:,:,1),st4m(:,:,1),  &
                dx,dy,radii,prthres(i),fss,timlab,unbiasopt)
  
  enddo
  
  
  
  
  !print*,"Finished"
   
end program st4vswrf_fss


subroutine ensfss(nx,ny,nen,data2d,st4obs,  &
                  dx,dy,radii,thres,fss,timlab,unbiasopt)
                     
implicit none

integer :: nx,ny,nen,unbiasopt
integer :: i,j,inen,itime,ii,jj,kk,iii,jjj,kkk
integer :: io,jo
real    :: dx,dy

real    :: data2d(nx,ny,nen)
real    :: st4obs(nx,ny)
integer :: temi01,temi02,temi03,temi04,temi05,temi06
real    :: temr01,temr02,temr03,temr04,temr05,temr06, &
           temr07,temr08,temr09,temr10
logical :: condition,condition1,condition2,condition3
integer :: icounta(nen)
real    :: fss(nen),fss_worst(nen)
real    :: radii,thres
integer :: is,ie,js,je
real    :: fbsm(nen),fbso(nen)
character(len=*) :: timlab
character(len=80) :: outfmt
real    :: rcount

real :: fcst_pctl(nx,ny,nen)
real :: obs_pctl(nx,ny)
real :: fss0d,fss0upper,fss0lower
real :: fss_accum(nen)
real :: fss_upper(nen)
real :: fss_lower(nen)

!print*,"entering ensfss"
fss=0
fss_worst=0
icounta=0
fss_accum=0
fss_upper=0
fss_lower=0

fcst_pctl=-1
obs_pctl=-1

if(unbiasopt==1.and.1==0) then
 do inen=1,nen

  temi01=1
  temi02=1
  temi03=0
  temi04=0
  do j=1,ny
   do i=1,nx
     temr01=-999999
     temr02=-999999   

     do jj=1,ny
       do ii=1,nx
         if(data2d(ii,jj,inen)>temr01.and.fcst_pctl(ii,jj,inen)<0) then
           temr01=data2d(ii,jj,inen)
           iii=ii
           jjj=jj  
         endif
       enddo
     enddo
     if(data2d(iii,jjj,inen)>1) temi03=temi03+1
     fcst_pctl(iii,jjj,inen)=1.0*temi01
     print*,'fsct:',iii,jjj,data2d(iii,jjj,inen), fcst_pctl(iii,jjj,inen)
     temi01=temi01+1
     
     if(inen==1) then
      do jj=1,ny
        do ii=1,nx
          if(st4obs(ii,jj)>temr02.and.obs_pctl(ii,jj)<0) then
            temr02=st4obs(ii,jj)
            iii=ii
            jjj=jj
          endif
         enddo
      enddo
      if(st4obs(iii,jjj)>=1) temi04=temi04+1
      obs_pctl(iii,jjj)=1.0*temi02
      temi02=temi02+1
      print*,'obs:',iii,jjj,st4obs(iii,jjj),obs_pctl(iii,jjj)
     endif

   enddo
  enddo
  if(temi03>=1) then
    fcst_pctl(:,:,inen)=fcst_pctl(:,:,inen)/temi03
  else
    fcst_pctl(:,:,inen)=0
  endif
  if(inen==1) then
    if(temi04>=1) then
      obs_pctl(:,:)=obs_pctl(:,:)/temi03
    else
      obs_pctl(:,:)=0
    endif
  endif
  print*,temi03,temi04
  !stop
 enddo
endif
if(unbiasopt==1.and.1==1) then
  obs_pctl(:,:)=st4obs(:,:)/maxval(maxval(st4obs,1),1)
  fcst_pctl(:,:,1)=data2d(:,:,1)/maxval(maxval(data2d(:,:,1),1),1)
endif

do jo=1,ny
  do io=1,nx
    !print*,jo,io
    fbsm=0
    fbso=0
    do inen=1,nen
	
	  if(unbiasopt==0)then
        condition=st4obs(io,jo)>=thres .or. data2d(io,jo,inen)>=thres
	  else if(unbiasopt==1) then 
        condition=obs_pctl(io,jo)>=thres .or. fcst_pctl(io,jo,inen)>=thres
	  endif
	  
      if(condition) then
        !print*,"inen,io,jo,st4obs,data2d",inen,io,jo,st4obs(io,jo),data2d(io,jo,inen)
      	!print*,"enter condition", jo,io
        is=max(io-int(radii/dx),1)
        ie=min(io+int(radii/dx),nx)
        js=max(jo-int(radii/dy),1)
        je=min(jo+int(radii/dy),ny)
      
        temr02=0
        temr03=0
        temi01=0

        do j=js,je
         do i=is,ie
            if(unbiasopt==1) then		 
              if(fcst_pctl(i,j,inen)>=thres) temr02=temr02+1
              if(obs_pctl(i,j)>=thres) temr03=temr03+1
			else if(unbiasopt==0) then
              if(data2d(i,j,inen)>=thres) temr02=temr02+1
              if(st4obs(i,j)>=thres) temr03=temr03+1			
			endif
            temi01=temi01+1
          enddo  
        enddo

        temr02=temr02/temi01 !*100
        temr03=temr03/temi01 !*100

  !      print*,temr02,temr03 
        fss(inen)=fss(inen)+ (temr02-temr03)**2 
        fss_worst(inen)=fss_worst(inen)+temr02**2+temr03**2  

        if(temr02**2+temr03**2>0) then
          fss0d=1-(temr02-temr03)**2/(temr02**2+temr03**2)
          fss0upper=(temr02-temr03)**2
          fss0lower=(temr02**2+temr03**2)
        else
          fss0d=0
        endif 
        fss_accum(inen)=fss_accum(inen)+fss0d
        fss_upper(inen)=fss_upper(inen)+fss0upper
        fss_lower(inen)=fss_lower(inen)+fss0lower
        !print*,'mem:',inen,fss(inen),fss_worst(inen),fss0d,fss_accum(inen)
        icounta(inen)=icounta(inen)+1
       
      endif        
    enddo
    
  enddo
enddo                     

temr01=0   
rcount=0
do inen=1,nen
  if(icounta(inen)>0) then    
!    fss(inen)=fss_accum(inen)/icounta(inen)  !1-fss(inen)/fss_worst(inen)
    fss(inen)=1-fss_upper(inen)/fss_lower(inen) 
 
    if(fss(inen)>=0.5) then
      rcount=rcount+1
    endif
  else
    fss(inen)=0
  endif  
  temr01=temr01+fss(inen)/nen

enddo
rcount=rcount/nen
write(outfmt,'(A,I2.2,A)') "(A,X,A,X,F8.3,",1,'F12.3)'

write(*,trim(outfmt)),"FSS:",timlab,thres,fss(1)

	                
end subroutine ensfss  



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
     !print*,iend
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
   !print*,"max,min:",maxval(maxval(maxval(vartem,1),1),1),minval(minval(minval(vartem,1),1),1)


   var=vartem
   deallocate(vartem)

end subroutine readnc





