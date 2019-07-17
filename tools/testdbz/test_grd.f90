program test_grd

implicit none

integer,parameter :: nx=3,ny=3
real,dimension(nx,ny) :: u,v,y
real,dimension(nx,ny) :: grdu,grdv
real    :: ds=1000
integer :: gropt

integer :: i,j
real    :: grd_lt,grd_rt
character(len=80)  :: wrtfmt


data u/7,2,2,4,2,4,3,2,8/
data v/4,6,3,2,3,2,6,4,5/

y=0
grdu=0
grdv=0
gropt=0
call fdiv(u,v,y,nx,ny,ds,gropt)
write(wrtfmt,'(A2,I1,A7)') "(",nx,"F16.8)"
print*,"u"
do j=ny,1,-1
  print(trim(wrtfmt)),(u(i,j),i=1,nx)
enddo
print*,"v"
do j=ny,1,-1
  print(trim(wrtfmt)),(v(i,j),i=1,nx)
enddo
print*,"Equation: Y=dU/dx+dV/dy"
print*,"dx=dy=ds=",ds
print*,"y"
do j=ny,1,-1
  print(trim(wrtfmt)),(y(i,j),i=1,nx)
enddo

gropt=1
call fdiv(grdu,grdv,y,nx,ny,ds,gropt)
print*,"grdu"
do j=ny,1,-1
  print(trim(wrtfmt)),(grdu(i,j),i=1,nx)
enddo
print*,"grdv"
do j=ny,1,-1
  print(trim(wrtfmt)),(grdv(i,j),i=1,nx)
enddo

grd_lt=0
grd_rt=0
do j=1,ny
 do i=1,nx
    grd_lt=grd_lt+y(i,j)*y(i,j)
    grd_rt=grd_rt+u(i,j)*grdu(i,j)+v(i,j)*grdv(i,j)
 enddo
enddo
print*,"------------------------------------------"
print*,"<y,y>,<x,A*y>",grd_lt,grd_rt
print*,"------------------------------------------"

end program test_grd


subroutine fdiv(u,v,y,nx,ny,ds,gropt)

implicit none

integer               :: gropt
integer               :: nx,ny
real,dimension(nx,ny) :: u,v,y
real                  :: ds

integer :: i,j

if(gropt==1) then
 u=0
 v=0
endif

! non boundary area
do j=2,ny-1
  do i=2,nx-1  

    if(gropt==0) then
      y(i,j)=0.5*(u(i+1,j)-u(i-1,j))/ds+0.5*(v(i,j+1)-v(i,j-1))/ds
    endif
    if(gropt==1) then
      u(i+1,j)=u(i+1,j)+0.5*y(i,j)/ds
      u(i-1,j)=u(i-1,j)-0.5*y(i,j)/ds
      v(i,j+1)=v(i,j+1)+0.5*y(i,j)/ds
      v(i,j-1)=v(i,j-1)-0.5*y(i,j)/ds	  
    endif

  enddo  
enddo

! boundary area, exclude corners
j=1
do i=2,nx-1
  if(gropt==0) then
    y(i,j)=0.5*(u(i+1,j)-u(i-1,j))/ds+(v(i,j+1)-v(i,j))/ds
  endif
  if(gropt==1) then
    u(i+1,j)=u(i+1,j)+0.5*y(i,j)/ds
    u(i-1,j)=u(i-1,j)-0.5*y(i,j)/ds
    v(i,j+1)=v(i,j+1)+y(i,j)/ds
    v(i,j)  =v(i,j)-y(i,j)/ds
  endif  
enddo

j=ny
do i=2,nx-1
  if(gropt==0) then
    y(i,j)=0.5*(u(i+1,j)-u(i-1,j))/ds+(v(i,j)-v(i,j-1))/ds
  endif
  if(gropt==1) then
    u(i+1,j)=u(i+1,j)+0.5*y(i,j)/ds
    u(i-1,j)=u(i-1,j)-0.5*y(i,j)/ds
    v(i,j)  =v(i,j)+y(i,j)/ds
    v(i,j-1)=v(i,j-1)-y(i,j)/ds
  endif  
enddo


i=1
do j=2,ny-1
  if(gropt==0) then
    y(i,j)=(u(i+1,j)-u(i,j))/ds+0.5*(v(i,j+1)-v(i,j-1))/ds
  endif
  if(gropt==1) then
    u(i,j)  =u(i,j)-y(i,j)/ds
    u(i+1,j)=u(i+1,j)+y(i,j)/ds
    v(i,j+1)=v(i,j+1)+0.5*y(i,j)/ds
    v(i,j-1)=v(i,j-1)-0.5*y(i,j)/ds
  endif  
enddo

i=nx
do j=2,ny-1
  if(gropt==0) then
    y(i,j)=(u(i,j)-u(i-1,j))/ds+0.5*(v(i,j+1)-v(i,j-1))/ds
  endif
  if(gropt==1) then
    u(i,j)  =u(i,j)+y(i,j)/ds
    u(i-1,j)=u(i-1,j)-y(i,j)/ds
    v(i,j+1)=v(i,j+1)+0.5*y(i,j)/ds
    v(i,j-1)=v(i,j-1)-0.5*y(i,j)/ds
  endif  
enddo

! Four corners

i=1
j=1
if(gropt==0) then
  y(i,j)=(u(i+1,j)-u(i,j))/ds+(v(i,j+1)-v(i,j))/ds
endif
if(gropt==1) then
  u(i+1,j)=u(i+1,j)+y(i,j)/ds
  u(i,j)=u(i,j)-y(i,j)/ds
  v(i,j+1)=v(i,j+1)+y(i,j)/ds
  v(i,j)=v(i,j)-y(i,j)/ds
endif

i=1
j=ny
if(gropt==0) then
  y(i,j)=(u(i+1,j)-u(i,j))/ds+(v(i,j)-v(i,j-1))/ds
endif
if(gropt==1) then
  u(i+1,j)=u(i+1,j)+y(i,j)/ds
  u(i,j)=u(i,j)-y(i,j)/ds
  v(i,j)=v(i,j)+y(i,j)/ds
  v(i,j-1)=v(i,j-1)-y(i,j)/ds
endif


i=nx
j=1
if(gropt==0) then
  y(i,j)=(u(i,j)-u(i-1,j))/ds+(v(i,j+1)-v(i,j))/ds
endif
if(gropt==1) then
  u(i,j)=u(i,j)+y(i,j)/ds
  u(i-1,j)=u(i-1,j)-y(i,j)/ds
  v(i,j+1)=v(i,j+1)+y(i,j)/ds
  v(i,j)=v(i,j)-y(i,j)/ds
endif

i=nx
j=ny
if(gropt==0) then
  y(i,j)=(u(i,j)-u(i-1,j))/ds+(v(i,j)-v(i,j-1))/ds
endif
if(gropt==1) then
  u(i,j)=u(i,j)+y(i,j)/ds
  u(i-1,j)=u(i-1,j)-y(i,j)/ds
  v(i,j)=v(i,j)+y(i,j)/ds
  v(i,j-1)=v(i,j-1)-y(i,j)/ds
endif


end subroutine fdiv
