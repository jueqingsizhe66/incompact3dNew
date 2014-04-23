subroutine square(ycenter,zcenter,xthick,xlenght,ux,uy,uz,esp)
!
!####################################################

USE param
USE decomp_2d
USE variables
USE IBM


implicit none

real(mytype), dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) :: ux,uy,uz,esp
real(mytype) :: ycenter,zcenter,xthick,xlenght,slenght,xz1,xz2,xy1,xy2
integer :: j1,j2,z1,z2,iep,i,j,k,k1,k2,ilen

iep=int(xthick*ny/yly)
ilen=int(xlenght*ny/yly)

j1=int(ycenter*ny/yly)
k1=int(zcenter*nz/zlz)+ilen
k2=int(zcenter*nz/zlz)-ilen

do k=xstart(3),xend(3)
do j=xstart(2),xend(2)
do i=xstart(1),xend(1)
   if ((k.ge.(k2-iep)).and.(k.le.(k1+iep))) then
   if ((j.ge.(j1+ilen-iep)).and.(j.le.(j1+ilen+iep))) then
   if ((i.ge.25).and.(i.le.36)) then
      esp(i,j,k)=1.
      ux(i,j,k)=0.
      uy(i,j,k)=0.
      uz(i,j,k)=0.
   endif
   endif
   endif
   if ((k.ge.(k2-iep)).and.(k.le.(k1+iep))) then
   if ((j.ge.(j1-ilen-iep)).and.(j.le.(j1-ilen+iep))) then
   if ((i.ge.25).and.(i.le.36)) then
      esp(i,j,k)=1.
      ux(i,j,k)=0.
      uy(i,j,k)=0.
      uz(i,j,k)=0.
   endif
   endif
   endif
   if ((k.ge.(k1-iep)).and.(k.le.(k1+iep))) then
   if ((j.ge.(j1-ilen-iep)).and.(j.le.(j1+ilen+iep))) then
   if ((i.ge.25).and.(i.le.36)) then
      esp(i,j,k)=1.
      ux(i,j,k)=0.
      uy(i,j,k)=0.
      uz(i,j,k)=0.
   endif
   endif
   endif
   if ((k.ge.(k2-iep)).and.(k.le.(k2+iep))) then
   if ((j.ge.(j1-ilen-iep)).and.(j.le.(j1+ilen+iep))) then
   if ((i.ge.25).and.(i.le.36)) then
      esp(i,j,k)=1.
      ux(i,j,k)=0.
      uy(i,j,k)=0.
      uz(i,j,k)=0.
   endif
   endif
   endif
enddo
enddo
enddo


return
end subroutine square
