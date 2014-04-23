subroutine init (ux1,uy1,uz1,phi1,gx1,gy1,gz1,phis1,hx1,hy1,hz1,phiss1)
!
!********************************************************************

USE decomp_2d
USE param
USE IBM
USE variables

implicit none

real,dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,phi1
real,dimension(xsize(1),xsize(2),xsize(3)) :: gx1,gy1,gz1,phis1
real,dimension(xsize(1),xsize(2),xsize(3)) :: hx1,hy1,hz1,phiss1

real :: y,r,um,r1,r2,r3
integer :: k,j,i

if (iin.eq.1) then !generation of a random noise

   call random_number(ux1)
   call random_number(uy1)
   call random_number(uz1)

   do k=1,xsize(3)
   do j=1,xsize(2)
   do i=1,xsize(1)
      ux1(i,j,k)=noise*ux1(i,j,k)
      uy1(i,j,k)=noise*uy1(i,j,k)
      uz1(i,j,k)=noise*uz1(i,j,k)
   enddo
   enddo
   enddo

   if (iscalar==1) then
      do k=1,xsize(3)
      do j=1,xsize(2)
      do i=1,xsize(1)
         phi1(i,j,k)=1.
         phis1(i,j,k)=phi1(i,j,k)
         phiss1(i,j,k)=phis1(i,j,k)
      enddo
      enddo
      enddo
   endif



endif

if (iin.eq.2) then !read a correlated noise generated before
endif

!MEAN FLOW PROFILE
call ecoule(ux1,uy1,uz1)
!INIT FOR G AND U=MEAN FLOW + NOISE
do k=1,xsize(3)
do j=1,xsize(2)
do i=1,xsize(1)
   ux1(i,j,k)=ux1(i,j,k)+bxx1(j,k)
   uy1(i,j,k)=uy1(i,j,k)+bxy1(j,k)
   uz1(i,j,k)=uz1(i,j,k)+bxz1(j,k)
   gx1(i,j,k)=ux1(i,j,k)
   gy1(i,j,k)=uy1(i,j,k)
   gz1(i,j,k)=uz1(i,j,k)
   hx1(i,j,k)=gx1(i,j,k)
   hy1(i,j,k)=gy1(i,j,k)
   hz1(i,j,k)=gz1(i,j,k)
enddo
enddo
enddo

return
end subroutine init
