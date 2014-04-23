subroutine inflow (ux,uy,uz,phi)
!
!*********************************************************

USE param
USE IBM
USE variables
USE decomp_2d

implicit none

integer  :: k,j
real,dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz,phi
real :: r1,r2,r3,y,um

call ecoule(ux,uy,uz)

call random_number(bxo)
call random_number(byo)
call random_number(bzo)

if (iin.eq.1) then
   do k=1,xsize(3)
   do j=1,xsize(2)
      bxx1(j,k)=bxx1(j,k)+bxo(j,k)*noise1
      bxy1(j,k)=bxy1(j,k)+byo(j,k)*noise1
      bxz1(j,k)=bxz1(j,k)+bzo(j,k)*noise1
   enddo
   enddo
   if (iscalar==1) then
      do k=1,xsize(3)
      do j=1,xsize(2)
         phi(1,j,k)=1.
      enddo
      enddo
   endif
endif

return
end subroutine inflow
