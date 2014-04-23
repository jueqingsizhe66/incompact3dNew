subroutine body(ux,uy,uz,ep1)
!
!*******************************************************************

USE param
USE decomp_2d
USE variables
USE IBM

implicit none

real(mytype), dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) :: ux,uy,uz,ep1
integer :: i,j,k
real(mytype) :: xm,ym,r

ep1=0.
do k=xstart(3),xend(3)
do j=xstart(2),xend(2)
do i=xstart(1),xend(1)
   xm=(i-1)*dx
   ym=(j-1)*dy
   r=sqrt((xm-cex)*(xm-cex)+(ym-cey)*(ym-cey))
   if (r-ra >= 0.) cycle
   ux(i,j,k)=0.
   uy(i,j,k)=0.
   uz(i,j,k)=0.
   ep1(i,j,k)=1.
enddo
enddo
enddo


return
end subroutine body
