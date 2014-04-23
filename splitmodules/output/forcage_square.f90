subroutine forcage_square(ux,uy,uz,esp)
!
!*******************************************************************

USE param
USE decomp_2d
USE variables
USE IBM


implicit none

real(mytype), dimension(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)) :: ux,uy,uz,esp

integer :: j, i, k, np,i1
real(mytype) :: ep0,ep1,ep2,ep3

esp(:,:,:)=0.

ep2=0.4 !2.1=1
ep1=1.2
ep0=3.2

!####fisrt fractal iteration##########################
 call square(57.6,57.6,ep0,28.8,ux,uy,uz,esp)
!####Second fractal iteration#########################
 call square(86.4,86.4,ep1,14.4,ux,uy,uz,esp)
 call square(86.4,28.8,ep1,14.4,ux,uy,uz,esp)
 call square(28.8,86.4,ep1,14.4,ux,uy,uz,esp)
 call square(28.8,28.8,ep1,14.4,ux,uy,uz,esp)
!####third fractal iteration########################
 call square(43.2,43.2,ep2,7.2,ux,uy,uz,esp)
 call square(43.2,14.4,ep2,7.2,ux,uy,uz,esp)
 call square(14.4,14.4,ep2,7.2,ux,uy,uz,esp)
 call square(14.4,43.2,ep2,7.2,ux,uy,uz,esp)

 call square(100.8,100.8,ep2,7.2,ux,uy,uz,esp)
 call square(100.8,72.,ep2,7.2,ux,uy,uz,esp)
 call square(72.,72.,ep2,7.2,ux,uy,uz,esp)
 call square(72.,100.8,ep2,7.2,ux,uy,uz,esp)

 call square(100.8,14.4,ep2,7.2,ux,uy,uz,esp)
 call square(72.,14.4,ep2,7.2,ux,uy,uz,esp)
 call square(100.8,43.2,ep2,7.2,ux,uy,uz,esp)
 call square(72.,43.2,ep2,7.2,ux,uy,uz,esp)

 call square(43.2,72.,ep2,7.2,ux,uy,uz,esp)
 call square(43.2,100.8,ep2,7.2,ux,uy,uz,esp)
 call square(14.4,72.,ep2,7.2,ux,uy,uz,esp)
 call square(14.4,100.8,ep2,7.2,ux,uy,uz,esp)


return
end subroutine forcage_square
