subroutine prepare (b,c,f,s,w,n)
!
!*******************************************************************

use decomp_2d, only : mytype

implicit none

integer :: i,n
real(mytype), dimension(n) :: b,c,f,s,w

do i=1,n
   w(i)=c(i)
enddo
do i=2,n
   s(i)=b(i-1)/w(i-1)
   w(i)=w(i)-f(i-1)*s(i)
enddo
do i=1,n
   w(i)=1.0_mytype/w(i)
enddo

return
end subroutine prepare
