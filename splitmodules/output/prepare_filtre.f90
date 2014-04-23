subroutine prepare_filtre(aff,af,a,ab,abb,al,n)
!
!*********************************************************************

use decomp_2d, only : mytype

implicit none

integer :: n ,i,j,k,l
real(mytype), dimension(n) :: aff,af,a,ab,abb
real(mytype), dimension(n,2) :: al
real(mytype) :: tiny, dum
!*********************************************************************

tiny=1.E-10



aff(1)=a(1)
af(1)=ab(1)
a(1)=abb(1)
ab(1)=0.
abb(1)=0.

aff(2)=af(2)
af(2)=a(2)
a(2)=ab(2)
ab(2)=abb(2)
abb(2)=0.



l=2

do k=1,n
   dum=aff(k)
   i=k
   if (l<n) l=l+1
   do j=k+1,l
      if (abs(aff(j)) <= abs(dum)) cycle
      dum=aff(j)
      i=j
   enddo
   if (dum==0.) aff(k)=0.
   if (i/=k) then
      dum=aff(k)
      aff(k)=aff(i)
      aff(i)=dum
      dum=af(k)
      af(k)=af(i)
      af(i)=dum
      dum=a(k)
      a(k)=a(i)
      a(i)=dum
      dum=ab(k)
      ab(k)=ab(i)
      ab(i)=dum
      dum=abb(k)
      abb(k)=abb(i)
      abb(i)=dum
   endif
   do i=k+1,l
      dum=aff(i)/aff(k)
      al(k,i-k)=dum
      aff(i)=af(i)-dum*af(k)
      af(i)=a(i)-dum*a(k)
      a(i)=ab(i)-dum*ab(k)
      ab(i)=abb(i)-dum*abb(k)
      abb(i)=0.
   enddo
enddo

return
end subroutine prepare_filtre
