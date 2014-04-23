subroutine derzz(tz,uz,rz,sz,sfz,ssz,swz,nx,ny,nz,npaire)
!
!********************************************************************

USE param
USE derivZ

implicit none

integer :: nx,ny,nz,npaire,i,j,k
real(mytype), dimension(nx,ny,nz) :: tz,uz,rz
real(mytype), dimension(nx,ny) :: sz
real(mytype), dimension(nz) :: sfz,ssz,swz

if (nclz==0) then
   do j=1,ny
   do i=1,nx
      tz(i,j,1)=askz*(uz(i,j,2)-uz(i,j,1   )&
           -uz(i,j,1)+uz(i,j,nz  ))&
           +bskz*(uz(i,j,3)-uz(i,j,1   )&
           -uz(i,j,1)+uz(i,j,nz-1))&
           +cskz*(uz(i,j,4)-uz(i,j,1   )&
           -uz(i,j,1)+uz(i,j,nz-2))
      rz(i,j,1)=-1.
      tz(i,j,2)=askz*(uz(i,j,3)-uz(i,j,2 )&
           -uz(i,j,2)+uz(i,j,1 ))&
           +bskz*(uz(i,j,4)-uz(i,j,2 )&
           -uz(i,j,2)+uz(i,j,nz))&
           +cskz*(uz(i,j,5)-uz(i,j,2 )&
           -uz(i,j,2)+uz(i,j,nz-1))
      rz(i,j,2)=0.
      tz(i,j,3)=askz*(uz(i,j,4)-uz(i,j,3 )&
           -uz(i,j,3)+uz(i,j,2 ))&
           +bskz*(uz(i,j,5)-uz(i,j,3 )&
           -uz(i,j,3)+uz(i,j,1 ))&
           +cskz*(uz(i,j,6)-uz(i,j,3 )&
           -uz(i,j,3)+uz(i,j,nz))
      rz(i,j,3)=0.
   enddo
   enddo
   do k=4,nz-3
   do j=1,ny
   do i=1,nx
      tz(i,j,k)=askz*(uz(i,j,k+1)-uz(i,j,k  )&
           -uz(i,j,k  )+uz(i,j,k-1))&
           +bskz*(uz(i,j,k+2)-uz(i,j,k  )&
           -uz(i,j,k  )+uz(i,j,k-2))&
           +cskz*(uz(i,j,k+3)-uz(i,j,k  )&
           -uz(i,j,k  )+uz(i,j,k-3))
      rz(i,j,k)=0.
   enddo
   enddo
   enddo
   do j=1,ny
   do i=1,nx
      tz(i,j,nz-2)=askz*(uz(i,j,nz-1)-uz(i,j,nz-2)&
           -uz(i,j,nz-2)+uz(i,j,nz-3))&
           +bskz*(uz(i,j,nz  )-uz(i,j,nz-2)&
           -uz(i,j,nz-2)+uz(i,j,nz-4))&
           +cskz*(uz(i,j,1   )-uz(i,j,nz-2)&
           -uz(i,j,nz-2)+uz(i,j,nz-5))
      rz(i,j,nz-2)=0.
      tz(i,j,nz-1)=askz*(uz(i,j,nz  )-uz(i,j,nz-1)&
           -uz(i,j,nz-1)+uz(i,j,nz-2))&
           +bskz*(uz(i,j,1   )-uz(i,j,nz-1)&
           -uz(i,j,nz-1)+uz(i,j,nz-3))&
           +cskz*(uz(i,j,2   )-uz(i,j,nz-1)&
           -uz(i,j,nz-1)+uz(i,j,nz-4))
      rz(i,j,nz-1)=0.
      tz(i,j,nz  )=askz*(uz(i,j,1 )-uz(i,j,nz  )&
           -uz(i,j,nz)+uz(i,j,nz-1))&
           +bskz*(uz(i,j,2 )-uz(i,j,nz  )&
           -uz(i,j,nz)+uz(i,j,nz-2))&
           +cskz*(uz(i,j,3 )-uz(i,j,nz  )&
           -uz(i,j,nz)+uz(i,j,nz-3))
      rz(i,j,nz  )=alsakz
   enddo
   enddo
   do k=2,nz
   do j=1,ny
   do i=1,nx
      tz(i,j,k)=tz(i,j,k)-tz(i,j,k-1)*ssz(k)
      rz(i,j,k)=rz(i,j,k)-rz(i,j,k-1)*ssz(k)
   enddo
   enddo
   enddo
   do j=1,ny
   do i=1,nx
      tz(i,j,nz)=tz(i,j,nz)*swz(nz)
      rz(i,j,nz)=rz(i,j,nz)*swz(nz)
   enddo
   enddo
   do k=nz-1,1,-1
   do j=1,ny
   do i=1,nx
      tz(i,j,k)=(tz(i,j,k)-sfz(k)*tz(i,j,k+1))*swz(k)
      rz(i,j,k)=(rz(i,j,k)-sfz(k)*rz(i,j,k+1))*swz(k)
   enddo
   enddo
   enddo
   do j=1,ny
   do i=1,nx
      sz(i,j)=(   tz(i,j,1)-alsakz*tz(i,j,nz))/&
           (1.+rz(i,j,1)-alsakz*rz(i,j,nz))
   enddo
   enddo
   do k=1,nz
   do j=1,ny
   do i=1,nx
      tz(i,j,k)=tz(i,j,k)-sz(i,j)*rz(i,j,k)
   enddo
   enddo
   enddo
endif

if (nclz==1) then
   if (npaire==1) then
      do j=1,ny
      do i=1,nx
         tz(i,j,1)=askz*(uz(i,j,2)-uz(i,j,1)&
              -uz(i,j,1)+uz(i,j,2))&
              +bskz*(uz(i,j,3)-uz(i,j,1)&
              -uz(i,j,1)+uz(i,j,3))&
              +cskz*(uz(i,j,4)-uz(i,j,1)&
              -uz(i,j,1)+uz(i,j,4))
         tz(i,j,2)=askz*(uz(i,j,3)-uz(i,j,2)&
              -uz(i,j,2)+uz(i,j,1))&
              +bskz*(uz(i,j,4)-uz(i,j,2)&
              -uz(i,j,2)+uz(i,j,2))&
              +cskz*(uz(i,j,5)-uz(i,j,2)&
              -uz(i,j,2)+uz(i,j,3))
         tz(i,j,3)=askz*(uz(i,j,4)-uz(i,j,3)&
              -uz(i,j,3)+uz(i,j,2))&
              +bskz*(uz(i,j,5)-uz(i,j,3)&
              -uz(i,j,3)+uz(i,j,1))&
              +cskz*(uz(i,j,6)-uz(i,j,3)&
              -uz(i,j,3)+uz(i,j,2))
      enddo
      enddo
      do k=4,nz-3
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=askz*(uz(i,j,k+1)-uz(i,j,k  )&
              -uz(i,j,k  )+uz(i,j,k-1))&
              +bskz*(uz(i,j,k+2)-uz(i,j,k  )&
              -uz(i,j,k  )+uz(i,j,k-2))&
              +cskz*(uz(i,j,k+3)-uz(i,j,k  )&
              -uz(i,j,k  )+uz(i,j,k-3))
      enddo
      enddo
      enddo
      do j=1,ny
      do i=1,nx
         tz(i,j,nz-2)=askz*(uz(i,j,nz-1)-uz(i,j,nz-2)&
              -uz(i,j,nz-2)+uz(i,j,nz-3))&
              +bskz*(uz(i,j,nz  )-uz(i,j,nz-2)&
              -uz(i,j,nz-2)+uz(i,j,nz-4))&
              +cskz*(uz(i,j,nz-1)-uz(i,j,nz-2)&
              -uz(i,j,nz-2)+uz(i,j,nz-5))
         tz(i,j,nz-1)=askz*(uz(i,j,nz  )-uz(i,j,nz-1)&
              -uz(i,j,nz-1)+uz(i,j,nz-2))&
              +bskz*(uz(i,j,nz-1)-uz(i,j,nz-1)&
              -uz(i,j,nz-1)+uz(i,j,nz-3))&
              +cskz*(uz(i,j,nz-2)-uz(i,j,nz-1)&
              -uz(i,j,nz-1)+uz(i,j,nz-4))
         tz(i,j,nz  )=askz*(uz(i,j,nz-1)-uz(i,j,nz  )&
              -uz(i,j,nz  )+uz(i,j,nz-1))&
              +bskz*(uz(i,j,nz-2)-uz(i,j,nz  )&
              -uz(i,j,nz  )+uz(i,j,nz-2))&
              +cskz*(uz(i,j,nz-3)-uz(i,j,nz  )&
              -uz(i,j,nz  )+uz(i,j,nz-3))
      enddo
      enddo
      do k=2,nz
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=tz(i,j,k)-tz(i,j,k-1)*ssz(k)
      enddo
      enddo
      enddo
      do j=1,ny
      do i=1,nx
         tz(i,j,nz)=tz(i,j,nz)*swz(nz)
      enddo
      enddo
      do k=nz-1,1,-1
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=(tz(i,j,k)-sfz(k)*tz(i,j,k+1))*swz(k)
      enddo
      enddo
      enddo
   endif
   if (npaire==0) then
      do j=1,ny
      do i=1,nx
         tz(i,j,1)=0.
         tz(i,j,2)=askz*(uz(i,j,3)-uz(i,j,2)&
              -uz(i,j,2)+uz(i,j,1))&
              +bskz*(uz(i,j,4)-uz(i,j,2)&
              -uz(i,j,2)-uz(i,j,2))&
              +cskz*(uz(i,j,5)-uz(i,j,2)&
              -uz(i,j,2)-uz(i,j,3))
         tz(i,j,3)=askz*(uz(i,j,4)-uz(i,j,3)&
              -uz(i,j,3)+uz(i,j,2))&
              +bskz*(uz(i,j,5)-uz(i,j,3)&
              -uz(i,j,3)+uz(i,j,1))&
              +cskz*(uz(i,j,6)-uz(i,j,3)&
              -uz(i,j,3)-uz(i,j,2))
      enddo
      enddo
      do k=4,nz-3
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=askz*(uz(i,j,k+1)-uz(i,j,k  )&
              -uz(i,j,k  )+uz(i,j,k-1))&
              +bskz*(uz(i,j,k+2)-uz(i,j,k  )&
              -uz(i,j,k  )+uz(i,j,k-2))&
              +cskz*(uz(i,j,k+3)-uz(i,j,k  )&
              -uz(i,j,k  )+uz(i,j,k-3))
      enddo
      enddo
      enddo
      do j=1,ny
      do i=1,nx
         tz(i,j,nz-2)=askz*( uz(i,j,nz-1)-uz(i,j,nz-2)&
              -uz(i,j,nz-2)+uz(i,j,nz-3))&
              +bskz*( uz(i,j,nz  )-uz(i,j,nz-2)&
              -uz(i,j,nz-2)+uz(i,j,nz-4))&
              +cskz*(-uz(i,j,nz-1)-uz(i,j,nz-2)&
              -uz(i,j,nz-2)+uz(i,j,nz-5))
         tz(i,j,nz-1)=askz*( uz(i,j,nz  )-uz(i,j,nz-1)&
              -uz(i,j,nz-1)+uz(i,j,nz-2))&
              +bskz*(-uz(i,j,nz-1)-uz(i,j,nz-1)&
              -uz(i,j,nz-1)+uz(i,j,nz-3))&
              +cskz*(-uz(i,j,nz-2)-uz(i,j,nz-1)&
              -uz(i,j,nz-1)+uz(i,j,nz-4))
         tz(i,j,nz  )=0.
      enddo
      enddo
      do k=2,nz
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=tz(i,j,k)-tz(i,j,k-1)*ssz(k)
      enddo
      enddo
      enddo
      do j=1,ny
      do i=1,nx
         tz(i,j,nz)=tz(i,j,nz)*swz(nz)
      enddo
      enddo
      do k=nz-1,1,-1
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=(tz(i,j,k)-sfz(k)*tz(i,j,k+1))*swz(k)
      enddo
      enddo
      enddo
   endif
endif

if (nclz==2) then
   do j=1,ny
   do i=1,nx
      tz(i,j,1)=as1z*uz(i,j,1)+bs1z*uz(i,j,2)&
           +cs1z*uz(i,j,3)+ds1z*uz(i,j,4)
      tz(i,j,2)=as2z*(uz(i,j,3)-uz(i,j,2)&
           -uz(i,j,2)+uz(i,j,1))
      tz(i,j,3)=as3z*(uz(i,j,4)-uz(i,j,3)&
           -uz(i,j,3)+uz(i,j,2))&
           +bs3z*(uz(i,j,5)-uz(i,j,3)&
           -uz(i,j,3)+uz(i,j,1))
   enddo
   enddo
   do k=4,nz-3
   do j=1,ny
   do i=1,nx
      tz(i,j,k)=askz*(uz(i,j,k+1)-uz(i,j,k  )&
           -uz(i,j,k  )+uz(i,j,k-1))&
           +bskz*(uz(i,j,k+2)-uz(i,j,k  )&
           -uz(i,j,k  )+uz(i,j,k-2))&
           +cskz*(uz(i,j,k+3)-uz(i,j,k  )&
           -uz(i,j,k  )+uz(i,j,k-3))
   enddo
   enddo
   enddo
   do j=1,ny
   do i=1,nx
      tz(i,j,nz-2)=astz*(uz(i,j,nz-1)-uz(i,j,nz-2)&
           -uz(i,j,nz-2)+uz(i,j,nz-3))&
           +bstz*(uz(i,j,nz  )-uz(i,j,nz-2)&
           -uz(i,j,nz-2)+uz(i,j,nz-4))
      tz(i,j,nz-1)=asmz*(uz(i,j,nz  )-uz(i,j,nz-1)&
           -uz(i,j,nz-1)+uz(i,j,nz-2))
      tz(i,j,nz  )=asnz*uz(i,j,nz  )+bsnz*uz(i,j,nz-1)&
           +csnz*uz(i,j,nz-2)+dsnz*uz(i,j,nz-3)
   enddo
   enddo
   do k=2,nz
   do j=1,ny
   do i=1,nx
      tz(i,j,k)=tz(i,j,k)-tz(i,j,k-1)*ssz(k)
   enddo
   enddo
   enddo
   do j=1,ny
   do i=1,nx
      tz(i,j,nz)=tz(i,j,nz)*swz(nz)
   enddo
   enddo
   do k=nz-1,1,-1
   do j=1,ny
   do i=1,nx
      tz(i,j,k)=(tz(i,j,k)-sfz(k)*tz(i,j,k+1))*swz(k)
   enddo
   enddo
   enddo
endif

return
end subroutine derzz
