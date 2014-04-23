subroutine derz(tz,uz,rz,sz,ffz,fsz,fwz,nx,ny,nz,npaire)
!
!********************************************************************

USE param
USE derivZ

implicit none

integer :: nx,ny,nz,npaire,i,j,k
real(mytype), dimension(nx,ny,nz) :: tz,uz,rz
real(mytype), dimension(nx,ny) :: sz
real(mytype), dimension(nz) :: ffz,fsz,fwz

if (nclz==0) then
   do j=1,ny
   do i=1,nx
      tz(i,j,1)=afkz*(uz(i,j,2)-uz(i,j,nz  ))&
           +bfkz*(uz(i,j,3)-uz(i,j,nz-1))
      rz(i,j,1)=-1.
      tz(i,j,2)=afkz*(uz(i,j,3)-uz(i,j,1 ))&
           +bfkz*(uz(i,j,4)-uz(i,j,nz))
      rz(i,j,2)=0.
   enddo
   enddo
   do k=3,nz-2
   do j=1,ny
   do i=1,nx
      tz(i,j,k)=afkz*(uz(i,j,k+1)-uz(i,j,k-1))&
               +bfkz*(uz(i,j,k+2)-uz(i,j,k-2))
      rz(i,j,k)=0.
   enddo
   enddo
   enddo
   do j=1,ny
   do i=1,nx
      tz(i,j,nz-1)=afkz*(uz(i,j,nz)-uz(i,j,nz-2))&
           +bfkz*(uz(i,j,1 )-uz(i,j,nz-3))
      rz(i,j,nz-1)=0.
      tz(i,j,nz  )=afkz*(uz(i,j,1)-uz(i,j,nz-1))&
           +bfkz*(uz(i,j,2)-uz(i,j,nz-2))
      rz(i,j,nz  )=alfakz
   enddo
   enddo
   do k=2,nz
   do j=1,ny
   do i=1,nx
      tz(i,j,k)=tz(i,j,k)-tz(i,j,k-1)*fsz(k)
      rz(i,j,k)=rz(i,j,k)-rz(i,j,k-1)*fsz(k)
   enddo
   enddo
   enddo
   do j=1,ny
   do i=1,nx
      tz(i,j,nz)=tz(i,j,nz)*fwz(nz)
      rz(i,j,nz)=rz(i,j,nz)*fwz(nz)
   enddo
   enddo
   do k=nz-1,1,-1
   do j=1,ny
   do i=1,nx
      tz(i,j,k)=(tz(i,j,k)-ffz(k)*tz(i,j,k+1))*fwz(k)
      rz(i,j,k)=(rz(i,j,k)-ffz(k)*rz(i,j,k+1))*fwz(k)
   enddo
   enddo
   enddo
   do j=1,ny
   do i=1,nx
      sz(i,j)=(   tz(i,j,1)-alfakz*tz(i,j,nz))/&
           (1.+rz(i,j,1)-alfakz*rz(i,j,nz))
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
         tz(i,j,1)=0.
         tz(i,j,2)=afkz*(uz(i,j,3)-uz(i,j,1))&
              +bfkz*(uz(i,j,4)-uz(i,j,2))
      enddo
      enddo
      do k=3,nz-2
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=afkz*(uz(i,j,k+1)-uz(i,j,k-1))&
              +bfkz*(uz(i,j,k+2)-uz(i,j,k-2))
      enddo
      enddo
      enddo
      do j=1,ny
      do i=1,nx
         tz(i,j,nz-1)=afkz*(uz(i,j,nz  )-uz(i,j,nz-2))&
              +bfkz*(uz(i,j,nz-1)-uz(i,j,nz-3))
         tz(i,j,nz  )=0.
      enddo
      enddo
      do k=2,nz
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=tz(i,j,k)-tz(i,j,k-1)*fsz(k)
      enddo
      enddo
      enddo
      do j=1,ny
      do i=1,nx
         tz(i,j,nz)=tz(i,j,nz)*fwz(nz)
      enddo
      enddo
      do k=nz-1,1,-1
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=(tz(i,j,k)-ffz(k)*tz(i,j,k+1))*fwz(k)
      enddo
      enddo
      enddo
   endif
   if (npaire==0) then
      do j=1,ny
      do i=1,nx
         tz(i,j,1)=afkz*(uz(i,j,2)+uz(i,j,2))&
              +bfkz*(uz(i,j,3)+uz(i,j,3))
         tz(i,j,2)=afkz*(uz(i,j,3)-uz(i,j,1))&
              +bfkz*(uz(i,j,4)+uz(i,j,2))
      enddo
      enddo
      do k=3,nz-2
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=afkz*(uz(i,j,k+1)-uz(i,j,k-1))&
              +bfkz*(uz(i,j,k+2)-uz(i,j,k-2))
      enddo
      enddo
      enddo
      do j=1,ny
      do i=1,nx
         tz(i,j,nz-1)=afkz*( uz(i,j,nz  )-uz(i,j,nz-2))&
              +bfkz*(-uz(i,j,nz-1)-uz(i,j,nz-3))
         tz(i,j,nz  )=afkz*(-uz(i,j,nz-1)-uz(i,j,nz-1))&
              +bfkz*(-uz(i,j,nz-2)-uz(i,j,nz-2))
      enddo
      enddo
      do k=2,nz
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=tz(i,j,k)-tz(i,j,k-1)*fsz(k)
      enddo
      enddo
      enddo
      do j=1,ny
      do i=1,nx
         tz(i,j,nz)=tz(i,j,nz)*fwz(nz)
      enddo
      enddo
      do k=nz-1,1,-1
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=(tz(i,j,k)-ffz(k)*tz(i,j,k+1))*fwz(k)
      enddo
      enddo
      enddo
   endif
endif

if (nclz==2) then
   do j=1,ny
   do i=1,nx
      tz(i,j,1)=af1z*uz(i,j,1)+bf1z*uz(i,j,2)&
           +cf1z*uz(i,j,3)
      tz(i,j,2)=af2z*(uz(i,j,3)-uz(i,j,1))
   enddo
   enddo
   do k=3,nz-2
   do j=1,ny
   do i=1,nx
      tz(i,j,k)=afkz*(uz(i,j,k+1)-uz(i,j,k-1))&
           +bfkz*(uz(i,j,k+2)-uz(i,j,k-2))
   enddo
   enddo
   enddo
   do j=1,ny
   do i=1,nx
      tz(i,j,nz-1)= afmz*(uz(i,j,nz)-uz(i,j,nz-2))
      tz(i,j,nz  )=-afnz*uz(i,j,nz)-bfnz*uz(i,j,nz-1)&
           -cfnz*uz(i,j,nz-2)
   enddo
   enddo
   do k=2,nz
   do j=1,ny
   do i=1,nx
      tz(i,j,k)=tz(i,j,k)-tz(i,j,k-1)*fsz(k)
   enddo
   enddo
   enddo
   do j=1,ny
   do i=1,nx
      tz(i,j,nz)=tz(i,j,nz)*fwz(nz)
   enddo
   enddo
   do k=nz-1,1,-1
   do j=1,ny
   do i=1,nx
      tz(i,j,k)=(tz(i,j,k)-ffz(k)*tz(i,j,k+1))*fwz(k)
   enddo
   enddo
   enddo
endif

return
end subroutine derz
