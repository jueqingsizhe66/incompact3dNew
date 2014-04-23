subroutine decz6(tz,uz,rz,sz,cfz6,csz6,cwz6,nx,ny,nz,nzm,npaire)
!
!********************************************************************

USE param
USE derivZ

implicit none

integer :: nx,ny,nz,nzm,npaire
real(mytype), dimension(nx,ny,nzm) :: tz
real(mytype), dimension(nx,ny,nz) :: uz
real(mytype), dimension(nx,ny,nz) :: rz
real(mytype), dimension(nx,ny) :: sz
real(mytype), dimension(nzm) :: cfz6,csz6,cwz6
integer :: i,j,k

if (nclz==0) then
   do j=1,ny
   do i=1,nx
      tz(i,j,1)=aciz6*(uz(i,j,2)-uz(i,j,1))&
           +bciz6*(uz(i,j,3)-uz(i,j,nz))
      rz(i,j,1)=-1.
      tz(i,j,2)=aciz6*(uz(i,j,3)-uz(i,j,2))&
           +bciz6*(uz(i,j,4)-uz(i,j,1))
      rz(i,j,2)=0.
   enddo
   enddo
   do k=3,nz-2
   do j=1,ny
   do i=1,nx
      tz(i,j,k)=aciz6*(uz(i,j,k+1)-uz(i,j,k))&
           +bciz6*(uz(i,j,k+2)-uz(i,j,k-1))
      rz(i,j,k)=0.
   enddo
   enddo
   enddo
   do j=1,ny
   do i=1,nx
      tz(i,j,nz-1)=aciz6*(uz(i,j,nz)-uz(i,j,nz-1))&
           +bciz6*(uz(i,j,1)-uz(i,j,nz-2))
      rz(i,j,nz-1)=0.
      tz(i,j,nz)=aciz6*(uz(i,j,1)-uz(i,j,nz))&
           +bciz6*(uz(i,j,2)-uz(i,j,nz-1))
      rz(i  ,j,nz)=alcaiz6
   enddo
   enddo
   do k=2,nz
   do j=1,ny
   do i=1,nx
      tz(i,j,k)=tz(i,j,k)-tz(i,j,k-1)*csz6(k)
      rz(i,j,k)=rz(i,j,k)-rz(i,j,k-1)*csz6(k)
   enddo
   enddo
   enddo
   do i=1,nx
   do j=1,ny
      tz(i,j,nz)=tz(i,j,nz)*cwz6(nz)
      rz(i,j,nz)=rz(i,j,nz)*cwz6(nz)
   enddo
   enddo
   do k=nz-1,1,-1
   do j=1,ny
   do i=1,nx
      tz(i,j,k)=(tz(i,j,k)-cfz6(k)*tz(i,j,k+1))*cwz6(k)
      rz(i,j,k)=(rz(i,j,k)-cfz6(k)*rz(i,j,k+1))*cwz6(k)
   enddo
   enddo
   enddo
   do j=1,ny
   do i=1,nx
      sz(i,j)=(tz(i,j,1)-alcaiz6*tz(i,j,nz))/&
           (1.+rz(i,j,1)-alcaiz6*rz(i,j,nz))
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

if ((nclz==1).or.(nclz==2)) then
   if (npaire==1) then
      do j=1,ny
      do i=1,nx
         tz(i,j,1)=aciz6*(uz(i,j,2)-uz(i,j,1))&
              +bciz6*(uz(i,j,3)-uz(i,j,2))
         tz(i,j,2)=aciz6*(uz(i,j,3)-uz(i,j,2))&
              +bciz6*(uz(i,j,4)-uz(i,j,1))
      enddo
      enddo
      do k=3,nzm-2
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=aciz6*(uz(i,j,k+1)-uz(i,j,k))&
              +bciz6*(uz(i,j,k+2)-uz(i,j,k-1))
      enddo
      enddo
      enddo
      do j=1,ny
      do i=1,nx
         tz(i,j,nzm-1)=aciz6*(uz(i,j,nzm)-uz(i,j,nzm-1))&
              +bciz6*(uz(nz,j,k)-uz(nzm-2,j,k))
         tz(i,j,nzm)=aciz6*(uz(i,j,nz)-uz(i,j,nzm))&
              +bciz6*(uz(i,j,nzm)-uz(i,j,nzm-1))
      enddo
      enddo
      do k=2,nzm
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=tz(i,j,k)-tz(i,j,k-1)*csz6(k)
      enddo
      enddo
      enddo
      do j=1,ny
      do i=1,nx
         tz(i,j,nzm)=tz(i,j,nzm)*cwz6(nzm)
      enddo
      enddo
      do k=nzm-1,1,-1
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=(tz(i,j,k)-cfz6(k)*tz(i,j,k+1))*cwz6(k)
      enddo
      enddo
      enddo
   endif
   if (npaire==0) then
      do j=1,ny
      do i=1,nx
         tz(i,j,1)=aciz6*(uz(i,j,2)-uz(i,j,1))&
              +bciz6*(uz(i,j,3)-2.*uz(i,j,1)+uz(i,j,2))
         tz(i,j,2)=aciz6*(uz(i,j,3)-uz(i,j,2))&
              +bciz6*(uz(i,j,4)-uz(i,j,1))
      enddo
      enddo
      do k=3,nzm-2
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=aciz6*(uz(i,j,k+1)-uz(i,j,k))&
              +bciz6*(uz(i,j,k+2)-uz(i,j,k-1))
      enddo
      enddo
      enddo
      do j=1,ny
      do i=1,nx
         tz(i,j,nzm-1)=aciz6*(uz(i,j,nz-1)-uz(i,j,nz-2))&
              +bciz6*(uz(i,j,nz)-uz(i,j,nz-3))
         tz(i,j,nzm)=aciz6*(uz(i,j,nz)-uz(i,j,nz-1))&
              +bciz6*(2.*uz(i,j,nz)-uz(i,j,nz-1)-uz(i,j,nz-2))
      enddo
      enddo
      do k=2,nzm
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=tz(i,j,k)-tz(i,j,k-1)*csz6(k)
      enddo
      enddo
      enddo
      do j=1,ny
      do i=1,nx
         tz(i,j,nzm)=tz(i,j,nzm)*cwz6(nzm)
      enddo
      enddo
      do k=nzm-1,1,-1
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=(tz(i,j,k)-cfz6(k)*tz(i,j,k+1))*cwz6(k)
      enddo
      enddo
      enddo
   endif
endif

return
end subroutine decz6
