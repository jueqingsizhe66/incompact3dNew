subroutine deciz6(tz,uz,rz,sz,cfiz6,csiz6,cwiz6,cfz6,csz6,cwz6,&
     nx,ny,nzm,nz,npaire)
!
!********************************************************************

USE param
USE derivZ

implicit none

integer :: nx,nzm,ny,nz,npaire
real(mytype), dimension(nx,ny,nz) :: tz
real(mytype), dimension(nx,ny,nzm) :: uz,rz
real(mytype), dimension(nx,ny) :: sz
real(mytype), dimension(nz) :: cfiz6,csiz6,cwiz6
real(mytype), dimension(nz) :: cfz6,csz6,cwz6
integer :: i,j,k

if (nclz==0) then
   do j=1,ny
   do i=1,nx
      tz(i,j,1)=aciz6*(uz(i,j,1)-uz(i,j,nz))&
           +bciz6*(uz(i,j,2)-uz(i,j,nz-1))
      rz(i,j,1)=-1.
      tz(i,j,2)=aciz6*(uz(i,j,2)-uz(i,j,1))&
           +bciz6*(uz(i,j,3)-uz(i,j,nz))
      rz(i,j,2)=0.
   enddo
   enddo
   do k=3,nz-2
   do j=1,ny
   do i=1,nx
      tz(i,j,k)=aciz6*(uz(i,j,k)-uz(i,j,k-1))&
           +bciz6*(uz(i,j,k+1)-uz(i,j,k-2))
      rz(i,j,k)=0.
   enddo
   enddo
   enddo
   do j=1,ny
   do i=1,nx
      tz(i,j,nz-1)=aciz6*(uz(i,j,nz-1)-uz(i,j,nz-2))&
           +bciz6*(uz(i,j,nz)-uz(i,j,nz-3))
      rz(i,j,nz-1)=0.
      tz(i,j,nz)=aciz6*(uz(i,j,nz)-uz(i,j,nz-1))&
           +bciz6*(uz(i,j,1)-uz(i,j,nz-2))
      rz(i,j,nz)=alcaiz6
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
   do j=1,ny
   do i=1,nx
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
         tz(i,j,1)=0.
         tz(i,j,2)=aciz6*(uz(i,j,2)-uz(i,j,1))&
              +bciz6*(uz(i,j,3)-uz(i,j,1))
      enddo
      enddo
      do k=3,nz-2
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=aciz6*(uz(i,j,k)-uz(i,j,k-1))&
              +bciz6*(uz(i,j,k+1)-uz(i,j,k-2))
      enddo
      enddo
      enddo
      do j=1,ny
      do i=1,nx
         tz(i,j,nz-1)=aciz6*(uz(i,j,nz-1)-uz(i,j,nz-2))&
              +bciz6*(uz(i,j,nz-1)-uz(i,j,nz-3))
         tz(i,j,nz)=0.
      enddo
      enddo
      do k=2,nz
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=tz(i,j,k)-tz(i,j,k-1)*csiz6(k)
      enddo
      enddo
      enddo
      do j=1,ny
      do i=1,nx
         tz(i,j,nz)=tz(i,j,nz)*cwiz6(nz)
      enddo
      enddo
      do k=nz-1,1,-1
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=(tz(i,j,k)-cfiz6(k)*tz(i,j,k+1))*cwiz6(k)
      enddo
      enddo
      enddo
   endif
endif

return
end subroutine deciz6
