subroutine decy6(ty,uy,ry,sy,cfy6,csy6,cwy6,ppyi,nx,ny,nym,nz,npaire)
!
!********************************************************************

USE param
USE derivY

implicit none

integer :: nx,ny,nym,nz,npaire
real(mytype), dimension(nx,nym,nz) :: ty
real(mytype), dimension(nx,ny,nz) :: uy
real(mytype), dimension(nx,ny,nz) :: ry
real(mytype), dimension(nx,nz) :: sy
real(mytype), dimension(nym) :: cfy6,csy6,cwy6,ppyi
integer :: i,j,k

if (ncly==0) then
   do k=1,nz
   do i=1,nx
      ty(i,1,k)=aciy6*(uy(i,2,k)-uy(i,1,k))&
           +bciy6*(uy(i,3,k)-uy(i,ny,k))
      ry(i,1,k)=-1.
      ty(i,2,k)=aciy6*(uy(i,3,k)-uy(i,2,k))&
           +bciy6*(uy(i,4,k)-uy(i,1,k))
      ry(i,2,k)=0.
   enddo
   enddo
   do k=1,nz
   do j=3,ny-2
   do i=1,nx
      ty(i,j,k)=aciy6*(uy(i,j+1,k)-uy(i,j,k))&
           +bciy6*(uy(i,j+2,k)-uy(i,j-1,k))
      ry(i,j,k)=0.
   enddo
   enddo
   enddo
   do k=1,nz
   do i=1,nx
      ty(i,ny-1,k)=aciy6*(uy(i,ny,k)-uy(i,ny-1,k))&
           +bciy6*(uy(i,1,k)-uy(i,ny-2,k))
      ry(i,ny-1,k)=0.
      ty(i,ny,k)=aciy6*(uy(i,1,k)-uy(i,ny,k))&
           +bciy6*(uy(i,2,k)-uy(i,ny-1,k))
      ry(i,ny,k)=alcaiy6
   enddo
   enddo
   do k=1,nz
   do j=2,ny
   do i=1,nx
      ty(i,j,k)=ty(i,j,k)-ty(i,j-1,k)*csy6(j)
      ry(i,j,k)=ry(i,j,k)-ry(i,j-1,k)*csy6(j)
   enddo
   enddo
   enddo
   do k=1,nz
   do i=1,nx
      ty(i,ny,k)=ty(i,ny,k)*cwy6(ny)
      ry(i,ny,k)=ry(i,ny,k)*cwy6(ny)
   enddo
   enddo
   do k=1,nz
   do j=ny-1,1,-1
   do i=1,nx
      ty(i,j,k)=(ty(i,j,k)-cfy6(j)*ty(i,j+1,k))*cwy6(j)
      ry(i,j,k)=(ry(i,j,k)-cfy6(j)*ry(i,j+1,k))*cwy6(j)
   enddo
   enddo
   enddo
   do k=1,nz
   do i=1,nx
      sy(i,k)=(ty(i,1,k)-alcaiy6*ty(i,ny,k))&
           /(1.+ry(i,1,k)-alcaiy6*ry(i,ny,k))
   enddo
   enddo
   do k=1,nz
   do j=1,ny
   do i=1,nx
      ty(i,j,k)=ty(i,j,k)-sy(i,k)*ry(i,j,k)
   enddo
   enddo
   enddo
endif
if ((ncly==1).or.(ncly==2)) then
   if (npaire==0) then
      do k=1,nz
      do i=1,nx
         ty(i,1,k)=aciy6*(uy(i,2,k)-uy(i,1,k))&
              +bciy6*(uy(i,3,k)-2.*uy(i,1,k)+uy(i,2,k))
         ty(i,2,k)=aciy6*(uy(i,3,k)-uy(i,2,k))&
              +bciy6*(uy(i,4,k)-uy(i,1,k))
      enddo
      enddo
      do k=1,nz
      do j=3,nym-2
      do i=1,nx
         ty(i,j,k)=aciy6*(uy(i,j+1,k)-uy(i,j,k))&
              +bciy6*(uy(i,j+2,k)-uy(i,j-1,k))
      enddo
      enddo
      enddo
      do k=1,nz
      do i=1,nx
         ty(i,nym-1,k)=aciy6*(uy(i,nym,k)-uy(i,nym-1,k))&
              +bciy6*(uy(i,ny,k)-uy(i,nym-2,k))
         ty(i,nym,k)=aciy6*(uy(i,ny,k)-uy(i,nym,k))&
              +bciy6*(2.*uy(i,ny,k)-uy(i,nym,k)-uy(i,nym-1,k))
      enddo
      enddo
      do k=1,nz
      do j=2,nym
      do i=1,nx
         ty(i,j,k)=ty(i,j,k)-ty(i,j-1,k)*csy6(j)
      enddo
      enddo
      enddo
      do k=1,nz
      do i=1,nx
         ty(i,nym,k)=ty(i,nym,k)*cwy6(nym)
      enddo
      enddo
      do k=1,nz
      do j=nym-1,1,-1
      do i=1,nx
         ty(i,j,k)=(ty(i,j,k)-cfy6(j)*ty(i,j+1,k))*cwy6(j)
      enddo
      enddo
      enddo
   endif
endif

if (istret.ne.0) then
   do k=1,nz
   do j=1,nym
   do i=1,nx
      ty(i,j,k)=ty(i,j,k)*ppyi(j)
   enddo
   enddo
   enddo
endif

return
end subroutine decy6
