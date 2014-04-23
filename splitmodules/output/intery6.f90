subroutine intery6(ty,uy,ry,sy,cify6,cisy6,ciwy6,nx,ny,nym,nz,npaire)
!
!********************************************************************

USE param
USE derivY

implicit none

integer :: nx,ny,nym,nz,npaire
real(mytype), dimension(nx,nym,nz) :: ty
real(mytype), dimension(nx,ny,nz) :: uy,ry
real(mytype), dimension(nx,nz) :: sy
real(mytype), dimension(nym) :: cify6,cisy6,ciwy6
integer :: i,j,k

if (ncly==0) then
   do k=1,nz
   do i=1,nx
      ty(i,1,k)=aiciy6*(uy(i,2,k)+uy(i,1,k))&
           +biciy6*(uy(i,3,k)+uy(i,ny,k))&
           +ciciy6*(uy(i,4,k)+uy(i,ny-1,k))
      ry(i,1,k)=-1.
      ty(i,2,k)=aiciy6*(uy(i,3,k)+uy(i,2,k))&
           +biciy6*(uy(i,4,k)+uy(i,1,k))&
           +ciciy6*(uy(i,5,k)+uy(i,ny,k))
      ry(i,2,k)=0.
   enddo
   enddo
   do k=1,nz
   do j=3,ny-3
   do i=1,nx
      ty(i,j,k)=aiciy6*(uy(i,j+1,k)+uy(i,j,k))&
           +biciy6*(uy(i,j+2,k)+uy(i,j-1,k))&
           +ciciy6*(uy(i,j+3,k)+uy(i,j-2,k))
      ry(i,j,k)=0.
   enddo
   enddo
   enddo
   do k=1,nz
   do i=1,nx
      ty(i,ny-2,k)=aiciy6*(uy(i,ny-1,k)+uy(i,ny-2,k))&
           +biciy6*(uy(i,ny,k)+uy(i,ny-3,k))&
           +ciciy6*(uy(i,1,k)+uy(i,ny-4,k))
      ry(i,ny-2,k)=0.
      ty(i,ny-1,k)=aiciy6*(uy(i,ny,k)+uy(i,ny-1,k))&
           +biciy6*(uy(i,1,k)+uy(i,ny-2,k))&
           +ciciy6*(uy(i,2,k)+uy(i,ny-3,k))
      ry(i,ny-1,k)=0.
      ty(i,ny,k)=aiciy6*(uy(i,1,k)+uy(i,ny,k))&
           +biciy6*(uy(i,2,k)+uy(i,ny-1,k))&
           +ciciy6*(uy(i,3,k)+uy(i,ny-2,k))
      ry(i,ny,k)=ailcaiy6
   enddo
   enddo
   do k=1,nz
   do j=2,ny
   do i=1,nx
      ty(i,j,k)=ty(i,j,k)-ty(i,j-1,k)*cisy6(j)
      ry(i,j,k)=ry(i,j,k)-ry(i,j-1,k)*cisy6(j)
   enddo
   enddo
   enddo
   do k=1,nz
   do i=1,nx
      ty(i,ny,k)=ty(i,ny,k)*ciwy6(ny)
      ry(i,ny,k)=ry(i,ny,k)*ciwy6(ny)
   enddo
   enddo
   do k=1,nz
   do j=ny-1,1,-1
   do i=1,nx
      ty(i,j,k)=(ty(i,j,k)-cify6(j)*ty(i,j+1,k))*ciwy6(j)
      ry(i,j,k)=(ry(i,j,k)-cify6(j)*ry(i,j+1,k))*ciwy6(j)
   enddo
   enddo
   enddo
   do k=1,nz
   do i=1,nx
      sy(i,k)=(ty(i,1,k)-ailcaiy6*ty(i,ny,k))&
           /(1.+ry(i,1,k)-ailcaiy6*ry(i,ny,k))
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
   if (npaire==1) then
      do k=1,nz
      do i=1,nx
         ty(i,1,k)=aiciy6*(uy(i,2,k)+uy(i,1,k))&
              +biciy6*(uy(i,3,k)+uy(i,2,k))&
              +ciciy6*(uy(i,4,k)+uy(i,3,k))
         ty(i,2,k)=aiciy6*(uy(i,3,k)+uy(i,2,k))&
              +biciy6*(uy(i,4,k)+uy(i,1,k))&
              +ciciy6*(uy(i,5,k)+uy(i,2,k))
      enddo
      enddo
      do k=1,nz
      do j=3,nym-2
      do i=1,nx
         ty(i,j,k)=aiciy6*(uy(i,j+1,k)+uy(i,j,k))&
              +biciy6*(uy(i,j+2,k)+uy(i,j-1,k))&
              +ciciy6*(uy(i,j+3,k)+uy(i,j-2,k))
      enddo
      enddo
      enddo
      do k=1,nz
      do i=1,nx
         ty(i,nym-1,k)=aiciy6*(uy(i,nym,k)+uy(i,nym-1,k))&
              +biciy6*(uy(i,ny,k)+uy(i,nym-2,k))&
              +ciciy6*(uy(i,nym,k)+uy(i,nym-3,k))
         ty(i,nym,k)=aiciy6*(uy(i,ny,k)+uy(i,nym,k))&
              +biciy6*(uy(i,nym,k)+uy(i,nym-1,k))&
              +ciciy6*(uy(i,nym-1,k)+uy(i,nym-2,k))
      enddo
      enddo
      do k=1,nz
      do j=2,nym
      do i=1,nx
         ty(i,j,k)=ty(i,j,k)-ty(i,j-1,k)*cisy6(j)
      enddo
      enddo
      enddo
      do k=1,nz
      do i=1,nx
         ty(i,nym,k)=ty(i,nym,k)*ciwy6(nym)
      enddo
      enddo
      do k=1,nz
      do j=nym-1,1,-1
      do i=1,nx
         ty(i,j,k)=(ty(i,j,k)-cify6(j)*ty(i,j+1,k))*ciwy6(j)
      enddo
      enddo
      enddo
   endif
endif

return
end subroutine intery6
