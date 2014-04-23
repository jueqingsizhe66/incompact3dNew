subroutine deciy6(ty,uy,ry,sy,cfi6y,csi6y,cwi6y,cfy6,csy6,cwy6,&
     ppy,nx,nym,ny,nz,npaire)
!
!********************************************************************

USE param
USE derivY

implicit none

integer :: nx,ny,nym,nz,npaire
real(mytype), dimension(nx,ny,nz) :: ty
real(mytype), dimension(nx,nym,nz) :: uy
real(mytype), dimension(nx,ny,nz) :: ry
real(mytype), dimension(nx,nz) :: sy
real(mytype), dimension(ny) :: cfi6y,csi6y,cwi6y,ppy
real(mytype), dimension(nym) :: cfy6,csy6,cwy6
integer :: i,j,k

if (ncly==0) then
   do k=1,nz
   do i=1,nx
      ty(i,1,k)=aciy6*(uy(i,1,k)-uy(i,ny,k))&
           +bciy6*(uy(i,2,k)-uy(i,ny-1,k))
      ry(i,1,k)=-1.
      ty(i,2,k)=aciy6*(uy(i,2,k)-uy(i,1,k))&
           +bciy6*(uy(i,3,k)-uy(i,ny,k))
      ry(i,2,k)=0.
   enddo
   enddo
   do j=3,ny-2
   do k=1,nz
   do i=1,nx
      ty(i,j,k)=aciy6*(uy(i,j,k)-uy(i,j-1,k))&
           +bciy6*(uy(i,j+1,k)-uy(i,j-2,k))
      ry(i,j,k)=0.
   enddo
   enddo
   enddo
   do k=1,nz
   do i=1,nx
      ty(i,ny-1,k)=aciy6*(uy(i,ny-1,k)-uy(i,ny-2,k))&
           +bciy6*(uy(i,ny,k)-uy(i,ny-3,k))
      ry(i,ny-1,k)=0.
      ty(i,ny,k)=aciy6*(uy(i,ny,k)-uy(i,ny-1,k))&
           +bciy6*(uy(i,1,k)-uy(i,ny-2,k))
      ry(i,ny,k)=alcaiy6
   enddo
   enddo
   do j=2,ny
   do k=1,nz
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
   do j=ny-1,1,-1
   do k=1,nz
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
   do j=1,ny
   do k=1,nz
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
         ty(i,1,k)=0.
         ty(i,2,k)=aciy6*(uy(i,2,k)-uy(i,1,k))&
              +bciy6*(uy(i,3,k)-uy(i,1,k))
      enddo
      enddo
      do j=3,ny-2
      do k=1,nz
      do i=1,nx
         ty(i,j,k)=aciy6*(uy(i,j,k)-uy(i,j-1,k))&
              +bciy6*(uy(i,j+1,k)-uy(i,j-2,k))
      enddo
      enddo
      enddo
      do k=1,nz
      do i=1,nx
         ty(i,ny-1,k)=aciy6*(uy(i,ny-1,k)-uy(i,ny-2,k))&
              +bciy6*(uy(i,ny-1,k)-uy(i,ny-3,k))
         ty(i,ny,k)=0.
      enddo
      enddo
      do j=2,ny
      do k=1,nz
      do i=1,nx
         ty(i,j,k)=ty(i,j,k)-ty(i,j-1,k)*csi6y(j)
      enddo
      enddo
      enddo
      do k=1,nz
      do i=1,nx
         ty(i,ny,k)=ty(i,ny,k)*cwi6y(ny)
      enddo
      enddo
      do j=ny-1,1,-1
      do k=1,nz
      do i=1,nx
         ty(i,j,k)=(ty(i,j,k)-cfi6y(j)*ty(i,j+1,k))*cwi6y(j)
      enddo
      enddo
      enddo
   endif
endif

if (istret.ne.0) then
   do k=1,nz
   do j=1,ny
   do i=1,nx
      ty(i,j,k)=ty(i,j,k)*ppy(j)
   enddo
   enddo
   enddo
endif

return
end subroutine deciy6
