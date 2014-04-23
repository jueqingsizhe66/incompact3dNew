subroutine dery(ty,uy,ry,sy,ffy,fsy,fwy,ppy,nx,ny,nz,npaire)
!
!********************************************************************

USE param
USE derivY

implicit none

integer :: nx,ny,nz,i,j,k,npaire
real(mytype), dimension(nx,ny,nz) :: ty,uy
real(mytype), dimension(nx,ny,nz) :: ry
real(mytype), dimension(nx,nz)  :: sy
real(mytype), dimension(ny) :: ffy,fsy,fwy,ppy

if (ncly==0) then
   do k=1,nz
   do i=1,nx
      ty(i,1,k)=afjy*(uy(i,2,k)-uy(i,ny,k))&
           +bfjy*(uy(i,3,k)-uy(i,ny-1,k))
      ry(i,1,k)=-1.
      ty(i,2,k)=afjy*(uy(i,3,k)-uy(i,1,k))&
           +bfjy*(uy(i,4,k)-uy(i,ny,k))
      ry(i,2,k)=0.
   enddo
   enddo
   do k=1,nz
   do j=3,ny-2
   do i=1,nx
      ty(i,j,k)=afjy*(uy(i,j+1,k)-uy(i,j-1,k))&
           +bfjy*(uy(i,j+2,k)-uy(i,j-2,k))
      ry(i,j,k)=0.
   enddo
   enddo
   enddo
   do k=1,nz
   do i=1,nx
      ty(i,ny-1,k)=afjy*(uy(i,ny,k)-uy(i,ny-2,k))&
           +bfjy*(uy(i,1,k)-uy(i,ny-3,k))
      ry(i,ny-1,k)=0.
      ty(i,ny,k)=afjy*(uy(i,1,k)-uy(i,ny-1,k))&
           +bfjy*(uy(i,2,k)-uy(i,ny-2,k))
      ry(i,ny,k)=alfajy
   enddo
   enddo
   do k=1,nz
   do j=2,ny
   do i=1,nx
      ty(i,j,k)=ty(i,j,k)-ty(i,j-1,k)*fsy(j)
      ry(i,j,k)=ry(i,j,k)-ry(i,j-1,k)*fsy(j)
   enddo
   enddo
   enddo
   do k=1,nz
   do i=1,nx
      ty(i,ny,k)=ty(i,ny,k)*fwy(ny)
      ry(i,ny,k)=ry(i,ny,k)*fwy(ny)
   enddo
   enddo
   do k=1,nz
   do j=ny-1,1,-1
   do i=1,nx
      ty(i,j,k)=(ty(i,j,k)-ffy(j)*ty(i,j+1,k))*fwy(j)
      ry(i,j,k)=(ry(i,j,k)-ffy(j)*ry(i,j+1,k))*fwy(j)
   enddo
   enddo
   enddo
   do k=1,nz
   do i=1,nx
      sy(i,k)=(ty(i,1,k)-alfajy*ty(i,ny,k))&
           /(1.+ry(i,1,k)-alfajy*ry(i,ny,k))
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

if (ncly==1) then
   if (npaire==1) then
      do k=1,nz
      do i=1,nx
         ty(i,1,k)=0.
         ty(i,2,k)=afjy*(uy(i,3,k)-uy(i,1,k))&
              +bfjy*(uy(i,4,k)-uy(i,2,k))
      enddo
      enddo
      do k=1,nz
      do j=3,ny-2
      do i=1,nx
         ty(i,j,k)=afjy*(uy(i,j+1,k)-uy(i,j-1,k))&
              +bfjy*(uy(i,j+2,k)-uy(i,j-2,k))
      enddo
      enddo
      enddo
      do k=1,nz
      do i=1,nx
         ty(i,ny-1,k)=afjy*(uy(i,ny,k)-uy(i,ny-2,k))&
              +bfjy*(uy(i,ny-1,k)-uy(i,ny-3,k))
         ty(i,ny,k)=0.
      enddo
      enddo
      do k=1,nz
      do j=2,ny
      do i=1,nx
         ty(i,j,k)=ty(i,j,k)-ty(i,j-1,k)*fsy(j)
      enddo
      enddo
      enddo
      do k=1,nz
      do i=1,nx
         ty(i,ny,k)=ty(i,ny,k)*fwy(ny)
      enddo
      enddo
      do k=1,nz
      do j=ny-1,1,-1
      do i=1,nx
         ty(i,j,k)=(ty(i,j,k)-ffy(j)*ty(i,j+1,k))*fwy(j)
      enddo
      enddo
      enddo
   endif
   if (npaire==0) then
      do k=1,nz
      do i=1,nx
         ty(i,1,k)=afjy*(uy(i,2,k)+uy(i,2,k))&
              +bfjy*(uy(i,3,k)+uy(i,3,k))
         ty(i,2,k)=afjy*(uy(i,3,k)-uy(i,1,k))&
              +bfjy*(uy(i,4,k)+uy(i,2,k))
      enddo
      enddo
      do k=1,nz
      do j=3,ny-2
      do i=1,nx
         ty(i,j,k)=afjy*(uy(i,j+1,k)-uy(i,j-1,k))&
              +bfjy*(uy(i,j+2,k)-uy(i,j-2,k))
      enddo
      enddo
      enddo
      do k=1,nz
      do i=1,nx
         ty(i,ny-1,k)=afjy*(uy(i,ny,k)-uy(i,ny-2,k))&
              +bfjy*((-uy(i,ny-1,k))-uy(i,ny-3,k))
         ty(i,ny,k)=afjy*((-uy(i,ny-1,k))-uy(i,ny-1,k))&
              +bfjy*((-uy(i,ny-2,k))-uy(i,ny-2,k))
      enddo
      enddo
      do k=1,nz
      do j=2,ny
      do i=1,nx
         ty(i,j,k)=ty(i,j,k)-ty(i,j-1,k)*fsy(j)
      enddo
      enddo
      enddo
      do k=1,nz
      do i=1,nx
         ty(i,ny,k)=ty(i,ny,k)*fwy(ny)
      enddo
      enddo
      do k=1,nz
      do j=ny-1,1,-1
      do i=1,nx
         ty(i,j,k)=(ty(i,j,k)-ffy(j)*ty(i,j+1,k))*fwy(j)
      enddo
      enddo
      enddo
   endif
endif

if (ncly==2) then
   do k=1,nz
   do i=1,nx
      ty(i,1,k)=af1y*uy(i,1,k)+bf1y*uy(i,2,k)+cf1y*uy(i,3,k)
      ty(i,2,k)=af2y*(uy(i,3,k)-uy(i,1,k))
   enddo
   enddo
   do k=1,nz
   do j=3,ny-2
   do i=1,nx
      ty(i,j,k)=afjy*(uy(i,j+1,k)-uy(i,j-1,k))&
           +bfjy*(uy(i,j+2,k)-uy(i,j-2,k))
   enddo
   enddo
   enddo
   do k=1,nz
   do i=1,nx
      ty(i,ny-1,k)=afmy*(uy(i,ny,k)-uy(i,ny-2,k))
      ty(i,ny,k)=-afny*uy(i,ny,k)-bfny*uy(i,ny-1,k)-cfny*uy(i,ny-2,k)
   enddo
   enddo
   do k=1,nz
   do j=2,ny
   do i=1,nx
      ty(i,j,k)=ty(i,j,k)-ty(i,j-1,k)*fsy(j)
   enddo
   enddo
   enddo
   do k=1,nz
   do i=1,nx
      ty(i,ny,k)=ty(i,ny,k)*fwy(ny)
   enddo
   enddo
   do k=1,nz
   do j=ny-1,1,-1
   do i=1,nx
      ty(i,j,k)=(ty(i,j,k)-ffy(j)*ty(i,j+1,k))*fwy(j)
   enddo
   enddo
   enddo
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
end subroutine dery
