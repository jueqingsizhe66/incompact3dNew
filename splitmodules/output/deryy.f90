subroutine deryy(ty,uy,ry,sy,sfy,ssy,swy,nx,ny,nz,npaire)
!
!********************************************************************

USE param
USE derivY

implicit none

integer :: nx,ny,nz,npaire,i,j,k
real(mytype), dimension(nx,ny,nz) :: ty,uy,ry
real(mytype), dimension(nx,nz) :: sy
real(mytype), dimension(ny) :: sfy,ssy,swy

if (ncly==0) then
   do k=1,nz
   do i=1,nx
      ty(i,1,k)=asjy*(uy(i,2,k)-uy(i,1,k)&
           -uy(i,1,k)+uy(i,ny,k))&
           +bsjy*(uy(i,3,k)-uy(i,1,k)&
           -uy(i,1,k)+uy(i,ny-1,k))&
           +csjy*(uy(i,4,k)-uy(i,1,k)&
           -uy(i,1,k)+uy(i,ny-2,k))
      ry(i,1,k)=-1.
      ty(i,2,k)=asjy*(uy(i,3,k)-uy(i,2,k)&
           -uy(i,2,k)+uy(i,1,k))&
           +bsjy*(uy(i,4,k)-uy(i,2,k)&
           -uy(i,2,k)+uy(i,ny,k))&
           +csjy*(uy(i,5,k)-uy(i,2,k)&
           -uy(i,2,k)+uy(i,ny-1,k))
      ry(i,2,k)=0.
      ty(i,3,k)=asjy*(uy(i,4,k)-uy(i,3,k)&
           -uy(i,3,k)+uy(i,2,k))&
           +bsjy*(uy(i,5,k)-uy(i,3,k)&
           -uy(i,3,k)+uy(i,1,k))&
           +csjy*(uy(i,6,k)-uy(i,3,k)&
           -uy(i,3,k)+uy(i,ny,k))
      ry(i,3,k)=0.
   enddo
   enddo
   do k=1,nz
   do j=4,ny-3
   do i=1,nx
      ty(i,j,k)=asjy*(uy(i,j+1,k)-uy(i,j,k)&
           -uy(i,j,k)+uy(i,j-1,k))&
           +bsjy*(uy(i,j+2,k)-uy(i,j,k)&
           -uy(i,j,k)+uy(i,j-2,k))&
           +csjy*(uy(i,j+3,k)-uy(i,j,k)&
           -uy(i,j,k)+uy(i,j-3,k))
      ry(i,j,k)=0.
   enddo
   enddo
   enddo
   do k=1,nz
   do i=1,nx
      ty(i,ny-2,k)=asjy*(uy(i,ny-1,k)-uy(i,ny-2,k)&
           -uy(i,ny-2,k)+uy(i,ny-3,k))&
           +bsjy*(uy(i,ny  ,k)-uy(i,ny-2,k)&
           -uy(i,ny-2,k)+uy(i,ny-4,k))&
           +csjy*(uy(i,1   ,k)-uy(i,ny-2,k)&
           -uy(i,ny-2,k)+uy(i,ny-5,k))
      ry(i,ny-2,k)=0.
      ty(i,ny-1,k)=asjy*(uy(i,ny  ,k)-uy(i,ny-1,k)&
           -uy(i,ny-1,k)+uy(i,ny-2,k))&
           +bsjy*(uy(i,1   ,k)-uy(i,ny-1,k)&
           -uy(i,ny-1,k)+uy(i,ny-3,k))&
           +csjy*(uy(i,2   ,k)-uy(i,ny-1,k)&
           -uy(i,ny-1,k)+uy(i,ny-4,k))
      ry(i,ny-1,k)=0.
      ty(i,ny  ,k)=asjy*(uy(i,1 ,k)-uy(i,ny  ,k)&
           -uy(i,ny,k)+uy(i,ny-1,k))&
           +bsjy*(uy(i,2 ,k)-uy(i,ny  ,k)&
           -uy(i,ny,k)+uy(i,ny-2,k))&
           +csjy*(uy(i,3 ,k)-uy(i,ny  ,k)&
           -uy(i,ny,k)+uy(i,ny-3,k))
      ry(i,ny  ,k)=alsajy
   enddo
   enddo
   do k=1,nz
   do j=2,ny
   do i=1,nx
      ty(i,j,k)=ty(i,j,k)-ty(i,j-1,k)*ssy(j)
      ry(i,j,k)=ry(i,j,k)-ry(i,j-1,k)*ssy(j)
   enddo
   enddo
   enddo
   do k=1,nz
   do i=1,nx
      ty(i,ny,k)=ty(i,ny,k)*swy(ny)
      ry(i,ny,k)=ry(i,ny,k)*swy(ny)
   enddo
   enddo
   do k=1,nz
   do j=ny-1,1,-1
   do i=1,nx
      ty(i,j,k)=(ty(i,j,k)-sfy(j)*ty(i,j+1,k))*swy(j)
      ry(i,j,k)=(ry(i,j,k)-sfy(j)*ry(i,j+1,k))*swy(j)
   enddo
   enddo
   enddo
   do k=1,nz
   do i=1,nx
      sy(i,k)=(   ty(i,1,k)-alsajy*ty(i,ny,k))/&
           (1.+ry(i,1,k)-alsajy*ry(i,ny,k))
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
         ty(i,1,k)=asjy*(uy(i,2,k)-uy(i,1,k)&
              -uy(i,1,k)+uy(i,2,k))&
              +bsjy*(uy(i,3,k)-uy(i,1,k)&
              -uy(i,1,k)+uy(i,3,k))&
              +csjy*(uy(i,4,k)-uy(i,1,k)&
              -uy(i,1,k)+uy(i,4,k))
         ty(i,2,k)=asjy*(uy(i,3,k)-uy(i,2,k)&
              -uy(i,2,k)+uy(i,1,k))&
              +bsjy*(uy(i,4,k)-uy(i,2,k)&
              -uy(i,2,k)+uy(i,2,k))&
              +csjy*(uy(i,5,k)-uy(i,2,k)&
              -uy(i,2,k)+uy(i,3,k))
         ty(i,3,k)=asjy*(uy(i,4,k)-uy(i,3,k)&
              -uy(i,3,k)+uy(i,2,k))&
              +bsjy*(uy(i,5,k)-uy(i,3,k)&
              -uy(i,3,k)+uy(i,1,k))&
              +csjy*(uy(i,6,k)-uy(i,3,k)&
              -uy(i,3,k)+uy(i,2,k))
      enddo
      enddo
      do k=1,nz
      do j=4,ny-3
      do i=1,nx
         ty(i,j,k)=asjy*(uy(i,j+1,k)-uy(i,j  ,k)&
              -uy(i,j  ,k)+uy(i,j-1,k))&
              +bsjy*(uy(i,j+2,k)-uy(i,j  ,k)&
              -uy(i,j  ,k)+uy(i,j-2,k))&
              +csjy*(uy(i,j+3,k)-uy(i,j  ,k)&
              -uy(i,j  ,k)+uy(i,j-3,k))
      enddo
      enddo
      enddo
      do k=1,nz
      do i=1,nx
         ty(i,ny-2,k)=asjy*(uy(i,ny-1,k)-uy(i,ny-2,k)&
              -uy(i,ny-2,k)+uy(i,ny-3,k))&
              +bsjy*(uy(i,ny  ,k)-uy(i,ny-2,k)&
              -uy(i,ny-2,k)+uy(i,ny-4,k))&
              +csjy*(uy(i,ny-1,k)-uy(i,ny-2,k)&
              -uy(i,ny-2,k)+uy(i,ny-5,k))
         ty(i,ny-1,k)=asjy*(uy(i,ny  ,k)-uy(i,ny-1,k)&
              -uy(i,ny-1,k)+uy(i,ny-2,k))&
              +bsjy*(uy(i,ny-1,k)-uy(i,ny-1,k)&
              -uy(i,ny-1,k)+uy(i,ny-3,k))&
              +csjy*(uy(i,ny-2,k)-uy(i,ny-1,k)&
              -uy(i,ny-1,k)+uy(i,ny-4,k))
         ty(i,ny  ,k)=asjy*(uy(i,ny-1,k)-uy(i,ny  ,k)&
              -uy(i,ny  ,k)+uy(i,ny-1,k))&
              +bsjy*(uy(i,ny-2,k)-uy(i,ny  ,k)&
              -uy(i,ny  ,k)+uy(i,ny-2,k))&
              +csjy*(uy(i,ny-3,k)-uy(i,ny  ,k)&
              -uy(i,ny  ,k)+uy(i,ny-3,k))
      enddo
      enddo
      do k=1,nz
      do j=2,ny
      do i=1,nx
         ty(i,j,k)=ty(i,j,k)-ty(i,j-1,k)*ssy(j)
      enddo
      enddo
      enddo
      do k=1,nz
      do i=1,nx
         ty(i,ny,k)=ty(i,ny,k)*swy(ny)
      enddo
      enddo
      do k=1,nz
      do j=ny-1,1,-1
      do i=1,nx
         ty(i,j,k)=(ty(i,j,k)-sfy(j)*ty(i,j+1,k))*swy(j)
      enddo
      enddo
      enddo
   endif
   if (npaire==0) then
      do k=1,nz
      do i=1,nx
         ty(i,1,k)=0.
         ty(i,2,k)=asjy*(uy(i,3,k)-uy(i,2,k)&
              -uy(i,2,k)+uy(i,1,k))&
              +bsjy*(uy(i,4,k)-uy(i,2,k)&
              -uy(i,2,k)-uy(i,2,k))&
              +csjy*(uy(i,5,k)-uy(i,2,k)&
              -uy(i,2,k)-uy(i,3,k))
         ty(i,3,k)=asjy*(uy(i,4,k)-uy(i,3,k)&
              -uy(i,3,k)+uy(i,2,k))&
              +bsjy*(uy(i,5,k)-uy(i,3,k)&
              -uy(i,3,k)+uy(i,1,k))&
              +csjy*(uy(i,6,k)-uy(i,3,k)&
              -uy(i,3,k)-uy(i,2,k))
      enddo
      enddo
      do k=1,nz
      do j=4,ny-3
      do i=1,nx
         ty(i,j,k)=asjy*(uy(i,j+1,k)-uy(i,j  ,k)&
              -uy(i,j  ,k)+uy(i,j-1,k))&
              +bsjy*(uy(i,j+2,k)-uy(i,j  ,k)&
              -uy(i,j  ,k)+uy(i,j-2,k))&
              +csjy*(uy(i,j+3,k)-uy(i,j  ,k)&
              -uy(i,j  ,k)+uy(i,j-3,k))
      enddo
      enddo
      enddo
      do k=1,nz
      do i=1,nx
         ty(i,ny-2,k)=asjy*( uy(i,ny-1,k)-uy(i,ny-2,k)&
              -uy(i,ny-2,k)+uy(i,ny-3,k))&
              +bsjy*( uy(i,ny  ,k)-uy(i,ny-2,k)&
              -uy(i,ny-2,k)+uy(i,ny-4,k))&
              +csjy*(-uy(i,ny-1,k)-uy(i,ny-2,k)&
              -uy(i,ny-2,k)+uy(i,ny-5,k))
         ty(i,ny-1,k)=asjy*( uy(i,ny  ,k)-uy(i,ny-1,k)&
              -uy(i,ny-1,k)+uy(i,ny-2,k))&
              +bsjy*(-uy(i,ny-1,k)-uy(i,ny-1,k)&
              -uy(i,ny-1,k)+uy(i,ny-3,k))&
              +csjy*(-uy(i,ny-2,k)-uy(i,ny-1,k)&
              -uy(i,ny-1,k)+uy(i,ny-4,k))
         ty(i,ny  ,k)=0.
      enddo
      enddo
      do k=1,nz
      do j=2,ny
      do i=1,nx
         ty(i,j,k)=ty(i,j,k)-ty(i,j-1,k)*ssy(j)
      enddo
      enddo
      enddo
      do k=1,nz
      do i=1,nx
         ty(i,ny,k)=ty(i,ny,k)*swy(ny)
      enddo
      enddo
      do k=1,nz
      do j=ny-1,1,-1
      do i=1,nx
         ty(i,j,k)=(ty(i,j,k)-sfy(j)*ty(i,j+1,k))*swy(j)
      enddo
      enddo
      enddo
   endif
endif

if (ncly==2) then
   do k=1,nz
   do i=1,nx
      ty(i,1,k)=as1y*uy(i,1,k)+bs1y*uy(i,2,k)&
           +cs1y*uy(i,3,k)+ds1y*uy(i,4,k)
      ty(i,2,k)=as2y*(uy(i,3,k)-uy(i,2,k)&
           -uy(i,2,k)+uy(i,1,k))
      ty(i,3,k)=as3y*(uy(i,4,k)-uy(i,3,k)&
           -uy(i,3,k)+uy(i,2,k))&
           +bs3y*(uy(i,5,k)-uy(i,3,k)&
           -uy(i,3,k)+uy(i,1,k))
   enddo
   enddo
   do k=1,nz
   do j=4,ny-3
   do i=1,nx
      ty(i,j,k)=asjy*(uy(i,j+1,k)-uy(i,j  ,k)&
           -uy(i,j  ,k)+uy(i,j-1,k))&
           +bsjy*(uy(i,j+2,k)-uy(i,j  ,k)&
           -uy(i,j  ,k)+uy(i,j-2,k))&
           +csjy*(uy(i,j+3,k)-uy(i,j  ,k)&
           -uy(i,j  ,k)+uy(i,j-3,k))
   enddo
   enddo
   enddo
   do k=1,nz
   do i=1,nx
      ty(i,ny-2,k)=asty*(uy(i,ny-1,k)-uy(i,ny-2,k)&
           -uy(i,ny-2,k)+uy(i,ny-3,k))&
           +bsty*(uy(i,ny  ,k)-uy(i,ny-2,k)&
           -uy(i,ny-2,k)+uy(i,ny-4,k))
      ty(i,ny-1,k)=asmy*(uy(i,ny  ,k)-uy(i,ny-1,k)&
           -uy(i,ny-1,k)+uy(i,ny-2,k))
      ty(i,ny  ,k)=asny*uy(i,ny  ,k)+bsny*uy(i,ny-1,k)&
           +csny*uy(i,ny-2,k)+dsny*uy(i,ny-3,k)
   enddo
   enddo
   do k=1,nz
   do j=2,ny
   do i=1,nx
      ty(i,j,k)=ty(i,j,k)-ty(i,j-1,k)*ssy(j)
   enddo
   enddo
   enddo
   do k=1,nz
   do i=1,nx
      ty(i,ny,k)=ty(i,ny,k)*swy(ny)
   enddo
   enddo
   do k=1,nz
   do j=ny-1,1,-1
   do i=1,nx
      ty(i,j,k)=(ty(i,j,k)-sfy(j)*ty(i,j+1,k))*swy(j)
   enddo
   enddo
   enddo
endif

return
end subroutine deryy
