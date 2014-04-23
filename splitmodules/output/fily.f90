subroutine fily(ty,uy,ry,sy,vy,fiffy,fify,ficy,fiby,fibby,filay,fiz1y,&
     fiz2y,nx,ny,nz,npaire)
!
!*********************************************************************

USE param
USE parfiY

implicit none

integer :: nx,ny,nz,i,j,k,npaire
real(mytype), dimension(nx,ny,nz) :: ty,uy,ry
real(mytype), dimension(nx,nz)  :: sy,vy
real(mytype), dimension(ny) :: fiffy,fify,ficy,fiby,fibby,fiz1y,fiz2y
real(mytype), dimension(ny,2) :: filay
real(mytype) :: xcoef

if (ncly==0) then
   do k=1,nz
   do i=1,nx
      ry(i,1,k)=fiaiy*uy(i,1,k)+&
           fibiy*(uy(i,2,k)+uy(i,ny,k))+&
           ficiy*(uy(i,3,k)+uy(i,ny-1,k))+&
           fidiy*(uy(i,4,k)+uy(i,ny-2,k))
      ry(i,2,k)=fiaiy*uy(i,2,k)+&
           fibiy*(uy(i,3,k)+uy(i,1,k))+&
           ficiy*(uy(i,4,k)+uy(i,ny,k))+ &
           fidiy*(uy(i,5,k)+uy(i,ny-1,k))
      ry(i,3,k)=fiaiy*uy(i,3,k)+&
           fibiy*(uy(i,4,k)+uy(i,2,k))+&
           ficiy*(uy(i,5,k)+uy(i,1,k))+&
           fidiy*(uy(i,6,k)+uy(i,ny,k))
   enddo
   enddo
   do j=4,ny-3
   do k=1,nz
   do i=1,nx
      ry(i,j,k)=fiaiy*uy(i,j,k)+&
           fibiy*(uy(i,j+1,k)+uy(i,j-1,k))+&
           ficiy*(uy(i,j+2,k)+uy(i,j-2,k))+&
           fidiy*(uy(i,j+3,k)+uy(i,j-3,k))
   enddo
   enddo
   enddo
   do k=1,nz
   do i=1,nx
      ry(i,ny,k)=fiaiy*uy(i,ny,k)+&
           fibiy*(uy(i,1,k)+uy(i,ny-1,k))+&
           ficiy*(uy(i,2,k)+uy(i,ny-2,k))+&
           fidiy*(uy(i,3,k)+uy(i,ny-3,k))
      ry(i,ny-1,k)=fiaiy*uy(i,ny-1,k)+&
           fibiy*(uy(i,ny,k)+uy(i,ny-2,k))+&
           ficiy*(uy(i,1,k)+uy(i,ny-3,k))+&
           fidiy*(uy(i,2,k)+uy(i,ny-4,k))
      ry(i,ny-2,k)=fiaiy*uy(i,ny-2,k)+&
           fibiy*(uy(i,ny-1,k)+uy(i,ny-3,k))+&
           ficiy*(uy(i,ny,k)+uy(i,ny-4,k))+&
           fidiy*(uy(i,1,k)+uy(i,ny-5,k))
   enddo
   enddo
   do k=1,nz
   do i=1,nx
   do j=1,ny-2
      ry(i,j+1,k)=ry(i,j+1,k)-filay(j,1)*ry(i,j,k)
      ry(i,j+2,k)=ry(i,j+2,k)-filay(j,2)*ry(i,j,k)
   enddo
   enddo
   enddo
   do k=1,nz
   do i=1,nx
      ry(i,ny,k)=ry(i,ny,k)-filay(ny-1,1)*ry(i,ny-1,k)
   enddo
   enddo
   do k=1,nz
   do i=1,nx
      ry(i,ny,k)=ry(i,ny,k)*fiffy(ny)
      ry(i,ny-1,k)=(ry(i,ny-1,k)-fify(ny-1)*ry(i,ny,k))*&
           fiffy(ny-1)
      ry(i,ny-2,k)=(ry(i,ny-2,k)-fify(ny-2)*ry(i,ny-1,k)-&
           ficy(ny-2)*ry(i,ny,k))*fiffy(ny-2)
      ry(i,ny-3,k)=(ry(i,ny-3,k)-fify(ny-3)*ry(i,ny-2,k)-&
           ficy(ny-3)*ry(i,ny-1,k)-&
           fiby(ny-3)*ry(i,ny,k))*fiffy(ny-3)
   enddo
   enddo
   do k=1,nz
   do i=1,nx
   do j=ny-4,1,-1
      ry(i,j,k)=(ry(i,j,k)-fify(j)*ry(i,j+1,k)-&
           ficy(j)*ry(i,j+2,k)-&
           fiby(j)*ry(i,j+3,k)-&
           fibby(j)*ry(i,j+4,k))*fiffy(j)
   enddo
   enddo
   enddo
   xcoef=1./2.
   do k=1,nz
   do i=1,nx
      sy(i,k)=fih1y*(-fibey*ry(i,1,k)+fibey*ry(i,ny-1,k)*xcoef+&
           fialy*ry(i,ny,k)*xcoef)+&
           fih2y*(fialy*ry(i,1,k)*xcoef+fibey*ry(i,2,k)*xcoef-&
           fibey*ry(i,ny,k))
      vy(i,k)=fih3y*(-fibey*ry(i,1,k)+fibey*ry(i,ny-1,k)*xcoef+&
           fialy*ry(i,ny,k)*xcoef)+&
           fih4y*(fialy*ry(i,1,k)*xcoef+fibey*ry(i,2,k)*xcoef-&
           fibey*ry(i,ny,k))
   enddo
   enddo
   do k=1,nz
   do j=1,ny
   do i=1,nx
      ty(i,j,k)=ry(i,j,k)-fiz1y(j)*sy(i,k)-fiz2y(j)*vy(i,k)
   enddo
   enddo
   enddo
endif

if (ncly==1) then
   if (npaire==1) then
      do k=1,nz
      do i=1,nx
         ty(i,1,k)=fiaiy*uy(i,1,k)+&
              fibiy*(uy(i,2,k)+uy(i,2,k))+&
              ficiy*(uy(i,3,k)+uy(i,3,k))+&
              fidiy*(uy(i,4,k)+uy(i,4,k))
         ty(i,2,k)=fiaiy*uy(i,2,k)+&
              fibiy*(uy(i,3,k)+uy(i,1,k))+&
              ficiy*(uy(i,4,k)+uy(i,2,k))+&
              fidiy*(uy(i,5,k)+uy(i,3,k))
         ty(i,3,k)=fiaiy*uy(i,3,k)+&
              fibiy*(uy(i,4,k)+uy(i,2,k))+&
              ficiy*(uy(i,5,k)+uy(i,1,k))+&
              fidiy*(uy(i,6,k)+uy(i,2,k))
      enddo
      enddo
      do j=4,ny-3
      do k=1,nz
      do i=1,nx
         ty(i,j,k)=fiaiy*uy(i,j,k)+&
              fibiy*(uy(i,j+1,k)+uy(i,j-1,k))+&
              ficiy*(uy(i,j+2,k)+uy(i,j-2,k))+&
              fidiy*(uy(i,j+3,k)+uy(i,j-3,k))
      enddo
      enddo
      enddo
      do k=1,nz
      do i=1,nx
         ty(i,ny,k)=fiaiy*uy(i,ny,k)+&
              fibiy*(uy(i,ny-1,k)+uy(i,ny-1,k))+&
              ficiy*(uy(i,ny-2,k)+uy(i,ny-2,k))+&
              fidiy*(uy(i,ny-3,k)+uy(i,ny-3,k))
         ty(i,ny-1,k)=fiaiy*uy(i,ny-1,k)+&
              fibiy*(uy(i,ny,k)+uy(i,ny-2,k))+&
              ficiy*(uy(i,ny-1,k)+uy(i,ny-3,k))+&
              fidiy*(uy(i,ny-2,k)+uy(i,ny-4,k))
         ty(i,ny-2,k)=fiaiy*uy(i,ny-2,k)+&
              fibiy*(uy(i,ny-1,k)+uy(i,ny-3,k))+&
              ficiy*(uy(i,ny,k)+uy(i,ny-4,k))+&
              fidiy*(uy(i,ny-1,k)+uy(i,ny-5,k))
      enddo
      enddo
      do k=1,nz
      do i=1,nx
      do j=1,ny-2
         ty(i,j+1,k)=ty(i,j+1,k)-filay(j,1)*ty(i,j,k)
         ty(i,j+2,k)=ty(i,j+2,k)-filay(j,2)*ty(i,j,k)
      enddo
      ty(i,ny,k)=ty(i,ny,k)-filay(ny-1,1)*ty(i,ny-1,k)
      enddo
      enddo
      do k=1,nz
      do i=1,nx
         ty(i,ny,k)=ty(i,ny,k)*fiffy(ny)
         ty(i,ny-1,k)=(ty(i,ny-1,k)-fify(ny-1)*ty(i,ny,k))*&
              fiffy(ny-1)
         ty(i,ny-2,k)=(ty(i,ny-2,k)-fify(ny-2)*ty(i,ny-1,k)-&
              ficy(ny-2)*ty(i,ny,k))*fiffy(ny-2)
         ty(i,ny-3,k)=(ty(i,ny-3,k)-fify(ny-3)*ty(i,ny-2,k)-&
              ficy(ny-3)*ty(i,ny-1,k)-&
              fiby(ny-3)*ty(i,ny,k))*fiffy(ny-3)
         do j=ny-4,1,-1
            ty(i,j,k)=(ty(i,j,k)-fify(j)*ty(i,j+1,k)-&
                 ficy(j)*ty(i,j+2,k)-&
                 fiby(j)*ty(i,j+3,k)-&
                 fibby(j)*ty(i,j+4,k))*fiffy(j)
         enddo
      enddo
      enddo
   endif
   if (npaire==0) then
      do k=1,nz
      do i=1,nx
         ty(i,1,k)=fiaiy*uy(i,1,k)+&
              fibiy*(uy(i,2,k)-uy(i,2,k))+&
              ficiy*(uy(i,3,k)-uy(i,3,k))+&
              fidiy*(uy(i,4,k)-uy(i,4,k))
         ty(i,2,k)=fiaiy*uy(i,2,k)+&
              fibiy*(uy(i,3,k)+uy(i,1,k))+&
              ficiy*(uy(i,4,k)-uy(i,2,k))+&
              fidiy*(uy(i,5,k)-uy(i,3,k))
         ty(i,3,k)=fiaiy*uy(i,3,k)+&
              fibiy*(uy(i,4,k)+uy(i,2,k))+&
              ficiy*(uy(i,5,k)+uy(i,1,k))+&
              fidiy*(uy(i,6,k)-uy(i,2,k))
      enddo
      enddo
      do j=4,ny-3
      do k=1,nz
      do i=1,nx
         ty(i,j,k)=fiaiy*uy(i,j,k)+&
              fibiy*(uy(i,j+1,k)+uy(i,j-1,k))+&
              ficiy*(uy(i,j+2,k)+uy(i,j-2,k))+&
              fidiy*(uy(i,j+3,k)+uy(i,j-3,k))
      enddo
      enddo
      enddo
      do k=1,nz
      do i=1,nx
         ty(i,ny,k)=fiaiy*uy(i,ny,k)+&
              fibiy*(uy(i,ny-1,k)-uy(i,ny-1,k))+&
              ficiy*(uy(i,ny-2,k)-uy(i,ny-2,k))+&
              fidiy*(uy(i,ny-3,k)-uy(i,ny-3,k))
         ty(i,ny-1,k)=fiaiy*uy(i,ny-1,k)+&
              fibiy*(uy(i,ny,k)+uy(i,ny-2,k))+&
              ficiy*(-uy(i,ny-1,k)+uy(i,ny-3,k))+&
              fidiy*(-uy(i,ny-2,k)+uy(i,ny-4,k))
         ty(i,ny-2,k)=fiaiy*uy(i,ny-2,k)+&
              fibiy*(uy(i,ny-1,k)+uy(i,ny-3,k))+&
              ficiy*(uy(i,ny,k)+uy(i,ny-4,k))+&
              fidiy*(-uy(i,ny-1,k)+uy(i,ny-5,k))
      enddo
      enddo
      do k=1,nz
      do i=1,nx
      do j=1,ny-2
         ty(i,j+1,k)=ty(i,j+1,k)-filay(j,1)*ty(i,j,k)
         ty(i,j+2,k)=ty(i,j+2,k)-filay(j,2)*ty(i,j,k)
      enddo
      ty(i,ny,k)=ty(i,ny,k)-filay(ny-1,1)*ty(i,ny-1,k)
      enddo
      enddo
      do k=1,nz
      do i=1,nx
         ty(i,ny,k)=ty(i,ny,k)*fiffy(ny)
         ty(i,ny-1,k)=(ty(i,ny-1,k)-fify(ny-1)*ty(i,ny,k))*&
              fiffy(ny-1)
         ty(i,ny-2,k)=(ty(i,ny-2,k)-fify(ny-2)*ty(i,ny-1,k)-&
              ficy(ny-2)*ty(i,ny,k))*fiffy(ny-2)
         ty(i,ny-3,k)=(ty(i,ny-3,k)-fify(ny-3)*ty(i,ny-2,k)-&
              ficy(ny-3)*ty(i,ny-1,k)-&
              fiby(ny-3)*ty(i,ny,k))*fiffy(ny-3)
         do j=ny-4,1,-1
            ty(i,j,k)=(ty(i,j,k)-fify(j)*ty(i,j+1,k)-&
                 ficy(j)*ty(i,j+2,k)-&
                 fiby(j)*ty(i,j+3,k)-&
                 fibby(j)*ty(i,j+4,k))*fiffy(j)
         enddo
      enddo
      enddo
   endif
endif

if (ncly.eq.2) then
   do k=1,nz
   do i=1,nx
      ty(i,1,k)=fia1y*uy(i,1,k)+fib1y*uy(i,2,k)+&
           fic1y*uy(i,3,k)+fid1y*uy(i,4,k)+&
           fie1y*uy(i,5,k)
      ty(i,2,k)=fia2y*uy(i,2,k)+fib2y*uy(i,1,k)+&
           fic2y*uy(i,3,k)+fid2y*uy(i,4,k)+&
           fie2y*uy(i,5,k)
      ty(i,3,k)=fia3y*uy(i,3,k)+fib3y*uy(i,1,k)+&
           fic3y*uy(i,2,k)+fid3y*uy(i,4,k)+&
           fie3y*uy(i,5,k)
   enddo
   enddo
   do j=4,ny-3
   do k=1,nz
   do i=1,nx
      ty(i,j,k)=fiaiy*uy(i,j,k)+&
           fibiy*(uy(i,j+1,k)+uy(i,j-1,k))+&
           ficiy*(uy(i,j+2,k)+uy(i,j-2,k))+&
           fidiy*(uy(i,j+3,k)+uy(i,j-3,k))
   enddo
   enddo
   enddo
   do k=1,nz
   do i=1,nx
      ty(i,ny,k)=fiany*uy(i,ny,k)+fibny*uy(i,ny-1,k)+&
           ficny*uy(i,ny-2,k)+fidny*uy(i,ny-3,k)+&
           fieny*uy(i,ny-4,k)
      ty(i,ny-1,k)=fiamy*uy(i,ny-1,k)+fibmy*uy(i,ny,k)+&
           ficmy*uy(i,ny-2,k)+fidmy*uy(i,ny-3,k)+&
           fiemy*uy(i,ny-4,k)
      ty(i,ny-2,k)=fiapy*uy(i,ny-2,k)+fibpy*uy(i,ny,k)+&
           ficpy*uy(i,ny-1,k)+fidpy*uy(i,ny-3,k)+&
           fiepy*uy(i,ny-4,k)
   enddo
   enddo
   do k=1,nz
   do i=1,nx
   do j=1,ny-2
      ty(i,j+1,k)=ty(i,j+1,k)-filay(j,1)*ty(i,j,k)
      ty(i,j+2,k)=ty(i,j+2,k)-filay(j,2)*ty(i,j,k)
   enddo
   ty(i,ny,k)=ty(i,ny,k)-filay(ny-1,1)*ty(i,ny-1,k)
   enddo
   enddo
   do k=1,nz
   do i=1,nx
      ty(i,ny,k)=ty(i,ny,k)*fiffy(ny)
      ty(i,ny-1,k)=(ty(i,ny-1,k)-fify(ny-1)*ty(i,ny,k))*&
           fiffy(ny-1)
      ty(i,ny-2,k)=(ty(i,ny-2,k)-fify(ny-2)*ty(i,ny-1,k)-&
           ficy(ny-2)*ty(i,ny,k))*fiffy(ny-2)
      ty(i,ny-3,k)=(ty(i,ny-3,k)-fify(ny-3)*ty(i,ny-2,k)-&
           ficy(ny-3)*ty(i,ny-1,k)-&
           fiby(ny-3)*ty(i,ny,k))*fiffy(ny-3)
      do j=ny-4,1,-1
         ty(i,j,k)=(ty(i,j,k)-fify(j)*ty(i,j+1,k)-&
              ficy(j)*ty(i,j+2,k)-&
              fiby(j)*ty(i,j+3,k)-&
              fibby(j)*ty(i,j+4,k))*fiffy(j)
      enddo
   enddo
   enddo
endif

return
end subroutine fily
