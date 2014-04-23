subroutine filx(tx,ux,rx,sx,vx,fiffx,fifx,ficx,fibx,fibbx,filax,fiz1x,&
     fiz2x,nx,ny,nz,npaire)
!
!*********************************************************************

USE param
USE parfiX

implicit none

integer :: nx,ny,nz,npaire,i,j,k
real(mytype), dimension(nx,ny,nz) :: tx,ux,rx
real(mytype), dimension(ny,nz) :: sx,vx
real(mytype), dimension(nx) :: fiffx, fifx,ficx,fibx,fibbx,fiz1x,fiz2x
real(mytype), dimension(nx,2) :: filax
real(mytype) :: xcoef

if (nclx==0) then
   do k=1,nz
   do j=1,ny
      rx(1,j,k)=fiaix*ux(1,j,k)+&
           fibix*(ux(2,j,k)+ux(nx,j,k))+&
           ficix*(ux(3,j,k)+ux(nx-1,j,k))+&
           fidix*(ux(4,j,k)+ux(nx-2,j,k))
      rx(2,j,k)=fiaix*ux(2,j,k)+&
           fibix*(ux(3,j,k)+ux(1,j,k))+&
           ficix*(ux(4,j,k)+ux(nx,j,k))+&
           fidix*(ux(5,j,k)+ux(nx-1,j,k))
      rx(3,j,k)=fiaix*ux(3,j,k)+&
           fibix*(ux(4,j,k)+ux(2,j,k))+&
           ficix*(ux(5,j,k)+ux(1,j,k))+&
           fidix*(ux(6,j,k)+ux(nx,j,k))
   enddo
   enddo
   do i=4,nx-3
   do k=1,nz
   do j=1,ny
      rx(i,j,k)=fiaix*ux(i,j,k)+&
           fibix*(ux(i+1,j,k)+ux(i-1,j,k))+&
           ficix*(ux(i+2,j,k)+ux(i-2,j,k))+&
           fidix*(ux(i+3,j,k)+ux(i-3,j,k))
   enddo
   enddo
   enddo
   do k=1,nz
   do j=1,ny
      rx(nx,j,k)=fiaix*ux(nx,j,k)+&
           fibix*(ux(1,j,k)+ux(nx-1,j,k))+&
           ficix*(ux(2,j,k)+ux(nx-2,j,k))+&
           fidix*(ux(3,j,k)+ux(nx-3,j,k))
      rx(nx-1,j,k)=fiaix*ux(nx-1,j,k)+&
           fibix*(ux(nx,j,k)+ux(nx-2,j,k))+&
           ficix*(ux(1,j,k)+ux(nx-3,j,k))+&
           fidix*(ux(2,j,k)+ux(nx-4,j,k))
      rx(nx-2,j,k)=fiaix*ux(nx-2,j,k)+&
           fibix*(ux(nx-1,j,k)+ux(nx-3,j,k))+&
           ficix*(ux(nx,j,k)+ux(nx-4,j,k))+&
           fidix*(ux(1,j,k)+ux(nx-5,j,k))
   enddo
   enddo
   do k=1,nz
   do j=1,ny
   do i=1,nx-2
      rx(i+1,j,k)=rx(i+1,j,k)-filax(i,1)*rx(i,j,k)
      rx(i+2,j,k)=rx(i+2,j,k)-filax(i,2)*rx(i,j,k)
   enddo
   enddo
   enddo
   do k=1,nz
   do j=1,ny
      rx(nx,j,k)=rx(nx,j,k)-filax(nx-1,1)*rx(nx-1,j,k)
   enddo
   enddo
   do k=1,nz
   do j=1,ny
      rx(nx,j,k)=rx(nx,j,k)*fiffx(nx)
      rx(nx-1,j,k)=(rx(nx-1,j,k)-fifx(nx-1)*rx(nx,j,k))*&
           fiffx(nx-1)
      rx(nx-2,j,k)=(rx(nx-2,j,k)-fifx(nx-2)*rx(nx-1,j,k)-&
           ficx(nx-2)*rx(nx,j,k))*fiffx(nx-2)
      rx(nx-3,j,k)=(rx(nx-3,j,k)-fifx(nx-3)*rx(nx-2,j,k)-&
           ficx(nx-3)*rx(nx-1,j,k)-&
           fibx(nx-3)*rx(nx,j,k))*fiffx(nx-3)
   enddo
   enddo
   do k=1,nz
   do j=1,ny
   do i=nx-4,1,-1
      rx(i,j,k)=(rx(i,j,k)-fifx(i)*rx(i+1,j,k)-&
           ficx(i)*rx(i+2,j,k)-&
           fibx(i)*rx(i+3,j,k)-&
           fibbx(i)*rx(i+4,j,k))*fiffx(i)
   enddo
   enddo
   enddo
   xcoef=1./2.
   do k=1,nz
   do j=1,ny
      sx(j,k)=fih1x*(-fibex*rx(1,j,k)+fibex*rx(nx-1,j,k)*xcoef+&
           fialx*rx(nx,j,k)*xcoef)+&
           fih2x*(fialx*rx(1,j,k)*xcoef+fibex*rx(2,j,k)*xcoef-&
           fibex*rx(nx,j,k))
      vx(j,k)=fih3x*(-fibex*rx(1,j,k)+fibex*rx(nx-1,j,k)*xcoef+&
           fialx*rx(nx,j,k)*xcoef)+&
           fih4x*(fialx*rx(1,j,k)*xcoef+fibex*rx(2,j,k)*xcoef-&
                    fibex*rx(nx,j,k))
   enddo
   enddo
   do k=1,nz
   do j=1,ny
   do i=1,nx
      tx(i,j,k)=rx(i,j,k)-fiz1x(i)*sx(j,k)-fiz2x(i)*vx(j,k)
   enddo
   enddo
   enddo
endif

if (nclx==1) then
   if (npaire==1) then
      do k=1,nz
      do j=1,ny
         tx(1,j,k)=fiaix*ux(1,j,k)+&
              fibix*(ux(2,j,k)+ux(2,j,k))+&
              ficix*(ux(3,j,k)+ux(3,j,k))+&
              fidix*(ux(4,j,k)+ux(4,j,k))
         tx(2,j,k)=fiaix*ux(2,j,k)+&
              fibix*(ux(3,j,k)+ux(1,j,k))+&
              ficix*(ux(4,j,k)+ux(2,j,k))+&
              fidix*(ux(5,j,k)+ux(3,j,k))
         tx(3,j,k)=fiaix*ux(3,j,k)+&
              fibix*(ux(4,j,k)+ux(2,j,k))+&
              ficix*(ux(5,j,k)+ux(1,j,k))+&
              fidix*(ux(6,j,k)+ux(2,j,k))
      enddo
      enddo
      do i=4,nx-3
      do k=1,nz
      do j=1,ny
         tx(i,j,k)=fiaix*ux(i,j,k)+&
              fibix*(ux(i+1,j,k)+ux(i-1,j,k))+&
              ficix*(ux(i+2,j,k)+ux(i-2,j,k))+&
              fidix*(ux(i+3,j,k)+ux(i-3,j,k))
      enddo
      enddo
      enddo
      do k=1,nz
      do j=1,ny
         tx(nx,j,k)=fiaix*ux(nx,j,k)+&
              fibix*(ux(nx-1,j,k)+ux(nx-1,j,k))+&
              ficix*(ux(nx-2,j,k)+ux(nx-2,j,k))+&
              fidix*(ux(nx-3,j,k)+ux(nx-3,j,k))
         tx(nx-1,j,k)=fiaix*ux(nx-1,j,k)+&
              fibix*(ux(nx,j,k)+ux(nx-2,j,k))+&
              ficix*(ux(nx-1,j,k)+ux(nx-3,j,k))+&
              fidix*(ux(nx-2,j,k)+ux(nx-4,j,k))
         tx(nx-2,j,k)=fiaix*ux(nx-2,j,k)+&
              fibix*(ux(nx-1,j,k)+ux(nx-3,j,k))+&
              ficix*(ux(nx,j,k)+ux(nx-4,j,k))+&
              fidix*(ux(nx-1,j,k)+ux(nx-5,j,k))
      enddo
      enddo
      do k=1,nz
      do j=1,ny
      do i=1,nx-2
         tx(i+1,j,k)=tx(i+1,j,k)-filax(i,1)*tx(i,j,k)
         tx(i+2,j,k)=tx(i+2,j,k)-filax(i,2)*tx(i,j,k)
      enddo
      tx(nx,j,k)=tx(nx,j,k)-filax(nx-1,1)*tx(nx-1,j,k)
      enddo
      enddo
      do k=1,nz
      do j=1,ny
         tx(nx,j,k)=tx(nx,j,k)*fiffx(nx)
         tx(nx-1,j,k)=(tx(nx-1,j,k)-fifx(nx-1)*tx(nx,j,k))*&
              fiffx(nx-1)
         tx(nx-2,j,k)=(tx(nx-2,j,k)-fifx(nx-2)*tx(nx-1,j,k)-&
              ficx(nx-2)*tx(nx,j,k))*fiffx(nx-2)
         tx(nx-3,j,k)=(tx(nx-3,j,k)-fifx(nx-3)*tx(nx-2,j,k)-&
              ficx(nx-3)*tx(nx-1,j,k)-&
              fibx(nx-3)*tx(nx,j,k))*fiffx(nx-3)
         do i=nx-4,1,-1
            tx(i,j,k)=(tx(i,j,k)-fifx(i)*tx(i+1,j,k)-&
                 ficx(i)*tx(i+2,j,k)-&
                 fibx(i)*tx(i+3,j,k)-&
                 fibbx(i)*tx(i+4,j,k))*fiffx(i)
         enddo
      enddo
      enddo
   endif
   if (npaire==0) then
      do k=1,nz
      do j=1,ny
         tx(1,j,k)=fiaix*ux(1,j,k)+&
              fibix*(ux(2,j,k)-ux(2,j,k))+&
              ficix*(ux(3,j,k)-ux(3,j,k))+&
              fidix*(ux(4,j,k)-ux(4,j,k))
         tx(2,j,k)=fiaix*ux(2,j,k)+&
              fibix*(ux(3,j,k)+ux(1,j,k))+&
              ficix*(ux(4,j,k)-ux(2,j,k))+&
              fidix*(ux(5,j,k)-ux(3,j,k))
         tx(3,j,k)=fiaix*ux(3,j,k)+&
              fibix*(ux(4,j,k)+ux(2,j,k))+&
              ficix*(ux(5,j,k)+ux(1,j,k))+&
              fidix*(ux(6,j,k)-ux(2,j,k))
      enddo
      enddo
      do i=4,nx-3
      do k=1,nz
      do j=1,ny
         tx(i,j,k)=fiaix*ux(i,j,k)+&
              fibix*(ux(i+1,j,k)+ux(i-1,j,k))+&
              ficix*(ux(i+2,j,k)+ux(i-2,j,k))+&
              fidix*(ux(i+3,j,k)+ux(i-3,j,k))
      enddo
      enddo
      enddo
      do k=1,nz
      do j=1,ny
         tx(nx,j,k)=fiaix*ux(nx,j,k)+&
              fibix*(ux(nx-1,j,k)-ux(nx-1,j,k))+&
              ficix*(ux(nx-2,j,k)-ux(nx-2,j,k))+&
              fidix*(ux(nx-3,j,k)-ux(nx-3,j,k))
         tx(nx-1,j,k)=fiaix*ux(nx-1,j,k)+&
              fibix*(ux(nx,j,k)+ux(nx-2,j,k))+&
              ficix*(-ux(nx-1,j,k)+ux(nx-3,j,k))+&
              fidix*(-ux(nx-2,j,k)+ux(nx-4,j,k))
         tx(nx-2,j,k)=fiaix*ux(nx-2,j,k)+&
              fibix*(ux(nx-1,j,k)+ux(nx-3,j,k))+&
              ficix*(ux(nx,j,k)+ux(nx-4,j,k))+&
              fidix*(-ux(nx-1,j,k)+ux(nx-5,j,k))
      enddo
      enddo
      do k=1,nz
      do j=1,ny
      do i=1,nx-2
         tx(i+1,j,k)=tx(i+1,j,k)-filax(i,1)*tx(i,j,k)
         tx(i+2,j,k)=tx(i+2,j,k)-filax(i,2)*tx(i,j,k)
      enddo
      tx(nx,j,k)=tx(nx,j,k)-filax(nx-1,1)*tx(nx-1,j,k)
      enddo
      enddo
      do k=1,nz
      do j=1,ny
         tx(nx,j,k)=tx(nx,j,k)*fiffx(nx)
         tx(nx-1,j,k)=(tx(nx-1,j,k)-fifx(nx-1)*tx(nx,j,k))*&
              fiffx(nx-1)
         tx(nx-2,j,k)=(tx(nx-2,j,k)-fifx(nx-2)*tx(nx-1,j,k)-&
              ficx(nx-2)*tx(nx,j,k))*fiffx(nx-2)
         tx(nx-3,j,k)=(tx(nx-3,j,k)-fifx(nx-3)*tx(nx-2,j,k)-&
              ficx(nx-3)*tx(nx-1,j,k)-&
              fibx(nx-3)*tx(nx,j,k))*fiffx(nx-3)
         do i=nx-4,1,-1
            tx(i,j,k)=(tx(i,j,k)-fifx(i)*tx(i+1,j,k)-&
                 ficx(i)*tx(i+2,j,k)-&
                 fibx(i)*tx(i+3,j,k)-&
                 fibbx(i)*tx(i+4,j,k))*fiffx(i)
         enddo
      enddo
      enddo
   endif
endif

if (nclx==2) then
   do k=1,nz
   do j=1,ny
      tx(1,j,k)=fia1x*ux(1,j,k)+fib1x*ux(2,j,k)+&
           fic1x*ux(3,j,k)+fid1x*ux(4,j,k)+&
           fie1x*ux(5,j,k)
      tx(2,j,k)=fia2x*ux(2,j,k)+fib2x*ux(1,j,k)+&
           fic2x*ux(3,j,k)+fid2x*ux(4,j,k)+&
           fie2x*ux(5,j,k)
      tx(3,j,k)=fia3x*ux(3,j,k)+fib3x*ux(1,j,k)+&
           fic3x*ux(2,j,k)+fid3x*ux(4,j,k)+&
           fie3x*ux(5,j,k)
   enddo
   enddo
   do i=4,nx-3
   do k=1,nz
   do j=1,ny
      tx(i,j,k)=fiaix*ux(i,j,k)+&
           fibix*(ux(i+1,j,k)+ux(i-1,j,k))+&
           ficix*(ux(i+2,j,k)+ux(i-2,j,k))+&
           fidix*(ux(i+3,j,k)+ux(i-3,j,k))
   enddo
   enddo
   enddo
   do k=1,nz
   do j=1,ny
      tx(nx,j,k)=fianx*ux(nx,j,k)+fibnx*ux(nx-1,j,k)+&
           ficnx*ux(nx-2,j,k)+fidnx*ux(nx-3,j,k)+&
           fienx*ux(nx-4,j,k)
      tx(nx-1,j,k)=fiamx*ux(nx-1,j,k)+fibmx*ux(nx,j,k)+&
           ficmx*ux(nx-2,j,k)+fidmx*ux(nx-3,j,k)+&
           fiemx*ux(nx-4,j,k)
      tx(nx-2,j,k)=fiapx*ux(nx-2,j,k)+fibpx*ux(nx,j,k)+&
           ficpx*ux(nx-1,j,k)+fidpx*ux(nx-3,j,k)+&
           fiepx*ux(nx-4,j,k)
   enddo
   enddo
   do i=1,nx-2
   do k=1,nz
   do j=1,ny
      tx(i+1,j,k)=tx(i+1,j,k)-filax(i,1)*tx(i,j,k)
      tx(i+2,j,k)=tx(i+2,j,k)-filax(i,2)*tx(i,j,k)
   enddo
   enddo
   enddo
   do k=1,nz
   do j=1,ny
      tx(nx,j,k)=tx(nx,j,k)-filax(nx-1,1)*tx(nx-1,j,k)
   enddo
   enddo
   do k=1,nz
   do j=1,ny
      tx(nx,j,k)=tx(nx,j,k)*fiffx(nx)
      tx(nx-1,j,k)=(tx(nx-1,j,k)-fifx(nx-1)*tx(nx,j,k))*&
           fiffx(nx-1)
      tx(nx-2,j,k)=(tx(nx-2,j,k)-fifx(nx-2)*tx(nx-1,j,k)-&
           ficx(nx-2)*tx(nx,j,k))*fiffx(nx-2)
      tx(nx-3,j,k)=(tx(nx-3,j,k)-fifx(nx-3)*tx(nx-2,j,k)-&
           ficx(nx-3)*tx(nx-1,j,k)-&
           fibx(nx-3)*tx(nx,j,k))*fiffx(nx-3)
   enddo
   enddo
   do i=nx-4,1,-1
   do k=1,nz
   do j=1,ny
      tx(i,j,k)=(tx(i,j,k)-fifx(i)*tx(i+1,j,k)-&
           ficx(i)*tx(i+2,j,k)-&
           fibx(i)*tx(i+3,j,k)-&
           fibbx(i)*tx(i+4,j,k))*fiffx(i)
   enddo
   enddo
   enddo
endif

return
end subroutine filx
