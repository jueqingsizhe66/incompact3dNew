subroutine derx(tx,ux,rx,sx,ffx,fsx,fwx,nx,ny,nz,npaire)
!
!********************************************************************

USE param
USE derivX

implicit none

integer :: nx,ny,nz,npaire,i,j,k
real(mytype), dimension(nx,ny,nz) :: tx,ux,rx
real(mytype), dimension(ny,nz):: sx
real(mytype), dimension(nx):: ffx,fsx,fwx

if (nclx==0) then
   do k=1,nz
   do j=1,ny
      tx(1,j,k)=afix*(ux(2,j,k)-ux(nx,j,k))&
           +bfix*(ux(3,j,k)-ux(nx-1,j,k))
      rx(1,j,k)=-1.
      tx(2,j,k)=afix*(ux(3,j,k)-ux(1,j,k))&
           +bfix*(ux(4,j,k)-ux(nx,j,k))
      rx(2,j,k)=0.
      do i=3,nx-2
         tx(i,j,k)=afix*(ux(i+1,j,k)-ux(i-1,j,k))&
              +bfix*(ux(i+2,j,k)-ux(i-2,j,k))
         rx(i,j,k)=0.
      enddo
      tx(nx-1,j,k)=afix*(ux(nx,j,k)-ux(nx-2,j,k))&
           +bfix*(ux(1,j,k)-ux(nx-3,j,k))
      rx(nx-1,j,k)=0.
      tx(nx,j,k)=afix*(ux(1,j,k)-ux(nx-1,j,k))&
           +bfix*(ux(2,j,k)-ux(nx-2,j,k))
      rx(nx,j,k)=alfaix
      do i=2, nx
         tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*fsx(i)
         rx(i,j,k)=rx(i,j,k)-rx(i-1,j,k)*fsx(i)
      enddo
      tx(nx,j,k)=tx(nx,j,k)*fwx(nx)
      rx(nx,j,k)=rx(nx,j,k)*fwx(nx)
      do i=nx-1,1,-1
         tx(i,j,k)=(tx(i,j,k)-ffx(i)*tx(i+1,j,k))*fwx(i)
         rx(i,j,k)=(rx(i,j,k)-ffx(i)*rx(i+1,j,k))*fwx(i)
      enddo
      sx(j,k)=(tx(1,j,k)-alfaix*tx(nx,j,k))&
           /(1.+rx(1,j,k)-alfaix*rx(nx,j,k))
      do i=1,nx
         tx(i,j,k)=tx(i,j,k)-sx(j,k)*rx(i,j,k)
      enddo
   enddo
   enddo
endif

if (nclx==1) then
   if (npaire==1) then
      do k=1,nz
      do j=1,ny
         tx(1,j,k)=0.
         tx(2,j,k)=afix*(ux(3,j,k)-ux(1,j,k))&
              +bfix*(ux(4,j,k)-ux(2,j,k))
         do i=3,nx-2
            tx(i,j,k)=afix*(ux(i+1,j,k)-ux(i-1,j,k))&
                 +bfix*(ux(i+2,j,k)-ux(i-2,j,k))
         enddo
         tx(nx-1,j,k)=afix*(ux(nx,j,k)-ux(nx-2,j,k))&
              +bfix*(ux(nx-1,j,k)-ux(nx-3,j,k))
         tx(nx,j,k)=0.
         do i=2,nx
            tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*fsx(i)
         enddo
         tx(nx,j,k)=tx(nx,j,k)*fwx(nx)
         do i=nx-1,1,-1
            tx(i,j,k)=(tx(i,j,k)-ffx(i)*tx(i+1,j,k))*fwx(i)
         enddo
      enddo
      enddo
   endif
   if (npaire==0) then
      do k=1,nz
      do j=1,ny
         tx(1,j,k)=afix*(ux(2,j,k)+ux(2,j,k))&
              +bfix*(ux(3,j,k)+ux(3,j,k))
         tx(2,j,k)=afix*(ux(3,j,k)-ux(1,j,k))&
              +bfix*(ux(4,j,k)+ux(2,j,k))
         do i=3,nx-2
            tx(i,j,k)=afix*(ux(i+1,j,k)-ux(i-1,j,k))&
                 +bfix*(ux(i+2,j,k)-ux(i-2,j,k))
         enddo
         tx(nx-1,j,k)=afix*(ux(nx,j,k)-ux(nx-2,j,k))&
              +bfix*((-ux(nx-1,j,k))-ux(nx-3,j,k))
         tx(nx,j,k)=afix*((-ux(nx-1,j,k))-ux(nx-1,j,k))&
              +bfix*((-ux(nx-2,j,k))-ux(nx-2,j,k))
         do i=2,nx
            tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*fsx(i)
         enddo
         tx(nx,j,k)=tx(nx,j,k)*fwx(nx)
         do i=nx-1,1,-1
            tx(i,j,k)=(tx(i,j,k)-ffx(i)*tx(i+1,j,k))*fwx(i)
         enddo
      enddo
      enddo
   endif
endif

if (nclx==2) then
   do k=1,nz
   do j=1,ny
      tx(1,j,k)=af1x*ux(1,j,k)+bf1x*ux(2,j,k)+cf1x*ux(3,j,k)
      tx(2,j,k)=af2x*(ux(3,j,k)-ux(1,j,k))
      do i=3,nx-2
         tx(i,j,k)=afix*(ux(i+1,j,k)-ux(i-1,j,k))&
              +bfix*(ux(i+2,j,k)-ux(i-2,j,k))
      enddo
      tx(nx-1,j,k)=afmx*(ux(nx,j,k)-ux(nx-2,j,k))
      tx(nx,j,k)=(-afnx*ux(nx,j,k))-bfnx*ux(nx-1,j,k)-cfnx*ux(nx-2,j,k)
      do i=2,nx
         tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*fsx(i)
      enddo
      tx(nx,j,k)=tx(nx,j,k)*fwx(nx)
      do i=nx-1,1,-1
         tx(i,j,k)=(tx(i,j,k)-ffx(i)*tx(i+1,j,k))*fwx(i)
      enddo
   enddo
   enddo
endif

return
end subroutine derx
