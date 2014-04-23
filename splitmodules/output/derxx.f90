subroutine derxx(tx,ux,rx,sx,sfx,ssx,swx,nx,ny,nz,npaire)
!
!********************************************************************

USE param
USE derivX

implicit none

integer :: nx,ny,nz,npaire,i,j,k
real(mytype), dimension(nx,ny,nz) :: tx,ux,rx
real(mytype), dimension(ny,nz) :: sx
real(mytype),  dimension(nx):: sfx,ssx,swx

if (nclx==0) then
   do k=1,nz
   do j=1,ny
      tx(1,j,k)=asix*(ux(2,j,k)-ux(1   ,j,k)&
           -ux(1,j,k)+ux(nx  ,j,k))&
           +bsix*(ux(3,j,k)-ux(1   ,j,k)&
           -ux(1,j,k)+ux(nx-1,j,k))&
           +csix*(ux(4,j,k)-ux(1   ,j,k)&
           -ux(1,j,k)+ux(nx-2,j,k))
      rx(1,j,k)=-1.
      tx(2,j,k)=asix*(ux(3,j,k)-ux(2   ,j,k)&
           -ux(2,j,k)+ux(1   ,j,k))&
           +bsix*(ux(4,j,k)-ux(2   ,j,k)&
           -ux(2,j,k)+ux(nx  ,j,k))&
           +csix*(ux(5,j,k)-ux(2   ,j,k)&
           -ux(2,j,k)+ux(nx-1,j,k))
      rx(2,j,k)=0.
      tx(3,j,k)=asix*(ux(4,j,k)-ux(3 ,j,k)&
           -ux(3,j,k)+ux(2 ,j,k))&
           +bsix*(ux(5,j,k)-ux(3 ,j,k)&
           -ux(3,j,k)+ux(1 ,j,k))&
           +csix*(ux(6,j,k)-ux(3 ,j,k)&
           -ux(3,j,k)+ux(nx,j,k))
      rx(3,j,k)=0.
      do i=4,nx-3
         tx(i,j,k)=asix*(ux(i+1,j,k)-ux(i  ,j,k)&
              -ux(i  ,j,k)+ux(i-1,j,k))&
              +bsix*(ux(i+2,j,k)-ux(i  ,j,k)&
              -ux(i  ,j,k)+ux(i-2,j,k))&
              +csix*(ux(i+3,j,k)-ux(i  ,j,k)&
              -ux(i  ,j,k)+ux(i-3,j,k))
         rx(i,j,k)=0.
      enddo
      tx(nx-2,j,k)=asix*(ux(nx-1,j,k)-ux(nx-2,j,k)&
           -ux(nx-2,j,k)+ux(nx-3,j,k))&
           +bsix*(ux(nx  ,j,k)-ux(nx-2,j,k)&
           -ux(nx-2,j,k)+ux(nx-4,j,k))&
           +csix*(ux(1   ,j,k)-ux(nx-2,j,k)&
           -ux(nx-2,j,k)+ux(nx-5,j,k))
      rx(nx-2,j,k)=0.
      tx(nx-1,j,k)=asix*(ux(nx  ,j,k)-ux(nx-1,j,k)&
           -ux(nx-1,j,k)+ux(nx-2,j,k))&
           +bsix*(ux(1   ,j,k)-ux(nx-1,j,k)&
           -ux(nx-1,j,k)+ux(nx-3,j,k))&
           +csix*(ux(2   ,j,k)-ux(nx-1,j,k)&
           -ux(nx-1,j,k)+ux(nx-4,j,k))
      rx(nx-1,j,k)=0.
      tx(nx  ,j,k)=asix*(ux(1 ,j,k)-ux(nx  ,j,k)&
           -ux(nx,j,k)+ux(nx-1,j,k))&
           +bsix*(ux(2 ,j,k)-ux(nx  ,j,k)&
           -ux(nx,j,k)+ux(nx-2,j,k))&
           +csix*(ux(3 ,j,k)-ux(nx  ,j,k)&
           -ux(nx,j,k)+ux(nx-3,j,k))
      rx(nx  ,j,k)=alsaix
      do i=2,nx
         tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*ssx(i)
         rx(i,j,k)=rx(i,j,k)-rx(i-1,j,k)*ssx(i)
      enddo
         tx(nx,j,k)=tx(nx,j,k)*swx(nx)
         rx(nx,j,k)=rx(nx,j,k)*swx(nx)
      do i=nx-1,1,-1
         tx(i,j,k)=(tx(i,j,k)-sfx(i)*tx(i+1,j,k))*swx(i)
         rx(i,j,k)=(rx(i,j,k)-sfx(i)*rx(i+1,j,k))*swx(i)
      enddo
      sx(j,k)=(   tx(1,j,k)-alsaix*tx(nx,j,k))/&
           (1.+rx(1,j,k)-alsaix*rx(nx,j,k))
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
         tx(1,j,k)=asix*(ux(2,j,k)-ux(1,j,k)&
              -ux(1,j,k)+ux(2,j,k))&
              +bsix*(ux(3,j,k)-ux(1,j,k)&
              -ux(1,j,k)+ux(3,j,k))&
              +csix*(ux(4,j,k)-ux(1,j,k)&
              -ux(1,j,k)+ux(4,j,k))
         tx(2,j,k)=asix*(ux(3,j,k)-ux(2,j,k)&
              -ux(2,j,k)+ux(1,j,k))&
              +bsix*(ux(4,j,k)-ux(2,j,k)&
              -ux(2,j,k)+ux(2,j,k))&
              +csix*(ux(5,j,k)-ux(2,j,k)&
              -ux(2,j,k)+ux(3,j,k))
         tx(3,j,k)=asix*(ux(4,j,k)-ux(3,j,k)&
              -ux(3,j,k)+ux(2,j,k))&
              +bsix*(ux(5,j,k)-ux(3,j,k)&
              -ux(3,j,k)+ux(1,j,k))&
              +csix*(ux(6,j,k)-ux(3,j,k)&
              -ux(3,j,k)+ux(2,j,k))
         do i=4,nx-3
            tx(i,j,k)=asix*(ux(i+1,j,k)-ux(i  ,j,k)&
                 -ux(i  ,j,k)+ux(i-1,j,k))&
                 +bsix*(ux(i+2,j,k)-ux(i  ,j,k)&
                 -ux(i  ,j,k)+ux(i-2,j,k))&
                 +csix*(ux(i+3,j,k)-ux(i  ,j,k)&
                 -ux(i  ,j,k)+ux(i-3,j,k))
         enddo
         tx(nx-2,j,k)=asix*(ux(nx-1,j,k)-ux(nx-2,j,k)&
              -ux(nx-2,j,k)+ux(nx-3,j,k))&
              +bsix*(ux(nx  ,j,k)-ux(nx-2,j,k)&
              -ux(nx-2,j,k)+ux(nx-4,j,k))&
              +csix*(ux(nx-1,j,k)-ux(nx-2,j,k)&
              -ux(nx-2,j,k)+ux(nx-5,j,k))
         tx(nx-1,j,k)=asix*(ux(nx  ,j,k)-ux(nx-1,j,k)&
              -ux(nx-1,j,k)+ux(nx-2,j,k))&
              +bsix*(ux(nx-1,j,k)-ux(nx-1,j,k)&
              -ux(nx-1,j,k)+ux(nx-3,j,k))&
              +csix*(ux(nx-2,j,k)-ux(nx-1,j,k)&
              -ux(nx-1,j,k)+ux(nx-4,j,k))
         tx(nx  ,j,k)=asix*(ux(nx-1,j,k)-ux(nx  ,j,k)&
              -ux(nx  ,j,k)+ux(nx-1,j,k))&
              +bsix*(ux(nx-2,j,k)-ux(nx  ,j,k)&
              -ux(nx  ,j,k)+ux(nx-2,j,k))&
              +csix*(ux(nx-3,j,k)-ux(nx  ,j,k)&
              -ux(nx  ,j,k)+ux(nx-3,j,k))
         do i=2,nx
            tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*ssx(i)
         enddo
         tx(nx,j,k)=tx(nx,j,k)*swx(nx)
         do i=nx-1,1,-1
            tx(i,j,k)=(tx(i,j,k)-sfx(i)*tx(i+1,j,k))*swx(i)
         enddo
      enddo
      enddo
   endif
   if (npaire==0) then
      do k=1,nz
      do j=1,ny
         tx(1,j,k)=0.
         tx(2,j,k)=asix*(ux(3,j,k)-ux(2,j,k)&
              -ux(2,j,k)+ux(1,j,k))&
              +bsix*(ux(4,j,k)-ux(2,j,k)&
              -ux(2,j,k)-ux(2,j,k))&
              +csix*(ux(5,j,k)-ux(2,j,k)&
              -ux(2,j,k)-ux(3,j,k))
         tx(3,j,k)=asix*(ux(4,j,k)-ux(3,j,k)&
              -ux(3,j,k)+ux(2,j,k))&
              +bsix*(ux(5,j,k)-ux(3,j,k)&
              -ux(3,j,k)+ux(1,j,k))&
              +csix*(ux(6,j,k)-ux(3,j,k)&
              -ux(3,j,k)-ux(2,j,k))
         do i=4,nx-3
            tx(i,j,k)=asix*(ux(i+1,j,k)-ux(i  ,j,k)&
                 -ux(i  ,j,k)+ux(i-1,j,k))&
                 +bsix*(ux(i+2,j,k)-ux(i  ,j,k)&
                 -ux(i  ,j,k)+ux(i-2,j,k))&
                 +csix*(ux(i+3,j,k)-ux(i  ,j,k)&
                 -ux(i  ,j,k)+ux(i-3,j,k))
         enddo
         tx(nx-2,j,k)=asix*( ux(nx-1,j,k)-ux(nx-2,j,k)&
              -ux(nx-2,j,k)+ux(nx-3,j,k))&
              +bsix*( ux(nx  ,j,k)-ux(nx-2,j,k)&
              -ux(nx-2,j,k)+ux(nx-4,j,k))&
              +csix*(-ux(nx-1,j,k)-ux(nx-2,j,k)&
              -ux(nx-2,j,k)+ux(nx-5,j,k))
         tx(nx-1,j,k)=asix*( ux(nx  ,j,k)-ux(nx-1,j,k)&
              -ux(nx-1,j,k)+ux(nx-2,j,k))&
              +bsix*(-ux(nx-1,j,k)-ux(nx-1,j,k)&
              -ux(nx-1,j,k)+ux(nx-3,j,k))&
              +csix*(-ux(nx-2,j,k)-ux(nx-1,j,k)&
              -ux(nx-1,j,k)+ux(nx-4,j,k))
         tx(nx  ,j,k)=0.
         do i=2,nx
            tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*ssx(i)
         enddo
         tx(nx,j,k)=tx(nx,j,k)*swx(nx)
         do i=nx-1,1,-1
            tx(i,j,k)=(tx(i,j,k)-sfx(i)*tx(i+1,j,k))*swx(i)
         enddo
      enddo
      enddo
   endif
endif

if (nclx==2) then
   do k=1,nz
   do j=1,ny
      tx(1,j,k)=as1x*ux(1,j,k)+bs1x*ux(2,j,k)&
           +cs1x*ux(3,j,k)+ds1x*ux(4,j,k)
      tx(2,j,k)=as2x*(ux(3,j,k)-ux(2,j,k)&
           -ux(2,j,k)+ux(1,j,k))
      tx(3,j,k)=as3x*(ux(4,j,k)-ux(3,j,k)&
           -ux(3,j,k)+ux(2,j,k))&
           +bs3x*(ux(5,j,k)-ux(3,j,k)&
           -ux(3,j,k)+ux(1,j,k))
      do i=4,nx-3
         tx(i,j,k)=asix*(ux(i+1,j,k)-ux(i  ,j,k)&
              -ux(i  ,j,k)+ux(i-1,j,k))&
              +bsix*(ux(i+2,j,k)-ux(i  ,j,k)&
              -ux(i  ,j,k)+ux(i-2,j,k))&
              +csix*(ux(i+3,j,k)-ux(i  ,j,k)&
              -ux(i  ,j,k)+ux(i-3,j,k))
      enddo
      tx(nx-2,j,k)=astx*(ux(nx-1,j,k)-ux(nx-2,j,k)&
           -ux(nx-2,j,k)+ux(nx-3,j,k))&
           +bstx*(ux(nx  ,j,k)-ux(nx-2,j,k)&
           -ux(nx-2,j,k)+ux(nx-4,j,k))
      tx(nx-1,j,k)=asmx*(ux(nx  ,j,k)-ux(nx-1,j,k)&
           -ux(nx-1,j,k)+ux(nx-2,j,k))
      tx(nx  ,j,k)=asnx*ux(nx  ,j,k)+bsnx*ux(nx-1,j,k)&
           +csnx*ux(nx-2,j,k)+dsnx*ux(nx-3,j,k)
      do i=2,nx
         tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*ssx(i)
      enddo
      tx(nx,j,k)=tx(nx,j,k)*swx(nx)
      do i=nx-1,1,-1
         tx(i,j,k)=(tx(i,j,k)-sfx(i)*tx(i+1,j,k))*swx(i)
      enddo
   enddo
   enddo
endif

return
end subroutine derxx
