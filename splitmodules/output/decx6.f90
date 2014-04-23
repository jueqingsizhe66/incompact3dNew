subroutine decx6(tx,ux,rx,sx,cfx6,csx6,cwx6,nx,nxm,ny,nz,npaire)
!
!********************************************************************

USE param
USE derivX

implicit none

integer :: nx,nxm,ny,nz,npaire
real(mytype), dimension(nxm,ny,nz) :: tx
real(mytype), dimension(nx,ny,nz) :: ux,rx
real(mytype), dimension(ny,nz) :: sx
real(mytype), dimension(nxm) :: cfx6,csx6,cwx6
integer :: i,j,k,nyz

nyz=ny*nz

if (nclx==0) then
   do k=1,nz
   do j=1,ny
      tx(1,j,k)=acix6*(ux(2,j,k)-ux(1  ,j,k))&
           +bcix6*(ux(3,j,k)-ux(nx,j,k))
      rx(1,j,k)=-1.
      tx(2,j,k)=acix6*(ux(3,j,k)-ux(2 ,j,k))&
           +bcix6*(ux(4,j,k)-ux(1,j,k))
      rx(2,j,k)=0.
      do i=3,nx-2
         tx(i,j,k)=acix6*(ux(i+1,j,k)-ux(i,j,k))&
              +bcix6*(ux(i+2,j,k)-ux(i-1,j,k))
         rx(i,j,k)=0.
      enddo
      tx(nx-1,j,k)=acix6*(ux(nx,j,k)-ux(nx-1,j,k))&
           +bcix6*(ux(1 ,j,k)-ux(nx-2,j,k))
      rx(nx-1,j,k)=0.
      tx(nx  ,j,k)=acix6*(ux(1,j,k)-ux(nx,j,k))&
           +bcix6*(ux(2,j,k)-ux(nx-1,j,k))
      rx(nx  ,j,k)=alcaix6
      do i=2,nx
         tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*csx6(i)
         rx(i,j,k)=rx(i,j,k)-rx(i-1,j,k)*csx6(i)
      enddo
      tx(nx,j,k)=tx(nx,j,k)*cwx6(nx)
      rx(nx,j,k)=rx(nx,j,k)*cwx6(nx)
      do i=nx-1,1,-1
         tx(i,j,k)=(tx(i,j,k)-cfx6(i)*tx(i+1,j,k))*cwx6(i)
         rx(i,j,k)=(rx(i,j,k)-cfx6(i)*rx(i+1,j,k))*cwx6(i)
      enddo
      sx(j,k)=(tx(1,j,k)-alcaix6*tx(nx,j,k))/&
           (1.+rx(1,j,k)-alcaix6*rx(nx,j,k))
      do i=1,nx
         tx(i,j,k)=tx(i,j,k)-sx(j,k)*rx(i,j,k)
      enddo
   enddo
   enddo
endif

if ((nclx==1).or.(nclx==2)) then
   if (npaire==1) then
      do k=1,nz
      do j=1,ny
         tx(1,j,k)=acix6*(ux(2,j,k)-ux(1,j,k))&
              +bcix6*(ux(3,j,k)-ux(2,j,k))
         tx(2,j,k)=acix6*(ux(3,j,k)-ux(2,j,k))&
              +bcix6*(ux(4,j,k)-ux(1,j,k))
         do i=3,nxm-2
            tx(i,j,k)=acix6*(ux(i+1,j,k)-ux(i,j,k))&
                 +bcix6*(ux(i+2,j,k)-ux(i-1,j,k))
         enddo
         tx(nxm-1,j,k)=acix6*(ux(nxm,j,k)-ux(nxm-1,j,k))&
              +bcix6*(ux(nx,j,k)-ux(nxm-2,j,k))
         tx(nxm,j,k)=acix6*(ux(nx,j,k)-ux(nxm,j,k))&
              +bcix6*(ux(nxm,j,k)-ux(nxm-1,j,k))
         do i=2,nxm
            tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*csx6(i)
         enddo
         tx(nxm,j,k)=tx(nxm,j,k)*cwx6(nxm)
         do i=nxm-1,1,-1
            tx(i,j,k)=(tx(i,j,k)-cfx6(i)*tx(i+1,j,k))*cwx6(i)
         enddo
      enddo
      enddo
   endif
   if (npaire==0) then
      do k=1,nz
      do j=1,ny
         tx(1,j,k)=acix6*(ux(2,j,k)-ux(1,j,k))&
              +bcix6*(ux(3,j,k)-2.*ux(1,j,k)+ux(2,j,k))
         tx(2,j,k)=acix6*(ux(3,j,k)-ux(2,j,k))&
              +bcix6*(ux(4,j,k)-ux(1,j,k))
         do i=3,nxm-2
            tx(i,j,k)=acix6*(ux(i+1,j,k)-ux(i,j,k))&
                 +bcix6*(ux(i+2,j,k)-ux(i-1,j,k))
         enddo
         tx(nxm-1,j,k)=acix6*(ux(nxm,j,k)-ux(nxm-1,j,k))&
              +bcix6*(ux(nx,j,k)-ux(nxm-2,j,k))
         tx(nxm,j,k)=acix6*(ux(nx,j,k)-ux(nxm,j,k))&
              +bcix6*(2.*ux(nx,j,k)-ux(nxm,j,k)-ux(nxm-1,j,k))
         do i=2,nxm
            tx(i,j,k)=tx(i,j,k)-tx(i-1,j,k)*csx6(i)
         enddo
         tx(nxm,j,k)=tx(nxm,j,k)*cwx6(nxm)
         do i=nxm-1,1,-1
            tx(i,j,k)=(tx(i,j,k)-cfx6(i)*tx(i+1,j,k))*cwx6(i)
         enddo
      enddo
      enddo
   endif
endif

return
end subroutine decx6
