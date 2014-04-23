subroutine interiz6(tz,uz,rz,sz,cifiz6,cisiz6,ciwiz6,cifz6,cisz6,ciwz6,&
     nx,ny,nzm,nz,npaire)
!
!********************************************************************

USE param
USE derivZ

implicit none

integer :: nx,ny,nz,nzm,npaire
real(mytype), dimension(nx,ny,nz) :: tz
real(mytype), dimension(nx,ny,nzm) :: uz
real(mytype), dimension(nx,ny,nz) :: rz
real(mytype), dimension(nx,ny) :: sz
real(mytype), dimension(nz) :: cifiz6,cisiz6,ciwiz6
real(mytype), dimension(nz) :: cifz6,cisz6,ciwz6
integer :: i,j,k

if (nclz==0) then
   do j=1,ny
   do i=1,nx
      tz(i,j,1)=aiciz6*(uz(i,j,1)+uz(i,j,nz))&
           +biciz6*(uz(i,j,2)+uz(i,j,nz-1))&
           +ciciz6*(uz(i,j,3)+uz(i,j,nz-2))
      rz(i,j,1)=-1.
      tz(i,j,2)=aiciz6*(uz(i,j,2)+uz(i,j,1))&
           +biciz6*(uz(i,j,3)+uz(i,j,nz))&
           +ciciz6*(uz(i,j,4)+uz(i,j,nz-1))
      rz(i,j,2)=0.
      tz(i,j,3)=aiciz6*(uz(i,j,3)+uz(i,j,2))&
           +biciz6*(uz(i,j,4)+uz(i,j,1))&
           +ciciz6*(uz(i,j,5)+uz(i,j,nz))
      rz(i,j,3)=0.
   enddo
   enddo
   do k=4,nz-2
   do j=1,ny
   do i=1,nx
      tz(i,j,k)=aiciz6*(uz(i,j,k)+uz(i,j,k-1))&
           +biciz6*(uz(i,j,k+1)+uz(i,j,k-2))&
           +ciciz6*(uz(i,j,k+2)+uz(i,j,k-3))
      rz(i,j,k)=0.
   enddo
   enddo
   enddo
   do j=1,ny
   do i=1,nx
      tz(i,j,nz-1)=aiciz6*(uz(i,j,nz-1)+uz(i,j,nz-2))&
           +biciz6*(uz(i,j,nz)+uz(i,j,nz-3))&
           +ciciz6*(uz(i,j,1)+uz(i,j,nz-4))
      rz(i,j,nz-1)=0.
      tz(i,j,nz)=aiciz6*(uz(i,j,nz)+uz(i,j,nz-1))&
           +biciz6*(uz(i,j,1)+uz(i,j,nz-2))&
           +ciciz6*(uz(i,j,2)+uz(i,j,nz-3))
      rz(i,j,nz)=ailcaiz6
   enddo
   enddo
   do k=2,nz
   do j=1,ny
   do i=1,nx
      tz(i,j,k)=tz(i,j,k)-tz(i,j,k-1)*cisz6(k)
      rz(i,j,k)=rz(i,j,k)-rz(i,j,k-1)*cisz6(k)
   enddo
   enddo
   enddo
   do j=1,ny
   do i=1,nx
      tz(i,j,nz)=tz(i,j,nz)*ciwz6(nz)
      rz(i,j,nz)=rz(i,j,nz)*ciwz6(nz)
   enddo
   enddo
   do k=nz-1,1,-1
   do j=1,ny
   do i=1,nx
      tz(i,j,k)=(tz(i,j,k)-cifz6(k)*tz(i,j,k+1))*ciwz6(k)
      rz(i,j,k)=(rz(i,j,k)-cifz6(k)*rz(i,j,k+1))*ciwz6(k)
   enddo
   enddo
   enddo
   do j=1,ny
   do i=1,nx
      sz(i,j)=(tz(i,j,1)-ailcaiz6*tz(i,j,nz))/&
           (1.+rz(i,j,1)-ailcaiz6*rz(i,j,nz))
   enddo
   enddo
   do k=1,nz
   do j=1,ny
   do i=1,nx
      tz(i,j,k)=tz(i,j,k)-sz(i,j)*rz(i,j,k)
   enddo
   enddo
   enddo
endif

if ((nclz==1).or.(nclz==2)) then
   if (npaire==1) then
      do j=1,ny
      do i=1,nx
         tz(i,j,1)=aiciz6*(uz(i,j,1)+uz(i,j,1))&
              +biciz6*(uz(i,j,2)+uz(i,j,2))&
              +ciciz6*(uz(i,j,3)+uz(i,j,3))
         tz(i,j,2)=aiciz6*(uz(i,j,2)+uz(i,j,1))&
              +biciz6*(uz(i,j,3)+uz(i,j,1))&
              +ciciz6*(uz(i,j,4)+uz(i,j,2))
         tz(i,j,3)=aiciz6*(uz(i,j,3)+uz(i,j,2))&
              +biciz6*(uz(i,j,4)+uz(i,j,1))&
              +ciciz6*(uz(i,j,5)+uz(i,j,1))
      enddo
      enddo
      do k=4,nz-3
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=aiciz6*(uz(i,j,k)+uz(i,j,k-1))&
              +biciz6*(uz(i,j,k+1)+uz(i,j,k-2))&
              +ciciz6*(uz(i,j,k+2)+uz(i,j,k-3))
      enddo
      enddo
      enddo
      do j=1,ny
      do i=1,nx
         tz(i,j,nz-2)=aiciz6*(uz(i,j,nz-2)+uz(i,j,nz-3))&
              +biciz6*(uz(i,j,nz-1)+uz(i,j,nz-4))&
              +ciciz6*(uz(i,j,nz-1)+uz(i,j,nz-5))
         tz(i,j,nz-1)=aiciz6*(uz(i,j,nz-1)+uz(i,j,nz-2))&
              +biciz6*(uz(i,j,nz-1)+uz(i,j,nz-3))&
              +ciciz6*(uz(i,j,nz-2)+uz(i,j,nz-4))
         tz(i,j,nz)=aiciz6*(uz(i,j,nz-1)+uz(i,j,nz-1))&
              +biciz6*(uz(i,j,nz-2)+uz(i,j,nz-2))&
              +ciciz6*(uz(i,j,nz-3)+uz(i,j,nz-3))
      enddo
      enddo
      do k=2,nz
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=tz(i,j,k)-tz(i,j,k-1)*cisiz6(k)
      enddo
      enddo
      enddo
      do j=1,ny
      do i=1,nx
         tz(i,j,nz)=tz(i,j,nz)*ciwiz6(nz)
      enddo
      enddo
      do k=nz-1,1,-1
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=(tz(i,j,k)-cifiz6(k)*tz(i,j,k+1))*ciwiz6(k)
      enddo
      enddo
      enddo
   endif
endif

return
end subroutine interiz6
