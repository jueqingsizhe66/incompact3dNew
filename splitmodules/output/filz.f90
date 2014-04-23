subroutine filz(tz,uz,rz,sz,vz,fiffz,fifz,ficz,fibz,fibbz,filaz,fiz1z,&
     fiz2z,nx,ny,nz,npaire)
!
!*********************************************************************

USE param
USE parfiZ

implicit none

integer :: nx,ny,nz,npaire,i,j,k
real(mytype), dimension(nx,ny,nz) :: tz,uz,rz
real(mytype), dimension(nx,ny) :: sz,vz
real(mytype), dimension(nz) :: fiffz,fifz,ficz,fibz,fibbz,fiz1z,fiz2z
real(mytype), dimension(nz,2) :: filaz
real(mytype) :: xcoef

if (nclz==0) then
   do j=1,ny
   do i=1,nx
      rz(i,j,1)=fiaiz*uz(i,j,1)+&
           fibiz*(uz(i,j,2)+uz(i,j,nz))+&
           ficiz*(uz(i,j,3)+uz(i,j,nz-1))+&
           fidiz*(uz(i,j,4)+uz(i,j,nz-2))
      rz(i,j,2)=fiaiz*uz(i,j,2)+&
           fibiz*(uz(i,j,3)+uz(i,j,1))+ &
           ficiz*(uz(i,j,4)+uz(i,j,nz))+&
           fidiz*(uz(i,j,5)+uz(i,j,nz-1))
      rz(i,j,3)=fiaiz*uz(i,j,3)+&
           fibiz*(uz(i,j,4)+uz(i,j,2))+&
           ficiz*(uz(i,j,5)+uz(i,j,1))+&
           fidiz*(uz(i,j,6)+uz(i,j,nz))
   enddo
   enddo
   do k=4,nz-3
   do j=1,ny
   do i=1,nx
      rz(i,j,k)=fiaiz*uz(i,j,k)+&
           fibiz*(uz(i,j,k+1)+uz(i,j,k-1))+&
           ficiz*(uz(i,j,k+2)+uz(i,j,k-2))+&
           fidiz*(uz(i,j,k+3)+uz(i,j,k-3))
   enddo
   enddo
   enddo
   do j=1,ny
   do i=1,nx
      rz(i,j,nz)=fiaiz*uz(i,j,nz)+&
           fibiz*(uz(i,j,1)+uz(i,j,nz-1))+&
           ficiz*(uz(i,j,2)+uz(i,j,nz-2))+&
           fidiz*(uz(i,j,3)+uz(i,j,nz-3))
      rz(i,j,nz-1)=fiaiz*uz(i,j,nz-1)+&
           fibiz*(uz(i,j,nz)+uz(i,j,nz-2))+&
           ficiz*(uz(i,j,1)+uz(i,j,nz-3))+&
           fidiz*(uz(i,j,2)+uz(i,j,nz-4))
      rz(i,j,nz-2)=fiaiz*uz(i,j,nz-2)+&
           fibiz*(uz(i,j,nz-1)+uz(i,j,nz-3))+&
           ficiz*(uz(i,j,nz)+uz(i,j,nz-4))+&
           fidiz*(uz(i,j,1)+uz(i,j,nz-5))
   enddo
   enddo
   do j=1,ny
   do i=1,nx
   do k=1,nz-2
      rz(i,j,k+1)=rz(i,j,k+1)-filaz(k,1)*rz(i,j,k)
      rz(i,j,k+2)=rz(i,j,k+2)-filaz(k,2)*rz(i,j,k)
   enddo
   rz(i,j,nz)=rz(i,j,nz)-filaz(nz-1,1)*rz(i,j,nz-1)
   enddo
   enddo
   do j=1,ny
   do i=1,nx
      rz(i,j,nz)=rz(i,j,nz)*fiffz(nz)
      rz(i,j,nz-1)=(rz(i,j,nz-1)-fifz(nz-1)*rz(i,j,nz))*&
           fiffz(nz-1)
      rz(i,j,nz-2)=(rz(i,j,nz-2)-fifz(nz-2)*rz(i,j,nz-1)-&
           ficz(nz-2)*rz(i,j,nz))*fiffz(nz-2)
      rz(i,j,nz-3)=(rz(i,j,nz-3)-fifz(nz-3)*rz(i,j,nz-2)-&
           ficz(nz-3)*rz(i,j,nz-1)-&
           fibz(nz-3)*rz(i,j,nz))*fiffz(nz-3)
   enddo
   enddo
   do j=1,ny
   do i=1,nx
   do k=nz-4,1,-1
      rz(i,j,k)=(rz(i,j,k)-fifz(k)*rz(i,j,k+1)-&
           ficz(k)*rz(i,j,k+2)-&
           fibz(k)*rz(i,j,k+3)-&
           fibbz(k)*rz(i,j,k+4))*fiffz(k)
   enddo
   enddo
   enddo
   xcoef=1./2.
   do j=1,ny
   do i=1,nx
      sz(i,j)=fih1z*(-fibez*rz(i,j,1)+fibez*rz(i,j,nz-1)*xcoef+&
           fialz*rz(i,j,nz)*xcoef)+&
           fih2z*(fialz*rz(i,j,1)*xcoef+fibez*rz(i,j,2)*xcoef-&
           fibez*rz(i,j,nz))
      vz(i,j)=fih3z*(-fibez*rz(i,j,1)+fibez*rz(i,j,nz-1)*xcoef+&
           fialz*rz(i,j,nz)*xcoef)+&
           fih4z*(fialz*rz(i,j,1)*xcoef+fibez*rz(i,j,2)*xcoef-&
           fibez*rz(i,j,nz))
   enddo
   enddo
   do k=1,nz
   do j=1,ny
   do i=1,nx
      tz(i,j,k)=rz(i,j,k)-fiz1z(k)*sz(i,j)-fiz2z(k)*vz(i,j)
   enddo
   enddo
   enddo
endif

if (nclz==1) then
   if (npaire==1) then
      do j=1,ny
      do i=1,nx
         tz(i,j,1)=fiaiz*uz(i,j,1)+&
              fibiz*(uz(i,j,2)+uz(i,j,2))+&
              ficiz*(uz(i,j,3)+uz(i,j,3))+&
              fidiz*(uz(i,j,4)+uz(i,j,4))
         tz(i,j,2)=fiaiz*uz(i,j,2)+&
              fibiz*(uz(i,j,3)+uz(i,j,1))+&
              ficiz*(uz(i,j,4)+uz(i,j,2))+&
              fidiz*(uz(i,j,5)+uz(i,j,3))
         tz(i,j,3)=fiaiz*uz(i,j,3)+&
              fibiz*(uz(i,j,4)+uz(i,j,2))+&
              ficiz*(uz(i,j,5)+uz(i,j,1))+&
              fidiz*(uz(i,j,6)+uz(i,j,2))
      enddo
      enddo
      do k=4,nz-3
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=fiaiz*uz(i,j,k)+&
              fibiz*(uz(i,j,k+1)+uz(i,j,k-1))+&
              ficiz*(uz(i,j,k+2)+uz(i,j,k-2))+&
              fidiz*(uz(i,j,k+3)+uz(i,j,k-3))
      enddo
      enddo
      enddo
      do j=1,ny
      do i=1,nx
         tz(i,j,nz)=fiaiz*uz(i,j,nz)+&
              fibiz*(uz(i,j,nz-1)+uz(i,j,nz-1))+&
              ficiz*(uz(i,j,nz-2)+uz(i,j,nz-2))+&
              fidiz*(uz(i,j,nz-3)+uz(i,j,nz-3))
         tz(i,j,nz-1)=fiaiz*uz(i,j,nz-1)+&
              fibiz*(uz(i,j,nz)+uz(i,j,nz-2))+&
              ficiz*(uz(i,j,nz-1)+uz(i,j,nz-3))+&
              fidiz*(uz(i,j,nz-2)+uz(i,j,nz-4))
         tz(i,j,nz-2)=fiaiz*uz(i,j,nz-2)+&
              fibiz*(uz(i,j,nz-1)+uz(i,j,nz-3))+&
              ficiz*(uz(i,j,nz)+uz(i,j,nz-4))+&
              fidiz*(uz(i,j,nz-1)+uz(i,j,nz-5))
      enddo
      enddo
      do j=1,ny
      do i=1,nx
      do k=1,nz-2
         tz(i,j,k+1)=tz(i,j,k+1)-filaz(k,1)*tz(i,j,k)
         tz(i,j,k+2)=tz(i,j,k+2)-filaz(k,2)*tz(i,j,k)
      enddo
      tz(i,j,nz)=tz(i,j,nz)-filaz(nz-1,1)*tz(i,j,nz-1)
      enddo
      enddo
      do j=1,ny
      do i=1,nx
         tz(i,j,nz)=tz(i,j,nz)*fiffz(nz)
         tz(i,j,nz-1)=(tz(i,j,nz-1)-fifz(nz-1)*tz(i,j,nz))*&
              fiffz(nz-1)
         tz(i,j,nz-2)=(tz(i,j,nz-2)-fifz(nz-2)*tz(i,j,nz-1)-&
              ficz(nz-2)*tz(i,j,nz))*fiffz(nz-2)
         tz(i,j,nz-3)=(tz(i,j,nz-3)-fifz(nz-3)*tz(i,j,nz-2)-&
              ficz(nz-3)*tz(i,j,nz-1)-&
              fibz(nz-3)*tz(i,j,nz))*fiffz(nz-3)
         do k=nz-4,1,-1
            tz(i,j,k)=(tz(i,j,k)-fifz(k)*tz(i,j,k+1)-&
                 ficz(k)*tz(i,j,k+2)-&
                 fibz(k)*tz(i,j,k+3)-&
                 fibbz(k)*tz(i,j,k+4))*fiffz(k)
         enddo
      enddo
      enddo
   endif
   if (npaire==0) then
      do j=1,ny
      do i=1,nx
         tz(i,j,1)=fiaiz*uz(i,j,1)+&
              fibiz*(uz(i,j,2)-uz(i,j,2))+&
              ficiz*(uz(i,j,3)-uz(i,j,3))+&
              fidiz*(uz(i,j,4)-uz(i,j,4))
         tz(i,j,2)=fiaiz*uz(i,j,2)+&
              fibiz*(uz(i,j,3)+uz(i,j,1))+&
              ficiz*(uz(i,j,4)-uz(i,j,2))+&
              fidiz*(uz(i,j,5)-uz(i,j,3))
         tz(i,j,3)=fiaiz*uz(i,j,3)+&
              fibiz*(uz(i,j,4)+uz(i,j,2))+&
              ficiz*(uz(i,j,5)+uz(i,j,1))+&
              fidiz*(uz(i,j,6)-uz(i,j,2))
      enddo
      enddo
      do k=4,nz-3
      do j=1,ny
      do i=1,nx
         tz(i,j,k)=fiaiz*uz(i,j,k)+&
              fibiz*(uz(i,j,k+1)+uz(i,j,k-1))+&
              ficiz*(uz(i,j,k+2)+uz(i,j,k-2))+&
              fidiz*(uz(i,j,k+3)+uz(i,j,k-3))
      enddo
      enddo
      enddo
      do j=1,ny
      do i=1,nx
         tz(i,j,nz)=fiaiz*uz(i,j,nz)+&
              fibiz*(uz(i,j,nz-1)-uz(i,j,nz-1))+&
              ficiz*(uz(i,j,nz-2)-uz(i,j,nz-2))+&
              fidiz*(uz(i,j,nz-3)-uz(i,j,nz-3))
         tz(i,j,nz-1)=fiaiz*uz(i,j,nz-1)+&
              fibiz*(uz(i,j,nz)+uz(i,j,nz-2))+&
              ficiz*(-uz(i,j,nz-1)+uz(i,j,nz-3))+&
              fidiz*(-uz(i,j,nz-2)+uz(i,j,nz-4))
         tz(i,j,nz-2)=fiaiz*uz(i,j,nz-2)+&
              fibiz*(uz(i,j,nz-1)+uz(i,j,nz-3))+&
              ficiz*(uz(i,j,nz)+uz(i,j,nz-4))+&
              fidiz*(-uz(i,j,nz-1)+uz(i,j,nz-5))
      enddo
      enddo
      do j=1,ny
      do i=1,nx
      do k=1,nz-2
         tz(i,j,k+1)=tz(i,j,k+1)-filaz(k,1)*tz(i,j,k)
         tz(i,j,k+2)=tz(i,j,k+2)-filaz(k,2)*tz(i,j,k)
      enddo
      tz(i,j,nz)=tz(i,j,nz)-filaz(nz-1,1)*tz(i,j,nz-1)
      enddo
      enddo
      do j=1,ny
      do i=1,nx
         tz(i,j,nz)=tz(i,j,nz)*fiffz(nz)
         tz(i,j,nz-1)=(tz(i,j,nz-1)-fifz(nz-1)*tz(i,j,nz))*&
              fiffz(nz-1)
         tz(i,j,nz-2)=(tz(i,j,nz-2)-fifz(nz-2)*tz(i,j,nz-1)-&
              ficz(nz-2)*tz(i,j,nz))*fiffz(nz-2)
         tz(i,j,nz-3)=(tz(i,j,nz-3)-fifz(nz-3)*tz(i,j,nz-2)-&
              ficz(nz-3)*tz(i,j,nz-1)-&
              fibz(nz-3)*tz(i,j,nz))*fiffz(nz-3)
         do k=nz-4,1,-1
            tz(i,j,k)=(tz(i,j,k)-fifz(k)*tz(i,j,k+1)-&
                 ficz(k)*tz(i,j,k+2)-&
                 fibz(k)*tz(i,j,k+3)-&
                 fibbz(k)*tz(i,j,k+4))*fiffz(k)
         enddo
      enddo
      enddo
   endif
endif

if (nclz==2) then
   do j=1,ny
   do i=1,nx
      tz(i,j,1)=fia1z*uz(i,j,1)+fib1z*uz(i,j,2)+&
           fic1z*uz(i,j,3)+fid1z*uz(i,j,4)+&
           fie1z*uz(i,j,5)
      tz(i,j,2)=fia2z*uz(i,j,2)+fib2z*uz(i,j,1)+&
           fic2z*uz(i,j,3)+fid2z*uz(i,j,4)+&
           fie2z*uz(i,j,5)
      tz(i,j,3)=fia3z*uz(i,j,3)+fib3z*uz(i,j,1)+&
           fic3z*uz(i,j,2)+fid3z*uz(i,j,4)+&
           fie3z*uz(i,j,5)
   enddo
   enddo
   do k=4,nz-3
   do j=1,ny
   do i=1,nx
      tz(i,j,k)=fiaiz*uz(i,j,k)+&
           fibiz*(uz(i,j,k+1)+uz(i,j,k-1))+&
           ficiz*(uz(i,j,k+2)+uz(i,j,k-2))+&
           fidiz*(uz(i,j,k+3)+uz(i,j,k-3))
   enddo
   enddo
   enddo
   do j=1,ny
   do i=1,nx
      tz(i,j,nz)=fianz*uz(i,j,nz)+fibnz*uz(i,j,nz-1)+&
           ficnz*uz(i,j,nz-2)+fidnz*uz(i,j,nz-3)+&
           fienz*uz(i,j,nz-4)
      tz(i,j,nz-1)=fiamz*uz(i,j,nz-1)+fibmz*uz(i,j,nz)+&
           ficmz*uz(i,j,nz-2)+fidmz*uz(i,j,nz-3)+&
           fiemz*uz(i,j,nz-4)
      tz(i,j,nz-2)=fiapz*uz(i,j,nz-2)+fibpz*uz(i,j,nz)+&
           ficpz*uz(i,j,nz-1)+fidpz*uz(i,j,nz-3)+&
           fiepz*uz(i,j,nz-4)
   enddo
   enddo
   do j=1,ny
   do i=1,nx
   do k=1,nz-2
      tz(i,j,k+1)=tz(i,j,k+1)-filaz(k,1)*tz(i,j,k)
      tz(i,j,k+2)=tz(i,j,k+2)-filaz(k,2)*tz(i,j,k)
   enddo
   tz(i,j,nz)=tz(i,j,nz)-filaz(nz-1,1)*tz(i,j,nz-1)
   enddo
   enddo
   do j=1,ny
   do i=1,nx
      tz(i,j,nz)=tz(i,j,nz)*fiffz(nz)
      tz(i,j,nz-1)=(tz(i,j,nz-1)-fifz(nz-1)*tz(i,j,nz))*&
           fiffz(nz-1)
      tz(i,j,nz-2)=(tz(i,j,nz-2)-fifz(nz-2)*tz(i,j,nz-1)-&
           ficz(nz-2)*tz(i,j,nz))*fiffz(nz-2)
      tz(i,j,nz-3)=(tz(i,j,nz-3)-fifz(nz-3)*tz(i,j,nz-2)-&
           ficz(nz-3)*tz(i,j,nz-1)-&
           fibz(nz-3)*tz(i,j,nz))*fiffz(nz-3)
      do k=nz-4,1,-1
         tz(i,j,k)=(tz(i,j,k)-fifz(k)*tz(i,j,k+1)-&
              ficz(k)*tz(i,j,k+2)-&
              fibz(k)*tz(i,j,k+3)-&
              fibbz(k)*tz(i,j,k+4))*fiffz(k)
      enddo
   enddo
   enddo
endif

return
end subroutine filz
