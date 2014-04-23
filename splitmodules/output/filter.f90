subroutine filter()
!
!*********************************************************************

USE param
USE parfiX
USE parfiY
USE parfiZ
USE variables

implicit none

integer  :: i,j,k
real(mytype) :: xcoef

call coefficients()

if (nx>1) then
   xcoef=1./2.
   if (nclx==0) then
      fiffx(1)=0.
      fifx(1)=0.
      ficx(1)=1.+fibex/xcoef
      fibx(1)=fialx
      fibbx(1)=fibex
      fiffx(2)=0.
      fifx(2)=fialx+fialx*xcoef
      ficx(2)=1.+fibex*xcoef
      fibx(2)=fialx
      fibbx(2)=fibex
      do i=3,nx-2
         fiffx(i)=fibex
         fifx(i)=fialx
         ficx(i)=1.
         fibx(i)=fialx
         fibbx(i)=fibex
      enddo
      fiffx(nx-1)=fibex
      fifx(nx-1)=fialx
      ficx(nx-1)=1.+fibex*xcoef
      fibx(nx-1)=fialx+fialx*xcoef
      fibbx(nx-1)=0.
      fiffx(nx)=fibex
      fifx(nx)=fialx
      ficx(nx)=1.+fibex/xcoef
      fibx(nx)=0.
      fibbx(nx)=0.
   endif
   if (nclx==1) then
      fiffx(1)=0.
      fifx(1)=0.
      ficx(1)=1.
      fibx(1)=fialx+fialx
      fibbx(1)=fibex+fibex
      fiffx(2)=0.
      fifx(2)=fialx
      ficx(2)=1.+fibex
      fibx(2)=fialx
      fibbx(2)=fibex
      do i=3,nx-2
         fiffx(i)=fibex
         fifx(i)=fialx
         ficx(i)=1.
         fibx(i)=fialx
         fibbx(i)=fibex
      enddo
      fiffx(nx-1)=fibex
      fifx(nx-1)=fialx
      ficx(nx-1)=1.+fibex
      fibx(nx-1)=fialx
      fibbx(nx-1)=0.
      fiffx(nx)=fibex+fibex
      fifx(nx)=fialx+fialx
      ficx(nx)=1.
      fibx(nx)=0.
      fibbx(nx)=0.
      do i=1,nx
         fiffxp(i)=fiffx(i)
         fifxp(i)=fifx(i)
         ficxp(i)=ficx(i)
         fibxp(i)=fibx(i)
         fibbxp(i)=fibbx(i)
      enddo
      fibx(1)=0.
      fibbx(1)=0.
      ficx(2)=1.-fibex
      ficx(nx-1)=1.-fibex
      fifx(nx)=0.
      fiffx(nx)=0.
   endif
   if (nclx==2) then
      fiffx(1)=0.
      fifx(1)=0.
      ficx(1)=1.
      fibx(1)=0.
      fibbx(1)=0.
      fiffx(2)=0.
      fifx(2)=0.
      ficx(2)=1.
      fibx(2)=0.
      fibbx(2)=0.
      fiffx(3)=0.
      fifx(3)=0.
      ficx(3)=1.
      fibx(3)=0.
      fibbx(3)=0.
      do i=4,nx-3
         fiffx(i)=fibex
         fifx(i)=fialx
         ficx(i)=1.
         fibx(i)=fialx
         fibbx(i)=fibex
      enddo
      fiffx(nx-2)=0.
      fifx(nx-2)=0.
      ficx(nx-2)=1.
      fibx(nx-2)=0.
      fibbx(nx-2)=0.
      fiffx(nx-1)=0.
      fifx(nx-1)=0.
      ficx(nx-1)=1.
      fibx(nx-1)=0.
      fibbx(nx-1)=0.
      fiffx(nx)=0.
      fifx(nx)=0.
      ficx(nx)=1.
      fibx(nx)=0.
      fibbx(nx)=0.
   endif
   call prepare_filtre(fiffx,fifx,ficx,fibx,fibbx,filax,nx)
   if (nclx==1) then
      call prepare_filtre(fiffxp,fifxp,ficxp,fibxp,&
           fibbxp,filaxp,nx)
      do i=1,nx
         fiffxp(i)=1./fiffxp(i)
      enddo
   endif
   if (nclx==0) then
      call cyclix(fiffx,fifx,ficx,fibx,fibbx,&
           filax,fiz1x,fiz2x,nx)
   endif
   do i=1,nx
      fiffx(i)=1./fiffx(i)
   enddo
endif

if (ny>1) then
   xcoef=1./2.
   if (ncly==0) then
      fiffy(1)=0.
      fify(1)=0.
      ficy(1)=1.+fibey/xcoef
      fiby(1)=fialy
      fibby(1)=fibey
      fiffy(2)=0.
      fify(2)=fialy+xcoef*fialy
      ficy(2)=1.+xcoef*fibey
      fiby(2)=fialy
      fibby(2)=fibey
      do j=3,ny-2
         fiffy(j)=fibey
         fify(j)=fialy
         ficy(j)=1.
         fiby(j)=fialy
         fibby(j)=fibey
      enddo
      fiffy(ny-1)=fibey
      fify(ny-1)=fialy
      ficy(ny-1)=1.+xcoef*fibey
      fiby(ny-1)=fialy+xcoef*fialy
      fibby(ny-1)=0.
      fiffy(ny)=fibey
      fify(ny)=fialy
      ficy(ny)=1.+fibey/xcoef
      fiby(ny)=0.
      fibby(ny)=0.
      do j=1,ny
         fiffyp(j)=fiffy(j)
         fifyp(j)=fify(j)
         ficyp(j)=ficy(j)
         fibyp(j)=fiby(j)
         fibbyp(j)=fibby(j)
      enddo
   endif
   if (ncly==1) then
      fiffy(1)=0.
      fify(1)=0.
      ficy(1)=1.
      fiby(1)=fialy+fialy
      fibby(1)=fibey+fibey
      fiffy(2)=0.
      fify(2)=fialy
      ficy(2)=1.+fibey
      fiby(2)=fialy
      fibby(2)=fibey
      do j=3,ny-2
         fiffy(j)=fibey
         fify(j)=fialy
         ficy(j)=1.
         fiby(j)=fialy
         fibby(j)=fibey
      enddo
      fiffy(ny-1)=fibey
      fify(ny-1)=fialy
      ficy(ny-1)=1.+fibey
      fiby(ny-1)=fialy
      fibby(ny-1)=0.
      fiffy(ny)=fibey+fibey
      fify(ny)=fialy+fialy
      ficy(ny)=1.
      fiby(ny)=0.
      fibby(ny)=0.
      do j=1,ny
         fiffyp(j)=fiffy(j)
         fifyp(j)=fify(j)
         ficyp(j)=ficy(j)
         fibyp(j)=fiby(j)
         fibbyp(j)=fibby(j)
      enddo
      fiby(1)=0.
      fibby(1)=0.
      ficy(2)=1.-fibey
      ficy(ny-1)=1.-fibey
      fify(ny)=0.
      fiffy(ny)=0.
   endif
   if (ncly==2) then
      fiffy(1)=0.
      fify(1)=0.
      ficy(1)=1.
      fiby(1)=0.
      fibby(1)=0.
      fiffy(2)=0.
      fify(2)=0.
      ficy(2)=1.
      fiby(2)=0.
      fibby(2)=0.
      fiffy(3)=0.
      fify(3)=0.
      ficy(3)=1.
      fiby(3)=0.
      fibby(3)=0.
      do j=4,ny-3
         fiffy(j)=fibey
         fify(j)=fialy
         ficy(j)=1.
         fiby(j)=fialy
         fibby(j)=fibey
      enddo
      fiffy(ny-2)=0.
      fify(ny-2)=0.
      ficy(ny-2)=1.
      fiby(ny-2)=0.
      fibby(ny-2)=0.
      fiffy(ny-1)=0.
      fify(ny-1)=0.
      ficy(ny-1)=1.
      fiby(ny-1)=0.
      fibby(ny-1)=0.
      fiffy(ny)=0.
      fify(ny)=0.
      ficy(ny)=1.
      fiby(ny)=0.
      fibby(ny)=0.
   endif
   call prepare_filtre(fiffy,fify,ficy,fiby,fibby,filay,ny)
   if (ncly==1) then
      call prepare_filtre(fiffyp,fifyp,ficyp,fibyp,&
           fibbyp,filayp,ny)
      do j=1,ny
         fiffyp(j)=1./fiffyp(j)
      enddo
   endif
   if (ncly==0) then
      call cycliy(fiffy,fify,ficy,fiby,fibby,&
           filay,fiz1y,fiz2y,ny)
   endif
   do j=1,ny
      fiffy(j)=1./fiffy(j)
   enddo
endif

#ifndef TWOD
   xcoef=1./2.
   if (nclz==0) then
      fiffz(1)=0.
      fifz(1)=0.
      ficz(1)=1.+fibez/xcoef
      fibz(1)=fialz
      fibbz(1)=fibez
      fiffz(2)=0.
      fifz(2)=fialz+fialz*xcoef
      ficz(2)=1.+fibez*xcoef
      fibz(2)=fialz
      fibbz(2)=fibez
      do k=3,nz-2
         fiffz(k)=fibez
         fifz(k)=fialz
         ficz(k)=1.
         fibz(k)=fialz
         fibbz(k)=fibez
      enddo
      fiffz(nz-1)=fibez
      fifz(nz-1)=fialz
      ficz(nz-1)=1.+fibez*xcoef
      fibz(nz-1)=fialz+fialz*xcoef
      fibbz(nz-1)=0.
      fiffz(nz)=fibez
      fifz(nz)=fialz
      ficz(nz)=1.+fibez/xcoef
      fibz(nz)=0.
      fibbz(nz)=0.
      do k=1,nz
         fiffzp(k)=fiffz(k)
         fifzp(k)=fifz(k)
         ficzp(k)=ficz(k)
         fibzp(k)=fibz(k)
         fibbzp(k)=fibbz(k)
      enddo
   endif
   if (nclz==1) then
      fiffz(1)=0.
      fifz(1)=0.
      ficz(1)=1.
      fibz(1)=fialz+fialz
      fibbz(1)=fibez+fibez
      fiffz(2)=0.
      fifz(2)=fialz
      ficz(2)=1.+fibez
      fibz(2)=fialz
      fibbz(2)=fibez
      do k=3,nz-2
         fiffz(k)=fibez
         fifz(k)=fialz
         ficz(k)=1.
         fibz(k)=fialz
         fibbz(k)=fibez
      enddo
      fiffz(nz-1)=fibez
      fifz(nz-1)=fialz
      ficz(nz-1)=1.+fibez
      fibz(nz-1)=fialz
      fibbz(nz-1)=0.
      fiffz(nz)=fibez+fibez
      fifz(nz)=fialz+fialz
      ficz(nz)=1.
      fibz(nz)=0.
      fibbz(nz)=0.
      do k=1,nz
         fiffzp(k)=fiffz(k)
         fifzp(k)=fifz(k)
         ficzp(k)=ficz(k)
         fibzp(k)=fibz(k)
         fibbzp(k)=fibbz(k)
      enddo
      fibz(1)=0.
      fibbz(1)=0.
      ficz(2)=1.-fibez
      ficz(nz-1)=1.-fibez
      fifz(nz)=0.
      fiffz(nz)=0.
   endif
   if (nclz==2) then
      fiffz(1)=0.
      fifz(1)=0.
      ficz(1)=1.
      fibz(1)=0.
      fibbz(1)=0.
      fiffz(2)=0.
      fifz(2)=0.
      ficz(2)=1.
      fibz(2)=0.
      fibbz(2)=0.
      fiffz(3)=0.
      fifz(3)=0.
      ficz(3)=1.
      fibz(3)=0.
      fibbz(3)=0.
      do k=4,nz-3
         fiffz(k)=fibez
         fifz(k)=fialz
         ficz(k)=1.
         fibz(k)=fialz
         fibbz(k)=fibez
      enddo
      fiffz(nz-2)=0.
      fifz(nz-2)=0.
      ficz(nz-2)=1.
      fibz(nz-2)=0.
      fibbz(nz-2)=0.
      fiffz(nz-1)=0.
      fifz(nz-1)=0.
      ficz(nz-1)=1.
      fibz(nz-1)=0.
      fibbz(nz-1)=0.
      fiffz(nz)=0.
      fifz(nz)=0.
      ficz(nz)=1.
      fibz(nz)=0.
      fibbz(nz)=0.
   endif
   call prepare_filtre(fiffz,fifz,ficz,fibz,fibbz,&
        filaz,nz)
   if (nclz==1) then
      call prepare_filtre(fiffzp,fifzp,ficzp,fibzp,&
           fibbzp,filazp,nz)
      do k=1,nz
         fiffzp(k)=1./fiffzp(k)
      enddo
   endif
   if (nclz==0) then
      call cycliz(fiffz,fifz,ficz,fibz,fibbz,&
           filaz,fiz1z,fiz2z,nz)
   endif
   do k=1,nz
      fiffz(k)=1./fiffz(k)
   enddo
#endif

return
end subroutine filter
