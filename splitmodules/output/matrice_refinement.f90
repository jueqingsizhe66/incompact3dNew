subroutine matrice_refinement()
!
!**************************************************************************

USE decomp_2d
USE variables
USE param
USE var
USE MPI
USE derivX
USE derivY
USE derivZ

implicit none

integer :: i,j,k

complex(mytype),dimension(sp%yst(1):sp%yen(1)) :: transx
complex(mytype),dimension(sp%yst(2):sp%yen(2)) :: transy
complex(mytype),dimension(sp%yst(3):sp%yen(3)) :: transz
real(mytype) :: xa0,xa1
complex(mytype) :: ytt,xtt,ztt,yt1,xt1,yt2,xt2
complex(mytype) :: xtt1,ytt1,ztt1,zt1,zt2,tmp1,tmp2,tmp3

do i = sp%yst(1),sp%yen(1)
   xtt=cmplx((bicix6*2.*cos(real(exs(i), kind=mytype)*3.*dx/2.)+&
        cicix6*2.*cos(real(exs(i), kind=mytype)*5.*dx/2.)),&
        (bicix6*2.*cos(real(exs(i), kind=mytype)*3.*dx/2.)+&
        cicix6*2.*cos(real(exs(i), kind=mytype)*5.*dx/2.)), kind=mytype)
   xtt1=cmplx((aicix6*2.*cos(real(exs(i), kind=mytype)*dx/2.)),&
        (aicix6*2.*cos(real(exs(i), kind=mytype)*dx/2.)), kind=mytype)
   xt1=cmplx((1.+2.*ailcaix6*cos(real(exs(i), kind=mytype)*dx)),&
        (1.+2.*ailcaix6*cos(real(exs(i), kind=mytype)*dx)), kind=mytype)
   transx(i)=cmplx(real(xtt1+xtt, kind=mytype)/real(xt1, kind=mytype),&
        real(xtt1+xtt, kind=mytype)/real(xt1, kind=mytype), kind=mytype)!(xtt+xtt)/xt1
enddo
do j = sp%yst(2),sp%yen(2)
   ytt=cmplx((biciy6*2.*cos(real(eys(j), kind=mytype)*3.*dy/2.)+&
        ciciy6*2.*cos(real(eys(j), kind=mytype)*5.*dy/2.)),&
        (biciy6*2.*cos(real(eys(j), kind=mytype)*3.*dy/2.)+&
        ciciy6*2.*cos(real(eys(j), kind=mytype)*5.*dy/2.)), kind=mytype)
   ytt1=cmplx((aiciy6*2.*cos(real(eys(j), kind=mytype)*dy/2.)),&
        (aiciy6*2.*cos(real(eys(j), kind=mytype)*dy/2.)), kind=mytype)
   yt1=cmplx((1.+2.*ailcaiy6*cos(real(eys(j), kind=mytype)*dy)),&
        (1.+2.*ailcaiy6*cos(real(eys(j), kind=mytype)*dy)), kind=mytype)
   transy(j)=cmplx(real(ytt1+ytt, kind=mytype)/real(yt1, kind=mytype),&
        real(ytt1+ytt, kind=mytype)/real(yt1, kind=mytype), kind=mytype)!(ytt+ytt)/yt1
enddo
if (nclz==0) then
   do k = sp%yst(3),sp%yen(3)
      ztt=cmplx((biciz6*2.*cos(real(ezs(k), kind=mytype)*3.*dz/2.)+&
           ciciz6*2.*cos(real(ezs(k), kind=mytype)*5.*dz/2.)),&
           (biciz6*2.*cos(real(ezs(k), kind=mytype)*3.*dz/2.)+&
           ciciz6*2.*cos(real(ezs(k), kind=mytype)*5.*dz/2.)), kind=mytype)
      ztt1=cmplx((aiciz6*2.*cos(real(ezs(k), kind=mytype)*dz/2.)),&
           (aiciz6*2.*cos(real(ezs(k), kind=mytype)*dz/2.)), kind=mytype)
      zt1=cmplx((1.+2.*ailcaiz6*cos(real(ezs(k), kind=mytype)*dz)),&
           (1.+2.*ailcaiz6*cos(real(ezs(k), kind=mytype)*dz)), kind=mytype)
      transz(k)=cmplx(real(ztt1+ztt, kind=mytype)/real(zt1, kind=mytype),&
           aimag(ztt1+ztt)/aimag(zt1), kind=mytype)!(ztt+ztt)/zt1
   enddo
else
   do k = sp%yst(3),sp%yen(3)
      ztt=cmplx((biciz6*2.*cos(real(ezs(k), kind=mytype)*3.*dz/2.)+&
           ciciz6*2.*cos(real(ezs(k), kind=mytype)*5.*dz/2.)),&
           (biciz6*2.*cos(aimag(ezs(k))*3.*dz/2.)+&
           ciciz6*2.*cos(aimag(ezs(k))*5.*dz/2.)), kind=mytype)
      ztt1=cmplx((aiciz6*2.*cos(real(ezs(k), kind=mytype)*dz/2.)),&
           (aiciz6*2.*cos(aimag(ezs(k))*dz/2.)), kind=mytype)
      zt1=cmplx((1.+2.*ailcaiz6*cos(real(ezs(k), kind=mytype)*dz)),&
           (1.+2.*ailcaiz6*cos(aimag(ezs(k))*dz)), kind=mytype)
      transz(k)=cmplx(real(ztt1+ztt, kind=mytype)/real(zt1, kind=mytype),&
           aimag(ztt1+ztt)/aimag(zt1), kind=mytype)!(ztt+ztt)/zt1
   enddo
endif

if ((istret==1).or.(istret==2)) then


   xa0=alpha/pi+1./2./beta/pi
   if (istret==1) xa1=1./4./beta/pi
   if (istret==2) xa1=-1./4./beta/pi
!
!construction of the pentadiagonal matrice
!
   do k=sp%yst(3),sp%yen(3)
   do j=1,ny/2
   do i=sp%yst(1),sp%yen(1)
      cw22(i,j,k)=cmplx(real(yky(2*j-1), kind=mytype)*real(transx(i), kind=mytype)*real(transz(k), k
           aimag(yky(2*j-1))*aimag(transx(i))*aimag(transz(k)), kind=mytype)
      cw2(i,j,k)=cmplx(real(yky(2*j), kind=mytype)*real(transx(i), kind=mytype)*real(transz(k), kind
           aimag(yky(2*j))*aimag(transx(i))*aimag(transz(k)), kind=mytype)
   enddo
   enddo
   enddo




!main diagonal
   do k=sp%yst(3),sp%yen(3)
   do j=2,ny/2-1
   do i=sp%yst(1),sp%yen(1)
      a(i,j,k,3)=cmplx(-(real(xk2(i), kind=mytype)*real(transy(2*j-1), kind=mytype)*real(transy(2*j-
           *real(transz(k), kind=mytype)*real(transz(k), kind=mytype))&
           -(real(zk2(k), kind=mytype)*real(transy(2*j-1), kind=mytype)*real(transy(2*j-1), kind=myt
           real(transx(i), kind=mytype)*real(transx(i), kind=mytype))&
           -real(cw22(i,j,k), kind=mytype)*real(cw22(i,j,k), kind=mytype)*xa0*xa0-&
           xa1*xa1*(real(cw22(i,j,k), kind=mytype)*real(cw22(i,j-1,k), kind=mytype)+real(cw22(i,j,k)
           real(cw22(i,j+1,k), kind=mytype)),&
           -(aimag(xk2(i))*aimag(transy(2*j-1))*aimag(transy(2*j-1))*aimag(transz(k))*aimag(transz(k
           -(aimag(zk2(k))*aimag(transy(2*j-1))*aimag(transy(2*j-1))*aimag(transx(i))*aimag(transx(i
           -aimag(cw22(i,j,k))*aimag(cw22(i,j,k))*xa0*xa0-&
           xa1*xa1*(aimag(cw22(i,j,k))*aimag(cw22(i,j-1,k))+aimag(cw22(i,j,k))*aimag(cw22(i,j+1,k)))
      a2(i,j,k,3)=cmplx(-(real(xk2(i), kind=mytype)*real(transy(2*j), kind=mytype)*real(transy(2*j),
           real(transz(k), kind=mytype)*real(transz(k), kind=mytype))&
           -(real(zk2(k), kind=mytype)*real(transy(2*j), kind=mytype)*real(transy(2*j), kind=mytype)
           real(transx(i), kind=mytype)*real(transx(i), kind=mytype))&
           -real(cw2(i,j,k), kind=mytype)*real(cw2(i,j,k), kind=mytype)*xa0*xa0-&
           xa1*xa1*(real(cw2(i,j,k), kind=mytype)*real(cw2(i,j-1,k), kind=mytype)+real(cw2(i,j,k), k
           real(cw2(i,j+1,k), kind=mytype)),&
           -(aimag(xk2(i))*aimag(transy(2*j))*aimag(transy(2*j))*aimag(transz(k))*aimag(transz(k)))&
           -(aimag(zk2(k))*aimag(transy(2*j))*aimag(transy(2*j))*aimag(transx(i))*aimag(transx(i)))&
           -aimag(cw2(i,j,k))*aimag(cw2(i,j,k))*xa0*xa0-&
           xa1*xa1*(aimag(cw2(i,j,k))*aimag(cw2(i,j-1,k))+aimag(cw2(i,j,k))*aimag(cw2(i,j+1,k))), ki
   enddo
   enddo
   do i=sp%yst(1),sp%yen(1)
      a(i,1,k,3)=cmplx(-(real(xk2(i), kind=mytype)*real(transy(1), kind=mytype)*real(transy(1), kind
           real(transz(k), kind=mytype)*real(transz(k), kind=mytype))&
           -(real(zk2(k), kind=mytype)*real(transy(1), kind=mytype)*real(transy(1), kind=mytype)*&
           real(transx(i), kind=mytype)*real(transx(i), kind=mytype))&
           -real(cw22(i,1,k), kind=mytype)*real(cw22(i,1,k), kind=mytype)*xa0*xa0-xa1*xa1*(real(cw22
           real(cw22(i,2,k), kind=mytype)),&
           -(aimag(xk2(i))*aimag(transy(1))*aimag(transy(1))*aimag(transz(k))*aimag(transz(k)))&
           -(aimag(zk2(k))*aimag(transy(1))*aimag(transy(1))*aimag(transx(i))*aimag(transx(i)))&
           -aimag(cw22(i,1,k))*aimag(cw22(i,1,k))*xa0*xa0-xa1*xa1*(aimag(cw22(i,1,k))*aimag(cw22(i,2
      a(i,ny/2,k,3)=cmplx(-(real(xk2(i), kind=mytype)*real(transy(ny-2), kind=mytype)*real(transy(ny
           *real(transz(k), kind=mytype)*real(transz(k), kind=mytype))&
           -(real(zk2(k), kind=mytype)*real(transy(ny-2), kind=mytype)*real(transy(ny-2), kind=mytyp
           real(transx(i), kind=mytype)*real(transx(i), kind=mytype))&
           -real(cw22(i,ny/2,k), kind=mytype)*real(cw22(i,ny/2,k), kind=mytype)*xa0*xa0-&
           xa1*xa1*(real(cw22(i,ny/2,k), kind=mytype)*real(cw22(i,ny/2-1,k), kind=mytype)),&
           -(aimag(xk2(i))*aimag(transy(ny-2))*aimag(transy(ny-2))*aimag(transz(k))*aimag(transz(k))
           -(aimag(zk2(k))*aimag(transy(ny-2))*aimag(transy(ny-2))*aimag(transx(i))*aimag(transx(i))
           -aimag(cw22(i,ny/2,k))*aimag(cw22(i,ny/2,k))*xa0*xa0-&
           xa1*xa1*(aimag(cw22(i,ny/2,k))*aimag(cw22(i,ny/2-1,k))), kind=mytype)
      a2(i,1,k,3)=cmplx(-(real(xk2(i), kind=mytype)*real(transy(2), kind=mytype)*real(transy(2), kin
           real(transz(k), kind=mytype)*real(transz(k), kind=mytype))&
           -(real(zk2(k), kind=mytype)*real(transy(2), kind=mytype)*real(transy(2), kind=mytype)*&
           real(transx(i), kind=mytype)*real(transx(i), kind=mytype))&
           -real(cw2(i,1,k), kind=mytype)*real(cw2(i,1,k), kind=mytype)*(xa0-xa1)*(xa0+xa1)-xa1*xa1*
           real(cw2(i,2,k), kind=mytype)),&
           -(aimag(xk2(i))*aimag(transy(2))*aimag(transy(2))*aimag(transz(k))*aimag(transz(k)))&
           -(aimag(zk2(k))*aimag(transy(2))*aimag(transy(2))*aimag(transx(i))*aimag(transx(i)))&
           -aimag(cw2(i,1,k))*aimag(cw2(i,1,k))*(xa0-xa1)*(xa0+xa1)-xa1*xa1*(aimag(cw2(i,1,k))*aimag
      a2(i,ny/2,k,3)=cmplx(-(real(xk2(i), kind=mytype)*real(transy(ny-1), kind=mytype)*real(transy(n
           real(transz(k), kind=mytype)*real(transz(k), kind=mytype))&
           -(real(zk2(k), kind=mytype)*real(transy(ny-1), kind=mytype)*real(transy(ny-1), kind=mytyp
           real(transx(i), kind=mytype)*real(transx(i), kind=mytype))&
           -real(cw2(i,ny/2,k), kind=mytype)*real(cw2(i,ny/2,k), kind=mytype)*(xa0+xa1)*(xa0+xa1)-&
           xa1*xa1*(real(cw2(i,ny/2,k), kind=mytype)*real(cw2(i,ny/2-1,k), kind=mytype)),&
           -(aimag(xk2(i))*aimag(transy(ny-1))*aimag(transy(ny-1))*aimag(transz(k))*aimag(transz(k))
           -(aimag(zk2(k))*aimag(transy(ny-1))*aimag(transy(ny-1))*aimag(transx(i))*aimag(transx(i))
           -aimag(cw2(i,ny/2,k))*aimag(cw2(i,ny/2,k))*(xa0+xa1)*(xa0+xa1)-&
           xa1*xa1*(aimag(cw2(i,ny/2,k))*aimag(cw2(i,ny/2-1,k))), kind=mytype)
   enddo
   enddo





!sup diag +1
   do k=sp%yst(3),sp%yen(3)
   do j=2,ny/2-1
   do i=sp%yst(1),sp%yen(1)
      a(i,j,k,4)=cmplx(xa0*xa1*(real(cw22(i,j,k), kind=mytype)*real(cw22(i,j+1,k), kind=mytype)+real
           real(cw22(i,j+1,k), kind=mytype)),&
           xa0*xa1*(aimag(cw22(i,j,k))*aimag(cw22(i,j+1,k))+aimag(cw22(i,j+1,k))*aimag(cw22(i,j+1,k)
      a2(i,j,k,4)=cmplx(xa0*xa1*(real(cw2(i,j,k), kind=mytype)*real(cw2(i,j+1,k), kind=mytype)+real(
           real(cw2(i,j+1,k), kind=mytype)),&
           xa0*xa1*(aimag(cw2(i,j,k))*aimag(cw2(i,j+1,k))+aimag(cw2(i,j+1,k))*aimag(cw2(i,j+1,k))),
   enddo
   enddo
   do i=sp%yst(1),sp%yen(1)
      a(i,1,k,4)=2.*cmplx((xa0*xa1*(real(cw22(i,1,k), kind=mytype)*real(cw22(i,2,k), kind=mytype)+re
           real(cw22(i,2,k), kind=mytype))),&
           (xa0*xa1*(aimag(cw22(i,1,k))*aimag(cw22(i,2,k))+aimag(cw22(i,2,k))*aimag(cw22(i,2,k)))),
      a2(i,1,k,4)=cmplx((xa0-xa1)*xa1*(real(cw2(i,1,k), kind=mytype)*real(cw2(i,2,k), kind=mytype))&
           +xa0*xa1*(real(cw2(i,2,k), kind=mytype)*real(cw2(i,2,k), kind=mytype)),&
           (xa0-xa1)*xa1*(aimag(cw2(i,1,k))*aimag(cw2(i,2,k)))&
           +xa0*xa1*(aimag(cw2(i,2,k))*aimag(cw2(i,2,k))), kind=mytype)
      a2(i,ny/2-1,k,4)=cmplx(xa0*xa1*(real(cw2(i,ny/2-1,k), kind=mytype)*real(cw2(i,ny/2,k), kind=my
           +(xa0+xa1)*xa1*(real(cw2(i,ny/2,k), kind=mytype)*real(cw2(i,ny/2,k), kind=mytype)),&
           xa0*xa1*(aimag(cw2(i,ny/2-1,k))*aimag(cw2(i,ny/2,k)))&
           +(xa0+xa1)*xa1*(aimag(cw2(i,ny/2,k))*aimag(cw2(i,ny/2,k))), kind=mytype)
      a2(i,ny/2,k,4)=0.
   enddo
   enddo



!sup diag +2
   do k=sp%yst(3),sp%yen(3)
   do i=sp%yst(1),sp%yen(1)
   do j=1,ny/2-2
      a(i,j,k,5)=cmplx(-real(cw22(i,j+1,k), kind=mytype)*real(cw22(i,j+2,k), kind=mytype)*xa1*xa1,&
           -aimag(cw22(i,j+1,k))*aimag(cw22(i,j+2,k))*xa1*xa1, kind=mytype)
      a2(i,j,k,5)=cmplx(-real(cw2(i,j+1,k), kind=mytype)*real(cw2(i,j+2,k), kind=mytype)*xa1*xa1,&
           -aimag(cw2(i,j+1,k))*aimag(cw2(i,j+2,k))*xa1*xa1, kind=mytype)
   enddo
   a(i,1,k,5)=cmplx(real(a(i,1,k,5), kind=mytype)*2.,aimag(a(i,1,k,5))*2., kind=mytype)
   a(i,ny/2-1,k,5)=0.
   a(i,ny/2,k,5)=0.
   a2(i,ny/2-1,k,5)=0.
   a2(i,ny/2,k,5)=0.
   enddo
   enddo



!inf diag -1
   do k=sp%yst(3),sp%yen(3)
   do i=sp%yst(1),sp%yen(1)
   do j=2,ny/2
      a(i,j,k,2)=cmplx(xa0*xa1*(real(cw22(i,j,k), kind=mytype)*real(cw22(i,j-1,k), kind=mytype)+real
           real(cw22(i,j-1,k), kind=mytype)),&
           xa0*xa1*(aimag(cw22(i,j,k))*aimag(cw22(i,j-1,k))+aimag(cw22(i,j-1,k))*aimag(cw22(i,j-1,k)
      a2(i,j,k,2)=cmplx(xa0*xa1*(real(cw2(i,j,k), kind=mytype)*real(cw2(i,j-1,k), kind=mytype)+real(
           real(cw2(i,j-1,k), kind=mytype)),&
           xa0*xa1*(aimag(cw2(i,j,k))*aimag(cw2(i,j-1,k))+aimag(cw2(i,j-1,k))*aimag(cw2(i,j-1,k))),
   enddo
   a(i,1,k,2)=0.
   a2(i,1,k,2)=0.
   a2(i,2,k,2)=cmplx(xa0*xa1*(real(cw2(i,2,k), kind=mytype)*real(cw2(i,1,k), kind=mytype))&
        +(xa0+xa1)*xa1*(real(cw2(i,1,k), kind=mytype)*real(cw2(i,1,k), kind=mytype)),&
        xa0*xa1*(aimag(cw2(i,2,k))*aimag(cw2(i,1,k)))&
        +(xa0+xa1)*xa1*(aimag(cw2(i,1,k))*aimag(cw2(i,1,k))), kind=mytype)
   a2(i,ny/2,k,2)=cmplx((xa0+xa1)*xa1*(real(cw2(i,ny/2,k), kind=mytype)*real(cw2(i,ny/2-1,k), kind=m
        +xa0*xa1*(real(cw2(i,ny/2-1,k), kind=mytype)*real(cw2(i,ny/2-1,k), kind=mytype)),&
        (xa0+xa1)*xa1*(aimag(cw2(i,ny/2,k))*aimag(cw2(i,ny/2-1,k)))&
        +xa0*xa1*(aimag(cw2(i,ny/2-1,k))*aimag(cw2(i,ny/2-1,k))), kind=mytype)
   enddo
   enddo
!inf diag -2
   do k=sp%yst(3),sp%yen(3)
   do i=sp%yst(1),sp%yen(1)
   do j=3,ny/2
      a(i,j,k,1)=cmplx(-real(cw22(i,j-1,k), kind=mytype)*real(cw22(i,j-2,k), kind=mytype)*xa1*xa1,&
           -aimag(cw22(i,j-1,k))*aimag(cw22(i,j-2,k))*xa1*xa1, kind=mytype)
      a2(i,j,k,1)=cmplx(-real(cw2(i,j-1,k), kind=mytype)*real(cw2(i,j-2,k), kind=mytype)*xa1*xa1,&
           -aimag(cw2(i,j-1,k))*aimag(cw2(i,j-2,k))*xa1*xa1, kind=mytype)
   enddo
   a(i,1,k,1)=0.
   a(i,2,k,1)=0.
   a2(i,1,k,1)=0.
   a2(i,2,k,1)=0.
   enddo
   enddo
!not to have a singular matrice
   do k=sp%yst(3),sp%yen(3)
   do i=sp%yst(1),sp%yen(1)
      if ((real(xk2(i), kind=mytype)==0.).and.(real(zk2(k), kind=mytype)==0)) then
         a(i,1,k,3)=cmplx(1.,1., kind=mytype)
         a(i,1,k,4)=0.
         a(i,1,k,5)=0.
      endif
   enddo
   enddo

else
   xa0=alpha/pi+1./2./beta/pi
   xa1=-1./4./beta/pi
!
!construction of the pentadiagonal matrice
!
   do k=sp%yst(3),sp%yen(3)
   do j=1,ny
   do i=sp%yst(1),sp%yen(1)
      cw22(i,j,k)=cmplx(real(yky(j), kind=mytype)*real(transx(i), kind=mytype)*real(transz(k), kind=
           aimag(yky(j))*aimag(transx(i))*aimag(transz(k)), kind=mytype)
   enddo
   enddo
   enddo

!main diagonal
   do k=sp%yst(3),sp%yen(3)
   do j=2,ny-1
   do i=sp%yst(1),sp%yen(1)
      a3(i,j,k,3)=cmplx(-(real(xk2(i), kind=mytype)*real(transy(j), kind=mytype)*real(transy(j), kin
           real(transz(k), kind=mytype)*real(transz(k), kind=mytype))&
           -(real(zk2(k), kind=mytype)*real(transy(j), kind=mytype)*real(transy(j), kind=mytype)*&
           real(transx(i), kind=mytype)*real(transx(i), kind=mytype))&
           -real(cw22(i,j,k), kind=mytype)*real(cw22(i,j,k), kind=mytype)*xa0*xa0-&
           xa1*xa1*(real(cw22(i,j,k), kind=mytype)*real(cw22(i,j-1,k), kind=mytype)+real(cw22(i,j,k)
           real(cw22(i,j+1,k), kind=mytype)),&
           -(aimag(xk2(i))*aimag(transy(j))*aimag(transy(j))*aimag(transz(k))*aimag(transz(k)))&
           -(aimag(zk2(k))*aimag(transy(j))*aimag(transy(j))*aimag(transx(i))*aimag(transx(i)))&
           -aimag(cw22(i,j,k))*aimag(cw22(i,j,k))*xa0*xa0-&
           xa1*xa1*(aimag(cw22(i,j,k))*aimag(cw22(i,j-1,k))+aimag(cw22(i,j,k))*aimag(cw22(i,j+1,k)))
   enddo
   enddo
   enddo
   do k=sp%yst(3),sp%yen(3)
   do i=sp%yst(1),sp%yen(1)
      a3(i,1,k,3)=cmplx(-(real(xk2(i), kind=mytype)*real(transy(1), kind=mytype)*real(transy(1), kin
           real(transz(k), kind=mytype)*real(transz(k), kind=mytype))&
           -(real(zk2(k), kind=mytype)*real(transy(1), kind=mytype)*real(transy(1), kind=mytype)*rea
           real(transx(i), kind=mytype))&
           -real(cw22(i,1,k), kind=mytype)*real(cw22(i,1,k), kind=mytype)*xa0*xa0-xa1*xa1*(real(cw22
           real(cw22(i,2,k), kind=mytype)),&
           -(aimag(xk2(i))*aimag(transy(1))*aimag(transy(1))*aimag(transz(k))*aimag(transz(k)))&
           -(aimag(zk2(k))*aimag(transy(1))*aimag(transy(1))*aimag(transx(i))*aimag(transx(i)))&
           -aimag(cw22(i,1,k))*aimag(cw22(i,1,k))*xa0*xa0-xa1*xa1*(aimag(cw22(i,1,k))*aimag(cw22(i,2
      a3(i,ny,k,3)=cmplx(-(real(xk2(i), kind=mytype)*real(transy(ny), kind=mytype)*real(transy(ny),
           real(transz(k), kind=mytype)*real(transz(k), kind=mytype))&
           -(real(zk2(k), kind=mytype)*real(transy(ny), kind=mytype)*real(transy(ny), kind=mytype)*r
           real(transx(i), kind=mytype))&
           -real(cw22(i,ny,k), kind=mytype)*real(cw22(i,ny,k), kind=mytype)*xa0*xa0-&
           xa1*xa1*(real(cw22(i,ny,k), kind=mytype)*real(cw22(i,ny-1,k), kind=mytype)),&
           -(aimag(xk2(i))*aimag(transy(ny))*aimag(transy(ny))*aimag(transz(k))*aimag(transz(k)))&
           -(aimag(zk2(k))*aimag(transy(ny))*aimag(transy(ny))*aimag(transx(i))*aimag(transx(i)))&
           -aimag(cw22(i,ny,k))*aimag(cw22(i,ny,k))*xa0*xa0-&
           xa1*xa1*(aimag(cw22(i,ny,k))*aimag(cw22(i,ny-1,k))), kind=mytype)
   enddo
   enddo
!sup diag +1
   do k=sp%yst(3),sp%yen(3)
   do i=sp%yst(1),sp%yen(1)
   do j=2,ny-1
      a3(i,j,k,4)=cmplx(xa0*xa1*(real(cw22(i,j,k), kind=mytype)*real(cw22(i,j+1,k), kind=mytype)+rea
           real(cw22(i,j+1,k), kind=mytype)),&
           xa0*xa1*(aimag(cw22(i,j,k))*aimag(cw22(i,j+1,k))+aimag(cw22(i,j+1,k))*aimag(cw22(i,j+1,k)
   enddo
   a3(i,1,k,4)=cmplx((xa0*xa1*(real(cw22(i,1,k), kind=mytype)*real(cw22(i,2,k), kind=mytype)+real(cw
        real(cw22(i,2,k), kind=mytype))),&
        (xa0*xa1*(aimag(cw22(i,1,k))*aimag(cw22(i,2,k))+aimag(cw22(i,2,k))*aimag(cw22(i,2,k)))), kin
   enddo
   enddo
!sup diag +2
   do k=sp%yst(3),sp%yen(3)
   do i=sp%yst(1),sp%yen(1)
   do j=1,ny-2
      a3(i,j,k,5)=cmplx(-real(cw22(i,j+1,k), kind=mytype)*real(cw22(i,j+2,k), kind=mytype)*xa1*xa1,&
           -aimag(cw22(i,j+1,k))*aimag(cw22(i,j+2,k))*xa1*xa1, kind=mytype)
   enddo
   !a3(i,1,k,5)=a3(i,1,k,5)*2.
   !a3(i,1,k,5)=0.
   a3(i,ny-1,k,5)=0.
   a3(i,ny,k,5)=0.
   enddo
   enddo
!inf diag -1
   do k=sp%yst(3),sp%yen(3)
   do i=sp%yst(1),sp%yen(1)
   do j=2,ny
      a3(i,j,k,2)=cmplx(xa0*xa1*(real(cw22(i,j,k), kind=mytype)*real(cw22(i,j-1,k), kind=mytype)+rea
           real(cw22(i,j-1,k), kind=mytype)),&
           xa0*xa1*(aimag(cw22(i,j,k))*aimag(cw22(i,j-1,k))+aimag(cw22(i,j-1,k))*aimag(cw22(i,j-1,k)
   enddo
   a3(i,1,k,2)=0.
   enddo
   enddo
!inf diag -2
   do k=sp%yst(3),sp%yen(3)
   do i=sp%yst(1),sp%yen(1)
   do j=3,ny
      a3(i,j,k,1)=cmplx(-real(cw22(i,j-1,k), kind=mytype)*real(cw22(i,j-2,k), kind=mytype)*xa1*xa1,&
           -aimag(cw22(i,j-1,k))*aimag(cw22(i,j-2,k))*xa1*xa1, kind=mytype)
   enddo
   a3(i,1,k,1)=0.
   a3(i,2,k,1)=0.
   enddo
   enddo
!not to have a singular matrice
!   do k=sp%yst(3),sp%yen(3)
!   do i=sp%yst(1),sp%yen(1)
!      if ((xkx(i)==0.).and.(zkz(k)==0)) then
         a3(1,1,1,3)=cmplx(1.,1., kind=mytype)
         a3(1,1,1,4)=0.
         a3(1,1,1,5)=0.
!      endif
!   enddo
!   enddo
endif

return
end subroutine matrice_refinement
