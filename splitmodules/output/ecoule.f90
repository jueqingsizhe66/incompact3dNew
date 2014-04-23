subroutine ecoule(ux1,uy1,uz1)
!
!**********************************************************************

USE param
USE IBM
USE variables
USE decomp_2d

implicit none

integer  :: i,j,k,jj1,jj2
real, dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1
real :: x,y,z,ym
real :: r1,r2,r3,r
real :: uh,ud,um,xv,bruit1

bxx1=0.;bxy1=0.;bxz1=0.
byx1=0.;byy1=0.;byz1=0.
bzx1=0.;bzy1=0.;bzz1=0.

!ITYPE=1 --> Constant flow field
!ITYPE=2 --> Channel flow
!ITYPE=3 --> Wake flow
!ITYPE=4 --> Mixing layer with splitter plate
!ITYPE=5 --> Channel flow
!ITYPE=6 --> Taylor Green vortices
!ITYPE=7 --> Cavity flow
!ITYPE=8 --> Flat plate Boundary layer
!ITYPE=9 --> Tank

if (itype.eq.1) then
   um=0.5*(u1+u2)
   do k=1,xsize(3)
   do j=1,xsize(2)
      bxx1(j,k)=um
      bxy1(j,k)=0.
      bxz1(j,k)=0.
   enddo
   enddo
endif

if (itype.eq.2) then
   do k=1,xsize(3)
   do j=1,xsize(2)
      if (istret.eq.0) y=(j+xstart(2)-1-1)*dy-yly/2.
      if (istret.ne.0) y=yp(j+xstart(2)-1)-yly/2.
!      print *,nrank,j+xstart(2)-1,yp(j+xstart(2)-1),1.-y*y
      do i=1,xsize(1)
         ux1(i,j,k)=ux1(i,j,k)+1.-y*y
      enddo
   enddo

   enddo
endif

if (itype.eq.3) then

endif

if (itype.eq.4) then

endif

if (itype.eq.5) then
   if (nclx.ne.0) then
      print *,'NOT POSSIBLE'
      stop
   endif
   if (nclz.ne.0) then
      print *,'NOT POSSIBLE'
      stop
   endif
   if (ncly==0) then
      print *,'NOT POSSIBLE'
      stop
   endif
   do k=1,xsize(3)
   do j=1,xsize(2)
      bxx1(j,k)=1.
      bxy1(j,k)=0.
      bxz1(j,k)=0.
   enddo
   enddo
!   do k=1,xsize(3)
!   do j=1,xsize(2)
!      if (istret.eq.0) y=(j+xstart(2)-1-1)*dy-yly/2.
!      if (istret.ne.0) y=yp(j)-yly/2.
!      do i=1,xsize(1)
!         ux1(i,j,k)=ux1(i,j,k)+1.-y*y
!      enddo
!   enddo
!   enddo
endif

if (itype.eq.6) then
   t=0.
   xv=1./100.
   xxk1=twopi/xlx
   xxk2=twopi/yly
   do k=1,xsize(3)
      z=(k+xstart(3)-1-1)*dz
   do j=1,xsize(2)
      y=(j+xstart(2)-1-1)*dy
      do i=1,xsize(1)
         x=(i-1)*dx
         ux1(i,j,k)=sin(2.*pi*x)*cos(2.*pi*y)*cos(2.*pi*z)
!sin(xxk1*x)*cos(xxk2*y)*exp(-(xxk1*xxk1+xxk2*xxk2)*xv*t)
         uy1(i,j,k)=sin(2.*pi*y)*cos(2.*pi*x)*cos(2.*pi*z)
!-xxk1/xxk2*sin(xxk2*y)*cos(xxk1*x)*exp(-(xxk1*xxk1+xxk2*xxk2)*xv*t)
         uz1(i,j,k)=sin(2.*pi*z)*cos(2.*pi*x)*cos(2.*pi*y)
!0.
         bxx1(j,k)=0.
         bxy1(j,k)=0.
         bxz1(j,k)=0.
      enddo
   enddo
   enddo
endif

if (itype.eq.7) then

endif

if (itype.eq.8) then

endif

if (itype.eq.9) then

endif


if (itype.eq.10) then
   do k=1,xsize(3)
   do j=1,xsize(2)
      bxx1(j,k)=0.
      bxy1(j,k)=0.
      bxz1(j,k)=0.
   enddo
   enddo
endif

return
end subroutine ecoule
