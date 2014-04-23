c
c   ********************* subr fftqua
c  this subroutine perform the calculation of trig* 
c   requested for temperton fft
c   in the 3D case real fft in x3 and then complex fft for x1
c   in the 2D case real fft in x1
c
      subroutine fftqua
      include 'param.f' 
      common/fftcm1/ifx1(13),trigx1(m12),nx3fft,nx1fft
      common/fftcm3/ifx3(13),trigx3(3*m3m/2+1)
      common/fftm11/ifx11(13),trigx11(3*m1m/2+1)
      pi=2.*asin(1.)
      n3mh=n3m/2
      n1mh=n1m/2
      n3mp=n3mh+1
      n1mp=n1mh+1
c
c     wave number definition
c
      do 16 k=1,n3mh
   16 an(k)=(k-1)*2.*pi
      do 17 k=n3mp,n3m
   17 an(k)=-(n3m-k+1)*2.*pi
      do 18 i=1,n1mh
   18 ap(i)=(i-1)*2.*pi
      do 19 i=n1mp,n1m
   19 ap(i)=-(n1m-i+1)*2.*pi
      nx3fft=n3m
      nx1fft=n1m
      if(n3m.ne.1) then
      call fftfax(nx3fft,ifx3,trigx3)
      if(n1m.ne.1) 
     %call cftfax(nx1fft,ifx1,trigx1)
                   endif
      if(n3m.eq.1) then
      call fftfax(nx1fft,ifx11,trigx11)
                   endif
c
c   modified wave number
c
      do 26 k=1,n3m
      ak3(k)=2.*(1.-cos(an(k)/n3m))*dx3q
   26 continue
      do 28 i=1,n1m
      ak1(i)=2.*(1.-cos(ap(i)/n1m))*dx1q
   28 continue
      return
      end
c
c   ********************* subr phcaij
c  this subroutine perform the calculation of dph in the r-theta planes
c        2D simulations
c
      subroutine phcaij
      include 'param.f'
      common/fftm11/ifx1(13),trigx1(3*m1m/2+1)
      dimension xr(m1m+2,m2m),work(m1m+1,m2m)
      n1mh=n1m/2+1
      n1md=n1m+2
      k=1
      do 11 j=1,n2m
      xr(1,j)=qcap(n1m,j,k)
      do 12 i=1,n1m
      is=i+1
      xr(is,j)=qcap(i,j,k)
   12 continue
      xr(n1md,j)=qcap(1,j,k)
   11 continue
c
c   2-d   fft applied to the divg(qhat) by cfft99
c   from physical to wave number space
c
      call fft99(xr,work,trigx1,ifx1,1,m1+1,n1m,n2m,-1)
      do 43 j=1,n2m
      do 43 i=1,n1mh
      ip=2*i
      id=2*i-1
      qcap(i,j,k)=xr(id,j)
      dq(i,j,k)=xr(ip,j)
   43 continue
c
c   solution of poisson equation real part
c
      call dsolv(qcap)
c
c   solution of poisson equation immag. part
c
      call dsolv(dq)
c
c   phi in wavenumber space
c
      do 44 j=1,n2m
      do 44 i=1,n1mh
      ip=2*i
      id=2*i-1
      xr(id,j)=qcap(i,j,k)
      xr(ip,j)=dq(i,j,k)
  44  continue
c
      call fft99(xr,work,trigx1,ifx1,1,m1+1,n1m,n2m,+1)
      do 22 i=1,n1m
      is=i+1
      do 22 j=1,n2m
   22 dph(i,j,k)=xr(is,j)
      return
      end
c
c   ********************* subr phcalc
c  this subroutine perform the calculation of dph , periodic direction
c  along x3 and x1 real fourier transform in x3
c      3D simulations
c
      subroutine phcalc
      include 'param.f' 
      common/fftcm1/ifx1(13),trigx1(m12),nx3fft,nx1fft
      dimension xr(m3m+2,m1m),work(m3m+1,m1m)
      common/fftcm3/ifx3(13),trigx3(3*m3m/2+1)
      complex xa(m1m,m3m),wor(m1m,m3m)
      n3mh=n3m/2+1
      n3md=n3m+2
      do 1 j=1,n2m
      do 11 i=1,n1m
      xr(1,i)=qcap(i,j,n3m)
      do 12 k=1,n3m
      ks=k+1
      xr(ks,i)=qcap(i,j,k)
   12 continue
      xr(n3md,i)=qcap(i,j,1)
   11 continue
c
c   2-d   fft applied to the divg(qhat) by cfft99
c   from physical to wave number space
c
      call fft99(xr,work,trigx3,ifx3,1,m3+1,n3m,n1m,-1)
      do 33 k=1,n3mh
      kp=2*k
      kd=2*k-1
      do 32 i=1,n1m
      xa(i,k)=cmplx(xr(kd,i),xr(kp,i))
  32  continue
  33  continue
c
c
      if(n1m.ne.1) 
     %call cfft99(xa,wor,trigx1,ifx1,1,m1m,n1m,n3mh,-1)
      do 13 i=1,n1m
      do 13 k=1,n3mh
      qcap(i,j,k)=real(xa(i,k)/(n1m))
      dq(i,j,k)=aimag(xa(i,k)/(n1m))
   13 continue
    1 continue
      qcap(1,1,1)=0.
      dq(1,1,1)=0.
c
c   solution of poisson equation real part
c
      call dsolv(qcap)
c
c   solution of poisson equation immag. part
c
      call dsolv(dq)
c
c   phi in wavenumber space
c
      do 2 j=1,n2m
      do 21 i=1,n1m
      do 21 k=1,n3mh
      xa(i,k)=cmplx(qcap(i,j,k),dq(i,j,k))
   21 continue
c
c   2-d fft applied to the phi by cfft99
c   from wave number space to physical space
c
      if(n1m.ne.1) 
     % call cfft99(xa,wor,trigx1,ifx1,1,m1m,n1m,n3mh,+1)
      do 34 k=1,n3mh
      kp=2*k
      kd=2*k-1
      do 35 i=1,n1m
      xr(kd,i)=real(xa(i,k))
      xr(kp,i)=aimag(xa(i,k))
  35  continue
  34  continue
c
      call fft99(xr,work,trigx3,ifx3,1,m3+1,n3m,n1m,+1)
      do 22 k=1,n3m
      ks=k+1
      do 22 i=1,n1m
   22 dph(i,j,k)=xr(ks,i)
    2 continue
      return
      end
c
c
c  ****************************** subrout tribj  **********************
c
      subroutine tribj(n,m)
      include 'param.f' 
      do 10 i=1,m
      bet(i)=acphj(i,1)
      qsbph(i,1)=fphj(i,1)/bet(i)
   10 continue
      do 11 j=2,n
      do 21 i=1,m
      gm(i,j)=apphj(j-1)/bet(i)
      bet(i)=acphj(i,j)-amphj(j)*gm(i,j)
      qsbph(i,j)=(fphj(i,j)-amphj(j)*qsbph(i,j-1))/bet(i)
   21 continue
   11 continue
      do 12 j=n-1,1,-1
      do 22 i=1,m
      qsbph(i,j)=qsbph(i,j)-gm(i,j+1)*qsbph(i,j+1)
   22 continue
   12 continue
      return
      end
c
c   ********************* subr dsolv
c
c   this subroutine performs the solution of the poisson equation
c   by solving a tridigonal matrix at each wave number k1 and k3
c   wave numbers
c
      subroutine dsolv(qkk)
      include 'param.f' 
      dimension qkk(m1,m2,m3)
      do 11 k=1,n3m/2+1
      do 15 j=1,n2m
      do 15 i=1,n1m
      fphj(i,j)=qkk(i,j,k)
      acphj(i,j)=acph(j)-ak1(i)-ak3(k)*rm(j)*rm(j)
      apphj(j)=apph(j)
      amphj(j)=amph(j)
   15 continue
      if(k.eq.1) then
      acphj(1,1)=1.
      apphj(1)=0.
      amphj(1)=0.
      endif
c
c   tridiagonal inversion
c
      call tribj(n2m,n1m)
c
c  solution
c
      do 14 j=1,n2m
      do 14 i=1,n1m
      qkk(i,j,k)=qsbph(i,j)
   14 continue
   11 continue
      return
      end
c
c  ****************************** subrout phini  **********************
c
c   in this subr the coefficients of the poisson eq. for dph
c   are calculated this subr. is called only at the beginning
c
      subroutine phini
      include 'param.f' 
      call fftqua
c
c   tridiagonal matrix coefficients from the discretization
c   along radial direction
c
      do 1 jc=1,n2m
c
c
      jm=jmv(jc)
      jp=jpv(jc)
      a22icc=rc(jc)*jmc(jc)*dx2q/g2rc(jc)
      a22icp=rc(jp)*jpc(jc)*dx2q/g2rc(jp)
      ac2=-(a22icc+a22icp)
      ugmmm=rm(jc)/g2rm(jc)
      amph(jc)=a22icc*ugmmm
      apph(jc)=a22icp*ugmmm
      acph(jc)=ac2*ugmmm
      write(66,166)jc,amph(jc),acph(jc),apph(jc)
  166 format(3x,i4,2x,3e12.5)
    1 continue
c
      return
      end
