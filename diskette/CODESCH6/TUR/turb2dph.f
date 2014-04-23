c
c   ********************* subr fftqua
c  this subroutine perform the calculation of trigz for temperton fft
c
      subroutine fftqua
      include 'param.f'
      common/dim/n1,n1m,n2,n2m
      common/fftcom/ifxz(13),trigz(3*m1m/2+1)
      call fftfax(n1m,ifxz,trigz)
      return
      end
c
c   ********************* subr dsolv
c
c   this subroutine performs the solution of  tridiagonal
c   matrix efficiently vectorized
c   the trid. matrix is periodic
c   this routine can be eliminated, in fact here only tripvv is called
c
      subroutine dsolv(n2m,qk,n1mh)
      include 'param.f'
      dimension qi(m1mh,m2),si(m1mh,m2),qei(m1mh,m2)
      dimension qk(m1mh,m2)
      call tripvv(qk,qi,si,qei,1,n2m,1,n1mh)
      return
      end
c
c   ********************* subr phcalc
c  this subroutine perform the calculation of stream function
c  using fft in the x1 direction and periodic tridiagonal
c  solver in x2
c
      subroutine phcal(qcap,psi)
      include 'param.f'
      common/dim/n1,n1m,n2,n2m
      dimension qkr(m1mh,m2),qki(m1mh,m2)
     1         ,qcap(m1,m2),psi(m1,m2)
      common/fftcom/ifxz(13),trigz(3*m1m/2+1)
      real xr(m1md,m2),work(m1,m2)
      n1mh=n1m/2+1
      n1md=n1m+2
      n1mdu=n1md+1
      n1mu=n1m-1
      do 1 j=1,n2m
      xr(1,j)=-qcap(n1m,j)
c
c  real fft-1 applied to the vor(i,j)
c  how to store the quantities is explained in the ncar temperton
c  routine
c
      do 2 k=1,n1m
      ks=k+1
      xr(ks,j)=-qcap(k,j)
    2 continue
      xr(n1md,j)=-qcap(1,j)
    1 continue
      call fft99(xr,work,trigz,ifxz,1,m1md,n1m,n2m,-1)
      do 4 kk=1,n1mh
      kd=2*kk-1
      kp=2*kk
c
c   storage as in temperton  qkr real qki immag.
c
      do 3 j=1,n2m
      qkr(kk,j)=xr(kd,j)
      qki(kk,j)=xr(kp,j)
    3 continue
    4 continue
c
c   solution of poisson equation real part
c
      call dsolv(n2m,qkr,n1mh)
c
c   solution of poisson equation immag. part
c
      call dsolv(n2m,qki,n1mh)
c
c   store the psi in the wave number space as in temperton
c
      do 11 kk=1,n1mh
      kd=2*kk-1
      kp=2*kk
      do 13 j=1,n2m
      xr(kd,j)=qkr(kk,j)
      xr(kp,j)=qki(kk,j)
   13 continue
   11 continue
c
c    real fft +1 to have dph in phisycal space
c
      call fft99(xr,work,trigz,ifxz,1,m1md,n1m,n2m,1)
c
c    the dph output fro fft99 as in temperton
c
   22 continue
      do 41 j=1,n2m
      psi(n1m,j)=xr(1,j)
      do 42 k=1,n1mu
      ks=k+1
      psi(k,j)=xr(ks,j)
   42 continue
   41 continue
      return
      end
c     ********************** phini ******************************
c
      subroutine phini
      include 'param.f'
      common/cophes/ap(m2),ac(m1mh,m2),am(m2)
      common/dim/n1,n1m,n2,n2m
      common/mesh/dx1,dx1q,dx2,dx2q
      common/phco/iph
      common/inph/n1mh
      iph=0
      pig=2.*asin(1.)
      n1mh=n1m/2+1
      af=2*pig/n1m
      do 588 k=1,n1mh
c
c     definition of modified wave numer
c
      akap=2.*(-1+cos((k-1)*af))*dx1q
          acc=(akap-2.*dx2q)
c
c  coefficients for stream function calculation
c
          do 588 j=1,n2m
          ap(j)=dx2q
          ac(k,j)=acc
          am(j)=dx2q
 588      continue
      call fftqua
      return
      end
c
c  ****************************** subrout tripvv **********************
c
      subroutine tripvv( f,q,s,qe,j1,j2,k1,k2)
c
      include 'param.f'
      common/cophes/c(m2),b(m1mh,m2),a(m2)
      dimension f(m1mh,1),q(m1mh,1)
     1         ,s(m1mh,1),qe(m1mh,1),fn(m1mh),p(m1mh)
      ja = j1 + 1
      jj = j1 + j2
      do 20 k=k1,k2
      q(k,j1) = -c(j1)/b(k,j1)
      s(k,j1) = - a(j1)/b(k,j1)
      fn(k) = f(k,j2)
      f(k,j1) = f(k,j1)/b(k,j1)
   20 continue
c
c     forward elimination sweep
c
      do 10 j=ja,j2
      do 21 k=k1,k2
      p(k) =1./( b(k,j) + a(j)*q(k,j-1))
      q(k,j) = - c(j)*p(k)
      s(k,j) = - a(j)*s(k,j-1)*p(k)
      f(k,j) = ( f(k,j) - a(j)*f(k,j-1))*p(k)
   21 continue
   10 continue
c
c     backward pass
c
      do 22 k=k1,k2
      s(k,j2) = 1.
      qe(k,j2) = 0.
   22 continue
      do 11 i=ja,j2
      j = jj - i
      do 23 k=k1,k2
      s(k,j) = s(k,j) + q(k,j)*s(k,j+1)
      qe(k,j) = f(k,j) + q(k,j)*qe(k,j+1)
   23 continue
   11 continue
      do 24 k=k1,k2
      f(k,j2)=(fn(k) - c(j2)*qe(k,j1) - a(j2)*qe(k,j2-1))
     &       /(c(j2)*s(k,j1) + a(j2)*s(k,j2-1)  +b(k,j2))
   24 continue
c
c     backward elimination pass
c
      do 12 i=ja,j2
      j = jj -i
      do 25 k=k1,k2
      f(k,j) = f(k,j2)*s(k,j) + qe(k,j)
   25 continue
   12 continue
      return
      end
c
c   ********************* subr phwa
c  this subroutine perform  transforms from
c  physical to wave number in two dimensions (ifft=1)
c  or from wave number to physical (ifft=-1)
c  the input or output in wave numbers is a
c  complex and  transposed matrix
c
      subroutine phwam(qcap,y)
      include 'param.f'
      common/dim/n1,n1m,n2,n2m
      common/fftrc/ifr(13),trigr(3*m2m/2+1)
      common/fftcc/ifc(13),trigc(2*m1m)
      dimension  qcap(m1,m2)
      real xr(m2m+2,m1m),workr(m2m+1,m1m)
      complex y(m1m,m2mh),workc(m1m,m2mh)
      common/n2mhi/n2mh
      ifft=-1
c  perform the calculation of one quantity from
c  physical to wave number in two dimensions
c  storage as in temperton for fft99
      do 10 i=1,n1m
      xr(1,i)=qcap(i,n2m)
      xr(n2m+2,i)=qcap(i,1)
      do 10 j=1,n2m
      xr(j+1,i)=qcap(i,j)
  10  continue
c
c  real  fft applied to transform the first direction
c
      call fft99(xr,workr,trigr,ifr,1,m2m+2,n2m,n1m,ifft)
      do 20 j=1,n2mh
      jp=2*j
      jd=2*j-1
      do 20 i=1,n1m
      y(i,j)=cmplx(xr(jd,i),xr(jp,i))
   20 continue
c
c complex  fft applied to transform the second direction
c
       call cfft99(y,workc,trigc,ifc,1,m1m,n1m,n2mh,ifft)
      sca=1.0/float(n1m)
      do 21 j=1,n2mh
      do 21 i=1,n1m
      y(i,j)=y(i,j)*sca
   21 continue
      return
      end
c
c   ********************* subr phwa
c  this subroutine perform  transforms from
c  physical to wave number in two dimensions (ifft=1)
c  or from wave number to physical (ifft=-1)
c  the input or output in wave numbers is a
c  complex and  transposed matrix
c
      subroutine phwap(qcap,y)
      include 'param.f'
      common/dim/n1,n1m,n2,n2m
      common/fftrc/ifr(13),trigr(3*m2m/2+1)
      common/fftcc/ifc(13),trigc(2*m1m)
      dimension  qcap(m1,m2)
      real xr(m2m+2,m1m),workr(m2m+1,m1m)
      complex y(m1m,m2mh),workc(m1m,m2mh)
      common/n2mhi/n2mh
      ifft=+1
c  perform the calculation of one quantity from
c  wave number to physical in two dimensions
c
c
c complex  fft applied to transform the second direction
c
      call cfft99(y,workc,trigc,ifc,1,m1m,n1m,n2mh,ifft)
      do 32 j=1,n2mh
      jp=2*j
      jd=2*j-1
      do 30 i=1,n1m
      xr(jd,i)=real(y(i,j))
      xr(jp,i)=aimag(y(i,j))
   30 continue
   32 continue
      call fft99(xr,workr,trigr,ifr,1,m2m+2,n2m,n1m,ifft)
      do 40 i=1,n1m
      qcap(i,n2m)=xr(1,i)
      do 40 j=1,n2m-1
      qcap(i,j)=xr(j+1,i)
  40  continue
      return
      end
