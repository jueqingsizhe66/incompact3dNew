c
c   ********************* subr fftqua
c  this subroutine perform the calculation of trigz for temperton fft
c
      subroutine fftqup
      include 'param.f'
      common/dim/n1,n1m,n2,n2m
      common/fftcop/ifxz(13),trigz(3*m1m/2+1)
      call fftfax(n1m,ifxz,trigz)
      return
      end
c
c   ********************* subr dsolv
c
c   this subroutine performs the solution of  tridiagonal
c   matrix efficiently vectorized
c   the trid. matrix is periodic
c
      subroutine dsolv(n2m,qk,n1mh)
      include 'param.f'
      dimension qi(m1mh,m2),si(m1mh,m2),qei(m1mh,m2)
      dimension qk(m1mh,m2)
      call tripvv(qk,qi,si,qei,1,n2m,1,n1mh)
      return
      end
c
c   ********************* subr psftpe
c  this subroutine perform the calculation of stream function
c  using fft in the x1 direction and periodic tridiagonal 
c  solver in x2
c
      subroutine psftpe(qcap,psi)
      include 'param.f'
      common/dim/n1,n1m,n2,n2m
      dimension qkr(m1mh,m2),qki(m1mh,m2)
     1         ,qcap(m1,m2),psi(m1,m2)
      common/fftcop/ifxz(13),trigz(3*m1m/2+1)
      real xr(m1md,m2),work(m1,m2)
      n1mh=n1m/2+1
      n1md=n1m+2
      n1mdu=n1md+1
      n1mu=n1m-1
      q1ma=0.
      q2ma=0.
      q3ma=0.
      q4ma=0.
c
c  fft-1 applied to the vor(i,j) RHS of the Laplacian
c  how to store the quantities is explained in the 
c  NCAR Temperton  routine
c
      do 1 j=1,n2m
      xr(1,j)=qcap(n1m,j)
      do 2 k=1,n1m
      ks=k+1
      xr(ks,j)=qcap(k,j)
      q1ma=max(abs(qcap(k,j)),q1ma)
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
      q2ma=max(abs(xr(kd,j)),q2ma)
    3 continue
    4 continue
      qkr(1,1)=0.
      qki(1,1)=0.
c
c    here the tridiagonal matrix with periodic conditions
c    are inverted
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
      q3ma=max(abs(xr(kd,j)),q3ma)
   13 continue
   11 continue
c
c   fft +1 to have dph in physical space
c
      call fft99(xr,work,trigz,ifxz,1,m1md,n1m,n2m,1)
c
c    the dph output fro fft99 as in temperton
c
   22 continue
c
c   the solution is stored as explained in the NCARFFT routines
c
      do 41 j=1,n2m
      psi(n1m,j)=xr(1,j)
      do 42 k=1,n1mu
      ks=k+1
      psi(k,j)=xr(ks,j)
      q4ma=max(abs(psi(k,j)),q4ma)
   42 continue
   41 continue
      write(6,*)'ppeft   ',q1ma,q2ma,q3ma,q4ma
      return
      end
c     ********************** phini ******************************
c
      subroutine phinip
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
c  coefficients of the tridiasgonal matrices   
c
          do 588 j=1,n2m
          ap(j)=dx2q
          ac(k,j)=acc
          am(j)=dx2q
 588      continue
      call fftqup
      return
      end
c
c  ****************************** subrout tripvv **********************
c
      subroutine tripvv( f,q,s,qe,j1,j2,k1,k2)
      include 'param.f'
c
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
