c
c   ********************* subr fftqua
c  this subroutine perform the calculation of trigz for temperton fft
c
      subroutine fftqua
      include 'param.f'
      parameter (m3m=m3-1,m1m=m1-1,m12=2*m1m,m32=2*m3m)
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/fftcm1/ifx1(13),trigx1(3*m1m/2+1)
      common/fftcm3/ifx3(13),trigx3(m32)
      common/waves/an(m3),ap(m1),ak3(m3),ak1(m1)
      common/d13/alx1,alx3
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      pi=2.*asin(1.)
      n3mh=n3m/2
      n1mh=n1m/2
      n3mp=n3mh+1
      n1mp=n1mh+1
c
c     wave number definition
c
      if(n3m.gt.1) then
      do 16 k=1,n3mh
   16 an(k)=(k-1)*2.*pi
      do 17 k=n3mp,n3m
   17 an(k)=-(n3m-k+1)*2.*pi
      nx3fft=n3m
      call cftfax(nx3fft,ifx3,trigx3)
c
c   modified wave number x3 direct.
c
      do 26 k=1,n3m
      ak3(k)=2.*(1.-cos(an(k)/n3m))*dx3q
   26 continue
                   else 
      ak3(1)=0.
                   endif
      do 18 i=1,n1mh
   18 ap(i)=(i-1)*2.*pi
      do 19 i=n1mp,n1m
   19 ap(i)=-(n1m-i+1)*2.*pi
      nx1fft=n1m
      call fftfax(nx1fft,ifx1,trigx1)
c
c   modified wave number x1 direct.
c
      do 28 i=1,n1m
      ak1(i)=2.*(1.-cos(ap(i)/n1m))*dx1q
   28 continue
c     write(6,*) 'ak1 ',(ak1(i),i=1,n1m)
c     write(6,*) 'ak3 ',(ak3(k),k=1,n3m)
      return
      end
c
c   ********************* subr phcalc
c  this subroutine perform the calculation of dph , periodic direction
c  along x3 and x1 real fourier transform along x1
c
      subroutine phcalc(qcap,dph)
      include 'param.f'
      parameter (m3m=m3-1,m1m=m1-1,m12=2*m1m,m32=2*m3m)
      dimension qcap(m1,m2,m3)
      dimension dph(m1,m2,m3)
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/fftcm1/ifx1(13),trigx1(3*m1m/2+1)
      common/fftcm3/ifx3(13),trigx3(m32)
      common/waves/an(m3),ap(m1),ak3(m3),ak1(m1)
      complex xa(m3m,m1m),wor(m3m,m1m)
      real xr(m1m+2,m3m),work(m1m+1,m3m)
      common/rhsc/rhs(m1,m2,m3)
      common/fftar/xr,work,xa,wor
      n1mh=n1m/2+1
      n1md=n1m+2
      do 1 j=1,n2m
      do 11 k=1,n3m
      xr(1,k)=qcap(n1m,j,k)
      xr(n1md,k)=qcap(1,j,k)
   11 continue
      do 12 i=1,n1m
      is=i+1
      do 12 k=1,n3m
      xr(is,k)=qcap(i,j,k)
   12 continue
c
c   2-d real  fft applied to the divg(qhat) by fft99
c   from physical to wave number space
c
      call fft99(xr,work,trigx1,ifx1,1,m1+1,n1m,n3m,-1)
      if(n3m.gt.1) then
      do 33 i=1,n1mh
      ip=2*i
      id=2*i-1
      do 33 k=1,n3m
      xa(k,i)=cmplx(xr(id,k),xr(ip,k))
  33  continue
c
c  complex fft
c
      call cfft99(xa,wor,trigx3,ifx3,1,m3m,n3m,n1mh,-1)
      do 13 k=1,n3m
      do 13 i=1,n1mh
      qcap(i,j,k)=real(xa(k,i)/(n3m))
      rhs(i,j,k)=aimag(xa(k,i)/(n3m))
   13 continue
                   else 
      do 43 i=1,n1mh
      ip=2*i
      id=2*i-1
      do 43 k=1,n3m
      qcap(i,j,k)=xr(id,k)
      rhs(i,j,k)=xr(ip,k)
   43 continue
                   endif
    1 continue
      qcap(1,1,1)=0.
      rhs(1,1,1)=0.
c
c   solution of poisson equation real part
c
      call dsolv(qcap)
c
c   solution of poisson equation immag. part
c
      call dsolv(rhs)
c
c   phi in wavenumber space
c
      do 2 j=1,n2m
      if(n3m.gt.1) then
      do 21 k=1,n3m
      do 21 i=1,n1mh
      xa(k,i)=cmplx(qcap(i,j,k),rhs(i,j,k))
   21 continue
c
c   2-d fft applied to the phi by cfft99
c   from wave number space to physical space
c
      call cfft99(xa,wor,trigx3,ifx3,1,m3m,n3m,n1mh,+1)
      do 34 i=1,n1mh
      ip=2*i
      id=2*i-1
      do 34 k=1,n3m
      xr(id,k)=real(xa(k,i))
      xr(ip,k)=aimag(xa(k,i))
  34  continue
                   else
      do 44 i=1,n1mh
      ip=2*i
      id=2*i-1
      do 44 k=1,n3m
      xr(id,k)=qcap(i,j,k)
      xr(ip,k)=rhs(i,j,k) 
  44  continue
                   endif
c
      call fft99(xr,work,trigx1,ifx1,1,m1+1,n1m,n3m,+1)
      do 22 i=1,n1m
      is=i+1
      do 22 k=1,n3m
      dph(i,j,k)=xr(is,k)
   22 continue
    2 continue
      return
      end
c
c   ********************* subr dsolv
c
c   this subroutine performs the solution of the poisson equation
c   by solving a tridigonal matrix at each wave number k1 and k3
      subroutine dsolv(qk)
      include 'param.f'
      parameter (m1m=m1-1,m2m=m2-1,m3m=m3-1)
      dimension qk(m1,m2,m3)
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/waves/an(m3),ap(m1),ak3(m3),ak1(m1)
      dimension amj(m2),acj(m1,m2),apj(m2),qsb(m1,m2),fj(m1,m2)
      common/ctrdph/amph(m2),acph(m2),apph(m2)
      n1mh=n1m/2+1
      do 11 k=1,n3m
      do 15 i=1,n1mh 
      do 16 j=1,n2m
      fj(i,j)=qk(i,j,k)
      acj(i,j)=acph(j)-ak1(i)-ak3(k)
      apj(j)=apph(j)
      amj(j)=amph(j)
   16 continue
c
c   zero wave number
c
      if(i.eq.1.and.k.eq.1) then
      acj(1,1)=1.
      apj(1)=0.
      amj(1)=0.
      endif
   15 continue
c
c   tridiagonal inversion
c
      call tribj(amj,acj,apj,fj,n2m,qsb,n1mh)
c
c  solution
c
      do 14 j=1,n2m
      do 14 i=1,n1mh
      qk(i,j,k)=qsb(i,j)
   14 continue
   11 continue
      return
      end
c
c
c  ****************************** subrout tribj  **********************
c
      subroutine tribj(a,b,c,r,n,u,m)
      include 'param.f'
      dimension gam(m1,m2),a(m2),b(m1,m2),c(m2),r(m1,m2)
      dimension bet(m1),u(m1,m2)
      do 10 i=1,m
      bet(i)=b(i,1)
      u(i,1)=r(i,1)/bet(i)
   10 continue
      do 11 j=2,n
      do 21 i=1,m
      gam(i,j)=c(j-1)/bet(i)
      bet(i)=b(i,j)-a(j)*gam(i,j)
      u(i,j)=(r(i,j)-a(j)*u(i,j-1))/bet(i)
   21 continue
   11 continue
      do 12 j=n-1,1,-1
      do 22 i=1,m
      u(i,j)=u(i,j)-gam(i,j+1)*u(i,j+1)
   22 continue
   12 continue
      return
      end
c  ****************************** subrout phini  **********************
c
c   in this subr the coefficients of the poisson eq. for dph
c   are calculated this subr. is called only at the beginning
c
      subroutine phini(qcap,dph)
      include 'param.f'
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      common/tstep/dt,beta,ren
      common/indbo/imv(m1),ipv(m1),jmmv(m2),jppv(m2),kmv(m3),kpv(m3)
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/ctrdph/amph(m2),acph(m2),apph(m2)
      common/metria/caj(m2),cac(m2)
      dimension qcap(m1,m2,m3)
      dimension dph(m1,m2,m3)
      call fftqua
c
c   tridiagonal matrix coefficients at each k and i
c   x1 and x3 cartesian coordinates
c
      do 1 jc=1,n2m
      jm=jmmv(jc)
      jp=jppv(jc)
      a22icc=float(jc-jm)/cac(jc)
      a22icp=float(jp-jc)/cac(jp)
      ac2=-(a22icc+a22icp)
      anc=a22icp
      asc=a22icc
      ugmmm=dx2q/caj(jc)
      amph(jc)=(asc)*ugmmm
      apph(jc)=(anc)*ugmmm
      acph(jc)=ac2*ugmmm
    1 continue
c
c  at the point i=1,j=1 a zero value for dph is assumed
c  this hold only for ak3=0. dph is calculated at least of a
c  constant due to the implicit newman conditions for dph eq.
c
      return
      end
