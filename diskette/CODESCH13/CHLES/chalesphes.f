c
c   ********************* subr fftqua
c  this subroutine perform the initialization for
c  the "pressure" solver here the ESSL libraries are
c  used for the FFT
c
      subroutine fftqua(qcap,dph)
      include 'param.f'
      parameter (m3m=m3-1,m1m=m1-1,m12=2*m1m,m32=2*m3m)
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/waves/an(m3),ap(m1),ak3(m3),ak1(m1)
      common/d13/alx1,alx3
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      common/ifftin/iftin
      dimension qcap(m1,m2,m3)
      dimension dph(m1,m2,m3)
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
c
c   modified wave number
c
      do 26 k=1,n3m
      ak3(k)=2.*(1.-cos(an(k)/n3m))*dx3q
   26 continue
      do 28 i=1,n1m
      ak1(i)=2.*(1.-cos(ap(i)/n1m))*dx1q
   28 continue
      iftin=1
      call phcalc(qcap,dph)
      iftin=0
      return
      end
c
c   ********************* subr phcalc
c  this subroutine perform the calculation of dph , periodic direction
c  along x3 and x1 the appropriate routine of ESSL is used
c
      subroutine phcalc(qcap,dph)
      include 'param.f'
      parameter (m3m=m3-1,m1m=m1-1,m12=2*m1m,m32=2*m3m)
      parameter (m3mh=m3m/2+1,m3p=m3+1)
      dimension qcap(m1,m2,m3)
      dimension dph(m1,m2,m3)
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/waves/an(m3),ap(m1),ak3(m3),ak1(m1)
      common/rhsc/rhs(m1,m2,m3)
      real*8 aux1(22000),aux2(10000),aux3(m3m*m1m/2)
      real*8 aaux1(22000),aaux2(10000),aaux3(m3m*m1m/2)
      real*4 xa(m3m,m1m),xb(m3m+2,m1m)
      complex*8 y(m3m/2+1,m1m)
      common/ifftin/iftin
      n3mh=n3m/2+1
      n3fft=n3m
      n1fft=n1m
      do j=1,n2m
      do i=1,n1m
      do k=1,n3m
      xa(k,i)=qcap(i,j,k)
      enddo
      enddo
c
c   2-d   fft applied to the divg(qhat) by essl IBM
c   from physical to wave number space
c
c     isign=+1 from physical to wave numbers
c     isign=-1 from wave to physical
c     scale=1. from phys to wav
c     scale=1./n1m from  wav to phys
c
c
      if(iftin.eq.1) then
      naux1=22000
      naux2=10000
      naux3=n3m*n1m/2
      call srcft2 (1,xa,m3m,y,m3mh,n3m,n1m,1,1.0,aux1,naux1,
     #aux2,naux2,aux3,naux3)
                      else
C
      call srcft2 (0,xa,m3m,y,m3mh,n3m,n1m,1,1.0,aux1,naux1,
     #aux2,naux2,aux3,naux3)
                      endif
      do  k=1,n3mh
      do  i=1,n1m
      qcap(i,j,k)=real(y(k,i))
      rhs(i,j,k)=aimag(y(k,i))
      enddo
      enddo
      enddo
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
 
      do  j=1,n2m
      do  i=1,n1m
      do  k=1,n3mh
      y(k,i)=cmplx(qcap(i,j,k),rhs(i,j,k))
      enddo
      enddo
c
c   2-d fft applied to the phi by ESSL IBM
c   from wave number space to physical space
c
       sca=1.0/float(n1m*n3m)
      if(iftin.eq.1) then
C
       call scrft2 (1,y,m3mh,xb,m3p,n3m,n1m,-1,sca,
     # aaux1,naux1,aaux2,naux2,aaux3,naux3)
                     else
C
       call scrft2 (0,y,m3mh,xb,m3p,n3m,n1m,-1,sca,
     # aaux1,naux1,aaux2,naux2,aaux3,naux3)
C
                     endif
      do  k=1,n3m
      do  i=1,n1m
      dph(i,j,k)=xb(k,i)
      enddo
      enddo
      enddo
      return
      end
c
c   ********************* subr dsolv
c
c   this subroutine performs the solution of the poisson equation
c   by solving a tridigonal matrix at each wave number n and p
      subroutine dsolv(qk)
      include 'param.f'
      parameter (m1m=m1-1,m2m=m2-1,m3m=m3-1)
      dimension qk(m1,m2,m3)
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/waves/an(m3),ap(m1),ak3(m3),ak1(m1)
      dimension amj(m2),acj(m1,m2),apj(m2),qsb(m1,m2),fj(m1,m2)
      common/ctrdph/amph(m2),acph(m2),apph(m2)
      n1mh=n1m/2+1
      do 11 k=1,n3m/2+1
      do 15 i=1,n1m
      do 16 j=1,n2m
      fj(i,j)=qk(i,j,k)
      acj(i,j)=acph(j)-ak1(i)-ak3(k)
      apj(j)=apph(j)
      amj(j)=amph(j)
   16 continue
      if(i.eq.1.and.k.eq.1) then
c     fj(1,1)=0.
      acj(1,1)=1.
      apj(1)=0.
      amj(1)=0.
      endif
   15 continue
c
c   tridiagonal inversion
c
      call tribj(amj,acj,apj,fj,n2m,qsb,n1m)
c
c  solution
c
      do 14 j=1,n2m
      do 14 i=1,n1m
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
      call fftqua(qcap,dph)
      return
      end
