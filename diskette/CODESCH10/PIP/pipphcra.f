c
c   ********************* subr fftqua
c  this subroutine perform the calculation of trigz for CRAY FFTs    
c
      subroutine fftqua
      include 'param.f' 
      common/fftcm1/ifx1(19),trigx1(2*m1m)
      common/fftcm33/ifx33(19),trigx33(2*m3m)
      common/fftcm3/ifx3(19),trigx3(2*m3m)
      common/fftcm11/ifx11(19),trigx11(2*m1m)
      common/idimsp/idsp
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
      call fftfax(n3m,ifx3,trigx3)
      call cftfax(n1m,ifx1,trigx1)
      call fftfax(n1m,ifx11,trigx11)
      call cftfax(n3m,ifx33,trigx33)
c
c   modified wave number
c
      do 26 k=1,n3m
      ak3(k)=2.*(1.-cos(an(k)/n3m))*dx3q
   26 continue
      do 28 i=1,n1m
      ak1(i)=2.*(1.-cos(ap(i)/n1m))*dx1q
   28 continue
      k1max=n1mh
      k3max=n3mh
      do i=1,n1mh
      apik2(i)=ap(i)/(2.*pi)
      write(16,*)i,apik2(i)
      enddo
      do k=1,n3mh
      ankk2(k)=sqrt((an(k)/n3m)**2*dx3q)
      write(16,*)k,ankk2(k)
      enddo
      return
      end
c
c   ********************* subr phcalc
c  this subroutine perform the calculation of dph , periodic direction
c  along x3 and x1 the real fourier transform along x3
c
      subroutine phcalc
      include 'param.f' 
      parameter (m3mh=m3m/2+1)
      dimension    xr(m3+2,m1m),work(2*m3,m1m)
      dimension    xrr(m1,m3mh),xri(m1,m3mh),wor(4*m1m,m3mh)
      common/fftcm1/ifx1(19),trigx1(2*m1m)
      common/fftcm3/ifx3(19),trigx3(2*m3m)
      common/phqua/xr,work,xrr,xri,wor
      n3mh=n3m/2+1
      n3md=n3m+2
      do j=1,n2m
      do i=1,n1m
      do k=1,n3m
      xr(k,i)=qcap(i,j,k)
      enddo
      enddo
c
c   2-d   real fft applied to the divg(qhat) 
c   from physical to wave number space
c
      call rfftmlt(xr,work,trigx3,ifx3,1,m3+2,n3m,n1m,-1)
      do  k=1,n3mh
      kp=2*k
      kd=2*k-1
      do  i=1,n1m
      xrr(i,k)=xr(kd,i)
      xri(i,k)=xr(kp,i)
      enddo
      enddo
c
c   2-d   complex fft 
c
      call cfftmlt(xrr,xri,wor,trigx1,ifx1,1,m1,n1m,n3mh,-1)
      do  k=1,n3mh
      do  i=1,n1m
      qcap(i,j,k)=xrr(i,k)/(n1m)
      dq(i,j,k)=xri(i,k)/(n1m)
      enddo
      enddo
      enddo
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
      do  j=1,n2m
      do  k=1,n3mh
      do  i=1,n1m
      xrr(i,k)=qcap(i,j,k)
      xri(i,k)=dq(i,j,k)
      enddo
      enddo
c
c   2-d complex fft 
c   from wave number space to physical space
c
      call cfftmlt(xrr,xri,wor,trigx1,ifx1,1,m1,n1m,n3mh,+1)
      do  k=1,n3mh
      kp=2*k
      kd=2*k-1
      do  i=1,n1m
      xr(kd,i)=xrr(i,k)
      xr(kp,i)=xri(i,k)
      enddo
      enddo
c
c  2-d real fft to get dph in physical space
c
      call rfftmlt(xr,work,trigx3,ifx3,1,m3+2,n3m,n1m,+1)
      do  k=1,n3m
      do  i=1,n1m
      dph(i,j,k)=xr(k,i)
      enddo
      enddo
      enddo
      return
      end
c
c   ********************* subr dsolv
c
c   this subroutine performs the solution of the poisson equation
c   by solving a tridigonal matrix at each wave number k1 and k3
c
      subroutine dsolv(qkk)
      include 'param.f' 
      dimension qkk(m1,m2,m3)
      do 11 k=1,n3m/2+1
      do 15 j=1,n2m
      do 15 i=1,n1m
      fphj(i,j)=qkk(i,j,k)
      acphj(i,j)=acph(j)-ak1(i)/rm(j)**2-ak3(k)
      apphj(j)=apph(j)
      amphj(j)=amph(j)
   15 continue
      if(k.eq.1) then
      fphj(1,1)=0.
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
c  ****************************** subrout phini  **********************
c
c   in this subr the coefficients of the poisson eq. for dph
c   are calculated this subr. is called only at the beginning
c
      subroutine phini
      include 'param.f' 
      call fftqua
c
c   tridiagonal matrix coefficients due to the discretization
c   along the radial direction
c
      do 1 jc=1,n2m
c
c
      jm=jmv(jc)
      jp=jpv(jc)
      a22icc=rc(jc)*jmc(jc)*dx2q/g2rc(jc)
      a22icp=rc(jp)*jpc(jc)*dx2q/g2rc(jp)
      ac2=-(a22icc+a22icp)
      ugmmm=1./rm(jc)/g2rm(jc)
      amph(jc)=a22icc*ugmmm
      apph(jc)=a22icp*ugmmm
      acph(jc)=ac2*ugmmm
c     write(66,166)jc,amph(jc),acph(jc),apph(jc)
  166 format(3x,i4,2x,3e12.5)
    1 continue
c
      write(6,*)'  end of phini'
      return
      end
