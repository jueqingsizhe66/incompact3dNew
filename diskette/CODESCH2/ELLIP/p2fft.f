c
c   ********************* subr fftqua
c  this subroutine perform the calculation of trigz for temperton fft
c
      subroutine fftq2p
      include 'param.f' 
      common/dim/n1,n1m,n2,n2m
      common/mesh/dx1,dx1q,dx2,dx2q
      common/fftcm1/ifx1(13),trigx1(m1f)
      common/fftcm2/ifx2(13),trigx2(3*m2m/2+1)
      common/wav2f/an(m2),ap(m1),ak2(m2),ak1(m1)
      pi=2.*asin(1.)
      n2mh=n2m/2
      n1mh=n1m/2
      n2mp=n2mh+1
      n1mp=n1mh+1
c
c     wave number definition
c
      do 16 k=1,n2mh
   16 an(k)=(k-1)*2.*pi
      do 17 k=n2mp,n2m
   17 an(k)=-(n2m-k+1)*2.*pi
      do 18 i=1,n1mh
   18 ap(i)=(i-1)*2.*pi
      do 19 i=n1mp,n1m
   19 ap(i)=-(n1m-i+1)*2.*pi
      nx2fft=n2m
      nx1fft=n1m
      call fftfax(nx2fft,ifx2,trigx2)
      call cftfax(nx1fft,ifx1,trigx1)
c
c   modified wave number
c
      do 26 k=1,n2m
      ak2(k)=2.*(1.-cos(an(k)/n2m))*dx2q
   26 continue
      do 28 i=1,n1m
      ak1(i)=2.*(1.-cos(ap(i)/n1m))*dx1q
   28 continue
      return
      end
c
c   ********************* subr phcalc
c  this subroutine perform the calculation of dph , periodic direction
c  along x3 and x1to use the real fourier transform
c
      subroutine p2fft(qcap,psi)
      include 'param.f' 
      complex xa(m1m,m2m),wor(m1m,m2m)
      dimension xr(m2m+2,m1m),work(m2m+1,m1m)
      dimension qcap(m1,m2),psi(m1,m2),rhs(m1,m2)
      common/dim/n1,n1m,n2,n2m
      common/fftcm1/ifx1(13),trigx1(m1f)
      common/fftcm2/ifx2(13),trigx2(3*m2m/2+1)
      common/wav2f/an(m2),ap(m1),ak2(m2),ak1(m1)
      n2mh=n2m/2+1
      n2md=n2m+2
      q1ma=0.
      q2ma=0.
      q3ma=0.
      q4ma=0.
      q5ma=0.
c
c   the RHS of the Laplacian is stored as requested
c   by the NCARFFT routines
c
      do 11 i=1,n1m
      xr(1,i)=qcap(i,n2m)
      do 12 k=1,n2m
      ks=k+1
      xr(ks,i)=qcap(i,k)
      q1ma=max(abs(qcap(i,k)),q1ma)
   12 continue
      xr(n2md,i)=qcap(i,1)
   11 continue
c
c   2-d   fft applied to the RHS  to go 
c   from physical to wave number space
c   isign=-1
c     real FFT fft99
c
      call fft99(xr,work,trigx2,ifx2,1,m2+1,n2m,n1m,-1)
      do 33 k=1,n2mh
      kp=2*k
      kd=2*k-1
      do 32 i=1,n1m
      xa(i,k)=cmplx(xr(kd,i),xr(kp,i))
      q2ma=max(abs(xr(kd,i)),q2ma)
  32  continue
  33  continue
c
c     complex FFT cfft99
c
      call cfft99(xa,wor,trigx1,ifx1,1,m1m,n1m,n2mh,-1)
      do 13 i=2,n1m
      do 13 k=1,n2mh
      den=-1./(float(n1m)*(ak1(i)+ak2(k)))
      rhs(i,k)=real(xa(i,k)*den)
      psi(i,k)=aimag(xa(i,k)*den)
   13 continue
      rhs(1,1)=0.
      psi(1,1)=0.
      do 21 i=1,n1m
      do 21 k=1,n2mh
      q5ma=max(abs(rhs(i,k)),q5ma)
      xa(i,k)=cmplx(rhs(i,k),psi(i,k))
   21 continue
c
c   2-d fft applied to the psi 
c   from wave number space to physical space
c   isign=+1
c    complex FFT by cfft99
c
      call cfft99(xa,wor,trigx1,ifx1,1,m1m,n1m,n2mh,+1)
      do 34 k=1,n2mh
      kp=2*k
      kd=2*k-1
      do 35 i=1,n1m
      xr(kd,i)=real(xa(i,k))
      xr(kp,i)=aimag(xa(i,k))
      q3ma=max(abs(xr(kd,i)),q3ma)
  35  continue
  34  continue
c
c    real FFT by fft99
c
      call fft99(xr,work,trigx2,ifx2,1,m2+1,n2m,n1m,+1)
c
c   solution stored as suggested in the NCARFFT routines
c
      do k=1,n2m
      ks=k+1
      do i=1,n1
      psi(i,k)=xr(ks,i)
      q4ma=max(abs(psi(i,k)),q4ma)
      enddo
      enddo
      write(6,*)'p2fft   ',q1ma,q2ma,q5ma,q3ma,q4ma
      return
      end
c
c  ****************************** subrout phini  **********************
c
c   in this subr the coefficients of the Poisson eq. for psi
c   are calculated this subr. is called only at the beginning
c
      subroutine phin2p
      include 'param.f' 
      call fftq2p
      return
      end
