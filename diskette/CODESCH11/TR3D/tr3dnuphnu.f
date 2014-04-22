c
c   ********************* subr fftqua
c  this subroutine perform the calculation of trigx1 for temperton fft
c   and evaluates reduced wave number
c
      subroutine fftqua
      include 'param.f' 
      common /fftcm1/ ifx1(13),trigx1(m12)
      pi=2.*asin(1.)
      n1mh=n1m/2
      n1mp=n1mh+1
      if(n1m.gt.1) then
      do 18 i=1,n1mh
   18 ap(i)=(i-1)*2.*pi
      do 19 i=n1mp,n1m
   19 ap(i)=-(n1m-i+1)*2.*pi
      nx1fft=n1m
      call fftfax(nx1fft,ifx1,trigx1)
c
c   modified wave number
c
      do 28 i=1,n1m
      ak1(i)=2.*(1.-cos(ap(i)/n1m))*dx1q
   28 continue
                   else
      ak1(1)=0.
      ap(1)=0.
                   endif
      return
      end
c
c   ********************* subr phcalc
c  this subroutine perform the calculation of dph , periodic direction
c   x1 azimuthal direct. to use the real fourier transform
c
      subroutine phcalc
      include 'param.f' 
      common /fftcm1/ ifx1(13),trigx1(m12)
      real xr(m1+1,m2m),work(m1,m2m)
      dimension acphjj(m2)

      my=m2
      pi=2.*asin(1.)
      n1mh=n1m/2+1
      if(n1m.gt.1) then
c
c   fft applied to the x1 direction to the divg of q^
c   from physical to wave number space
c    Temperton storage
c
      do k=1,n3m
      do j=1,n2m
      xr(1,j)=qcap(n1m,j,k) 
      do i=1,n1m
      is=i+1
      xr(is,j)=qcap(i,j,k)
      enddo
      enddo
c
c   2-d   fft applied to the divg(qhat) by fft99
c   from physical to wave number space
c
      call fft99(xr,work,trigx1,ifx1,1,m1+1,n1m,n2m,-1)
c
c   real and imaginary parts
c
      do j=1,n2m
      do i=1,n1mh
      ip=2*i
      id=2*i-1
      qcap(i,j,k)=xr(id,j)
      dq(i,j,k)=xr(ip,j)
      enddo
      enddo
      enddo
                 endif
c
c   end  FFT of divg qhat from phys to wave number 
c
c
c    solution of real part by fishpack
c
            do i=1,n1mh
      do j=1,n2m
      acphjj(j)=acphj(i,j)
      enddo
      do j=1,n2m
      do k=1,n3m
      yfis(j,k)=qcap(i,j,k)
      enddo
      enddo
      call blktri(1,np,n3m,amphk,acphk,apphk
     1             ,mp,n2m,amphj,acphjj,apphj,my,yfis,ierror,w)
      do j=1,n2m
      do k=1,n3m
      qcap(i,j,k)=yfis(j,k)
      enddo
      enddo
      if(n1m.gt.1) then
c
c    solution of imaginary part by fishpack
c
      do j=1,n2m
      do k=1,n3m
      yfis(j,k)=dq(i,j,k)
      enddo
      enddo
      call blktri(1,np,n3m,amphk,acphk,apphk
     1             ,mp,n2m,amphj,acphjj,apphj,my,yfis,ierror,w)
      do j=1,n2m
      do k=1,n3m
      dq(i,j,k)=yfis(j,k)
      enddo
      enddo
                else
      do j=1,n2m
      do k=1,n3m
      dph(i,j,k)=qcap(i,j,k)
      enddo
      enddo
                endif
                enddo
      if(n1m.gt.1) then
c
c   phi in wavenumber space by direct real fft
c
      n1mu=n1m-1
      do k=1,n3m
      do j=1,n2m
      do i=1,n1mh
      ip=2*i
      id=2*i-1
      xr(id,j)=qcap(i,j,k)
      xr(ip,j)=dq(i,j,k)
      enddo
      enddo
      call fft99(xr,work,trigx1,ifx1,1,m1+1,n1m,n2m,+1)
      do j=1,n2m
      dph(n1m,j,k)=xr(1,j)
      do i=1,n1mu
      is=i+1
      dph(i,j,k)=xr(is,j)
      enddo
      enddo
      enddo
                endif
      return
      end
c
c  ****************************** subrout phini  **********************
c
c   in this subr the coefficients of the poisson eq. for dph
c   are calculated this subr. is called only at the beginning
c   these are the coefficients entering in the fishpack
c
      subroutine phini
      include 'param.f' 
      dimension acphjj(m2)
      n1mh=n1m/2+1
c
c   tridiagonal matrix coefficients at each k1      
c    for the radial direction
c
      do jc=1,n2m
      jm=jmv(jc)
      jp=jpv(jc)
      a22icc=rc(jc)*jmc(jc)*dx2q/g2rc(jc)
      a22icp=rc(jp)*jpc(jc)*dx2q/g2rc(jp)
      ac2=-(a22icc+a22icp)
      ugmmm=1./rm(jc)/g2rm(jc)
      amphj(jc)=a22icc*ugmmm
      apphj(jc)=a22icp*ugmmm
      acphjj(jc)=-(amphj(jc)+apphj(jc))
      enddo
      do ic=1,n1mh
      do jc=1,n2m
      acphj(ic,jc)=acphjj(jc)-ak1(ic)/rm(jc)**2
      enddo
      enddo

      do ic=1,n1mh
c     write(6,*)' end phini coeff j=1,n2m'
      jc=1
c     write(6,101)jc,amphj(jc),acphj(ic,jc),apphj(jc)
      jc=n2m
c     write(6,101)jc,amphj(jc),acphj(ic,jc),apphj(jc)
  101 format(3x,i4,3x,3e12.5)
      enddo
c
c   tridiagonal matrix coefficients 
c    for the axial direction
c
      do kc=1,n3m
      km=kmv(kc)
      kp=kpv(kc)
      a33icc=kmc(kc)*dx3q/g3rc(kc)
      a33icp=kpc(kc)*dx3q/g3rc(kp)
      ugmmm=1./g3rm(kc)
      amphk(kc)=a33icc*ugmmm
      apphk(kc)=a33icp*ugmmm
      acphk(kc)=-(amphk(kc)+apphk(kc))
      enddo
c     write(6,*)' end phini coeff k=1,n3m'
      kc=1
c     write(6,101)kc,amphk(kc),acphk(kc),apphk(kc)
      kc=n3m
c     write(6,101)kc,amphk(kc),acphk(kc),apphk(kc)
c
c   initialization for fishpack
c
      an3=n3m
      ax=alog(an3)/alog(2.)
      k=ax+1
      L=2**(K+1)
      mw=(K-2)*L+K+5+MAX(2*n3m,6*n2m)
      np=1
      mp=1
      do j=1,n2m
      do k=1,n3m
      yfis(j,k)=1.
      enddo
      enddo
      my=m2
      call blktri(0,np,n3m,amphk,acphk,apphk
     1             ,mp,n2m,amphj,acphjj,apphj,my,yfis,ierror,w)
      return
      end
 
