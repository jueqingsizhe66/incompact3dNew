cc
c************************************************************
c****************  spectre   **********************************
c************************************************************
c
c
c   calculation of three-dimensional energy spectra
c
      subroutine spectre(q,qtil)
c     energy spectrum 
      include 'param.f'
      dimension q(ndv,m1,m2,m3)
      complex qtil(m1,(m2-1)/2+1,m3)
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/fft001/ifxx1(13),trigxx1(2*(m1-1))
      common/fft002/ifxx2(13),trigxx2(3*(m2-1)/2+1)
      common/fft003/ifxx3(13),trigxx3(2*(m3-1))
      common/numond/kx(m1),ky(m2),kz(m3),kkmax
      common/speene/e(ndv,0:m3)
      complex xa(m3-1,m1-1),wor(m3-1,m1-1)
      complex xa2(m1-1,m3-1),wor2(m1-1,m3-1)
      real xr(m2+1,m1-1),work(m2,m1-1)
c
      
      n2mh=n2m/2+1
c
      do 1 n=1,3
      do  10 k=1,n3m
c
        do i=1,n1m
         xr(1,i)=q(n,i,n2m,k)
         xr(n2m+2,i)=q(n,i,1,k)
         do j=1,n2m
          js=j+1
          xr(js,i)=q(n,i,j,k)
         enddo
        enddo
c
c   real fft applied from  physical space to wave number 
c   in y direction
c
        call fft99(xr,work,trigxx2,ifxx2,1,m2+1,n2m,n1m,-1)
c
        do j=1,n2mh
         jp=2*j   
         jd=2*j-1
         do i=1,n1m
          qtil(i,j,k)=cmplx(xr(jd,i),xr(jp,i))
         enddo
        enddo
c
 10   continue
c
c   2-d  cfft applied (twice) from
c   physical space to wave number
c
      do 20 j=1,n2mh
c
        do k=1,n3m
         do i=1,n1m
          xa(k,i)=qtil(i,j,k)
         enddo
        enddo
c
        call cfft99(xa,wor,trigxx3,ifxx3,1,m3-1,n3m,n1m,-1)
c
        do k=1,n3m
         do i=1,n1m
          xa2(i,k)=xa(k,i)/float(n3m)
         enddo
        enddo
c
        call cfft99(xa2,wor2,trigxx1,ifxx1,1,m1-1,n1m,n3m,-1)
c
        do i=1,n1m
         do k=1,n3m
          qtil(i,j,k)=xa2(i,k)/float(n1m)
         enddo
        enddo
c
  20   continue
       do kk=1,kkmax
       e(n,kk)=0.0
       enddo
       do k=1,n3m
        do j=2,n2mh
         do i=1,n1m
         akk=sqrt(float(kx(i)*kx(i)+ky(j)*ky(j)+kz(k)*kz(k)))
         kk=nint(akk)
         urea=real(qtil(i,j,k))
         uimm=aimag(qtil(i,j,k))
         e(n,kk)=e(n,kk)+2.*(urea*urea+uimm*uimm)
         end do
        end do
       end do
       j=1
       do k=1,n3m
         do i=1,n1m
         akk=sqrt(float(kx(i)*kx(i)+ky(j)*ky(j)+kz(k)*kz(k)))
         kk=nint(akk)
         urea=real(qtil(i,j,k))
         uimm=aimag(qtil(i,j,k))
         e(n,kk)=e(n,kk)+(urea*urea+uimm*uimm)
         end do
       end do
  1    continue
       return
       end 
c
c************************************************************
c****************  speini   **********************************
c************************************************************
c***********************************************************
c   here arrays for the FFT are calculated before to
c   perform the calculation
c
      subroutine speini
      include 'param.f'
      common/fft001/ifxx1(13),trigxx1(2*(m1-1))
      common/fft002/ifxx2(13),trigxx2(3*(m2-1)/2+1)
      common/fft003/ifxx3(13),trigxx3(2*(m3-1))
      common/ffc003/ifxc3(13),trigxc3(3*(m3-1)/2+1)
      common/ffc002/ifxc2(13),trigxc2(2*(m2-1))
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/numond/kx(m1),ky(m2),kz(m3),kkmax
      common/waves/an(m3),ap(m1),ak3(m3),ak1(m1)
      nx3fft=n3m
      call cftfax(nx3fft,ifxx3,trigxx3)
      call fftfax(nx3fft,ifxc3,trigxc3)
      nx1fft=n1m
      call cftfax(nx1fft,ifxx1,trigxx1)
      nx2fft=n2m
      call fftfax(nx2fft,ifxx2,trigxx2)
      call cftfax(nx2fft,ifxc2,trigxc2)
      pi=2.*asin(1.)
      n1mh=n1m/2
      n2mh=n2m/2
      n3mh=n3m/2
      n1mp=n1mh+1
      n2mp=n2mh+1
      n3mp=n3mh+1
c
c     wave number definition
c
      do k=1,n3mh
       kz(k)=(k-1)
      enddo
      do k=n3mp,n3m
       kz(k)=-(n3m-k+1)
      enddo
      do i=1,n1mh
      kx(i)=(i-1)
      enddo
      do i=n1mp,n1m
      kx(i)=-(n1m-i+1)
      enddo
      do j=1,n2mh
        ky(j)=(j-1)
      end do
      do j=n2mp,n2m
       ky(j)=-(n2m-j+1)
      end do
      akk=sqrt(float(n1mh*n1mh+n2mh*n2mh+n3mh*n3mh))
      kkmax=nint(akk)
      return 
      end
c************************************************************
c****************  correl   **********************************
c************************************************************
c***********************************************************
c   the correlation are calculated by using the FFT
c   this is similar to the routine for the spectrum 
c
      subroutine correl(q,qtil)
      include 'param.f'
      dimension q(ndv,m1,m2,m3)
      complex qtil(m1,m2,(m3-1)/2+1)
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/fft001/ifxx1(13),trigxx1(2*(m1-1))
      common/ffc003/ifxc3(13),trigxc3(3*(m3-1)/2+1)
      common/ffc002/ifxc2(13),trigxc2(2*(m2-1))
      common/numond/kx(m1),ky(m2),kz(m3),kkmax
      common/corre/corf(ndv,m3)
      complex xa(m2-1,m1-1),wor(m2-1,m1-1)
      complex xa2(m1-1,m2-1),wor2(m1-1,m2-1)
      real xr(m3+1,m1-1),work(m3,m1-1)
      real cor(m3+1,1),worcok(m3,1)
c
      
      n3mh=n3m/2+1
c
      do 1 n=1,3
      do  10 j=1,n2m
c
        do i=1,n1m
         xr(1,i)=q(n,i,j,n3m)
         xr(n3m+2,i)=q(n,i,j,1)
         do k=1,n3m
          ks=k+1
          xr(ks,i)=q(n,i,j,k)
         enddo
        enddo
c
c   real fft applied from  physical space to wave number 
c   in y direction
c
        call fft99(xr,work,trigxc3,ifxc3,1,m3+1,n3m,n1m,-1)
c
        do k=1,n3mh
         kp=2*k   
         kd=2*k-1
         do i=1,n1m
          qtil(i,j,k)=cmplx(xr(kd,i),xr(kp,i))
         enddo
        enddo
c
 10   continue
c
c   2-d  cfft applied (twice) from
c   physical space to wave number
c
      do 20 k=1,n3mh
c
        do j=1,n2m
         do i=1,n1m
          xa(j,i)=qtil(i,j,k)
         enddo
        enddo
c
        call cfft99(xa,wor,trigxc2,ifxc2,1,m2-1,n2m,n1m,-1)
c
        do j=1,n2m
         do i=1,n1m
          xa2(i,j)=xa(j,i)/float(n2m)
         enddo
        enddo
c
        call cfft99(xa2,wor2,trigxx1,ifxx1,1,m1-1,n1m,n2m,-1)
c
        do i=1,n1m
         do j=1,n2m
          qtil(i,j,k)=xa2(i,j)/float(n1m)
         enddo
        enddo
c
  20   continue
c
c   the sum of the k1 k2 wave numbers is done
c   the results is a real array to be transformed
c   from wave number k3  to x3
c
       anpo=1./float(n1m*n2m)
       do k=1,n3mh
       kp=2*k
       kd=2*k-1
       cor(kd,1)=0.0
       cor(kp,1)=0.0
       do i=1,n1m
      do  j=1,n2m
         urea=real(qtil(i,j,k))
         uimm=aimag(qtil(i,j,k))
         cor(kp,1)=0.
         cor(kd,1)=cor(kd,1)+(urea*urea+uimm*uimm)
         end do
        end do
         cor(kd,1)=cor(kd,1)*anpo
       end do
c
c   real fft applied from  physical space to wave number 
c   in x3 direction
c
        call fft99(cor,worcok,trigxc3,ifxc3,1,m3+1,n3m,1,+1)
c
       do k=1,n3m
         ks=k+1
         corf(n,k)=cor(ks,1)
       end do
  1    continue
       return
       end 
