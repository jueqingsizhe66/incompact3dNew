      subroutine qphywa(q,qtil,n)
c
c   permits to transform a quantity from the physical
c   to the wave number space
c
c     energy spectrum 
      include 'param.f'
      dimension q(ndv,m1,m2,m3)
      complex qtil(m1,(m2-1)/2+1,m3)
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/fft001/ifxx1(13),trigxx1(2*(m1-1))
      common/fft002/ifxx2(13),trigxx2(3*(m2-1)/2+1)
      common/fft003/ifxx3(13),trigxx3(2*(m3-1))
      common/numond/kx(m1),ky(m2),kz(m3),kkmax
      common/spequa/xa,wor,xa2,wor2,xr,work
      complex xa(m3-1,m1-1),wor(m3-1,m1-1)
      complex xa2(m1-1,m3-1),wor2(m1-1,m3-1)
      real xr(m2+1,m1-1),work(m2,m1-1)
c
      
      n2mh=n2m/2+1
c
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
       return
       end 
      subroutine spectre(q,qtil)
c
c   evaluates the energy spectrum
c
      include 'param.f'
      dimension q(ndv,m1,m2,m3)
      complex qtil(m1,(m2-1)/2+1,m3)
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/numond/kx(m1),ky(m2),kz(m3),kkmax
      common/speene/e(ndv,0:m3)
      dimension ei(ndv,0:m3)
c
c  three-dimensional spectrum
c
      do 1 n=1,3
      call qphywa(q,qtil,n)
      call enekco(qtil,n,ei)
       do kk=1,kkmax
       e(n,kk)=ei(n,kk)
       enddo
  1    continue
       return
       end 
      subroutine enekco(qtil,n,e)
c
c     energy spectrum  of one velocity component
c
      include 'param.f'
      complex qtil(m1,(m2-1)/2+1,m3)
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/numond/kx(m1),ky(m2),kz(m3),kkmax
      dimension e(ndv,0:m3)
c
      
      n2mh=n2m/2+1
c
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
         e(n,kk)=e(n,kk)+(urea*urea+uimm*uimm)
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
         e(n,kk)=e(n,kk)+(urea*urea+uimm*uimm)*0.5
         end do
       end do
       return
       end 
