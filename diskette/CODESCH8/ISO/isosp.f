c
c   in this routine the energy spectra are calculated
c
      subroutine spectre(q,qtil)
c     energy spectrum 
      include 'param.f'
      dimension q(ndv,m1,m2,m3)
      complex qtil((m1-1)/2+1,m2,m3)
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/fft002/ifxx2(13),trigxx2(2*(m2-1))
      common/fft001/ifxx1(13),trigxx1(3*(m1-1)/2+1)
      common/fft003/ifxx3(13),trigxx3(2*(m3-1))
      common/numond/kx(m1),ky(m2),kz(m3),kkmax
      common/speene/e(ndv,0:m3)
      common/sparr/xa,wor,xa2,wor2,xr,work
      complex xa(m3-1,m2-1),wor(m3-1,m2-1)
      complex xa2(m2-1,m3-1),wor2(m2-1,m3-1)
      real xr(m1+1,m2-1),work(m1,m2-1)
c
      
      n1mh=n1m/2+1
c
      do 1 n=1,3
      do  10 k=1,n3m
c
        do j=1,n2m
         xr(1,j)=q(n,n1m,j,k)
         xr(n1+1,j)=q(n,1,j,k)
         do i=1,n1m
          is=i+1
          xr(is,j)=q(n,i,j,k)
         enddo
        enddo
c
c   real fft applied from  physical space to wave number 
c   in x1 direction
c
        call fft99(xr,work,trigxx1,ifxx1,1,m1+1,n1m,n2m,-1)
c
        do i=1,n1mh
         ip=2*i   
         id=2*i-1
         do j=1,n2m
          qtil(i,j,k)=cmplx(xr(id,j),xr(ip,j))
         enddo
        enddo
c
 10   continue
c
c   complex  cfft applied (twice  x3 and then x2) from
c   physical space to wave number
c
      do 20 i=1,n1mh
c
        do k=1,n3m
         do j=1,n2m
          xa(k,j)=qtil(i,j,k)
         enddo
        enddo
c
        call cfft99(xa,wor,trigxx3,ifxx3,1,m3-1,n3m,n2m,-1)
c
        do k=1,n3m
         do j=1,n2m
          xa2(j,k)=xa(k,j)/float(n3m)
         enddo
        enddo
c
        call cfft99(xa2,wor2,trigxx2,ifxx2,1,m2-1,n2m,n3m,-1)
c
        do j=1,n2m
         do k=1,n3m
          qtil(i,j,k)=xa2(j,k)/float(n2m)
         enddo
        enddo
c
  20   continue
c
c    3D energy spectrum
c
       do kk=1,kkmax
       e(n,kk)=0.0
       enddo
       do k=1,n3m
        do i=2,n1mh
         do j=1,n2m
         akk=sqrt(float(kx(i)*kx(i)+ky(j)*ky(j)+kz(k)*kz(k)))
         kk=nint(akk)
         urea=real(qtil(i,j,k))
         uimm=aimag(qtil(i,j,k))
         e(n,kk)=e(n,kk)+2.*(urea*urea+uimm*uimm)
         end do
        end do
       end do
       i=1
       do k=1,n3m
         do j=1,n2m
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
