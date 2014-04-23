c
c   ********************* subr ruuthe
c
c   CORRELATION IN THE AZIMUTHAL DIRECTION
c
c  this subroutine evaluates the energy spectrum 
c
      subroutine ruuthe(n2r)
      include 'param.f' 
      parameter (m1mh=m1m/2+1)
      dimension    xr(m1+1,m3m),work(m1m+1,m3m)
      dimension rkk(m1+1,m2),worj(m1,m2)
      dimension rip(m1m,m2)
      common/fftcm11/ifx11(13),trigx11(3*m1m/2+1)
      common/ifftin/iftin
      common/wavin/dlx1,dlx3,dkk1,dkk3
      n1mh=n1m/2+1
      n1md=n1m+2
      do j=1,n2r
      do k=1,n3m
      xr(1,k)=rhs(n1m,j,k)
      do i=1,n1m
      is=i+1
      xr(is,k)=rhs(i,j,k)
      enddo
      xr(n1md,k)=rhs(1,j,k)
      enddo
c
c   2-d   fft applied to the the generic fluctuating quantity
c   from physical to wave number space
c
      call fft99(xr,work,trigx11,ifx11,1,m1+1,n1m,n3m,-1)
c
c   azimuthal energy spectrum
c
            do ke=1,k1max
      enet(j,ke)=0.
            enddo
      do i=2,n1mh
      ip=2*i
      id=2*i-1
      ke=i
      rk=0.
                  do k=1,n3m
      enel=(xr(id,k)**2+xr(ip,k)**2)/n3m*2.
      enet(j,ke)=enet(j,ke)+enel*dlx1
c     rk=rk+enel  
                  enddo
c     rkk(id,j)=rk
c     rkk(ip,j)=0.
      enddo
c
c   rkkll  0 wave number   
c
      i=1
      ip=2*i
      id=2*i-1
      ke=i
      rk=0.
                  do k=1,n3m
      enel=(xr(id,k)**2+xr(ip,k)**2)/n3m
      enet(j,ke)=enet(j,ke)+enel*dlx1
                  enddo
      enddo
      return
      end
c
c   ********************* subr ruuaxi
c
c   CORRELATION IN THE AXIAL DIRECTION
c
c
c  this subroutine evaluates the energy spectrum and
c
      subroutine ruuaxi(n2r)
      include 'param.f' 
      parameter (m3mh=m3m/2+1,m3p=m3+1)
      dimension rkk(m3+1,m2),worj(m3,m2)
      dimension    xr(m3+1,m1m),work(m3m+1,m1m)
      dimension rip(m3m,m2)
      common/ifftin/iftin
      common/wavin/dlx1,dlx3,dkk1,dkk3
      common/fftcm3/ifx3(13),trigx3(3*m3m/2+1)
      n3mh=n3m/2+1
      n3md=n3m+2
      do j=1,n2r
      do i=1,n1m
      xr(1,i)=rhs(i,j,n3m)
      do k=1,n3m
      ks=k+1
      xr(ks,i)=rhs(i,j,k)
      enddo
      xr(n3md,i)=rhs(i,j,1)
      enddo

c
c   2-d   fft applied to the generic fluctuating quantity
c   from physical to wave number space
c
      call fft99(xr,work,trigx3,ifx3,1,m3+1,n3m,n1m,-1)
c
c
      do ke=1,k3max
      enez(j,ke)=0.
      enddo
      do k=2,n3mh
      kp=2*k
      kd=2*k-1
      ke=k
      do i=1,n1m
      enel=(xr(kd,i)**2+xr(kp,i)**2)/n1m*2.
      enez(j,ke)=enez(j,ke)+enel*dlx3
      enddo
      enddo
c
c    zero wave number
c
      k=1
      kp=2*k
      kd=2*k-1
      ke=k
      rk=0.
      do i=1,n1m
      enel=(xr(kd,i)**2+xr(kp,i)**2)/n1m
      enez(j,ke)=enez(j,ke)+enel*dlx3
      enddo
      enddo
      return
      end
c   ********************* subr fftqua
c  this subroutine perform the calculation of trigz for temperton fft
c
      subroutine fftqua
      include 'param.f'
      common/ifftin/iftin
      common/wavin/dlx1,dlx3,dkk1,dkk3
      common/fftcm3/ifx3(13),trigx3(3*m3m/2+1)
      common/fftcm11/ifx11(13),trigx11(3*m1m/2+1)
      pi=2.*asin(1.)
      n1mh=n1m/2
      n3mh=n3m/2
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
      dlx1=1.
      dlx3=alx3d/(2.*pi)
      dkk1=1.
      dkk3=2.*pi/alx3d
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
      iftin=1
      nx3fft=n3m
      nx1fft=n1m
      n2r=n2
      n2rm=n2m
      call fftfax(n3m,ifx3,trigx3)
      call fftfax(n1m,ifx11,trigx11)
      iftin=0
      return
      end
c
c   ********************* subr correl
c
c   CORRELATIONS  
c   these correlations are calculated in the physical
c   space and then require a large computational time
c   to speed up the FFT can be used
c   the user can write this more efficient way
c   by following what was done for the correlations
c   see routines cospx1 and cospx3 in
c   chatrasp.f in the directory 
c   CODESCH9/POST2
c
c   
      subroutine correl(n2r)
      include 'param.f'
      parameter (m1mh=m1m/2+1)
      parameter (m3mh=m3m/2+1)
      dimension riip(m1m,m2),rkkp(m3m,m2)
      dimension sum_j(m2)
      common/ruudi/ruuthd(m2,m1m),ruuaxd(m2,m3m)
c
        n1mh=n1m/2+1
        n3mh=n3m/2+1
c
      do   j=1,n2r
            sum_j(j)=0.
         do i=1,n1m
          riip(i,j)=0.
         end do
         do k=1,n3m
          rkkp(k,j)=0.
         end do
      end do
 
         
      do  1 j=1,n2r
c
c
      do k=1,n3m
       do i= 1,n1m
c
        sum_j(j)=sum_j(j)+ rhs(i,j,k)*rhs(i,j,k)
c
c correlation in x direction
c
        do jj= 0,n1m-1
         ii=i+jj
         if(ii.gt.n1m) ii=ii-n1m
         riip(jj+1,j)= riip(jj+1,j)+rhs(i,j,k)*rhs(ii,j,k)
        end do
c
c correlation in z direction
c
        do jj= 0,n3m-1
         kk=k+jj
         if(kk.gt.n3m) kk=kk-n3m
         rkkp(jj+1,j)= rkkp(jj+1,j)+rhs(i,j,k)*rhs(i,j,kk)
        end do
c
       end do
      end do

      sum_j(j)=sum_j(j)/float(n1m*n3m)

      do i= 1,n1m
       riip(i,j)= riip(i,j)/float(n1m*n3m)
      end do
      do k=1,n3m
       rkkp(k,j)= rkkp(k,j)/float(n1m*n3m)
      end do
c
c  average the correlations for positive 
c  and negative  separation
c
      do i= 2,n1m/2
        ii = n1m+2 -i
        riip(i,j)= .5*(riip(i,j)+ riip(ii,j))
      end do
      do k= 2,n3m/2
        kk = n3m+2 -k
        rkkp(k,j)= .5*(rkkp(k,j)+ rkkp(kk,j))
      end do
c
      do i= 1,n1m/2+1
        ruuthd(j,i)= riip(i,j)/riip(1,j)
      end do
      do k= 1,n3m/2+1
        ruuaxd(j,k)= rkkp(k,j)/rkkp(1,j)
      end do

   1  continue    
      return
      end
