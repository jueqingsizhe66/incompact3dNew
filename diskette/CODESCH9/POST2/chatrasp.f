c
c   ********************* subr ruux1
c
c   CORRELATION IN THE spanwise DIRECTION
c  this subroutine evaluates the inverse fft in x1 of a general
c  quantity rhs 
c  it is stored in en1ik that is used  in the routine cospx3
c
c
      subroutine ruux1(j,en1ik,l)
      include 'param.f' 
      dimension    xr(m1+1,m3m),work(m1m+1,m3m)
      common/fftcm1/ifx1(13),trigx1(3*m1m/2+1)
      dimension en1ik(7,m1+1,m3m)
      n1mh=n1m/2+1
      n1md=n1m+2
      do k=1,n3m
      xr(1,k)=rhs(n1m,j,k)
      do i=1,n1m
      is=i+1
      xr(is,k)=rhs(i,j,k)
      enddo
      xr(n1md,k)=rhs(1,j,k)
      enddo
c
c   2D real  fft applied to the rhs by fft99 along x1
c   from physical to wave number space
c
      call fft99(xr,work,trigx1,ifx1,1,m1+1,n1m,n3m,-1)
c
c
      do i=1,n1mh
      ip=2*i
      id=2*i-1
                  do k=1,n3m
      en1ik(l,id,k)=xr(id,k)
      en1ik(l,ip,k)=xr(ip,k)
                  enddo
      enddo
      return
      end
c
c   ********************* subr ruux3
c  this subroutine evaluates the inverse fft in x3 of a general
c  quantity rhs 
c  it is stored in en3ki that is used  in the routine cospx3
c
c
      subroutine ruux3(j,en3ki,l)
      include 'param.f' 
      dimension    xr(m3+1,m1m),work(m3m+1,m1m)
      common/fftcm3/ifx3(13),trigx3(3*m3m/2+1)
      dimension en3ki(7,m3+1,m1m)
      n3mh=n3m/2+1
      n3md=n3m+2
      do i=1,n1m
      xr(1,i)=rhs(i,j,n3m)
      do k=1,n3m
      ks=k+1
      xr(ks,i)=rhs(i,j,k)
      enddo
      xr(n3md,i)=rhs(i,j,1)
      enddo
c
c   2D real  fft applied to the rhs by fft99 along x3
c   from physical to wave number space
c   
c
c
      call fft99(xr,work,trigx3,ifx3,1,m3+1,n3m,n1m,-1)
c
c
      do k=1,n3mh
      kp=2*k
      kd=2*k-1
            do i=1,n1m
      en3ki(l,kd,i)=xr(kd,i)
      en3ki(l,kp,i)=xr(kp,i)
                  enddo
      enddo
      return
      end
c
c   ********************* subr cospx3
c
c   spectra and cospectra calculation
c
c  this subroutine evaluates the cospectra  in the
c       STREAMWISE direction
c
      subroutine cospx3(j,en3ki,l,m)
      include 'param.f' 
      parameter (m3mh=m3m/2+1,m3p=m3+1)
      dimension en3ki(7,m3+1,m1m)
      dimension corr3(m3+1,1),work(m3,1)
      common/fftcm3/ifx3(13),trigx3(3*m3m/2+1)
      n3mh=n3m/2+1
      do ke=1,n3
      corr3(ke,1)=0.
      enddo
      do ke=1,n3m
      ene3(j,ke)=0.
      enddo
      do k=2,n3mh
      kp=2*k
      kd=2*k-1
      ke=k
            do i=1,n1m
      enel=(en3ki(l,kd,i)*en3ki(m,kd,i)
     1     +en3ki(l,kp,i)*en3ki(m,kp,i))/n1m*2.
      ene3(j,ke)=ene3(j,ke)+enel
      corr3(kd,1)=corr3(kd,1)+enel
      corr3(kp,1)=0.
            enddo
      enddo
      k=1
      kp=2*k
      kd=2*k-1
      ke=k
            do i=1,n1m
      enel=(en3ki(l,kd,i)*en3ki(m,kd,i)
     1     +en3ki(l,kp,i)*en3ki(m,kp,i))/n1m
      ene3(j,ke)=ene3(j,ke)+enel
      corr3(kd,1)=corr3(kd,1)+enel
      corr3(kp,1)=0.
            enddo
c
c   fft from the wave number to physical space
c
      call fft99(corr3,work,trigx3,ifx3,1,m3+1,n3m,1,+1)
      do k=1,n3m
      ks=k+1
      corx3(j,k)=corr3(ks,1)/corr3(2,1)
      enddo
c     write(62,133)j,l,m,ene3(j,1),ene3(j,2),corx3(j,1),corx3(j,2)
  133 format(3x,'  x3  ',3i5,2(2x,2e12.5))
      return
      end
c
c   ********************* subr cospx1
c
c   spectra and cospectra calculation
c
c  this subroutine evaluates the cospectra  in the
c       SPANWISE direction
c
      subroutine cospx1(j,en1ik,l,m)
      include 'param.f' 
      parameter (m3mh=m3m/2+1,m3p=m3+1)
      dimension en1ik(7,m1+1,m3m)
      common/fftcm1/ifx1(13),trigx1(3*m1m/2+1)
      dimension corr1(m1+1,1),work(m1,1)
      n1mh=n1m/2+1
      do ke=1,n1
      corr1(ke,1)=0.
      enddo
      do ke=1,n1m
      ene1(j,ke)=0.
      enddo
      do i=2,n1mh
      ip=2*i
      id=2*i-1
      ke=i
            do k=1,n3m
      enel=(en1ik(l,id,k)*en1ik(m,id,k)
     1     +en1ik(l,ip,k)*en1ik(m,ip,k))/n3m*2.
      ene1(j,ke)=ene1(j,ke)+enel
      corr1(id,1)=corr1(id,1)+enel
      corr1(ip,1)=0.
            enddo
      enddo
      i=1
      ip=2*i
      id=2*i-1
      ke=i
            do k=1,n3m
      enel=(en1ik(l,id,k)*en1ik(m,id,k)
     1     +en1ik(l,ip,k)*en1ik(m,ip,k))/n3m
      ene1(j,ke)=ene1(j,ke)+enel
      corr1(id,1)=corr1(id,1)+enel
      corr1(ip,1)=0.
            enddo
      call fft99(corr1,work,trigx1,ifx1,1,m1+1,n1m,1,+1)
      do i=1,n1m
      is=i+1
      corx1(j,i)=corr1(is,1)/corr1(2,1)
      enddo
c     write(62,133)j,l,m,ene1(j,1),ene1(j,2),corx1(j,1),corx1(j,2)
  133 format(3x,'  x1  ',3i5,2(2x,2e12.5))
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
      common/fftcm1/ifx1(13),trigx1(3*m1m/2+1)
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
      dlx3=1.
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
      nx3fft=n3m
      nx1fft=n1m
      n2r=n2
      n2rm=n2m
      call fftfax(n1m,ifx1,trigx1)
      call fftfax(n3m,ifx3,trigx3)
      return
      end
