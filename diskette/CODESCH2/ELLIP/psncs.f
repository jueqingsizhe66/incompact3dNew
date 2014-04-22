c************************************************************************
c                                                                       *
c  ****************************** subrout tribi   ********************** *
c                                                                       *
c************************************************************************
c
c     ROUTINE FOR THE INVERSION OF BOUNDED TRIDIAGONAL MATRICES
c     by LU decomposition
c
      subroutine tribi(n,mi,mf,fk)
      include 'param.f'
      common/ctrdpi/amk(m1),ackj(m2,m1),apk(m1)
      dimension btj(m2),gmkj(m2,m1)
      dimension fk(m2,m1)
      do i=mi,mf
      btj(i)=ackj(i,1)
      fk(i,1)=fk(i,1)/btj(i)
      enddo
      do j=2,n
      do i=mi,mf
      gmkj(i,j)=apk(j-1)/btj(i)
      btj(i)=ackj(i,j)-amk(j)*gmkj(i,j)
      fk(i,j)=(fk(i,j)-amk(j)*fk(i,j-1))/btj(i)
      enddo
      enddo
      do j=n-1,1,-1
      do i=mi,mf
      fk(i,j)=fk(i,j)-gmkj(i,j+1)*fk(i,j+1)
      enddo
      enddo
      return
      end
c
c   ********************* subr fftquj
c  this subroutine perform the calculation of trigz for temperton fft
c
      subroutine fftquj
      include 'param.f'
      common/dim/n1,n1m,n2,n2m
      common/fftcoj/ifxz(13),trigz(2*m2f)
      common/wavej/anj(m2)
      common/pi/pi
      pi=2.*asin(1.)
      n2mh=n2m/2+1
      n2f=2*n2m
c     modified wave number definition
c
      do j=1,n2
      anj(j)=(j-1)*pi
      enddo
      call cftfax(n2f,ifxz,trigz)
      return
      end
c
c   ********************* subr pnssc
c  this subroutine perform the calculation of dph , periodic direction
c  along x3 allows to use the cfft99
c
      subroutine pnssc(qcap,psi)
      include 'param.f'
      common/dim/n1,n1m,n2,n2m
      dimension qki(m2,m1),qcap(m1,m2),psi(m1,m2)
      common/ctrdpi/amph(m1),acph(m2,m1),apph(m1)
      common/fftcoj/ifxz(13),trigz(2*m2f)
      common/boupa/n1b1,n1bn1,n2b1,n2bn2
      real xr(2,m2fd,m1),work(2,m2f,m1)
c
      n2f=n2m*2
c
c  fft-1 applied to  RHS of Laplacian by cfft99
c  how to store the quantities is explained in the ncar temperton
c  routine
c
      q1ma=0.
      q2ma=0.
      q3ma=0.
      q4ma=0.
      do i=1,n1
      do j=1,n2m
        xr(1,j,i)=qcap(i,j)
        xr(2,j,i)=0.
      q1ma=max(abs(qcap(i,j)),q1ma)
      enddo
      do j=n2,m2fd
          xr(2,j,i)=0.
          xr(1,j,i)=0.
      enddo
      enddo
      write(6,*)'n2b1=',n2b1,'q1ma= ',q1mA
      call cfft99(xr,work,trigz,ifxz,1,m2fd,n2f,n1,-1)
c
c
c   storage as in temperton  qki immag.
c
      if(n2b1.eq.1) then
      do i=2,n1m
      do k=1,n2
      qki(k,i)=-xr(2,k,i)*2./float(n2m)
      q2ma=max(abs(qki(k,i)),q2ma)
      enddo
      enddo
      write(6,*)'n2b1=',n2b1,q2mA
                    endif
      if(n2b1.eq.2) then
      do i=2,n1m
      qki(1,i)=-xr(1,1,i)/float(n2m)
      do k=2,n2
      qki(k,i)=-xr(2,k,i)*2./float(n2m)
      q2ma=max(abs(qki(k,i)),q2ma)
      enddo
      enddo
      write(6,*)'n2b1=',n2b1,q2mA
                    endif
      do i=1,n1,n1m
      do k=1,n2
      qki(k,i)=0.
      enddo
      enddo
c
c   solution of poisson equation immag. part
c
      ne=n2
      nt=n1
      call tribi(nt,1,ne,qki)
c
c   store the psi in the wave number space as in temperton cftt99
c
      if(n2b1.eq.1) then
      do i=1,n1
      do k=1,n2m
          xr(1,k,i)=0.
          xr(2,k,i)=-qki(k,i)
      q3ma=max(abs(qki(k,i)),q3ma)
      enddo
      enddo
      write(6,*)'n2b1=',n2b1,q3mA
                    endif
      if(n2b1.eq.2) then
      do i=1,n1
      do k=1,n2m
          xr(2,k,i)=0.
          xr(1,k,i)=-qki(k,i)
      q3ma=max(abs(qki(k,i)),q3ma)
      enddo
      enddo
      write(6,*)'n2b1=',n2b1,q3mA
                    endif
      do i=1,n1
      do k=n2,m2fd
          xr(1,k,i)=0.
          xr(2,k,i)=0.
      enddo
      enddo
c
c   fft +1 to have dph in phisycal space
c
      call cfft99(xr,work,trigz,ifxz,1,m2fd,n2f,n1,1)
c
c    the psi in physical space from cfft99 as in temperton
c
c    evaluation of psi  k from 2 to n2m-1
c
      do i=1,n1
        do k=1,n2
      q4ma=max(abs(xr(1,k,i)),q4ma)
         psi(i,k)=xr(1,k,i)
      enddo
      enddo
      write(6,*)'pnssc   ',q1ma,q2ma,q3ma,q4ma
      return
      end
c
c  ****************************** subrout matsbi **********************
c
c   in this subr the coefficients of the poisson eq. for psi
c   are calculated this subr. is called only at the beginning
c   inside phini.
      subroutine matsbi(k)
      include 'param.f'
      common/mesh/dx1,dx1q,dx2,dx2q
      common/dim/n1,n1m,n2,n2m
      common/wavej/anj(m2)
      common/ctrdpi/amph(m1),acph(m2,m1),apph(m1)
      common/boupa/n1b1,n1bn1,n2b1,n2bn2
c
c  ******** modified wave number
c
      ak2=2.*(1.-cos(anj(k)/float(n2m)))*dx2q
c
c   tridiagonal matrix coefficients at each k
c   and cartesian coordinates in x2
c
      do i=2,n1m
        acph(k,i)=-dx1q*2.-ak2
      enddo
        if(n1b1.eq.1) then
        acph(k,1)=1.
        acph(k,n1)=1.
                      endif
        if(n1b1.eq.2) then
        acph(k,1)=-dx1
        acph(k,n1)=dx1
                      endif
      return
      end
c
c   ********************* subr phinii
c  this subroutine calculate the l,u matrices to solve the poisson
c  eq. for dph this subroutine is called only once.
c
      subroutine phinii
      include 'param.f'
      common/dim/n1,n1m,n2,n2m
      common/mesh/dx1,dx1q,dx2,dx2q
      common/ctrdpi/amph(m1),acph(m2,m1),apph(m1)
      common/boupa/n1b1,n1bn1,n2b1,n2bn2
c
      call fftquj
      do ic=2,n1m
        amph(ic)=dx1q
        apph(ic)=dx1q
      enddo
        if(n1b1.eq.1) then
        amph(1)=0.
        apph(1)=0.
        amph(n1)=0.
        apph(n1)=0.
                      endif
        if(n1b1.eq.2) then
        amph(1)=0.
        apph(1)=dx1
        amph(n1)=-dx1
        apph(n1)=0.
                      endif
      do j=1,n2
        call matsbi(j)
      enddo
      i=1
      write(6,132)i,amph(i),acph(1,i),apph(i)
      i=(n1-1)/2+2
      write(6,132)i,amph(i),acph(1,i),apph(i)
      i=n1
      write(6,132)i,amph(i),acph(1,i),apph(i)
  132 format(1x,'psncs i=',i4,3x,4e12.5)
      return
      end
