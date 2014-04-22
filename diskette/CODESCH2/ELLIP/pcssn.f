c************************************************************************
c                                                                       *
c  ****************************** subrout tribj   ********************** *
c                                                                       *
c************************************************************************
c
c     ROUTINE FOR THE INVERSION OF BOUNDED TRIDIAGONAL MATRICES
c     by LU decomposition
c
      subroutine tribj(n,mi,mf,fk)
      include 'param.f'
      common/ctrdpj/amk(m2),ackj(m1,m2),apk(m2)
      dimension btj(m1),gmkj(m2,m1)
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
c   ********************* subr fftqui
c  this subroutine perform the calculation of trigz for temperton fft
c
      subroutine fftqui
      include 'param.f'
      common/dim/n1,n1m,n2,n2m
      common/fftcoi/ifxz(13),trigz(2*m1f)
      common/wavei/ani(m1)
      common/pi/pi
      pi=2.*asin(1.)
      n1f=2*n1m
c
c     modified wave number definition
c
      do j=1,n1
      ani(j)=(j-1)*pi
      enddo
      call cftfax(n1f,ifxz,trigz)
      return
      end
c
c   ********************* subr pscsn
c  this subroutine perform the calculation of dph , periodic direction
c  along x3 allows to use the cfft99
c
      subroutine pscns(qcap,psi)
      include 'param.f'
      common/dim/n1,n1m,n2,n2m
      dimension qki(m1,m2),qcap(m1,m2),psi(m1,m2)
      common/ctrdpj/amph(m2),acph(m1,m2),apph(m2)
      common/fftcoi/ifxz(13),trigz(2*m1f)
      real xr(2,m1fd,m2),work(2,m1f,m2)
      common/boupa/n1b1,n1bn1,n2b1,n2bn2
c
      n1f=n1m*2
c
c  fft-1 applied to  RHS of Laplacian by cfft99
c  how to store the quantities is explained in 
c  the NCAR Temperton routines
c
      q1ma=0.
      q2ma=0.
      q3ma=0.
      q4ma=0.
      do j=1,n2m
      do i=1,n1m
        xr(1,i,j)=qcap(i,j)
        xr(2,i,j)=0.
      q1ma=max(abs(qcap(i,j)),q1ma)
      enddo
      do i=n1,m1fd
          xr(2,i,j)=0.
          xr(1,i,j)=0.
      enddo
      enddo
      write(6,*)'n1b1=',n1b1,'q1ma= ',q1ma
      call cfft99(xr,work,trigz,ifxz,1,m1fd,n1f,n2,-1)
c
c   storage as in temperton  qki immag.
c
      if(n1b1.eq.2) then
      do j=2,n2m
      qki(1,j)=-xr(1,1,j)/float(n1m)
      do i=2,n1
      qki(i,j)=-xr(1,i,j)*2./float(n1m)
      q2ma=max(abs(qki(i,j)),q2ma)
      enddo
      enddo
      write(6,*)'n1b1=',n1b1,'q2ma= ',q2ma
                    endif
      if(n1b1.eq.1) then
      do j=2,n2m
      do i=1,n1
      qki(i,j)=-xr(2,i,j)*2./float(n1m)
      q2ma=max(abs(qki(i,j)),q2ma)
      enddo
      enddo
      write(6,*)'n1b1=',n1b1,'q2ma= ',q2ma
                    endif
      do j=1,n2,n2m
      do i=1,n1
      qki(i,j)=0.
      enddo
      enddo
c
c   solution of poisson equation immag. part
c
      ne=n1
      nt=n2
      call tribj(nt,1,ne,qki)
c
c   store the psi in the wave number space as in temperton cftt99
c
      if(n1b1.eq.2) then
      do j=1,n2
      do i=1,n1m
          xr(2,i,j)=0.
          xr(1,i,j)=-qki(i,j)
      q3ma=max(abs(qki(i,j)),q3ma)
      enddo
      enddo
      write(6,*)'n1b1=',n1b1,'q2ma= ',q3ma
                    endif
      if(n1b1.eq.1) then
      do j=1,n2
      do i=1,n1m
          xr(1,i,j)=0.
          xr(2,i,j)=-qki(i,j)
      q3ma=max(abs(qki(i,j)),q3ma)
      enddo
      enddo
      write(6,*)'n1b1=',n1b1,'q2ma= ',q3ma
                    endif
      do j=1,n2
      do i=n1,m1fd
          xr(1,i,j)=0.
          xr(2,i,j)=0.
      enddo
      enddo
c
c   fft +1 to have dph in phisycal space
c
      call cfft99(xr,work,trigz,ifxz,1,m1fd,n1f,n2,1)
c
c    the psi in physical space from cfft99 as in temperton
c
c    evaluation of psi  k from 2 to n2m-1
c
      do i=1,n1
        do  j=1,n2
        psi(i,j)=xr(1,i,j)
      q4ma=max(abs(psi(i,j)),q4ma)
      enddo
      enddo
      write(6,*)'pscns   ',q1ma,q2ma,q3ma,q4ma
      return
      end
c
c  ****************************** subrout matsbj **********************
c
c   in this subr the coefficients of the poisson eq. for psi
c   are calculated this subr. is called only at the beginning
c   inside phini.
c
      subroutine matsbj(k)
      include 'param.f'
      common/mesh/dx1,dx1q,dx2,dx2q
      common/dim/n1,n1m,n2,n2m
      common/wavei/ani(m1)
      common/ctrdpj/amph(m2),acph(m1,m2),apph(m2)
      common/boupa/n1b1,n1bn1,n2b1,n2bn2
c
c  ******** modified wave number
c
      ak1=2.*(1.-cos(ani(k)/float(n1m)))*dx1q
c
c   tridiagonal matrix coefficients at each k
c   Cartesian coordinates in x2
c
      do ic=2,n2m
        acph(k,ic)=-dx2q*2.-ak1
      enddo
        if(n2b1.eq.1) then
        acph(k,1)=1.
        acph(k,n2)=1.
                      endif
        if(n2b1.eq.2) then
        acph(k,1)=-dx2
        acph(k,n2)=dx2
                      endif
      return
      end
c
c   ********************* subr phinij
c  this subroutine calculate the l,u matrices to solve the Poisson
c  eq. for psi this subroutine is called only once.
c
      subroutine phinij
      include 'param.f'
      common/dim/n1,n1m,n2,n2m
      common/mesh/dx1,dx1q,dx2,dx2q
      common/ctrdpj/amph(m2),acph(m1,m2),apph(m2)
      common/boupa/n1b1,n1bn1,n2b1,n2bn2
c
      call fftqui
      do ic=2,n2m
        amph(ic)=dx2q
        apph(ic)=dx2q
      enddo
        if(n2b1.eq.1) then
        amph(1)=0.
        apph(1)=0.
        amph(n2)=0.
        apph(n2)=0.
                      endif
        if(n2b1.eq.2) then
        amph(1)=0.
        apph(1)=dx2
        amph(n2)=-dx2
        apph(n2)=0.
                      endif
      do  j=1,n1
        call matsbj(j)
      enddo
      j=1
      write(6,133)j,amph(j),acph(1,j),apph(j)
      j=(n2-1)/2+1
      write(6,133)j,amph(j),acph(1,j),apph(j)
      j=n2
      write(6,133)j,amph(j),acph(1,j),apph(j)
  133 format(1x,'pcssn  j=',i4,3x,4e12.5)
      return
      end
