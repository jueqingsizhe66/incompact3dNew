c
c
c  ****************************** subrout trib  **********************
c   this routine evaluates the tridiagonal matrix
c   for the psi calculation. Can be substituted by the other 
c   tridiagonal solvers in  psomgbcnn.f 
      subroutine trib(a,b,c,r,n,u)
      include 'param.f'
      dimension gam(m1),a(m1),b(m1),c(m1),r(m1),u(m1)
      bet=b(1)
      u(1)=r(1)/bet
      do 11 j=2,n
      gam(j)=c(j-1)/bet
      bet=b(j)-a(j)*gam(j)
c      if(bet.eq.0) write(6,781)
  781 format(3x,'error in tridiag')
      u(j)=(r(j)-a(j)*u(j-1))/bet
   11 continue
      do 12 j=n-1,1,-1
      u(j)=u(j)-gam(j+1)*u(j+1)
   12 continue
      return
      end
c
c   ********************* subr dsolc
c
c   this subroutine performs the solution of the poisson equation
c   by solving l*u*ph=divg(q)
c   qk is divg(qk) stored as a vector
c   the calculation must be done at each modified wave number
c   the l, and u matrices have been stored in the asbn2 matrix
c   some memory could be furtherly reduced due to the large number of
c   zeros.
c
      subroutine dsolc(qk,n2mh)
      include 'param.f'
      parameter (m2m=m2-1,m2mh=m2m/2+1)
      dimension qsb(m1),am(m1),ac(m1),ap(m1),f(m1)
      dimension qk(m1,m2)
      common/ctrdpc/amph(m1),acph(m1,m2),apph(m1)
      common/dim/n1,n1m,n2,n2m
c
c
      do 11 j=1,n2mh
        do 15 ir=1,n1
          f(ir) = qk(ir,j)
          ac(ir)= acph(ir,j)
          ap(ir)= apph(ir)
          am(ir)= amph(ir)
   15   continue
c
c   tridiagonal inversion
c
        call trib(am,ac,ap,f,n1,qsb)
c
c  solution
c
        do 14 ij=1,n1
          qk(ij,j)=qsb(ij)
   14   continue
   11 continue
      return
      end
c
c   ********************* subr dsolr
c
c   this subroutine performs the solution of the poisson equation
c   by solving l*u*ph=divg(q)
c   qk is divg(qk) stored as a vector
c   the calculation must be done at each modified wave number
c   the l, and u matrices have been stored in the asbn2 matrix
c   some memory could be furtherly reduced due to the large number of
c   zeros.
c
      subroutine dsolr(qk,n2mh)

      include 'param.f'
      parameter (m2m=m2-1,m2mh=m2m/2+1)
      dimension qsb(m1),am(m1),ac(m1),ap(m1),f(m1)
      dimension qk(m1,m2mh)
      common/ctrdpr/amph(m1),acph(m1,m2mh),apph(m1)
      common/dim/n1,n1m,n2,n2m
c
c
      do 11 j=1,n2mh
        do 15 ir=1,n1
          f(ir) = qk(ir,j)
          ac(ir)= acph(ir,j)
          ap(ir)= apph(ir)
          am(ir)= amph(ir)
   15   continue
c
c   tridiagonal inversion
c
        call trib(am,ac,ap,f,n1,qsb)
c
c  solution
c
        do 14 ij=1,n1
          qk(ij,j)=qsb(ij)
   14   continue
   11 continue
      return
      end
c
c   ********************* subr fftqua
c  this subroutine perform the calculation of trigz for temperton fft
c
      subroutine fftqua
      include 'param.f'
      parameter (m2m=m2-1)
      parameter (m2f=m2m*2)
      common/dim/n1,n1m,n2,n2m
      common/fftcoc/ifxc(13),trigc(2*m2f)
      common/fftcor/ifxr(13),trigr(3*m2m/2+1)
      common/iperi/ib2per
      common/waver/anr(m2)
      common/wavec/anc(m2)
      common/pi/pi
      n2mh=n2m/2+1
      n2ff=2*n2m
c
c     modified wave number definition
c
c     write(6,764) (an(j),j=1,n2)
  764 format(1x,'an',2x,11f8.3)
      if(ib2per.eq.1) then
      do j=1,n2mh
      anr(j)=(j-1)*2.*pi
      enddo
      call fftfax(n2m,ifxr,trigr)
                      else
      do j=1,n2
      anc(j)=(j-1)*pi
      enddo
      call cftfax(n2ff,ifxc,trigc)
                      endif
      return
      end
c
c
c   ********************* subr phcalc
c  this subroutine perform the calculation of dph , periodic direction
c  along x3 allows to use the cfft99
c
      subroutine phcalc(qcap,psi)
      include 'param.f'
      parameter (m1m=m1-1,m2m=m2-1)
      parameter (m2f=m2m*2,m2fd=m2f+2)
      dimension qki(m1,m2),qcap(m1,m2),psi(m1,m2)
      common/dim/n1,n1m,n2,n2m
      common/fftcoc/ifxc(13),trigc(2*m2f)
      common/tstep/dt
      common/bvedph/dpest(m2),dpwes(m2)
      common/bchodp/dpsud(m1),dpnor(m1)
c
c
      real xr(2,m2fd,m1),work(2,m2f,m1)
c
      n2ff=n2m*2
c
c  fft-1 applied to  vor(i,j) by cfft99
c  how to store the quantities is explained in the ncar temperton
c  routine
c
      do 1 i=1,n1 
      do 2 j=2,n2m
        xr(1,j,i)=-qcap(i,j)
        xr(2,j,i)=0.
    2 continue
      do 52 j=n2,m2fd
          xr(2,j,i)=0.
          xr(1,j,i)=0.
   52   continue
    1 continue
c
c   boundary conditions inlet radiative or free-slip
c
        do 53 j=2,n2m
          xr(2,j,1)=0.
          xr(1,j,1)=dpwes(j)
   53   continue
c
c   boundary conditions outflow radiative or zero
c
        do 54 j=2,n2m
          xr(2,j,n1)=0.
          xr(1,j,n1)=dpest(j)
   54   continue
        do 55 j=n2,m2fd
          xr(2,j,1)=0.
          xr(1,j,1)=0.
          xr(2,j,n1)=0.
          xr(1,j,n1)=0.
   55   continue
c
c   boundary conditions wall or symmetry psi=0.
c
      do 56 i=1,n1
C        xr(1,1,i)=dpsud(i)
         xr(1,1,i)=0.
         xr(2,1,i)=0.
C        xr(1,n2,i)=dpnor(i)
         xr(1,n2,i)=0.
         xr(2,n2,i)=0.
   56   continue
      call cfft99(xr,work,trigc,ifxc,1,m2fd,n2ff,n1,-1)
c
c
c   storage as in temperton  qki immag.
c
      do 4 k=1,n2
      do 3 i=1,n1
      qki(i,k)=-xr(2,k,i)*2./float(n2m)
    3 continue
    4 continue
c
c   solution of poisson equation immag. part
c
      call dsolc(qki,n2)
c
c   store the psi in the wave number space as in temperton cftt99
c
      do 13 i=1,n1
      do 11 k=1,n2
          xr(1,k,i)=0.
          xr(2,k,i)=-qki(i,k)
   11 continue
      do 31 k=n2,m2fd
          xr(1,k,i)=0.
          xr(2,k,i)=0.
   31 continue
   13 continue
c
c   fft +1 to have dph in phisycal space
c
      call cfft99(xr,work,trigc,ifxc,1,m2fd,n2ff,n1,1)
c
c    the psi in physical space from cfft99 as in temperton
c
c    evaluation of psi  k from 2 to n2m-1
c
      do i=1,n1
c        psi(i,n2)=0.
c        psi(i,1)=0.
        do k=2,n2m
       psi(i,k)=xr(1,k,i)
      enddo
      enddo
c
c   boundary conditions 
c
c       do j=2,n2m
c         psi(1,j)  = dpwes(j)
c         psi(n1,j) = dpest(j)
c       end do
      return
      end
c
c
c   ********************* subr phcalp
c  this subroutine perform the calculation of dph , periodic direction
c  along x3 allows to use the cfft99
c
      subroutine phcalp(qcap,psi)
      include 'param.f'
      parameter (m1m=m1-1,m2m=m2-1)
      parameter (m2md=m2m+2,m2mh=m2m/2+1)
      dimension qki(m1,m2mh),qkr(m1,m2mh)
      dimension qcap(m1,m2),psi(m1,m2)
      common/dim/n1,n1m,n2,n2m
      common/fftcor/ifxr(13),trigr(3*m2m/2+1)
      common/tstep/dt
      common/bvedph/dphest(m2),dphwes(m2)
c
c
      real xr(m2md,m1),work(m2,m1)
c
c
c  fft-1 applied to  vor(i,j) by cfft99
c  how to store the quantities is explained in the ncar temperton
c  routine
c
      do i=2,n1m
        xr(1,i)=-qcap(i,n2m)
      do j=1,n2m
        js=j+1
        xr(js,i)=-qcap(i,j)
      enddo
      enddo
      do i=1,n1
      do j=n2+1,m2md
          xr(j,i)=0.
      enddo
      enddo
c
c   boundary conditions inlet radiative or free-slip
c
          xr(1,1)=dphwes(n2m)
        do j=1,n2m
        js=j+1
          xr(js,1)=dphwes(j)
        enddo
c
c   boundary conditions outflow radiative or zero
c
          xr(1,n1)=dphest(n2m)
        do j=1,n2m
        js=j+1
          xr(js,n1)=dphest(j)
        enddo
      call fft99(xr,work,trigr,ifxr,1,m2md,n2m,n1,-1)
c
c
c   storage as in temperton  qki immag.
c
      n2mh=n2m/2+1
      do kk=1,n2mh
      kd=2*kk-1
      kp=2*kk
      do i=1,n1
      qki(i,kk)=xr(kp,i)
      qkr(i,kk)=xr(kd,i)
      enddo
      enddo
c
c   solution of poisson equation immag. part
c
      call dsolr(qki,n2mh)
      call dsolr(qkr,n2mh)
c
c   store the psi in the wave number space as in temperton cftt99
c
      do kk=1,n2mh
      kd=2*kk-1
      kp=2*kk
      do i=1,n1
      xr(kp,i)=qki(i,kk)
      xr(kd,i)=qkr(i,kk)
      enddo
      enddo
c
c   fft +1 to have dph in phisycal space
c
      call fft99(xr,work,trigr,ifxr,1,m2md,n2m,n1,1)
c
c    the psi in physical space from cfft99 as in temperton
c
c    evaluation of psi  k from 2 to n2m-1
c
      do i=1,n1
          psi(i,n2m)=xr(1,i)
        do k=1,n2m-1
          ks=k+1
          psi(i,k)=xr(ks,i)
        enddo
      enddo
      return
      end
c
c  ****************************** subrout matsc  **********************
c
c   in this subr the coefficients of the poisson eq. for psi
c   are calculated this subr. is called only at the beginning
c   inside phini.
c
      subroutine matsc(k)
      include 'param.f'
      parameter (m2m=m2-1)
      common/mesh/dx1,dx1q,dx2,dx2q
      common/metri/g2m(m1),g2c(m1)
      common/dim/n1,n1m,n2,n2m
      common/wavec/an(m2)
      common/ctrdpc/amph(m1),acph(m1,m2),apph(m1)
      if(k.eq.1) then
      do ic=2,n1m
        amph(ic)=dx1q/g2m(ic-1)/g2c(ic)
        apph(ic)=dx1q/g2m(ic)/g2c(ic)
      enddo
      do ic=1,n1,n1m
        amph(ic)=0.
        apph(ic)=0.
      enddo
                  endif
c
c  ******** modified wave number
c
      ak2=2.*(1.-cos(an(k)/float(n2m)))*dx2q
c
c   tridiagonal matrix coefficients at each k
c    cartesian coordinates in x2
c
      do 1 ic=2,n1m
        acph(ic,k)=-dx1q*(1./g2m(ic-1)+1./g2m(ic))/g2c(ic)-ak2
    1 continue
      do 2 ic=1,n1,n1m
        acph(ic,k)=1.
    2 continue
      return
      end
c
c  ****************************** subrout matsr  **********************
c
c   in this subr the coefficients of the poisson eq. for psi
c   are calculated this subr. is called only at the beginning
c   inside phini.
c
      subroutine matsr(k)
      include 'param.f'
      parameter (m2m=m2-1,m2mh=m2m/2+1)
      common/mesh/dx1,dx1q,dx2,dx2q
      common/metri/g2m(m1),g2c(m1)
      common/dim/n1,n1m,n2,n2m
      common/waver/an(m2)
      common/ctrdpr/amph(m1),acph(m1,m2mh),apph(m1)
      if(k.eq.1) then
      do ic=2,n1m
        amph(ic)=dx1q/g2m(ic-1)/g2c(ic)
        apph(ic)=dx1q/g2m(ic)/g2c(ic)
      enddo
      do ic=1,n1,n1m
        amph(ic)=0.
        apph(ic)=0.
      enddo
                  endif
c
c  ******** modified wave number
c
      ak2=2.*(1.-cos(an(k)/float(n2m)))*dx2q
c
c   tridiagonal matrix coefficients at each k
c    cartesian coordinates in x2
c
      do 1 ic=2,n1m
        acph(ic,k)=-dx1q*(1./g2m(ic-1)+1./g2m(ic))/g2c(ic)-ak2
    1 continue
      do 2 ic=1,n1,n1m
        acph(ic,k)=1.
    2 continue
      return
      end
c
c   ********************* subr phini
c  this subroutine calculate the l,u matrices to solve the poisson
c  eq. for dph this subroutine is called only once.
c  correspond to the operation ijob=1 in the leqt1b imsl routine.
c
      subroutine phini
      include 'param.f'
      common/dim/n1,n1m,n2,n2m
      common/mesh/dx1,dx1q,dx2,dx2q
      common/indpe2/n2i,n2f
      common/iperi/ib2per
c
      call fftqua
c
c
c
      if(ib2per.eq.0) then
      do j=1,n2
        call matsc(j)
      enddo
                      else
      n2ft=n2m/2+1
      do j=1,n2ft
        call matsr(j)
      enddo
                      endif
      return
      end
