c
c  ****************************** subrout initur **********************
c
c   initial  conditions for the simulation of isotropic
c   turbulence.
c   The spectrum is assigned and from this the velocity field
c   is calculated See the Rogallo report.
c
      subroutine initur(q,pr,qcap,dph,time)
c
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      dimension q(ndv,m1,m2,m3),pr(m1,m2,m3)
      dimension qcap(m1,m2,m3),dph(m1,m2,m3)
      common/corpri/yp1(m1),yp2(m2),yp3(m3)
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/tstep/dt,beta,ren
      common/rhs3p/dp3ns
      common/speene/e(ndv,0:m3)
      common/numond/kx(m1),ky(m2),kz(m3),kkmax
      common/ispec/imic
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      pi=4.*atan(1.)
C
c call mic to make initial conditions. Use qcap,dph and pr as work files
C
      call turbis(qcap,dph,pr,q,time)
      open(58,file='spe3dcheco.out',status='unknown')
c
c    here the spectrum of the field is evaluated
c
      call spectre(q,qcap)
       do k=1,kkmax
       write (58,*)k,e(1,k),e(2,k),e(3,k)
       end do
      close(58)
        al=1.
c
c   the velocity field in physical space is not
c   free divergent
c
      call divg(qcap,q,al)
      call divgck(q,qmax)
      write(6,*)'  the divergence of the field is qmax=',qmax
c
c qcap now contains the divergence of the velocity field q
c
c  ********* calculation of the "pressure" dph by fft in two
c            directions and tridiag in vertical
      call phcalc(qcap,dph)
c
c  ********* calculation of solenoidal vel field
c
      call updvp(dph,q,1.)
      vm1m=0.
      vm2m=0.
      vm3m=0.
      vl123=1./float(n1m*n2m*n3m)
      do 411 k=1,n3m
      do 411 j=1,n2m
      do 411 i=1,n1m
       vm1m=q(1,i,j,k)+vm1m
       vm2m=q(2,i,j,k)+vm2m
       vm3m=q(3,i,j,k)+vm3m
  411 continue
      vm1=vm1m*vl123
      vm2=vm2m*vl123
      vm3=vm3m*vl123
      enej=0.
      print*,'vm1,vm2,vm3'
      print*,vm1,vm2,vm3
      call divgck(q,qmax)
      print *,'divergenza massima=',qmax
      open(57,file='spe3dinto.out',status='unknown')
      open(58,file='spe3dinco.out',status='unknown')
c
c   spectrum of the free-divergent field
c
      call spectre(q,qcap)
       do k=1,kkmax
       write (58,*)k,e(1,k),e(2,k),e(3,k)
       write (57,*)k,e(1,k)+e(2,k)+e(3,k)
       end do
      close(57)
      close(58)
      return
      end
*************************************************************************
      subroutine turbis(qtil1,qtil2,qtil3,q,time)
c
c make intitial conditions
c from assigned spectra depending on imic
c
      include 'param.f'
      dimension q(ndv,m1,m2,m3)
      complex qtil1((m1-1)/2+1,m2,m3)
      complex qtil2((m1-1)/2+1,m2,m3)
      complex qtil3((m1-1)/2+1,m2,m3)
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/numond/kx(m1),ky(m2),kz(m3),kkmax
      common/rot/f0
      common/strat1/istrat,rho0,g,schm
      common/strat2/igrad,bvais
      common/d123/alx1,alx2,alx3
      complex alp,bet,cexp1,cexp2,cexp3,cexp4
      dimension ekk(2*m1),e(3,2*m1)
      common/ispec/imic
      common/spemw/akkpp,qq,sig
c
      n2mh=n2m/2+1
      n1mh=n1m/2+1

      pi=4.*atan(1.)

      if(imic.eq.1)  then
c
c++++++++++COMPTE BELOLOT - Corrsin DATA +++++++++++++++
c      ek(0)=0.
       ekk(1)=.007
       ekk(2)=.048
       ekk(3)=.068
       ekk(4)=.099
       ekk(5)=.09
       ekk(6)=.068
       ekk(7)=.059
       ekk(8)=.05
       ekk(9)=.052
       ekk(10)=.041
       ekk(11)=.036
       ekk(12)=.034
       ekk(13)=.030
       ekk(14)=.0275
       ekk(15)=.0225
       ekk(16)=.024
       ekk(17)=.0230
       ekk(18)=.021
       ekk(19)=.018
       ekk(20)=.0175
       ekk(21)=.017 
       ekk(22)=.0155
       ekk(23)=.0140
       ekk(24)=.013
       ekk(25)=.012
       ekk(26)=.0115
       ekk(27)=.011
       ekk(28)=.0105
       ekk(29)=.0095
       ekk(30)=.009
       ekk(31)=.008
       ekk(32)=.0075
                      endif
       if(imic.eq.3)  then
       open(21,file='speread.dat',status='unknown')
      do kk=1,500
      read(21,*,end=211)k,e(1,k),e(2,k),e(3,k)
      ekk(kk)=e(1,kk)+e(2,kk)+e(3,kk)
      enddo
  211 continue
                      endif
       if(imic.eq.2)  then
c      
c  here we try  the initial spectrum of Mansour and Wray
c  see their Physics of Fluids paper
c               DNS
       aa=0.
       do kk=1,kkmax
       ak=kk
       rk=ak
       aint=rk**sig*exp(-sig*rk**2/2.)
       aa=aa+aint
       enddo
       aa1=0.
       eto=0.
       do kk=1,kkmax
       ak=kk
       rk=ak/akkpp
       ek0kl=qq/(2.*aa)/akkpp
       ek1kl=rk**sig*exp(-sig*rk**2/2.)
       aa1=aa1+ek1kl
       ekk(kk)=ek0kl*ek1kl
       eto=eto+ekk(kk)
       enddo
       write(6,*)'initial DNS Mans Wray spectrum'
       write(6,*)'     akkpp,qq,sig',akkpp,qq,sig
       write(6,*)'     aa,aa1,eto',aa,aa1,eto
                       endif
       if(imic.ne.4) then
       sig=10.
       sig2=sig**2
       isd=17
       call srand(isd)
c      call ranset(isd)
       qtima1=0.
       qtima2=0.
       qtima3=0.
c
c   velocity field in wave number space with assigned modul
c   to get the above spectra
c   the phase are random
c   attention to the random number generator, machine dependent
c
       do k=1,n3m
            akz=kz(k)
        do i=2,n1mh
            akx=kx(i)
         do j=1,n2m
            aky=ky(j)
            fk2=float(kx(i)*kx(i)+ky(j)*ky(j)+kz(k)*kz(k))
            fk=sqrt(fk2)
            kk=fk
            fkhor2=float(kx(i)*kx(i)+ky(j)*ky(j))
            fkhor=sqrt(fk2)
            ampk=sqrt(ekk(kk)/(4.*pi*fk2))
c           x1=rnd()
c           x2=rnd()
c           x3=rnd()
            x1=rand()
            x2=rand()
            x3=rand()
c           x1=ranf()
c           x2=ranf()
c           x3=ranf()
            phi=2.*pi*x1
            the1=2.*pi*x2
            the2=2.*pi*x3
            cexp1=cmplx(cos(the1),sin(the1))
            cexp2=cmplx(cos(the2),sin(the2))
            alp=ampk*cexp1*cos(phi)
            bet=ampk*cexp2*sin(phi)
            qtil1(i,j,k)=(alp*fk*akx+bet*akx*aky)/(fk*fkhor)
            qtil2(i,j,k)=(-alp*fk*aky+bet*akz*akx)/(fk*fkhor)
            qtil3(i,j,k)=-(bet*fkhor)/fk
           qtima1=max(abs(qtil1(i,j,k)),qtima1)
           qtima2=max(abs(qtil2(i,j,k)),qtima2)
           qtima3=max(abs(qtil3(i,j,k)),qtima3)
         end do
        end do
       end do
c
c    zero wave number
c
       i=1
            akx=kx(i)
       do k=1,n3m
            akz=kz(k)
         do j=2,n2m
            aky=ky(j)
            fk2=float(kx(i)*kx(i)+ky(j)*ky(j)+kz(k)*kz(k))
            fk=sqrt(fk2)
            kk=fk
            fkhor2=float(kx(i)*kx(i)+ky(j)*ky(j))
            fkhor=sqrt(fk2)
            ampk=sqrt(ekk(kk)/(4.*pi*fk2))
c           x1=rnd()
c           x2=rnd()
c           x3=rnd()
            x1=rand()
            x2=rand()
            x3=rand()
c           x1=ranf()
c           x2=ranf()
c           x3=ranf()
            phi=2.*pi*x1
            the1=2.*pi*x2
            the2=2.*pi*x3
            cexp1=cmplx(cos(the1),sin(the1))
            cexp2=cmplx(cos(the2),sin(the2))
            alp=ampk*cexp1*cos(phi)
            bet=ampk*cexp2*sin(phi)
            qtil1(i,j,k)=(alp*fk*akx+bet*akx*aky)/(fk*fkhor)
            qtil2(i,j,k)=(-alp*fk*aky+bet*akz*akx)/(fk*fkhor)
            qtil3(i,j,k)=-(bet*fkhor)/fk
           qtima1=max(abs(qtil1(i,j,k)),qtima1)
           qtima2=max(abs(qtil2(i,j,k)),qtima2)
           qtima3=max(abs(qtil3(i,j,k)),qtima3)
         end do
       end do
      write(6,*)'max qtil in wave numbers',qtima1,qtima2,qtima3
      open(57,file='spmanwrint.out',status='unknown')
      open(58,file='spmanwrinc.out',status='unknown')
                     endif

      if(imic.eq.4)  then
c
c   in the following part the initail field in wave
c   number space is read from a file produced by a pseudospectral
c   simulation. In here a field evaluated by the Rogallo code
c   by Wray is read.
c
      write(6,*)'  imic=',imic
      n1mh=n1m/2
      open(17,file='fieldinik.res',form='unformatted',status='unknown')
      read(17)time,n2l,n1l,n3l,ns
      write(6,*)time,n2l,n1l,n3l,ns
      write(6,*)' initila field read by wave number field generated'
     1         ,' by a pseudospectral codee at t=',time
      read(17)((((qtil1(i,j,k),qtil2(i,j,k),qtil3(i,j,k))
     1            ,i=1,n1mh),j=1,n2m),k=1,n3m) 
      close(17)
      write(6,*)' the file was read kkmax=',kkmax
      open(57,file='spalaninto.out',status='unknown')
      open(58,file='spalaninco.out',status='unknown')
                     endif
c
c   3D spectrum is calculated
c
       do kk=1,kkmax
       do n=1,3
       e(n,kk)=0.0
       enddo
       enddo
       do k=1,n3m
        do i=2,n1mh
         do j=1,n2m
         akk=sqrt(float(kx(i)*kx(i)+ky(j)*ky(j)+kz(k)*kz(k)))
         kk=nint(akk)
         urea=real(qtil1(i,j,k))
         uimm=aimag(qtil1(i,j,k))
         e(1,kk)=e(1,kk)+2.*(urea*urea+uimm*uimm)
         urea=real(qtil2(i,j,k))
         uimm=aimag(qtil2(i,j,k))
         e(2,kk)=e(2,kk)+2.*(urea*urea+uimm*uimm)
         urea=real(qtil3(i,j,k))
         uimm=aimag(qtil3(i,j,k))
         e(3,kk)=e(3,kk)+2.*(urea*urea+uimm*uimm)
         end do
        end do
       end do
       i=1
       do k=1,n3m
         do j=1,n2m
         akk=sqrt(float(kx(i)*kx(i)+ky(j)*ky(j)+kz(k)*kz(k)))
         kk=nint(akk)
         urea=real(qtil1(i,j,k))
         uimm=aimag(qtil1(i,j,k))
         e(1,kk)=e(1,kk)+(urea*urea+uimm*uimm)
         urea=real(qtil2(i,j,k))
         uimm=aimag(qtil2(i,j,k))
         e(2,kk)=e(2,kk)+(urea*urea+uimm*uimm)
         urea=real(qtil3(i,j,k))
         uimm=aimag(qtil3(i,j,k))
         e(3,kk)=e(3,kk)+(urea*urea+uimm*uimm)
         end do
       end do
       do k=1,kkmax
       write (58,*)k,e(1,k),e(2,k),e(3,k)
       write (57,*)k,e(1,k)+e(2,k)+e(3,k)
       end do
      close(57)
      close(58)

212    continue
c
c   in the fft3d routine the field is transformed from wave number to
c   physical space.
c
      call fft3d(q,1,qtil1,+1)
      call fft3d(q,2,qtil2,+1)
      call fft3d(q,3,qtil3,+1)
      return
      end
c***********************************************************   
      subroutine fft3d(q,ll,qtil,ifft)
c   
c   in this routine a field
c    in wave number space is transformed in 
c   physical space  ifft=+1
c   the viceversa for ifft=-1
c   Real FFT in the direction x1
c
      include 'param.f'
      dimension q(3,m1,m2,m3)
      complex qtil((m1-1)/2+1,m2,m3)
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/numond/kx(m1),ky(m2),kz(m3),kkmax
      common/speene/e(ndv,0:m3)
      common/fft002/ifxx2(13),trigxx2(2*(m2-1))
      common/fft001/ifxx1(13),trigxx1(3*(m1-1)/2+1)
      common/fft003/ifxx3(13),trigxx3(2*(m3-1))
      common/sparr/xa,wor,xa2,wor2,xr,work
      complex xa(m3-1,m2-1),wor(m3-1,m2-1)
      complex xa2(m2-1,m3-1),wor2(m2-1,m3-1)
      real xr(m1+1,m2-1),work(m1,m2-1)
c
      n1mh=n1m/2+1
c
      if(ifft.eq.-1) then
c
c  from  PHYSICAL to WAVE NUMBER space
c
      do  10 k=1,n3m
c
        do j=1,n2m
         xr(1,j)=q(ll,n1m,j,k)
         xr(n1+1,j)=q(ll,1,j,k)
         do i=1,n1m
          is=i+1
          xr(is,j)=q(ll,i,j,k)
         enddo
        enddo
c
c   real fft applied from  physical space to wave number
c   in x1 direction
c
        call fft99(xr,work,trigxx1,ifxx1,1,m1+1,n1m,n2m,ifft)
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
c   2-d  cfft applied (twice) from
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
        call cfft99(xa,wor,trigxx3,ifxx3,1,m3-1,n3m,n2m,ifft)
c
        do k=1,n3m
         do j=1,n2m
          xa2(j,k)=xa(k,j)/float(n3m)
         enddo
        enddo
c
        call cfft99(xa2,wor2,trigxx2,ifxx2,1,m2-1,n2m,n3m,ifft)
c
        do j=1,n2m
         do k=1,n3m
          qtil(i,j,k)=xa2(j,k)/float(n2m)
         enddo
        enddo
c
  20   continue
         else
c
c   inverse  fft  ifft+1
c
c  from  WAVE NUMBER to PHYSICAL space
c
c
      do 41 i=1,n1mh
        do j=1,n2m
         do k=1,n3m
          xa2(j,k)=qtil(i,j,k)
         enddo
        enddo
        call cfft99(xa2,wor2,trigxx2,ifxx2,1,m2-1,n2m,n3m,ifft)
        do k=1,n3m
         do j=1,n2m
          xa(k,j)=xa2(j,k)
         enddo
        enddo
c
        call cfft99(xa,wor,trigxx3,ifxx3,1,m3-1,n3m,n2m,ifft)
c
c
        do j=1,n2m
         do k=1,n3m
          qtil(i,j,k)=xa(k,j)
         enddo
        enddo
c
c
  41   continue
      do  40 k=1,n3m
c
        do i=1,n1mh
         ip=2*i   
         id=2*i-1
         do j=1,n2m
          xr(id,j)=real(qtil(i,j,k)) 
          xr(ip,j)=aimag(qtil(i,j,k)) 
         enddo
        enddo
c
c   real fft applied from  physical space to wave number 
c   in y direction
c
        call fft99(xr,work,trigxx1,ifxx1,1,m1+1,n1m,n2m,ifft)
        do j=1,n2m
         do i=1,n1m
          is=i+1
          q(ll,i,j,k)=xr(is,j)
         enddo
        enddo
c
c
 40   continue
             endif
       return
       end 
