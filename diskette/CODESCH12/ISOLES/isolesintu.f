c
c  ****************************** subrout initur **********************
c
c   initial conditions 
c   from the asigned spectrum a velocity field with
c   random phases is given
c   this is non-solenoidal then it is projected
c   in a solenoidal field by evaluating the divergence from this
c   dph and finally the solenoidal field
c
      subroutine initur(q,pr,qcap,dph)
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
      pi=4.*atan(1.)
c
c call mic to make initial conditions. Use qcap,dph and pr as work files
c
      call turbis(qcap,dph,pr,q)
        al=1.
      call divg(qcap,q,al)
c
c qcap now contains the divergence of the velocity field q
c
c  ********* calculation of the pressure dph by fft in two
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
c
c   the mean velocity is subctracted
c
      do 414 k=1,n3m
      do 414 j=1,n2m
      do 414 i=1,n1m
       q(1,i,j,k) =q(1,i,j,k)-vm1
       q(2,i,j,k) =q(2,i,j,k)-vm2
       q(3,i,j,k) =q(3,i,j,k)-vm3
414   continue
      vm1m=0.
      vm2m=0.
      vm3m=0.
      vl123=1./float(n1m*n2m*n3m)
      do 511 k=1,n3m
      do 511 j=1,n2m
      do 511 i=1,n1m
       vm1m=q(1,i,j,k)+vm1m
       vm2m=q(2,i,j,k)+vm2m
       vm3m=q(3,i,j,k)+vm3m
  511 continue
      vm1=vm1m*vl123
      vm2=vm2m*vl123
      vm3=vm3m*vl123
      print*,'vm1,vm2,vm3',vm1,vm2,vm3
      call divgck(q,qmax)
      print *,'divergenza massima=',qmax
      return
      end
*************************************************************************
c
c   here the initial spectra is chosen
c
      subroutine turbis(qtil1,qtil2,qtil3,q)
c make intitial conditions
      include 'param.f'
      dimension q(ndv,m1,m2,m3)
      complex qtil1(m1,(m2-1)/2+1,m3)
      complex qtil2(m1,(m2-1)/2+1,m3)
      complex qtil3(m1,(m2-1)/2+1,m3)
      common/rhsc/rhs(m1,m2,m3)
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/numond/kx(m1),ky(m2),kz(m3),kkmax
      common/rot/f0
      common/d123/alx1,alx2,alx3
      complex alp,bet,cexp1,cexp2,cexp3,cexp4
      dimension e(3,m1)
      dimension ee(3,0:m1)
      common/eispe/ekk(m1)
      common/ispec/imic
      common/spemw/akkpp,qq,sig
      common/teddy/akedd,eeto
      common/ledat/ics0,cvisc,ifiltr,ibox
      common/tstep/dt,beta,ren
      common/resca/iresca,tresca,rlamas
      common/eneas/engive
      common/kpmasp/kmax
c
      n2mh=n2m/2+1

      pi=4.*atan(1.)

      if(imic.eq.1)  then
c++++++++++COMPTE BELOLOT - Corrsin+++++++++++++++++++++
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
       eneco=0.
       disco=0.
       open(21,file='Compte.plo')
       do kk=1,32
       eneco=eneco+ekk(kk)
       write(21,*)kk,ekk(kk)
       disco=disco+ekk(kk)*kk**2
       enddo
       close(21)
       write(6,*)' Compte Bellot init. spectrum',
     1      cvisc,eneco,disco
       cvisco=cvisc
       reno=ren
       rlamco=reno*eneco*sqrt(20./(3.*disco))
       ren=reno*(rlamas/rlamco)
       write(6,*)' eneco, disco ',eneco,disco
       cvisc=1./ren
       rlam=ren*eneco*sqrt(20./(3.*disco))
       write(6,*)'  reno,ren',reno,ren,rlamco,rlam
                      endif
       if(imic.eq.0)  then
       open(21,file='speread.dat')
       do kk=1,500
       read(21,*,end=211)k,e(1,k),e(2,k),e(3,k)
       ekk(k)=e(1,k)+e(2,k)+e(3,k)
       enddo
211    continue
                      endif
       if(imic.eq.2)  then
c      
c  here we try  the initial spectrum of Mansour and Wray
c               DNS
       aa=0.
       do kk=1,kkmax
       ak=kk
c      rk=ak/akkpp
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
       if(imic.eq.3)  then
c      
c  here we try  the initial spectrum of Lesieur  and Rogallo
c               
       aa=0.
       do kk=1,kkmax
       ak=kk
       rk=ak/akkpp
       aint=ak**sig*exp(-sig*rk**2/2.)
       aa=aa+aint
       enddo
       aa1=qq/aa
       eto=0.
       do kk=1,kkmax
       ak=kk
       rk=ak/akkpp
       ek1kl=aa1*ak**sig*exp(-sig*rk**2/2.)
       ekk(kk)=ek1kl
       eto=eto+ekk(kk)
       enddo
       write(6,*)'   akpp=',akkpp,'   aa1=',aa1
       write(6,*)'   qq=',qq,'   eto=',eto
                    endif
       if(imic.eq.4)  then
c      
c  here we try  the initial spectrum of  Chasnov (PFA)
c               
       aa1=qq/akkpp*(2./pi)**0.5*256/35
       eto=0.
       do kk=1,kkmax
       ak=kk
       rk=ak/akkpp
       ek1kl=aa1*rk**8*exp(-2*rk**2)
       ekk(kk)=ek1kl
       eto=eto+ekk(kk)
       enddo
       write(6,*)'   akpp=',akkpp,'   aa1=',aa1
       write(6,*)'   qq=',qq,'   eto=',eto
                    endif
       akedd=akkpp
       eeto=eto
       isd=17
c
c   velocity field with random phase and magnitude
c   derived from the energy spectrum
c   Here the field in wave number space
c   Random number generator depends on the computer
c
       call srand(isd)
       do k=1,n3m
            akz=kz(k)
        do j=2,n2mh
            aky=ky(j)
         do i=1,n1m
            akx=kx(i)
            fk2=float(kx(i)*kx(i)+ky(j)*ky(j)+kz(k)*kz(k))
            fk=sqrt(fk2)
            kk=fk
            fkhor2=float(kx(i)*kx(i)+ky(j)*ky(j))
            fkhor=sqrt(fkhor2)
            ampk=sqrt(ekk(kk)/(2.*pi*fk2))
c           x1=rnd()
c           x2=rnd()
c           x3=rnd()
            x1=rand()
            x2=rand()
            x3=rand()
c           x1=rand()
c           x2=rand()
c           x3=rand()
            phi=2.*pi*x1
            the1=2.*pi*x2
            the2=2.*pi*x3
            cexp1=cmplx(cos(the1),sin(the1))
            cexp2=cmplx(cos(the2),sin(the2))
            alp=ampk*cexp1*cos(phi)
            bet=ampk*cexp2*sin(phi)
            qtil1(i,j,k)=(alp*fk*aky+bet*akx*akz)/(fk*fkhor)
            qtil2(i,j,k)=(-alp*fk*akx+bet*akz*aky)/(fk*fkhor)
            qtil3(i,j,k)=-(bet*fkhor)/fk
         end do
        end do
       end do
       j=1
            aky=ky(j)
       do k=1,n3m
            akz=kz(k)
         do i=2,n1m
            akx=kx(i)
            fk2=float(kx(i)*kx(i)+ky(j)*ky(j)+kz(k)*kz(k))
            fk=sqrt(fk2)
            kk=fk
            fkhor2=float(kx(i)*kx(i)+ky(j)*ky(j))
            fkhor=sqrt(fkhor2)
            ampk=sqrt(ekk(kk)/(2.*pi*fk2))
c           x1=rnd()
c           x2=rnd()
c           x3=rnd()
            x1=rand()
            x2=rand()
            x3=rand()
c           x1=rand()
c           x2=rand()
c           x3=rand()
            phi=2.*pi*x1
            the1=2.*pi*x2
            the2=2.*pi*x3
            cexp1=cmplx(cos(the1),sin(the1))
            cexp2=cmplx(cos(the2),sin(the2))
            alp=ampk*cexp1*cos(phi)
            bet=ampk*cexp2*sin(phi)
            qtil1(i,j,k)=(alp*fk*aky+bet*akx*akz)/(fk*fkhor)
            qtil2(i,j,k)=(-alp*fk*akx+bet*akz*aky)/(fk*fkhor)
            qtil3(i,j,k)=-(bet*fkhor)/fk
         end do
       end do
c
c   the spectrum is calculated to see how it differs
c   from the assigned one
c
      do 1 n=1,3
       do kk=1,kkmax
       ee(n,kk)=0.0
       enddo
       do k=1,n3m
        do j=2,n2mh
         do i=1,n1m
         akk=sqrt(float(kx(i)*kx(i)+ky(j)*ky(j)+kz(k)*kz(k)))
         kk=nint(akk)
         if(n.eq.1) then
         urea=real(qtil1(i,j,k))
         uimm=aimag(qtil1(i,j,k))
                    endif
         if(n.eq.2) then
         urea=real(qtil2(i,j,k))
         uimm=aimag(qtil2(i,j,k))
                    endif
         if(n.eq.3) then
         urea=real(qtil3(i,j,k))
         uimm=aimag(qtil3(i,j,k))
                    endif
         ee(n,kk)=ee(n,kk)+(urea*urea+uimm*uimm)
         end do
        end do
       end do
       j=1
       do k=1,n3m
         do i=1,n1m
         akk=sqrt(float(kx(i)*kx(i)+ky(j)*ky(j)+kz(k)*kz(k)))
         kk=nint(akk)
         if(n.eq.1) then
         urea=real(qtil1(i,j,k))
         uimm=aimag(qtil1(i,j,k))
                    endif
         if(n.eq.2) then
         urea=real(qtil2(i,j,k))
         uimm=aimag(qtil2(i,j,k))
                    endif
         if(n.eq.3) then
         urea=real(qtil3(i,j,k))
         uimm=aimag(qtil3(i,j,k))
                    endif
         ee(n,kk)=ee(n,kk)+(urea*urea+uimm*uimm)*0.5
         end do
       end do
  1    continue
       open(51,file='inispe.plo')
       eneco=0.
       do kk=1,kmax
       eneco=eneco+ee(1,kk)+ee(2,kk)+ee(3,kk)
       write(51,*)kk,ekk(kk), ee(1,kk)+ee(2,kk)+ee(3,kk)
       enddo
212    continue
       engive=eneco
       write(6,*)' engive=',kmax,engive
c
c    the field is transformed from the wave number to
c    the physical space
c
      call fft3d(rhs,qtil1,+1)
      do kc=1,n3m
      do jc=1,n2m
      do ic=1,n1m
      q(1,ic,jc,kc)=rhs(ic,jc,kc)
      enddo
      enddo
      enddo
      call fft3d(rhs,qtil2,+1)
      do kc=1,n3m
      do jc=1,n2m
      do ic=1,n1m
      q(2,ic,jc,kc)=rhs(ic,jc,kc)
      enddo
      enddo
      enddo
      call fft3d(rhs,qtil3,+1)
      do kc=1,n3m
      do jc=1,n2m
      do ic=1,n1m
      q(3,ic,jc,kc)=rhs(ic,jc,kc)
      enddo
      enddo
      enddo
      return
      end
c***********************************************************   
      subroutine fft3d(a,b,ifft)
c
c   routine to transform a quantity given in wave
c   number space to that in physical space  ifft=-1
c    VICEVERSA   with ifft=+1
c
c
      include 'param.f'
      dimension a(m1,m2,m3)
      complex b(m1,(m2-1)/2+1,m3)
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/fft001/ifxx1(13),trigxx1(2*(m1-1))
      common/fft002/ifxx2(13),trigxx2(3*(m2-1)/2+1)
      common/fft003/ifxx3(13),trigxx3(2*(m3-1))
      common/numond/kx(m1),ky(m2),kz(m3),kkmax
      common/speene/e(ndv,0:m3)
      common/spequa/xa,wor,xa2,wor2,xr,work
      complex xa(m3-1,m1-1),wor(m3-1,m1-1)
      complex xa2(m1-1,m3-1),wor2(m1-1,m3-1)
      real xr(m2+1,m1-1),work(m2,m1-1)
c
      n2mh=n2m/2+1
c
      if(ifft.eq.-1) then
c     do 1 n=1,3
      do  10 k=1,n3m
c
        do i=1,n1m
         xr(1,i)=a(i,n2m,k)
         xr(n2m+2,i)=a(i,1,k)
         do j=1,n2m
          js=j+1
          xr(js,i)=a(i,j,k)
         enddo
        enddo
c
c   real fft applied from  physical space to wave number 
c   in y direction
c
        call fft99(xr,work,trigxx2,ifxx2,1,m2+1,n2m,n1m,ifft)
c
        do j=1,n2mh
         jp=2*j   
         jd=2*j-1
         do i=1,n1m
          b(i,j,k)=cmplx(xr(jd,i),xr(jp,i))
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
          xa(k,i)=b(i,j,k)
         enddo
        enddo
c
        call cfft99(xa,wor,trigxx3,ifxx3,1,m3-1,n3m,n1m,ifft)
c
        do k=1,n3m
         do i=1,n1m
          xa2(i,k)=xa(k,i)/float(n3m)
         enddo
        enddo
c
        call cfft99(xa2,wor2,trigxx1,ifxx1,1,m1-1,n1m,n3m,ifft)
c
        do i=1,n1m
         do k=1,n3m
          b(i,j,k)=xa2(i,k)/float(n1m)
         enddo
        enddo
c
  20   continue
c
c     b is the velocity component in Fourier space.
c
      
     

      else

      do 41 j=1,n2mh
c
        do i=1,n1m
         do k=1,n3m
          xa2(i,k)=b(i,j,k)
         enddo
        enddo
        call cfft99(xa2,wor2,trigxx1,ifxx1,1,m1-1,n1m,n3m,ifft)
        do k=1,n3m
         do i=1,n1m
          xa(k,i)=xa2(i,k)
         enddo
        enddo
c
        call cfft99(xa,wor,trigxx3,ifxx3,1,m3-1,n3m,n1m,ifft)
c
c
        do i=1,n1m
         do k=1,n3m
          b(i,j,k)=xa(k,i)
         enddo
        enddo
c
c
  41   continue
      do  40 k=1,n3m
c
        do j=1,n2mh
         jp=2*j   
         jd=2*j-1
         do i=1,n1m
          xr(jd,i)=real(b(i,j,k)) 
          xr(jp,i)=aimag(b(i,j,k)) 
         enddo
        enddo
c
c   real fft applied from  physical space to wave number 
c   in y direction
c
        call fft99(xr,work,trigxx2,ifxx2,1,m2+1,n2m,n1m,ifft)
        do i=1,n1m
         do j=1,n2m
          js=j+1
          a(i,j,k)=xr(js,i)
         enddo
        enddo
 40   continue
       endif
       return
       end 
*************************************************************************
      subroutine rescal(q,qtil1,qtil2,qtil3)
c
c   here there is the eventuality to rescal a field
c   to have the desired  R-lambda
c
      include 'param.f'
      dimension q(ndv,m1,m2,m3)
      complex qtil1(m1,(m2-1)/2+1,m3)
      complex qtil2(m1,(m2-1)/2+1,m3)
      complex qtil3(m1,(m2-1)/2+1,m3)
      common/rhsc/rhs(m1,m2,m3)
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/numond/kx(m1),ky(m2),kz(m3),kkmax
      complex alp,bet,cexp1,cexp2,cexp3,cexp4
      dimension ei(3,0:m1)
      dimension ee(3,0:m1)
      common/eispe/ekk(m1)
      common/ispec/imic
      common/spemw/akkpp,qq,sig
      common/teddy/akedd,eeto
      common/ledat/ics0,cvisc,ifiltr,ibox
      common/tstep/dt,beta,ren
      common/kpmasp/kmax
      common/resca/iresca,tresca,rlamas
      common/eneas/engive
c
c   evaluates the velocities in wave number space
c
      n2mh=n2m/2+1
            call qphywa(q,qtil1,1)
         call enekco(qtil1,1,ei)
            call qphywa(q,qtil2,2)
         call enekco(qtil2,2,ei)
            call qphywa(q,qtil3,3)
         call enekco(qtil3,3,ei)
       eeto=0.
       disk=0.
       do k=1,kmax
       eeto=eeto+ei(1,k)+ei(2,k)+ei(3,k)
       disk=disk+(ei(1,k)+ei(2,k)+ei(3,k))*k**2
       enddo
       rlamfi=ren*eeto*sqrt(20./(3.*disk))
       write(6,*)' eeto,disk,rlamas,rlamfi',eeto,disk,rlamas,rlamfi
       cvrerl=rlamas/rlamfi
       cvresc=sqrt(engive/eeto)
       write(6,*)' rescale the velocities rlam',cvrerl
     1       ,' energy =',cvresc
       do k=1,n3m
        do j=1,n2m
         do i=1,n1m
            q(1,i,j,k)=q(1,i,j,k)*cvresc
            q(2,i,j,k)=q(2,i,j,k)*cvresc
            q(3,i,j,k)=q(3,i,j,k)*cvresc
         end do
        end do
       end do
      call divgck(q,qmax)
      print *,'max diverg of rescaled field =',qmax
            call qphywa(q,qtil1,1)
         call enekco(qtil1,1,ei)
            call qphywa(q,qtil2,2)
         call enekco(qtil2,2,ei)
            call qphywa(q,qtil3,3)
         call enekco(qtil3,3,ei)
       eeto=0.
       disk=0.
       do k=1,kmax
       eeto=eeto+ei(1,k)+ei(2,k)+ei(3,k)
       disk=disk+(ei(1,k)+ei(2,k)+ei(3,k))*k**2
       enddo
       rlamne=ren*eeto*sqrt(20./(3.*disk))
       write(6,*)' eeto,disk,rlamne',eeto,disk,rlamne
       return
       end 

