c
c  ****************************** subrout updvp  **********************
c
c  this subroutine calculate the solenoidal vel field
c       q(n+1)=qhat-grad(dph)*dt ,  pr=dph
c  third order runge-kutta is used.
c
      subroutine updvp(dph,q,al)
c
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      common/metria/caj(m2),cac(m2)
      dimension dph(m1,m2,m3),q(ndv,m1,m2,m3)
      common/tstep/dt,beta,ren
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/qwallo/q1s(m1,m3),q2s(m1,m3),q3s(m1,m3)
      common/qwalup/q1n(m1,m3),q2n(m1,m3),q3n(m1,m3)

c
c  ***********  compute the q1 velocity component
c               dfx11=component 1 of grad(dph)
      do 1 kc=1,n3m
      do 1 jc=1,n2m
      do 1 ic=1,n1m
      im=imv(ic)
      dfx11=(dph(ic,jc,kc)-dph(im,jc,kc))*dx1
      q(1,ic,jc,kc)=q(1,ic,jc,kc)-dfx11*dt*al
    1 continue
c
c  ***********  compute the q2 velocity component
c               dfx22=component 2 of grad(dph)
      do 2 kc=1,n3m
      do 2 jc=2,n2m
      sucac=1./cac(jc)
      do 2 ic=1,n1m
      jm=jc-1
      dfx22=(dph(ic,jc,kc)-dph(ic,jm,kc))*dx2
      q(2,ic,jc,kc)=q(2,ic,jc,kc)-dfx22*dt*al*sucac
    2 continue
      do  kc=1,n3m
      do  ic=1,n1m
       q(2,ic,1,kc)=q2s(ic,kc)
       q(2,ic,n2,kc)=q2n(ic,kc)
      end do
      end do
c
c  ***********  compute the q3 velocity component
c               q3 is the cartesian component
c               dfx33=component 3 of grad(dph)
      do 5 kc=1,n3m
      km=kmv(kc)
      do 5 jc=1,n2m
      do 5 ic=1,n1m
      dfx33=(dph(ic,jc,kc)-dph(ic,jc,km))*dx3
      q(3,ic,jc,kc)=q(3,ic,jc,kc)-dfx33*al*dt
    5 continue
      return
      end

c
c  ***********  compute the coefficients for tridiagonal 
c               in the 2 direction
c
      subroutine coeinv
c  
      include 'param.f'
      common/amcpj1/amj1(m2),acj1(m2),apj1(m2)
      common/amcpj2/amj2(m2),acj2(m2),apj2(m2)
      common/amcpj3/amj3(m2),acj3(m2),apj3(m2)
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/metria/caj(m2),cac(m2)
      common/islwal/islv1s,islv1n,islv3s,islv3n
c
c   set up the coefficients apj1, acj1, amj1 at interior points
c
      do jc=2,n2m-1
       jp=jc+1
       jm=jc-1
       ucaj=1./caj(jc)
       apj1(jc)=1./cac(jp)*ucaj
       acj1(jc)=(1./cac(jp)+1./cac(jc))*ucaj
       amj1(jc)=1./cac(jc)*ucaj
       apj3(jc)=apj1(jc)
       acj3(jc)=acj1(jc)
       amj3(jc)=amj1(jc)
      end do
c
c   set up the cefficients apj1, acj1, amj1 at the boundaries
      if(islv1s.eq.0) then
c
c  jc=1  spanwise velocity
c   in the case of q1  assigned
c
      ucaj=4./(2.*cac(2)+caj(1))
      apj1(1)=1./cac(2)*ucaj
      acj1(1)=(1./cac(2)+2./cac(1))*ucaj
      amj1(1)=0.
                       else
c
c   in the case of q1 shear free  
c
      ucaj=1./caj(1)
      apj1(1)=1./cac(2)*ucaj
      acj1(1)=1./cac(2)*ucaj
      amj1(1)=0.
                       endif
      if(islv3s.eq.0) then
c
c  jc=1  streamwise velocity
c   in the case of  q3 assigned
c
      ucaj=4./(2.*cac(2)+caj(1))
      apj3(1)=1./cac(2)*ucaj
      acj3(1)=(1./cac(2)+2./cac(1))*ucaj
      amj3(1)=0.
                       else
c
c   in the case of q3  shear free
c
      ucaj=1./caj(1)
      apj3(1)=1./cac(2)*ucaj
      acj3(1)=1./cac(2)*ucaj
      amj3(1)=0.
                       endif
c
c  jc=n2m streamwise velocity
c   in the case of  q1 assigned
c
      if(islv1n.eq.0) then
      ucaj=4./(2.*cac(n2m)+caj(n2m))
      apj1(n2m)=0.
      acj1(n2m)=(2./cac(n2)+1./cac(n2m))*ucaj
      amj1(n2m)=1./cac(n2m)*ucaj
                       else
c
c   in the case of q1  shear free
c
      ucaj=1./caj(n2m)
      apj1(n2m)=0.
      acj1(n2m)=1./cac(n2m)*ucaj
      amj1(n2m)=1./cac(n2m)*ucaj
                       endif
c
c   in the case of  q3 assigned
c
      if(islv3n.eq.0) then
      ucaj=4./(2.*cac(n2m)+caj(n2m))
      apj3(n2m)=0.
      acj3(n2m)=(2./cac(n2)+1./cac(n2m))*ucaj
      amj3(n2m)=1./cac(n2m)*ucaj
                       else
c
c   in the case of q3  shear free
c
      ucaj=1./caj(n2m)
      apj3(n2m)=0.
      acj3(n2m)=1./cac(n2m)*ucaj
      amj3(n2m)=1./cac(n2m)*ucaj
                       endif
c
c   set up the coefficients apj2, acj2, amj2 at interior points
c
      do jc=2,n2m
      jm=jc-1
      ucac=1./cac(jc)
       apj2(jc)=1./caj(jc)*ucac
       acj2(jc)=(1./caj(jc)+1./caj(jm))*ucac
       amj2(jc)=1./caj(jm)*ucac
      end do
c
c   the coefficients at jc=1  and n2 are evaluated in the
c   invtr_i   routines
c
      open(17,file='coeff')
      do j=2,n2m
      write(17,*) j,apj1(j),acj1(j),amj1(j),
     1              apj2(j),acj2(j),amj2(j)
      end do
      close(17)
      return
      end
c
c  **************  subrout tschem ***************************************
c  time advancement 
c
      subroutine tschem(q,pr,ru,qcap,dph,time)
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      dimension q(ndv,m1,m2,m3),qcap(m1,m2,m3)
     1         ,dph(m1,m2,m3),pr(m1,m2,m3)
      dimension ru(ndv,m1,m2,m3)
      common/tstep/dt,beta,ren
      common/tscoe/gam(3),rom(3),nsst
      common/movwal/tosc,uosc
      common/slotfl/flowq2,tau2
      common/slotii/tim0sl
c
c     open(61,file='ck.out')
      rewind(61)
      tloc=time
      do 2000 ns=1,nsst
c
c   time integration implicit viscous 
c
c  if nsst=1 adams bashford if nsst=3 runge-kutta
c
      al=(gam(ns)+rom(ns))
      ga=gam(ns)
      ro=rom(ns)
c
c    wall velocity b.conditions
c
      tloc=tloc+al*dt
      tslo=(tloc-tim0sl)/tau2
      ft=ftir(tslo)
      call boucdq(time,ft)
c
c  *****  computation of nonlinear terms
c
c 
      call hdnl1(q)
c 
      call hdnl2(q,dph)
c
      call hdnl3(q,qcap)
c
c  *****  solve the dqhat=qhat-q(n) momentum equations
c  
      call invtr1(q,al,ga,ro,pr,ru)
c
c
c 
      call invtr2(q,dph,al,ga,ro,pr,ru)
c
c
c 
      call invtr3(q,qcap,al,ga,ro,pr,ru)
c
      do k=1,n3m
      do j=1,n2
      do i=1,n1m
      q(2,i,j,k)=dph(i,j,k)+q(2,i,j,k)
      enddo
      enddo
      enddo
      do k=1,n3m
      do j=1,n2m
      do i=1,n1m
      q(3,i,j,k)=qcap(i,j,k)+q(3,i,j,k)
      enddo
      enddo
      enddo
c
c  ********* calculation of divg(dqhat)
c
c 
      call divg(qcap,q,al)
c
c  ********* calculation of the pressure dph by fft in two
c            directions and tridiag in vertical
c 
      call phcalc(qcap,dph)
c
c  ********* calculation of solenoidal vel field
c
c 
      call updvp(dph,q,al)
c
c  ********* calculation of pressure field
c
c 
      call prcalc(pr,dph,al)
c     call divgck(q,qmax)
c  
c
c
 2000 continue
c
      return
      end
c  ****************************** subrout invtr1  **********************
c
c   this subroutine performs the inversion of the q1 momentum equation
c   by a factored implicit scheme, the derivatives 11,22,33 of q1
c   are treated implicitly
      subroutine invtr1(q,al,ga,ro,pr,ru)
c
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      common/amcpj1/amj1(m2),acj1(m2),apj1(m2)
      common/trici/ami(m2,m1),aci(m2,m1),api(m2,m1),
     1           fi(m2,m1),fei(m2,m1),qq(m2,m1),ss(m2,m1)
      common/tvbkj/amj(m3,m2),acj(m3,m2),apj(m3,m2),fj(m3,m2)
      common/tvpjk/amk(m2,m3),ack(m2,m3),apk(m2,m3),fk(m2,m3)
     1         ,qek(m2,m3),qk(m2,m3),sk(m2,m3)
      common/metria/caj(m2),cac(m2)
      dimension q(ndv,m1,m2,m3)
      common/rhsc/rhs(m1,m2,m3)
      dimension pr(m1,m2,m3)
      dimension ru(ndv,m1,m2,m3)
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/tstep/dt,beta,ren
      common/qwallo/q1s(m1,m3),q2s(m1,m3),q3s(m1,m3)
      common/qwalup/q1n(m1,m3),q2n(m1,m3),q3n(m1,m3)
      common/dqwalo/dq1s(m1,m3),dq2s(m1,m3),dq3s(m1,m3)
      common/dqwalu/dq1n(m1,m3),dq2n(m1,m3),dq3n(m1,m3)
      common/islwal/islv1s,islv1n,islv3s,islv3n

c
c  ********* compute the rhs of the factored equation
c  the rhs include h(n), h(n-1)=ru, grad(pr) component
c  and the 11, 22,3 second derivat. of q1(n)
c  everything at i,j+1/2,k+1/2
c  dq=qhat-q(n)
c
      alre=al/ren
      n2mm=n2m-1
      do kc=1,n3m
      km=kmv(kc)
      kp=kpv(kc)
      do jc=1,n2m
      do ic=1,n1m
      ip=ipv(ic)
      im=imv(ic)
c
c   11 second deriv. of q1(n)
c
      d11q1=(q(1,ip,jc,kc)-2.*q(1,ic,jc,kc)+q(1,im,jc,kc))*dx1q
c
c   add the 33 derivat.
c
      d33q1=(q(1,ic,jc,kp)-2.*q(1,ic,jc,kc)+q(1,ic,jc,km))*dx3q
      dcq1=d11q1+d33q1
c
c   grad(pr) along 1
c
      dpx11=(pr(ic,jc,kc)-pr(im,jc,kc))*dx1
      gradp=dpx11*al
      rhsccc=(ga*rhs(ic,jc,kc)+ro*ru(1,ic,jc,kc)-gradp
     1             +alre*dcq1)*dt
      ru(1,ic,jc,kc)=rhs(ic,jc,kc) 
      rhs(ic,jc,kc)=rhsccc
      enddo
      enddo
      enddo
c
c   add the second derivative d22q1 in the inner points
c
      do kc=1,n3m
      do jc=2,n2mm
      jmm=jmv(jc)
      jpp=jpv(jc)
      do ic=1,n1m
c
c   22 second deriv. of q1(n)
c
      d22q1=(apj1(jc)*q(1,ic,jpp,kc)
     1      -acj1(jc)*q(1,ic,jc,kc)
     1      +amj1(jc)*q(1,ic,jmm,kc))*dx2q
      rhs(ic,jc,kc)=rhs(ic,jc,kc)+alre*d22q1*dt
      enddo
      enddo
      enddo
c
c   add the second derivative d22q1 at the lower wall 
c
      jc=1
      jp=jc+1
      do kc=1,n3m
      do ic=1,n1m
c
c   22 second deriv. of q1(n) O(dx2^2)
c
      am=4./(2.*cac(2)+caj(1))*2./cac(jc)
      d22q1=(apj1(jc)*q(1,ic,jp,kc)
     1      -acj1(jc)*q(1,ic,jc,kc)
     1      +am*q1s(ic,kc)*(1-islv1s) )*dx2q
      rhs(ic,jc,kc)=rhs(ic,jc,kc)+alre*d22q1*dt
      enddo
      enddo
c
c   add the second derivative d22q1 aqt the upper wall 
c
      jc=n2m
      jm=jc-1
      do kc=1,n3m
      do ic=1,n1m
c
c   22 second deriv. of q1(n) O(dx2^2)
c
      ap=4./(2.*cac(n2m)+caj(n2m))*2./cac(n2)
      d22q1=(amj1(jc)*q(1,ic,jm,kc)
     1      -acj1(jc)*q(1,ic,jc,kc)
     1      +ap*q1n(ic,kc)*(1-islv1n) )*dx2q
      rhs(ic,jc,kc)=rhs(ic,jc,kc)+alre*d22q1*dt
      enddo
      enddo
c
c  ********* compute dq1*  sweeping in the x1 direction
c            periodic
      betadx=dx1q*al*beta
      do 1 kc=1,n3m
      do 9 jc=1,n2m
      do 9 ic=1,n1m
      api(jc,ic)=-betadx
      aci(jc,ic)=1.+betadx*2.
      ami(jc,ic)=-betadx
      fi(jc,ic)=rhs(ic,jc,kc)
    9 continue
c
c   periodic tridiagonal solver
c
      call tripvi(1,n1m,1,n2m)
      do 3 jc=1,n2m
      do 3 ic=1,n1m
      rhs(ic,jc,kc)=fi(jc,ic)
    3 continue
    1 continue
c
c  ************ compute  from dq** sweeping along the x3 direction
c               periodic
c
      betadz=dx3q*al*beta
      do 6 ic=1,n1m
      do 8 kc=1,n3m
      do 8 jc=1,n2m
      apk(jc,kc)=-betadz
      ack(jc,kc)=1.+betadz*2.
      amk(jc,kc)=-betadz
      fk(jc,kc)=rhs(ic,jc,kc)
    8 continue
c
c   periodic tridiagonal solver
c
      call trvpjk(1,n3m,1,n2m)
      do 5 kc=1,n3m
      do 5 jc=1,n2m
      rhs(ic,jc,kc)=fk(jc,kc)
    5 continue
    6 continue
c
c  ************ compute dq1  sweeping along the x2 direction
c
      aldt=dx2q*al*beta
      ucaj=4./(2.*cac(2)+caj(1))
      am=ucaj/cac(1)*2.
      ucaj=4./(2.*cac(n2m)+caj(n2m))
      ap=ucaj/cac(n2)*2.
      do 110 ic=1,n1m
      do 120 jc=1,n2m
      do 120 kc=1,n3m
      apj(kc,jc)=-aldt*apj1(jc)
      acj(kc,jc)=1.+aldt*acj1(jc)
      amj(kc,jc)=-aldt*amj1(jc)
  120 continue
      do 135 kc=1,n3m
      fj(kc,1)=rhs(ic,1,kc)+aldt*am*dq1s(ic,kc)*(1-islv1s)
      fj(kc,n2m)=rhs(ic,n2m,kc)+aldt*ap*dq1n(ic,kc)*(1-islv1n)
      do 136 jc=2,n2mm
      fj(kc,jc)=rhs(ic,jc,kc)
  136 continue
  135 continue
c
c   tridiagonal solver
c
      call trvbkj(n2m,1,n3m)
c
      do 30 kc=1,n3m
      do 30 jc=1,n2m
      q(1,ic,jc,kc)=fj(kc,jc)+q(1,ic,jc,kc)
   30 continue
  110 continue
      return
      end
c
c  ****************************** subrout invtr2  **********************
c   this subroutine performs the inversion of the q2 momentum equation
c   by a factored implicit scheme, the derivatives 11,22,33 of q2
c   are treated implicitly
c
      subroutine invtr2(q,dq,al,ga,ro,pr,ru)
c
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      common/trici/ami(m2,m1),aci(m2,m1),api(m2,m1),
     1           fi(m2,m1),fei(m2,m1),qq(m2,m1),ss(m2,m1)
      common/tvbkj/amj(m3,m2),acj(m3,m2),apj(m3,m2),fj(m3,m2)
      common/tvpjk/amk(m2,m3),ack(m2,m3),apk(m2,m3),fk(m2,m3)
     1         ,qek(m2,m3),qk(m2,m3),sk(m2,m3)
      common/amcpj2/amj2(m2),acj2(m2),apj2(m2)
      common/metria/caj(m2),cac(m2)
      dimension ru(ndv,m1,m2,m3)
      common/rhsc/rhs(m1,m2,m3)
      dimension q(ndv,m1,m2,m3),pr(m1,m2,m3)
      dimension dq(m1,m2,m3)
      common/tstep/dt,beta,ren
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/qwallo/q1s(m1,m3),q2s(m1,m3),q3s(m1,m3)
      common/qwalup/q1n(m1,m3),q2n(m1,m3),q3n(m1,m3)
      common/dqwalo/dq1s(m1,m3),dq2s(m1,m3),dq3s(m1,m3)
      common/dqwalu/dq1n(m1,m3),dq2n(m1,m3),dq3n(m1,m3)
      common/islwal/islv1s,islv1n,islv3s,islv3n

c
c  ********* compute the rhs of the factored equation
c  the rhs include h(n), h(n-1)=ru, grad(pr) component
c  and the 11, 22,3 second derivat. of q2(n)
c  everything at i+1/2,j,k+1/2
c
      alre=al/ren
      do 16 kc=1,n3m
      km=kmv(kc)
      kp=kpv(kc)
      do 16 jc=2,n2m
      jm=jc-1
      jp=jc+1
      sucac=1./cac(jc)
      do 16 ic=1,n1m
      im=imv(ic)
      ip=ipv(ic)
c
c   11 second derivative of q2
c
      d11q2=(q(2,ip,jc,kc)-2.*q(2,ic,jc,kc)+q(2,im,jc,kc))*dx1q
c
c   22 second derivative of q2
c
      d22q2=( apj2(jc)*q(2,ic,jp,kc) 
     1       -acj2(jc)*q(2,ic,jc,kc)+
     1        amj2(jc)*q(2,ic,jm,kc))*dx2q
c
c   add 33 second derivative of q2
c
      d33q2=(q(2,ic,jc,kp)-2.*q(2,ic,jc,kc)+q(2,ic,jc,km))*dx3q
      dcq2=d11q2+d33q2+d22q2
c
c   component of grap(pr) along 2 direction
c
      dpx22=(pr(ic,jc,kc)-pr(ic,jm,kc))*dx2
      gradp=dpx22*al*sucac
      rhs(ic,jc,kc)=(ga*dq(ic,jc,kc)+ro*ru(2,ic,jc,kc)-gradp
     1             +alre*dcq2)*dt
      ru(2,ic,jc,kc)=dq(ic,jc,kc) 
   16 continue
c
c  ************ compute dq2** sweeping along the x1 direction
c               periodic
c
      betadx=dx1q*al*beta
      do 10 kc=1,n3m
      do 21 jc=2,n2m
      do 21 ic=1,n1m
      api(jc,ic)=-betadx
      aci(jc,ic)=1.+betadx*2.
      ami(jc,ic)=-betadx
      fi(jc,ic)=rhs(ic,jc,kc)
   21 continue
c
c  periodic tridiagonal solver
c
      call tripvi(1,n1m,2,n2m)
      do 30 jc=2,n2m
      do 30 ic=1,n1m
      rhs(ic,jc,kc)=fi(jc,ic)
   30 continue
   10 continue
c
      betadx=dx3q*al*beta
      do 5 ic=1,n1m
      do 8 kc=1,n3m
      do 8 jc=2,n2m
      apk(jc,kc)=-betadx
      ack(jc,kc)=1.+betadx*2.
      amk(jc,kc)=-betadx
      fk(jc,kc)=rhs(ic,jc,kc)
    8 continue
c
c  periodic tridiagonal solver
c
      call trvpjk(1,n3m,2,n2m)
      do 7 kc=1,n3m
      do 7 jc=2,n2m
      rhs(ic,jc,kc)=fk(jc,kc)
    7 continue
    5 continue


c
c  ********* compute the dq2* sweeping in the x2 direction
c            wall boundaries direction
c
      aldt=dx2q*al*beta
      do 1 ic=1,n1m
      do 11 kc=1,n3m
      amj(kc,1)=0.
      apj(kc,1)=0.
      acj(kc,1)=1.
      amj(kc,n2)=0.
      apj(kc,n2)=0.
      acj(kc,n2)=1.
      fj(kc,1)=dq2s(ic,kc)
      fj(kc,n2)=dq2n(ic,kc)
      do 11 jc=2,n2m
      apj(kc,jc)=-aldt*apj2(jc)
      acj(kc,jc)=1.+aldt*acj2(jc)
      amj(kc,jc)=-aldt*amj2(jc)
      fj(kc,jc)=rhs(ic,jc,kc)
   11 continue
c
c   tridiagonal solver
c
      call trvbkj(n2,1,n3m)
c
      do 3 kc=1,n3m
      do 3 jc=1,n2
      dq(ic,jc,kc)=fj(kc,jc)
    3 continue
    1 continue
c
      do 9 kc=1,n3m
      do 9 ic=1,n1m
      dq(ic,1,kc)=dq2s(ic,kc)
      dq(ic,n2,kc)=dq2n(ic,kc)
    9 continue
      return
      end
c
c  ****************************** subrout invtr3  **********************
c   this subroutine performs the inversion of the q3 momentum equation
c   by a factored implicit scheme, the derivatives 11,22,33 of q3
c   are treated implicitly
c
      subroutine invtr3(q,dq,al,ga,ro,pr,ru)
c

      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      common/amcpj3/amj3(m2),acj3(m2),apj3(m2)
      common/trici/ami(m2,m1),aci(m2,m1),api(m2,m1),
     1           fi(m2,m1),fei(m2,m1),qq(m2,m1),ss(m2,m1)
      common/tvbkj/amj(m3,m2),acj(m3,m2),apj(m3,m2),fj(m3,m2)
      common/tvpjk/amk(m2,m3),ack(m2,m3),apk(m2,m3),fk(m2,m3)
     1         ,qek(m2,m3),qk(m2,m3),sk(m2,m3)
      dimension ru(ndv,m1,m2,m3)
      common/rhsc/rhs(m1,m2,m3)
      dimension q(ndv,m1,m2,m3),pr(m1,m2,m3)
      dimension dq(m1,m2,m3)
      common/metria/caj(m2),cac(m2)
      common/tstep/dt,beta,ren
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/rhs3p/dp3ns
      common/neno/dpnew
      common/s3t/s3tot
c
      common/qwallo/q1s(m1,m3),q2s(m1,m3),q3s(m1,m3)
      common/qwalup/q1n(m1,m3),q2n(m1,m3),q3n(m1,m3)
      common/dqwalo/dq1s(m1,m3),dq2s(m1,m3),dq3s(m1,m3)
      common/dqwalu/dq1n(m1,m3),dq2n(m1,m3),dq3n(m1,m3)
      common/islwal/islv1s,islv1n,islv3s,islv3n
c
c
c  ********* compute the rhs of the factored equation
c  the rhs include h(n), h(n-1)=ru, grad(pr) component
c  and the 11, 22,3 second derivat. of q3(n)
c  everything at i,j+1/2,k=1/2
c
      alre=al/ren
      s3tot=0.
      n2mm=n2m-1
      do kc=1,n3m
      km=kmv(kc)
      kp=kpv(kc)
      do jc=1,n2m
      do ic=1,n1m
      im=imv(ic)
      ip=ipv(ic)
c
c   11 second derivatives of q3
c
      dq31=(q(3,ip,jc,kc)-2.*q(3,ic,jc,kc)+q(3,im,jc,kc))*dx1q
c
c   33 second derivatives of q3
c
      dq33=(q(3,ic,jc,kp)-2.*q(3,ic,jc,kc)+q(3,ic,jc,km))*dx3q
      dcq3=dq31+dq33
      rhs(ic,jc,kc)=(ga*dq(ic,jc,kc)+ro*ru(3,ic,jc,kc)+
     1               alre*dcq3)*dt
      s3tot=s3tot+dcq3*caj(jc)/ren
      ru(3,ic,jc,kc)=dq(ic,jc,kc) 
      enddo
      enddo
      enddo
c
c  add 22 second derivatives of q3 inner field
c
      do kc=1,n3m
      do jc=2,n2mm
      jmm=jmv(jc)
      jpp=jpv(jc)
      do ic=1,n1m
c
c   22 second derivatives of q3
c
      dq32=(apj3(jc)*q(3,ic,jpp,kc)
     1     -acj3(jc)*q(3,ic,jc,kc)
     1     +amj3(jc)*q(3,ic,jmm,kc))*dx2q
      rhs(ic,jc,kc)=rhs(ic,jc,kc)+alre*dq32*dt
      s3tot=s3tot+dq32*caj(jc)/ren
      enddo
      enddo
      enddo
c
c  add 22 second derivatives of q3 lower wall   
c
      jc=1
      jp=jc+1
      do kc=1,n3m
      do ic=1,n1m
c
c   22 second deriv. of q3(n) O(dx2^2)
c
      am=4./(2.*cac(2)+caj(1))*2./cac(jc)
      dq32=(apj3(jc)*q(3,ic,jp,kc)
     1     -acj3(jc)*q(3,ic,jc,kc)
     1     +am*q3s(ic,kc)*(1-islv3s) )*dx2q
      s3tot=s3tot+dq32*caj(jc)/ren
      rhs(ic,jc,kc)=rhs(ic,jc,kc)+alre*dq32*dt
      enddo
      enddo
c
c   add the second derivative d22q3 aqt the upper wall
c
      jc=n2m
      jm=jc-1
      do kc=1,n3m
      do ic=1,n1m
c
c   22 second deriv. of q3(n) O(dx2^2)
c
      ap=4./(2.*cac(n2m)+caj(n2m))*2./cac(n2)
      dq32=(amj3(jc)*q(3,ic,jm,kc)
     1     -acj3(jc)*q(3,ic,jc,kc)
     1     +ap*q3n(ic,kc)*(1-islv3n) )*dx2q
      s3tot=s3tot+dq32*caj(jc)/ren
      rhs(ic,jc,kc)=rhs(ic,jc,kc)+alre*dq32*dt
      enddo
      enddo
c
c   in dp3ns there is the mean pressure gradient to keep
c   constant mass
c
      dp3ns=s3tot/(2.*n1m*n2m*n3m)
c
      do 38 kc=1,n3m
      km=kmv(kc)
      do 38 jc=1,n2m
      do 38 ic=1,n1m
c
c  component of grad(pr) along x3 direction
c
      dpx33=(pr(ic,jc,kc)-pr(ic,jc,km))*dx3
      gradp=(dpx33+dp3ns)*al
      rhs(ic,jc,kc)=rhs(ic,jc,kc)-gradp*dt
   38 continue
c
c  ********* compute the dq3* sweeping in the x3 direction
c
      betadz=dx3q*al*beta
      do 1 ic=1,n1m
      do 15 kc=1,n3m
      do 15 jc=1,n2m
      apk(jc,kc)=-betadz
      ack(jc,kc)=1.+betadz*2.
      amk(jc,kc)=-betadz
      fk(jc,kc)=rhs(ic,jc,kc)
   15 continue
c
c    periodic tridiagonal solver
c
      call trvpjk(1,n3m,1,n2m)
c
      do 3 kc=1,n3m
      do 3 jc=1,n2m
      rhs(ic,jc,kc)=fk(jc,kc)
    3 continue
    1 continue
c
c  ************ compute dq3** sweeping along the x1 direction
c
      betadx=dx1q*al*beta
      do 4 kc=1,n3m
      do 14 jc=1,n2m
      do 14 ic=1,n1m
      api(jc,ic)=-betadx
      aci(jc,ic)=1.+betadx*2.
      ami(jc,ic)=-betadx
      fi(jc,ic)=rhs(ic,jc,kc)
   14 continue
c
c    periodic tridiagonal solver
c
      call tripvi(1,n1m,1,n2m)
c
      do 6 jc=1,n2m
      do 6 ic=1,n1m
      rhs(ic,jc,kc)=fi(jc,ic)
    6 continue
    4 continue
c
c  ************ compute dq3 sweeping along the x2 direction
c
      aldt=dx2q*al*beta
      ucaj=4./(2.*cac(2)+caj(1))
      am=ucaj/cac(1)*2.
      ucaj=4./(2.*cac(n2m)+caj(n2m))
      ap=ucaj/cac(n2)*2.
      do 10 ic=1,n1m
      do 11 kc=1,n3m
      do 11 jc=1,n2m
      apj(kc,jc)=-aldt*apj3(jc)
      acj(kc,jc)=1.+aldt*acj3(jc)
      amj(kc,jc)=-aldt*amj3(jc)
   11 continue
      do 21 kc=1,n3m
      fj(kc,1)=rhs(ic,1,kc)+am*aldt*dq3s(ic,kc)*(1-islv3s)
      fj(kc,n2m)=rhs(ic,n2m,kc)+aldt*ap*dq3n(ic,kc)*(1-islv3n)
      do 22 jc=2,n2m
      fj(kc,jc)=rhs(ic,jc,kc)
   22 continue
   21 continue
c
c    tridiagonal solver
c
      call trvbkj(n2m,1,n3m)
c
      do 30 kc=1,n3m
      do 30 jc=1,n2m
      dq(ic,jc,kc)=fj(kc,jc)
   30 continue
   10 continue
      return
      end
c****************************subrout prcalc ***************************
      subroutine prcalc(pre,pr,al)
c
c    the real pressure is evaluated by the phi
c
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      dimension pre(m1,m2,m3),pr(m1,m2,m3)
      common/tstep/dt,beta,ren
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/metria/caj(m2),cac(m2)
c
c    the pressure is evaluated at the center of the box
c    at the near boundary cells newman b.c. for dph are assumed
c    pre=pressure pr=dph .
c
      be=al*beta
      do 1 kc=1,n3m
      kp=kpv(kc)
      km=kmv(kc)
      do 1 jc=1,n2m
      jm=jmv(jc)
      jp=jpv(jc)
      ac2=1./caj(jc)*dx2q
      co8=float(jp-jc)/cac(jp)
      co4=float(jc-jm)/cac(jc)
      do 1 ic=1,n1m
      ip=ipv(ic)
      im=imv(ic)
      pre(ic,jc,kc)=pre(ic,jc,kc)+pr(ic,jc,kc)-be*(
     1  dx1q*(pr(ip,jc,kc)-2.*pr(ic,jc,kc)+pr(im,jc,kc))+
     1  ac2*(co8*pr(ic,jp,kc)-(co8+co4)*pr(ic,jc,kc)+co4*pr(ic,jm,kc))+
     1  dx3q*(pr(ic,jc,kp)-2.*pr(ic,jc,kc)+pr(ic,jc,km))
     1                                            )
    1 continue
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine boucdq(time,tloc)
      include 'param.f'
      common/qwaloo/q1sold(m1,m3),q2sold(m1,m3),q3sold(m1,m3)
      common/qwalou/q1nold(m1,m3),q2nold(m1,m3),q3nold(m1,m3)
      common/qwallo/q1s(m1,m3),q2s(m1,m3),q3s(m1,m3)
      common/qwalup/q1n(m1,m3),q2n(m1,m3),q3n(m1,m3)
      common/dqwalo/dq1s(m1,m3),dq2s(m1,m3),dq3s(m1,m3)
      common/dqwalu/dq1n(m1,m3),dq2n(m1,m3),dq3n(m1,m3)
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/movwal/tosc,uosc
      common/slotfl/flowq2,tau2
      common/slotii/tim0sl
      common/vslot/v2inf(m1,m3)
      pi=2.*asin(1.)
      omeg1=2.*pi/tosc
      omeg2=2.*pi/tau2
      tpert=tim0sl+tau2
      dq1sma=0.
      dq2sma=0.
      dq3sma=0.
c
c   lower wall
c
      do kc=1,n3m
       do ic=1,n1m
       q1sold(ic,kc)=q1s(ic,kc)
       q2sold(ic,kc)=q2s(ic,kc)
       q3sold(ic,kc)=q3s(ic,kc)
       q1s(ic,kc)=0.
c
c   eventual periodic transpiration
c
c      q2s(ic,kc)=v2inf(ic,kc)*sin(omeg2*tloc)
       q2s(ic,kc)=v2inf(ic,kc)*tloc
       q3s(ic,kc)=uosc*tloc
       dq1s(ic,kc)=q1s(ic,kc)-q1sold(ic,kc)
       dq2s(ic,kc)=q2s(ic,kc)-q2sold(ic,kc)
       dq3s(ic,kc)=q3s(ic,kc)-q3sold(ic,kc)
       dq1sma=max(abs(dq1s(ic,kc)),dq1sma)
       dq2sma=max(abs(dq2s(ic,kc)),dq2sma)
       dq3sma=max(abs(dq3s(ic,kc)),dq3sma)
       end do
      end do
c
c   upper wall at the moment set as no-slip wall
c   slip conditions for q1 and q3 are directly assigned in invtr1 and invtr3
c
      dq1nma=0.
      dq2nma=0.
      dq3nma=0.
      do kc=1,n3m
       do ic=1,n1m
       q1nold(ic,kc)=q1n(ic,kc)
       q2nold(ic,kc)=q2n(ic,kc)
       q3nold(ic,kc)=q3n(ic,kc)
       q1n(ic,kc)=0.
       q2n(ic,kc)=0.
       q3n(ic,kc)=uosc*tloc
       dq1n(ic,kc)=0.
       dq2n(ic,kc)=0.
       dq3n(ic,kc)=q3n(ic,kc)-q3nold(ic,kc)
       dq1nma=max(abs(dq1n(ic,kc)),dq1nma)
       dq2nma=max(abs(dq2n(ic,kc)),dq2nma)
       dq3nma=max(abs(dq3n(ic,kc)),dq3nma)
       end do
      end do
      if(uosc.gt.0.and.time.lt.tpert) then
      write(95,*)'t,dqisma  ',time,tloc,dq1sma,dq2sma,dq3sma
     1           ,dq1nma,dq2nma,dq3nma
c     print *,'q1 parete inferiore',tloc,sin(omega*tloc)
c     print *,'dq1 parete inferiore',dq1s(8,8)
      endif
      return
      end
