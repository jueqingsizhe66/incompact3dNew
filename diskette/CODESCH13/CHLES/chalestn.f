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
c               v1dgf=component 1 of grad(dph)
c
      do 1 kc=1,n3m
      do 1 jc=1,n2m
      do 1 ic=1,n1m
      dfx11=(dph(ic,jc,kc)-dph(imv(ic),jc,kc))*dx1
      q(1,ic,jc,kc)=q(1,ic,jc,kc)-dfx11*dt*al
    1 continue
c
c  ***********  compute the q2 velocity component
c               v2dgf=component 2 of grad(dph)
c
      do 2 kc=1,n3m
      do 2 jc=2,n2m
      sucac=1./cac(jc)
      do 2 ic=1,n1m
      dfx22=(dph(ic,jc,kc)-dph(ic,jmv(jc),kc))*dx2
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
c               v3dgf=component 3 of grad(dph)
c
      do 5 kc=1,n3m
      do 5 jc=1,n2m
      do 5 ic=1,n1m
      dfx33=(dph(ic,jc,kc)-dph(ic,jc,kmv(kc)))*dx3
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
c  jc=1  spanwise velocity
c   in the case of q1  assigned
      ucaj=4./(2.*cac(2)+caj(1))
      apj1(1)=1./cac(2)*ucaj
      acj1(1)=(1./cac(2)+2./cac(1))*ucaj
      amj1(1)=0.
                       else
c   in the case of q1 shear free  
      ucaj=1./caj(1)
      apj1(1)=1./cac(2)*ucaj
      acj1(1)=1./cac(2)*ucaj
      amj1(1)=0.
                       endif
      if(islv3s.eq.0) then
c  jc=1  streamwise velocity
c   in the case of  q3 assigned
      ucaj=4./(2.*cac(2)+caj(1))
      apj3(1)=1./cac(2)*ucaj
      acj3(1)=(1./cac(2)+2./cac(1))*ucaj
      amj3(1)=0.
                       else
c   in the case of q3  shear free
      ucaj=1./caj(1)
      apj3(1)=1./cac(2)*ucaj
      acj3(1)=1./cac(2)*ucaj
      amj3(1)=0.
                       endif
c  jc=n2m streamwise velocity
c   in the case of  q1 assigned
      if(islv1n.eq.0) then
      ucaj=4./(2.*cac(n2m)+caj(n2m))
      apj1(n2m)=0.
      acj1(n2m)=(2./cac(n2)+1./cac(n2m))*ucaj
      amj1(n2m)=1./cac(n2m)*ucaj
                       else
c   in the case of q1  shear free
      ucaj=1./caj(n2m)
      apj1(n2m)=0.
      acj1(n2m)=1./cac(n2m)*ucaj
      amj1(n2m)=1./cac(n2m)*ucaj
                       endif
c   in the case of  q3 assigned
      if(islv3n.eq.0) then
      ucaj=4./(2.*cac(n2m)+caj(n2m))
      apj3(n2m)=0.
      acj3(n2m)=(2./cac(n2)+1./cac(n2m))*ucaj
      amj3(n2m)=1./cac(n2m)*ucaj
                       else
c   in the case of q3  shear free
      ucaj=1./caj(n2m)
      apj3(n2m)=0.
      acj3(n2m)=1./cac(n2m)*ucaj
      amj3(n2m)=1./cac(n2m)*ucaj
                       endif
c
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
      open(17,file='coeff')
      do j=1,n2m
      write(17,*) j,apj1(j),acj1(j),amj1(j),
     1              apj2(j),acj2(j),amj2(j)
      end do
      close(17)
      return
      end
c
c  **************  subrout tschem ***************************************
c
c
      subroutine tschem(q,pr,ru,qcap,dph,time,ncount)
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
      common/ledat/csma,cvisc,iles
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
cccccccccccccccccccccccccccccccccccccccccccccccc
c
c Dynamic model:  Lilly version 
c
cccccccccccccccccccccccccccccccccccccccccccccccc
      if(ns.eq.1.and.iles.eq.1)  then
        ncount=ncount+1
        write(61,*)ncount,iles
        call dynamic(q,dph,ncount)
                                 end if
cccccccccccccccccccccccccccccccccccccccccccccccc
c
c Smagorinsky model
c
cccccccccccccccccccccccccccccccccccccccccccccccc
       if(ns.eq.1.and.iles.eq.3)  then
        ncount=ncount+1
        write(61,*)ncount,iles
        call smago(q,dph)
                                  endif
cccccccccccccccccccccccccccccccccccccccccccccccc
c
c Constant nu turb. =<nu_T(i,j,k)>
c
cccccccccccccccccccccccccccccccccccccccccccccccc
       if(ns.eq.1.and.iles.eq.5)  then
        ncount=ncount+1
        write(61,*)ncount,iles
        call consma(q,dph)
                                  endif
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
c  *****  computation of nonlinear terms
c
c 
      call hdnl1(q)
c 
      call hdnl2(q,dph)
c
      call hdnl3(q,qcap)
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
c   by a factored implicit scheme, only the derivatives 11,22,33 of q1
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
      common/sma/vis(m1,m2,m3)
      common/ledat/csma,cvisc,iles

c
c  ********* compute the rhs of the factored equation
c  the rhs include h(n), h(n-1)=ru, grad(pr) component
c  and the 11, 22,3 second derivat. of q1(n)
c  everything at i,j+1/2,k=1/2
c  dq=qhat-q(n)
c
      n2mm=n2m-1
      do kc=1,n3m
      do jc=1,n2m
      do ic=1,n1m
c
c   11 second deriv. of q1(n)
c
      d11q1=(vis(ic,jc,kc)*q(1,ipv(ic),jc,kc)
     1     -(vis(ic,jc,kc)+vis(imv(ic),jc,kc))*
     1      q(1,ic,jc,kc)+vis(imv(ic),jc,kc)*q(1,imv(ic),jc,kc))*dx1q
c
c   add the 33 derivat.
c
      visb=(vis(ic,jc,kc)+vis(imv(ic),jc,kc)
     1     +vis(ic,jc,kpv(kc))+vis(imv(ic),jc,kpv(kc)))*.25
      visa=(vis(ic,jc,kc)+vis(imv(ic),jc,kc)
     1     +vis(ic,jc,kmv(kc))+vis(imv(ic),jc,kmv(kc)))*.25
      d33q1=(visb*q(1,ic,jc,kpv(kc))-(visb+visa)*q(1,ic,jc,kc)+
     1      visa*q(1,ic,jc,kmv(kc)))*dx3q
      dcq1=d11q1+d33q1
c
c   grad(pr) along 1
c
      dpx11=(pr(ic,jc,kc)-pr(imv(ic),jc,kc))*dx1
      gradp=dpx11*al
      rhsccc=(ga*rhs(ic,jc,kc)+ro*ru(1,ic,jc,kc)-gradp
     1             +al*dcq1)*dt
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
      do ic=1,n1m
      visc=(vis(ic,jc,kc)+vis(imv(ic),jc,kc)
     1     +vis(ic,jmv(jc),kc)+vis(imv(ic),jmv(jc),kc))*.25
      visd=(vis(ic,jc,kc)+vis(imv(ic),jc,kc)
     1     +vis(ic,jpv(jc),kc)+vis(imv(ic),jpv(jc),kc))*.25
c
c   22 second deriv. of q1(n)
c
      d22q1=(apj1(jc)*q(1,ic,jpv(jc),kc)*visd
     1      -(apj1(jc)*visd+amj1(jc)*visc)*q(1,ic,jc,kc)
     1      +amj1(jc)*q(1,ic,jmv(jc),kc)*visc)*dx2q
      rhs(ic,jc,kc)=rhs(ic,jc,kc)+al*d22q1*dt
      enddo
      enddo
      enddo
c
c   add the second derivative d22q1 at the lower wall 
c
      jc=1
      am=4./(2.*cac(2)+caj(1))*2./cac(jc)
      do kc=1,n3m
      do ic=1,n1m
      visc=cvisc
      visd=(vis(ic,jc,kc)+vis(imv(ic),jc,kc)
     1     +vis(ic,jc+1,kc)+vis(imv(ic),jc+1,kc))*.25
c
c   22 second deriv. of q1(n) O(dx2^2)
c
      d22q1=(apj1(jc)*q(1,ic,jc+1,kc)*visd
     1      -(apj1(jc)*visd+am*visc)*q(1,ic,jc,kc)
     1      +am*q1s(ic,kc)*(1-islv1s)*visc )*dx2q
      rhs(ic,jc,kc)=rhs(ic,jc,kc)+al*d22q1*dt
      enddo
      enddo
c
c   add the second derivative d22q1 aqt the upper wall 
c
      jc=n2m
      ap=4./(2.*cac(n2m)+caj(n2m))*2./cac(n2)
      do kc=1,n3m
      do ic=1,n1m
      visc=(vis(ic,jc,kc)+vis(imv(ic),jc,kc)
     1     +vis(ic,jc-1,kc)+vis(imv(ic),jc-1,kc))*.25
      visd=cvisc

c
c   22 second deriv. of q1(n) O(dx2^2)
c
      d22q1=(amj1(jc)*q(1,ic,jc-1,kc)*visc
     1      -(amj1(jc)*visc+ap*visd)*q(1,ic,jc,kc)
     1      +ap*q1n(ic,kc)*(1-islv1n)*visd )*dx2q
      rhs(ic,jc,kc)=rhs(ic,jc,kc)+al*d22q1*dt
      enddo
      enddo
c
c  ********* compute dq1*  sweeping in the x1 direction
c            periodic
      betadx=dx1q*al*beta
      do 1 kc=1,n3m
      do 9 jc=1,n2m
      do 9 ic=1,n1m
      api(jc,ic)=-betadx*vis(ic,jc,kc)
      aci(jc,ic)=1.+betadx*(vis(ic,jc,kc)+vis(imv(ic),jc,kc))
      ami(jc,ic)=-betadx*vis(imv(ic),jc,kc)
      fi(jc,ic)=rhs(ic,jc,kc)
    9 continue
c
      call tripvi(1,n1m,1,n2m)
c
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
      visb=(vis(ic,jc,kc)+vis(imv(ic),jc,kc)
     1     +vis(ic,jc,kpv(kc))+vis(imv(ic),jc,kpv(kc)))*.25
      visa=(vis(ic,jc,kc)+vis(imv(ic),jc,kc)
     1     +vis(ic,jc,kmv(kc))+vis(imv(ic),jc,kmv(kc)))*.25
      apk(jc,kc)=-betadz*visb
      ack(jc,kc)=1.+betadz*(visb+visa)
      amk(jc,kc)=-betadz*visa
      fk(jc,kc)=rhs(ic,jc,kc)
    8 continue
c
      call trvpjk(1,n3m,1,n2m)
c
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
      do 120 jc=2,n2m-1
      do 120 kc=1,n3m
      visc=(vis(ic,jc,kc)+vis(imv(ic),jc,kc)
     1     +vis(ic,jmv(jc),kc)+vis(imv(ic),jmv(jc),kc))*.25
      visd=(vis(ic,jc,kc)+vis(imv(ic),jc,kc)
     1     +vis(ic,jpv(jc),kc)+vis(imv(ic),jpv(jc),kc))*.25
      apj(kc,jc)=-aldt*apj1(jc)*visd
      acj(kc,jc)=1.+aldt*(apj1(jc)*visd+amj1(jc)*visc)
      amj(kc,jc)=-aldt*amj1(jc)*visc
  120 continue
      do 130 kc=1,n3m
      visd=(vis(ic,1,kc)+vis(imv(ic),1,kc)
     1     +vis(ic,2,kc)+vis(imv(ic),2,kc))*.25
      visc=(vis(ic,n2m,kc)+vis(imv(ic),n2m,kc)
     1     +vis(ic,n2m-1,kc)+vis(imv(ic),n2m-1,kc))*.25
      apj(kc,1)=-aldt*visd*apj1(1)
      acj(kc,1)=1.+aldt*(visd*apj1(1)+cvisc*am)
      amj(kc,1)=0.
      apj(kc,n2m)=0.
      acj(kc,n2m)=1.+aldt*(cvisc*ap+visc*amj1(n2m))
      amj(kc,n2m)=-aldt*visc*amj1(n2m)
  130 continue

      do 135 kc=1,n3m
      fj(kc,1)=rhs(ic,1,kc)+aldt*am*dq1s(ic,kc)*(1-islv1s)*cvisc
      fj(kc,n2m)=rhs(ic,n2m,kc)+aldt*ap*dq1n(ic,kc)*(1-islv1n)*cvisc
      do 136 jc=2,n2mm
      fj(kc,jc)=rhs(ic,jc,kc)
  136 continue
  135 continue
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
c   by a factored implicit scheme, only the derivatives 11,22,33 of q2
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
      common/sma/vis(m1,m2,m3)
      common/ledat/csma,cvisc,iles
c
c  ********* compute the rhs of the factored equation
c  the rhs include h(n), h(n-1)=ru, grad(pr) component
c  and the 11, 22,3 second derivat. of q2(n)
c  everything at i,j+1/2,k=1/2
c
      do 16 kc=1,n3m
      do 16 jc=2,n2m
      sucac=1./cac(jc)
      do 16 ic=1,n1m
c
c   11 second derivative of q2
c
      visb=(vis(ic,jc,kc)+vis(ipv(ic),jc,kc)
     1     +vis(ic,jmv(jc),kc)+vis(ipv(ic),jmv(jc),kc))*.25
      visa=(vis(ic,jc,kc)+vis(imv(ic),jc,kc)
     1     +vis(ic,jmv(jc),kc)+vis(imv(ic),jmv(jc),kc))*.25
      d11q2=(visb*q(2,ipv(ic),jc,kc)-(visb+visa)*q(2,ic,jc,kc)
     1      +visa*q(2,imv(ic),jc,kc))*dx1q
c
c   22 second derivative of q2
c
      d22q2=(apj2(jc)*q(2,ic,jpv(jc),kc)*vis(ic,jc,kc) 
     1       -(apj2(jc)*vis(ic,jc,kc)
     1        +amj2(jc)*vis(ic,jmv(jc),kc))*q(2,ic,jc,kc)+
     1       amj2(jc)*q(2,ic,jmv(jc),kc)*vis(ic,jmv(jc),kc))*dx2q
c
c   add 33 second derivative of q2
c
      visb=(vis(ic,jc,kc)+vis(ic,jc,kpv(kc))
     1     +vis(ic,jmv(jc),kc)+vis(ic,jmv(jc),kpv(kc)))*.25
      visa=(vis(ic,jc,kc)+vis(ic,jc,kmv(kc))
     1     +vis(ic,jmv(jc),kc)+vis(ic,jmv(jc),kmv(kc)))*.25
      d33q2=(visb*q(2,ic,jc,kpv(kc))-(visb+visa)*q(2,ic,jc,kc)+
     1       visa*q(2,ic,jc,kmv(kc)))*dx3q
      dcq2=d11q2+d33q2+d22q2
c
c   component of grap(pr) along 2 direction
c
      dpx22=(pr(ic,jc,kc)-pr(ic,jmv(jc),kc))*dx2
      gradp=dpx22*al*sucac
      rhs(ic,jc,kc)=(ga*dq(ic,jc,kc)+ro*ru(2,ic,jc,kc)-gradp
     1             +al*dcq2)*dt
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
      visb=(vis(ic,jc,kc)+vis(ipv(ic),jc,kc)
     1     +vis(ic,jmv(jc),kc)+vis(ipv(ic),jmv(jc),kc))*.25
      visa=(vis(ic,jc,kc)+vis(imv(ic),jc,kc)
     1     +vis(ic,jmv(jc),kc)+vis(imv(ic),jmv(jc),kc))*.25
      api(jc,ic)=-betadx*visb
      aci(jc,ic)=1.+betadx*(visb+visa)
      ami(jc,ic)=-betadx*visa
      fi(jc,ic)=rhs(ic,jc,kc)
      fi(jc,ic)=rhs(ic,jc,kc)
   21 continue
c
      call tripvi(1,n1m,2,n2m)
c
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
      visb=(vis(ic,jc,kc)+vis(ic,jc,kpv(kc))
     1     +vis(ic,jmv(jc),kc)+vis(ic,jmv(jc),kpv(kc)))*.25
      visa=(vis(ic,jc,kc)+vis(ic,jc,kmv(kc))
     1     +vis(ic,jmv(jc),kc)+vis(ic,jmv(jc),kmv(kc)))*.25
      apk(jc,kc)=-betadx*visb
      ack(jc,kc)=1.+betadx*(visb+visa)
      amk(jc,kc)=-betadx*visa
      fk(jc,kc)=rhs(ic,jc,kc)
    8 continue
c
      call trvpjk(1,n3m,2,n2m)
c
      do 7 kc=1,n3m
      do 7 jc=2,n2m
      rhs(ic,jc,kc)=fk(jc,kc)
    7 continue
    5 continue


c
c  ********* compute the dq2* sweeping in the x2 direction
c            wall boundaries direction
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
      apj(kc,jc)=-aldt*apj2(jc)*vis(ic,jc,kc)
      acj(kc,jc)=1.+aldt*(vis(ic,jc,kc)*apj2(jc)
     1                   +vis(ic,jmv(jc),kc)*amj2(jc))
      amj(kc,jc)=-aldt*amj2(jc)*vis(ic,jmv(jc),kc)
      fj(kc,jc)=rhs(ic,jc,kc)
   11 continue
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
c   by a factored implicit scheme, only the derivatives 11,22,33 of q3
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
      common/sma/vis(m1,m2,m3)
      common/ledat/csma,cvisc,iles
c
c
c  ********* compute the rhs of the factored equation
c  the rhs include h(n), h(n-1)=ru, grad(pr) component
c  and the 11, 22,3 second derivat. of q3(n)
c  everything at i,j+1/2,k=1/2
c
      s3tot=0.
      n2mm=n2m-1
      do kc=1,n3m
      do jc=1,n2m
      do ic=1,n1m
c
c   11 second derivatives of q3
c
      visb=(vis(ic,jc,kmv(kc))+vis(ic,jc,kc)
     1     +vis(ipv(ic),jc,kc)+vis(ipv(ic),jc,kmv(kc)))*.25
      visa=(vis(imv(ic),jc,kmv(kc))+vis(imv(ic),jc,kc)
     1     +vis(ic,jc,kc)+vis(ic,jc,kmv(kc)))*.25
      dq31=(visb*q(3,ipv(ic),jc,kc)-(visb+visa)*q(3,ic,jc,kc)+
     1      visa*q(3,imv(ic),jc,kc))*dx1q
c
c   33 second derivatives of q3
c
      dq33=(vis(ic,jc,kc)*q(3,ic,jc,kpv(kc))
     1    -(vis(ic,jc,kc)+vis(ic,jc,kmv(kc)))*q(3,ic,jc,kc)
     1     +vis(ic,jc,kmv(kc))*q(3,ic,jc,kmv(kc)))*dx3q
      dcq3=dq31+dq33
      rhs(ic,jc,kc)=(ga*dq(ic,jc,kc)+ro*ru(3,ic,jc,kc)+
     1               al*dcq3)*dt
      s3tot=s3tot+dcq3*caj(jc)
      ru(3,ic,jc,kc)=dq(ic,jc,kc) 
      enddo
      enddo
      enddo
c
c  add 22 second derivatives of q3 inner field
c
      do kc=1,n3m
      do jc=2,n2mm
      do ic=1,n1m
c
c   22 second derivatives of q3
c
      visc=(vis(ic,jc,kmv(kc))+vis(ic,jc,kc)
     1     +vis(ic,jmv(jc),kmv(kc))+vis(ic,jmv(jc),kc))*.25
      visd=(vis(ic,jc,kmv(kc))+vis(ic,jc,kc)
     1     +vis(ic,jpv(jc),kmv(kc))+vis(ic,jpv(jc),kc))*.25
      dq32=(apj3(jc)*q(3,ic,jpv(jc),kc)*visd
     1     -(apj3(jc)*visd+amj3(jc)*visc)*q(3,ic,jc,kc)
     1     +amj3(jc)*q(3,ic,jmv(jc),kc)*visc)*dx2q
      rhs(ic,jc,kc)=rhs(ic,jc,kc)+al*dq32*dt
      s3tot=s3tot+dq32*caj(jc)
      enddo
      enddo
      enddo
c
c  add 22 second derivatives of q3 lower wall   
c
      jc=1
      am=4./(2.*cac(2)+caj(1))*2./cac(jc)
      do kc=1,n3m
      do ic=1,n1m
c
c   22 second deriv. of q3(n) O(dx2^2)
c
      visc=cvisc
      visd=(vis(ic,jc,kmv(kc))+vis(ic,jc,kc)
     1     +vis(ic,jc+1,kmv(kc))+vis(ic,jc+1,kc))*.25
      dq32=(apj3(jc)*q(3,ic,jc+1,kc)*visd
     1     -(apj3(jc)*visd+am*visc)*q(3,ic,jc,kc)
     1     +am*q3s(ic,kc)*(1-islv3s)*visc )*dx2q
      s3tot=s3tot+dq32*caj(jc)
      rhs(ic,jc,kc)=rhs(ic,jc,kc)+al*dq32*dt
      enddo
      enddo
c
c   add the second derivative d22q3 aqt the upper wall
c
      jc=n2m
      ap=4./(2.*cac(n2m)+caj(n2m))*2./cac(n2)
      do kc=1,n3m
      do ic=1,n1m
c
c   22 second deriv. of q3(n) O(dx2^2)
c
      visc=(vis(ic,jc,kmv(kc))+vis(ic,jc,kc)
     1     +vis(ic,jc-1,kmv(kc))+vis(ic,jc-1,kc))*.25
      visd=cvisc
      dq32=(amj3(jc)*q(3,ic,jc-1,kc)*visc
     1     -(amj3(jc)*visc+ap*visd)*q(3,ic,jc,kc)
     1     +ap*q3n(ic,kc)*(1-islv3n)*visd )*dx2q
      s3tot=s3tot+dq32*caj(jc)
      rhs(ic,jc,kc)=rhs(ic,jc,kc)+al*dq32*dt
      enddo
      enddo
      dp3ns=s3tot/(2.*n1m*n2m*n3m)
c
c
      do 38 kc=1,n3m
      do 38 jc=1,n2m
      do 38 ic=1,n1m
c
c  component of grad(pr) along x3 direction
c
      dpx33=(pr(ic,jc,kc)-pr(ic,jc,kmv(kc)))*dx3
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
      apk(jc,kc)=-betadz*vis(ic,jc,kc)
      ack(jc,kc)=1.+betadz*(vis(ic,jc,kc)+vis(ic,jc,kmv(kc)))
      amk(jc,kc)=-betadz*vis(ic,jc,kmv(kc))
      fk(jc,kc)=rhs(ic,jc,kc)
   15 continue
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
      visb=(vis(ic,jc,kmv(kc))+vis(ic,jc,kc)
     1     +vis(ipv(ic),jc,kc)+vis(ipv(ic),jc,kmv(kc)))*.25
      visa=(vis(imv(ic),jc,kmv(kc))+vis(imv(ic),jc,kc)
     1     +vis(ic,jc,kc)+vis(ic,jc,kmv(kc)))*.25
      api(jc,ic)=-betadx*visb
      aci(jc,ic)=1.+betadx*(visb+visa)
      ami(jc,ic)=-betadx*visa
      fi(jc,ic)=rhs(ic,jc,kc)
   14 continue
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
      do 11 jc=2,n2m-1
      visc=(vis(ic,jc,kmv(kc))+vis(ic,jc,kc)
     1     +vis(ic,jmv(jc),kmv(kc))+vis(ic,jmv(jc),kc))*.25
      visd=(vis(ic,jc,kmv(kc))+vis(ic,jc,kc)
     1     +vis(ic,jpv(jc),kmv(kc))+vis(ic,jpv(jc),kc))*.25
      apj(kc,jc)=-aldt*apj3(jc)*visd
      acj(kc,jc)=1.+aldt*(apj3(jc)*visd+amj3(jc)*visc)
      amj(kc,jc)=-aldt*amj3(jc)*visc
   11 continue
      do 21 kc=1,n3m
      visc=(vis(ic,n2m,kmv(kc))+vis(ic,n2m,kc)
     1     +vis(ic,n2m-1,kmv(kc))+vis(ic,n2m-1,kc))*.25
      visd=(vis(ic,1,kmv(kc))+vis(ic,1,kc)
     1     +vis(ic,2,kmv(kc))+vis(ic,2,kc))*.25
      apj(kc,1)=-aldt*visd*apj3(1)
      acj(kc,1)=1.+aldt*(visd*apj3(1)+cvisc*am)
      amj(kc,1)=0.
      apj(kc,n2m)=0.
      acj(kc,n2m)=1.+aldt*(cvisc*ap+visc*amj3(n2m))
      amj(kc,n2m)=-aldt*visc*amj3(n2m)
      fj(kc,1)=rhs(ic,1,kc)+am*aldt*dq3s(ic,kc)*(1-islv3s)*cvisc
      fj(kc,n2m)=rhs(ic,n2m,kc)+aldt*ap*dq3n(ic,kc)*(1-islv3n)*cvisc
      do 22 jc=2,n2m-1
      fj(kc,jc)=rhs(ic,jc,kc)
   22 continue
   21 continue
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
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      dimension pre(m1,m2,m3),pr(m1,m2,m3)
      common/tstep/dt,beta,ren
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/metria/caj(m2),cac(m2)
      common/sma/vis(m1,m2,m3)
c
c    the pressure is evaluated at the center of the box
c    at the near boundary cells newman b.c. for dph are assumed
c    pre=pressure pr=dph .
c
      be=al*beta
      ac1=dx1q
      do 1 kc=1,n3m
      kp=kpv(kc)
      km=kmv(kc)
      do 1 jc=1,n2m
      jm=jmv(jc)
      jp=jpv(jc)
      ac2=1./caj(jc)*dx2q
      do 1 ic=1,n1m
      ip=ipv(ic)
      im=imv(ic)
      vb=.5*(vis(ic,jc,kc)+vis(ip,jc,kc))
      va=.5*(vis(ic,jc,kc)+vis(im,jc,kc))
      vd=.5*(vis(ic,jc,kc)+vis(ic,jp,kc))*float(jp-jc)/cac(jp)
      vc=.5*(vis(ic,jc,kc)+vis(ic,jm,kc))*float(jc-jm)/cac(jc)
      vf=.5*(vis(ic,jc,kc)+vis(ic,jc,kp))
      ve=.5*(vis(ic,jc,kc)+vis(ic,jc,km))
      pre(ic,jc,kc)=pre(ic,jc,kc)+pr(ic,jc,kc)-be*(
     1     ac1*(vb*pr(ip,jc,kc)-(vb+va)*pr(ic,jc,kc)+va*pr(im,jc,kc))+
     1     ac2*(vd*pr(ic,jp,kc)-(vd+vc)*pr(ic,jc,kc)+vc*pr(ic,jm,kc))+
     1    dx3q*(vf*pr(ic,jc,kp)-(vf+ve)*pr(ic,jc,kc)+ve*pr(ic,jc,km)))
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
