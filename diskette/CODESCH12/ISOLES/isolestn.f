c  **************  subrout tschem ***************************************
c   here the time integration of the transport
c   equations is performed
c
c
      subroutine tschem(q,pr,al,ga,ro,ru,qcap,dph,rho,ru4)
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      dimension q(ndv,m1,m2,m3),qcap(m1,m2,m3)
     1         ,dph(m1,m2,m3),pr(m1,m2,m3)
      dimension ru(ndv,m1,m2,m3)
      common/rot/f0
      common/velmax/vmax(ndv),vmaxo(ndv)
      common/tstep/dt,beta,ren
      dimension rho(m1,m2,m3),ru4(m1,m2,m3)
      dimension hro(m1,m2,m3)
      common/strat1/istrat,rho0,g,schm
      common/strat3/grbar(m3)
c
c   time integration implicit viscous adams baschfort nonlinear
c   or runge-kutta depending on the value of nsst before this
c   subroutine is called
c use rhs=dq1, dph=dq2, qcap=dq3 as temporary storage locations for the dq's
      if(n1m.gt.1) then
      call hdnl1(q)
                   endif
      call hdnl2(q,dph)
c
c  *****  computation of nonlinear terms
c use  qcap=dq3 as temporary storage location
c calculate all terms that are explicit for the time-step in hdnl...
      if(n3m.gt.1) then
      call hdnl3(q,qcap,rho)
                   endif
      if(istrat.eq.1) then
c
c  *****  computation of nonlinear terms  for density equation
c
       call hdnl4(rho,q,hro)
                      endif
c
c add rotation -fxu=-f(0,0,1)x(u1,u2,u3)=f(u2,-u1,0)  to the momentum equations
c  *****  solve the dqhat=qhat-q(n) momentum equations
c
      if(ren.gt..1e+15)     then
      call invisc(q,al,ga,ro,dph,qcap,pr,ru)
                            else
      if(n1m.gt.1) then
      call invtr1(q,al,ga,ro,pr,ru)
                   endif
      call invtr2(q,dph,al,ga,ro,pr,ru)
      if(n3m.gt.1) then
      call invtr3(q,qcap,al,ga,ro,pr,ru)
                   endif
                            endif
c
c
c dq=dq1, dph=dq2, qcap=dq3 now contain the 'divergenceful'
c corrections to the velocity field
c q now contains the 'divergenceful' velocity field
c
c  ********* calculation of divg(dqhat)
c
      call divg(qcap,q,al)
c qcap now contains the divergence of the velocity field q
c
c  ********* calculation of the pressure dph by fft in two
c            directions and tridiag in vertical
      call phcalc(qcap,dph)
c
c  ********* calculation of solenoidal vel field
c
      call updvp(dph,q,al)
c     call vmaxv(q)
c     print *,'dopo updvp',vmax(1),vmax(2),vmax(3)
c
c  ********* calculation of pressure field
c
      call prcalc(pr,dph,al)
c
      if(istrat.eq.1) then
        if(ren.lt..1e+12) then
         call invtr4(rho,hro,al,ga,ro,ru4)
                          endif
        if(ren.gt..1e+07) then
         call invti4(rho,hro,al,ga,ro,ru4)
                          endif
                      endif
      return
      end
c  ****************************** subrout invtr1  **********************
c
c   this subroutine performs the inversion of the q1 momentum equation
c   by a factored implicit scheme, only the derivatives 11,22,33 of q1
c   are treated implicitly
c
      subroutine invtr1(q,al,ga,ro,pr,ru)
c
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      common/amcpj1/amj1(m2),acj1(m2),apj1(m2)
      common/trici/ami(m2,m1),aci(m2,m1),api(m2,m1),fi(m2,m1),
     1             fei(m2,m1),qq(m2,m1),ss(m2,m1)
      common/tvpkj/amj(m3,m2),acj(m3,m2),apj(m3,m2),fj(m3,m2),
     1             fej(m3,m2),qj(m3,m2),sj(m3,m2)
      common/tvpjk/amk(m2,m3),ack(m2,m3),apk(m2,m3),fk(m2,m3),
     1             qek(m2,m3),qk(m2,m3),sk(m2,m3)
      dimension q(ndv,m1,m2,m3)
      common/rhsc/rhs(m1,m2,m3)
      dimension pr(m1,m2,m3)
      dimension ru(ndv,m1,m2,m3)
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/tstep/dt,beta,ren
      common/sma/vis(m1,m2,m3)
      common/ledat/ics0,cvisc,ifiltr,ibox
c
c  ********* compute the rhs of the factored equation
c  the rhs include h(n), h(n-1)=ru, grad(pr) component
c  and the 11, 22,3 second derivat. of q1(n)
c  everything at i,j+1/2,k=1/2
c  dq=qhat-q(n)
c
      do 15 kc=1,n3m
      km=kmv(kc)
      kp=kpv(kc)
      do 15 jc=1,n2m
      jm=jmv(jc)
      jp=jpv(jc)
      do 15 ic=1,n1m
      ip=ipv(ic)
      im=imv(ic)
      visc=(vis(ic,jc,kc)+vis(im,jc,kc)+vis(ic,jm,kc)+vis(im,jm,kc))*.25
      visd=(vis(ic,jc,kc)+vis(im,jc,kc)+vis(ic,jp,kc)+vis(im,jp,kc))*.25
c
c   11 second deriv. of q1(n)
c
      d11q1=(vis(ic,jc,kc)*q(1,ip,jc,kc)
     1      -(vis(ic,jc,kc)+vis(im,jc,kc))*q(1,ic,jc,kc)
     1      +vis(im,jc,kc)*q(1,im,jc,kc))*dx1q
c
c   22 second deriv. of q1(n)
c
      d22q1=(visd*q(1,ic,jp,kc)-(visd+visc)*q(1,ic,jc,kc)+
     1       visc*q(1,ic,jm,kc))*dx2q
c
c   add the 33 derivat.
c
      visb=(vis(ic,jc,kc)+vis(im,jc,kc)+vis(ic,jc,kp)+vis(im,jc,kp))*.25
      visa=(vis(ic,jc,kc)+vis(im,jc,kc)+vis(ic,jc,km)+vis(im,jc,km))*.25
      d33q1=(visb*q(1,ic,jc,kp)-(visb+visa)*q(1,ic,jc,kc)+
     1      visa*q(1,ic,jc,km))*dx3q
      dcq1=d11q1+d33q1+d22q1
c
c   grad(pr) along 1
c
      dpx11=(pr(ic,jc,kc)-pr(im,jc,kc))*dx1
      gradp=dpx11*al
      rhsccc=(ga*rhs(ic,jc,kc)+ro*ru(1,ic,jc,kc)-gradp
     1             +al*dcq1)*dt
      ru(1,ic,jc,kc)=rhs(ic,jc,kc)
      rhs(ic,jc,kc)=rhsccc
   15 continue
c
c  ********* compute dq1*  sweeping in the x1 direction
c            periodic
      if(n1m.gt.1) then
      betadx=.5*dx1q*al*dt
      do 1 kc=1,n3m
      do 9 jc=1,n2m
      do 9 ic=1,n1m
      im=imv(ic)
      api(jc,ic)=-betadx*vis(ic,jc,kc)
      aci(jc,ic)=1.+betadx*(vis(ic,jc,kc)+vis(im,jc,kc))
      ami(jc,ic)=-betadx*vis(im,jc,kc)
      fi(jc,ic)=rhs(ic,jc,kc)
    9 continue
c
      call tripvi(1,n1m,1,n2m)
c
      do 3 jc=1,n2m
      do 3 ic=1,n3m
      rhs(ic,jc,kc)=fi(jc,ic)
    3 continue
    1 continue
                   endif
c
c  ************ compute  from dq** sweeping along the x3 direction
c               periodic
c
      if(n3m.gt.1) then
      betadz=.5*dx3q*al*dt
      do 6 ic=1,n1m
      im=imv(ic)
      do 8 kc=1,n3m
      km=kmv(kc)
      kp=kpv(kc)
      do 8 jc=1,n2m
      visb=(vis(ic,jc,kc)+vis(im,jc,kc)+vis(ic,jc,kp)+vis(im,jc,kp))
     1      *.25
      visa=(vis(ic,jc,kc)+vis(im,jc,kc)+vis(ic,jc,km)+vis(im,jc,km))
     1      *.25
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
                   endif
c
c  ************ compute dq1  sweeping along the x2 direction
c
      aldt=al*dt*.5*dx2q
      do 110 ic=1,n1m
      im=imv(ic)
      ip=ipv(ic)
      do 120 jc=1,n2m
      jm=jmv(jc)
      jp=jpv(jc)
      do 120 kc=1,n3m
      visc=(vis(ic,jc,kc)+vis(im,jc,kc)+vis(ic,jm,kc)+vis(im,jm,kc))
     1      *.25
      visd=(vis(ic,jc,kc)+vis(im,jc,kc)+vis(ic,jp,kc)+vis(im,jp,kc))
     1      *.25
      apj(kc,jc)=-aldt*visd
      acj(kc,jc)=1.+aldt*(visd+visc)
      amj(kc,jc)=-aldt*visc
      fj(kc,jc)=rhs(ic,jc,kc)
  120 continue
c
      call trvpkj(1,n2m,1,n3m)
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
      common/tvpkj/amj(m3,m2),acj(m3,m2),apj(m3,m2),fj(m3,m2),
     1             fej(m3,m2),qj(m3,m2),sj(m3,m2)
      common/tvpjk/amk(m2,m3),ack(m2,m3),apk(m2,m3),fk(m2,m3)
     1         ,qek(m2,m3),qk(m2,m3),sk(m2,m3)
      dimension ru(ndv,m1,m2,m3)
      common/rhsc/rhs(m1,m2,m3)
      dimension q(ndv,m1,m2,m3),pr(m1,m2,m3)
      dimension dq(m1,m2,m3)
      common/tstep/dt,beta,ren
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/sma/vis(m1,m2,m3)
      common/ledat/ics0,cvisc,ifiltr,ibox
c
c  ********* compute the rhs of the factored equation
c  the rhs include h(n), h(n-1)=ru, grad(pr) component
c  and the 11, 22,3 second derivat. of q2(n)
c  everything at i,j+1/2,k=1/2
c
      do 16 kc=1,n3m
      km=kmv(kc)
      kp=kpv(kc)
      do 16 jc=1,n2m
      jm=jmv(jc)
      jp=jpv(jc)
      do 16 ic=1,n1m
      im=imv(ic)
      ip=ipv(ic)
c
c   11 second derivative of q2
c
      visb=(vis(ic,jc,kc)+vis(ip,jc,kc)+vis(ic,jm,kc)+vis(ip,jm,kc))
     1     *.25
      visa=(vis(ic,jc,kc)+vis(im,jc,kc)+vis(ic,jm,kc)+vis(im,jm,kc))
     1     *.25
      d11q2=(visb*q(2,ip,jc,kc)-(visb+visa)*q(2,ic,jc,kc)
     1      +visa*q(2,im,jc,kc))*dx1q
c
c   22 second derivative of q2
c
      d22q2=( vis(ic,jc,kc)*q(2,ic,jp,kc)-
     1     (vis(ic,jc,kc)+vis(ic,jm,kc))*q(2,ic,jc,kc)+
     1        vis(ic,jm,kc)*q(2,ic,jm,kc))*dx2q
c
c   add 33 second derivative of q2
c
      visb=(vis(ic,jc,kc)+vis(ic,jc,kp)+vis(ic,jm,kc)+vis(ic,jm,kp))*.25
      visa=(vis(ic,jc,kc)+vis(ic,jc,km)+vis(ic,jm,kc)+vis(ic,jm,km))*.25
      d33q2=(visb*q(2,ic,jc,kp)-(visb+visa)*q(2,ic,jc,kc)+
     1      visa*q(2,ic,jc,km))*dx3q
      dcq2=d11q2+d33q2+d22q2
c
c   component of grap(pr) along 2 direction
c
      dpx22=(pr(ic,jc,kc)-pr(ic,jm,kc))*dx2
      gradp=dpx22*al
      rhs(ic,jc,kc)=(ga*dq(ic,jc,kc)+ro*ru(2,ic,jc,kc)-gradp
     1             +al*dcq2)*dt
      ru(2,ic,jc,kc)=dq(ic,jc,kc) 
   16 continue
c
c  ************ compute dq2** sweeping along the x1 direction
c               periodic
c
      if(n1m.gt.1) then
      betadx=dt*.5*dx1q*al
      do 10 kc=1,n3m
      do 21 jc=1,n2m
      jm=jmv(jc)
      do 21 ic=1,n1m
      ip=ipv(ic)
      im=imv(ic)
      visb=(vis(ic,jc,kc)+vis(ip,jc,kc)+vis(ic,jm,kc)+vis(ip,jm,kc))*.25
      visa=(vis(ic,jc,kc)+vis(im,jc,kc)+vis(ic,jm,kc)+vis(im,jm,kc))*.25
      api(jc,ic)=-betadx*visb
      aci(jc,ic)=1.+betadx*(visb+visa)
      ami(jc,ic)=-betadx*visa
      fi(jc,ic)=rhs(ic,jc,kc)
   21 continue
c
      call tripvi(1,n1m,2,n2m)
c
      do 30 jc=1,n2m
      do 30 ic=1,n1m
      rhs(ic,jc,kc)=fi(jc,ic)
   30 continue
   10 continue
                   endif
c
c  ************ compute dq2 sweeping along the x3 direction
c               periodic
c
      if(n3m.gt.1) then
      betadx=dt*.5*dx3q*al
      do 5 ic=1,n1m
      do 8 kc=1,n3m
      kp=kpv(kc)
      km=kmv(kc)
      do 8 jc=1,n2m
      jm=jmv(jc)
      visb=(vis(ic,jc,kc)+vis(ic,jc,kp)+vis(ic,jm,kc)+vis(ic,jm,kp))*.25
      visa=(vis(ic,jc,kc)+vis(ic,jc,km)+vis(ic,jm,kc)+vis(ic,jm,km))*.25
      apk(jc,kc)=-betadx*visb
      ack(jc,kc)=1.+betadx*(visb+visa)
      amk(jc,kc)=-betadx*visa
      fk(jc,kc)=rhs(ic,jc,kc)
    8 continue
c
      call trvpjk(1,n3m,1,n2m)
c
      do 7 kc=1,n3m
      do 7 jc=1,n2m
      rhs(ic,jc,kc)=fk(jc,kc)
    7 continue
    5 continue
                  endif
c
c  ********* compute the dq2* sweeping in the x2 direction
c            wall boundaries direction
      aldt=al*dt*.5*dx2q
      do 1 ic=1,n1m
      do 11 kc=1,n3m
      do 11 jc=1,n2m
      jm=jmv(jc)
      visu=vis(ic,jc,kc)
      visl=vis(ic,jm,kc)
      apj(kc,jc)=-aldt*visu
      acj(kc,jc)=1.+aldt*(visu+visl)
      amj(kc,jc)=-aldt*visl
      fj(kc,jc)=rhs(ic,jc,kc)
   11 continue
c
      call trvpkj(1,n2m,1,n3m)
c
      do 3 kc=1,n3m
      do 3 jc=1,n2
      q(2,ic,jc,kc)=fj(kc,jc)+q(2,ic,jc,kc)
    3 continue
    1 continue
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
      common/trici/ami(m2,m1),aci(m2,m1),api(m2,m1),fi(m2,m1),
     1             fei(m2,m1),qq(m2,m1),ss(m2,m1)
      common/tvpkj/amj(m3,m2),acj(m3,m2),apj(m3,m2),fj(m3,m2),
     1             fej(m3,m2),qj(m3,m2),sj(m3,m2)
      common/tvpjk/amk(m2,m3),ack(m2,m3),apk(m2,m3),fk(m2,m3),
     1             qek(m2,m3),qk(m2,m3),sk(m2,m3)
      dimension ru(ndv,m1,m2,m3)
      common/rhsc/rhs(m1,m2,m3)
      dimension q(ndv,m1,m2,m3),pr(m1,m2,m3)
      dimension dq(m1,m2,m3)
      common/tstep/dt,beta,ren
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/sma/vis(m1,m2,m3)
      common/ledat/ics0,cvisc,ifiltr,ibox
c
c
c  ********* compute the rhs of the factored equation
c  the rhs include h(n), h(n-1)=ru, grad(pr) component
c  and the 11, 22,3 second derivat. of q3(n)
c  everything at i,j+1/2,k=1/2
c
      do 18 kc=1,n3m
      km=kmv(kc)
      kp=kpv(kc)
      do 18 jc=1,n2m
      jm=jmv(jc)
      jp=jpv(jc)
      do 18 ic=1,n1m
      im=imv(ic)
      ip=ipv(ic)
c
c   11 second derivatives of q3
c
      visb=(vis(ic,jc,km)+vis(ic,jc,kc)+vis(ip,jc,kc)+vis(ip,jc,km))*.25
      visa=(vis(im,jc,km)+vis(im,jc,kc)+vis(ic,jc,kc)+vis(ic,jc,km))*.25
      dq31=(visb*q(3,ip,jc,kc)-(visb+visa)*q(3,ic,jc,kc)+
     1      visa*q(3,im,jc,kc))*dx1q
c
c   22 second derivatives of q3
c

      visc=(vis(ic,jc,km)+vis(ic,jc,kc)+vis(ic,jm,km)+vis(ic,jm,kc))
     1     *.25
      visd=(vis(ic,jc,km)+vis(ic,jc,kc)+vis(ic,jp,km)+vis(ic,jp,kc))
     1     *.25
      dq32=(q(3,ic,jp,kc)*visd
     1     -q(3,ic,jc,kc)*(visd+visc)
     1     +q(3,ic,jm,kc)*visc)*dx2q
c
c   33 second derivatives of q3
c
      dq33=(vis(ic,jc,kc)*q(3,ic,jc,kp)-(vis(ic,jc,kc)+vis(ic,jc,km))*
     1     q(3,ic,jc,kc)+vis(ic,jc,km)*q(3,ic,jc,km))*dx3q
      dcq3=dq31+dq33+dq32
c
c  component of grad(pr) along x3 direction
c
      dpx33=(pr(ic,jc,kc)-pr(ic,jc,km))*dx3
      rhs(ic,jc,kc)=(ga*dq(ic,jc,kc)+ro*ru(3,ic,jc,kc)+al*dcq3
     1              -al*dpx33)*dt
      ru(3,ic,jc,kc)=dq(ic,jc,kc) 
   18 continue
c
c  ************ compute dq3** sweeping along the x1 direction
c
      if(n1m.gt.1) then
      betadx=dt*.5*dx1q*al
      do 4 kc=1,n3m
      km=kmv(kc)
      do 14 jc=1,n2m
      do 14 ic=1,n1m
      im=imv(ic)
      ip=ipv(ic)
      visb=(vis(ic,jc,km)+vis(ic,jc,kc)+vis(ip,jc,kc)+vis(ip,jc,km))*.25
      visa=(vis(im,jc,km)+vis(im,jc,kc)+vis(ic,jc,kc)+vis(ic,jc,km))*.25
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
                   endif
c
c
c  ********* compute the dq3* sweeping in the x3 direction
c
      if(n3m.gt.1) then
      betadx=dt*.5*dx3q*al
      do 1 ic=1,n1m
      do 15 kc=1,n3m
      km=kmv(kc)
      do 15 jc=1,n2m
      apk(jc,kc)=-betadx*vis(ic,jc,kc)
      ack(jc,kc)=1.+betadx*(vis(ic,jc,kc)+vis(ic,jc,km))
      amk(jc,kc)=-betadx*vis(ic,jc,km)
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
                   endif
c
c  ************ compute dq3 sweeping along the x2 direction
c
      aldt=al*dt*.5*dx2q
      do 10 ic=1,n1m
      do 11 kc=1,n3m
      km=kmv(kc)
      do 11 jc=1,n2m
      jm=jmv(jc)
      jp=jpv(jc)
      visc=(vis(ic,jc,km)+vis(ic,jc,kc)+vis(ic,jm,km)+vis(ic,jm,kc))
     1      *.25
      visd=(vis(ic,jc,km)+vis(ic,jc,kc)+vis(ic,jp,km)+vis(ic,jp,kc))
     1      *.25
      apj(kc,jc)=-aldt*visd
      acj(kc,jc)=1.+aldt*(visd+visc)
      amj(kc,jc)=-aldt*visc
      fj(kc,jc)=rhs(ic,jc,kc)
   11 continue
c
      call trvpkj(1,n2m,1,n3m)
c
      do 30 kc=1,n3m
      do 30 jc=1,n2m
      q(3,ic,jc,kc)=fj(kc,jc)+q(3,ic,jc,kc)
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
      common/sma/vis(m1,m2,m3)
c
c    the pressure is evaluated at the center of the box
c    at the near boundary cells newman b.c. for dph are assumed
c    pre=pressure pr=dph .
c
      be=al*dt*.5
      ac1=dx1q
      ac2=dx2q
      do 1 kc=1,n3m
      kp=kpv(kc)
      km=kmv(kc)
      do 1 jc=1,n2m
      jm=jmv(jc)
      jp=jpv(jc)
      do 1 ic=1,n1m
      ip=ipv(ic)
      im=imv(ic)
      vb=.5*(vis(ic,jc,kc)+vis(ip,jc,kc))
      va=.5*(vis(ic,jc,kc)+vis(im,jc,kc))
      vd=.5*(vis(ic,jc,kc)+vis(ic,jp,kc))
      vc=.5*(vis(ic,jc,kc)+vis(ic,jm,kc))
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
c  ****************************** subrout invisc  **********************
c
c   this subroutine performs the inversion of the 
c   three momentum equations in the inviscid case
c   Also the term due to the s.g.s models are discarded
c   this case is done only to check the energy conservation
c   by this run was understood to locate the rho at the
c   same location of q3 to permit the potential
c   energy to become kinetic energy and viceversa
c
      subroutine invisc(q,al,ga,ro,hnl2,hnl3,pr,ru)
c
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      dimension q(ndv,m1,m2,m3),hnl2(m1,m2,m3),hnl3(m1,m2,m3)
      common/rhsc/rhs(m1,m2,m3)
      dimension pr(m1,m2,m3)
      dimension ru(ndv,m1,m2,m3)
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/tstep/dt,beta,ren
c
c  ********* compute the rhs of the q1 equation
c  the rhs include h(n), h(n-1)=ru, grad(pr) component
c  everything at i,j+1/2,k+1/2
c  dq=qhat-q(n)
c
      do kc=1,n3m
            do jc=1,n2m
                  do ic=1,n1m
      im=imv(ic)
c
c   grad(pr) along 1
c

      dpx11=(pr(ic,jc,kc)-pr(im,jc,kc))*dx1
      gradp=dpx11*al
      rhsccc=(ga*rhs(ic,jc,kc)+ro*ru(1,ic,jc,kc)-gradp)*dt
      ru(1,ic,jc,kc)=rhs(ic,jc,kc)
      q(1,ic,jc,kc)=q(1,ic,jc,kc)+rhsccc
                  enddo
            enddo
      enddo
c
c  ********* compute the rhs of the q2 equation
c  the rhs include h(n), h(n-1)=ru, grad(pr) component
c  everything at i+1/2,j,k+1/2
c
      do kc=1,n3m
            do jc=1,n2m
      jm=jmv(jc)
                  do ic=1,n1m
c
c   component of grap(pr) along 2 direction
c
      dpx22=(pr(ic,jc,kc)-pr(ic,jm,kc))*dx2
      gradp=dpx22*al
      rhsccc=(ga*hnl2(ic,jc,kc)+ro*ru(2,ic,jc,kc)-gradp)*dt
      ru(2,ic,jc,kc)=hnl2(ic,jc,kc)
      q(2,ic,jc,kc)=q(2,ic,jc,kc)+rhsccc
                  enddo
            enddo
      enddo
c
c
c  ********* compute the rhs of the q3 equation
c  the rhs include h(n), h(n-1)=ru, grad(pr) component
c  everything at i+1/2,j+1/2,k
c
      do kc=1,n3m
      km=kmv(kc)
            do jc=1,n2m
                  do ic=1,n1m
c
c  component of grad(pr) along x3 direction
c
      dpx33=(pr(ic,jc,kc)-pr(ic,jc,km))*dx3
      rhsccc=(ga*hnl3(ic,jc,kc)+ro*ru(3,ic,jc,kc)-al*dpx33)*dt
      ru(3,ic,jc,kc)=hnl3(ic,jc,kc)
      q(3,ic,jc,kc)=q(3,ic,jc,kc)+rhsccc
                  enddo
            enddo
      enddo
      return
      end

