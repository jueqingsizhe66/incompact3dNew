c
c  **************  subrout tsinv ***************************************
c  in this routine the time marching is performed in the inviscid case
c  this permits to check the energy conservation properties of
c  the numerical scheme.
c  nsst=3  for Runge-Kutta
c  nsst=1  for Adams Bashfort
c
      subroutine tsinv(q,pr,ru,qcap,dph,time)
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      dimension q(ndv,m1,m2,m3),qcap(m1,m2,m3)
     1         ,dph(m1,m2,m3),pr(m1,m2,m3)
      dimension ru(ndv,m1,m2,m3)
      common/tstep/dt,beta,ren
      common/tscoe/gam(3),rom(3),nsst
      do 2000 ns=1,nsst
c
c   time integration implicit viscous 
c
      al=(gam(ns)+rom(ns))
      ga=gam(ns)
      ro=rom(ns)
c
c    routine for the boudary conditions including
c    the injection in the walls
c
      call boucdq(time,ft)

c  *****  computation of nonlinear terms
c
      call hdnl1(q)
c 
      call hdnl2(q,dph)
c
      call hdnl3(q,qcap)
c
c  *****  solve the dqhat=qhat-q(n) momentum equations
c  
      call invti1(q,al,ga,ro,pr,ru)
c
c
c 
      call invti2(q,dph,al,ga,ro,pr,ru)
c
c
c 
      call invti3(q,qcap,al,ga,ro,pr,ru)
c
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
c  ****************************** subrout invti1  **********************
c
c   this subroutine performs the calculation of the q1 momentum equation
c   in the inviscid case 
      subroutine invti1(q,al,ga,ro,pr,ru)
c
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/metria/caj(m2),cac(m2)
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      dimension q(ndv,m1,m2,m3)
      common/rhsc/rhs(m1,m2,m3)
      dimension pr(m1,m2,m3)
      dimension ru(ndv,m1,m2,m3)
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/tstep/dt,beta,ren

c
c  ********* compute the rhs of the factored equation
c  the rhs include h(n), h(n-1)=ru, grad(pr) component
c  everything at i,j+1/2,k=1/2
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
      rhs(ic,jc,kc)=rhsccc
      q(1,ic,jc,kc)=rhsccc+q(1,ic,jc,kc)
      enddo
      enddo
      enddo
      return
      end
c
c  ****************************** subrout invti2  **********************
c   this subroutine performs the calculation of the q2 momentum equation
c   in the inviscid case 
c
      subroutine invti2(q,dq,al,ga,ro,pr,ru)
c
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      common/metria/caj(m2),cac(m2)
      dimension ru(ndv,m1,m2,m3)
      common/rhsc/rhs(m1,m2,m3)
      dimension q(ndv,m1,m2,m3),pr(m1,m2,m3)
      dimension dq(m1,m2,m3)
      common/tstep/dt,beta,ren
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/qwallo/q1s(m1,m3),q2s(m1,m3),q3s(m1,m3)
      common/qwalup/q1n(m1,m3),q2n(m1,m3),q3n(m1,m3)
c
c  ********* compute the rhs of the factored equation
c  the rhs include h(n), h(n-1)=ru, grad(pr) component
c  everything at i,j+1/2,k=1/2
c
      alre=al/ren
      do 16 kc=1,n3m
      do 16 jc=2,n2m
      jm=jc-1
      sucac=1./cac(jc)
      do 16 ic=1,n1m
c
c   component of grap(pr) along 2 direction
c
      dpx22=(pr(ic,jc,kc)-pr(ic,jm,kc))*dx2
      gradp=dpx22*al*sucac
      rhs(ic,jc,kc)=(ga*dq(ic,jc,kc)+ro*ru(2,ic,jc,kc)-gradp)*dt
      ru(2,ic,jc,kc)=dq(ic,jc,kc) 
      q(2,ic,jc,kc)=rhs(ic,jc,kc)+q(2,ic,jc,kc)
   16 continue
c
      do 9 kc=1,n3m
      do 9 ic=1,n1m
      q(2,ic,1,kc)=q2s(ic,kc)
      q(2,ic,n2,kc)=q2n(ic,kc)
    9 continue
      return
      end
c
c  ****************************** subrout invti3  **********************
c   this subroutine performs the calculation of the q3 momentum equation
c   in the inviscid case 
c
      subroutine invti3(q,dq,al,ga,ro,pr,ru)
c

      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
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
      common/islwal/islv1s,islv1n,islv3s,islv3n
c
c
c  ********* compute the rhs of the factored equation
c  the rhs include h(n), h(n-1)=ru, grad(pr) component
c  everything at i,j+1/2,k=1/2
c
      s3tot=0.
      do kc=1,n3m
      do jc=1,n2m
      do ic=1,n1m
      rhs(ic,jc,kc)=(ga*dq(ic,jc,kc)+ro*ru(3,ic,jc,kc))*dt
      s3tot=s3tot+ru(3,ic,jc,kc)*caj(jc)
      ru(3,ic,jc,kc)=dq(ic,jc,kc) 
      enddo
      enddo
      enddo
      dp3ns=s3tot/(2.*n1m*n2m*n3m)
c
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
      q(3,ic,jc,kc)=rhs(ic,jc,kc)+q(3,ic,jc,kc)
   38 continue
      return
      end
