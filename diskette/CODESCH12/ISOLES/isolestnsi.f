c
c
c  ****************************** subrout hdnl4  **********************
c
c calculate nonlinear advection and diffusion of rho
c in this subroutine are calculated the non-linear terms and sub-grid terms.
c
      subroutine hdnl4(rho,q,hro)
c
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      dimension rho(m1,m2,m3)
      dimension hro(m1,m2,m3)
      dimension q(ndv,m1,m2,m3)
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/smapra/vis(m1,m2,m3)
      common/ledat/ics0,cvisc,ifiltr,ibox
      common/strat1/istrat,rho0,g,schm   
      common/strat3/grbar(m3)
      common/strat2/igrad,bvais
c
c
c
c   **********  compute the non-linear terms by centered differences
c    rho is defined at the same position of q3  for global conservation
c    properties
c                       i+1/2,j+1/2,k
c
      qdx1=dx1*0.25
      qdx2=dx2*0.25
      qdx3=dx3*0.25
      do 10 kc=1,n3m
      km=kmv(kc)
      kp=kpv(kc)
      do 10 jc=1,n2m
      jm=jmv(jc)
      jp=jpv(jc)
      do 10 ic=1,n1m
      ip=ipv(ic)
      im=imv(ic)
c
c
c    rho q1 term
c
      h32=( (q(2,ic,jp,kc)+q(2,ic,jp,km))*
     1      (rho(ic,jp,kc)+rho(ic,jc,kc))
     1     -(q(2,ic,jc,kc)+q(2,ic,jc,km))*
     1      (rho(ic,jc,kc)+rho(ic,jm,kc))
     1    )*qdx2
c
c    rho q2 term
c
      h31=((q(1,ip,jc,kc)+q(1,ip,jc,km))
     1    *(rho(ip,jc,kc)+rho(ic,jc,kc))
     3    -(q(1,ic,jc,kc)+q(1,ic,jc,km))
     1    *(rho(ic,jc,kc)+rho(im,jc,kc))
     1    )*qdx1
c
c    rho q3 term
c
      h33=((q(3,ic,jc,kp)+q(3,ic,jc,kc))
     1    *(rho(ic,jc,kp)+rho(ic,jc,kc))
     3    -(q(3,ic,jc,kc)+q(3,ic,jc,km))
     1    *(rho(ic,jc,kc)+rho(ic,jc,km))
     1    )*qdx3
      hq=h31+h32+h33
      hro(ic,jc,kc)=-hq
   10 continue
c
c compute all viscous terms in invtr4 (implicit)
c
      if(igrad.ne.0) then
      do kc=1,n3m
      do jc=1,n2m
      do ic=1,n1m
      stra=(q(3,ic,jc,kc))*grbar(kc)
      hro(ic,jc,kc)= hro(ic,jc,kc)-stra
      enddo
      enddo
      enddo
                     endif
      return
      end
c
c
c
c  ****************************** subrout mgrbar  **********************
c
c make the gradient of the density  
c make d/dz rhobar
c
c
      subroutine mgrbar
c
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/strat1/istrat,rho0,g,schm   
      common/strat2/igrad,bvais   
      common/strat3/grbar(m3)
      common/stratf/rich
c
c
c grbar=d(rhobar)/dz is vector that is either specified in the file grbar.in
c  or calculated  from the constant brunt-vaisala frequency bvais. 
c  the sign - is since grbar is the stratif. term on the right hand side
c
      if(igrad.eq.1) then
      do 10 kc=1,n3m
      grbar(kc)=-rich
   10 continue
      print*,'grbar(1) (d rhobar/dz)'
      write(6,*) grbar(1)
      endif
      if(igrad.eq.-1) then
      open(17,file='grbar.in')
      do 20 kc=1,n3m
      read(17,*)k,grbar(k)
   20 continue
      close(17)
      endif
      return
      end

c
c   this routine is used to check global conservation
c   properties
      subroutine invti4(rho,drho,al,ga,ro,ru4)
c
      include 'param.f'
      common/strat1/istrat,rho0,g,schm   
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      dimension rho(m1,m2,m3),drho(m1,m2,m3),ru4(m1,m2,m3)
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/tstep/dt,beta,ren
      common/ledat/ics0,cvisc,ifiltr,ibox
c
c  ********* compute the rhs of the factored equation
c  the rhs include h(n), h(n-1)=ru, grad(pr) component
c  everything at i+1/2,j+1/2,k
c  dq=qhat-q(n)
c
      do 100 kc=1,n3m
      km=kmv(kc)
      kp=kpv(kc)
      do 100 jc=1,n2m
      jp=jpv(jc)
      jm=jmv(jc)
      do 100 ic=1,n1m
      im=imv(ic)
      ip=ipv(ic) 
      rhs=(ga*drho(ic,jc,kc)+ro*ru4(ic,jc,kc)
     1             )*dt
      ru4(ic,jc,kc)=drho(ic,jc,kc) 
      rho(ic,jc,kc)=rho(ic,jc,kc)+rhs
 100  continue
      return
      end
c
c  ****************************** subrout invtr4  **********************
c
c   this subroutine performs the inversion of the rho evolution equation
c   by a factored implicit scheme, only the diffusive derivatives 
c   11,22,33 of rho 
c   are treated implicitly
c   the other linear terms are treated explicitly like the nonlinear 
c   terms
c   Thus we take the equation
c   r'=rn+ga*dt*NL+ro*drold+al*dt*L*(rn+r')
c   This is a mixed version of 3rd order runge-kutta and crank-nicklson
c   Here L is the diffusive operator d/dxi(kappa*d/dxi(.))   which
c   is handled implicitly for the sake of stability. 
c   Thus we must solve
c   (1-al*dt*L/2)*(r'-r)=rhs=ga*dt*NL+ro*drold+2*al*dt*L*rn/2
c   Note we have solved in terms of DELTA rho=(r-r'), and that
c   accounts for the extra factor of 2 in the last term on the rhs.
c   Thus we must invert the matrix (1-al*dt*L/2), and this is done in
c   tripvi.  First we rewrite (1-al*dt*L/2)=(1-A1)(1-A2)(1-A3)
c   which is accurate to order dt, where Ai is proportional to 
c   d/dxi(kappa*d/dxi(.)).
c   The matrix Ai is tridaigonal so we must give tripvi 3 vectors
c   corresponding to the diagonal and two off diagonals. These vectors
c   are called ap, ac and am. For more efficiency, these vectors are made
c   into 2 dimensional arrays. They not only carry the matrix structure 
c   apporpriat to the d/dxi(kappa*d/dxi(.)) derivative, but also
c   one additional index corresponging to one aditional spatial dimension.
c   For example, if i=1 correspoinging to the x direction, then we also
c   index the matrices with the y index and loop on the k index. The
c   subroutine tripvi is desiged to vectorize for most efficient operation
c   in this manner.   
      subroutine invtr4(rho,drho,al,ga,ro,ru4)
c
      include 'param.f'
      common/strat1/istrat,rho0,g,schm   
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      common/amcpj1/amj1(m2),acj1(m2),apj1(m2)
      common/trici/ami(m2,m1),aci(m2,m1),api(m2,m1),fi(m2,m1),
     1             fei(m2,m1),qq(m2,m1),ss(m2,m1)
      common/tvpkj/amj(m3,m2),acj(m3,m2),apj(m3,m2),fj(m3,m2),
     1             fej(m3,m2),qj(m3,m2),sj(m3,m2)
      common/tvpjk/amk(m2,m3),ack(m2,m3),apk(m2,m3),fk(m2,m3),
     1             qek(m2,m3),qk(m2,m3),sk(m2,m3)
      common/rhsc/rhs(m1,m2,m3)
      dimension rho(m1,m2,m3),drho(m1,m2,m3),ru4(m1,m2,m3)
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/tstep/dt,beta,ren
      common/smapra/vis(m1,m2,m3)
      common/ledat/ics0,cvisc,ifiltr,ibox
c
c  ********* compute the rhs of the factored equation
c  the rhs include h(n), h(n-1)=ru, grad(pr) component
c  and the 11, 22,3 second derivat. of q1(n)
c  everything at i,j+1/2,k=1/2
c  dq=qhat-q(n)
c
c compute all viscoud terms in invtr4 (implicit)
c******    large eddies simulation term    *********************
c******    assuming nu=schm*kappa  with schm an assigned constant ***********
c    remember that rho is defined in the same location as q3
c
      do 100 kc=1,n3m
      km=kmv(kc)
      kp=kpv(kc)
      do 100 jc=1,n2m
      jp=jpv(jc)
      jm=jmv(jc)
      do 100 ic=1,n1m
      im=imv(ic)
      ip=ipv(ic) 
c
c   11 second derivatives of  rho
c

      visb=(vis(ic,jc,km)+vis(ic,jc,kc)+vis(ip,jc,kc)+vis(ip,jc,km))*.25
      visa=(vis(im,jc,km)+vis(im,jc,kc)+vis(ic,jc,kc)+vis(ic,jc,km))*.25
      drho1=(visb*rho(ip,jc,kc)-(visb+visa)*rho(ic,jc,kc)+
     1      visa*rho(im,jc,kc))*dx1q
c
c   22 second derivatives of rho
c

      visc=(vis(ic,jc,km)+vis(ic,jc,kc)+vis(ic,jm,km)+vis(ic,jm,kc))*.25
      visd=(vis(ic,jc,km)+vis(ic,jc,kc)+vis(ic,jp,km)+vis(ic,jp,kc))*.25
      drho2=(rho(ic,jp,kc)*visd
     1     -rho(ic,jc,kc)*(visd+visc)
     1     +rho(ic,jm,kc)*visc)*dx2q
c
c   33 second derivatives of rho
c
      drho3=(vis(ic,jc,kc)*rho(ic,jc,kp)-(vis(ic,jc,kc)+vis(ic,jc,km))*
     1     rho(ic,jc,kc)+vis(ic,jc,km)*rho(ic,jc,km))*dx3q
      dcr=+(drho1+drho2+drho3)
      rhs(ic,jc,kc)=(ga*drho(ic,jc,kc)+ro*ru4(ic,jc,kc)
     1             +al*dcr)*dt
      ru4(ic,jc,kc)=drho(ic,jc,kc) 
 100  continue
c
c  ********* compute dq1*  sweeping in the x1 direction
c            periodic
      betadx=.5*dx1q*al*dt
      do 1 kc=1,n3m
      do 9 ic=1,n1m
      im=imv(ic)
      ip=ipv(ic) 
      do 9 jc=1,n2m
      visb=(vis(ic,jc,km)+vis(ic,jc,kc)+vis(ip,jc,kc)+vis(ip,jc,km))*.25
      visa=(vis(im,jc,km)+vis(im,jc,kc)+vis(ic,jc,kc)+vis(ic,jc,km))*.25
      api(jc,ic)=-betadx*visa
      aci(jc,ic)=1.+betadx*(visa+visb)
      ami(jc,ic)=-betadx*visb
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
c
c  ************ compute  from drho** sweeping along the x3 direction
c               periodic
      betadx=dt*.5*dx3q*al
      do 6 ic=1,n1m
      do 8 kc=1,n3m
      km=kmv(kc)
      do 8 jc=1,n2m
      apk(jc,kc)=-betadx*vis(ic,jc,kc)
      ack(jc,kc)=1.+betadx*(vis(ic,jc,kc)+vis(ic,jc,km))
      amk(jc,kc)=-betadx*vis(ic,jc,km)
      fk(jc,kc)=rhs(ic,jc,kc)
    8 continue
 
      call trvpjk(1,n3m,1,n2m)
c
      do 5 kc=1,n3m
      do 5 jc=1,n2m
      rhs(ic,jc,kc)=fk(jc,kc)
    5 continue
    6 continue
c
c  ************ compute drho  sweeping along the x2 direction
c
      aldt=al*dt*.5*dx2q
      do 110 ic=1,n1m
      do 120 kc=1,n3m
      km=kmv(kc)
      do 120 jc=1,n2m
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
  120 continue
c
      call trvpkj(1,n2m,1,n3m)
c
      do 30 kc=1,n3m
      do 30 jc=1,n2m
      rho(ic,jc,kc)=rho(ic,jc,kc)+fj(kc,jc)
   30 continue
  110 continue
      return
      end
