c  ****************************** subrout hdnl1  **********************
c
c  in this subroutine are calculated the non-linear terms and sub-grid terms.
c
      subroutine hdnl1(q)
c
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      dimension q(ndv,m1,m2,m3)
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/sma/vis(m1,m2,m3)
      common/ledat/ics0,cvisc,ifiltr,ibox
      common/rhsc/rhs(m1,m2,m3)
c     common/ss1/sum1,sue1
      common/rot/f0
c
c
c
c   **********  compute the non-linear terms by centered differences
c
c   h term for the q1 momentum equation at i,j+1/2,k=1/2
c
      qdx1=dx1*0.25
      qdx2=dx2*0.25
      qdx3=dx3*0.25
      if(n1m.gt.1) then
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
c    q1 q1 term
c
      h11=( (q(1,ip,jc,kc)+q(1,ic,jc,kc))*
     1      (q(1,ip,jc,kc)+q(1,ic,jc,kc))
     1    - (q(1,im,jc,kc)+q(1,ic,jc,kc))*
     1      (q(1,im,jc,kc)+q(1,ic,jc,kc)) )*qdx1
c
c   q1 q2 term
c
      h12=( (q(2,ic,jp,kc)+q(2,im,jp,kc))*
     1      (q(1,ic,jp,kc)+q(1,ic,jc,kc))
     1    - (q(2,ic,jc,kc)+q(2,im,jc,kc))*
     1      (q(1,ic,jc,kc)+q(1,ic,jm,kc)) )*qdx2
c
c   q1 q3 term
c
      h13=((q(3,ic,jc,kp)+q(3,im,jc,kp))
     1    *(q(1,ic,jc,kp)+q(1,ic,jc,kc))
     3    -(q(3,ic,jc,kc)+q(3,im,jc,kc))
     1    *(q(1,ic,jc,kc)+q(1,ic,jc,km))  )*qdx3
      hq1=h11+h12+h13
      rhs(ic,jc,kc)=-hq1
c
c    add rotation 
c   -fxu=-f(0,0,1)x(u1,u2,u3)=f(u2,-u1,0) 
c   to rhs of the momentum equations
c
      rhs(ic,jc,kc)=rhs(ic,jc,kc)+
     & 0.25*f0*(q(2,ic,jc,kc)+q(2,ic,jp,kc)+q(2,im,jc,kc)+q(2,im,jp,kc))
   10 continue
      if(ics0.ne.0) then
c
c******    large eddies simulation term    *********************
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
      visb=.25*(vis(ic,jc,kc)+vis(ic,jp,kc)+vis(im,jc,kc)+vis(im,jp,kc))
      visa=.25*(vis(ic,jc,kc)+vis(ic,jm,kc)+vis(im,jc,kc)+vis(im,jm,kc))
      visd=.25*(vis(ic,jc,kc)+vis(ic,jc,kp)+vis(im,jc,kc)+vis(im,jc,kp))
      visc=.25*(vis(ic,jc,kc)+vis(ic,jc,km)+vis(im,jc,kc)+vis(im,jc,km))
      sgs1=dx1q*(vis(ic,jc,kc)*(q(1,ip,jc,kc)-q(1,ic,jc,kc))-
     1           vis(im,jc,kc)*(q(1,ic,jc,kc)-q(1,im,jc,kc)))
      sgs2=dx1*dx2*(
     1    visb*(q(2,ic,jp,kc)-q(2,im,jp,kc))-
     1    visa*(q(2,ic,jc,kc)-q(2,im,jc,kc)))
      sgs3=dx3*dx1*(
     1     visd*(q(3,ic,jc,kp)-q(3,im,jc,kp))- 
     1     visc*(q(3,ic,jc,kc)-q(3,im,jc,kc))) 
      rhs(ic,jc,kc)=rhs(ic,jc,kc)+sgs1+sgs2+sgs3
 100  continue
      end if
      endif
      return
      end
c
c  ****************************** subrout hdnl2  **********************
c
c  in this subroutine are calculated the non-linear terms and sub-grid terms.
c
c
      subroutine hdnl2(q,h)
c
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      dimension h(m1,m2,m3)
      dimension q(ndv,m1,m2,m3)
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/sma/vis(m1,m2,m3)
      common/ledat/ics0,cvisc,ifiltr,ibox
      common/rot/f0
c
c  **********  compute the non-linear terms by centered differences
c
c   h term for the q2 momentum equation at i+1/2,j,k+1/2
c
      qdx1=dx1*0.25
      qdx2=dx2*0.25
      qdx3=dx3*0.25
      do 20 kc=1,n3m
      km=kmv(kc)
      kp=kpv(kc)
      do 20 jc=1,n2m
      jm=jmv(jc)
      jp=jpv(jc)
      do 20 ic=1,n1m
      ip=ipv(ic)
      im=imv(ic)
c
c     q2 q1 term
c
      h21=((q(1,ip,jc,kc)+q(1,ip,jm,kc))*
     1     (q(2,ip,jc,kc)+q(2,ic,jc,kc))
     1    -(q(1,ic,jc,kc)+q(1,ic,jm,kc))*
     1     (q(2,ic,jc,kc)+q(2,im,jc,kc)) )*qdx1
c
c     q2 q2 term
c
      h22=((q(2,ic,jp,kc)+q(2,ic,jc,kc))*
     1     (q(2,ic,jp,kc)+q(2,ic,jc,kc))
     1    -(q(2,ic,jc,kc)+q(2,ic,jm,kc))*
     1     (q(2,ic,jc,kc)+q(2,ic,jm,kc)) )*qdx2
c
c     q2 q3 term
c
      h23=((q(3,ic,jc,kp)+q(3,ic,jm,kp))
     1    *(q(2,ic,jc,kp)+q(2,ic,jc,kc))
     3    -(q(3,ic,jc,kc)+q(3,ic,jm,kc))
     1    *(q(2,ic,jc,kc)+q(2,ic,jc,km)) )*qdx3
      hq2=h21+h22+h23
      h(ic,jc,kc)=-hq2
c
c  add rotation -fxu=-f(0,0,1)x(u1,u2,u3)=f(u2,-u1,0)  
c  to the momentum equations
c
      h(ic,jc,kc)=h(ic,jc,kc)-
     & 0.25*f0*(q(1,ic,jc,kc)+q(1,ip,jc,kc)+q(1,ic,jm,kc)+q(1,ip,jm,kc))
   20 continue
C
c******    large eddies simulation term    *********************
c
      if(ics0.ne.0) then
      do 100 kc=1,n3m
      km=kmv(kc)
      kp=kpv(kc)
      do 100 jc=1,n2m
      jp=jpv(jc)
      jm=jmv(jc)
      do 100 ic=1,n1m
      im=imv(ic)
      ip=ipv(ic) 
      visb=.25*(vis(ic,jc,kc)+vis(ic,jm,kc)+vis(ip,jc,kc)+vis(ip,jm,kc))
      visa=.25*(vis(ic,jc,kc)+vis(ic,jm,kc)+vis(im,jc,kc)+vis(im,jm,kc))
      visd=.25*(vis(ic,jc,kc)+vis(ic,jm,kc)+vis(ic,jc,kp)+vis(ic,jm,kp))
      visc=.25*(vis(ic,jc,kc)+vis(ic,jm,kc)+vis(ic,jc,km)+vis(ic,jm,km))
      sgs1=dx1*dx2*(
     1     visb*(q(1,ip,jc,kc)-q(1,ip,jm,kc))-
     1     visa*(q(1,ic,jc,kc)-q(1,ic,jm,kc)))
      sgs2=dx2q*(vis(ic,jc,kc)*(q(2,ic,jp,kc)-q(2,ic,jc,kc))
     1          -vis(ic,jm,kc)*(q(2,ic,jc,kc)-q(2,ic,jm,kc)))
      sgs3=dx3*dx2*(
     1     visd*(q(3,ic,jc,kp)-q(3,ic,jm,kp))-
     1     visc*(q(3,ic,jc,kc)-q(3,ic,jm,kc)))
      h(ic,jc,kc)=h(ic,jc,kc)+sgs1+sgs2+sgs3
 100  continue
      end if
      return
      end
c
c  ****************************** subrout hdnl3  **********************
c
c  in this subroutine are calculated the non-linear terms and sub-grid terms.
c
c
      subroutine hdnl3(q,h,rho)
c
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      dimension h(m1,m2,m3)
      dimension q(ndv,m1,m2,m3)
      dimension rho(m1,m2,m3)
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/sma/vis(m1,m2,m3)
      common/ledat/ics0,cvisc,ifiltr,ibox
      common/strat1/istrat,rho0,g,schm
      common/strat2/igrad,bvais
c
c
c  **********  compute the non-linear terms by centered differences
c   h term for the q3 momentum equation at i+1/2,j+1/2,k
c
      qdx1=dx1*0.25
      qdx2=dx2*0.25
      qdx3=dx3*0.25
      do 30 kc=1,n3m
      km=kmv(kc)
      kp=kpv(kc)
      do 30 jc=1,n2m
      jm=jmv(jc)
      jp=jpv(jc)
      do 30 ic=1,n1m
      ip=ipv(ic)
      im=imv(ic)
c
c    q3 q1 term
c
      h32=( (q(2,ic,jp,kc)+q(2,ic,jp,km))*
     1      (q(3,ic,jp,kc)+q(3,ic,jc,kc))
     1     -(q(2,ic,jc,kc)+q(2,ic,jc,km))*
     1      (q(3,ic,jc,kc)+q(3,ic,jm,kc))
     1    )*qdx2
c
c    q3 q2 term
c
      h31=((q(1,ip,jc,kc)+q(1,ip,jc,km))
     1    *(q(3,ip,jc,kc)+q(3,ic,jc,kc))
     3    -(q(1,ic,jc,kc)+q(1,ic,jc,km))
     1    *(q(3,ic,jc,kc)+q(3,im,jc,kc))
     1    )*qdx1
c
c    q3 q3 term
c
      h33=((q(3,ic,jc,kp)+q(3,ic,jc,kc))
     1    *(q(3,ic,jc,kp)+q(3,ic,jc,kc))
     3    -(q(3,ic,jc,kc)+q(3,ic,jc,km))
     1    *(q(3,ic,jc,kc)+q(3,ic,jc,km))
     1    )*qdx3
      hq3=h31+h32+h33
      h(ic,jc,kc)=-hq3
   30 continue
      if(istrat.eq.1.and.igrad.ne.0) then 
c
c add stratification term to rhs of z momentum eq
c du3/dt= ... -g*rho'/rho0 the dimensionalization
c gives du3/dt= ... -rho
c
      do kc=1,n3m
      do jc=1,n2m
      do ic=1,n1m
      buoya=-rho(ic,jc,kc)
      h(ic,jc,kc)=h(ic,jc,kc)+buoya
      enddo
      enddo
      enddo
                      endif
c
c******    large eddies simulation term    *********************
c
      if (ics0.ne.0) then
      do 100 kc=1,n3m
      km=kmv(kc)
      kp=kpv(kc)
      do 100 jc=1,n2m
      jp=jpv(jc)
      jm=jmv(jc)
      do 100 ic=1,n1m
      im=imv(ic)
      ip=ipv(ic) 
      visb=.25*(vis(ic,jc,km)+vis(ic,jc,kc)+vis(ip,jc,km)+vis(ip,jc,kc))
      visa=.25*(vis(ic,jc,km)+vis(ic,jc,kc)+vis(im,jc,km)+vis(im,jc,kc))
      visd=.25*(vis(ic,jc,km)+vis(ic,jp,kc)+vis(ic,jp,km)+vis(ic,jc,kc))
      visc=.25*(vis(ic,jc,km)+vis(ic,jm,kc)+vis(ic,jm,km)+vis(ic,jc,kc))
      sgs1=dx1*dx3*(
     1     visb*(q(1,ip,jc,kc)-q(1,ip,jc,km))-
     1     visa*(q(1,ic,jc,kc)-q(1,ic,jc,km)))
      sgs2=dx2*dx3*(
     1     visd*(q(2,ic,jp,kc)-q(2,ic,jp,km))-
     1     visc*(q(2,ic,jc,kc)-q(2,ic,jc,km)))
      sgs3=dx3q*(vis(ic,jc,kc)*q(3,ic,jc,kp)-
     1     (vis(ic,jc,kc)+vis(ic,jc,km))*q(3,ic,jc,kc)+
     1     vis(ic,jc,km)*q(3,ic,jc,km))
      h(ic,jc,kc)=h(ic,jc,kc)+sgs1+sgs2+sgs3
 100  continue
      end if
      return
      end
