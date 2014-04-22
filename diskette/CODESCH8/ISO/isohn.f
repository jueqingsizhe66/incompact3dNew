c  ****************************** subrout hdnl1  **********************
c
c  in this subroutine are calculated the non-linear terms 
c  with staggered velocities the schemes conserves energy
c
      subroutine hdnl1(q)
c
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      dimension h(m1,m2,m3)
      dimension q(ndv,m1,m2,m3)
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/rhsc/rhs(m1,m2,m3)
      common/ledat/cvisc
      common/rot/f0
c
c   **********  compute the non-linear terms by centered differences
c   bilinear interpolation is used 
c   non-linear terms
c
c   h term for the q1 momentum equation at i,j+1/2,k+1/2
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
c  add rotation in the direction x3to rhs:
c   -fxu=-f(0,0,1)x(u1,u2,u3)=f(u2,-u1,0) 
c
      rhs(ic,jc,kc)=rhs(ic,jc,kc)+
     & 0.25*f0*(q(2,ic,jc,kc)+q(2,ic,jp,kc)+q(2,im,jc,kc)+q(2,im,jp,kc))
   10 continue
c************
      endif
      return
      end
c
c  ****************************** subrout hdnl2  **********************
c
c  in this subroutine are calculated the non-linear terms,
c  for the q2 velocity component
c
      subroutine hdnl2(q,h)
c
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      dimension h(m1,m2,m3)
      dimension q(ndv,m1,m2,m3)
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/ledat/cvisc
      common/rot/f0
c
c  **********  compute the non-linear terms by centered differences
c   bilinear interpolations are used 
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
c  add rotation in the x3 direction to rhs: 
c  -fxu=-f(0,0,1)x(u1,u2,u3)=f(u2,-u1,0) 
c
      h(ic,jc,kc)=h(ic,jc,kc)-
     & 0.25*f0*(q(1,ic,jc,kc)+q(1,ip,jc,kc)+q(1,ic,jm,kc)+q(1,ip,jm,kc))
   20 continue
c*********************** 
      return
      end
c
c  ****************************** subrout hdnl3  **********************
c
c  in this subroutine are calculated the non-linear terms,
c  of the q3 velocity component
c
      subroutine hdnl3(q,h)
c
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      dimension h(m1,m2,m3)
      dimension q(ndv,m1,m2,m3)
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/ledat/cvisc
c
c
c  **********  compute the non-linear terms by centered differences
c   bilinear interpolations are used 
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
      return
      end
