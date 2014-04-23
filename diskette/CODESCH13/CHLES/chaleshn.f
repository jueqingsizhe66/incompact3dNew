c  ****************************** subrout hdnl1  **********************
c
c  in this subroutine are calculated the non-linear terms 
c  cross derivatives and terms where appear the other vel. components
c  e.g. in the eqution for q1 appear terms with second derivatives
c  of q2 e.g. LES  . It does not occur in DNS   
c
c
      subroutine hdnl1(q)
c
c
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      common/metria/caj(m2),cac(m2)
      common/qwallo/q1s(m1,m3),q2s(m1,m3),q3s(m1,m3)
      common/qwalup/q1n(m1,m3),q2n(m1,m3),q3n(m1,m3)
      common/rhsc/rhs(m1,m2,m3)
      common/sma/vis(m1,m2,m3)
      common/ledat/csma,cvisc,iles
      dimension q(ndv,m1,m2,m3)
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/islwal/islv1s,islv1n,islv3s,islv3n
c
c   **********  compute the non-linear terms by centered differences
c
c   h term for the q1 momentum equation at i,j+1/2,k=1/2
c
      qdx1=dx1*0.25
      qdx2=dx2*0.25
      qdx3=dx3*0.25
      n2mm=n2m-1
      do kc=1,n3m
      km=kmv(kc)
      kp=kpv(kc)
      do jc=1,n2m
      do ic=1,n1m
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
c   q1 q3 term
c
      h13=((q(3,ic,jc,kp)+q(3,im,jc,kp))
     1    *(q(1,ic,jc,kp)+q(1,ic,jc,kc))
     3    -(q(3,ic,jc,kc)+q(3,im,jc,kc))
     1    *(q(1,ic,jc,kc)+q(1,ic,jc,km))  )*qdx3
      hq1=h11+h13
      rhs(ic,jc,kc)=-hq1
      enddo
      enddo
      enddo
      do kc=1,n3m
      do jc=2,n2mm
      jmm=jmv(jc)
      jpp=jpv(jc)
      jp=jc+1
      uqdx2=qdx2/caj(jc)
      do ic=1,n1m
      ip=ipv(ic)
      im=imv(ic)
c
c   q1 q2 term  inner field
c
      h12=( (q(2,ic,jp,kc)+q(2,im,jp,kc))*
     1      (q(1,ic,jpp,kc)+q(1,ic,jc,kc))
     1    - (q(2,ic,jc,kc)+q(2,im,jc,kc))*
     1      (q(1,ic,jc,kc)+q(1,ic,jmm,kc)) )*uqdx2
      hq1=h12
      rhs(ic,jc,kc)=rhs(ic,jc,kc)-hq1
      enddo
      enddo
      enddo
c
c   lower wall
c
      jc=1
      jp=jc+1
      jpp=jpv(jc)
      jmm=jmv(jc)
      uqdx2=dx2/caj(jc)
      if(islv1s.eq.0) then
c
c     q1 q2  q1 no-slip
c
      do kc=1,n3m
      do ic=1,n1m
      ip=ipv(ic)
      im=imv(ic)
      h12=( (q(2,ic,jp,kc)+q(2,im,jp,kc))*
     1      (q(1,ic,jp,kc)+q(1,ic,jc,kc))*.25
     1    - (q(2,ic,jc,kc)+q(2,im,jc,kc))*.5*
     1      q1s(ic,kc)   ) *uqdx2
c
      hq1=h12
      rhs(ic,jc,kc)=rhs(ic,jc,kc)-hq1
      enddo
      enddo
                       else
c
c     q1 q2 q1   indices of jc account for shear-free 
c
      do kc=1,n3m
      uqdx2=dx2/caj(jc)
      do ic=1,n1m
      im=imv(ic)
      h12=( (q(2,ic,jp,kc)+q(2,im,jp,kc))*
     1      (q(1,ic,jpp,kc)+q(1,ic,jc,kc))
     1    - (q(2,ic,jc,kc)+q(2,im,jc,kc))*
     1      (q(1,ic,jc,kc)+q(1,ic,jmm,kc)) )*uqdx2
c
      hq1=h12
      rhs(ic,jc,kc)=rhs(ic,jc,kc)-hq1
      enddo
      enddo
                       endif
c
c   upper wall
c
      jc=n2m
      jp=jc+1
      jm=jc-1
      jpp=jpv(jc)
      jmm=jmv(jc)
      uqdx2=dx2/caj(jc)
      if(islv1n.eq.0) then
c
c     q1 q2 q1 no-slip
c
      do kc=1,n3m
      do ic=1,n1m
      im=imv(ic)
      h12=( (q(2,ic,jp,kc)+q(2,im,jp,kc))*0.5*
     1      q1n(ic,kc)   
     1    - (q(2,ic,jc,kc)+q(2,im,jc,kc))*
     1      (q(1,ic,jm,kc)+q(1,ic,jc,kc))*.25   ) *uqdx2
c
      hq1=h12
      rhs(ic,jc,kc)=rhs(ic,jc,kc)-hq1
      enddo
      enddo
                       else
c
c     q1 q2 q1   indices of jc account for shear-free
c
      do kc=1,n3m
      do ic=1,n1m
      im=imv(ic)
      h12=( (q(2,ic,jp,kc)+q(2,im,jp,kc))*
     1      (q(1,ic,jpp,kc)+q(1,ic,jc,kc))
     1    - (q(2,ic,jc,kc)+q(2,im,jc,kc))*
     1      (q(1,ic,jc,kc)+q(1,ic,jmm,kc)) )*uqdx2
      hq1=h12
      rhs(ic,jc,kc)=rhs(ic,jc,kc)-hq1
      enddo
      enddo
                        endif
      if(iles.ne.0)  then
c     
c******    large eddies simulation term    *********************
c    inner field
c
      do kc=1,n3m
      km=kmv(kc)
      kp=kpv(kc)
      do jc=2,n2m-1
      jpp=jc+1
      jp=jpv(jc)
      jm=jmv(jc)
      do ic=1,n1m
      im=imv(ic)
      ip=ipv(ic)
      visb=.25*(vis(ic,jc,kc)+vis(ic,jp,kc)+vis(im,jc,kc)+vis(im,jp,kc))
      if(jc.eq.n2m) visb=cvisc
      visa=.25*(vis(ic,jc,kc)+vis(ic,jm,kc)+vis(im,jc,kc)+vis(im,jm,kc))
      if(jc.eq.1) visa=cvisc
      visd=.25*(vis(ic,jc,kc)+vis(ic,jc,kp)+vis(im,jc,kc)+vis(im,jc,kp))
      visc=.25*(vis(ic,jc,kc)+vis(ic,jc,km)+vis(im,jc,kc)+vis(im,jc,km))
      sgs1=dx1q*(vis(ic,jc,kc)*(q(1,ip,jc,kc)-q(1,ic,jc,kc))-
     1           vis(im,jc,kc)*(q(1,ic,jc,kc)-q(1,im,jc,kc)))
      sgs2=dx1*dx2/caj(jc)*(
     1    visb*(q(2,ic,jpp,kc)-q(2,im,jpp,kc))-
     1    visa*(q(2,ic,jc,kc)-q(2,im,jc,kc)))
      sgs3=dx3*dx1*(
     1     visd*(q(3,ic,jc,kp)-q(3,im,jc,kp))-
     1     visc*(q(3,ic,jc,kc)-q(3,im,jc,kc)))
      rhs(ic,jc,kc)=rhs(ic,jc,kc)+sgs1+sgs2+sgs3
      enddo
      enddo
      enddo
c     
c******    large eddies simulation term    *********************
c    lower wall  
c
      jc=1
      jpp=jc+1
      jp=jpv(jc)
      jm=jmv(jc)
      do kc=1,n3m
      km=kmv(kc)
      kp=kpv(kc)
      do ic=1,n1m
      im=imv(ic)
      ip=ipv(ic)
      visb=.25*(vis(ic,jc,kc)+vis(ic,jp,kc)+vis(im,jc,kc)+vis(im,jp,kc))
      visa=cvisc
      visd=.25*(vis(ic,jc,kc)+vis(ic,jc,kp)+vis(im,jc,kc)+vis(im,jc,kp))
      visc=.25*(vis(ic,jc,kc)+vis(ic,jc,km)+vis(im,jc,kc)+vis(im,jc,km))
      sgs1=dx1q*(vis(ic,jc,kc)*(q(1,ip,jc,kc)-q(1,ic,jc,kc))-
     1           vis(im,jc,kc)*(q(1,ic,jc,kc)-q(1,im,jc,kc)))
      sgs2=dx1*dx2/caj(jc)*(
     1    visb*(q(2,ic,jpp,kc)-q(2,im,jpp,kc))-
     1    visa*(q(2,ic,jc,kc)-q(2,im,jc,kc)))
      sgs3=dx3*dx1*(
     1     visd*(q(3,ic,jc,kp)-q(3,im,jc,kp))-
     1     visc*(q(3,ic,jc,kc)-q(3,im,jc,kc)))
      rhs(ic,jc,kc)=rhs(ic,jc,kc)+sgs1+sgs2+sgs3
      enddo
      enddo
c     
c******    large eddies simulation term    *********************
c    upper wall 
c
      jc=n2m
      jpp=jc+1
      jp=jpv(jc)
      jm=jmv(jc)
      do kc=1,n3m
      km=kmv(kc)
      kp=kpv(kc)
      do ic=1,n1m
      im=imv(ic)
      ip=ipv(ic)
      visb=cvisc
      visa=.25*(vis(ic,jc,kc)+vis(ic,jm,kc)+vis(im,jc,kc)+vis(im,jm,kc))
      visd=.25*(vis(ic,jc,kc)+vis(ic,jc,kp)+vis(im,jc,kc)+vis(im,jc,kp))
      visc=.25*(vis(ic,jc,kc)+vis(ic,jc,km)+vis(im,jc,kc)+vis(im,jc,km))
      sgs1=dx1q*(vis(ic,jc,kc)*(q(1,ip,jc,kc)-q(1,ic,jc,kc))-
     1           vis(im,jc,kc)*(q(1,ic,jc,kc)-q(1,im,jc,kc)))
      sgs2=dx1*dx2/caj(jc)*(
     1    visb*(q(2,ic,jpp,kc)-q(2,im,jpp,kc))-
     1    visa*(q(2,ic,jc,kc)-q(2,im,jc,kc)))
      sgs3=dx3*dx1*(
     1     visd*(q(3,ic,jc,kp)-q(3,im,jc,kp))-
     1     visc*(q(3,ic,jc,kc)-q(3,im,jc,kc)))
      rhs(ic,jc,kc)=rhs(ic,jc,kc)+sgs1+sgs2+sgs3
      enddo
      enddo
                      endif
      return
      end
c
c  ****************************** subrout hdnl2  **********************
c
c  in this subroutine are calculated the non-linear terms,
c  cross derivatives and terms where appear the other vel. components
c  e.g. in the eqution for q1 appear terms with second derivatives
c  of q2 e.g. LES  . It does not occur in DNS   
c
      subroutine hdnl2(q,h)
c
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      dimension h(m1,m2,m3)
      dimension q(ndv,m1,m2,m3)
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/metria/caj(m2),cac(m2)
      common/sma/vis(m1,m2,m3)
      common/ledat/csma,cvisc,iles
c
c  **********  compute the non-linear terms by centered differences
c
c   h term for the q2 momentum equation at i+1/2,j,k=1/2
c
      qdx1=dx1*0.25
      qdx2=dx2*0.25
      qdx3=dx3*0.25
      do kc=1,n3m
      km=kmv(kc)
      kp=kpv(kc)
      do jc=2,n2m
      jm=jc-1
      jp=jc+1
      sucac=1./cac(jc)
      do ic=1,n1m
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
     1     (q(2,ic,jc,kc)+q(2,ic,jm,kc)) )*qdx2*sucac
c
c     q2 q3 term
c
      h23=((q(3,ic,jc,kp)+q(3,ic,jm,kp))
     1    *(q(2,ic,jc,kp)+q(2,ic,jc,kc))
     3    -(q(3,ic,jc,kc)+q(3,ic,jm,kc))
     1    *(q(2,ic,jc,kc)+q(2,ic,jc,km))
     1    )*qdx3
      hq2=h21+h22+h23
      h(ic,jc,kc)=-hq2
      enddo
      enddo
      enddo
      if(iles.ne.0)  then
c     
c******    large eddies simulation term    *********************
c    inner field
c
      do kc=1,n3m
      km=kmv(kc)
      kp=kpv(kc)
      do jc=2,n2m
      jp=jc+1
      jm=jc-1
      suc1=1./cac(jc)
      do ic=1,n1m
      im=imv(ic)
      ip=ipv(ic)
      visb=.25*(vis(ic,jc,kc)+vis(ic,jm,kc)+vis(ip,jc,kc)+vis(ip,jm,kc))
      visa=.25*(vis(ic,jc,kc)+vis(ic,jm,kc)+vis(im,jc,kc)+vis(im,jm,kc))
      visd=.25*(vis(ic,jc,kc)+vis(ic,jm,kc)+vis(ic,jc,kp)+vis(ic,jm,kp))
      visc=.25*(vis(ic,jc,kc)+vis(ic,jm,kc)+vis(ic,jc,km)+vis(ic,jm,km))
      sgs1=dx1*dx2*(
     1     visb*(q(1,ip,jc,kc)-q(1,ip,jm,kc))-
     1     visa*(q(1,ic,jc,kc)-q(1,ic,jm,kc)))
      sgs2=dx2q*(
     1     vis(ic,jc,kc)/caj(jc)*(q(2,ic,jp,kc)-q(2,ic,jc,kc))-
     1     vis(ic,jm,kc)/caj(jm)*(q(2,ic,jc,kc)-q(2,ic,jm,kc)))
      sgs3=dx3*dx2*(
     1     visd*(q(3,ic,jc,kp)-q(3,ic,jm,kp))-
     1     visc*(q(3,ic,jc,kc)-q(3,ic,jm,kc)))
      sgs=suc1*(sgs1+sgs2+sgs3)
      h(ic,jc,kc)=h(ic,jc,kc)+sgs
      enddo
      enddo
      enddo
               endif
      return
      end
c
c  ****************************** subrout hdnl3  **********************
c
c  in this subroutine are calculated the non-linear terms,
c  cross derivatives and terms where appear the other vel. components
c  e.g. in the eqution for q1 appear terms with second derivatives
c  of q2 e.g. in LES . It does not occur in DNS
c
      subroutine hdnl3(q,h)
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      common/qwallo/q1s(m1,m3),q2s(m1,m3),q3s(m1,m3)
      common/qwalup/q1n(m1,m3),q2n(m1,m3),q3n(m1,m3)
      dimension h(m1,m2,m3)
      dimension q(ndv,m1,m2,m3)
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/metria/caj(m2),cac(m2)
      common/sma/vis(m1,m2,m3)
      common/ledat/csma,cvisc,iles
      common/islwal/islv1s,islv1n,islv3s,islv3n
c
c
      qdx1=dx1*0.25
      qdx2=dx2*0.25
      qdx3=dx3*0.25
      n2mm=n2m-1
      do kc=1,n3m
      km=kmv(kc)
      kp=kpv(kc)
      do jc=1,n2m
      sucaj=1./caj(jc)
      do ic=1,n1m
      ip=ipv(ic)
      im=imv(ic)
c
c    q3 q1 term
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
      hq3=h31+h33
      h(ic,jc,kc)=-hq3
      enddo
      enddo
      enddo
      do kc=1,n3m
      km=kmv(kc)
      do jc=2,n2mm
      jp=jc+1
      jmm=jmv(jc)
      jpp=jpv(jc)
      sucaj=1./caj(jc)
      do ic=1,n1m
      ip=ipv(ic)
      im=imv(ic)
c
c    q3 q2 term
c
      h32=( (q(2,ic,jp,kc)+q(2,ic,jp,km))*
     1      (q(3,ic,jpp,kc)+q(3,ic,jc,kc))
     1     -(q(2,ic,jc,kc)+q(2,ic,jc,km))*
     1      (q(3,ic,jc,kc)+q(3,ic,jmm,kc))
     1    )*qdx2*sucaj
      hq3=h32
      h(ic,jc,kc)=h(ic,jc,kc)-hq3
      enddo
      enddo
      enddo
c    lower wall
      jc = 1
      jp=jc+1
      jmm=jmv(jc)
      jpp=jpv(jc)
      if(islv3s.eq.0)  then
      sucaj=1./caj(jc)
         do kc=1,n3m
         km=kmv(kc)
         do ic=1,n1m
c   q3 no-slip
          h32=( (q(2,ic,jp,kc)+q(2,ic,jp,km))*
     1            (q(3,ic,jp,kc)+q(3,ic,jc,kc))*0.25
     1           -(q(2,ic,jc,kc)+q(2,ic,jc,km))*0.5*
     1             q3s(ic,kc)
     1          )*dx2*sucaj
           hq3=h32
      h(ic,jc,kc)=h(ic,jc,kc)-hq3
      enddo
      enddo
                        else
         do kc=1,n3m
         km=kmv(kc)
         do ic=1,n1m
c   q3 indices in jc account for shear free
c
      h32=( (q(2,ic,jp,kc)+q(2,ic,jp,km))*
     1      (q(3,ic,jpp,kc)+q(3,ic,jc,kc))
     1     -(q(2,ic,jc,kc)+q(2,ic,jc,km))*
     1      (q(3,ic,jc,kc)+q(3,ic,jmm,kc))
     1    )*qdx2*sucaj
           hq3=h32
      h(ic,jc,kc)=h(ic,jc,kc)-hq3
      enddo
      enddo
                        endif
c    upper wall
      jc = n2m
      jp=jc+1
      jm=jc-1
      jmm=jmv(jc)
      jpp=jpv(jc)
      sucaj=1./caj(jc)
      if(islv3n.eq.0)  then
         do kc=1,n3m
         km=kmv(kc)
         do ic=1,n1m
c   q3 no-slip
c
          h32=( (q(2,ic,jp,kc)+q(2,ic,jp,km))*0.5*
     1             q3n(ic,kc)
     1           -(q(2,ic,jc,kc)+q(2,ic,jc,km))*
     1            (q(3,ic,jm,kc)+q(3,ic,jc,kc))*0.25
     1          )*dx2*sucaj
           hq3=h32
      h(ic,jc,kc)=h(ic,jc,kc)-hq3
      enddo
      enddo
                        else
         do kc=1,n3m
         km=kmv(kc)
         do ic=1,n1m
c   q3 indices in jc account for shear free
c
      h32=( (q(2,ic,jp,kc)+q(2,ic,jp,km))*
     1      (q(3,ic,jpp,kc)+q(3,ic,jc,kc))
     1     -(q(2,ic,jc,kc)+q(2,ic,jc,km))*
     1      (q(3,ic,jc,kc)+q(3,ic,jmm,kc))
     1    )*qdx2*sucaj
           hq3=h32
      h(ic,jc,kc)=h(ic,jc,kc)-hq3
      enddo
      enddo
                        endif
      if(iles.ne.0)  then
c     
c******    large eddies simulation term    *********************
c    inner field
c
      do kc=1,n3m
      km=kmv(kc)
      kp=kpv(kc)
      do jc=2,n2m-1
      jpp=jc+1
      jp=jpv(jc)
      jm=jmv(jc)
      do ic=1,n1m
      im=imv(ic)
      ip=ipv(ic)
      visb=.25*(vis(ic,jc,km)+vis(ic,jc,kc)+vis(ip,jc,km)+vis(ip,jc,kc))
      visa=.25*(vis(ic,jc,km)+vis(ic,jc,kc)+vis(im,jc,km)+vis(im,jc,kc))
      visd=.25*(vis(ic,jc,km)+vis(ic,jp,kc)+vis(ic,jp,km)+vis(ic,jc,kc))
      visc=.25*(vis(ic,jc,km)+vis(ic,jm,kc)+vis(ic,jm,km)+vis(ic,jc,kc))
      sgs1=dx1*dx3*(
     1     visb*(q(1,ip,jc,kc)-q(1,ip,jc,km))-
     1     visa*(q(1,ic,jc,kc)-q(1,ic,jc,km)))
      sgs2=dx2*dx3/caj(jc)*(
     1     visd*(q(2,ic,jpp,kc)-q(2,ic,jpp,km))-
     1     visc*(q(2,ic,jc,kc)-q(2,ic,jc,km)))
      sgs3=dx3q*(
     1     vis(ic,jc,kc)*(q(3,ic,jc,kp)-q(3,ic,jc,kc))-
     1     vis(ic,jc,km)*(q(3,ic,jc,kc)-q(3,ic,jc,km)))
      h(ic,jc,kc)=h(ic,jc,kc)+sgs1+sgs2+sgs3
      enddo
      enddo
      enddo
c     
c******    large eddies simulation term    *********************
c    lower wall 
c
      jc=1
      jpp=jc+1
      jp=jpv(jc)
      jm=jmv(jc)
      do kc=1,n3m
      km=kmv(kc)
      kp=kpv(kc)
      do ic=1,n1m
      im=imv(ic)
      ip=ipv(ic)
      visb=.25*(vis(ic,jc,km)+vis(ic,jc,kc)+vis(ip,jc,km)+vis(ip,jc,kc))
      visa=.25*(vis(ic,jc,km)+vis(ic,jc,kc)+vis(im,jc,km)+vis(im,jc,kc))
      visd=.25*(vis(ic,jc,km)+vis(ic,jp,kc)+vis(ic,jp,km)+vis(ic,jc,kc))
      visc=cvisc
      sgs1=dx1*dx3*(
     1     visb*(q(1,ip,jc,kc)-q(1,ip,jc,km))-
     1     visa*(q(1,ic,jc,kc)-q(1,ic,jc,km)))
      sgs2=dx2*dx3/caj(jc)*(
     1     visd*(q(2,ic,jpp,kc)-q(2,ic,jpp,km))-
     1     visc*(q(2,ic,jc,kc)-q(2,ic,jc,km)))
      sgs3=dx3q*(
     1     vis(ic,jc,kc)*(q(3,ic,jc,kp)-q(3,ic,jc,kc))-
     1     vis(ic,jc,km)*(q(3,ic,jc,kc)-q(3,ic,jc,km)))
      h(ic,jc,kc)=h(ic,jc,kc)+sgs1+sgs2+sgs3
      enddo
      enddo
c     
c******    large eddies simulation term    *********************
c    upper wall 
c
      jc=n2m
      jpp=jc+1
      jp=jpv(jc)
      jm=jmv(jc)
      do kc=1,n3m
      km=kmv(kc)
      kp=kpv(kc)
      do ic=1,n1m
      im=imv(ic)
      ip=ipv(ic)
      visb=.25*(vis(ic,jc,km)+vis(ic,jc,kc)+vis(ip,jc,km)+vis(ip,jc,kc))
      visa=.25*(vis(ic,jc,km)+vis(ic,jc,kc)+vis(im,jc,km)+vis(im,jc,kc))
      visd=cvisc
      visc=.25*(vis(ic,jc,km)+vis(ic,jm,kc)+vis(ic,jm,km)+vis(ic,jc,kc))
      sgs1=dx1*dx3*(
     1     visb*(q(1,ip,jc,kc)-q(1,ip,jc,km))-
     1     visa*(q(1,ic,jc,kc)-q(1,ic,jc,km)))
      sgs2=dx2*dx3/caj(jc)*(
     1     visd*(q(2,ic,jpp,kc)-q(2,ic,jpp,km))-
     1     visc*(q(2,ic,jc,kc)-q(2,ic,jc,km)))
      sgs3=dx3q*(
     1     vis(ic,jc,kc)*(q(3,ic,jc,kp)-q(3,ic,jc,kc))-
     1     vis(ic,jc,km)*(q(3,ic,jc,kc)-q(3,ic,jc,km)))
      h(ic,jc,kc)=h(ic,jc,kc)+sgs1+sgs2+sgs3
      enddo
      enddo
               endif

      return
      end
