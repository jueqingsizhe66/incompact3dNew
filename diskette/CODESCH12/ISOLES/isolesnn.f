c
c
c  ****************************** subrout cfl  **********************
c
c  in this subroutine is calculated the cfl
c
      subroutine cfl(q,cflm)
c
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      common/tstep/dt,beta,ren
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      dimension q(ndv,m1,m2,m3)
c
      cflm=0.
      cflma=0.
      do 7 k=1,n3m
      kp=kpv(k)
      do 7 j=1,n2m
      jp=jpv(j)
      do 7 i=1,n1m
      ip=ipv(i)
      qcf=(abs((q(1,i,j,k)+q(1,ip,j,k))*dx1)+
     1     abs((q(2,i,j,k)+q(2,i,jp,k))*dx2)+
     1     abs((q(3,i,j,k)+q(3,i,j,kp))*dx3))*.5
    7 cflma=amax1(cflma,qcf)
      cflm=cflma
      return
      end
c
c  ****************************** subrout coordi  **********************
c    this subroutine assign the physical coordinates 
c    these are uniform 
c
      subroutine coordi(y)
c
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/d123/alx1,alx2,alx3
      common/corpri/yp1(m1),yp2(m2),yp3(m3)
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      pi=2.*asin(1.)
c
      do 61 i=1,n1
      x1=float(i-1)/float(n1m)
      yp1(i)=-alx1*0.5+x1*alx1
   61 continue
c
      do 62 j=1,n2
      x2=float(j-1)/float(n2m)
      yp2(j)=-alx2*0.5+x2*alx2
   62 continue
c
      do 63 k=1,n3
      x3=float(k-1)/float(n3m)
      yp3(k)=-alx3*0.5+x3*alx3
   63 continue
      return
      end
c
c  ****************************** subrout divg  **********************
c
c  this subroutine perform the calculation of divg(q)
c
      subroutine divg(qcap,vq,al)
c
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      dimension qcap(m1,m2,m3),vq(ndv,m1,m2,m3)
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/tstep/dt,beta,ren
c
c  ***** compute the divg(u) 
      do 11 kc=1,n3m
      kp=kpv(kc)
      do 11 jc=1,n2m
      jp=jpv(jc)
      do 11 ic=1,n1m
      ip=ipv(ic)
      dqcap=(vq(1,ip,jc,kc)-vq(1,ic,jc,kc))*dx1
     1     +(vq(2,ic,jp,kc)-vq(2,ic,jc,kc))*dx2
     1     +(vq(3,ic,jc,kp)-vq(3,ic,jc,kc))*dx3
      qcap(ic,jc,kc)=dqcap/(dt*al)
   11 continue
      return
      end
c  ****************************** subrout divgck  **********************
c
c  this subroutine perform the calculation of divg(q)
c
      subroutine divgck(vq,qmax)
c
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      dimension vq(ndv,m1,m2,m3)
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
c
c  ***** compute the divg(q)
      qmax=0.
      do 11 kc=1,n3m
      kp=kpv(kc)
      do 11 jc=1,n2m
      jp=jpv(jc)
      do 11 ic=1,n1m
      ip=ipv(ic)
      dqcap=(vq(1,ip,jc,kc)-vq(1,ic,jc,kc))*dx1
     1     +(vq(2,ic,jp,kc)-vq(2,ic,jc,kc))*dx2
     1     +(vq(3,ic,jc,kp)-vq(3,ic,jc,kc))*dx3
      qmax=amax1(abs(dqcap),qmax)
   11 continue
      return
      end
c
c  ****************************** subrout indic **********************
c
c  in this subroutine the indices ip,im,jp,jm,kp,km are calculated
c  these are necessary when the equation are solved for the periodic conditions.
c
      subroutine indic
c
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
c
c
c   periodic directions
c
      do 1 ic=1,n1m
      ipv(ic)=ic+1
      if(ic.eq.n1m) ipv(ic)=1
      imv(ic)=ic-1
      if(ic.eq.1) imv(ic)=n1m
    1 continue
      do 2 kc=1,n3m
      kmv(kc)=kc-1
      kpv(kc)=kc+1
      if(kc.eq.1) kmv(kc)=n3m
      if(kc.eq.n3m) kpv(kc)=1
    2 continue
      do 3 jc=1,n2m
      jpv(jc)=jc+1
      if(jc.eq.n2m) jpv(jc)=1
      jmv(jc)=jc-1
      if(jc.eq.1) jmv(jc)=n2m
    3 continue
      return
      end
c
c
c  ****************************** subrout meshes **********************
c
c  generates the mesh inverse of the spatial steps dx1, dx2, dx3,
c  and the squares of the dx1 and dx2.
c
      subroutine meshes
c
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      common/d123/alx1,alx2,alx3
      common/meshu/d1x1,d1x2,d1x3
c
      d1x1=alx1/float(n1m)
      d1x2=alx2/float(n2m)
      d1x3=alx3/float(n3m)
      dx1=1./d1x1
      dx2=1./d1x2
      dx3=1./d1x3
      dx1q=dx1*dx1
      dx2q=dx2*dx2
      dx3q=dx3*dx3
      return
      end
c
c  ****************************** subrout oldqua **********************
c
c   transfer of old time mean quantities to ??old is performed only
c   after intia or inirea.
c   In the present version is not used it works for statistical
c   steady state flows.
c
      subroutine oldqua
c
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/vmeao/vmo(ndv),vrmso(4)
      common/qmean/qm(ndv),qrms(4)
      common/skfl/ske(4),fla(4)
      common/prrm/prm,prms
      common/skflo/skeo(4),flao(4)
      common/prrmo/prmo,prmso
c
      enavo=enav
      do 410 l=1,3
      vmo(l)=qm(l)
  410 continue
      prmo=prm
      prmso=prms
      do 413 l=1,4
      skeo(l)=ske(l)
      flao(l)=fla(l)
      vrmso(l)=qrms(l)
  413 continue
      return
      end
c
c  ****************************** subrout taver **********************
c
c   calculation of time averaged mean quantities
c   veloc. and stresses
c   In the present version is not used it works for statistical
c   steady state flows.
c
      subroutine taver
c
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/vmeao/vmo(ndv),vrmso(4)
      common/qmean/qm(ndv),qrms(4)
      common/eneav/enav,enavo
      common/skfl/ske(4),fla(4)
      common/prrm/prm,prms
      common/skflo/skeo(4),flao(4)
      common/prrmo/prmo,prmso
c
      enavo=enav+enavo
      do 410 l=1,3
      vmo(l)=(vmo(l)+qm(l))
  410 continue
      prmo=(prmo+prm)
      prmso=(prmso+prms)
  411 continue
      do 413 l=1,4
      vrmso(l)=(vrmso(l)+qrms(l))
      skeo(l)=(skeo(l)+ske(l))
      flao(l)=(flao(l)+fla(l))
  413 continue
      return
      end
c
c  ******************** subrout tavwri *********************
c
c   write quantities to evaluate averages in time
c   In the present version is not used it works for statistical
c   steady state flows.
c
      subroutine tavwri(nav)
c
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/vmeao/vmo(ndv),vrmso(4)
      common/eneav/enav,enavo
      common/skflo/skeo(4),flao(4)
      common/prrmo/prmo,prmso
      character*10 tawrfi
      tawrfi='tavwri.out'
      open(70,file=tawrfi,form='unformatted')
      rewind(70)
      write(70)enavo,nav
      write(70)(vmo(l),l=1,3),prmo,prmso,
     1     (vrmso(l),l=1,4),(skeo(l),l=1,4),(flao(l),l=1,4)
   72 continue
      close(70)
      return
      end
c
c  ******************** subrout tavrea *********************
c
c   read quantities to evaluate averages in time
c   In the present version is not used it works for statistical
c   steady state flows.
c
      subroutine tavrea(nav)
c
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/vmeao/vmo(ndv),vrmso(4)
      common/eneav/enav,enavo
      common/skflo/skeo(4),flao(4)
      common/prrmo/prmo,prmso
      character*10 tawrfi
      tawrfi='tavwri.out'
      open(70,file=tawrfi,form='unformatted')
      read(70)enavo,nav
      read(70)(vmo(l),l=1,3),prmo,prmso,
     1     (vrmso(l),l=1,4),(skeo(l),l=1,4),(flao(l),l=1,4)
      close(70)
      return
      end
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
      dimension dph(m1,m2,m3),q(ndv,m1,m2,m3)
      common/tstep/dt,beta,ren
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
c
c  ***********  compute the q1 velocity component
c               v1dgf=component 1 of grad(dph)
      do 1 kc=1,n3m
      do 1 jc=1,n2m
      do 1 ic=1,n1m
      im=imv(ic)
      dfx11=(dph(ic,jc,kc)-dph(im,jc,kc))*dx1
      q(1,ic,jc,kc)=q(1,ic,jc,kc)-dfx11*dt*al
    1 continue
c
c  ***********  compute the q2 velocity component
c               v2dgf=component 2 of grad(dph)
      do 2 kc=1,n3m
      do 2 jc=1,n2m
      do 2 ic=1,n1m
      jm=jmv(jc)
      dfx22=(dph(ic,jc,kc)-dph(ic,jm,kc))*dx2
      q(2,ic,jc,kc)=q(2,ic,jc,kc)-dfx22*dt*al
    2 continue
c
c  ***********  compute the q3 velocity component
c               q3 is the cartesian component
c               v3dgf=component 3 of grad(dph)
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
c  ****************************** subrout vmaxv **********************
c
      subroutine vmaxv(q)
c
      include 'param.f'
      dimension q(ndv,m1,m2,m3)
      common/velmax/vmax(ndv),vmaxo(ndv)
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/metria/caj(m2),cac(m2)
c
c  ***** calculation of the maximum velocities
c  in order to check convergency or for stability calculations
c  to derive stability conditions.
c
      do 311 l=1,ndv
      vmaxo(l)=vmax(l)
      vmax(l)=0.
      do 310 k=1,n3m
      do 310 j=1,n2m
      do 310 i=1,n1m
      if(l.eq.1) vca=q(1,i,j,k)
      if(l.eq.2) vca=q(2,i,j,k)
      if(l.eq.3) vca=q(3,i,j,k)
      vfm=abs(vca)
      vm=vmax(l)
      vmax(l)=amax1(vm,vfm)
  310 continue
  311 continue
      return
      end
c
