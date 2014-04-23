c***********************************************************************
c                                                                       *
c  ****************************** subrout indic **********************  *
c                                                                       *
c     in this subroutine the indices ip,im,jp,jm,kp,km are calculated.  *
c                                                               c                                                                       *
c************************************************************************
      subroutine indic
      include 'param.f'
c
c
c   azimuthal periodic direction

c
      do 1 ic=1,n1m
      imv(ic)=ic-1
      if(ic.eq.1) imv(ic)=n1m
      ipv(ic)=ic+1
      if(ic.eq.n1m) ipv(ic)=1
    1 continue
c
c   vertical periodic direction

c
      do 4 kc=1,n3m
      kmv(kc)=kc-1
      kpv(kc)=kc+1
        if(kc.eq.1) kmv(kc)=kc
        if(kc.eq.n3m) kpv(kc)=kc
    4 continue
c
c     direction normal to the radial boundary

      do 3 jc=1,n2m
      jmv(jc)=jc-1
      jpv(jc)=jc+1
      if(jc.eq.1) jmv(jc)=jc
      if(jc.eq.n2m) jpv(jc)=jc
    3 continue
c
c   indices for the axis of symmetry  in the hdnl routines
c
      do i=1,n1m
      isym(i) = i + n1m/2
      if(isym(i).gt.n1m) isym(i) = isym(i) - n1m
      enddo
c
c   indices for radial direction
c
      do jc=1,n2m
        jpc(jc)=jpv(jc)-jc
        jmc(jc)=jc-jmv(jc)
      enddo
c
c   indices for axial direction
c
      do kc=1,n3m
        kpc(kc)=kpv(kc)-kc
        kmc(kc)=kc-kmv(kc)
      enddo
      return
      end
c************************************************************************
c                                                                       *
c  ****************************** subrout divg  **********************  *
c                                                                       *
c     this subroutine calculates divg(q).                               *
c                                                                       *
c************************************************************************
c
c     Here it is performed the calculation of the divergence of
c     the intermediate velocity field. This term is the source in the
c     Poisson equation for the pressure correction and for the
c     projection of the non-solenoidal velocity field into a solenoidal one
c
      subroutine divg(al)
      include 'param.f'
      qmax=0.
      qtot=0.
      do 11 kc=1,n3m
      usrnu3=dx3/g3rm(kc)
      kp=kc+1
      do 11 jc=1,n2m
      jp=jc+1
      usrnu1=dx1/rm(jc)
      usrnu2=dx2/g2rm(jc)/rm(jc)
      udt=1./(dt*al)
      do 11 ic=1,n1m
      ip=ipv(ic)
      dqcap=+(q1(ip,jc,kc)-q1(ic,jc,kc))*usrnu1
     1      +(q2(ic,jp,kc)-q2(ic,jc,kc))*usrnu2
     1      +(q3(ic,jc,kp)-q3(ic,jc,kc))*usrnu3
      qcap(ic,jc,kc)=dqcap*udt 
      qmax=max(abs(qcap(ic,jc,kc)),qmax)
      qtot=qtot+dqcap*rm(jc)*g2rm(jc)*g3rm(kc)
   11 continue

      return
      end
c                                                                       *
c  ****************************** subrout divgck  ***********************
c                                                                       *
c     this subroutine checks divg(q) to see the local residue.          *
c                                                                       *
c************************************************************************
      subroutine divgck(qmax,qtot)
      include 'param.f'
c
c     This is just a check on the divergence of the free-divergent
c     velocity field. The calculation is stopped if QMAX > RESID
c
      qtot=0.
      qmax=0.
      do 11 kc=1,n3m
      usrnu3=dx3/g3rm(kc)
      kp=kc+1
      do 11 jc=1,n2m
      jp=jc+1
      usrnu1=dx1/rm(jc)
      usrnu2=dx2/g2rm(jc)/rm(jc)
      do 11 ic=1,n1m
      ip=ipv(ic)
      dqca1= (q1(ip,jc,kc)-q1(ic,jc,kc))*usrnu1
      dqca2= (q2(ic,jp,kc)-q2(ic,jc,kc))*usrnu2
      dqca3= (q3(ic,jc,kp)-q3(ic,jc,kc))*usrnu3
      dqcap= (dqca1+dqca2+dqca3)
      qtot=qtot+dqcap*rm(jc)*g2rm(jc)*g3rm(kc)
      if(abs(dqcap).ge.qmax) then
      qmax=abs(dqcap)
      imxq=ic
      jmxq=jc
      kmxq=kc
                             endif
   11 continue
      return
      end
c************************************************************************
c                                                                       *
c  ****************************** subrout cfl  **********************   *
c                                                                       *
c************************************************************************
      subroutine cfl(cflm)
      include 'param.f'
c
c     in this routine the COURANT number is calculated.
c     This parameter determines the stability condition
c     (CFL < 1 for Adams-Bashfort  CFL < 1.7 for third order Runge-Kutta)
c     for the calculation at dt = const.
c     If a calculation at variable dt is used, this parameter determines
c     the dt itself.
c
      cflm=0.
c
      do j=1,n2m
      jp=jpv(j)
      usrnu1=dx1/rm(j)
      usrnu=dx2/g2rm(j)/rm(j)
      do k=1,n3m
      usrnu3=dx3/g3rm(k)
      kp=k+1
      do i=1,n1m
      ip=ipv(i)
      qcf1=abs((q1(i,j,k)+q1(ip,j,k))*0.5*usrnu1)
      qcf2=abs((q2(i,j,k)+q2(i,jp,k))*0.5*usrnu)
      qcf3=abs((q3(i,j,k)+q3(i,j,kp))*0.5*usrnu3)
      qcf=qcf1+qcf2+qcf3
      cflm=max(cflm,qcf)
      enddo
      enddo
      enddo
      return
      end
