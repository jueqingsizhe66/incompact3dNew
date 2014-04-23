c
      subroutine openfi
      common/newdat/itimsq,timeav
      common/names/filcnw,filcnr,filth,filou,filuu,filsf
      character*27 filcnw,filcnr
      character*17 filth,filou,filuu,filsf
      character*1 pit
      open(46,file='nfchcapn')
      read(46,'(a)')filcnw
      read(46,'(a)')filcnr
      read(46,'(a)')filth
      read(46,'(a)')filou
      read(46,'(a)')filuu
      read(46,'(a)')filsf
      open(13,file=filcnw,form='unformatted')
      open(23,file=filcnr,form='unformatted')
      open(32,file=filth)
      open(17,file=filou)
      open(18,file=filuu)
      open(19,file=filsf)
      open(20,file='chaoo.out')
      open(95,file='chabou.out')
      rewind 14
      rewind 13
      rewind 23
      rewind 32
      rewind 17
      rewind 18
      rewind 19
      if(itimsq.eq.1) then
      do l=1,5
      write(pit,81)l
   81 format(i1.1)
      nfil=70+l
      open(nfil,file='timseq.dat'//pit,form='unformatted')
      enddo
                      endif
      return
      end
c
c  ****************************** subrout timseq **********************
c
c   write the time sequence of the velocity at
c   some desired points to anlise the evolution in time
c
      subroutine timseq(time,q,pr)
c
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      dimension pr(m1,m2,m3)
      dimension q(ndv,m1,m2,m3)
      common/rhs3p/dp3ns
      common/jprts/ipri(5),kpri(5)
      jfi=n2m/2
      do l=1,5
      jin=1
      nfil=70+l
      iin=ipri(l)-3
      ifi=ipri(l)+3
      kin=kpri(l)-3
      kfi=kpri(l)+3
      write(nfil)time,dp3ns
      do i =iin,ifi
      do j =jin,jfi
      do k =kin,kfi
      write(nfil) q(1,i,j,k),
     1            q(2,i,j,k),q(3,i,j,k),pr(i,j,k)
      enddo
      enddo
      enddo
      enddo
      return
      end
c
c  ****************************** subrout contwr **********************
c   write the restarting file to be considered as
c   initial conditions when the calculation start from a previous
c   calculation
c
      subroutine contwr(ntime,time,q,ru,pr,enen,ncount)
c
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      dimension pr(m1,m2,m3)
      dimension q(ndv,m1,m2,m3)
      common/velmax/vmax(ndv),vmaxo(ndv)
      dimension ru(ndv,m1,m2,m3)
      common/tstep/dt,beta,ren
      common/wallst/cflw,cfuw
      common/veltot/vit(ndv)
      common/omegas/omr,omi
      common/timref/tref
      common/inener/ene0
      common/rhs3p/dp3ns
      common/eneav/enav,enavo,cfn,cfo
      common/qwallo/q1s(m1,m3),q2s(m1,m3),q3s(m1,m3)
      common/qwalup/q1n(m1,m3),q2n(m1,m3),q3n(m1,m3)
      common/iprf/iprfi
      common/averou/iav
      real*4 lijmij,mijmij
      common/outd1/cdyn(m2),alij(m2),amij(m2)
      common/outd2/cs(m2),lijmij(m2),mijmij(m2)
      common/outd3/accdyn(m2),acclij(m2),accmij(m2)
      common/ledat/csma,cvisc,iles
      common/names/filcnw,filcnr,filth,filou,filuu,filsf
      character*27 filcnw,filcnr
      character*17 filth,filou,filuu,filsf
      character*80 namfil
      character*4 ipfi

      itime=nint(time)
      write(ipfi,82)itime
   82 format(i4.4)
      if(iprfi.eq.1.and.iav.eq.1) then
      namfil='field'//ipfi//'.dat'
      open(13,file=namfil,form='unformatted')
      if(iles.eq.1) then
      namfil='dynsta'//ipfi//'.dat'
      open(29,file=namfil,form='unformatted')
                     endif
                     else
c     write(6,*)' wrote ',namfil
      open(13,file=filcnw,form='unformatted')
c     write(6,*)' wrote ',filcnw
      if(iles.eq.1) then
      namfil='dyncnw.dat'
      open(29,file=namfil,form='unformatted')
                     endif
                     endif
c
      nfil=13
      rewind(nfil)
      write(nfil)n1m,n2,n3m
      write(nfil)ntime,time,ene0,dp3ns,enav,cfn,dt
      write(nfil) (((q(1,i,j,k),i=1,n1m),j=1,n2),k=1,n3m),
     1            (((q(2,i,j,k),i=1,n1m),j=1,n2),k=1,n3m),
     1            (((q(3,i,j,k),i=1,n1m),j=1,n2),k=1,n3m),
     1            (((pr(i,j,k),i=1,n1m),j=1,n2),k=1,n3m)
      write(nfil) ((q1s(i,k),i=1,n1m),k=1,n3m),
     1            ((q2s(i,k),i=1,n1m),k=1,n3m),
     1            ((q3s(i,k),i=1,n1m),k=1,n3m),
     1            ((q1n(i,k),i=1,n1m),k=1,n3m),
     1            ((q2n(i,k),i=1,n1m),k=1,n3m),
     1            ((q3n(i,k),i=1,n1m),k=1,n3m)
      write(nfil)ntime,time,vit(1),vit(2),vit(3)
     1         ,dp3ns,cfuw,enen,vmax(2),vmax(3),cflw
      close(nfil)
      if(iles.eq.1) then
      nfil=29
      rewind(nfil)
      do j=1,n2m
      write(29)ncount,accmij(j),acclij(j),accdyn(j)
      enddo
      close(nfil)
                    endif
      return
      end
c
c  ****************************** subrout enerca **********************
c
c   calculation of total energy = (v1**2+v2**2+v3**2)*0.5
c   and other quantities derived by fluctuating
c   velocity and vorticity
c
      subroutine enerca(q,enej,pr)
c
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      dimension q(ndv,m1,m2,m3),pr(m1,m2,m3)
      common/sma/vis(m1,m2,m3)
      common/strain/st(m1,m2,m3,6)
      dimension vmc2(m2)
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      common/meshu/d1x1,d1x2,d1x3
      common/metria/caj(m2),cac(m2)
      common/vmean/vm(ndv,m2),vrms(6,m2)
      common/qmean/qm(ndv,m2),qrms(6,m2)
      common/tstep/dt,beta,re
      common/d13/alx1,alx3
      common/veltot/vit(ndv)
      common/vini0/v30(m2)
      common/eneav/enav,enavo,cfn,cfo
      common/eneto/enet,enstro
      common/wallst/cflw,cfuw
      common/skfl/ske(4,m2),fla(4,m2)
      common/prrm/prm(m2),prms(m2)
      common/skflo/skeo(4,m2),flao(4,m2)
      common/prrmo/prmo(m2),prmso(m2)
      common/qwallo/q1s(m1,m3),q2s(m1,m3),q3s(m1,m3)
      common/qwalup/q1n(m1,m3),q2n(m1,m3),q3n(m1,m3)
      common/islwal/islv1s,islv1n,islv3s,islv3n
      common/y2sta/y2s(m2)
      common/corpri/yp1(m1),yp2(m2),yp3(m3)
      common/visme/vismed(m2),visrms(m2)
      common/vismo/vismeo(m2),visrmo(m2)
      common/sgsp/sgs(6,m2),sgso(6,m2)
      common/ledat/csma,cvisc,iles

c
      vol=2.
      vl13=1./float(n1m*n3m)
      vl123=1./float(n1m*n2m*n3m)
      vit(1)=0.
      vit(2)=0.
      vit(3)=0.
      do 410 j=1,n2m
      jp=j+1
      vm1m=0.
      vm2m=0.
      vm2c=0.
      vm3m=0.
      pmm=0.
      vismm=0.
      do 411 k=1,n3m
      kp=kpv(k)
      do 411 i=1,n1m
      ip=ipv(i)
      vm1m=(q(1,ip,j,k)+q(1,i,j,k))*0.5+vm1m
      vm2m=(q(2,i,jp,k)+q(2,i,j,k))*0.5+vm2m
      vm2c=q(2,i,j,k)+vm2c
      vm3m=(q(3,i,j,kp)+q(3,i,j,k))*0.5+vm3m
      pmm=pmm+pr(i,j,k) 
      vismm=vismm+vis(i,j,k) 
      vit(1)=(q(1,ip,j,k)+q(1,i,j,k))*0.5*caj(j)+vit(1)
      vit(2)=(q(2,i,jp,k)+q(2,i,j,k))*0.5*caj(j)+vit(2)
      vit(3)=(q(3,i,j,kp)+q(3,i,j,k))*0.5*caj(j)+vit(3)
  411 continue
      vm(1,j)=vm1m*vl13
      vm(2,j)=vm2m*vl13
      vmc2(j)=vm2c*vl13
      vm(3,j)=vm3m*vl13
      prm(j)=pmm*vl13
      vismed(j)=vismm*vl13
  410 continue
      vit(1)=vit(1)*vl123/vol
      vit(2)=vit(2)*vl123/vol
      vit(3)=vit(3)*vl123/vol
      enej=0.
      enet=0.
      do j=1,n2m
      jp=j+1
      u1rmm=0.
      u2rmm=0.
      u3rmm=0.
      u13mm=0.
      u23mm=0.
      u12mm=0.
      ppmm=0.1e-09
      vism=0.
      do k=1,n3m
      kp=kpv(k)
      do i=1,n1m
      ip=ipv(i)
      v1t=(q(1,ip,j,k)+q(1,i,j,k))*0.5
      v2t=(q(2,i,jp,k)+q(2,i,j,k))*0.5
      v3t=(q(3,i,j,kp)+q(3,i,j,k))*0.5
      v1m=v1t-vm(1,j)
      v2m=v2t-vm(2,j)
      v3m=v3t-vm(3,j)  
      visl=vis(i,j,k)-vismed(j)
      ppm=pr(i,j,k)-prm(j)
      u1rmm=u1rmm+v1m**2
      u2rmm=u2rmm+v2m**2
      u3rmm=u3rmm+v3m**2
      u23mm=u23mm+v2m*v3m
      u13mm=u13mm+v1m*v3m
      u12mm=u12mm+v1m*v2m
      ppmm=ppmm+ppm**2
      vism=vism+visl**2
      enejp=(v1m**2+v2m**2+v3m**2)*0.25*vl123*caj(j)
      enetp=(v1t**2+v2t**2+v3t**2)*0.25*vl123*caj(j)
      enej=enej+enejp
      enet=enet+enetp
      enddo
      enddo
      visrms(j)=vism*vl13
      vrms(1,j)=u1rmm*vl13
      vrms(2,j)=u2rmm*vl13
      vrms(3,j)=u3rmm*vl13
      vrms(4,j)=u23mm*vl13
      vrms(5,j)=u13mm*vl13
      vrms(6,j)=u12mm*vl13
      prms(j)=ppmm*vl13
      enddo
      do j=1,n2m
      jp=j+1
      sk1m=0.
      sk2m=0.
      sk3m=0.
      sk4m=0.
      fl1m=0.
      fl2m=0.
      fl3m=0.
      fl4m=0.
      do k=1,n3m
      kp=kpv(k)
      do i=1,n1m
      ip=ipv(i)
      v1t=(q(1,ip,j,k)+q(1,i,j,k))*0.5
      v2t=(q(2,i,jp,k)+q(2,i,j,k))*0.5
      v3t=(q(3,i,j,kp)+q(3,i,j,k))*0.5
      v1m=v1t-vm(1,j)
      v2m=v2t-vm(2,j)
      v3m=v3t-vm(3,j)  
      sk1m=sk1m+(v1m/sqrt(vrms(1,j)))**3
      sk2m=sk2m+(v2m/sqrt(vrms(2,j)))**3
      sk3m=sk3m+(v3m/sqrt(vrms(3,j)))**3
      sk4m=sk4m+(ppm/sqrt(prms(j)))**3
      fl1m=fl1m+(v1m/sqrt(vrms(1,j)))**4
      fl2m=fl2m+(v2m/sqrt(vrms(2,j)))**4
      fl3m=fl3m+(v3m/sqrt(vrms(3,j)))**4
      fl4m=fl4m+(ppm/sqrt(prms(j)))**4
      enddo
      enddo
      ske(1,j)=sk1m*vl13
      ske(2,j)=sk2m*vl13
      ske(3,j)=sk3m*vl13
      ske(4,j)=sk4m*vl13
      fla(1,j)=fl1m*vl13   
      fla(2,j)=fl2m*vl13   
      fla(3,j)=fl3m*vl13   
      fla(4,j)=fl4m*vl13  
      enddo
      cfm=(cflw-cfuw)*0.5
      cfn=cfm
      enav=enej
      do 427 j=1,n2m
      do 428 l=1,3
      qrms(l,j)=vrms(l,j)
      qm(l,j)=vm(l,j)
  428 continue
      qrms(4,j)=vrms(4,j)
      qrms(5,j)=vrms(5,j)
      qrms(6,j)=vrms(6,j)
  427 continue
      enstro=0.
      do kc=1,n3m
      km=kmv(kc)
      do jc=1,n2m
      do ic=1,n1m
      im=imv(ic)
c
c      OM_y COMPONENT
c
      omy=+(q(1,ic,jc,kc)-q(1,ic,jc,km))*dx3
     1    -(q(3,ic,jc,kc)-q(3,im,jc,kc))*dx1
      enstro=enstro+.5*(omy*omy)*vl123*caj(jc)
      enddo
      enddo
      enddo
      do kc=1,n3m
      km=kmv(kc)
      do jc=2,n2m-1
      jm=jmv(jc)
      jp=jpv(jc)
      do ic=1,n1m
      im=imv(ic)
c
c      OM_x COMPONENT
c
      omxp=+(q(3,ic,jp,kc)-q(3,ic,jc,kc))*dx2/cac(jp)
     1     -(q(2,ic,jp,kc)-q(2,ic,jp,km))*dx3
      omxm=+(q(3,ic,jc,kc)-q(3,ic,jm,kc))*dx2/cac(jc)
     1     -(q(2,ic,jc,kc)-q(2,ic,jc,km))*dx3
      omx=(omxp+omxm)*0.5
c
c      OM_z COMPONENT
c
      omzp=-(q(1,ic,jp,kc)-q(1,ic,jc,kc))*dx2/cac(jp)
     1     +(q(2,ic,jp,kc)-q(2,im,jp,kc))*dx1
      omzm=-(q(1,ic,jc,kc)-q(1,ic,jm,kc))*dx2/cac(jc)
     1     +(q(2,ic,jc,kc)-q(2,im,jc,kc))*dx1
      omz=(omzp+omzm)*0.5
      enstro=enstro+.5*(omx*omx+omz*omz)*vl123*caj(jc)
      enddo
      enddo
      enddo
c
c   lower wall
c
      jc=1
      jp=jpv(jc)
      do kc=1,n3m
      km=kmv(kc)
      do ic=1,n1m
      im=imv(ic)
c
c      OM_x COMPONENT
c
      omxp=+(q(3,ic,jp,kc)-q(3,ic,jc,kc))*dx2/cac(jp)
     1     -(q(2,ic,jp,kc)-q(2,ic,jp,km))*dx3
      omxm=+(q(3,ic,jc,kc)-q3s(ic,kc))/(y2s(jc)-yp2(jc))*(1-islv3s)
      omx=(omxp+omxm)*0.5
c
c      OM_z COMPONENT
c
      omzp=-(q(1,ic,jp,kc)-q(1,ic,jc,kc))*dx2/cac(jp)
     1     +(q(2,ic,jp,kc)-q(2,im,jp,kc))*dx1
      omzm=-(q(1,ic,jc,kc)-q1s(ic,kc))/(y2s(jc)-yp2(jc))*(1-islv1s)
      omz=(omzp+omzm)*0.5
      enstro=enstro+.5*(omx*omx+omz*omz)*vl123*caj(jc)
      enddo
      enddo
c
c   upper wall
c
      jc=n2m
      jm=jmv(jc)
      do kc=1,n3m
      km=kmv(kc)
      do ic=1,n1m
      im=imv(ic)
c
c      OM_x COMPONENT
c
      omxp=+(q3s(ic,kc)-q(3,ic,jc,kc))/(yp2(n2)-y2s(jc))*(1-islv3n)
      omxm=+(q(3,ic,jc,kc)-q(3,ic,jm,kc))*dx2/cac(jc)
     1     -(q(2,ic,jc,kc)-q(2,ic,jc,km))*dx3
      omx=(omxp+omxm)*0.5
c
c      OM_z COMPONENT
c
      omzp=-(q1s(ic,kc)-q(1,ic,jc,kc))/(yp2(n2)-y2s(jc))*(1-islv1n)
      omzm=-(q(1,ic,jc,kc)-q(1,ic,jm,kc))*dx2/cac(jc)
     1     +(q(2,ic,jc,kc)-q(2,im,jc,kc))*dx1
      omz=(omzp+omzm)*0.5
      enstro=enstro+.5*(omx*omx+omz*omz)*vl123*caj(jc)
      enddo
      enddo
      do j=1,n2m
      do l=1,6
      sgs(l,j)=0.
      enddo
      enddo
      if(iles.ne.0) then
      call strai(q)
      do j=1,n2m
      do l=1,6
      sgs(l,j)=0.
      enddo
      do k=1,n3m
      do i=1,n1m
      sgs11=2.*vis(i,j,k)*st(i,j,k,1)
      sgs22=2.*vis(i,j,k)*st(i,j,k,2)
      sgs33=2.*vis(i,j,k)*st(i,j,k,3)
      sgs12=2.*vis(i,j,k)*st(i,j,k,4)
      sgs13=2.*vis(i,j,k)*st(i,j,k,5)
      sgs23=2.*vis(i,j,k)*st(i,j,k,6)
      sgs(1,j)=sgs(1,j)+sgs11 
      sgs(2,j)=sgs(2,j)+sgs22 
      sgs(3,j)=sgs(3,j)+sgs33 
      sgs(4,j)=sgs(4,j)+sgs23 
      sgs(5,j)=sgs(5,j)+sgs13 
      sgs(6,j)=sgs(6,j)+sgs12 
      enddo
      enddo
      do l=1,6
      sgs(l,j)=sgs(l,j)*vl13
      enddo
      enddo
                     endif
      return
      end
c
c  ****************************** subrout inirea **********************
c
c   reading initial  conditions from file
c   when calculation are continued from a previous
c   calculation.
c   This reading should be consistent with contwr
c
      subroutine inirea(ntii,time,q,ru,pr)
c
      include 'param.f'
      dimension pr(m1,m2,m3)
      common/velmax/vmax(ndv),vmaxo(ndv)
      common/dim/n1,n1m,n2,n2m,n3,n3m
      dimension q(ndv,m1,m2,m3)
      dimension ru(ndv,m1,m2,m3)
      common/tstep/dt,beta,ren
      common/wallst/cflw,cfuw
      common/veltoo/vit(ndv)
      common/omegas/omr,omi
      common/timref/tref
      common/inener/ene0
      common/rhs3p/dp3ns
      common/eneav/enav,enavo,cfn,cfo
      common/metria/caj(m2),cac(m2)
      common/d13/alx1,alx3
      common/oldso/ifield
      common/qwallo/q1s(m1,m3),q2s(m1,m3),q3s(m1,m3)
      common/qwalup/q1n(m1,m3),q2n(m1,m3),q3n(m1,m3)
      real*4 lijmij,mijmij
      common/outd1/cdyn(m2),alij(m2),amij(m2)
      common/outd2/cs(m2),lijmij(m2),mijmij(m2)
      common/outd3/accdyn(m2),acclij(m2),accmij(m2)
      common/ledat/csma,cvisc,iles
      dimension y2sol(2*m2)
      common/codat/inicos
      common/cordo/ichy2
      common/dimo/n2o,str2o
      character*80 namfil
c
      nfil=23
      rewind(nfil)
      read(nfil)n1lm,n2o,n3lm
      read(nfil)ntime,time,ene0,dp3ns,enavo,cfo,dtold
      read(nfil) (((q(1,i,j,k),i=1,n1lm),j=1,n2o),k=1,n3lm),
     1            (((q(2,i,j,k),i=1,n1lm),j=1,n2o),k=1,n3lm),
     1            (((q(3,i,j,k),i=1,n1lm),j=1,n2o),k=1,n3lm),
     1            (((pr(i,j,k),i=1,n1lm),j=1,n2o),k=1,n3lm)
      read(nfil)  ((q1s(i,k),i=1,n1lm),k=1,n3lm),
     1            ((q2s(i,k),i=1,n1lm),k=1,n3lm),
     1            ((q3s(i,k),i=1,n1lm),k=1,n3lm),
     1            ((q1n(i,k),i=1,n1lm),k=1,n3lm),
     1            ((q2n(i,k),i=1,n1lm),k=1,n3lm),
     1            ((q3n(i,k),i=1,n1lm),k=1,n3lm)
      read(nfil)ntime,time,vit(1),vit(2),vit(3)
     1         ,dp3ns,cfuw,enen,vmax(2),vmax(3),cflm
c     write(6,158)ntime,time,vit(1),vit(2),vit(3)
c    1         ,dp3ns,cfuw,enen,vmax(2),vmax(3),cflm
  158 format(3x,i4,4(1x,e11.5),2x,8(1x,e11.5))
      close(nfil)
      if(iles.eq.1.and.inicos.ne.0) then
      namfil='dyncnw.dat'
      open(29,file=namfil,form='unformatted')
      nfil=29
      rewind(nfil)
      do j=1,n2o-1
      read(29)ncount,accmij(j),acclij(j),accdyn(j)
      enddo
      close(nfil)
                    endif
      if(n1lm.eq.2*n1m.and.n3lm.eq.2*n3m) then
      write(6,*)' inx1x3 fine grid n1lm=',n1lm,' n3lm=',n3lm,
     1          ' to coarse grid n1m=',n1m,' n3m=',n3m
      call inx1x3(n1lm,n3lm,q,ru,pr)
                                                    endif
      if(n1lm.eq.2*n1m.and.n3lm.eq.n3m) then
      write(6,*)' inx1x1 fine grid n1lm=',n1lm,' n3lm=',n3lm,
     1          ' to coarse grid n1m=',n1m,' n3m=',n3m
      call inx1x1(n1lm,q,ru,pr)
                                            endif
      if(n1lm.eq.n1m.and.n3lm.eq.2*n3m) then
      write(6,*)' inx3x3 fine grid n1lm=',n1lm,' n3lm=',n3lm,
     1          ' to coarse grid n1m=',n1m,' n3m=',n3m
      call inx3x3(n3lm,q,ru,pr)
                                            endif
      if(n1m.eq.2*n1lm.and.n3m.eq.n3lm) then
      write(6,*)' exx1x1 coarse grid n1lm=',n1lm,' n3lm=',n3lm,
     1          ' to fine grid n1m=',n1m,' n3m=',n3m
      call exx1x1(n1lm,n3lm,q,ru,pr)
                                                    endif
      if(n1m.eq.n1lm.and.n3m.eq.2*n3lm) then
      write(6,*)' exx3x3 coarse grid n1lm=',n1lm,' n3lm=',n3lm,
     1          ' to fine grid n1m=',n1m,' n3m=',n3m
      call exx3x3(n1lm,n3lm,q,ru,pr)
                                                    endif
      if(n1m.eq.2*n1lm.and.n3m.eq.2*n3lm) then
      write(6,*)' exx1x3 coarse grid n1lm=',n1lm,' n3lm=',n3lm,
     1          ' to fine grid n1m=',n1m,' n3m=',n3m
      call exx1x3(n1lm,n3lm,q,ru,pr)
                                                    endif
      if(n2o.ne.n2) ichy2=1
      if(ichy2.eq.1) then 
      call cordio(y2sol)
      q1mao=0.
      q2mao=0.
      q3mao=0.
      do j=1,n2m
      do k=1,n3m
      do i=1,n1m
      q1mao=max(abs(q(1,i,j,k)),q1mao)
      q2mao=max(abs(q(2,i,j,k)),q2mao)
      q3mao=max(abs(q(3,i,j,k)),q3mao)
      enddo
      enddo
      enddo
      call intrpr(q,ru,y2sol)
      q1man=0.
      q2man=0.
      q3man=0.
      do j=1,n2m
      do k=1,n3m
      do i=1,n1m
      q1man=max(abs(q(1,i,j,k)),q1man)
      q2man=max(abs(q(2,i,j,k)),q2man)
      q3man=max(abs(q(3,i,j,k)),q3man)
      enddo
      enddo
      enddo
      write(6,*)'   q max old ',q1mao,q2mao,q3mao
      write(6,*)'   q max new ',q1man,q2man,q3man
                     endif
      call divgck(q,qmax)
      write(6,*)' in inirea after read divg max=',qmax
      return
      end
c
c  ****************************** subrout coordi  **********************
c    this subroutine assign the physical coordinates as function
c    of the computational ones x2
c    coordinate clustering near the wallstretchings in x2 direction can be   assumed
c    This is the coordinate of a previous  calculation
c   The coordinate in the new simulation is different and
c   the velocity field is interpolated to restart
c
      subroutine cordio(y2sol)
c
      include 'param.f'
      common/dimo/n2o,str2o
      dimension y2sol(2*m2),eta(2*m2),yp2ol(2*m2)
      common/y2sta/y2s(m2)
c
      n2om=n2o-1
      dx2o=float(n2om)
      tstr2=tanh(str2o*0.5)
      do 63 j=1,n2o
      x2=(j-1)/float(n2om)
      eta(j)=0.5*(1.+tanh(str2o*(x2-0.5))/tstr2)
   63 continue
      do 65 j=1,n2o
      yp2ol(j)=(-0.5+eta(j))*2.
   65 continue
c     print *,'y al centro',y(n2/2+1)
      open(67,file='coordneol.out')
      do 67 j=1,n2om
      y2sol(j)=(yp2ol(j)+yp2ol(j+1))*0.5
      x2=(j-1)/float(n2om)
      write(67,*)x2,y2s(j),y2sol(j)
   67 continue
      close(67)
      return
      end
c
c  *************  subroutine intrpr  ************
c
      subroutine intrpr(q,ru,y2sol)
      include'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/dimo/n2o,str2o
      dimension fo(2*m2),yn(2*m2),fni(2*m2),yu(2*m2)
      dimension y2sol(2*m2)
      common/y2sta/y2s(m2)
      dimension q(ndv,m1,m2,m3)
      dimension ru(ndv,m1,m2,m3)
      common/veltoo/vito(ndv)
      dimension vitn(ndv),vit(ndv)
      common/metria/caj(m2),cac(m2)

c
c in this subr evaluate the velocity, read from a continuation
c file, by an cubic spline interpolation in the new grid
c
c
c
      n2om=n2o-1
      write(6,*)'n2o=',n2o
      do 1 k=1,n3m
      do 1 i=1,n1m
      do 10 j=1,n2om
      fo(j)=q(1,i,j,k)
   10 continue
      dpn2=0.
      dp1=0.
c
c
      call spline(y2sol,fo,n2om,dp1,dpn2,yu)
c
      do 11 j=1,n2m
      xx=y2s(j)
      call splint(y2sol,fo,yu,n2om,xx,yy)
      fni(j)=yy
   11 continue
      do 12 j=1,n2m
      ru(1,i,j,k)=fni(j)
   12 continue
    1 continue
c
c
c   velocity q3 b.c. islv3s=0 old derivative islv3s=1 zero derivative
c   velocity q3 b.c. islv3n=0 old derivative islv3n=1 zero derivative
c
      do 2 k=1,n3m
      do 2 i=1,n1m
      do 20 j=1,n2om
      fo(j)=q(3,i,j,k)
   20 continue
      dpn2=0.
      dp1=0.
c
c
      call spline(y2sol,fo,n2om,dp1,dpn2,yu)
c
      do 21 j=1,n2m
      xx=y2s(j)
      call splint(y2sol,fo,yu,n2om,xx,yy)
      fni(j)=yy
   21 continue
      do 22 j=1,n2m
      ru(3,i,j,k)=fni(j)
   22 continue
    2 continue
      vol=2.
      vl123=1./float(n1m*n2m*n3m)
      vitn(1)=0.
      vitn(3)=0.
      do j=1,n2m
      do k=1,n3m
      kp=kpv(k)
      do i=1,n1m
      ip=ipv(i)
      vitn(1)=(ru(1,ip,j,k)+ru(1,i,j,k))*0.5*caj(j)+vitn(1)
      vitn(3)=(ru(3,i,j,kp)+ru(3,i,j,k))*0.5*caj(j)+vitn(3)
      enddo
      enddo
      enddo
      vitn(1)=vitn(1)*vl123/vol
      vitn(3)=vitn(3)*vl123/vol
      write(6,*)' old  vit(1,vit(3 ',vito(1),vito(3)
      write(6,*)' new  vit(1,vit(3 ',vitn(1),vitn(3)
      do k=1,n3m
      do i=1,n1m
      do j=1,n2om
      ru(3,i,j,k)=ru(3,i,j,k)*vito(3)/vitn(3)
      ru(1,i,j,k)=ru(1,i,j,k)*vito(1)/vitn(1)
      enddo
      enddo
      enddo
c
c
c
      call q2div(q,ru)
      return
      end
      subroutine spline(xin,y,n,rp1,rpn,y2)
      include'param.f'
      dimension xin(2*m2),y(2*m2),y2(2*m2),u(2*m2)
      if (rp1.gt..99e30) then
      y2(1)=0.
      u(1)=0.
      else
      y2(1)=-0.5
      u(1)=(3./(xin(2)-xin(1)))*((y(2)-y(1))
     *     /(xin(2)-xin(1))-rp1)
      endif
      do 11 i=2,n-1
      sig=(xin(i)-xin(i-1))/(xin(i+1)-xin(i-1))
      pnn=sig*y2(i-1)+2.
      y2(i)=(sig-1.)/pnn
      u(i)=(6.*((y(i+1)-y(i))
     *   /(xin(i+1)-xin(i))-(y(i)-y(i-1))
     *   /(xin(i)-xin(i-1)))
     *   /(xin(i+1)-xin(i-1))-sig*u(i-1))/pnn
11    continue
      if (rpn.gt..99e30) then
      qn=0.
      un=0.
      else
      qn=0.5
      un=(3./(xin(n)-xin(n-1)))*(rpn-(y(n)-y(n-1))
     *   /(xin(n)-xin(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do 12 k=n-1,1,-1
      y2(k)=y2(k)*y2(k+1)+u(k)
12    continue
      return
      end
      subroutine splint(xain,ya,y2a,n,x,y)
      dimension xain(1),ya(1),y2a(1)
      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xain(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xain(khi)-xain(klo)
      if (h.eq.0.) pause 'bad xa input.'
      a=(xain(khi)-x)/h
      b=(x-xain(klo))/h
      y=a*ya(klo)+b*ya(khi)+
     *      ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
      return
      end
c************************************************************************
c                                                                       *
c  ****************************** subrout inx1x3 ********************** *
c     interpolates in coarser grid both in x1 and x3 direction     *
c                                                                       *
c************************************************************************
      subroutine inx1x3(n1lm,n3lm,q,ru,pr) 
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      dimension pr(m1,m2,m3)
      dimension q(ndv,m1,m2,m3)
      dimension ru(ndv,m1,m2,m3)
      do j=1,n2m
      do i=1,n1m
      inm=2*(i-1)
      if(i.eq.1) inm=2*n1m
      inc=2*(i-1)+1
      inp=2*i
      do k=1,n3m
      knm=2*(k-1)+1
      knp=2*k
      ru(1,i,j,k)=(q(1,inp,j,knp)+q(1,inp,j,knm)
     1           +q(1,inm,j,knp)+q(1,inm,j,knm)
     1       +2.*(q(1,inc,j,knp)+q(1,inc,j,knm)) )/8.
      pr(i,j,k)=0.
      enddo
      enddo
      enddo
      do j=1,n2m
      do k=1,n3m
      knm=2*(k-1)
      if(k.eq.1) knm=2*n3m
      knc=2*(k-1)+1
      knp=2*k
      do i=1,n1m
      inm=2*(i-1)+1
      inp=2*i
      ru(3,i,j,k)=(q(3,inp,j,knp)+q(3,inp,j,knm)
     1           +q(3,inm,j,knp)+q(3,inm,j,knm)
     1       +2.*(q(3,inm,j,knc)+q(3,inp,j,knc)) )/8.
      enddo
      enddo
      enddo
      call q2div(q,ru)
      return                                                            
      end                                                               
c************************************************************************
c                                                                       *
c  ****************************** subrout inx1x1 ********************** *
c     interpolates in coarser grid  in x1  direction                 *
c                                                                       *
c************************************************************************
      subroutine inx1x1(n1lm,q,ru,pr) 
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      dimension pr(m1,m2,m3)
      dimension q(ndv,m1,m2,m3)
      dimension ru(ndv,m1,m2,m3)
      do j=1,n2m
      do i=1,n1m
      inm=2*(i-1)
      if(i.eq.1) inm=2*n1m
      inc=2*(i-1)+1
      inp=2*i
      do k=1,n3m
      ru(1,i,j,k)=(q(1,inp,j,k)+q(1,inm,j,k)+2.*q(1,inc,j,k) )/4.
      pr(i,j,k)=0.
      enddo
      enddo
      enddo
      do j=1,n2m
      do k=1,n3m
      do i=1,n1m
      inm=2*(i-1)+1
      inp=2*i
      ru(3,i,j,k)=(q(3,inp,j,k)+q(3,inp,j,k))/2.
      enddo
      enddo
      enddo
      call q2div(q,ru)
      return                                                            
      end                                                               
c************************************************************************
c                                                                       *
c  ****************************** subrout exx1x3 ********************** *
c     extrapolates fro coarse grid in x1 and x3 direction          *
c                                                                       *
c************************************************************************
      subroutine exx1x3(n1lm,n3lm,q,ru,pr) 
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      dimension pr(m1,m2,m3)
      dimension q(ndv,m1,m2,m3)
      dimension ru(ndv,m1,m2,m3)
      do j=1,n2m
      do in=1,n1lm
      ifn=2*in-1
      ifp=2*in
      ip=in+1
      if(in.eq.n1lm) ip=1
      do kn=1,n3lm
      kp=kn+1
      kfn=2*kn
      kfp=2*kn+1
      if(kn.eq.n3lm) then
      kfp=1
      kp=1
                    endif
      ru(1,ifn,j,kfn)=(3.*q(1,in,j,kn)+1.*q(1,in,j,kp))/4.
      ru(1,ifn,j,kfp)=(1.*q(1,in,j,kn)+3.*q(1,in,j,kp))/4.
      ru(1,ifp,j,kfn)=(3.*q(1,in,j,kn)+1.*q(1,in,j,kp))/8.
     1              +(3.*q(1,ip,j,kn)+1.*q(1,ip,j,kp))/8.
      ru(1,ifp,j,kfp)=(1.*q(1,in,j,kn)+3.*q(1,in,j,kp))/8.
     1              +(1.*q(1,ip,j,kn)+3.*q(1,ip,j,kp))/8.
      pr(ifn,j,kfn)=0.
      pr(ifn,j,kfp)=0.
      pr(ifp,j,kfn)=0.
      pr(ifp,j,kfp)=0.
      enddo
      enddo
      enddo
      do j=1,n2m
      do kn=1,n3lm
      kfn=2*kn-1
      kfp=2*kn
      kp=kn+1
      if(kn.eq.n3lm) kp=1
      do in=1,n1lm
      ip=in+1
      ifn=2*in
      ifp=2*in+1
      if(in.eq.n1lm) then
      ifp=1
      ip=1
                    endif
      ru(3,ifn,j,kfn)=(3.*q(3,in,j,kn)+1.*q(3,ip,j,kn))/4.
      ru(3,ifp,j,kfn)=(1.*q(3,in,j,kn)+3.*q(3,ip,j,kn))/4.
      ru(3,ifn,j,kfp)=(3.*q(3,in,j,kn)+1.*q(3,ip,j,kn))/8.
     1              +(3.*q(3,in,j,kp)+1.*q(3,ip,j,kp))/8.
      ru(3,ifp,j,kfp)=(1.*q(3,in,j,kn)+3.*q(3,ip,j,kn))/8.
     1              +(1.*q(3,in,j,kp)+3.*q(3,ip,j,kp))/8.
      enddo
      enddo
      enddo
      call q2div(q,ru)
      return                                                            
      end                                                               
c************************************************************************
c                                                                       *
c  ****************************** subrout exx1x1 ********************** *
c     extrapolates fro coarse grid in x1  direction          *
c                                                                       *
c************************************************************************
      subroutine exx1x1(n1lm,n3lm,q,ru,pr) 
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      dimension pr(m1,m2,m3)
      dimension q(ndv,m1,m2,m3)
      dimension ru(ndv,m1,m2,m3)
      do j=1,n2m
      do in=1,n1lm
      ifn=2*in-1
      ifp=2*in
      ip=in+1
      if(in.eq.n1lm) ip=1
      do k=1,n3m
      ru(1,ifn,j,k)=q(1,in,j,k)
      ru(1,ifp,j,k)=(q(1,in,j,k)+q(1,ip,j,k))/2.
      pr(ifn,j,k)=0.
      pr(ifp,j,k)=0.
      enddo
      enddo
      enddo
      do j=1,n2m
      do k=1,n3m
      do in=1,n1lm
      ip=in+1
      ifn=2*in
      ifp=2*in+1
      if(in.eq.n1lm) then
      ifp=1
      ip=1
                    endif
      ru(3,ifn,j,k)=(3.*q(3,in,j,k)+1.*q(3,ip,j,k))/4.
      ru(3,ifp,j,k)=(1.*q(3,in,j,k)+3.*q(3,ip,j,k))/4.
      enddo
      enddo
      enddo
      call q2div(q,ru)
      return                                                            
      end                                                               
c************************************************************************
c                                                                       *
c  ****************************** subrout exx3x3 ********************** *
c     extrapolates fro coarse grid in x3  direction          *
c     values of q1 and q3 components in a fine grid in  Zeta   *
c                                                                       *
c************************************************************************
      subroutine exx3x3(n1lm,n3lm,q,ru,pr) 
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      dimension pr(m1,m2,m3)
      dimension q(ndv,m1,m2,m3)
      dimension ru(ndv,m1,m2,m3)
      do j=1,n2m
      do i=1,n1m
      do kn=1,n3lm
      kp=kn+1
      kfn=2*kn
      kfp=2*kn+1
      if(kn.eq.n3lm) then
      kfp=1
      kp=1
                    endif
      ru(1,i,j,kfn)=(3.*q(1,i,j,kn)+1.*q(1,i,j,kp))/4.
      ru(1,i,j,kfp)=(1.*q(1,i,j,kn)+3.*q(1,i,j,kp))/4.
      pr(i,j,kfn)=0.
      pr(i,j,kfp)=0.
      enddo
      enddo
      enddo
      do j=1,n2m
      do kn=1,n3lm
      kfn=2*kn-1
      kfp=2*kn
      kp=kn+1
      if(kn.eq.n3lm) kp=1
      do i=1,n1m
      ru(3,i,j,kfn)=q(3,i,j,kn)
      ru(3,i,j,kfp)=(q(3,i,j,kn)+q(3,i,j,kp))/2.
      enddo
      enddo
      enddo
      call q2div(q,ru)
      return                                                            
      end                                                               
c************************************************************************
c                                                                       *
c  ****************************** subrout q2div ********************** *
c    evaluates q2 from div q=0 and put dq1,dq2,dq3 in q1,q2,q3          *
c                                                                       *
c************************************************************************
      subroutine q2div(q,ru)
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      dimension q(ndv,m1,m2,m3)
      dimension ru(ndv,m1,m2,m3)
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/metria/caj(m2),cac(m2)
c
c  ***** compute the divg(u)
      do kc=1,n3m
      kp=kpv(kc)
      do ic=1,n1m
      ip=ipv(ic)
      ru(2,ic,1,kc)=0.
      do jc=1,n2m
      jp=jc+1
      sucaj=dx2/caj(jc)
      dqcap=(ru(1,ip,jc,kc)-ru(1,ic,jc,kc))*dx1
     1     +(ru(3,ic,jc,kp)-ru(3,ic,jc,kc))*dx3                          
      ru(2,ic,jp,kc)=ru(2,ic,jc,kc)-dqcap/sucaj
      enddo
      enddo
      enddo
      do kc=1,n3m                                                    
      do ic=1,n1m                                                   
      do jc=1,n2                                                     
      q(2,ic,jc,kc)=ru(2,ic,jc,kc)
      enddo
      enddo
      enddo
      do kc=1,n3m                                                    
      do ic=1,n1m                                                   
      do jc=1,n2m
      q(1,ic,jc,kc)=ru(1,ic,jc,kc)
      q(3,ic,jc,kc)=ru(3,ic,jc,kc)
      enddo
      enddo
      enddo
      return                                                            
      end                                                               
c************************************************************************
c                                                                       *
c  ****************************** subrout inx3x3 ********************** *
c     interpolates in coarser grid in  x3 direction                   *
c                                                                       *
c***************************************************************c************************************************************************
      subroutine inx3x3(n3lm,q,ru,pr)
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      dimension pr(m1,m2,m3)
      dimension q(ndv,m1,m2,m3)
      dimension ru(ndv,m1,m2,m3)
      do j=1,n2m
      do i=1,n1m
      do k=1,n3m
      knm=2*(k-1)+1
      knp=2*k
      ru(1,i,j,k)=(q(1,i,j,knp)+q(1,i,j,knm))/2.
      pr(i,j,k)=0.
      enddo
      enddo
      enddo
      do j=1,n2m
      do k=1,n3m
      knm=2*(k-1)
      if(k.eq.1) knm=2*n3m
      knc=2*(k-1)+1
      knp=2*k
      do i=1,n1m
      ru(3,i,j,k)=(q(3,i,j,knp)+q(3,i,j,knm)+2.*q(3,i,j,knc))/4.
      enddo
      enddo
      enddo
      call q2div(q,ru)
      return
      end

c
c  ****************************** subrout oldqua **********************
c
c   transfer of old time mean quantities to ??old is performed only
c   after intia or inirea.
c
      subroutine oldqua
c
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/vmeao/vmo(ndv,m2),vrmso(6,m2)
      common/qmean/qm(ndv,m2),qrms(6,m2)
      common/wallst/cflw,cfuw
      common/skfl/ske(4,m2),fla(4,m2)
      common/prrm/prm(m2),prms(m2)
      common/skflo/skeo(4,m2),flao(4,m2)
      common/prrmo/prmo(m2),prmso(m2)
c
      cfo=cfn
      enavo=enav
      do 410 l=1,3
      do 410 j=1,n2m
      vmo(l,j)=qm(l,j)
  410 continue
      do 411 j=1,n2m
      prmo(j)=prm(j)
      prmso(j)=prms(j)
  411 continue
      do 413 l=1,4
      do 413 j=1,n2m
      skeo(l,j)=ske(l,j)
      flao(l,j)=fla(l,j)
      vrmso(l,j)=qrms(l,j)
  413 continue
      do j=1,n2m
       vrmso(5,j)=qrms(5,j)
       vrmso(6,j)=qrms(6,j)
      end do
      return
      end
c
c   *****************  subrout outh  ********************
c
c  writes output to file
c
      subroutine outh(ntime,time,cflm,enen)
c
      include 'param.f'
      parameter (m3m=m3-1)
      common/velmax/vmax(ndv),vmaxo(ndv)
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/tstep/dt,beta,ren
      common/wallst/cflw,cfuw
      common/veltot/vit(ndv)
      common/omegas/omr,omi
      common/eneto/enet,enstro
      common/timref/tref
      common/inener/ene0
      common/rhs3p/dp3ns
      common/neno/dpnew
c
c   write on file thist
c
      nfilt=32
      cflwp=-cflw
      write(6,158)ntime,time,vit(1),vit(2),vit(3)
     1         ,dp3ns,dpnew,enen,enet,cflm
      write(nfilt,159)time,vit(1),vit(2),vit(3)
     1         ,dp3ns,cfuw,cflw,enen,enet
  158 format(3x,i4,4(1x,e11.5),2x,8(1x,e11.5))
      write(17,159) time,cflwp,cfuw,dpnew,dp3ns
      write(18,159) time,enen,enet
      write(20,159) time,enstro
      write(19,159) time,dp3ns
  159 format(3x,11(1x,e11.5))
      return
      end
c
c  ****************************** subrout outpf  **********************
c
      subroutine outpf(time,enen,nav,q)
c
      include 'param.f'
      dimension q(ndv,m1,m2,m3)
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      common/tstep/dt,beta,ren
      dimension vprms(6),vmp(3),skp(4),flp(4),vistzr(m2)
      dimension psgs(6)
      common/filep/ifilp
      common/wallst/cflw,cfuw
      common/wallsto/cflwo,cfuwo
      common/vmean/vm(ndv,m2),vrms(6,m2)
      common/y2sta/y2s(m2)
      common/y13sta/y1s(m1),y3s(m3)
      common/corpri/yp1(m1),yp2(m2),yp3(m3)
      common/eneav/enav,enavo,cfn,cfo
      common/rhs3p/dp3ns
      common/qmean/qm(ndv,m2),qrms(6,m2)
      common/skfl/ske(4,m2),fla(4,m2)
      common/prrm/prm(m2),prms(m2)
      common/skflo/skeo(4,m2),flao(4,m2)
      common/prrmo/prmo(m2),prmso(m2)
      common/vmeao/vmo(ndv,m2),vrmso(6,m2)
      common/vismo/vismeo(m2),visrmo(m2)
      common/sgsp/sgs(6,m2),sgso(6,m2)
      common/names/filcnw,filcnr,filth,filou,filuu,filsf
      common/islwal/islv1s,islv1n,islv3s,islv3n
      common/ledat/csma,cvisc,iles
      real*4 lijmij,mijmij
      common/outd1/cdyn(m2),alij(m2),amij(m2)
      common/outd2/cs(m2),lijmij(m2),mijmij(m2)
      common/outd3/accdyn(m2),acclij(m2),accmij(m2)
      common/qles1/deltax1,deltay1,deltaz1,ell1c(m2)
c      
      character*27 filcnw,filcnr
      character*17 filth,filou,filuu,filsf
      character*17 filcos,filrms,filcfl,filvp,filrp,filtrs
      character*17 fildyn,filvis,filsgs
      character*4 pntim
      character*6 slp
      character*1 is1s,is3s
      character*1 is1n,is3n
c
      write(is1s,61) islv1s
      write(is1n,61) islv1n
      write(is3s,61) islv3s
      write(is3n,61) islv3n
   61 format(i1)
      slp='s'//is1s//is3s//'n'//is1n//is3n
      enavp=enavo/nav
      cfnp=cfo/nav
      uta=sqrt(abs(cfnp))
      cfm=(cflw-cfuw)*0.5
      utal=sqrt(abs(cflw))
      utau=sqrt(abs(cfuw))
      itim=time+0.3
      write(pntim,77) itim
   77 format(i4.4)
      filou='cpou.'//pntim
      filvp='cpvp.'//pntim
      filrp='cprp.'//pntim
      filsf='cpsf.'//pntim
      filrms='cprm.'//pntim
      filtrs='cpts.'//pntim
      filvis='cpvis.'//pntim
      open(39,file=filvis)
      if(iles.eq.1) then
      fildyn='cpdyn.'//pntim
      open(29,file=fildyn)
      fildyn='cpsgsn.'//pntim
      open(28,file=filsgs)
                    endif
      open(42,file=filou)
      open(52,file=filvp)
      open(53,file=filrp)
      open(43,file=filrms)
      open(62,file=filsf)
      open(63,file=filtrs)
      rewind(29)
      rewind(39)
      rewind(52)
      rewind(53)
      rewind(42)
      rewind(43)
      rewind(62)
      rewind(63)
c     write(52,783)time,utal,utau
      write(42,783)time,utal,utau,nav,cfnp,uta
c     write(43,783)time,utal,utau,nav,cfnp,uta
c     write(62,783)time,utal,utau,nav,cfnp,uta
  783 format(2x,e9.3,2x,e12.4,2x,e12.4,2x,i4,2x,2e14.7)
c     write(6,781)time,enavp,utal,utau,dp3ns
  781 format(3x,2x,'t=',e9.3,2x,'ener=',e10.4,2x,'utal=',e10.4,
     1  'utau=',e11.4,2x,'dp3ns=',e11.4
     1 /,3x,'j',8x,'y',12x,'u1',10x,'u2',10x,'u3',10x,'u1',10x,'u2',
     1  10x,'u3',8x,'u2u3',10x,'y+')
      vistzr(1)=vmo(3,2)/(y2s(2)-yp2(1))/nav/ren
      vistzr(n2)=-vmo(3,n2-1)/(yp2(n2)-y2s(n2m))/nav/ren
      do j=2,n2m
      vistzr(j)=(vmo(3,j)-vmo(3,j-1))/(y2s(j)-y2s(j-1))/nav/ren
      enddo
      do 611 j=1,n2m
      do l=1,6
      psgs(l)=sgso(l,j)/nav
      enddo
      do l=1,3
      skp(l)=skeo(l,j)/nav
      flp(l)=flao(l,j)/nav
      vmp(l)=vmo(l,j)/nav
      qrms(l,j)=sqrt(qrms(l,j))
      vprms(l)=sqrt(vrmso(l,j)/nav)
      enddo
      skp(4)=skeo(4,j)/nav
      flp(4)=flao(4,j)/nav
      vprms(4)=vrmso(4,j)/nav
      vprms(5)=vrmso(5,j)/nav
      vprms(6)=vrmso(6,j)/nav
      ppmp=prmso(j)  
      pvisme=vismeo(j)/nav
      pvisrm=visrmo(j)/nav
      toszrv=(vistzr(j)+vistzr(j+1))*0.5/uta/uta
      toszrt=vprms(4)/uta/uta
      toszrg=psgs(4)/uta/uta
      totszr=(toszrv-toszrt-toszrg)
      if(y2s(j).le.0) ypl=(1.+y2s(j))*ren*sqrt(abs(cflwo/nav))
      if(y2s(j).ge.0) ypl=(1.-y2s(j))*ren*sqrt(abs(cfuwo/nav))
      upl=vmp(3)/uta
      nsp1=vprms(1)/uta
      nsp2=vprms(2)/uta
      nsp3=vprms(3)/uta
      write(42,612)j,y2s(j),(vmp(l),l=1,3),(vprms(l),l=1,3),vprms(4)
     1            ,ypl,ppmp
      write(43,613)y2s(j),(vprms(l),l=1,6)
      write(39,613)y2s(j),pvisme,pvisrm    
      write(52,613)ypl,upl
      write(53,613)ypl,nsp1,nsp2,nsp3
      write(63,613)y2s(j),toszrv,toszrt,toszrg,totszr          
      write(62,613)y2s(j),(skp(l),l=1,4),(flp(l),l=1,4)
      if(iles.eq.1) then
      write(28,613)y2s(j),(psgs(l),l=1,6)
      write(29,613)y2s(j),cdyn(j),cs(j),sqrt(abs(cdyn(j)*ell1c(j))),
     1            lijmij(j),mijmij(j)
                    endif
  611 continue
  612 format(1x,i3,2x,e14.7,2x,9(1x,e15.7)) 
  613 format(1x,e14.7,2x,9(1x,e15.7)) 
      close(28)
      close(29)
      close(39)
      close(52)
      close(53)
      close(42)
      close(43)
      close(62)
      close(63)
      vmp3=vmo(3,n2m/2)/nav
      write(26,*)' nav =',nav,' vmo=',vmo(3,n2m/2),' vm=',vmp3,
     1             sqrt(abs(cflwo/nav)),sqrt(abs(cfuwo/nav))
      filcfl='cpcfl.'//pntim
      open(82,file=filcfl)
      dyl=(y2s(1)-yp2(1))
      dyu=-(y2s(n2m)-yp2(n2))
      do kc=1,n3m
      dql=0.
      dqu=0.
      do ic=1,n1m
      dql=q(3,ic,1,kc)/dyl+dql
      dqu=-q(3,ic,n2m,kc)/dyu+dqu
      enddo
      cfllk=dql/(ren*n1m)
      cfulk=dqu/(ren*n1m)
      write(82,613)y3s(kc),cfllk,cfulk
      enddo
      close(82)
      return
      end
c
c  ****************************** subrout tave0 **********************
c
c   calculation of time averaged mean quantities
c   veloc. and stresses
c
      subroutine tave0
c
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/vmeao/vmo(ndv,m2),vrmso(6,m2)
      common/eneav/enav,enavo,cfn,cfo
      common/wallsto/cflwo,cfuwo
      common/skflo/skeo(4,m2),flao(4,m2)
      common/prrmo/prmo(m2),prmso(m2)
      common/vismo/vismeo(m2),visrmo(m2)
      common/sgsp/sgs(6,m2),sgso(6,m2)
c
      enavo=0.
      cfo=0.
      cflwo=0.
      cfuwo=0.
      do l=1,3
      do j=1,n2m
      vmo(l,j)=0.
      enddo
      enddo
      do j=1,n2m
      prmo(j)=0.
      prmso(j)=0.
      enddo
      do l=1,4
      do j=1,n2m
      vrmso(l,j)=0.
      skeo(l,j)=0.
      flao(l,j)=0.
      enddo
      enddo
      do j=1,n2m
       vrmso(5,j)=0.
       vrmso(6,j)=0.
      end do
      do j=1,n2m
      vismeo(j)=0.
      visrmo(j)=0.
      end do
      do l=1,6
      do j=1,n2m
      sgso(l,j)=0.
      enddo
      enddo
      return
      end
c
c  ****************************** subrout taver **********************
c
c   calculation of time averaged mean quantities
c   veloc. and stresses
c
      subroutine taver
c
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/vmeao/vmo(ndv,m2),vrmso(6,m2)
      common/qmean/qm(ndv,m2),qrms(6,m2)
      common/eneav/enav,enavo,cfn,cfo
      common/wallst/cflw,cfuw
      common/wallsto/cflwo,cfuwo
      common/skfl/ske(4,m2),fla(4,m2)
      common/prrm/prm(m2),prms(m2)
      common/skflo/skeo(4,m2),flao(4,m2)
      common/prrmo/prmo(m2),prmso(m2)
      common/visme/vismed(m2),visrms(m2)
      common/vismo/vismeo(m2),visrmo(m2)
      common/sgsp/sgs(6,m2),sgso(6,m2)
      common/tauwal/cfnp
      common/rhs3p/dp3ns

c
c     cfm=(cflw-cfuw)*.5
      cfm=dp3ns
      enavo=enav+enavo
      cfo=cfm+cfo
      cfnp=cfo
      cflwo=cflwo+cflw
      cfuwo=cfuwo+cfuw
      do l=1,3
      do j=1,n2m
      vmo(l,j)=(vmo(l,j)+qm(l,j))
      enddo
      enddo
      do j=1,n2m
      prmo(j)=(prmo(j)+prm(j))
      prmso(j)=(prmso(j)+prms(j))
      enddo
      do l=1,4
      do j=1,n2m
      vrmso(l,j)=(vrmso(l,j)+qrms(l,j))
      skeo(l,j)=(skeo(l,j)+ske(l,j))
      flao(l,j)=(flao(l,j)+fla(l,j))
      enddo
      enddo
      do j=1,n2m
       vrmso(5,j)=(vrmso(5,j)+qrms(5,j))
       vrmso(6,j)=(vrmso(6,j)+qrms(6,j))
      end do
      do j=1,n2m
      vismeo(j)=vismeo(j)+vismed(j)
      visrmo(j)=visrmo(j)+visrms(j)
      end do
      do l=1,6
      do j=1,n2m
      sgso(l,j)=(sgso(l,j)+sgs(l,j))
      enddo
      enddo
      return
      end
c
c  ******************** subrout tavwri *********************
c
c   write quantities to evaluate averages in time
c
      subroutine tavwri(nav)
c
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/vmeao/vmo(ndv,m2),vrmso(6,m2)
      common/eneav/enav,enavo,cfn,cfo
      common/skflo/skeo(4,m2),flao(4,m2)
      common/prrmo/prmo(m2),prmso(m2)
      common/vismo/vismeo(m2),visrmo(m2)
      common/sgsp/sgs(6,m2),sgso(6,m2)
      character*10 tawrfi
      tawrfi='tavwri.out'
      open(70,file=tawrfi,form='unformatted')
      rewind(70)
      write(70)enavo,cfo,nav
      do 72 j=1,n2m
      write(70)(vmo(l,j),l=1,3),prmo(j),prmso(j),
     1     (vrmso(l,j),l=1,6),(skeo(l,j),l=1,4),(flao(l,j),l=1,4)
     1     ,(sgso(l,j),l=1,6),vismeo(j),visrmo(j)
   72 continue
      close(70)
      return
      end
c
c  ******************** subrout tavrea *********************
c
c   read quantities to evaluate averages in time
c
      subroutine tavrea(nav)
c
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/vmeao/vmo(ndv,m2),vrmso(6,m2)
      common/eneav/enav,enavo,cfn,cfo
      common/skflo/skeo(4,m2),flao(4,m2)
      common/prrmo/prmo(m2),prmso(m2)
      common/vismo/vismeo(m2),visrmo(m2)
      common/sgsp/sgs(6,m2),sgso(6,m2)
      character*10 tawrfi
      tawrfi='tavwri.out'
      open(70,file=tawrfi,form='unformatted')
      read(70)enavo,cfo,nav
      do 72 j=1,n2m
      read(70)(vmo(l,j),l=1,3),prmo(j),prmso(j),
     1     (vrmso(l,j),l=1,6),(skeo(l,j),l=1,4),(flao(l,j),l=1,4)
     1     ,(sgso(l,j),l=1,6),vismeo(j),visrmo(j)
   72 continue
      close(70)
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
      vmax(l)=max(vm,vfm)
  310 continue
  311 continue
      return
      end
c
c  ****************************** subrout wstre **********************
c
      subroutine wstre(q)
c
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      common/meshu/d1x1,d1x2,d1x3
      dimension q(ndv,m1,m2,m3)
      common/wallst/cflw,cfuw
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/tstep/dt,beta,re
      common/d13/alx1,alx3
      common/neno/dpnew
      common/y2sta/y2s(m2)
      common/corpri/yp1(m1),yp2(m2),yp3(m3)
c
c     lower and upper walls shear
ccccccccccccccccccccccccccccccccccccccccccccccccc
      dql=0.
      dqu=0.
      dyl=(y2s(1)-yp2(1))
      dyu=-(y2s(n2m)-yp2(n2))
      do 450 kc=1,n3m
      do 450 ic=1,n1m
      dql=(q(3,ic,1,kc)+q(3,ic,1,kpv(kc)))*.5/dyl+dql
      dqu=-(q(3,ic,n2m,kc)+q(3,ic,n2m,kpv(kc)))*.5/dyu+dqu
  450 continue
      cfl1=dql/(re*n1m*n3m)
      cfu1=dqu/(re*n1m*n3m)
      cfuw=cfu1
      cflw=cfl1
      dpnew=.5*(cfuw-cflw)
c     print *,'cflold=',cflw,'cfuold=',cfuw
      return
      end
