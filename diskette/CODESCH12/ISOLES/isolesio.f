c
c  ****************************** subrout contwr **********************
c
c   the field is writen for restarting or for visualizations
c
      subroutine contwr(ntime,time,q,rho,pr,enen)
c
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      dimension pr(m1,m2,m3)
      dimension q(ndv,m1,m2,m3)
      common/velmax/vmax(ndv),vmaxo(ndv)
      dimension rho(m1,m2,m3)
      common/tstep/dt,beta,ren
      common/wallst/cflw,cfuw
      common/veltot/vit(ndv)
      common/omegas/omr,omi
      common/timref/tref
      common/inener/ene0
      common/rhs3p/dp3ns
      common/eneav/enav,enavo
      common/strat1/istrat,rho0,g,schm
      common/names/filcnw,filcnr,filth,filou,filuu,filsf
      character*17 filcnw,filcnr,filth,filou,filuu,filsf
      open(13,file=filcnw,form='unformatted')
c
      nfil=13
      rewind(nfil)
      write(nfil) n1,n2,n3
      write(nfil) time,ren,time,time
      if(istrat.eq.0) then
      write(nfil) (((q(1,i,j,k),i=1,n1),j=1,n2),k=1,n3),
     1            (((q(2,i,j,k),i=1,n1),j=1,n2),k=1,n3),
     1            (((q(3,i,j,k),i=1,n1),j=1,n2),k=1,n3)
     1           ,(((pr(i,j,k),i=1,n1),j=1,n2),k=1,n3)
                      else
      write(nfil) (((q(1,i,j,k),i=1,n1),j=1,n2),k=1,n3),
     1            (((q(2,i,j,k),i=1,n1),j=1,n2),k=1,n3),
     1            (((q(3,i,j,k),i=1,n1),j=1,n2),k=1,n3),
     1            (((pr(i,j,k),i=1,n1),j=1,n2),k=1,n3),
     1            (((rho(i,j,k),i=1,n1),j=1,n2),k=1,n3)
                      endif
      close(nfil)
      return
      end
c
c  ****************************** subrout enerca **********************
c
c   calculation of total energy = (v1**2+v2**2+v3**2)*0.5
c   vi cartesian components
c
      subroutine enerca(q,pr,enej,time)
c
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      dimension q(ndv,m1,m2,m3),pr(m1,m2,m3)
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      common/meshu/d1x1,d1x2,d1x3
      common/qmean/qm(ndv),qrms(4)
      common/d1/re
      common/d123/alx1,alx2,alx3
      common/veltot/vit(ndv)
      common/eneav/enav,enavo
      common/vmean/vm(ndv),vrms(6)
      common/skfl/ske(4),fla(4)
      common/prrm/prm,prms
      common/skflo/skeo(4),flao(4)
      common/prrmo/prmo,prmso
      common/vormo/vo1rms,vo2rms,vo3rms
      common/skder/skd1,skd2,skd3,fld1,fld2,fld3
      common/strain/st(m1,m2,m3,6)
      common/ledat/ics0,cvisc,ifiltr,ibox
      common/sma/vis(m1,m2,m3)
      common/smapra/dif(m1,m2,m3)
      common/vista/vimed,virms,viske,vifla
      common/dista/dimed,dirms,diske,difla
      common/enssta/ensmed,ensrms,ensske,ensfla
      common/kolqua/diss,eta
      common/tstep/dt,beta,ren
      common/sij2po/sijji
      common/oooo/enstro,diss1,visma,vismi
      common/rhsc/rhs(m1,m2,m3)
c
      pi=2.*asin(1.)
      vol=(2.*pi)**3.
      vl13=1./float(n1m*n3m)
      vl123=1./float(n1m*n2m*n3m)
c     vit(1)=0.
c     vit(2)=0.
c     vit(3)=0.
      vm1m=0.
      vm2m=0.
      vm3m=0.
      prm=0.
      vimed=0.
      dimed=0.
      do 411 k=1,n3m
      kp=kpv(k)
      do 411 j=1,n2m
      jp=jpv(j)
      do 411 i=1,n1m
      ip=ipv(i)
       vm1m=q(1,i,j,k)+vm1m
       vm2m=q(2,i,j,k)+vm2m
       vm3m=q(3,i,j,k)+vm3m
      prm=prm+pr(i,j,k) 
      vimed=vimed+vis(i,j,k) 
      dimed=dimed+dif(i,j,k) 
c     vit(1)=q(1,i,j,k)+vit(1)
c     vit(2)=q(2,i,j,k)+vit(2)
c     vit(3)=q(3,i,j,k)+vit(3)
  411 continue
      vm(1)=vm1m*vl123
      vm(2)=vm2m*vl123
      vm(3)=vm3m*vl123
      prm=prm*vl123
      vimed=vimed*vl123
      dimed=dimed*vl123
  410 continue
c     vit(1)=vit(1)*vl123/vol
c     vit(2)=vit(2)*vl123/vol
c     vit(3)=vit(3)*vl123/vol
c
c    here the velocity rms 
c    and the stresses are evaluated
c    even the rms of the eddy viscosity and
c    eddy diffusivity are calculated
c
      enej=0.
      u1rmm=0.
      u2rmm=0.
      u3rmm=0.
      u23mm=0.
      u31mm=0.
      u12mm=0.
      ppmm=0.
      virms=0.
      dirms=0.
      do 414 k=1,n3m
      kp=kpv(k)
      do 414 j=1,n2m
      jp=jpv(j)
      do 414 i=1,n1m
      ip=ipv(i)
       v1m=q(1,i,j,k)-vm(1)
       v2m=q(2,i,j,k)-vm(2)
       v3m=q(3,i,j,k)-vm(3)
       v1c=(q(1,i,j,k)+q(1,ip,j,k))*0.5
       v2c=(q(2,i,j,k)+q(2,i,jp,k))*0.5
       v3c=(q(3,i,j,k)+q(3,i,j,kp))*0.5
      ppm=pr(i,j,k)-prm
      vipm=vis(i,j,k)-vimed
      dipm=dif(i,j,k)-dimed
      u1rmm=u1rmm+v1m**2
      u2rmm=u2rmm+v2m**2
      u3rmm=u3rmm+v3m**2
      u23mm=u23mm+v2m*v3m
      u12mm=u12mm+v1c*v2c
      u31mm=u31mm+v1c*v3c
      ppmm=ppmm+ppm**2
      virms=virms+vipm**2
      dirms=dirms+dipm**2
      enejp=(v1m**2+v2m**2+v3m**2)*0.5*vl123
      enej=enej+enejp
  414 continue
      vrms(1)=u1rmm*vl123
      vrms(2)=u2rmm*vl123
      vrms(3)=u3rmm*vl123
      vrms(4)=u23mm*vl123
      vrms(5)=u12mm*vl123
      vrms(6)=u31mm*vl123
      prms=ppmm*vl123
      virms=virms*vl123
      dirms=dirms*vl123
      if(time.eq.0.) virms=1.
      if(time.eq.0.) dirms=1.
c
c
c    here the skewness and flatness of the fluctuating velocities
c    are evaluated
c    even the skewness and flatness  of the eddy viscosity and
c    eddy diffusivity are calculated
c
      sk1m=0.
      sk2m=0.
      sk3m=0.
      sk4m=0.
      fl1m=0.
      fl2m=0.
      fl3m=0.
      fl4m=0.
      viske=0.
      vifla=0.
      diske=0.
      difla=0.
      do 514 k=1,n3m
      kp=kpv(k)
      do 514 j=1,n2m
      jp=jpv(j)
      do 514 i=1,n1m
      ip=ipv(i)
       v1m=q(1,i,j,k)-vm(1)
       v2m=q(2,i,j,k)-vm(2)
       v3m=q(3,i,j,k)-vm(3)
      ppm=pr(i,j,k)-prm
      vipm=vis(i,j,k)-vimed
      dipm=dif(i,j,k)-dimed
      sk1m=sk1m+(v1m/sqrt(vrms(1)))**3
      sk2m=sk2m+(v2m/sqrt(vrms(2)))**3
      sk3m=sk3m+(v3m/sqrt(vrms(3)))**3
      sk4m=sk4m+(ppm/sqrt(prms))**3
      fl1m=fl1m+(v1m/sqrt(vrms(1)))**4
      fl2m=fl2m+(v2m/sqrt(vrms(2)))**4
      fl3m=fl3m+(v3m/sqrt(vrms(3)))**4
      fl4m=fl4m+(ppm/sqrt(prms))**4
      viske=viske+(vipm/sqrt(virms))**3
      vifla=vifla+(vipm/sqrt(virms))**4
      diske=diske+(dipm/sqrt(dirms))**3
      difla=difla+(dipm/sqrt(dirms))**4
  514 continue
      ske(1)=sk1m*vl123
      ske(2)=sk2m*vl123
      ske(3)=sk3m*vl123
      ske(4)=sk4m*vl123
      fla(1)=fl1m*vl123   
      fla(2)=fl2m*vl123   
      fla(3)=fl3m*vl123   
      fla(4)=fl4m*vl123  
      viske=viske*vl123
      vifla=vifla*vl123
      diske=diske*vl123
      difla=difla*vl123
  513 continue
      enav=enej
      do 428 l=1,3
      qrms(l)=vrms(l)
      qm(l)=vm(l)
  428 continue
      qrms(4)=vrms(4)
c     write(*,*) 'vm=',vm(1),vm(2),vm(3),prm
c     write(*,*) 'vrms=',vrms(1),vrms(2),vrms(3),vrms(4)
c     write(*,*) 'ske=',ske(1),ske(2),ske(3),ske(4)
c     write(*,*) 'fla=',fla(1),fla(2),fla(3),fla(4)
c
c     skewness della derivata della velocita'
c
      rmsdux=0.
      rmsduy=0.
      rmsduz=0.
      do 714 k=1,n3m
      kp=kpv(k)
      do 714 j=1,n2m
      jp=jpv(j)
      do 714 i=1,n1m
      ip=ipv(i)
      dux=(q(1,ip,j,k)-q(1,i,j,k))*dx1
      duy=(q(2,i,jp,k)-q(2,i,j,k))*dx2
      duz=(q(3,i,j,kp)-q(3,i,j,k))*dx3
      rmsdux=rmsdux+dux*dux
      rmsduy=rmsduy+duy*duy
      rmsduz=rmsduz+duz*duz
 714  continue
      rmsdux=rmsdux*vl123
      rmsduy=rmsduy*vl123
      rmsduz=rmsduz*vl123
c
c     write(*,*) 'rmsder=',rmsdux,rmsduy,rmsduz
c
      skd1=0.
      skd2=0.
      skd3=0.
      fld1=0.
      fld2=0.
      fld3=0.
      do 724 k=1,n3m
      kp=kpv(k)
      do 724 j=1,n2m
      jp=jpv(j)
      do 724 i=1,n1m
      ip=ipv(i)
      dux=(q(1,ip,j,k)-q(1,i,j,k))*dx1
      duy=(q(2,i,jp,k)-q(2,i,j,k))*dx2
      duz=(q(3,i,j,kp)-q(3,i,j,k))*dx3
      skd1=skd1+(dux/sqrt(rmsdux))**3
      skd2=skd2+(duy/sqrt(rmsduy))**3
      skd3=skd3+(duz/sqrt(rmsduz))**3
      fld1=fld1+(dux/sqrt(rmsdux))**4
      fld2=fld2+(duy/sqrt(rmsduy))**4
      fld3=fld3+(duz/sqrt(rmsduz))**4
 724  continue
      skd1=skd1*vl123
      skd2=skd2*vl123
      skd3=skd3*vl123
      fld1=fld1*vl123
      fld2=fld2*vl123
      fld3=fld3*vl123
c     write(*,*) 'skeder=',skd1,skd2,skd3
c     write(*,*) 'flader=',fld1,fld2,fld3
c
c     compute the dissiapation 
c
      call strper(q)
      diss=0.
      diss1=0.
      visma=-1000.
      vismi=+1000.
      do k=1,n3m
       do j=1,n2m
        do i=1,n1m
        app1=2.*
     1   (st(i,j,k,1)*st(i,j,k,1)+
     1    st(i,j,k,2)*st(i,j,k,2)+
     1    st(i,j,k,3)*st(i,j,k,3)+2.*
     1   (st(i,j,k,4)*st(i,j,k,4)+
     1    st(i,j,k,5)*st(i,j,k,5)+
     1    st(i,j,k,6)*st(i,j,k,6)))*vl123
        diss=diss+app1*cvisc
        diss1=diss1+app1*(vis(i,j,k)-cvisc)
        visma=max(visma,(vis(i,j,k)-cvisc))
        vismi=min(vismi,(vis(i,j,k)-cvisc))
        end do
       end do
      end do
c
c               compute the  enstrophy
c                                                                       
      vo1rms=0.
      vo2rms=0.
      vo3rms=0.
      enstro=0.
      do kc=1,n3m 
      km=kmv(kc)                                                    
      do jc=1,n2m                                                     
      jm=jmv(jc)                                                        
      do ic=1,n1m                                                     
      im=imv(ic)                                                        
c
c      OM_x COMPONENT
c
      omx=+(q(3,ic,jc,kc)-q(3,ic,jm,kc))*dx2
     1    -(q(2,ic,jc,kc)-q(2,ic,jc,km))*dx3
c
c      OM_y COMPONENT
c
      omy=+(q(1,ic,jc,kc)-q(1,ic,jc,km))*dx3
     1    -(q(3,ic,jc,kc)-q(3,im,jc,kc))*dx1
c
c      OM_z COMPONENT
c
      omz=-(q(1,ic,jc,kc)-q(1,ic,jm,kc))*dx2
     1    +(q(2,ic,jc,kc)-q(2,im,jc,kc))*dx1
      rhs(ic,jc,kc)=(omx*omx+omy*omy+omz*omz)*0.5
      vo1rms=vo1rms+omx*omx
      vo2rms=vo2rms+omy*omy
      vo3rms=vo3rms+omz*omz
      enstro=enstro+(omx*omx+omy*omy+omz*omz)
      enddo
      enddo
      enddo
      vo1rms=vo1rms*vl123
      vo2rms=vo2rms*vl123
      vo3rms=vo3rms*vl123
      enstro=enstro*vl123*0.5
      ensmed=enstro
      eta=(cvisc**3./diss)**(1./4.)
      ensrms=0.
      ensske=0.
      ensfla=0.
      do kc=1,n3m 
      do jc=1,n2m                                                     
      do ic=1,n1m                                                     
      enspm=rhs(ic,jc,kc)-ensmed
      ensrms=ensrms+enspm**2
      enddo
      enddo
      enddo
      ensrms=ensrms*vl123
      if(ensrms.eq.0.) ensrms=1.
      do kc=1,n3m 
      do jc=1,n2m                                                     
      do ic=1,n1m                                                     
      enspm=rhs(ic,jc,kc)-ensmed
      ensske=ensske+(enspm/sqrt(ensrms))**3
      ensfla=ensfla+(enspm/sqrt(ensrms))**4
      ensrms=ensrms
      enddo
      enddo
      enddo
      return
      end
c
c  ****************************** subrout inirea **********************
c
c   reading initial  conditions from file
c   when calculation are continued from a previous
c   calculation.
c
      subroutine inirea(ntii,time,q,rho,pr)
c
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      dimension pr(m1,m2,m3)
      dimension q(ndv,m1,m2,m3)
      common/velmax/vmax(ndv),vmaxo(ndv)
      dimension rho(m1,m2,m3)
      common/tstep/dt,beta,ren
      common/wallst/cflw,cfuw
      common/veltot/vit(ndv)
      common/strat1/istrat,rho0,g,schm
      common/omegas/omr,omi
      common/timref/tref
      common/inener/ene0
      common/eneav/enav,enavo
      common/inior/indrea,icosma,icont
      common/d123/alx1,alx2,alx3
c
      nfil=23
      rewind(nfil)
      read(nfil) n1,n2,n3
      read(nfil) time,ren,time,dum
      if(istrat.eq.0) then
      read(nfil) (((q(1,i,j,k),i=1,n1),j=1,n2),k=1,n3),
     1            (((q(2,i,j,k),i=1,n1),j=1,n2),k=1,n3),
     1            (((q(3,i,j,k),i=1,n1),j=1,n2),k=1,n3)
     1           ,(((pr(i,j,k),i=1,n1),j=1,n2),k=1,n3)
                      else
      read(nfil) (((q(1,i,j,k),i=1,n1),j=1,n2),k=1,n3),
     1            (((q(2,i,j,k),i=1,n1),j=1,n2),k=1,n3),
     1            (((q(3,i,j,k),i=1,n1),j=1,n2),k=1,n3),
     1            (((pr(i,j,k),i=1,n1),j=1,n2),k=1,n3),
     1            (((rho(i,j,k),i=1,n1),j=1,n2),k=1,n3)
                      endif
      close(nfil)
      return
      end
c
c   *****************  subrout outh  ********************
c
c  writes time histories of some global quantity
c
      subroutine outh(ntime,time,cflm,enen)
c
      include 'param.f'
      parameter (m3m=m3-1)
      common/velmax/vmax(ndv),vmaxo(ndv)
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/tstep/dt,beta,ren
      common/veltot/vit(ndv)
      common/omegas/omr,omi
      common/timref/tref
      common/inener/ene0
      common/vmean/vm(ndv),vrms(6)
      common/skfl/ske(4),fla(4)
      common/prrm/prm,prms
      common/vormo/vo1rms,vo2rms,vo3rms
      common/skder/skd1,skd2,skd3,fld1,fld2,fld3
      common/vista/vimed,virms,viske,vifla
      common/dista/dimed,dirms,diske,difla
      common/enssta/ensmed,ensrms,ensske,ensfla
      common/kolqua/diss,eta
      common/ledat/ics0,cvisc,ifiltr,ibox
      common/resca/iresca,tresca,rlamas
      common/dycopr/cdynav,visad,avg,avmij,nif,cist,flback,flint
      common/icomp/irunpc
      common/sij2po/sijji
      common/oooo/enstro,diss1,visma,vismi
      common/cosma/csma,pratu
      common/csmac/cosmac
      common/qwolm/sij2,sij2fi,raddis,denvis,anuvis,
     1             dendif,anudif,visav,difav
      common/names/filcnw,filcnr,filth,filou,filuu,filsf
      character*17 filcnw,filcnr,filth,filou,filuu,filsf
c
c   write on file thist
c
      if(irunpc.eq.0) then
      call openfi(time)
                      endif
      rlamb=ren*enen*sqrt(20./(3.*sijji))
      sijsij=diss/cvisc
      write(34,135)time,enen,enstro,sijji,sijsij,diss,rlamb
      write(35,*)time,eta,diss,diss1,visma,vismi
      nfilt=32
      write(6,158)ntime,time,enen,diss,rlamb,eta
     1           ,(vrms(1)+vrms(2)+vrms(3))/3.
     1           ,(ske(1)+ske(2)+ske(3))/3.
     1           ,(fla(1)+fla(2)+fla(3))/3.
     1           ,(skd1+skd2+skd3)/3.,(fld1+fld2+fld3)/3.
      write(nfilt,158)ntime,time,enen,vmax(1),vmax(2),vmax(3)
     1           ,(skd1+skd2+skd3)/3.,(fld1+fld2+fld3)/3.
  158 format(3x,i4,10(1x,e11.5))
      write(42,258)time,vrms(1),vrms(2),vrms(3)
      write(52,258)time,ske(1),ske(2),ske(3),ske(4)
      write(62,258)time,fla(1),fla(2),fla(3),fla(4)
      write(72,258)time,skd1,skd2,skd3,fld1,fld2,fld3
      write(82,258)time,vo1rms,vo2rms,vo3rms         
      write(92,258)time,vrms(4),vrms(5),vrms(6)
      write(93,258)time,vimed,virms,viske,vifla
      write(95,258)time,dimed,dirms,diske,difla
      write(94,258)time,ensmed,ensrms,ensske,ensfla
  258 format(7(e13.6,1x))
       if(ics0.eq.1) then
       if(ibox.eq.1) then
       write(33,134) time,sqrt(cdynav),visad,-avg,avmij
     1           ,nif
  134 format(5e12.5,3x,i5)
       else
       write(33,135) time,visad,cist,-avg,avmij
     1           ,flback,flinf
  135 format(8e12.5)
       end if
                      endif
       if(ics0.eq.2) then
      if(time.gt.0.05) then
      write(33,161)time,sij2,sij2fi,raddis,visav,difav,cosmac,pratu
                      endif
  161 format(11e12.4)
                     endif
      if(irunpc.eq.0) then
      call closefi
                      endif
      return
      end
c
c  ****************************** subrout openfi  **********************
c
      subroutine openfi(time)
      common/resca/iresca,tresca,rlamas
      if(iresca.eq.1.and.time.le.tresca) then
      open(32,file='isolesth.tran')
      open(34,file='enstr.tran')
      open(35,file='visctu.tran')
      open(42,file='rms.tran')
      open(52,file='ske.tran')
      open(62,file='fla.tran')
      open(82,file='vorms.tran')
      open(92,file='tarms.tran')
      open(93,file='virms.tran')
      open(95,file='dirms.tran')
      open(94,file='ensrms.tran')
      open(72,file='momder.tran')
      open(33,file='dyncost.tra')
                                         else
      open(32,file='isolesth.out')
      open(34,file='enstr.out')
      open(35,file='visctu.out')
      open(42,file='rms.out')
      open(52,file='ske.out')
      open(62,file='fla.out')
      open(82,file='vorms.out')
      open(92,file='tarms.out')
      open(93,file='virms.out')
      open(95,file='dirms.out')
      open(94,file='ensrms.out')
      open(72,file='momder.out')
      open(33,file='dyncost.out')
                                         endif
      return
      end
c
c  ****************************** subrout closefi  **********************
c
      subroutine closefi
      close(34)
      close(35)
      close(33)
      close(32)
      close(42)
      close(52)
      close(62)
      close(72)
      close(82)
      close(92)
      close(93)
      close(94)
      close(95)
      return
      end
c
c  ****************************** subrout outpf  **********************
c
      subroutine outpf(time,enen,nav,q,qcap)
c
c    here writes the spectra of velocity density and
c    also the spectra kolmogorov scaled 
c
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      dimension q(ndv,m1,m2,m3),qcap(m1,m2,m3)
      common/tstep/dt,beta,ren
      dimension vprms(4),vmp(3),skp(4),flp(4)
      common/filep/ifilp
      common/vmean/vm(ndv),vrms(6)
      common/eneav/enav,enavo
      common/qmean/qm(ndv),qrms(4)
      common/skfl/ske(4),fla(4)
      common/prrm/prm,prms
      common/skflo/skeo(4),flao(4)
      common/prrmo/prmo,prmso
      common/vmeao/vmo(ndv),vrmso(4)
      common/names/filcnw,filcnr,filth,filou,filuu,filsf
      common/outd1/cdyn(m2),ell1(m2),alij(m2),amij(m2)
      common/outd2/cs(m2),flmij(m2),fmmij(m2)
      common/corpri/yp1(m1),yp2(m2),yp3(m3)
      common/ledat/ics0,cvisc,ifiltr,ibox
      common/speene/e(ndv,0:m3)
      common/viscsp/visk(0:m3)
      common/numond/kx(m1),ky(m2),kz(m3),kkmax
      common/kolqua/diss,eta
      common/kpmasp/kmax
      common/mtime/multim
      common/teddy/akedd,eeto
      common/strat2/igrad,bvais
      common/strat1/istrat,rho0,g,schm   
      common/spemw/akkpp,qq,sig
      common/resca/iresca,tresca,rlamas
      common/ispec/imic
      common/stratf/rich
c      
      character*17 filcnw,filcnr,filth,filou,filuu,filsf
      character*17 filcos
      character*4 pntim
      character*6 slp
      character*1 is1s,is3s
      character*1 is1n,is3n
c
      enavp=enavo/nav
      itim=(multim*time+0.3)
      write(pntim,77) itim
   77 format(i4.4)
  611 continue
  612 format(1x,i3,2x,e14.7,2x,9(1x,e11.4)) 
  613 format(e14.7,1x,9(1x,e11.4)) 
      if(iresca.eq.1.and.time.le.tresca) then
      filcos='spec.'//pntim//'tra'
      open(57,file=filcos)
      rewind(57)
      filcos='spto.'//pntim//'tra'
      open(58,file=filcos)
      rewind(58)
      filcos='spkoto.'//pntim//'tra'
      open(56,file=filcos)
      filcos='speko.'//pntim//'tra'
      open(59,file=filcos)
      rewind(59)
      call spectre(q,qcap)
       if(ics0.ne.0) then
       filcos='edvdyn.'//pntim//'tra'
                     endif
                                         else
      filcos='spec.'//pntim
      open(57,file=filcos)
      rewind(57)
      filcos='spto.'//pntim
      open(58,file=filcos)
      rewind(58)
      filcos='spkoto.'//pntim
      open(56,file=filcos)
      filcos='speko.'//pntim
      open(59,file=filcos)
      rewind(59)
      call spectre(q,qcap)
       if(ics0.ne.0) then
       filcos='edvdyn.'//pntim
                     endif
                                        endif
      speska=(1./(cvisc**5.*diss))**(1./4.)
       akrap=1.
       if(igrad.eq.1.and.istrat.eq.1) akrap=akkpp
       ee1to=0.
       ee2to=0.
       ee3to=0.
       disk=0.
       do k=1,kmax
       ak=k/akrap
       write (58,*)ak,e(1,k)+e(2,k)+e(3,k)
       write (57,*)k,e(1,k),e(2,k),e(3,k)
       ee1to=ee1to+e(1,k)
       ee2to=ee2to+e(2,k)
       ee3to=ee3to+e(3,k)
       eeto=ee1to+ee2to+ee3to
       disk=disk+(e(1,k)+e(2,k)+e(3,k))*k**2
       akol=k*eta
       se1=e(1,k)*speska
       se2=e(2,k)*speska
       se3=e(3,k)*speska
       write (56,*)akol,se1+se2+se3
       write (59,*)akol,se1,se2,se3
       end do
      close(58)
      close(57)
      close(56)
      close(59)
      vsca=sqrt(2.*eeto)
c     ted=time/(vsca*akedd)
      ted=time
      if(imic.eq.10) then
      ted=time/sqrt(rich)
                    endif
      sijsij=diss/cvisc
      write(6,133)ted,ee1to,ee2to,ee3to,eeto,enen,disk,sijsij
  133 format('te=',e12.4,3x,3e12.4,3x,'En sp,phy',2e12.5
     1       ,'Dis sp,ph',2e12.5)
       if(ics0.ne.0) then
       call spevis(qcap)
       filcos='edvdyn.'//pntim
       open(29,file=filcos,form='formatted')
       do kk=0,kmax
       akk=kk/akrap
       ekk=(e(1,kk)+e(2,kk)+e(3,kk))
       write(29,*) akk,ekk,visk(kk)
        end do
                     endif
      close(29)
      return
      end
      subroutine spevis(qtil)
c
c    here the spectrum of the eddy viscosity is calculated
c
c     energy spectrum
      include 'param.f'
      common/sma/vis(m1,m2,m3)
      complex qtil(m1,(m2-1)/2+1,m3)
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/fft001/ifxx1(13),trigxx1(2*(m1-1))
      common/fft002/ifxx2(13),trigxx2(3*(m2-1)/2+1)
      common/fft003/ifxx3(13),trigxx3(2*(m3-1))
      common/numond/kx(m1),ky(m2),kz(m3),kkmax
      common/viscsp/visk(0:m3)
      common/spequa/xa,wor,xa2,wor2,xr,work
      complex xa(m3-1,m1-1),wor(m3-1,m1-1)
      complex xa2(m1-1,m3-1),wor2(m1-1,m3-1)
      real xr(m2+1,m1-1),work(m2,m1-1)
c
c     print *,'calcolo dello spettro viscosita'
      n2mh=n2m/2+1
c
      do  10 k=1,n3m
c
        do i=1,n1m
         xr(1,i)=vis(i,n2m,k)
         xr(n2m+2,i)=vis(i,1,k)
         do j=1,n2m
          js=j+1
          xr(js,i)=vis(i,j,k)
         enddo
        enddo
c
c   real fft applied from  physical space to wave number
c   in x2 direction
c
        call fft99(xr,work,trigxx2,ifxx2,1,m2+1,n2m,n1m,-1)
c
        do j=1,n2mh
         jp=2*j
         jd=2*j-1
         do i=1,n1m
          qtil(i,j,k)=cmplx(xr(jd,i),xr(jp,i))
         enddo
        enddo
c
 10   continue
c
c   2-d  cfft applied (twice) from
c   physical space to wave number
c
      do 20 j=1,n2mh
c
        do k=1,n3m
         do i=1,n1m
          xa(k,i)=qtil(i,j,k)
         enddo
        enddo
c
        call cfft99(xa,wor,trigxx3,ifxx3,1,m3-1,n3m,n1m,-1)
c
        do k=1,n3m
         do i=1,n1m
          xa2(i,k)=xa(k,i)/float(n3m)
         enddo
        enddo
c
        call cfft99(xa2,wor2,trigxx1,ifxx1,1,m1-1,n1m,n3m,-1)
c
        do i=1,n1m
         do k=1,n3m
          qtil(i,j,k)=xa2(i,k)/float(n1m)
         enddo
        enddo
c
  20   continue
       do k=0,kkmax
       visk(k)=0.
       end do
c
c     qtil is the velocity component in Fourier space.
c
       do k=1,n3m
        do j=1,n2mh
         do i=1,n1m
         akk=sqrt(float(kx(i)*kx(i)+ky(j)*ky(j)+kz(k)*kz(k)))
         kk=nint(akk)
         urea=real(qtil(i,j,k))
         uimm=aimag(qtil(i,j,k))
         visk(kk)=visk(kk)+(urea*urea+uimm*uimm)
         end do
        end do
       end do
  1    continue
       return
       end

