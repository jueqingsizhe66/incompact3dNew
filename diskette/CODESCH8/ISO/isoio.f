c
c  ****************************** subrout contwr **********************
c
c   the velocity field is saved for postprocessing or as a
c   restarting file
c   initial conditions when the calculation start from a previous
c   calculation
c   only two velocity components are necessary together 
c   with the third component in a single plane.
c   From this the third component is evaluated in the whole
c   box by the divergence.
c
      subroutine contwr(ntime,time,q,enen)
c
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      dimension q(ndv,m1,m2,m3)
      common/velmax/vmax(ndv),vmaxo(ndv)
      common/tstep/dt,beta,ren
      common/wallst/cflw,cfuw
      common/veltot/vit(ndv)
      common/omegas/omr,omi
      common/timref/tref
      common/inener/ene0
      common/rhs3p/dp3ns
      common/eneav/enav,enavo
      common/names/filcnw,filcnr,filth,filou,filuu,filsf
      character*17 filcnw,filcnr,filth,filou,filuu,filsf
      character*17 filcos
      character*4 pntim
      common/ispec/imic
      common/cflco/icfl,cflc,tpin,tprin,tfin,twrfi
      if(imic.ge.0) then
      if(amod(time,twrfi).lt.dt) then
      itim=time+0.3
      write(pntim,77) itim
      filcos='field3d'//pntim//'.dat'
      open(13,file=filcos,form='unformatted',status='unknown')
                                 else
      open(13,file=filcnw,form='unformatted',status='unknown')
                                 endif
c
                    else
      itim=time+0.3
      write(pntim,77) itim
   77 format(i4.4)
      filcos='field3d'//pntim//'.dat'
      open(13,file=filcos,form='unformatted',status='unknown')
                    endif
      nfil=13
      rewind(nfil)
      write(nfil) n1,n2,n3
      write(nfil) time,ren,time,time
      write(nfil) (((q(1,i,j,k),i=1,n1),j=1,n2),k=1,n3),
     1            (((q(2,i,j,k),i=1,n1),j=1,n2),k=1,n3),
     1            ((q(3,i,j,1),i=1,n1),j=1,n2)
      close(nfil)
      return
      end
c
c  ****************************** subrout enerca **********************
c
c   calculation of total energy = (v1**2+v2**2+v3**2)*0.5
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
      common/skder/skd1,skd2,skd3,fld1,fld2,fld3
      common/ledat/cvisc
c
      pi=2.*asin(1.)
      vol=(2.*pi)**3.
      vl123=1./float(n1m*n2m*n3m)
      vm1m=0.
      vm2m=0.
      vm3m=0.
      prm=0.
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
  411 continue
      vm(1)=vm1m*vl123
      vm(2)=vm2m*vl123
      vm(3)=vm3m*vl123
      prm=prm*vl123
  410 continue
c
c   the mean value has been evaluated 
c   it should be zero in isotropic turbulence
c   if at t=0 the mean velocities were zero
c
c
c  the fluctuating quantities are calculted at the cell
c  centre together with the rms values of the one-point
c  correlations
c
      enej=0.
      u1rmm=0.
      u2rmm=0.
      u3rmm=0.
      u31mm=0.
      u12mm=0.
      u23mm=0.
      ppmm=0.
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
      u1rmm=u1rmm+v1m**2
      u2rmm=u2rmm+v2m**2
      u3rmm=u3rmm+v3m**2
      u23mm=u23mm+v2c*v3c
      u12mm=u12mm+v1c*v2c
      u31mm=u31mm+v1c*v3c
      ppmm=ppmm+ppm**2
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
c
c   here skewness and flatness of V and p are evaluated
c
      if(ppm.eq.0.) prms=1.
      sk1m=0.
      sk2m=0.
      sk3m=0.
      sk4m=0.
      fl1m=0.
      fl2m=0.
      fl3m=0.
      fl4m=0.
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
      sk1m=sk1m+(v1m/sqrt(vrms(1)))**3
      sk2m=sk2m+(v2m/sqrt(vrms(2)))**3
      sk3m=sk3m+(v3m/sqrt(vrms(3)))**3
      sk4m=sk4m+(ppm/sqrt(prms))**3
      fl1m=fl1m+(v1m/sqrt(vrms(1)))**4
      fl2m=fl2m+(v2m/sqrt(vrms(2)))**4
      fl3m=fl3m+(v3m/sqrt(vrms(3)))**4
      fl4m=fl4m+(ppm/sqrt(prms))**4
  514 continue
      ske(1)=sk1m*vl123
      ske(2)=sk2m*vl123
      ske(3)=sk3m*vl123
      ske(4)=sk4m*vl123
      fla(1)=fl1m*vl123   
      fla(2)=fl2m*vl123   
      fla(3)=fl3m*vl123   
      fla(4)=fl4m*vl123  
  513 continue
      enav=enej
      do 428 l=1,3
      qrms(l)=vrms(l)
      qm(l)=vm(l)
  428 continue
      qrms(4)=vrms(4)
c
c     rms of velocity derivative 
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
c
c     skewness and flatness of velocity derivative 
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
      return
      end
c
c  ****************************** subrout inirea **********************
c
c   initial conditions when the calculation start from a previous
c   calculation
c   only two velocity components are necessary together 
c   with the third component in a single plane.
c   From this the third component is evaluated in the whole
c   box by the divergence.
c   calculation.
c
      subroutine inirea(ntii,time,q)
c
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      dimension q(ndv,m1,m2,m3)
      common/velmax/vmax(ndv),vmaxo(ndv)
      common/tstep/dt,beta,ren
      common/wallst/cflw,cfuw
      common/veltot/vit(ndv)
      common/omegas/omr,omi
      common/timref/tref
      common/inener/ene0
      common/eneav/enav,enavo
      common/d123/alx1,alx2,alx3
c
      nfil=23
      rewind(nfil)
      read(nfil) n1,n2,n3
      read(nfil) time,ren,time,dum
      read(nfil) (((q(1,i,j,k),i=1,n1),j=1,n2),k=1,n3),
     1            (((q(2,i,j,k),i=1,n1),j=1,n2),k=1,n3),
     1            ((q(3,i,j,1),i=1,n1),j=1,n2)
      close(nfil)
c
c   computes q3 from continuity
c
      do 11 kc=1,n3m-1
      kp=kpv(kc)
      do 11 jc=1,n2m
      jp=jpv(jc)
      do 11 ic=1,n1m
      ip=ipv(ic)
      dqcap=(q(1,ip,jc,kc)-q(1,ic,jc,kc))*dx1
     1     +(q(2,ic,jp,kc)-q(2,ic,jc,kc))*dx2
      q(3,ic,jc,kp)=q(3,ic,jc,kc)-dqcap/dx3
   11 continue
      call divgck(q,qmax)
      call vmaxv(q)
      write(6,*)'****************  read from field'
      write(6,*)time,(vmax(l),l=1,3),qmax,dt
      return
      end
c
c   *****************  subrout outh  ********************
c
c  writes time history of global quantities
c
      subroutine outh(ntime,time,cflm,enen,q,pr,vor)
c
      include 'param.f'
      parameter (m3m=m3-1)
      dimension q(ndv,m1,m2,m3),pr(m1,m2,m3)
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      common/corpri/yp1(m1),yp2(m2),yp3(m3)
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
      common/skder/skd1,skd2,skd3,fld1,fld2,fld3
      common/kolqua/diss,eta
      common/prenst/enstro
      common/ledat/cvisc
      dimension vor(m1,m2),vorz(m1,m2)
      common/ispec/imic
      common/names/filcnw,filcnr,filth,filou,filuu,filsf
      character*17 filcnw,filcnr,filth,filou,filuu,filsf
      pi=2.*asin(1.)
      vl123=1./float(n1m*n2m*n3m)
c                
      vo1max=0.
      vo2max=0.
      vo3max=0.
      vo1min=0.
      vo2min=0.
      vo3min=0.
      vo1rms=0.
      vo2rms=0.
      vo3rms=0.
      enstro=0.
      do 6 kc=1,n3m
      km=kmv(kc)
      do 6 jc=1,n2m
      jm=jmv(jc)
      do 6 ic=1,n1m
      im=imv(ic) 
c
c      OM_x COMPONENT
c
      omx=+(q(3,ic,jc,kc)-q(3,ic,jm,kc))*dx2
     1    -(q(2,ic,jc,kc)-q(2,ic,jc,km))*dx3
      vo1max=max(omx,vo1max)
      vo1min=min(omx,vo1min)
c
c      OM_y COMPONENT
c
      omy=+(q(1,ic,jc,kc)-q(1,ic,jc,km))*dx3
     1    -(q(3,ic,jc,kc)-q(3,im,jc,kc))*dx1
      vo2max=max(omy,vo2max)
      vo2min=min(omy,vo2min)
c
c      OM_z COMPONENT
c
      omz=-(q(1,ic,jc,kc)-q(1,ic,jm,kc))*dx2
     1    +(q(2,ic,jc,kc)-q(2,im,jc,kc))*dx1
      vo3max=max(omz,vo3max)
      vo3min=min(omz,vo3min)
      enstro=enstro+.5*(omx*omx+omy*omy+omz*omz)*vl123
      vo1rms=vo1rms+omx*omx*vl123
      vo2rms=vo2rms+omy*omy*vl123
      vo3rms=vo3rms+omz*omz*vl123
    6 continue
      if(n3m.gt.1) then
c
c   write on file thist
c
c     print*,'in outh'
      call enerca(q,pr,enen,time)
      call pdiss(q,enen,diss,flamb,eta)
c
c               compute the  enstrophy
c
c   these files are opened when a particular compilation
c   is done on IBM RISC workstations to avoi to have
c   file empty is there is some problem on the machine
c   during the run
c
c     open(34,file='enstr.out')
c     open(32,file=filth)
c     open(42,file=filou)
c     open(52,file=filuu)
c     open(62,file=filsf)
c     open(82,file='vorms.out')
c     open(92,file='tarms.out')
c     open(72,file='momder.out')
      enstrp=enstro*2.
      write(34,*)time,enen,enstrp,diss,
     1         sqrt(10./3.)*enen/cvisc/sqrt(enstro)
      nfilt=32
      write(68,158)ntime,time,enen,vmax(1),vmax(2),vmax(3)
     1           ,(vrms(1)+vrms(2)+vrms(3))/3.
     1           ,(ske(1)+ske(2)+ske(3))/3.
     1           ,(fla(1)+fla(2)+fla(3))/3.
     1           ,(skd1+skd2+skd3)/3.,(fld1+fld2+fld3)/3.
  158 format(3x,i4,10(1x,e11.5))
      write(42,258)time,vrms(1),vrms(2),vrms(3),prms
      write(92,258)time,vrms(4),vrms(5),vrms(6)
      rlamb=enen*sqrt(20./(3.*cvisc*diss))
      write(32,258)time,enen,diss,flamb,eta,rlamb
      write(52,258)time,ske(1),ske(2),ske(3),ske(4)
      write(62,258)time,fla(1),fla(2),fla(3),fla(4)
      write(72,258)time,skd1,skd2,skd3,fld1,fld2,fld3
      write(82,258)time,vo1rms,vo2rms,vo3rms
  258 format(7(e11.5,2x))
c     close(34)
c     close(32)
c     close(42)
c     close(52)
c     close(62)
c     close(72)
c     close(82)
c     close(92)
                      endif
      xmax=0.
      ymax=0.
      vomax=-1000.
      xmin=0.
      ymin=0.
      vomin=+1000.
      v1max=0.
      v2max=0.
      do i=1,n1m
            do j=1,n2m
      vorzz=0.
                 do k=1,n3m
      v1max=max(q(1,i,j,k),v1max)
      v2max=max(q(2,i,j,k),v2max)
      vor(i,j)=(q(2,i,j,k)-q(2,imv(i),j,k))*dx1
     1        -(q(1,i,j,k)-q(1,i,jmv(j),k))*dx2
      vorzz=vorzz+vor(i,j)
      if(vor(i,j).gt.vomax) then
      xmax=yp1(i)
      ymax=yp2(j)
      vomax=vor(i,j)
                            endif
      if(vor(i,j).lt.vomin) then
      xmin=yp1(i)
      ymin=yp2(j)
      vomin=vor(i,j)
                            endif
                  enddo
      vorz(i,j)=vorzz/n3m
            enddo
      enddo
c     open(32,file=filth)
      write(6,158)ntime,time,vomax,xmax,ymax,vomin,xmin,ymin
     1        ,v1max,v2max
      write(33,258)time,vomax,xmax,ymax,vomin,xmin,ymin
c     close(32)
      n2mh=n2m/2+1
      avg=1./(dx1*dx2)
      gampla=0.
      gamnea=0.
      gamplb=0.
      gamneb=0.
      do i=1,n1m
            do j=1,n2mh
      if(vorz(i,j).ge.0.)then
      gampla=gampla+vorz(i,j)*avg
                         endif
      if(vorz(i,j).le.0.)then
      gamnea=gamnea+vorz(i,j)*avg
                         endif
             enddo
      enddo
      do i=1,n1m
            do j=n2mh,n2m
      if(vorz(i,j).ge.0.)then
      gamplb=gamplb+vorz(i,j)*avg
                         endif
      if(vorz(i,j).le.0.)then
      gamneb=gamneb+vorz(i,j)*avg
                         endif
             enddo
      enddo
c     open(53,file=filuu)
      write(53,258)time,gampla,gamnea,gamplb,gamneb
c     close(53)
c     open(83,file='vomax.out')
      write(83,258)time,vo1max,vo1min,vo2max,vo2min,vo3max,vo3min
c     close(83)
      return
      end
c
c  ****************************** subrout outpf  **********************
c
      subroutine outpf(time,enen,nav,q,qcap,pr)
c
c   in here some other quantities are evaluated and printed
c   to see how the simulation proceed
c   e.g the energy spectra
c
      include 'param.f'
      parameter (m1m=m1-1,m2m=m2-1)
      common/dim/n1,n1m,n2,n2m,n3,n3m
      dimension q(ndv,m1,m2,m3),qcap(m1,m2,m3),pr(m1,m2,m3)
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
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
      common/corpri/yp1(m1),yp2(m2),yp3(m3)
      common/speene/e(ndv,0:m3)
      common/numond/kx(m1),ky(m2),kz(m3),kkmax
      common/ledat/cvisc
      common/prenst/enstro
      common/kolqua/diss,eta
      common/ispec/imic
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      dimension vorzm(m1,m2),velzm(m1,m2),premm(m1,m2)
      dimension vorxm(m1,m2),vorym(m1,m2)
      dimension velxm(m1,m2),velym(m1,m2)
      dimension vorxp(m1,m2),voryp(m1,m2),vorzp(m1,m2)
c      
      character*17 filcnw,filcnr,filth,filou,filuu,filsf
      character*17 filcos
      character*4 pntim
      character*6 slp
      character*1 is1s,is3s
      character*1 is1n,is3n
      itim=time+0.3
      do 615 l=1,3
      skp(l)=skeo(l)/nav
      flp(l)=flao(l)/nav
      vmp(l)=vmo(l)/nav
      qrms(l)=sqrt(qrms(l))
  615 vprms(l)=sqrt(vrmso(l)/nav)
      skp(4)=skeo(4)/nav
      flp(4)=flao(4)/nav
      vprms(4)=vrmso(4)/nav
      ppmp=prmso
  611 continue
  612 format(1x,i3,2x,e14.7,2x,9(1x,e11.4)) 
  613 format(e14.7,1x,9(1x,e11.4)) 
c     close(42)
c     close(62)
      if(imic.gt.0) then
      rlapr=sqrt(10./3.)*enen/cvisc/sqrt(enstro)
      if(time.lt.2.) then
      rlapr=sqrt(10./3.)*enen/cvisc/sqrt(enstro)
      itim=rlapr   
                     endif
   77 format(i4.4)
      enavp=enavo/nav
      if(time.lt.2.) then
      filcos='spetrl.'//pntim
      open(58,file=filcos,status='unknown')
      rewind(58)
      filcos='specrl.'//pntim
      open(57,file=filcos,status='unknown')
      rewind(57)
                     else
      filcos='spet.'//pntim
      open(58,file=filcos,status='unknown')
      rewind(58)
      filcos='spec.'//pntim
      open(57,file=filcos,status='unknown')
      rewind(57)
                     endif
      call spectre(q,qcap)
       ee1to=0.
       ee2to=0.
       ee3to=0.
       ensto=0.
       n1mhp=n1m/2+1
       n2mhp=n2m/2+1
       n3mhp=n3m/2+1
       ikma=max(n1mhp,n3mhp)
       kmax=max(ikma,n2mhp)
       do k=1,kmax
       write (58,*)k,e(1,k)+e(2,k)+e(3,k)
       if(k.le.kmax) write (57,*)k,e(1,k),e(2,k),e(3,k)
       ee1to=ee1to+e(1,k)
       ee2to=ee2to+e(2,k)
       ee3to=ee3to+e(3,k)
       eeto=ee1to+ee2to+ee3to
       ensto=ensto+k**2*(e(1,k)+e(2,k)+e(3,k))
       end do
      dissk=2.*ensto*cvisc
      close(58)
      close(57)
      if(time.lt.2.) then
c
c   the spectra in kolmogorov units are calculated
c
      filcos='spekorl.'//pntim
                     else
      filcos='speko.'//pntim
                     endif
      open(59,file=filcos,status='unknown')
      rewind(59)
      speska=(1./(cvisc**5.*dissk))**(1./4.)
      etak=(cvisc**3./dissk)**(1./4.)
       do k=1,kmax
       akol=k*etak
       se1=e(1,k)*speska
       se2=e(2,k)*speska
       se3=e(3,k)*speska
       if(k.le.kmax) write (59,*)akol,se1,se2,se3          
       end do
      close(59)
      eetp=eeto*0.5
      rlamc=eetp*sqrt(10./(3.*cvisc*dissk))
      write(6,133)time,ee1to,ee2to,ee3to,eeto,enen
      write(66,*)time,eeto,ensto,rlamc,etak
  133 format('energy',e12.4,3x,3e12.4,3x,'etot sp,phy',2e12.5)
                     endif
c
c   here averaged quantities are calculated in the x1-x2 plane
c   this calculation makes sense for vortex dynamics simulations
c
      write(pntim,77) itim
c     write(6,*)' in outh print 2d'
      xmax=0.
      ymax=0.
      vomax=-100.
      xmin=0.
      ymin=0.
      vomin=+100.
      do ic=1,n1m
      im=imv(ic)                                                    
            do jc=1,n2m
      jm=jmv(jc)                                                    
      vorzp(ic,jc)=0.
      vorxp(ic,jc)=0.
      voryp(ic,jc)=0.
      premm(ic,jc)=0.
      velzm(ic,jc)=0.
      velxm(ic,jc)=0.
      velym(ic,jc)=0.
                  do kc=1,n3m
      km=kmv(kc)                                                    
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
      vmx=(q(1,ic,jc,kc)+q(1,im,jc,kc))*0.5
      vmy=(q(2,ic,jc,kc)+q(2,ic,jm,kc))*0.5
      vmz=(q(3,ic,jc,kc)+q(3,ic,jc,km))*0.5
      premm(ic,jc)=premm(ic,jc)+pr(ic,jc,kc)
      velxm(ic,jc)=velxm(ic,jc)+vmx
      velym(ic,jc)=velym(ic,jc)+vmy
      velzm(ic,jc)=velzm(ic,jc)+vmz
      vorxp(ic,jc)=vorxp(ic,jc)+omx
      voryp(ic,jc)=voryp(ic,jc)+omy
      vorzp(ic,jc)=vorzp(ic,jc)+omz
                  enddo
      velxm(ic,jc)=velxm(ic,jc)/n3m
      velym(ic,jc)=velym(ic,jc)/n3m
      velzm(ic,jc)=velzm(ic,jc)/n3m
      premm(ic,jc)=premm(ic,jc)/n3m
      vorxp(ic,jc)=vorxp(ic,jc)/n3m
      voryp(ic,jc)=voryp(ic,jc)/n3m
      vorzp(ic,jc)=vorzp(ic,jc)/n3m
            enddo
      enddo
c     print*,'evaluated averages outh ',nav,imic
      do ic=1,n1m
            do jc=1,n2m
      vorzm(ic,jc)=(vorzp(ic,jc)+vorzp(ipv(ic),jc)
     1             +vorzp(ic,jpv(jc))+vorzp(ipv(ic),jpv(jc)))*0.25
      vorxm(ic,jc)=(vorxp(ic,jc)+vorxp(ic,jpv(jc)))*0.5
      vorym(ic,jc)=(voryp(ic,jc)+voryp(ipv(ic),jc))*0.5
      if(vorzp(ic,jc).gt.vomax) then
      xmax=yp1(ic)
      ymax=yp2(jc)
      vomax=vorzm(ic,jc)
                            endif
      if(vorzp(ic,jc).lt.vomin) then
      xmin=yp1(ic)
      ymin=yp2(jc)
      vomin=vorzm(ic,jc)
                            endif
            enddo
      enddo
      write(6,158)ntime,time,vomax,xmax,ymax,vomin,xmin,ymin
  158 format(3x,i4,10(1x,e11.5))
c
c   here the mean velocity and vorticity are
c   written into a 2D file to perform flow
c   visualizations by TURB3D
c
      filcos='vevoavgz'//pntim//'.dat'
      open(57,file=filcos,form='formatted',status='unknown')
      write(6,*)' wrote ',filcos
      write(57,*) n1m,n2m,1
      write(57,*) time,time,time,time
      write(57,*) ((velxm(i,j),i=1,n1m),j=1,n2m)
     1         ,((velym(i,j),i=1,n1m),j=1,n2m)
     1         ,((velzm(i,j),i=1,n1m),j=1,n2m)
     1         ,((premm(i,j),i=1,n1m),j=1,n2m)
     1         ,((vorxm(i,j),i=1,n1m),j=1,n2m)
     1         ,((vorym(i,j),i=1,n1m),j=1,n2m)
     1         ,((vorzm(i,j),i=1,n1m),j=1,n2m)
      close(57)
      return
      end
c************************************************************
c****************  pdiss   **********************************
c************************************************************
c
c    calculate the dissipation. average over strain at 4 points 
c  diss  proportional to SijSij
c
      subroutine pdiss(q,enej,diss,flamb,eta)
      include 'param.f'
      dimension q(ndv,m1,m2,m3)
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      common/d123/alx1,alx2,alx3
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/ledat/cvisc
      common/d2/nstop,nprint,ntst,npin,npstf
      dimension visout(m2)
      character*20 titfil
      character*3 icount
c      
c     print*,'in pdiss'
      vl123=1./float(n1m*n2m*n3m)
      diss=0.
      diss1=0.
      do 4 kc=1,n3m 
      kp=kpv(kc)                                                    
      km=kmv(kc)                                                    
      do 4 jc=1,n2m                                                     
      jp=jpv(jc)
      jm=jmv(jc)
      do 4 ic=1,n1m                                                     
      ip=ipv(ic)                                                        
      im=imv(ic)                                                        
c
c   stl  =Sll  at i+1/2,j+1/2,k+1/2
c   st2  =S22  at i+1/2,j+1/2,k+1/2
c   st3  =S33  at i+1/2,j+1/2,k+1/2
c
      st1=(q(1,ip,jc,kc)-q(1,ic,jc,kc))*dx1
      st2=(q(2,ic,jp,kc)-q(2,ic,jc,kc))*dx2
      st3=(q(3,ic,jc,kp)-q(3,ic,jc,kc))*dx3
c
c   st4 = S12    at  i,j,k+1/2
c   st6 = S32    at  i+1/2,j,k
c   st5 = S13    at  i,j+1/2,k
c
c
      stccc4=(q(1,ic,jc,kc)-q(1,ic,jm,kc))*dx2
     1               +(q(2,ic,jc,kc)-q(2,im,jc,kc))*dx1
      stccc5=(q(1,ic,jc,kc)-q(1,ic,jc,km))*dx3
     1              +(q(3,ic,jc,kc)-q(3,im,jc,kc))*dx1
      stccc6=(q(3,ic,jc,kc)-q(3,ic,jm,kc))*dx2
     1               +(q(2,ic,jc,kc)-q(2,ic,jc,km))*dx3
c
      stcpc4=(q(1,ic,jp,kc)-q(1,ic,jc,kc))*dx2
     1               +(q(2,ic,jp,kc)-q(2,im,jp,kc))*dx1
      stccp5=(q(1,ic,jc,kp)-q(1,ic,jc,kc))*dx3
     1              +(q(3,ic,jc,kp)-q(3,im,jc,kp))*dx1
      stcpc6=(q(3,ic,jp,kc)-q(3,ic,jc,kc))*dx2
     1               +(q(2,ic,jp,kc)-q(2,ic,jp,km))*dx3
c
      stpcc4=(q(1,ip,jc,kc)-q(1,ip,jm,kc))*dx2
     1               +(q(2,ip,jc,kc)-q(2,ic,jc,kc))*dx1
      stpcc5=(q(1,ip,jc,kc)-q(1,ip,jc,km))*dx3
     1              +(q(3,ip,jc,kc)-q(3,ic,jc,kc))*dx1
      stccp6=(q(3,ic,jc,kp)-q(3,ic,jm,kp))*dx2
     1               +(q(2,ic,jc,kp)-q(2,ic,jc,kc))*dx3
c
      stppc4=(q(1,ip,jp,kc)-q(1,ip,jc,kc))*dx2
     1               +(q(2,ip,jp,kc)-q(2,ic,jp,kc))*dx1
      stpcp5=(q(1,ip,jc,kp)-q(1,ip,jc,kc))*dx3
     1              +(q(3,ip,jc,kp)-q(3,ic,jc,kp))*dx1
      stcpp6=(q(3,ic,jp,kp)-q(3,ic,jc,kp))*dx2
     1               +(q(2,ic,jp,kp)-q(2,ic,jp,kc))*dx3
c
c   st4 = S12    at  center
c   st5 = S13    at  center
c   st6 = S32    at  center
c
      st4=.5*0.25*(stccc4+stcpc4+stpcc4+stppc4)
      st5=.5*0.25*(stccc5+stccp5+stpcc5+stpcp5)
      st6=.5*0.25*(stccc6+stcpc6+stccp6+stcpp6)
        app1=2.*(st1*st1+st2*st2+ st3*st3+
     1   2.*(st4*st4+st5*st5+st6*st6))*vl123
        diss=diss+app1*cvisc
    4 continue                                                          
c      
c  Kolmogorov scale
c  and Taylor microscale
c
      eta=(cvisc**3./diss)**(1./4.)
      flamb=sqrt(20*cvisc*enej/diss)
      return
      end
