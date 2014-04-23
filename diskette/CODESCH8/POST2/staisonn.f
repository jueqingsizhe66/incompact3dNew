c
c  ****************************** subrout coordi  **********************
c    this subroutine assign the physical coordinates as function
c    of the computational ones x1,x2
c
      subroutine coordi
c
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/d123/alx1,alx2,alx3
      common/corpri/yp1(m1),yp2(m2),yp3(m3)
      dimension y1s(m1),y3s(m3),y2s(m2)
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      common/ispec/imic
      character*60 namfile
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
      do  i=1,n1m
          y1s(i) = 0.5*(yp1(i) + yp1(i+1) )
      enddo
      do  k=1,n3m
          y3s(k) = 0.5*(yp3(k) + yp3(k+1) )
      enddo
      do j=1,n2m
          y2s(j)=(yp2(j)+yp2(j+1))*0.5
      enddo

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
c  these are necessary when the equation are solved near the
c  walls and for the periodic conditions.
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
  311 continue
      do 310 k=1,n3m
      do 310 j=1,n2m
      do 310 i=1,n1m
      vca1=q(1,i,j,k)
      vmax(1)=max(vmax(1),abs(vca1))
      vca2=q(2,i,j,k)
      vmax(2)=max(vmax(2),abs(vca2))
      vca3=q(3,i,j,k)
      vmax(3)=max(vmax(3),abs(vca3))
  310 continue
      return
      end
c
c  ****************************** subrout inirea **********************
c
c   reading initial  conditions from file
c   this file is different from the restarting file.
c   It must be produced by running only one time step
c   of the ISO code with a very small Delta t.
c   The user should insert in the code ISO the routine
c   that print the field in the following format.
c   ************   IMPORTANT ******************************
c
      subroutine inirea(ntii,time,q,pr)
c
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      dimension pr(m1,m2,m3)
      dimension q(ndv,m1,m2,m3)
      common/velmax/vmax(ndv),vmaxo(ndv)
      common/tstep/dt,beta,ren
      common/wallst/cflw,cfuw
      common/veltot/vit(ndv)
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
      read(nfil) (((q(1,i,j,k),i=1,n1),j=1,n2),k=1,n3),
     1            (((q(2,i,j,k),i=1,n1),j=1,n2),k=1,n3),
     1            (((q(3,i,j,k),i=1,n1),j=1,n2),k=1,n3)
     1           ,(((pr(i,j,k),i=1,n1),j=1,n2),k=1,n3)
      close(nfil)
      return
      end
c
c  ****************************** subrout vort  **********************
c
      subroutine vort(q,ru)
c
c    The vorticity components are calculated
c
      include 'param.f'
      parameter (m1m=m1-1,m2m=m2-1)
      common/dim/n1,n1m,n2,n2m,n3,n3m
      dimension ru(ndv,m1,m2,m3)
      dimension q(ndv,m1,m2,m3)
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/omean/vo(ndv),vorms(3)
c
c  THE vorticity components are staggered  
C        o1 at the  POINT I+1/2,J,K
C        o2 at the  POINT I,J+1/2,K
C        o3 at the  POINT I,J,K+1/2
c
      o1max=0.
      o2max=0.
      o3max=0.
      vo1m=0.
      vo2m=0.
      vo3m=0.
      vl123=1./float(n1m*n2m*n3m)
      do k=1,n3m
      kp=kpv(k)
      km=kmv(k)
      do i=1,n1m
      ip=ipv(i)
      im=imv(i)
      do j=1,n2m
      jp=jpv(j)
      jm=jmv(j)
c
c      OM_1 COMPONENT
c
      om1=+(q(3,i,j,k)-q(3,i,jm,k))*dx2
     1    -(q(2,i,j,k)-q(2,i,j,km))*dx3
c
c      OM_2 COMPONENT
c
      om2=+(q(1,i,j,k)-q(1,i,j,km))*dx3
     1    -(q(3,i,j,k)-q(3,im,j,k))*dx1
c
c      OM_3 COMPONENT
c
      om3=-(q(1,i,j,k)-q(1,i,jm,k))*dx2
     1    +(q(2,i,j,k)-q(2,im,j,k))*dx1
      o1max=max(abs(om1),o1max)
      o2max=max(abs(om2),o2max)
      o3max=max(abs(om3),o3max)
      ru(1,i,j,k)=om1
      ru(2,i,j,k)=om2
      ru(3,i,j,k)=om3
       vo1m=ru(1,i,j,k)+vo1m
       vo2m=ru(2,i,j,k)+vo2m
       vo3m=ru(3,i,j,k)+vo3m
      enddo
      enddo
      enddo
      vo(1)=vo1m*vl123
      vo(2)=vo2m*vl123
      vo(3)=vo3m*vl123
      o1rms=0.
      o2rms=0.
      o3rms=0.
      do k=1,n3m
      do j=1,n2m
      do i=1,n1m
       o1m=ru(1,i,j,k)-vo(1)
       o2m=ru(2,i,j,k)-vo(2)
       o3m=ru(3,i,j,k)-vo(3)
      o1rms=o1rms+o1m**2
      o2rms=o2rms+o2m**2
      o3rms=o3rms+o3m**2
      enddo
      enddo
      enddo
      vorms(1)=o1rms*vl123
      vorms(2)=o2rms*vl123
      vorms(3)=o3rms*vl123
      write(6,*) 'vort max',o1max,o2max,o3max
      write(6,*) 'vo=',vo(1),vo(2),vo(3)
      write(6,*) 'vorms=',vorms(1),vorms(2),vorms(3)
      return
      end
c
c  ****************************** subrout enerca **********************
c
c   calculation of the one-point statistics
c   e.g.  total energy = (v1**2+v2**2+v3**2)*0.5
c
      subroutine enerca(q,pr,enej,time)
c
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      dimension q(ndv,m1,m2,m3),pr(m1,m2,m3)
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      common/meshu/d1x1,d1x2,d1x3
      common/d1/re
      common/d123/alx1,alx2,alx3
      common/veltot/vit(ndv)
      common/eneav/enav,enavo
      common/vmean/vm(ndv),vrms(4)
      common/skfl/ske(4),fla(4)
      common/prrm/prm,prms
      common/skder/skd1,skd2,skd3,fld1,fld2,fld3
      common/ledat/ics0,cvisc,ifiltr,ibox
c
      pi=2.*asin(1.)
      vol=(2.*pi)**3.
      vl13=1./float(n1m*n3m)
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
      enej=0.
      u1rmm=0.
      u2rmm=0.
      u3rmm=0.
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
      ppm=pr(i,j,k)-prm
      pr(i,j,k)=ppm
      u1rmm=u1rmm+v1m**2
      u2rmm=u2rmm+v2m**2
      u3rmm=u3rmm+v3m**2
      u23mm=u23mm+v2m*v3m
      ppmm=ppmm+ppm**2
      enej=enej+(v1m**2+v2m**2+v3m**2)
  414 continue
      enej=enej*vl123
      vrms(1)=u1rmm*vl123
      vrms(2)=u2rmm*vl123
      vrms(3)=u3rmm*vl123
      vrms(4)=u23mm*vl123
      prms=ppmm*vl123
c
c   skewness and flatness
c
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
      v1m=(q(1,ip,j,k)+q(1,i,j,k))*0.5-vm(1)
      v2m=(q(2,i,jp,k)+q(2,i,j,k))*0.5-vm(2)
      v3m=(q(3,i,j,kp)+q(3,i,j,k))*0.5-vm(3)
      ppm=pr(i,j,k)
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
      write(6,*) 'vm=',vm(1),vm(2),vm(3)
      write(6,*) 'pm=',prm,'   prms=',prms
      write(6,*) 'vrms=',vrms(1),vrms(2),vrms(3),vrms(4)
      write(6,*) 'ske=',ske(1),ske(2),ske(3),ske(4)
      write(6,*) 'fla=',fla(1),fla(2),fla(3),fla(4)
c
c     rms of the velocity derivative 
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
c     write(*,*) 'rmsder=',rmsdux,rmsduy,rmsduz
c
c
c     skewness and flatness of the velocity derivative 
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
      write(*,*) 'skeder=',skd1,skd2,skd3
      write(*,*) 'flader=',fld1,fld2,fld3
      return
      end
cc
c************************************************************
c****************  pdiss   **********************************
c************************************************************
c calculate the dissipation. average over strain at 4 points 
c
      subroutine pdiss(q,diss,sca)
      include 'param.f'
      dimension q(ndv,m1,m2,m3)
      dimension sca(m1,m2,m3)
      dimension ell2(m2),ell1(m2)
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      common/d123/alx1,alx2,alx3
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/ledat/ics0,cvisc,ifiltr,ibox
      common/d2/nstop,nprint,ntst,npin,npstf
      common/cosma/csma
      dimension visout(m2)
      character*20 titfil
      character*3 icount
c      
c
c      computation of Sij
c
c     this calculates the  strain rate tensor at the center of
c     the cell. Perodic box
c                                                                       
c     
c     print*,'in pdiss'
      dima=-100.
      dimi=+100.
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
        app1=.5*(st1*st1+st2*st2+ st3*st3+
     1   2.*(st4*st4+st5*st5+st6*st6))*4.
        diss=diss+app1
      dima=max(app1,dima)
      dimi=min(app1,dimi)
      sca(ic,jc,kc)=app1
    4 continue                                                          
      diss=diss*vl123
      write(6,*)'diss ',diss,dima,dimi
      return
      end
