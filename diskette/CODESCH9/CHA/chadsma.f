c
c  random nuber generator appearing as ranset(1) and ranf( ) in inqpr
c
      program main
      common/tstep/dt,beta,ren
      common/d2/nstop,nprint,ntst,npin,npstf
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/wrre/nwrit,nread
      common/strpar/str2
      common/tscoe/ga(3),ro(3),nsst
      common/d13/alx1,alx3
      common/vperin/vper
      common/averou/iav
      common/islwal/islv1s,islv1n,islv3s,islv3n
      common/jrrr/jri,jrf,djr,irejr,iruuca
      common/newdat/itimsq,timeav
      common/movwal/tosc,uosc
      common/slotfl/flowq2,tau2
      common/slotii/tim0sl
      common/slotpa/y1gsl,y1ssl,y3gsl,y3ssl
      common/slotdi/y1disl,y3disl
      common/oldso/ifield
      common/cflco/icfl,cflc,tpin,tprin,tfin                      
      common/iprf/iprfi
      open(15,file='chapn.d')
      read(15,*) n1,n2,n3,nsst
      read(15,*) nwrit,nread,iav,iprfi
      read(15,*) alx3d,alx1d,str2
      read(15,*) ren,vper
      read(15,*) dt,ntst,nprint,npin,npstf
      read(15,*) nstop
      read(15,*) icfl,cflc,tpin,tprin,tfin                      
      read(15,*) jri,jrf,djr,irejr,iruuca
      read(15,*) timeav,itimsq
      read(15,*) islv1s,islv1n,islv3s,islv3n
      read(15,*) tosc,uosc
      read(15,*) flowq2,tau2,tim0sl
      read(15,*) y1gsd,y1ssd,y3gsd,y3ssd
      read(15,*) y1disd,y3disd 
      read(15,*) ifield
c     ifield=0 legge solo la soluzione
c     ifield=1 legge anche le velocita' della parete inferiore
      pi=2.*asin(1.)
      n1m=n1-1
      n2m=n2-1
      n3m=n3-1
      alx1=alx1d*pi
      alx3=alx3d*pi
      y1gsl=y1gsd*pi
      y1ssl=y1ssd*pi
      y3gsl=y3gsd*pi
      y3ssl=y3ssd*pi
      y1disl=y1disd*pi
      y3disl=y3disd*pi
c     if nsst=3 Runge-Kutta scheme
      if(nsst.gt.1) then
      ga(1)=8./15.
      ga(2)=5./12.
      ga(3)=3./4.
      ro(1)=0.
      ro(2)=-17./60.
      ro(3)=-5./12.
      write(6,100) (ga(n),n=1,nsst),(ro(n),n=1,nsst)
  100 format(10x,'third order runge gutta=',4x,
     1       'gam=',3f8.3,4x,'ro',3f8.3)
      else
c     if nsst=1 Adams-Bashfort scheme
      ga(1)=1.5
      ga(2)=0.
      ga(3)=0.
      ro(1)=-.5
      ro(2)=0.
      ro(3)=0.
      write(6,110) ga(1),ro(1)
  110 format(10x,'adams bashfort=',4x,
     1       'gam=',f8.3,4x,'ro',f8.3)
      endif
c
      write(6,112)alx1d,alx3d
  112 format(10x,'chann. dimens ly=2',3x,'lx=',f4.2,'*pi'
     1        ,3x,'lz=',f4.2,'*pi')
      write(6,201)vper
      write(6,200)
  200 format(10x,'3d channel periodic in x3 and x1 ')
  201 format(3x,'programma con init. random perturbation'
     1   ,5x,'vper',e10.3)
      write(6,202) tosc,uosc,tim0sl,tau2
  202 format(10x,'parete inferiore, tosc=',e11.4,2x,
     1   'uosc=',e11.4,2x,'tim0sl=',e11.4,2x,'tau2=',e11.4)
      call openfi
      call solve
      stop
      end
c
c  ************************************************************
c  ************subrout solve **********************************
c  ************************************************************
c
      subroutine solve
c
c     code for computation of three-dimensional incompressible flows
c     in cartesian coordinates.
c
c     this code solves flow fields periodic in x1 and x3 and non-slip
c     conditions in x2.
c     non uniform coordiante only in x2 direction.
c     the geometry y(j) is given in the subroutine coordi.
c     the equations for yhe cartesian component u(i) i=1,3
c     are discretized by a finite difference scheme.
c ??  cij=dyi/dxj   are metric quantities in the subroutine metric. 
c     all the spatial derivatives are discretized by centered schemes.
c     including the non linear terms.
c     a factored scheme brings to the solut. of three tridiag. matrices
c     for each vel. component  subr invtr
c     in time a fractional step is used in the version of nagi mansour
c     introducing the press. in the first step.
c     the non-linear terms and the cross der. are discr. by explicit
c     adams-bashfort or third order range-kutta.
c     the scalar ph which gives the div. free vel. is solved by a direct
c     method. for the x3 direction a modified wave number is introduced
c     and a fft method is used to have the ph in the physical space.
c     no boundary conditions are necess. for the poisson equation.
c                a*ph=div(q)
c     in the subr. phini the coeff of the matrix a are calculated
c
c     all the second derivatives are discretized by a second order conse
c     rvative scheme.
c
      include 'param.f'
      parameter (m1m=m1-1)
c
      parameter (m2k=13)
      dimension xrko(3,m1m,m2k)
      common/wrre/nwrit,nread
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      dimension y(m2)
      dimension ru(ndv,m1,m2,m3)
      dimension q(ndv,m1,m2,m3),qcap(m1,m2,m3)
     1         ,dph(m1,m2,m3),pr(m1,m2,m3)
      common/velmax/vmax(ndv),vmaxo(ndv)
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/d13/alx1,alx3
      common/d2/nstop,nprint,ntst,npin,npstf
      common/tstep/dt,beta,ren
      common/inener/ene0
      common/averou/iav
      common/sc/sc,sc1
      common/names/filcnw,filcnr,filth,filou,filuu,filsf
      common/jrrr/jri,jrf,djr,irejr,iruuca
      common/newdat/itimsq,timeav
      common/cflco/icfl,cflc,tpin,tprin,tfin                      
      character*17 filcnw,filcnr,filth,filou,filuu,filsf
      character*4 pntim
      istop=0
      nat=0
c
c
      npfile=npin*npstf
      pi=2.*asin(1.)
c
c     step  and mesh sizes calculations
c
c
      call meshes
      call indic
      call coordi(y)
      call metric(y)
c
      write(6,754)n1,n2,n3
  754 format(10x,'centered velocities',2x
     1      ,5x,'n1=',i3,2x,'n2=',i3,2x,'n3=',i3)
      write(6,755) ren,dx1,dx2,dx3,ntst
  755 format(3x,'ren=',e10.3,3x,'dx1=',e10.3,3x,'dx2=',e10.3,3x
     1      ,3x,'dx3=',e10.3,3x,
     1      'ntst=',i5,3x)
      write(6,*) '  '
      write(32,*) '  '
      if(icfl.eq.1) then
      ntst=100000
      nstop=100000
      dtl=dt
      write(6,*)' calculation at dt variable by fixing cfl=',cflc
      write(32,*)' calculation at dt variable by fixing cfl=',cflc
                     else
      write(6,765) dt
      write(32,765) dt
  765 format(' time step dt= ',e10.3)
                     endif
      write(6,*) '  '
      write(32,*) '  '
      nti=0
      time=0.
      ntii=0
      n1mh=n1m/2+1
c   correlation
      jj=0
      do 20 j=jri,jrf,djr
      jj=jj+1
      write(26,*)jj
      do 20 i=1,n1mh
      do 20 l=1,3
      xrko(l,i,jj)=0.
   20 continue
c
c     initial conditions
c
      if(nread.eq.0) then
      do 22 l=1,ndv
      vmax(l)=0.1e-10
   22 continue
      ntii=0
      ntime=0
      time=0.
      nmedr=0
      cflm=0.
      nav=0
c
      call initia(q,y,pr)
c     call wstre(q)
c     call enerca(q,ene0,pr)
c     call oldqua
c     call outh(ntime,time,cflm,ene0)
      if(ene0.lt..1e-10) ene0=.1e-10
      else
c

      call slotin
      call inirea(ntii,time,q,ru,pr)
      if(icfl.eq.1) then
c
c   courant number calculation for concentric cilinders
c
      call cfl(q,cflmm)
      dt=cflc/cflmm
      cflm=cflc
      if(dt.gt.dtl) then
      dt=dtl
      beta=dt/ren*0.5
      cflm=cflmm*dt
                    endif
      write(6,*)' with inirea=1 dt=',dt
                   endif
      if(iav.eq.1) then
      call tavrea(nav)
                   endif
      if(irejr.eq.1.and.iruuca.eq.1) then
      call rearuu(xrko,nmedr)
                   endif
      ntime=ntii
      itim=time+0.3
      write(pntim,77)itim
   77 format(i4.4)
      filth='cpth.'//pntim
      close(32)
      open(32,file=filth)
c
      call wstre(q)
c
      call enerca(q,ene0,pr)
c
c     modifica 15/10/93 per retstart da campi interpolati
c     call divgck(q,qmax)
c
c     write(6,7698) qmax,ntime
c     if(qmax.gt..1e-03) go to 169
c      call outh(ntime,time,cflm,ene0)
      endif
c
c  ********* start of time dependent calculation ***
c
      do i=1,n1m
      do j=1,n2m
      do k=1,n3m
      qcap(i,j,k)=1.
      dph(i,j,k)=1.
      enddo
      enddo
      enddo
      call phini(qcap,dph)
c
      ntstf=ntii+ntst
      ntii=ntii+1
      write(6,711)nprint,ntii,ntstf,nstop,nti,dt
  711 format(3x,'check in cond',5i8,2x,e11.3)
  181 continue
c
c    the dipendent calculation start
c
c
      beta=dt/ren*0.5
c
      call coeinv
      ncount=0
      do 350 ntime=ntii,ntstf
c     if(time.gt.250.) iav=1
      if(time.gt.timeav) iav=1
      if(ntime.gt.nstop) go to 167
c     
c    
      if(icfl.eq.1) then
c
c   courant number calculation for running at cfl=constant
c
      call cfl(q,cflum)
      cflmm=cflum
      cflmp=cflum*dt
      dt=cflc/cflmm
      cflm=cflc
      write(56,*)ntime,dtl,dt,cflmp
      if(dt.gt.dtl) then
      dt=dtl
      beta=dt/ren*0.5
                    endif
      if(dt.lt..1e-02)      then
      write(6,368)dt
  368 format(3x,'end calcul for dt=',e12.4)
        call contwr(ntime,time,q,ru,pr,enen)
        if(iav.eq.1) then
         call tavwri(nav)
      call outpf(time,enen,nav,q)
                     endif
      stop
                            endif
                    endif
      if(ren.gt.1.e07) then
      call tsinv(q,pr,ru,qcap,dph,time) 
                       else
      call tschem(q,pr,ru,qcap,dph,time) 
                       endif
c
      time=time+dt
      if(time.gt.timeav.and.itimsq.eq.1) then
      call timseq(time,q,pr)
                         endif
c
c     print some quantity cf,max(divg),total energy,amx(veloc),
c     and theirs gradient in time.
c                                                                       
      if(icfl.eq.1) then
      if(amod(time,tpin).lt.dt) go to 306
      go to 305
                     else
      if(mod(ntime,npin).eq.0) then
      call cfl(q,cflum)
      cflmm=cflum
      cflm=cflmm*dt
           go to 306
                               endif
                     endif
  306 continue
c
      call vmaxv(q)
c
c     the calculation stop if the velocities are diverging for numer
c     stab conditions (courant number restrictions)
c
      if(vmax(1).gt.1000.and.vmax(2).gt.1000) go to 166
c
c    lower and upper skin friction
c
      call wstre(q)
c
c     total energy
c
      call enerca(q,enen,pr)
c
      if(iav.eq.1) then
       nav=nav+1
       call taver
      endif
      if(iruuca.eq.1) then
      nmedr=nmedr+1
      call ruucal(qcap,dph,xrko,nmedr,q,nav)
                   endif


c
      call divgck(q,qmax)
c
      ntick=ntime-ntii 
      if(qmax.gt..2e-03.and.ntick.ge.21) go to 169
      write(66,7698) qmax,ntime,cflm,dt
 7698 format(3x,'max div=',g10.4,2x,'iteration = ',i5
     1       ,3x,'cflm=',e12.5,3x,'dt=',e12.5)
c
      call outh(ntime,time,cflm,enen)
c
  305 continue
c  
c    write the flow field
c
c*******************************************************
      if(icfl.eq.1.and.time.ge.tpin) then
      if(amod(time,tprin).lt.dt) go to 301
            go to 300
                     else
      if(mod(ntime,nprint).eq.0) go to 301
            go to 300
                     endif
  301 continue
            if(nwrit.eq.1) then                                             
        call contwr(ntime,time,q,ru,pr,enen)
        if(iav.eq.1) then
         call tavwri(nav)
      call outpf(time,enen,nav,q)
                     endif
                           endif                                                           
  300 continue
       if(time.gt.tfin.and.icfl.eq.1) then
      write(6,367)time
  367 format(3x,'end calcul for t=',e12.4)
        call contwr(ntime,time,q,ru,pr,enen)
      stop
                                        endif
  350 continue                                                          
c
c
      go to 167
  166 continue
      write(6,168)
  168 format(10x,'calculation diverged in time for the vel field')
      write(6,*)'  vmax(1,2,3 ',(vmax(l),l=1,3)
      write(6,7698) qmax,ntime,cflm
      if(iav.eq.0) nav=1
      call outpf(time,enen,nav,q)
      go to 167
  169 continue
      write(6,178)qmax
  178 format(10x,'calculation diverged for  dmax='
     1 ,e12.5)
  167 continue
      return
      end
