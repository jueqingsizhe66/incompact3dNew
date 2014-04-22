      subroutine openfile
c
c  in this routine certain files to write the results are opened
      open(15,file='turb2d.d')
      open(18,file='co.dat',form='unformatted')
      open(20,file='thdite.out')
      open(22,file='hioinv.out')
      open(23,file='specin.out')
      open(24,file='specut.out')
      open(16,file='outsc.out')
      rewind(18)
      rewind(20)
      rewind(16)
      return
      end
c *********************************************************************************
c  In this code th Navier-Stokes equations in a Cartesian grid
c  written in vorticity streamfunction formulations are solved. 
c  The grid is uniform in both directions and only periodic
c  conditions are treated.
c  The code allows the solution of 2D turbulence strating from
c  a vorticity field relative to a certain preassigned energy
c  spectrum. There is the possibility to insert te initial 
c  field in a much larger domain.
c  The spoace discretization for the non-linear terms
c  is obtained by second or fourth order Arakawa scheme,
c  explained in the book in SEct.6.2
c  The time advancement can be by the Adams Bashfort or third
c  order Runge Kutta scheme
c
      program main
      common/d1/tfin,alx1f,alx1i,alx2i,alx2f
      common/d2/nstop,nprint,ntst,npin,nread,nwrit,nprspe
      common/dim/n1,n1m,n2,n2m
      common/tscoe/ga(3),ro(3),nsst
      common/chnlc/chal,chbe,chga
      common/visct/re
      common/tstep/dt
      common/icci/icut
      common/angmd/thet0
      common/spet0/ak0,en0,gam0
      common/idue/timrea
      common/inisp/ityp
      common/iark4/iforth
      common/icflp/icflm,cflma,dtl
      common/tfinst/tprin,tpin,tprspe
      common/dimrat/factpi
      common/ened/enede
c
      call openfile
c
      read(15,*) n1,n2,nsst,factpi
      read(15,*) dt,ntst,nprint,npin,nprspe
      read(15,*) nstop,nread,nwrit,nkt
      read(15,*) icflm,cflma,dtl
      read(15,*) tprin,tpin,tfin,tprspe
      read(15,*) re,timrea
      read(15,*) chal,chbe,chga
      read(15,*) icut
      read(15,*) ak0,en0,gam0,ityp,enede
      read(15,*) iforth
c
  33  format(a3)
      pi=2.*asin(1.)
      aldim=2.*pi*factpi
      alx1f=aldim
      alx2f=aldim
      alx1i=0.
      alx2i=0.
      n1m=n1-1
      n2m=n2-1
      if(nsst.eq.1) go to 10
      ga(1)=8./15.
      ga(2)=5./12.
      ga(3)=3./4.
      ro(1)=0.
      ro(2)=-17./60.
      ro(3)=-5./12.
      write(6,100) (ga(n),n=1,nsst),(ro(n),n=1,nsst)
  100 format(10x,'third order runge gutta=',4x,
     1       'gam=',3f8.3,4x,'ro',3f8.3)
      if(iforth.eq.0) 
     1     write(6,*)'second order Arakawa with J=J1'
      if(iforth.eq.1) 
     1     write(6,*)'forth order Arakawa with J=2J1-J2'
      if(iforth.eq.2) 
     1     write(6,*)'second order Arakawa with J=J2'
      go to 12
   10 continue
      ga(1)=1.5
      ga(2)=0.
      ga(3)=0.
      ro(1)=-.5
      ro(2)=0.
      ro(3)=0.
      write(6,110) ga(1),ro(1)
  110 format(10x,'adams bashfort=',4x,
     1       'gam=',f8.3,4x,'ro',f8.3)
   12 continue
      call solve
      stop
      end
c
c     subroutine solve  *******************************************
c
      subroutine solve
c     code for computation of two-dimensional incompressible flows
c     in cartesian coordinates.  with vorticity stream function.
c
c     this code solves flow fields periodic in x1 and periodic
c     in x2
c     This code permits to test arakawa schemes. 
c     The time discretization is the third order runge kutta.
c
      include 'param.f'
      dimension vor(m1,m2),psi(m1,m2)
      dimension ru(m1,m2)
      common/mesh/dx1,dx1q,dx2,dx2q
      common/coor/yp1(m1),yp2(m2)
      common/dim/n1,n1m,n2,n2m
      common/d1/tfin,alx1f,alx1i,alx2i,alx2f
      common/d2/nstop,nprint,ntst,npin,nread,nwrit,nprspe
      common/tstep/dt
      common/ntstep/nkt
      common/tscoe/ga(3),ro(3),nsst
      common/chnlc/chal,chbe,chga
      common/visct/re
      common/inmodo/inmod
      common/camo/enerv,vorip,vorim,ensti,cflm
      common/pale/palen
      common/hioo/oor3,oor4,oor6,oor8,oor10   
      common/hioo0/oo0r3,oo0r4,oo0r6,oo0r8,oo0r10   
      common/camo0/vmi0,ens0,ene0,circm0,pal0
      common/phco/iph
      common/icflp/icflm,cflma,dtl
      common/tfinst/tprin,tpin,tprspe
      iph=0
      pi=acos(-1.)
c
c     step  and messh sizes calculations
c
      call meshes
  758 format(1x,2i3)
      write(6,754)n1,n2
  754 format(10x,'n1=',i3,2x,'n2=',i3)
      write(6,755) dx1,dx2,dt,ntst,re
  755 format(3x,'dx1=',e10.3,3x,'dx2=',e10.3,3x
     1      ,'dt=',e10.3,3x,'ntst=',i4,3x,'re=',e10.3)
      call indic
      call coordi
      time=0.
      ntime=0
      ntii=0
c
c     read or creat initial conditions
c
      call phini
      call psini
      call initia(vor,psi,ru)
      ntstf=ntii+ntst
      ntii=ntii+1
      write(6,811) nprint,ntii,ntstf,dt
cc    write(nfil,811) nprint,ntii,ntstf,dt
  811 format (3x,'check initial conditions : nprint =',i8,
     1       '  ntii =',i8,'  ntstf =',i8,2x,' dt =',e11.3//)
c   evaluation of certain global quantities
c
      call vorqua(ntime,vor,psi,time,vormax,vormin)
      call cfield(vor,psi)
      circm0=vorip
      vmi0=vormax
      pal0=palen
      ens0=ensti
      ene0=enerv
      oo0r4=oor4
      oo0r6=oor6
      oo0r8=oor8
      oo0r10=oor10
      write(6,815) vorip,vormax,ensti,enerv
  815 format(3x,'in cir,vma,enst,ene',2x,4e12.4)
c
c   write the spectra and the fields
c
      call wrispe(ntime,time,vor,psi)
      call writfi(ntime,time,vor,psi,ru)
c
c
      write(6,711)nprint,ntst,nstop,dt,chal,chbe,chga
  711 format(1x,'in cond',3i8,2x,e11.3,/3x,'chal=',e8.2,3x,'chbe=',e8.2,
     1       3x,'chga=',e8.2)
  181 continue
      cflm=0.
      if(icflm.eq.1) then
      nstop=100000
                     endif
c
c    the time dependent calculation starts
c
      write(20,717)n1,n2,re,vorim
  717 format(1x,2i4,3x,2e15.4)
      call outth(ntime,time,vor,psi,vormax,vormin)
      do 350 ntime=ntii,ntstf
      if(ntime.gt.nstop) then
       write(6,*) ntime,' stop  >>>> ntime'
      stop
       endif
c
c  ********* calculation of the vorticity and stream function
c
      do 2000 ns=1,nsst
      alp=(ga(ns)+ro(ns))
      gam=ga(ns)
      rho =ro(ns)
      call tschem(vor,psi,ru,alp,gam,rho)
 2000 continue
      time=time+dt
      call vorqua(ntime,vor,psi,time,vormax,vormin)
      if(icflm.eq.1) then
c
c   courant number calculation for concentric cilinders
c
      call cfl(psi,cflmm)
c     print*,'88888******   dt=',dt,'  cflmm',cflmm
      dt=cflma/cflmm
      cflm=cflma
      if(dt.gt.dtl) dt=dtl
      if(dt.lt..1e-05)      then
      write(6,368)dt
  368 format(3x,'end calcul for dt=',e12.4)
      stop
                            endif
                      endif
c
c     print some quantity ,total energy,total enstrophy
c     and others
c
      if(icflm.eq.1) then
      if(amod(time,tpin).lt.dt) go to 306
      go to 305
                     else
      ntim=ntime-1
      if(mod(ntime,npin).eq.0) then
      go to 306       
                                endif
      go to 305
                     endif
  306 continue
      call cfield(vor,psi)
      write(96,252) enerv
  252 format(20x,'integrale energia =',e12.3)
      call outth(ntime,time,vor,psi,vormax,vormin)
  158 format(1x,i4,1x,e10.3,1x,10(1x,e10.4))
  305 continue
      if(icflm.eq.1) then
      if(amod(time,tprspe).lt.dt) go to 406
      go to 405
                     else
      ntim=ntime-1
      if(mod(ntime,nprspe).eq.0) then
      go to 406       
                                endif
      go to 405
                     endif
  406 continue
      call wrispe(ntime,time,vor,psi)
  405 continue
c
c
c    write the flow field
c
      if(icflm.eq.1) then
      if(amod(time,tprin).lt.dt) go to 301
            go to 300
                     else
      if(mod(ntime,nprint).eq.0) go to 301
            go to 300
                     endif
  301 continue
      call writfi(ntime,time,vor,psi,ru)
  300 continue
            if(time.gt.tfin.and.icflm.eq.1) then
      call wrispe(ntime,time,vor,psi)
      call writfi(ntime,time,vor,psi,ru)
      write(6,367)time
  367 format(3x,'end calcul for t=',e12.4)
      stop
                                            endif
  350 continue
  156 continue
      write(6,159) ntime
  159 format(/20x,'calculation stop for ntime=',i5)
      return
      end

