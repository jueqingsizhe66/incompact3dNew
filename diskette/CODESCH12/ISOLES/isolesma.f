c
      program main
      common/d1/re
      common/d2/nstop,nprint,ntst,npin,npstf
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/wrre/nwrit,nread
      common/tscoe/ga(3),ro(3),nsst
      common/d123/alx1,alx2,alx3
      common/averou/iav
      common/inior/indrea,icosma,icont
      common/ledat/ics0,cvisc,ifiltr,ibox
      common/cosma/csma,pratu
      common/tstep/dt,beta,ren
      common/cflco/icfl,cflc,tpin,tprin,tfin
      common/rot/f0
      common/ispec/imic
      common/vpeini/vper,omtres
      common/spemw/akkpp,qq,sig
      common/strat1/istrat,rho0,g,schm
      common/strat2/igrad,bvais
      common/kpmasp/kmax
      common/mtime/multim
      common/resca/iresca,tresca,rlamas 
      common/icomp/irunpc
      common/stratf/rich
c
      open(15,file='isoles.d')
      read(15,*) n1,n2,n3,nsst,alx1d,alx2d,alx3d
      read(15,*) nwrit,nread,iav,multim
      read(15,*) re,vper,dt,ntst,nprint,npin ,npstf
      read(15,*) nstop
      read(15,*) icfl,cflc,tpin,tprin,tfin
      read(15,*) ics0,ifiltr,csma,pratu,ibox
      read(15,*) f0
c add rotation to the momentum equation
c coriolis parameter f0, rotation about z-axis i.e.  axis 3
c try writing in inverse sec.
      read(15,*) imic
      read(15,*) akkpp,qq,sig
      read(15,*) istrat,rho0,g,schm
      read(15,*) igrad,bvais
      read(15,*) iresca,tresca,rlamas 
      read(15,*) irunpc
c
c it is conventient to read in bvais in cph. Thus we
c convert to inverse seconds.
c     bvais=bvais*2.*pi/(60.*60.)
c if istrat=1, evolve the density field rho
c (this is actually the perturbation to
c the mean density field which is rhobar
c --- here we only need use d/dz rhobar which
c we call grbar (but only after initial setup so that dimensions
c n1m etc are defined first!), and the vertical
c average int rhobar dz = rho0 in the equations)
c if igrad=1: the gradient is constant and we
c calculate grbar from the Brunt-Vaisala frequency (bvais)
c if igrad=-1: the gradient is not constant and we read in
c grbar from a file grbar.in
c
      cvisc=1./re
      rich=bvais
      pi=2.*asin(1.)
      n1m=n1-1
      n2m=n2-1
      n3m=n3-1
       n1mhp=n1m/2+1
       n2mhp=n2m/2+1
       n3mhp=n3m/2+1
       ikma=max(n1mhp,n3mhp)
       kmax=max(ikma,n2mhp)
      alx1=alx1d*pi
      alx2=alx2d*pi
      alx3=alx3d*pi
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
      write(6,112)alx1d,alx2d,alx3d
  112 format(10x,'box dimension',3x,'lx=',f4.2,'*pi'
     1        ,3x,'lz=',f4.2,'*pi',3x,'lz=',f4.2,'*pi')
      write(6,*) '**************************************'
      write(6,200)bvais
  200 format(10x,'3d isotropic turbulence. with strat,bvais=',e11.4)
      write(6,*) '**************************************'
      call openin
      if(igrad.ne.0) call mgrbar
      call solve
      stop
      end
c
      subroutine openin
      common/names/filcnw,filcnr,filth,filou,filuu,filsf
      character*17 filcnw,filcnr,filth,filou,filuu,filsf
      open(46,file='nfisodyn')
c     continuation file
      read(46,'(a)')filcnw
c     restart file
      read(46,'(a)')filcnr
c     time history
      read(46,'(a)')filth
c     rms file
      read(46,'(a)')filou
c     skewness
      read(46,'(a)')filuu
c     flatness
      read(46,'(a)')filsf
      open(23,file=filcnr,form='unformatted')
      rewind 13
      rewind 23
      return
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
c     this code solves flow fields periodic in x1 x2 and x3 
c     the equations for yhe cartesian component u(i) i=1,3
c     are discretized by a finite difference scheme.
c
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
c     THIS CODE PERMITS TO TEST SEVERAL s.g.s. MODELS AND THUS
C     PERFORM LES SIMULATIONS.
C     CASES WITHOUT STRATIFICATION OR BUOYANCY CAN BE CONSIDERED
C     THE INFLUENCE OF THE s.g.s MODELS CAN BE CHECKED
C     THE USER CAN EASILY INTRODUCE DIFFERENT s.g.s models
c
      include 'param.f'
      parameter (m1m=m1-1)
c
      common/wrre/nwrit,nread
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      dimension y(m2)
      dimension ru(ndv,m1,m2,m3)
      dimension q(ndv,m1,m2,m3),qcap(m1,m2,m3)
     1         ,dph(m1,m2,m3),pr(m1,m2,m3)
      dimension rho(m1,m2,m3),ru4(m1,m2,m3)
      dimension srho(0:m3)
      common/strain/st(m1,m2,m3,6)
      common/stsgs/stmed(6),strms(6)
      common/velmax/vmax(ndv),vmaxo(ndv)
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/d1/re
      common/d123/alx1,alx2,alx3
      common/d2/nstop,nprint,ntst,npin,npstf
      common/tstep/dt,beta,ren
      common/tscoe/ga(3),ro(3),nsst
      common/inener/ene0
      common/averou/iav
      common/names/filcnw,filcnr,filth,filou,filuu,filsf
      common/jrrr/jri,jrf,djr,irejr,iruuca
      common/ledat/ics0,cvisc,ifiltr,ibox
      common/sma/vis(m1,m2,m3)
      common/newdat/icost,timeav
      common/cflco/icfl,cflc,tpin,tprin,tfin
      common/rot/f0
      common/cosma/csma,pratu
      common/strat1/istrat,rho0,g,schm
      common/strat2/igrad,bvais
      common/resca/iresca,tresca,rlamas 
      common/ispec/imic
      common/icomp/irunpc
      common/numond/kx(m1),ky(m2),kz(m3),kkmax
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
c
      write(6,754)n1,n2,n3
  754 format(10x,'centered velocities',2x
     1      ,5x,'n1=',i2,2x,'n2=',i3,2x,'n3=',i2)
      write(6,755) re,dx1,dx2,dx3,dt,ntst
  755 format(3x,'re=',e10.3,3x,'dx1=',e10.3,3x,'dx2=',e10.3,3x
     1      ,3x,'dx3=',e10.3,3x,
     1      'dt=',e10.3,3x,'ntst=',i5,3x)
      write(6,756)  f0
  756 format(3x,'f0=',e10.3)
      if(istrat.eq.1) then
      write(6,757) istrat,rho0,g,schm
  757 format(3x,'istrat=',i3,3x,'rho0=',e14.7,
     1       3x,'g=',e10.4,3x,'schm=',e10.3)
      write(6,758) igrad,bvais
  758 format(3x,'igrad=',i3,3x,'bvais=',e14.7)
                      endif
      ren=re
      if(icfl.eq.1) then
      ntst=100000
      nstop=100000
      dtl=dt
      write(6,*)' calculation at dt variable by fixing cfl=',cflc
                     else
      write(6,765) dt
c     write(32,765) dt
  765 format(' time step dt= ',e10.3)
                     endif
c
c   calculates the quantities necessary for fft's
c
      call phini
      call speini
      ren=re
      nti=0
      time=0.
      ntii=0
      n1mh=n1m/2+1
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
      nav=1
c
      do i=1,n1m
      write(19,*)' wave numbers',i,kx(i),ky(i),kz(i)
      enddo
      call initur(q,pr,qcap,dph)
      call enerca(q,pr,ene0,time)
      call oldqua
      if(istrat.eq.1) then
      if(igrad.ne.0) then
      call inchas(rho,qcap,dph,q)
                     else
      call inlero(rho,qcap,dph)
                      endif
      call outrho(ntime,time,q,rho)
      call ouprho(time,nav,rho,qcap)
                      endif
      call outh(ntime,time,cflm,ene0)
      if(ene0.lt..1e-10) ene0=.1e-10
      call outpf(time,ene0,nav,q,qcap)
      else
c
      call inirea(ntii,time,q,rho,pr)
      print *,'leggo la soluzione'
      call divgck(q,qmax)
      print *,'divergenza della soluzione letta=',qmax
      call vmaxv(q)
      print *,'vmax(1)=',vmax(1)
      print *,'vmax(2)=',vmax(2)
      print *,'vmax(3)=',vmax(3)
      if(icfl.eq.1) then
c
c   courant number calculation 
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
c
c   this will be used only for statistical steady state simulations
c
      call tavrea(nav)
                   endif
c
      call enerca(q,pr,ene0,time)
      endif
c
c  ********* start of time dependent calculation ***
c
      ntstf=ntii+ntst
      ntii=ntii+1
      write(6,711)nprint,ntii,ntstf,nstop,nti,dt
  711 format(3x,'check in cond',5i8,2x,e11.3)
  181 continue
c
c    the dipendent calculation start
c
      beta=dt/re*0.5
c
      print *,'ics0',ics0
      if(ics0.eq.0) then
       print *,'calcolo con modello di large eddy inibito'
       do 320 k=1,n3m
       do 320 j=1,n2m
       do 320 i=1,n1m
       vis(i,j,k)=cvisc
 320  continue
      end if 
      if(ics0.eq.1) then
       print * ,' LES with dynamic model'
       if (ifiltr.eq.0) print * ,'physical filter '
       if (ifiltr.eq.1) print * ,'wave number filter '
      end if
      if(ics0.eq.2) then
       print * ,'LES  dynamic Wong lilly modified'
       print * ,'filtro fisico '
      end if
      if(ics0.eq.3) then
       print * ,'Smagorinsky  LES '
      end if
      if(ics0.eq.4) then
       print * ,' LES Structure function model'
       if(istrat.eq.1) then
       write(6,*)' program stops because this model has not'
       write(6,*)' been modified for the eddy diffusivity  '
       stop
                       endif
       call indstf
      end if
      if(ics0.eq.5) then
       print * ,' LES constant nu turb. =<nu_T(i,j,k)> '
      end if
      ichang=0
  351 continue
      if(irunpc.eq.1) then
      call openfi(time)
                      endif
      if(ichang.eq.1) then
      write(6,*)' velocity  rescaled and the simulation restart'
      write(6,*)' icfl,dt,ren ',icfl,dt,ren 
      call enerca(q,pr,ener,time)
      call outh(ntime,time,cflm,ener)
      call outpf(time,ener,1,q,qcap)
                      endif
      ncount=0
      do 350 ntime=ntii,ntstf
      if(ntime.gt.nstop) go to 167
      if(icfl.eq.1)                 then
c
c   courant number calculation for running at cfl=constant
c
      call cfl(q,cflum)
      cflmm=cflum
      cflmp=cflum*dt
      dt=cflc/cflmm
      cflm=cflc
      if(dt.gt.dtl) then
      dt=dtl
      beta=dt/ren*0.5
                    endif
      if(dt.lt..1e-02)      then
      write(6,368)dt
  368 format(3x,'end calcul for dt=',e12.4)
        call contwr(ntime,time,q,rho,pr,enen)
      call outpf(time,enen,nav,q,qcap)
      if(istrat.eq.1) then
      call ouprho(time,nav,rho,qcap)
      call outrho(ntime,time,q,rho)
                      endif
      stop
                            endif
                                   endif
      call quales
c
c  if nsst=1 adams bashford if nsst=3 runge-kutta
c
      do 2000 ns=1,nsst
c
c    computation of turbulent viscosity only at the 
c    first substep
cccccccccccccccccccccccccccccccccccccccccccccccc
c
c Dynamic model
c
cccccccccccccccccccccccccccccccccccccccccccccccc
      if(ns.eq.1.and.ics0.eq.1)  then 
        ncount=ncount+1
cccccccccccccccccccccccccccccccccccccccccccccccc
c
c Dynamic model:  Lilly version without density
c
cccccccccccccccccccccccccccccccccccccccccccccccc
      if(istrat.eq.0) then
        call dynamic(q,dph,ncount,time)
                      endif
      if(istrat.eq.1.and.igrad.eq.0) then
        call dynamic(q,dph,ncount,time)
                      end if
      if(istrat.eq.1.and.igrad.eq.1) then
      if(time.gt.0.5) then
        call dywolil(q,dph,rho,ncount,time)
                     else
        call smarho(q,dph,rho,ncount,time)
                     end if
                                     end if
      end if
      if(ns.eq.1.and.ics0.eq.2)  then 
        ncount=ncount+1
cccccccccccccccccccccccccccccccccccccccccccccccc
c
c Dynamic model:  Modified Lilly version 
c
cccccccccccccccccccccccccccccccccccccccccccccccc
      if(time.gt.0.05) then
        call dywolmo(q,dph,rho,ncount,time)
                     else
        call smarho(q,dph,rho,ncount,time)
                     end if
                                  end if
cccccccccccccccccccccccccccccccccccccccccccccccc
c
c Smagorinsky model
c
cccccccccccccccccccccccccccccccccccccccccccccccc
       if(ns.eq.1.and.ics0.eq.3)  then
        ncount=ncount+1
      if(istrat.eq.0) then
        call smago(q,dph,ncount,time)
                      endif
      if(istrat.eq.1.and.igrad.eq.0) then
        call smago(q,dph,ncount,time)
                      end if
      if(istrat.eq.1.and.igrad.eq.1) then
        call smarho(q,dph,rho,ncount,time)
                      end if
      end if
cccccccccccccccccccccccccccccccccccccccccccccccc
c
c Structure function Metais Lesieur
c
cccccccccccccccccccccccccccccccccccccccccccccccc
       if(ns.eq.1.and.ics0.eq.4)  then
        ncount=ncount+1
        call strfun(q,dph,ncount,time)
      end if
cccccccccccccccccccccccccccccccccccccccccccccccc
c
c Constant nu turb. =<nu_T(i,j,k)>
c
cccccccccccccccccccccccccccccccccccccccccccccccc
       if(ns.eq.1.and.ics0.eq.5)  then
        ncount=ncount+1
      if(istrat.eq.0) then
        call consma(q,dph,ncount,time)
                      endif
      if(istrat.eq.1.and.igrad.eq.0) then
        call consma(q,dph,ncount,time)
                      end if
      if(istrat.eq.1.and.igrad.eq.1) then
        call cosmrh(q,dph,rho,ncount,time)
                      end if
      end if
c
       if(ns.eq.1)  then
      call strper(q)
      do n=1,6
      stmed(n)=0.
      do k=1,n3m
      do j=1,n2m
      do i=1,n1m
      stmed(n)=st(i,j,k,n)*(vis(i,j,k)-cvisc)+stmed(n)
      enddo
      enddo
      enddo
      stmed(n)=stmed(n)/float(n1m*n2m*n3m)
      enddo
      write(68,118)time,(stmed(n),n=1,6)
  118 format('sma',e10.3,6e12.5)
                     endif
      alp=(ga(ns)+ro(ns))
      gam=ga(ns)
      rhh=ro(ns)
      call tschem(q,pr,alp,gam,rhh,ru,qcap,dph,rho,ru4) 
 2000 continue
c
      time=time+dt
c
c   in the Compte-Bellot case there is the possibility
c   to restart the simulation with a rescaled energy
c   to mach Rlambda
c
      if(imic.eq.1.and.iresca.eq.1.and.time.ge.tresca) then
      call rescal(q,pr,qcap,dph)
      if(irunpc.eq.1) then
      call closefi
                      endif
      iresca=0
      time=0.
      ntii=1
      ichang=1
      go to 351
                                                       endif
      if(icfl.eq.1)then
      if(cflm.gt.cflmm) then
      istop=1
      go to 306
                     endif
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
      go to 305
                     endif
  306 continue
c
c
c     print some quantity cf,max(divg),total energy,amx(veloc),
c     and theirs gradient in time.
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
c
c     total energy
c
      call enerca(q,pr,enen,time)
c
      if(iav.eq.1) then
       nav=nav+1
       call taver
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
      if(istrat.eq.1) then
      call outrho(ntime,time,q,rho)
                      endif
c
  305 continue
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
c     write(6,*)'ntime=',ntime,' prima  contwr'
        call contwr(ntime,time,q,rho,pr,enen)
c     write(6,*)'ntime=',ntime,' dopo  contwr'
        if(iav.eq.1) then
         call tavwri(nav)
c
                     endif
                           endif
      if(iav.eq.0) nav=1
      call outpf(time,enen,nav,q,qcap)
      if(istrat.eq.1) then
      call ouprho(time,nav,rho,qcap)
                      endif
  300 continue
       if(time.gt.tfin.and.icfl.eq.1) then
      write(6,367)time
  367 format(3x,'end calcul for t=',e12.4)
c       call contwr(ntime,time,q,rho,pr,enen)
      stop
      endif
  350 continue
      go to 167
  166 continue
      write(6,168)
  168 format(10x,'calculation diverged in time for the vel field')
      write(6,*)'  vmax(1,2,3 ',(vmax(l),l=1,3)
      write(6,7698) qmax,ntime,cflm
      if(iav.eq.0) nav=1
      call outpf(time,enen,nav,q,qcap)
      if(istrat.eq.1) then
      call ouprho(time,nav,rho,qcap)
                      endif
      go to 167
  169 continue
      write(6,178)qmax
  178 format(10x,'calculation diverged for  dmax='
     1 ,e12.5)
  167 continue
c     call spectre(q)
c     call spevis(vis)
      return
      end
c***********************************************************
      subroutine speini
      include 'param.f'
      common/fft001/ifxx1(13),trigxx1(2*(m1-1))
      common/fft002/ifxx2(13),trigxx2(3*(m2-1)/2+1)
      common/fft003/ifxx3(13),trigxx3(2*(m3-1))
      common/dim/n1,n1m,n2,n2m,n3,n3m
      nx3fft=n3m
      call cftfax(nx3fft,ifxx3,trigxx3)
      nx1fft=n1m
      call cftfax(nx1fft,ifxx1,trigxx1)
      nx2fft=n2m
      call fftfax(nx2fft,ifxx2,trigxx2)
      return
      end

