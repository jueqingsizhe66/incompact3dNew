
c
      program main
      common/d1/re
      common/d2/nstop,nprint,ntst,npin,npstf
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/wrre/nwrit,nread
      common/tscoe/ga(3),ro(3),nsst
      common/d123/alx1,alx2,alx3
      common/averou/iav
      common/ledat/cvisc
      common/cflco/icfl,cflc,tpin,tprin,tfin,twrfi
      common/rot/f0
      common/ispec/imic
      common/vpeini/vper,omtres
      common/spemw/akkpp,qq,sig
      common/tstep/dt,beta,ren
      common/pardip/thet0,vsi,yc1mo,yc2mo,akmo,velmo
      common/anmod/ar1,ar2
      common/pargam/gell1,gell2,bw1,bw2

      pi=2.*asin(1.)
c
      open(15,file='iso.d',status='unknown')
      read(15,*) n1,n2,n3,nsst,alx1d,alx2d,alx3d
      write(96,*) n1,n2,n3,nsst,alx1d,alx2d,alx3d
      read(15,*) nwrit,nread,iav
      write(96,*) nwrit,nread,iav
      read(15,*) re,vper,dt,ntst,nprint,npin ,npstf
      write(96,*) re,vper,dt,ntst,nprint,npin ,npstf
      read(15,*) nstop
      write(96,*) nstop
      read(15,*) icfl,cflc,tpin,tprin,tfin,twrfi
      write(96,*) icfl,cflc,tpin,tprin,tfin,twrfi
      read(15,*) f0
      write(96,*) f0
c add rotation to the momentum equation
c coriolis parameter f0, rotation about z-axis i.e.  axis 3
c try writing in inverse sec.
      read(15,*) imic
      write(96,*) imic
      if(imic.ge.0) then
      read(15,*) akkpp,qq,sig
      write(96,*) akkpp,qq,sig
                    else
      read(15,*) yc1mo,yc2mo,ramo,velmo,vsi
      write(96,*) yc1mo,yc2mo,ramo,velmo,vsi
      read(15,*) thet0,omtres
      write(96,*) thet0,omtres
      read(15,*) ar1,ar2
      write(96,*) ar1,ar2
      read(15,*) gell1,gell2,bw1,bw2
      write(96,*) gell1,gell2,bw1,bw2
      akmo=3.83711/ramo
      pi=2.*asin(1.)
      thet0=thet0*pi/180.
      yc1mo=yc1mo*pi
      yc2mo=yc2mo*pi
      write(6,201)thet0
  201 format(3x,'modone vort stream funct. tht0=',e10.4)
                    endif
      cvisc=1./re
      n1m=n1-1
      n2m=n2-1
      n3m=n3-1
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
c     if nsst=1 Adams-Bashforth scheme
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
      write(6,200)
  200 format(10x,'3d isotropic turbulence.')
      write(6,*) '**************************************'
      call openfi
      call solve
      stop
      end
c
      subroutine openfi
      common/names/filcnw,filcnr,filth,filou,filuu,filsf
      character*17 filcnw,filcnr,filth,filou,filuu,filsf
      open(46,file='nfisodyn',status='unknown')
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
      open(23,file=filcnr,form='unformatted',status='unknown')
      open(32,file=filth,status='unknown')
      open(42,file=filou,status='unknown')
      open(52,file=filuu,status='unknown')
      open(62,file=filsf,status='unknown')
      open(72,file='momder.out',status='unknown')
      open(34,file='enstr.out',status='unknown')
      open(82,file='vorms.out',status='unknown')
      open(92,file='tarms.out',status='unknown')
      open(72,file='momder.out',status='unknown')
      open(53,file='circu.out',status='unknown')
      open(33,file='traje.out',status='unknown')
      open(83,file='vomax.out',status='unknown')
      open(66,file='rlack.out',status='unknown')
      open(67,file='vmax.out',status='unknown')
      open(68,file='charith.out',status='unknown')
      rewind 13
      rewind 23
      rewind 32
      rewind 42
      rewind 52
      rewind 62
      rewind 72
c     close(32)
c     close(42)
c     close(52)
c     close(62)
c     close(72)
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
c
      include 'param.f'
      parameter (m1m=m1-1)
c
      common/cflco/icfl,cflc,tpin,tprin,tfin,twrfi
      common/wrre/nwrit,nread
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      dimension q(ndv,m1,m2,m3),qcap(m1,m2,m3)
     1         ,dph(m1,m2,m3),pr(m1,m2,m3)
      dimension ru(ndv,m1,m2,m3)
      common/velmax/vmax(ndv),vmaxo(ndv)
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/d1/re
      common/d123/alx1,alx2,alx3
      common/d2/nstop,nprint,ntst,npin,npstf
      common/tstep/dt,beta,ren
      common/tscoe/ga(3),ro(3),nsst
      common/inener/ene0
      common/averou/iav
      common/sc/sc,sc1
      common/names/filcnw,filcnr,filth,filou,filuu,filsf
      common/jrrr/jri,jrf,djr,irejr,iruuca
      common/ledat/cvisc
      common/newdat/icost,timeav
      common/rot/f0
      common/speene/e(ndv,0:m3)
      common/numond/kx(m1),ky(m2),kz(m3),kkmax
      common/ispec/imic
      dimension vor(m1,m2)
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
      call meshes
      call indic
      call coordi
c
      write(6,754)n1,n2,n3
  754 format(10x,'centered velocities',2x
     1      ,5x,'n1=',i2,2x,'n2=',i3,2x,'n3=',i2)
      write(6,755) re,dx1,dx2,dx3,dt,ntst
  755 format(3x,'re=',e10.3,3x,'dx1=',e10.3,3x,'dx2=',e10.3,3x
     1      ,3x,'dx3=',e10.3,3x,
     1      'dt=',e10.3,3x,'ntst=',i5)
      write(6,756)  f0
  756 format(3x,'f0=',e10.3)
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
      write(6,*) '  '
c     write(32,*) '  '
c     open(67,file='vma.out')
      nti=0
      time=0.
      ntii=0
      n1mh=n1m/2+1
c
c     initial conditions
c 
      if(n3m.gt.1) then
      call speini
                   endif
      call phinip
      call phini
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
c
      if(imic.gt.0) then
c
c    turbulence
c
      call initur(q,pr,qcap,dph,time)
      call spectre(q,qcap)
       open(77,file='spec.0000',status='unknown')
       do k=1,kkmax
       write (77,*)k,e(1,k),e(2,k),e(3,k)
       end do
       close(77)
        if(tprin.lt.0.and.tfin.ge.0) then
        call contwr(ntime,time,q,enen)
       stop
                       endif
      call oldqua
                    else
c
c   vortex dynamics
c
      call inidip(q,pr,qcap,dph)
        call contwr(0,time,q,enen)
                    endif
c     call outh(ntime,time,cflm,ene0,q,pr,vor)
c     call outpf(time,ene0,nav,q,qcap,pr)
      if(ene0.lt..1e-10) ene0=.1e-10
      else
c
c   from restarting file
c
      call inirea(ntii,time,q)
      print *,'leggo la soluzione'
      call divgck(q,qmax)
      print *,'divergenza della soluzione letta=',qmax
      call outh(ntime,time,cflm,ene0,q,pr,vor)
c     call outpf(time,enen,nav,q,qcap,pr)
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
      ntime=ntii
      itim=time+0.3
      write(pntim,77)itim
   77 format(i4.4)
c
      call outh(ntime,time,cflm,ene0,q,pr,vor)
      endif
      call vmaxv(q)
      write(6,*)'  vmax',(vmax(l),l=1,3)    
c     open(57,file='vma.out')
      write(67,*)time,(vmax(l),l=1,3)    
c     close(57)
c
c  ********* start of time dependent calculation ***
      ntstf=ntii+ntst
      ntii=ntii+1
      write(6,711)nprint,npin,ntii,ntstf,nstop,nti,dt
  711 format(3x,'check in cond',6i8,2x,e11.3)
  181 continue
c
c    the dipendent calculation start
c
c
       print *,'calcolo DI DIRECT SIMULATION             '
      ncount=0
      do 350 ntime=ntii,ntstf
c     if(time.gt.100.) iav=1
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
        call contwr(ntime,time,q,enen)
      call outpf(time,enen,nav,q,qcap,pr)
      stop
                            endif
                                   endif
c
c  if nsst=1 adams bashford if nsst=3 runge-kutta
c
      do 2000 ns=1,nsst
c
      alptmp=(ga(ns)+ro(ns))
      gamtmp=ga(ns)
      rhotmp=ro(ns)
      call tschem(q,pr,alptmp,gamtmp,rhotmp,
     &              ru,qcap,dph) 
 2000 continue
c
      time=time+dt
c     write(6,*)ntime,time
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
      call vmaxv(q)
      vmaa=0.
      do l=1,3
      vmaa=max(vmaa,vmax(l))
      enddo
      if(vmaa.gt.10.) go to 166
c
c     the calculation stop if the velocities are diverging for numer
c     stab conditions (courant number restrictions)
c
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
c     open(66,file='divma.out')
      write(66,7698) ntime,qmax,cflm,dt
c     close(66)
c     open(57,file='vma.out')
      write(67,*)time,(vmax(l),l=1,3),qmax,dt
c     close(57)
 7698 format(i5,6e13.5)
c
      call outh(ntime,time,cflm,enen,q,pr,vor)
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
        call contwr(ntime,time,q,enen)
        if(iav.eq.1) then
         call tavwri(nav)
c
                     endif
                           endif
      call outpf(time,enen,nav,q,qcap,pr)
  300 continue
       if(time.gt.tfin.and.icfl.eq.1) then
      write(6,367)time
  367 format(3x,'end calcul for t=',e12.4)
        call contwr(ntime,time,q,enen)
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
      call outpf(time,enen,nav,q,qcap,pr)
      go to 167
  169 continue
      write(6,178)qmax
  178 format(10x,'calculation diverged for  dmax='
     1 ,e12.5)
  167 continue
      return
      end
c***********************************************************
      subroutine speini
c
c   evaluates quantities necessary for FFT
c   before to initiate the simulations
c
      include 'param.f'
      common/fft002/ifxx2(13),trigxx2(2*(m2-1))
      common/fft001/ifxx1(13),trigxx1(3*(m1-1)/2+1)
      common/fft003/ifxx3(13),trigxx3(2*(m3-1))
      common/dim/n1,n1m,n2,n2m,n3,n3m
      nx3fft=n3m
      nx2fft=n2m
      nx1fft=n1m
      call cftfax(nx3fft,ifxx3,trigxx3)
      call cftfax(nx2fft,ifxx2,trigxx2)
      call fftfax(nx1fft,ifxx1,trigxx1)
      return 
      end
