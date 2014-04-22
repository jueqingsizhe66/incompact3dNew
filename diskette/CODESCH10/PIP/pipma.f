c************************************************************************
c                                                                       *
c                                                                       *
c     this cod is mad for simulating three-dimensional flows in polar *
c     ( cilindrical ) coordinates.                                      *
c     boundary condition are slip-walls in the radial direction and     *
c     periodic in the axial direction.                               *
c                                                                       *
c      navier-stokes equations are solved by a fractional step method   *
c     ( Kim and Moin ) with the pressure in the first step.             *
c                                                                       *
c     The time advancement of the solution is obtained by a             *
c     Runge-Kutta 3rd order low storage scheme (Wray) or a 2nd          *
c     order Adams-Bashfort scheme.                                      *
c                                                                       *
c     The Poisson  equation for the pressure is solved directly         *      
c     introducing FFT in the azimutal and vertical direction.           *
c                                                                       *
c                 roberto verzicco and paolo orlandi                    *  
c                 dipartimento di meccanica ed aeronautica              *
c                 universita' la sapienza di roma                       *
c       the method is explained in                                      *
c       Verzicco, R. \& Orlandi, P. \\1996a\\
c       A finite-difference scheme for the three-dimensional
c       incompressible flows in cylindrical coordinates\\
c       J. Comp. Phys\\123\\402--414 .
c                                                                       *
c                                                                       *
c     All variables are calculated in a staggered grid:                 *
c                                                                       *
c        Instead of velocities, flux variables are introduced to avoid  *
c        the singularity for r=0                                        *
c        q1=vtheta*r, q2= r*vr, q3= vz                                    *      
c                                                                       *
c        qcap :divergence of the  non free divergent velocity field     *
c                                                                       *
c        non linear terms:                                              *
c                                                                       *
c        ru1, ru2, ru3 : old step                                       *
c        h1,h2,h3 : new step                                            *
c        pr : pressure                                                  *
c        dph : pressure correction                                      *
c       pressure solver is in the package pipph*.f                       *
c       pipphcra.f     for CRAYS computers
c       pipphesl.f     for IBM workstations with ESSl library
c       pipphnc.f      for any computer Temperton's FFT
c       non linear terms are calculated in hdnl routines                *
c       the invertions of  equations are performed in invt....          *
c       routines                                                        *
c                                                                       *
c                                                                       *
c************************************************************************
      program main                                                      
      include 'param.f'
      common/npjet/n2t,n2v
      common/irety/ireq2,iwrq2
      open(15,file='pipe.d')
      open(5,file='iq2.d')
      read(15,301) dummy                                                
      read(15,*) n1,n2,n3,nsst,nwrit,nread
      read(15,301) dummy                                                
      read(15,*) n1p,n2p,n3p
      read(15,301) dummy                                                
      read(15,*) ntst,nprint,npin,npouth,ireset
      read(15,301) dummy                                                
      read(15,*) alx3p                                                  
      read(15,301) dummy                                                
      read(15,*) vper,dt
      read(15,301) dummy                                                
      read(15,*) cflc,icfl,tpin,tprin,dtfin,tpouth
      read(15,301) dummy                                                
      read(15,*) re,ros,ipipe,islip
c$$$$$ parameters for non uniform grid (r distribution) $$$$$$$$$$$$$$$$$
      read(15,301) dummy                                                
      read(15,*)strr,rext,rint,rmed1
      read(15,301) dummy                                                
      read(15,*)istr,rmed,etdp,strb,n2t
c$$$$$$$parameters for boundary conditions $$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      read(15,301) dummy
      read(15,*) iav,timav
      read(15,301) dummy
      read(15,*) ichrc
      read(15,301) dummy
      read(15,*)strro,rexto,rinto,rmed1o
      read(15,301) dummy
      read(15,*)istro,rmedo,etdpo,strbo,n2to
      read(15,301) dummy
      read(15,*)iwlop,njprs,(npjp(j),j=1,njprs)
c
c    by ipipe=0
c    perform the simulation of a dipole moving in
c    r-theta planes this to check axis accuracy
c
      if(ipipe.eq.0) then
      read(15,301) dummy
      read(15,*)yc1mo,yc2mo,velmo
                      endif
301   format(a4)                                                        
101   format(a60)                                                        
      r0=rmed
c
c    this quantities permit to write smaller restarting files
c
      write(6,*)' enter ireq2,iwrq2 for read and write q2'
      read(5,*) ireq2,iwrq2
c
      pi=2.*asin(1.)                                                    
      alx3d=2.*pi*alx3p
      n1m=n1-1                                                          
      n2m=n2-1                                                          
      n3m=n3-1                                                          
      n3mh=n3m/2+1                                                      
c                                                                       
c                                                                       
      alx3=alx3d                                                        
c                                                                       
      call openfi
c                                                                       
c     assign coefficients for time marching schemes                     
c
      if(nsst.gt.1) then                                                
      gam(1)=8./15.                                                      
      gam(2)=5./12.                                                      
      gam(3)=3./4.                                                       
      rom(1)=0.                                                          
      rom(2)=-17./60.                                                    
      rom(3)=-5./12.                                                     
      write(6,*)'*****************************************************'
      write(32,*)'*****************************************************'
      write(6,*)'*                                                    *'
      write(32,*) '*                                                  *'
      write(6,*)'*               TURBULENT PIPE WITH ROTATION         *'
      write(32,*)'*              TURBULENT PIPE WITH ROTATION         *'
      write(6,*)'*                                                    *'
      write(32,*) '*                                                  *'
      write(6,*)'*****************************************************'
      write(32,*)'*****************************************************'
      write(6,*)'  '
      write(32,*) '  '
      write(6,100) (gam(n),n=1,nsst),(rom(n),n=1,nsst)                    
      write(32,100) (gam(n),n=1,nsst),(rom(n),n=1,nsst)                    
  100 format(3x,'the time advancement of the solution is obtained by',/,
     1  'a third order Runge-Kutta scheme. Coefficients are:  ',/, 
     1       'gam=  ',3f8.3,6x,'ro=  ',3f8.3)                                
      else                                                              
      gam(1)=1.5                                                         
      gam(2)=0.                                                          
      gam(3)=0.                                                          
      rom(1)=-.5                                                         
      rom(2)=0.                                                          
      rom(3)=0.                                                          
      write(6,110) gam(1),rom(1)                                          
      write(32,110) gam(1),rom(1)                                          
  110 format(3x,'the time advancement of the solution is obtained by',/,
     1  'an Adams-Bashfort scheme. Coefficients are:  ',/, 
     1       'gam=  ',f8.3,6x,'ro=  ',f8.3)                                
      endif                                                             
c                                                                       
c                                                                       
      do 10 ns=1,nsst
      alm(ns)=(gam(ns)+rom(ns))
   10 continue
      write(6,*)'  '
      write(32,*) '  '
      write(6,112)alx3d,r0                                              
      write(32,112)alx3d,r0                                             
  112 format(10x,'domain dimensions L_z =',f8.5,' L_r = ',f8.5)                
      write(6,*)'  '
      write(32,*) '  '
      write(6,*)'  '
      write(32,*) '  '
c                                                                       
      call gcurv                                                        
                                                                        
      stop                                                              
      end                                                               
c************************************************************************
c************************************************************************
c                                                                       *
c                                                                       *
c    code for computation of three-dimensional incompressible flows     *
c    in polar cylindrical coordinates.                                  * 
c                                                                       * 
c    this code solves flow fields which are periodic in the             *
c    axial and azimuthal directions and with free-slip               *
c    or no-slip conditions in the radial direction                     *
c                                                                       * 
c    the non uniform grid in r  is given in the subroutine cordi        *
c                                                                       * 
c    the equations for q(i) i=1,2,3                                     *
c     q1=v(theta)*r  q2=v(r)*r     q3=v(zeta)                           *
c    are discretized by a finite difference scheme.                     *
c                                                                       *
c                                                                       *
c    all the spatial derivatives are discretized by centered schemes    *
c    including the non linear terms.                                    * 
c    the non-linear terms are computed explicitly the viscous terms     *    
c    implicitly.                                                        *
c    a factored scheme brings to the solution of three tridiagonal      *
c    matrices for each velocity component                               *
c    in time a fractional step is used in the version of nagi mansour   * 
c    introducing the pressure in the first step.                        *
c    the non-linear terms and the cross der. are discr. by explicit     *
c    adams-bashfort or 3rd order r-k wray method.                       *
c    the scalar ph which gives the div. free vel. is solved by a direct *
c    method. for the x3 and x1 directions modified wave numbers are used *
c    and a fft method gives  the ph in the physical space.     *
c    no boundary conditions are necess. for the poisson equation.       *
c               a*ph=div(q)                                             *
c                                                                       *
c    all the second derivatives are discretized by a second order conse *
c    rvative scheme.                                                     *
c                                                                       *
c                                                                       *
c************************************************************************
c************************************************************************
      subroutine gcurv                                                  
      include 'param.f'
      common/timw/timew
c
c                                                                       
c     grid definition, indices and mesh size calculation                          
c                                                                       
      call meshes
      call indic                                                        
      call cordin                                                       
      if(ipipe.eq.0) then
      call pricor
                     endif
c
c
c     print some informations on the run
c
      write(6,*) '  '
      write(32,*) '  '
      write(6,754)n1,n2,n3                                              
      write(32,754)n1,n2,n3                                             
  754 format(10x,'number of grid points :'/                             
     1      5x,'n1=',i4,2x,'n2=',i4,2x,'n3=',i4/)                       
      write(6,*) '  '
      write(32,*) '  '
      write(6,756)n1p,n2p,n3p                                              
      write(32,756)n1p,n2p,n3p                                             
  756 format(10x,'plotting stride  :'/                             
     1      5x,'n1p=',i4,2x,'n2p=',i4,2x,'n3p=',i4/)                       
      write(6,*) '  '
      write(32,*) '  '
      write(6,755) re,ros
      write(32,755) re,ros
  755 format(3x,' Parameters of the flow: ',/,
     1 ' Reynolds number = ',e10.3,3x,' Rossby number = ',e10.3) 
      write(6,*) '  '
      write(32,*) '  '
      if(icfl.eq.1) then
      ntst=100000
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
      ren=re                                                            
      time=0.                                                           
      ntii=0                                                            
      beta=dt/re*0.5                                                    
            nav=0                                                           
c*******************************************************
c
c    read or generate the initial fields
c
c*******************************************************
c                                                                       
c    create the initial fields                                       
c                                                                       
      call initia
c
c  evaluation of metric quantities for the inversion
c
      call phini    
      call coetar
      call csoq13i
      call csoq2i
      if(nread.eq.0) then                                               
            write(6,*)' nread=0 ==>  the calculation is starting from '      
            write(6,*)' generated initial conditions                  '      
            write(32,*)' nread=0 ==>  the calculation is starting from'      
            write(32,*)' generated initial conditions                 '      
            write(6,*) '  '
            write(32,*) '  '
            ntii=0                                                          
            ntt=0.       
            ntime=0                                                         
            ncount=0                                                         
            time=0.                                                         
            cflm=0.                                                         
            nav=0                                                           
            nap=0                                                           
            if(ipipe.eq.1) then
            write(6,*)'************ TURBULENT PIPE ************'
            call intur 
                           endif
            if(ipipe.eq.0) then
            write(6,*)'************ 2D DIPOLE      ************'
            call indipo 
c                                                                       
c   print the initial flow field                                        
c                                                                       
            call outhdi
            call oudip
                           endif
c                                                                       
      else                                                              
c*******************************************************
c                                                                       
c   read  the initial fields                                       
c                                                                       
          write(6,*)' nread=1 ==>  the calculation is starting from  '      
          write(6,*)' read initial conditions                   '      
          write(32,*)' nread=1 ==>  the calculation is starting from  '      
          write(32,*)' read initial conditions                   '      
          write(6,*) '  '
          write(32,*) '  '
          call inirea(ntii,time,ntt,ncount,nap)           
          call divuck(dmax,dtot)                                 
          write(6,900)dmax,dtot                 
          write(32,900)dmax,dtot                
  900   format(3x,'maxima local and global divergence of the read field',
     1     /,' dmax = ',e11.4,' dtot = ',e11.4)
       if(iav.eq.1) then
       call tavrea(nav,nap)
          write(6,*)' read the averaged quantity at the last time'
          write(6,*)' the averaged were done for nav=',nav
          write(6,*)' the averaged spectra and correl done for nap=',nap
                    endif
      if(icfl.eq.1) then
c
c   courant number calculation for the first step from a restarting file
c
      call cflu(cflum)
      cflmm=cflum
      dt=cflc/cflmm
      cflm=cflc
      if(dt.gt.dtl) then
      dt=dtl
      beta=dt/re*0.5
      cflm=cflmm*dt
                    endif
      write(6,*)' with inirea=1 dt=',dt
                   endif
      endif                                                             
c
c     print other informations
c
      timel=time
      tfin=time+dtfin
      ntstf=ntst                                                   
      write(6,*) '  '
      write(32,*) '  '
      write(6,711)nprint,ntii,npin,ntt,timel,tfin
      write(32,711)nprint,ntii,npin,ntt,timel,tfin
  711 format(3x,'check in conditions :',/,
     1'  nprint = ',i5,'  ntii = ',i5,'  npin = ',i5
     1,'   ntt=',i5,'  timel=',e11.4,'  tfin=',e11.4)
      write(6,*) '  '
      write(32,*) '  '
  159 format(1x,i4,2x,e10.4,3e10.3,3(1x,e10.4,1x,i3,1x,i3),e10.3)
c
c*******************************************************
c
      write(6,*)'           time ',
     1'         v1tot        v3tot         dp3ns',
     1'         utat        utap        cflm'
      write(32,*)'           time ',
     1'         v1tot        v3tot         dp3ns',
     1'         utat        utap        cflm'
c
c      write the time history
c
        ntii=ntii+1
c
c                                                                       
c  ********* starts the time dependent calculation ***                  
c                                                                       
      do 350 ntime=ntii,ntstf                                           
      timew=time                                                                  
      if(icfl.eq.1) then
c
c   courant number calculation 
c
      call cflu(cflum)
      cflmm=cflum
      dt=cflc/cflmm
      cflm=cflc
      if(dt.gt.dtl) then
      dt=dtl
      beta=dt/re*0.5
      call csoq13i
      call csoq2i
      cflm=cflmm*dt
                    endif
      cflummm=dt*cflum
      cflbmmm=dt*cflbm
c
c     write(6,*) cflummm,cflbmmm
c
      if(dt.lt..1e-05)      then
      write(6,368)dt
  368 format(3x,'end calcul for dt=',e12.4)
      stop
                            endif
                    endif
      time=time+dt                                                      
      call tschem(ntime,time)
c                                                                       
c  write velocity vorticity and pressure at the desired
c  locations at each time step.
c  The analisys is similar as that in the experiments
c
      if(iwlop.gt.0) then
      call outiqu(time)
                     endif
      if(time.gt.timav) iav=1
      if(icfl.eq.1) then
      if(amod(time,tpin).lt.dt) go to 306
      go to 305
                     else
      if(mod(ntime,npin).eq.0) then
      call cflu(cflum)
      cflmm=cflum
      cflm=cflmm*dt
           go to 306
                               endif
                     endif
      go to 305
  306 continue
c*******************************************************
c
c     print time history of   cf, max(vel), max(vor) etc.  etc 
c
c*******************************************************
        if(ipipe.eq.1) then
        call outhth(time,nap,ntime,cflm)            
                       endif
        if(ipipe.eq.0) then
        call outhdi
                       endif
        call vmaxv(n1m,n2m,n3m)
          call divuck(dmax,dtot)                                 
       write(34,112)time,dt,(vmax(l),l=1,3),dmax,dtot
  112 format(3x,10e12.5)
        ntt=ntt+1
c                                                                       
  305 continue                                                          
c*******************************************************
c                                                                       
c     write the flow field  as a restaring file for
c     postprocessing                                            
c                                                                       
c*******************************************************
      if(icfl.eq.1) then
      if(amod(time,tprin).lt.dt) go to 301
            go to 300
                     else
      if(mod(ntime,nprint).eq.0) go to 301
            go to 300
                     endif
  301 continue
            if(ipipe.eq.1)     then
            if(nwrit.eq.1)       then
             call contwr(ntime,time,ntt,nap)
                                 endif
           if(iav.eq.1)    then
c
c    here the time averages of statistics are done
c    and written
c
           call outh(time,nav,ntime,cflm)
           call tavwri(nav,time)
                           endif
                               endif
            if(ipipe.eq.0)     then
           call oudip
                               endif
  300 continue
       if(time.gt.tfin.and.icfl.eq.1) then
      write(6,367)time
  367 format(3x,'end calcul for t=',e12.4)
      stop
                                        endif
  350 continue                                                          
c*******************************************************
c                                                                       
c     end of the time dipendent calculation
c
c*******************************************************
      go to 167                                                         
  166 continue                                                          
      write(6,168)                                                      
      write(32,168)                                                     
  168 format(10x,'too large cfl number')                                
      go to 167                                                         
  266 continue                                                          
      write(6,268)                                                      
      write(32,268)                                                     
  268 format(10x,'velocities diverged')                                 
      go to 167                                                         
  169 continue                                                          
      write(6,178) dmax                                                 
      write(32,178) dmax                                                
  178 format(10x,'too large local residue for mass conservation : '     
     1       ,e12.5)                                                    
  167 continue                                                          
      return                                                            
      end                                                               
c                                                                       
c  **************  subrout totvel                                       
c   evaluates the integral of the velocity components
c   to check conservation properties
c                                                                       
      subroutine totvel
      include 'param.f'
      vit(1)=0.
      vit(2)=0.
      vit(3)=0.
      vit(4)=0.
      volto=1./(pi*alx3d)
      avgn=1./(float(n1m*n3m))
      do jc=1,n2m
      jp=jc+1
      do kc=1,n3m
      do ic=1,n1m
      vit(1)=q1(ic,jc,kc)/rm(jc)*volz(jc)+vit(1)
      vit(2)=(q2(ic,jp,kc)+q2(ic,jc,kc))/rm(jc)*volz(jc)+vit(2)
      vit(3)=q3(ic,jc,kc)*volz(jc)+vit(3)
      vit(4)=pr(ic,jc,kc)*volz(jc)+vit(4)
      enddo
      enddo
      enddo
      vit(1)=vit(1)*volto
      vit(2)=vit(2)*volto
      vit(3)=vit(3)*volto
      vit(4)=vit(4)*volto
      return
      end
c                                                                       
c  **************  subrout tschem                                       
c                                                                       
      subroutine tschem(ntime,time)
      include 'param.f'
      do 2000 ns=1,nsst                                                 
c
c
      al=alm(ns)
      ga=gam(ns)                                                        
      ro=rom(ns)                                                        
c                                                                       
c   time integration : implicit viscous, 3rd order rk (adams bashfort)  
c
      do  kc=1,n3m
      do  jc=1,n2m
      do  ic=1,n1m
      dq(ic,jc,kc)=0.
      qcap(ic,jc,kc)=0.
      dph(ic,jc,kc)=0.
         iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2m
      rhs(iadd)=0.
      enddo
      enddo
      enddo
c
c    non linear terms for momentum equation
c         uu usual term
c   from uu1   dq
      call hdnl1
c   from uu2   dph
      call hdnl2
c   from uu3   qcap
      if(n3m.gt.1) then
      call hdnl3
                   endif
      if(ren.gt..1e07) then
c
c   inviscid equations for energy conservation check
c
      call invinv(al,ga,ro,ns)                            
                       else
c     solve the dq1hat=q1hat-q1(n) momentum equation                    
      call invtu1(al,ga,ro,ns)                            
c     solve the dq2hat=q2hat-q2(n) momentum equation                    
      call invtu2(al,ga,ro,ns)                            
c     solve the dq3hat=q3hat-q3(n) momentum equation                    
      if(n3m.gt.1) then
      call invtu3(al,ga,ro,ns)                            
                   endif
                       endif
c
c     calculation of divg(dqhat)                                        
c***************************************** check 
      call divgu(al)                                         
c
c     calculation of the pressure dph by fft in two directions and      
c      tridiagonal in vertical                                          
      if(n3m.gt.1) then
c
c   3D simulations
c
      call phcalc
                   else
c
c   2D simulations in  r-theta planes
c
      call phcaij
                   endif
c     calculation of solenoidal vel field                               
      call updvp(al) 
c                                                                       
c
      call prcalc(al)                                                
c***************************************** check divg
 2000 continue                                                          
      return                                                            
      end                                                               
c************************************************************************
c                                                                       *
      subroutine openfi                                                 *
c                                                                       *
c************************************************************************
      include 'param.f'
c
c   open a bunch of files to write the outputs
c
      open(46,file='nfpipe')
      read(46,'(a)')filcnw
      read(46,'(a)')filcnr
      read(46,'(a)')filth
      read(46,'(a)')filvm
      read(46,'(a)')filpo
      read(46,'(a)')filen
      read(46,'(a)')filet
      read(46,'(a)')filer
      read(46,'(a)')filez
      read(46,'(a)')filed
      read(46,'(a)')filev
      open(32,file=filth)
      open(33,file=filvm)
      open(34,file=filpo)
      open(39,file=filen)
      open(40,file=filet)
      open(41,file=filer)
      open(42,file=filez)
      open(49,file=filed)
      open(50,file='piqm.out')
      open(59,file='piav.out')
      rewind 12
      rewind 33
      rewind 34
      rewind 32
      rewind 39
      rewind 40
      rewind 41
      rewind 42
      rewind 49
      rewind 48
      if(iwlop.gt.0) then
c
c   open the files where time evolution of velocity
c   vorticity and pressure are written
c
      numf=90
      kpr=n3m/12
      kpi=n3m/12+1
      ipr=n1m/4
      do jp=1,njprs
      j=npjp(jp)
      write(6,*) jp,j
c     do k=kpi,n3m,kpr
c     write(ipfk,82)k
      write(ipfj,82)j
   82 format(i3.3)
c     namfi3='fk'//ipfk//'j'//ipfj//'.dat'
      namfi3='fj'//ipfj//'.dat'
      open(numf,file=namfi3,form='unformatted')
      numf=numf+1
c     enddo
      enddo
      write(6,*)numf
                      endif
      return
      end   
