c************************************************************************
c                                                                       *
c                                                                       *
c     this cod is mad for simulating three-dimensional flows in polar *
c     ( cilindrical ) coordinates.                                      *
c     boundary condition can be free-slip on the boundaries in
c     the radial directions if inslws,inslwn are set = 0
c     no-slip if inslws,inslwn are set = 1
c     the boundary in the radial direction is free-slip if inslwr.eq.0
c     no-slip if inslwr.eq.1
c                                                                       *
c      navier-stokes equations are solved by a fractional step method   *
c     ( Kim and Moin ) with the pressure in the first step.             *
c                                                                       *
c     The time advancement of the solution is obtained by a             *
c     Runge-Kutta 3rd order low storage scheme (Wray) or a 2nd          *
c     order Adams-Bashfort scheme.                                      *
c                                                                       *
c     The Poisson  equation for the pressure is solved directly         *      
c     introducing FFT in the azimutal direction and the FISHPACK
c     package to solve the pentadioganal matrix in the r-z directions
c                                                                       *
c     All variables are calculated in a staggered grid:                 *
c                                                                       *
c        Flux variables are introduced to avoid  *
c        the singularity for r=0                                        *
c        in this case differently than in the code in the directory
c        CODESCH10/PIP    q1 = v_theta is used
c        the other variables are  q2= r*vr, q3= vz                                    *      
c                                                                       *
c        qcap :divergence of the  non free divergent velocity field     *
c                                                                       *
c        non linear terms:                                              *
c                                                                       *
c        ru1, ru2, ru3 : old step                                       *
c        h1,h2,h3 : new step                                            *
c        pr : pressure                                                  *
c        dph : pressure correction                                      *
c       pressure solver is in the file tr3dnuphnu.f                    *
c       non linear terms are calculated in hdnl routines in tr3dnuhn.f
c       the invertions of momentum equations in invtr in tr3dnutn.f
c
c       to run the different case the data files in CODESCH11/DAT
c       must be copied in the file tr3dfs.d
c                                                                       *
c                                                                       *
c************************************************************************
      program main                                                      
      include 'param.f'
      open(15,file='tr3dfs.d')
      read(15,301) dummy                                                
      read(15,*) n1,n2,n3,nsst,nwrit,nread
      read(15,301) dummy                                                
      read(15,*) n1p,n2p,n3p
      read(15,301) dummy                                                
      read(15,*) ntst,nprint,npin,nprde
      read(15,301) dummy                                                
      read(15,*) alx3p,istr3,str3,rmed31,etdp3,strb3
      read(15,301) dummy                                                
      read(15,*) re,vper,dt,pran,pscwal,irid
      read(15,301) dummy
      read(15,*) inslws,inslwn,inslwr,lamb
      read(15,301) dummy                                                
      read(15,*) cflc,icfl,tpin,tprin,tfin,tchpr,tprich
      read(15,301) dummy                                                
      read(15,*) ros                                                    
c$$$$$ parameters for non uniform grid (r distribution) $$$$$$$$$$$$$$$$$
      read(15,301) dummy                                                
      read(15,*)strr,rext,rint,rmed1
      read(15,301) dummy                                                
      read(15,*)istr,rmed,etdp,strb,n2t,n2v
c$$$$$$$parameters for tripol init. conditions $$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      read(15,301) dummy
      read(15,*) amp,alpha,crad,amp3,itrip
      read(15,301) dummy
      read(15,*) nwa,nwa3,iran1,iran3,sigr,rap,rapsc
c$$$$$$$parameters for vortex ring            
      read(15,301) dummy
      read(15,*) vmx,sig,yc2mo,yc3mo,iring
      read(15,301) dummy
      read(15,*) vper,epsil

301   format(a4)                                                        
      pi=2.*asin(1.)                                                    
      alx3d=2.*pi*alx3p
      n1m=n1-1                                                          
      n2m=n2-1                                                          
      n3m=n3-1                                                          
      n3mh=n3m/2+1                                                      
      if(n1m.eq.1) then
      iaxsy = 0
       else
      iaxsy=1
       endif
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
      write(6,*)'*               Tripol 3d      WITH ROTATION         *'
      write(32,*)'*              Tripol 3d      WITH ROTATION         *'
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
      write(6,112)alx3d,str3
      write(32,112)alx3d,str3
  112 format(10x,'domain dimensions L_z =',f8.5
     1          ,'     with str3=',f8.5)                
      write(6,*)'  '
      write(32,*) '  '
      write(6,*)'  '
      write(32,*) '  '
c                                                                       
      call solve                                                        
                                                                        
      stop                                                              
      end                                                               
c************************************************************************
c************************************************************************
c                                                                       *
c                                                                       *
c    code for computation of three-dimensional incompressible flows     *
c    in polar cylindrical coordinates.                                  * 
c                                                                       * 
c                                                                       * 
c    the geometry  is given in the subroutine cordi                     *
c                                                                       * 
c    the equations for q(i) i=1,2,3                                     *
c     q1=v(theta)    q2=v(r)*r     q3=v(zeta)                           *
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
c    method using fft routines and the fishpack package                  
c    no boundary conditions are necess. for the poisson equation.       *
c               a*ph=div(q)                                             *
c                                                                       *
c    all the second derivatives are discretized by a second order conse *
c    rvative scheme.                                                     *
c                                                                       *
c                                                                       *
c************************************************************************
c************************************************************************
      subroutine solve                                                  
      include 'param.f'
c
c*******************************************************
c                                                                       
c     grid definition, indices and mesh size calculation                          
c                                                                       
      call meshes
      call indic                                                        
      call cordin                                                       
c
c*******************************************************
c
c     print the grid
c
       ioldf=0
      call pricor
c
c*******************************************************
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
      ntst=10000
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
      call fftqua    
      call coetar
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
            if(itrip.eq.1) then
            write(6,*)'    initial vorticity of an isolated or Gaussian' 
     1               ,'vortex axial vorticity'
            call intrip
            write(6,*) ' initial velocity field generated'
            iring=0
                           endif
            if(iring.eq.1) then
            write(6,*)'    initial vorticity of a vortex ring ' 
     1               ,' azimuthal vorticity'
            call inring
                           endif
      else                                                              
c                                                                       
c   read  the initial fields   to continue a previous simulation   
c                                                                       
          write(6,*)' nread=1 ==>  the calculation is starting from  '      
          write(6,*)' read initial conditions                   '      
          write(32,*)' nread=1 ==>  the calculation is starting from  '      
          write(32,*)' read initial conditions                   '      
          write(6,*) '  '
          write(32,*) '  '
          call inirea(ntii,time,ntt,ncount,nap)           
          call divgck(dmax,dtot)                                 
          write(6,900)dmax,dtot                 
          write(32,900)dmax,dtot                
  900   format(3x,'maxima local and global divergence of the read field',
     1     /,' dmax = ',e11.4,' dtot = ',e11.4)
      endif                                                             
c
c     print other informations
c
      ntstf=ntst                                                   
      write(6,*) '  '
      write(32,*) '  '
      write(6,711)nprint,ntii,ntstf,npin,ntt
      write(32,711)nprint,ntii,ntstf,npin,ntt
  711 format(3x,'check in conditions :',/,
     1'  nprint = ',i5,'  ntii = ',i5,'  ntstf = ',i5,'  npin = ',i5
     1,'   ntt=',i5)            
      write(6,*) '  '
      write(32,*) '  '
  159 format(1x,i4,2x,e10.4,3e10.3,3(1x,e10.4,1x,i3,1x,i3),e10.3)
c
c*******************************************************
c
      write(6,*)'        time        vmath',
     1'         vmara        vmaax         omath',
     1'         omara        omaax      omiax',
     1'         cflm  '
      write(32,*)'        time        vmath',
     1'         vmara        vmaax         omath',
     1'         omara        omaax      omiax',
     1'         cflm  '
           call cfl(cflm)
           call outh(time,nav,ntime,cflm,nvv,navbu)            
            call outpf(time,nav)
            if(nprde.eq.1)       then
          write(6,*)'  nprde=',nprde,'   written init field'
             call contwr(ntime,time,ntt,nap)
                                 endif
c
c      write the time history
c
        ntii=ntii+1
c
c                                                                       
c  ********* starts the time dependent calculation ***                  
c                                                                       
      call phini    
      timew=time
      do 350 ntime=ntii,ntstf                                           
      if(icfl.eq.1) then
c
c   courant number calculation 
c
      call cfl(cflmm)
      dt=cflc/cflmm
      cflm=cflc
      if(dt.gt.dtl) then
      dt=dtl
      beta=dt/re*0.5                                                    
      cflm=cflmm*dt
                    endif
      timew=time
      if(dt.lt..1e-05)      then
      write(6,368)dt
  368 format(3x,'end calcul for dt=',e12.4)
      stop
                            endif
                    endif
      time=time+dt                                                      
c
c    time advancement
c
      call tschem(ntime,time)
c                                                                       
      if(icfl.eq.1) then
      if(amod(time,tpin).lt.dt) go to 306
      go to 305
                     else
      if(mod(ntime,npin).eq.0) then
      call cfl(cflmm)
      cflm=cflmm*dt
           go to 306                                
                               endif
                     endif
      go to 305                                                         
  306 continue                                                          
c*******************************************************
c
c     print some quantities  max(vel), max(vor) 
c
c*******************************************************
           call outh(time,nav,ntime,cflm,nvv,navbu)            
        ntt=ntt+1
c                                                                       
  305 continue                                                          
c*******************************************************
c                                                                       
c     write the flow field                                              
c                                                                       
c*******************************************************
      if(icfl.eq.1) then
      if(time.lt.tchpr)  then
      if(amod(time,tprin).lt.dt) go to 301
            go to 300
                         else
      if(amod(time,tprich).lt.dt) go to 301
            go to 300
                         endif
                     else
      if(mod(ntime,nprint).eq.0) go to 301
            go to 300
                     endif
  301 continue
            call outpf(time,ntime)
            if(nwrit.eq.1)       then
             call contwr(ntime,time,ntt,nap)
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
************************************************************************
c
c           SUBROUTINE  TSCHEM
c
c   This subroutine manages the whole integration scheme.
c   The following equations are solved:
c  
c    ~~     n
c   Q  -  Q                n         n       n-1   alp       2  ~~   n
c  --------- = -alp*grad (P ) + gam*H + rho*H   + ----- nabla ( Q + Q )
c    d t                                          2 Re
c
c          i                           i               i
c   where H  are the nonlinear terms, P  the pressure Q  the velocities
c       ~~
c   and Q  the provisional non solenoidal velocity field.
c   The superscripts (~~, n, n-1) indicate the time step level.
c                        n
c   The nonlinear terms H  are computed in the routines HDNL*, while
c   in the routines INVTR* are computed the remaining terms, updated
c   the non linear terms and inverted the equation to find the provisional
c   field at the new time step.
c       ~~
c   The Q  velocity field is projected onto a solenoidal field by a
c   scalar Phi computed through the equation
c
c                         2            1          ~~
c                    nabla (Phi ) =  ------ div ( Q  )
c                                    alp dt
c
c   The right hand side of this equation is computed in the routine
c   DIVG, while the equation is solved in PHCALC.
c
c   In the routine UPDVP the solenoidal velocity field at the new time
c   step is then computed through
c
c                n+1  ~~
c               Q   = Q  - alt*dt grad (Phi)
c
c   Finally in the routine PRCALC is updated the pressure field
c
c                n+1   n        alp dt      2
c               P   = P + Phi - ------ nabla (Phi)
c                                2 Re
c   When the scalar field is computed (density, concentration,
c   temperature) the routines HDNLps and INVTRps are used. The same
c   strategy at the velocity field is used, except that the scalar
c   field does not need any correction.
c
c   All variables are located on a staggered grid with the velocities
c   on the faces of the computational cell and all the scalars at the
c   centre. This is important when terms belonging to different equations
c   are avaluated.
c
c   Further details of the scheme can be found in the paper
c   "A finite-difference scheme for three-dimensional incompressible
c    flows in cylindrical coordinates" by R. Verzicco and P. Orlandi
c    J. of Comp. Phys. 1996.
c
c   REMEMBER that in this code q1(=v_theta
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
c    boundary conditions for passive scalar
c
      call pscbo
c                                                                       
c   time integration : implicit viscous, 3rd order rk (adams bashfort)  
c
c     solve the scalar   equation
c
      call hdnlps
      call invtrps(al,ga,ro,ns)
c

c
c
      do  kc=1,n3m
      do  jc=1,n2m
      do  ic=1,n1m
      dq(ic,jc,kc)=0.
      qcap(ic,jc,kc)=0.
      dph(ic,jc,kc)=0.
      enddo
      enddo
      enddo
      if(n1m.gt.1) then
      call hdnl1
                   endif
      call hdnl2
      call hdnl3
c     solve the dq1hat=q1hat-q1(n) momentum equation                    
      if(n1m.gt.1) then
      call invtr1(al,ga,ro,ns)                            
                   endif
c     solve the dq2hat=q2hat-q2(n) momentum equation                    
      call invtr2(al,ga,ro,ns)                            
                                                                        
c     solve the dq3hat=q3hat-q3(n) momentum equation                    
      call invtr3(al,ga,ro,ns)                            
c
c     calculation of divg(dqhat)                                        
c
      call divg(al)                                         
c
c     calculation of the pressure dph by fft and fishpack 
c      
      call phcalc
c
c     calculation of solenoidal vel field                               
c
      call updvp(al) 
c
c     calculation of real pressure                                      
c
      call prcalc(al)                                                
 2000 continue                                                          
      return                                                            
      end                                                               

c************************************************************************
c                                                                       *
      subroutine openfi                                                 *
c                                                                       *
c************************************************************************
      include 'param.f'
      open(46,file='nftr3dcy')
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
      return
      end

