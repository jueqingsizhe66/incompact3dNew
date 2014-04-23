c************************************************************************
c                                                                       *
c                                                                       *
c     this cod is mad for simulating three-dimensional flows in polar *
c     ( cilindrical ) coordinates.                                      *
c     boundary condition are slip-walls in the radial direction and     *
c     periodic in the vertical direction.                               *
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
c                                                                       *
c                                                                       *
c                                                                       *
c     All variables are calculated in a staggered grid:                 *
c                                                                       *
c        Instead of velocities, flux variables are introduced to avoid  *
c        the singularity for r=0                                        *
c        q1=vtheta*r, q2= r*vr, q3= vz                                    *      
c                                                                       *
c        dq1,dq2,dq3 : velocity correction                              *
c                                                                       *
c        qcap :divergence of the  non free divergent velocity field     *
c                                                                       *
c        non linear terms:                                              *
c                                                                       *
c        ru1, ru2, ru3 : old step                                       *
c        h1,h2,h3 : new step                                            *
c        pr : pressure                                                  *
c        dph : pressure correction                                      *
c       pressure solver is in the package phr3dn.f                      *
c       non linear terms are calculated in hdnl routines                *
c       the invertions of momentum equations is performed in invtr      *
c       routines                                                        *
c                                                                       *
c                                                                       *
c************************************************************************
      program main                                                      
      include 'param.f'
      common/timavg/timav
      common/timini/timei
      common/npjet/n2t
      common/nbal/nba,ibudg
      common/strpar/str2
      open(15,file='chabudg.d')
      read(15,*) n1,n2,n3,nsst
      read(15,*) nwrit,nread,iav,iprfi
      read(15,*) alx3d,alx1d,str2
      read(15,*) ren,vper
      read(15,*) dt,ntst,nprint,nprde,npouth
      read(15,*) nstop
      read(15,*) icfl,cflc,tpin,tprin,tfin                      
      read(15,*) jri,jrf,djr,irejr,iruuca
      read(15,*) timei
      read(15,*) islv1s,islv1n,islv3s,islv3n
      read(15,*) tosc,uosc
      read(15,*) flowq2,tau2,tim0sl
      read(15,*) y1gsd,y1ssd,y3gsd,y3ssd
      read(15,*) y1disd,y3disd 
      read(15,*) ifield
      read(15,*) itot,icorspe,icorr,timav,ibudg
301   format(a4)                                                        
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
c
      write(6,112)alx1d,alx3d
  112 format(10x,'chann. dimens ly=2',3x,'lx=',f4.2,'*pi'
     1        ,3x,'lz=',f4.2,'*pi')
      write(6,201)vper
      write(6,200)
  200 format(10x,'3d channel periodic in x3 and x1 ')
  201 format(3x,'programma con init. random perturbation'
     1   ,5x,'vper',e10.3)
      write(6,202) tosc,uosc
  202 format(10x,'parete inferiore, tosc=',e11.4,'uosc=',e11.4)
c                                                                       
      call openfi
c                                                                       
c     assign coefficients for time marching schemes                     
c
      write(6,*)'*************************************************'
      write(32,*)'*************************************************'
      write(6,*)'*                                                *'
      write(32,*) '*                                              *'
      write(6,*)'*      CHANNEL WITH WALL TRANSPIRATION           *'
      write(32,*)'*      CHANNEL WITH WALL TRANSPIRATION           *'
      write(6,*)'*      staistics from fields                     *'
      write(32,*)'*     staistics from fields                     *'
      write(6,*)'*************************************************'
      write(32,*)'*************************************************'
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
c    the geometry  is given in the subroutine cordi                     *
c                                                                       * 
c    the equations for q(i) i=1,2,3                                     *
c     q1=v(theta)*r  q2=v(r)*r     q3=v(zeta)                           *
c    are discretized by a finite difference scheme.                     *
c                                                                       *
c                                                                       *
c                                                                       *
c                                                                       *
c************************************************************************
c************************************************************************
      subroutine solve                                                  
      include 'param.f'
      common/tima3/tiax3d
      common/timavg/timav
      common/timini/timei
      common/nbal/nba,ibudg
      common/iolfil/iovco,iohel
      common/jspc/jspe,ifftco
      character*3 njpse
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
      write(6,755) ren,dt
      write(32,755) ren,dt
  755 format(3x,' Parameters of the flow: ',/,
     1 ' Reynolds number = ',e10.3,3x,'dt=',e10.3) 
      write(6,*) '  '
      write(32,*) '  '
      re=ren
      time=0.                                                           
      ntii=0                                                            
      beta=dt/re*0.5                                                    
            nap=0                                                           
            nav=0                                                           
            nvv=0                                                           
            nba=0
c
c  evaluation of metric quantities for the inversion
c
      call coetar
      ifftco=0
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
      write(32,*) '  '
  159 format(1x,i4,2x,e10.4,3e10.3,3(1x,e10.4,1x,i3,1x,i3),e10.3)
c
c*******************************************************
c
      write(6,*)'        time        pscto',
     1'         v1tot        v3tot         enert',
     1'         utat        utap        ttat',
     1'       ttav  cflm  '
      write(32,*)'       time        v1tot',
     1'         v2tot        v3tot         prto ',
     1'         pscto        enav         utat  ',
     1'         dp3ns qsour   '
c
c      write the time history
c
c
c*******************************************************
c                                                                       
c  ********* starts the time dependent calculation ***                  
c                                                                       
      ncount=0
      ntii=1
       time=timei
      write(6,*)'ntii,ntst,dt   ',ntii,ntst,dt
      do 350 ntime=ntii,ntst                                           
c     the calculation stops if the velocities are diverging for numerical
c     stability conditions (courant number restrictions)
      call tschem(ntime,time,ncount)
       ncount=ncount+1
c
        ntt=ntt+1
c                                                                       
c*******************************************************
c                                                                       
c     write the flow field                                              
c                                                                       
c*******************************************************
           call outh(time,nav,ntime,cflm,nvv,navbu)            
         if(ntime.eq.ntst) go to 351
        time=time+dt
  350 continue                                                          
  351 continue                                                          
      return                                                            
      end                                                               
c************************************************************************
c                                                                       *
c ****************************** subrout coetar  **********************  *
c                                                                       *
c    this subroutine calculates the coefficients for the              *
c    integration in the radial direction with non-uniform coor. trasf.  *
c   NEW VERSION FOR q\theta=v\theta*r
c    Bound. Cond. 1st oreder accurate
c                                                                       *
c************************************************************************
      subroutine coetar
      include 'param.f'
      common/cor3j/ap3j(m2),ac3j(m2),am3j(m2)
c
c  *******   coefficients in several eq. funct of j
c
      do jc=1,n2m
      udx1q(jc)=dx1q
      volz(jc)=caj(jc)/(dx1*dx2*dx3)
      enddo
c
c  ***********  coefficients for q3   inner points
c
      do jc=2,n2m-1
       jp=jc+1
       jm=jc-1
       ucaj=1./caj(jc)
       ap3j(jc)=1./cac(jp)*ucaj
       ac3j(jc)=(1./cac(jp)+1./cac(jc))*ucaj
       am3j(jc)=1./cac(jc)*ucaj
      end do
c
c   set up the cefficients apj1, acj1, amj1 at the boundaries
c
c  jc=1
      ucaj=4./(2.*cac(2)+caj(1))
      ap3j(1)=1./cac(2)*ucaj
      ac3j(1)=(1./cac(2)+2./cac(1))*ucaj
      am3j(1)=0.
c  jc=n2m
      ucaj=4./(2.*cac(n2m)+caj(n2m))
      ap3j(n2m)=0.
      ac3j(n2m)=(2./cac(n2)+1./cac(n2m))*ucaj
      am3j(n2m)=1./cac(n2m)*ucaj
      return
      end
c************************************************************************
c                                                                       *
c ****************************** subrout prgqso  ********************** *
c  this subroutine performs the calculation of pressure gradient and    *
c  heat flux.                                                           *
c                                                                       *
c************************************************************************
      subroutine prgqso              
      include 'param.f'
      common/cor3j/ap3j(m2),ac3j(m2),am3j(m2)
      pi=2.*asin(1.)
      alre=al/ren
c
c
c    Pressure gradient
c
      alre=al/ren
      s3tot=0.
      do 18 kc=1,n3m
      km=kmv(kc)
      kp=kpv(kc)
      do 18 jc=1,n2m
      jmm=jmv(jc)
      jpp=jpv(jc)
      do 18 ic=1,n1m
      im=imv(ic)
      ip=ipv(ic)
c
c   11 second derivatives of q3
c
      dq31=(q3(ip,jc,kc)-2.*q3(ic,jc,kc)+q3(im,jc,kc))*dx1q
c
c   22 second derivatives of q3
c
      dq32=(ap3j(jc)*q3(ic,jpp,kc)
     1     -ac3j(jc)*q3(ic,jc,kc)
     1     +am3j(jc)*q3(ic,jmm,kc))*dx2q
c
c   33 second derivatives of q3
c
      dq33=(q3(ic,jc,kp)-2.*q3(ic,jc,kc)+q3(ic,jc,km))*dx3q
      dcq3=dq31+dq33+dq32
      s3tot=s3tot+dcq3*caj(jc)/ren
   18 continue
      dp3ns=s3tot/(2.*n1m*n2m*n3m)
      return
      end
c                                                                       
c  **************  subrout tschem                                       
c                                                                       
      subroutine tschem(ntime,time,ncount)
      include 'param.f'
      dimension voz1(m1),vot1(m1)
      common/timavg/timav
      common/timini/timei
      common/nbal/nba,ibudg
      common/jspc/jspe,ifftco
      common/wavin/dlx1,dlx3,dkk1,dkk3
c
c   fine check spectra
c
          call inirea(ntime,time,ntt,ncount,nap)           
          call divgck(dmax,dtot)                                 
          write(6,900)dmax,dtot                 ,ntime
          write(32,900)dmax,dtot                ,ntime
  900   format(3x,'maxima local and global divergence of the read field',
     1     /,' dmax = ',e11.4,' dtot = ',e11.4,'ntime=',i5)
          call prgqso
           call velc
           call vorc
      n3mh=n3m/2+1
      n1mh=n1m/2+1
      return                                                            
      end                                                               
