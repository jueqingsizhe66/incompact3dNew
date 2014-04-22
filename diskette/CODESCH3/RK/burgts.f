c************************************************************************
c                                                                       *
c                                                                       *
c     this cod is mad for simulating burgers staggered                *
c     The time advancement of the solution is obtained by a             *
c     Runge-Kutta 3rd order low storage scheme (Wray) or a 2nd          *
c     order Adams-Bashfort scheme.                                      *
c     This code permits to understand how to solve the Navier-Stokes    *
c     Equations in 2D and 3D. The structure of the code is similar      *
c     the pressure solver do not appear.
c                                                                       *
c************************************************************************
      program main                                                      
      include 'param.f'
      character*4 nprep
      character*4 npn1m
      open(15,file='burgts.d')
      read(15,301) dummy                                                
      read(15,*) n1,nsst
      read(15,301) dummy                                                
      read(15,*) ntst,nprint,npin,npouth
      read(15,301) dummy                                                
      read(15,*) vper,dt
      read(15,301) dummy                                                
      read(15,*) cflc,icfl,ptpin,ptprin,ptfin,tpouth
      read(15,301) dummy                                                
      read(15,*) pnu
c
c$$$$$ parameters for non uniform grid  $$$$$$$$$$$$$$$$$
c
      read(15,301) dummy                                                
      read(15,*)strr,rext
      read(15,301) dummy                                                
      read(15,*)istr
301   format(a4)                                                        
      pi=2.*asin(1.)                                                    
      n1m=n1-1                                                          
      n1mh=n1m/2
      n1mp=n1m/2+1
      re=pi/(pnu)
      dt=dt/pi
      tpin=ptpin/pi
      tprin=ptprin/pi
      tfin=ptfin/pi
c
      if(nsst.gt.1) then                                                
      gam(1)=8./15.                                                      
      gam(2)=5./12.                                                      
      gam(3)=3./4.                                                       
      rom(1)=0.                                                          
      rom(2)=-17./60.                                                    
      rom(3)=-5./12.                                                     
      write(6,*)'*****************************************************'
      write(31,*)'*****************************************************'
      write(6,*)'*                                                    *'
      write(31,*) '*                                                  *'
      write(6,*)'* BURGERS STAGGERED                                  *'
      write(31,*)'* BURGERS STAGGERED                                 *'
      write(6,*)'*                                                    *'
      write(31,*) '*                                                  *'
      write(6,*)'*****************************************************'
      write(31,*)'*****************************************************'
      write(6,*)'  '
      write(31,*) '  '
      write(6,100) (gam(n),n=1,nsst),(rom(n),n=1,nsst)                    
      write(31,100) (gam(n),n=1,nsst),(rom(n),n=1,nsst)                    
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
      write(31,110) gam(1),rom(1)                                          
  110 format(3x,'the time advancement of the solution is obtained by',/,
     1  'an Adams-Bashfort scheme. Coefficients are:  ',/, 
     1       'gam=  ',f8.3,6x,'ro=  ',f8.3)                                
      endif                                                             
c                                                                       
c                                                                       
      do 10 ns=1,nsst
      alm(ns)=(gam(ns)+rom(ns))
   10 continue
c                                                                       
      irepn=nint(1./pnu)
      write(6,*)'pnu=',pnu,'    irepn',irepn
      write(nprep,82) irepn
      write(npn1m,83) n1m
   82 format(i4.4)
   83 format(i4.4)
      open(31,file='burgth'//nprep//'.'//npn1m)
      open(32,file='burg'//nprep//'.'//npn1m)
      call solve                                                        
                                                                        
      stop                                                              
      end                                                               
c************************************************************************
c************************************************************************
c                                                                       *
c                                                                       *
c                                                                       *
c                                                                       *
c************************************************************************
c************************************************************************
      subroutine solve                                                  
      include 'param.f'
c                                                                       
c     grid definition, indices and mesh size calculation                          
c                                                                       
      call meshes
      call indic                                                        

      call cordin                                                       
c
c*******************************************************
c
c
c*******************************************************
c
c     print some informations on the run
c
      write(6,*) '  '
      write(6,754)n1
  754 format(10x,'number of grid points :'/                             
     1      5x,'n1=',i4)                       
      write(6,*) '  '
      write(6,755) re,pnu
      write(31,755) re,pnu
  755 format(3x,' Parameters of the flow: ',/,
     1 ' Reynolds number = ',e10.3,3x,3x,'pi*nu=',e12.3) 
      if(icfl.eq.1) then
      ntst=100000
      dtl=dt
      write(6,*)' calculation at dt variable by fixing cfl=',cflc
      write(31,*)' calculation at dt variable by fixing cfl=',cflc
                     else
      write(6,765) dt
      write(31,765) dt
  765 format(' time step dt= ',e10.3)
                     endif
      write(6,*) '  '
      write(31,*) '  '
      ren=re                                                            
      time=0.                                                           
      ntii=1                                                            
      beta=dt/re*0.5                                                    
      devma=0.
            nav=0                                                           
c                                                                       
c    create the initial fields                                       
c                                                                       
      call initia
c
c     print other informations
c
      ntstf=ntst                                                   
      write(6,*) '  '
      write(31,*) '  '
      write(6,711)nprint,ntii,ntstf,npin
      write(31,711)nprint,ntii,ntstf,npin
  711 format(3x,'check in conditions :',/,
     1'  nprint = ',i5,'  ntii = ',i5,'  ntstf = ',i5,' npin = ',i5)
      write(6,*) '  '
      write(31,*) '  '
c
c                                                                       
c  ********* starts the time dependent calculation ***                  
c                                                                       
      call coetar
           time=0.
      duxc=(q1(n1mp+1)-q1(n1mp-1))/(yp1(n1mp+1)-yp1(n1mp-1))
           call outth(time,dt)
           call outp(time)
      ntt=0
      do 350 ntime=ntii,ntstf                                           
                                                                        
c     the calculation stops if the velocities are diverging for numerical
c     stability conditions (courant number restrictions)                
      timew=time                                                                  
      if(icfl.eq.1) then
c
c   Courant number calculation 
c
      call cflu(cflum,q1ma)
      cflmm=cflum
      cflm=cflmm*dt
      
      dt=cflc/cflmm
      cflm=cflc
      if(dt.gt.dtl) then
      dt=dtl
      beta=dt/re*0.5
      cflm=cflmm*dt
                    endif
      cflummm=dt*cflum
c     write(6,*) cflummm
      if(dt.lt..1e-05)      then
      write(6,368)ntime,dt
  368 format(3x,'after ntim=',i4,3x,'end calcul for dt=',e12.4)
      stop
                            endif
                    endif
      time=time+dt                                                      
      ntt=ntt+1
      call tschem(ntime,time)
      duxc=(q1(n1mp+1)-q1(n1mp-1))/(yp1(n1mp+1)-yp1(n1mp-1))
      if(abs(duxc).gt.devma) then
      devma=abs(duxc)
      tidem=time*pi
                            endif
c                                                                       
      if(icfl.eq.1) then
      if(amod(time,tpin).lt.dt) go to 306
      go to 305
                     else
      if(mod(ntime,npin).eq.0) then
      call cflu(cflum,q1ma)
      cflmm=cflum
      cflm=cflmm*dt
      if(q1ma.ge.1000.) go to 166
           go to 306
                               endif
                     endif
      go to 305
  306 continue
c*******************************************************
c
c     print some global quantities 
c
c*******************************************************
           call outth(time,cflm)
c                                                                       
  305 continue                                                          
c*******************************************************
c                                                                       
c     write the q1 field                                              
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
           call outp(time)
  300 continue
       if(time.gt.tfin.and.icfl.eq.1) then
      dymin=yp1(n1mp+1)
      write(31,367)time,ntt
      write(31,369)devma,tidem,dymin
      write(6,367)time,ntt
      write(6,369)devma,tidem,dymin
  367 format(3x,'end calcul for t=',e12.4,'  ntt=',i6)
  369 format(3x,'devma=',e15.6,3x,'at t=',e15.6,' dymin= ',e12.4)
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
      write(31,168)                                                     
  168 format(10x,'too large cfl number')                                
      go to 167                                                         
  167 continue                                                          
      dymin=yp1(n1mp+1)
      write(6,367)time,ntt
      write(6,369)devma,tidem,dymin
      write(31,367)time,ntt
      write(31,369)devma,tidem,dymin
      return                                                            
      end                                                               
c************************************************************************
c                                                                       *
c    *****   subro outp   ********************************              *
c                                                                       *
c************************************************************************
      subroutine outp(tim)
      include 'param.f'
      character*3 nptim
      character*4 npn1m
      enen=0.
      time=tim*pi
      ntim=nint(time*tpouth)
      write(npn1m,82) n1m
   82 format(i4.4)
      write(nptim,83) ntim
   83 format(i3.3)
      open(13,file='prof'//nptim//'.'//npn1m)
      do ic=1,n1
      write(13,133) yp1(ic),q1(ic)
  133 format(5e12.5)
      enddo
      return
      end
c************************************************************************
c                                                                       *
c    *****   subro outh   ********************************              *
c                                                                       *
c************************************************************************
      subroutine outth(tim,cflm)
      include 'param.f'
      call totvel
      time=tim*pi
      enen=0.
      q1ma=0.
      do ic=1,n1m
      ip=ipv(ic)
      q1p=(q1(ic)+q1(ip))*0.5
      enen=enen+q1p**2/dx1*g1m(ic)
      q1ma=max(abs(q1(ic)),q1ma)
      enddo
      write(31,133)time,q1ma,vit(1),enen,duxc,cflm
      write(6,133)time,q1ma,vit(1),enen,duxc,cflm
      dupr=abs(duxc*pnu)
      write(32,133)time,dupr
  133 format(8e12.5)
      return
      end
c                                                                       
c  **************  subrout totvel                                       
c                                                                       
      subroutine totvel
      include 'param.f'
      vit(1)=0.
      do ic=1,n1m
      ip=ipv(ic)
      vit(1)=(q1(ic)+q1(ip))*0.5+vit(1)
      enddo
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
c
c    non linear terms for momentum equation
c 
c   from uu1   dq
      call hdnuu1
      
c
c     solve the dq1hat=q1(n+1)-q1(n) momentum equation                    
c
      call invtu1(al,ga,ro)                            
      q1ma=0.
      dq1ma=0.
        do ic=1,n1
         q1ma=max(q1(ic),q1ma)
         dq1ma=max(rhs(ic),dq1ma)
         q1(ic)=q1(ic) + rhs(ic)
        end do
c***************************************** check divg
 2000 continue                                                          
      write(61,*)time,q1ma,dq1ma
      return                                                            
      end                                                               
c**********************************************************************ni*
c                                                                      *
c     this routine gives the initial condition                         *
c                                                                      * 
c***********************************************************************
      subroutine initia
      include 'param.f'
      write(6,*)  '  '
      write(31,*) '  '
      write(6,*) ' I N I T I A L   C O N D I T I O N '
      write(31,*) ' I N I T I A L   C O N D I T I O N '
      write(6,*)  '  '
      write(31,*) '  '
c*******************    q1 velocity          ************************
      do  i=1,n1
      q1(i)=-sin(pi*yp1(i))
      enddo
      return                                                            
      end                                                               
c***********************************************************************
c                                                                       *
c  ****************************** subrout indic **********************  *
c                                                                       *
c     in this subroutine the indices ip,im are evaluated in all points  *
c     in reality these are not used to account for periodicity or
c     boundary conditions as is done for the Navier-Stokes codes
c                                                                       *
c************************************************************************
      subroutine indic                                                  
      include 'param.f'
c                                                                       
c                                                                       
c   azimuthal periodic direction                                               
c                                                                       
      do 1 ic=1,n1m                                                     
      imv(ic)=ic-1                                                      
      ipv(ic)=ic+1                                                      
    1 continue                                                          
      return                                                            
      end                                                               
c************************************************************************
c                                                                       *
c  ****************************** subrout cflu  **********************   *
c                                                                       *
c************************************************************************
      subroutine cflu(cflm,q1ma)                                     
      include 'param.f'
c                                                                       
c     in this routine the COURANT number is calculated.
c     This parameter determines the stability condition
c     (CFL < 1 for Adams-Bashfort  CFL < 1.7 for third order Runge-Kutta)
c     for the calculation at dt = const.
c     If a calculation at variable dt is used, this parameter determines
c     the dt itself.
c
      cflm=0.                                                           
      q1ma=0.                                                           
c                                                                       
      do i=2,n1m                                                    
      ip=ipv(i)                                                         
      qcf=abs(q1(i)+q1(ip))*0.5*dx1/g1m(i)
      cflm=max(cflm,qcf)                                              
      q1ma=max(q1ma,q1(i))                                              
      enddo
      return                                                            
      end                                                               
c************************************************************************
c                                                                       *
c  ****************************** subrout meshes ********************** *
c                                                                       *
c************************************************************************
c
c     The calculation of the mesh is performed.
c     THe physical coordinate are evaluated in the routine CORDI
c
      subroutine meshes
      include 'param.f'
      dx1=1./float(n1m)
      write(6,100)dx1
 100  format(3x,'mesh size: d_theta= ',e9.4,' d_r= ',e9.4,' dz= ',e9.4)
      write(6,*)'  '
      write(31,*) '  '
      dx1=1./dx1                                                        
      dx1q=dx1*dx1                                                      
      return                                                            
      end                                                               
c************************************************************************
c                                                                       *
c   ************** subroutine cordin                                    *
c                                                                       *
c************************************************************************
      subroutine cordin                                                 
c
c     Physical coordinates are assigned
c
      include 'param.f'
c
      open(78,file='rg1c.out')
      open(79,file='rg1m.out')
      if (istr.lt.0) then
c
c   the physical coordinate has constant mesh
c

      do 111 i=1,n1
       x1=(i-1)/dx1
       yp1(i)=-rext+2.*rext*x1
  111 continue
       endif
      if (istr.eq.0) then
c
c   the physical coordinate is clustered near x=0
c
      tstr2=tanh(strr)
      do 112 i=1,n1
       x1=(i-1)/dx1
       if(x1.le.0.5) then
       arg=2.*strr*x1
       yp1(i)=(-1.+tanh(arg)/tstr2)*rext
                     else
       arg=2.*strr*(x1-1)
       yp1(i)=(+1.+tanh(arg)/tstr2)*rext
                     endif
  133 format(4e12.4)
  112 continue
       endif
      do 12 i=1,n1m                                                     
      ym1(i)=(yp1(i)+yp1(i+1))*0.5 
      g1m(i)=(yp1(i+1)-yp1(i))*dx1
   12 continue                                                          
      do 122 i=2,n1m
      g1c(i)=(yp1(i+1)-yp1(i-1))*dx1*0.5
122   continue
      g1c(1)=(ym1(1)-yp1(1))*dx1*2.
      g1c(n1)=(yp1(n1)-ym1(n1m))*dx1*2.
      do i=1,n1
       x1=(i-1)/dx1
      write(78,133)x1,yp1(i),g1c(i)
      enddo
      do i=1,n1m
       x1=(i-1+0.5)/dx1
      write(79,133)x1,ym1(i),g1m(i)
      enddo
      close(78)
      close(79)
      return                                                            
      end                                                               
c************************************************************************
c                                                                       *
c ****************************** subrout hdnuu1  **********************  *
c    the non-linear term  is calculated.                                *
c                                                                       *
c*************************************************************************
      subroutine hdnuu1
      include 'param.f'

c
c     compute the non-linear terms by centered difference.
c     H term for the q1 momentum equation at i+1/2
c
c
      do ic=2,n1m
      ip=ipv(ic)
      im=imv(ic)
      udx1c=dx1/g1c(ic)*0.25
c
c    q1 q1 term
c
c                   d (q_1)^2
c                  ---------
c                   d  x 
c
      htt=( (q1(ip)+q1(ic))*(q1(ip)+q1(ic))
     1     -(q1(im)+q1(ic))*(q1(im)+q1(ic))
     1    )*udx1c
      dq(ic)=-htt*0.5
      enddo
      return                                                            
      end                                                               
c************************************************************************
c                                                                       *
c  ************************ subrout coetar  **********************      *
c                                                                       *
c  this subroutine evaluates the coefficients of the tridiagonal matrix *
c  associated to the second derivative and to the Dirichlet b.c.        *
c                                                                       *
c************************************************************************
      subroutine coetar
      include 'param.f'
      common/cor1i/ap1i(m1),ac1i(m1),am1i(m1)
      do ic=2,n1m
      ip=ipv(ic)
      im=imv(ic)
      ucai=dx1q/g1c(ic)
      ap1i(ic)=1./g1m(ic)*ucai
      am1i(ic)=1./g1m(im)*ucai
      ac1i(ic)=-(1./g1m(im)+1./g1m(ic))*ucai
      enddo
      do ic=1,n1,n1m
      ap1i(ic)=0.
      am1i(ic)=0.
      ac1i(ic)=1.
      enddo
      return                                                            
      end                                                               
c************************************************************************
c                                                                       *
c  ************************ subrout invtu1  **********************      *
c                                                                       *
c  this subroutine performs the inversion of the q1 momentum equation   *
c                                                                       *
c************************************************************************
      subroutine invtu1(al,ga,ro)                      
      include 'param.f'
      common/cor1i/ap1i(m1),ac1i(m1),am1i(m1)
      alre=al/ren                                                       
c                                                                       
c  ********* compute the rhs of the factored equation                   
c  the rhs includes h(n), h(n-1)=ru
c  and the 11 second derivat. of q1(n)                            
c  everything at i+1/2                                          
c  dq=q(n+1)-q(n)                                                         
c                                                                       
      do ic=2,n1m                                                   
      ip=ipv(ic)                                                        
      im=imv(ic)                                                        
      d11q1=q1(ip)*ap1i(ic)
     1     +q1(ic)*ac1i(ic)
     1     +q1(im)*am1i(ic)
      rhs(ic)=(ga*dq(ic)+ro*ru1(ic)+alre*d11q1)*dt
      ru1(ic)=dq(ic)                                        
      enddo                                                             
      call solq1i(al) 
      return                                                            
      end                                                               
c************************************************************************
c                                                                       *
c ****************************** subrout solq1i  ********************** *
c                                                                       *
c  this subroutine performs the inversion of the q1 momentum equation   *
c  by a factored implicit scheme, only the derivatives 11,22,33 of q1   *
c  are treated implicitly                                               *
c                                                                       *
c************************************************************************
      subroutine solq1i(al)                                      
      include 'param.f'
      common/cor1i/ap1i(m1),ac1i(m1),am1i(m1)
      dimension api(m1),aci(m1),ami(m1)
c
c    Solve for the equation:
c
c
c      (       al * dt [ d     d    ] )
c      ( 1 -  -------- [ ---  ----  ] ) ( q_1) =  RHS
c      (       2 * Re  [ dx1   d x1 ] )
c                                                                       
c  ************ compute dq1  sweeping along the x1 direction            
c                                                                       
c
      betadx=beta*al
      do 21 ic=2,n1m
      api(ic)= -ap1i(ic)*betadx
      aci(ic)=1.-ac1i(ic)*betadx
      ami(ic)= -am1i(ic)*betadx
   21 continue
      do ic=1,n1,n1m
      api(ic)=0.
      ami(ic)=0.
      aci(ic)=1.
      rhs(ic)=0.
      enddo
      nx=n1
      call trib(ami,aci,api,rhs,nx)                                    
      return                                                            
      end                                                               
c
c
c  ****************************** subrout trib  **********************
c
      subroutine trib(a,b,c,f,n)
      include 'param.f'
      dimension c(m1),b(m1),a(m1),f(m1),gm(m1)
c
c   solution of a tridiagonal matrix
c

      bet=b(1)
      f(1)=f(1)/bet
      do j=2,n
      gm(j)=c(j-1)/bet
      bet=b(j)-a(j)*gm(j)
      f(j)=(f(j)-a(j)*f(j-1))/bet
      enddo
      do j=n-1,1,-1
      f(j)=f(j)-gm(j+1)*f(j+1)
      enddo
      return
      end
