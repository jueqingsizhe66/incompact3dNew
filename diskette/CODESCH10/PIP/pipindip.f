c  *************************************************************        
c     this routine evaluates the radial and azimuthal velocity 
c     distribution in the r-theta planes by the streamfunction
c     of a Lamb dipole evaluated by solving the nabla psi =-omega
c     with omega from the Lamb expression
c                                                                       
c  *************************************************************        
      subroutine indipo
      include 'param.f'
      call psini                                                        
      do 100 k=1,n3m                                                    
      call dipole
c     q3 comp.                                                          
c     q1 comp.                                                          
      q1max=0.                                                          
      do 810 i=1,n1m                                                    
      do 810 j=1,n2m                                                    
      q1(i,j,k)=-rm(j)*(psi(i,j+1)-psi(i,j))*dx2/g2rm(j)
      q1max=max(q1max,abs(q1(i,j,k)))                                 
  810 continue                                                          
c                                                                       
c     q2 comp.                                                          
      q2max=0.                                                          
      do 918 i=1,n1m                                                    
      q2(i,1,k)=0.                                                      
      q2(i,n2,k)=0.                                                     
  918 continue                                                          
      do 815 j=1,n2m                                                    
      do 815 i=1,n1m                                                    
      ip=ipv(i)                                                         
      q2(i,j,k)=(psi(ip,j)-psi(i,j))*dx1
      q2max=max(q2max,abs(q2(i,j,k)))                                 
  815 continue       
  100 continue                                                          
c                                                                       
c q3 is set zero                             
c                                                                       
      do 454 kc=1,n3m                                                   
      do 454 jc=1,n2m                                                   
      do 454 ic=1,n1m                                                 
      q3(ic,jc,kc)=0.                       
 454  continue                                                          
      call vmaxv(n1m,n2m,n3m)
      call divuck(qmax,qtot)
      write(6,705) qmax,qtot
      write(32,705) qmax,qtot
  705 format(3x,'from inqpr  qmax and qtot  =',2e11.4)                   
      write(6,700) (vmax(l),l=1,3)                                      
      write(32,700) (vmax(l),l=1,3)                                     
  700 format(1x,'vmx1=',e11.4,2x,'vmx2=',e11.4,2x,'vmx3=',e11.4)        
      return                                                            
      end                                                               
c  ****************************** subrout tribjj ********************** 
c  triagonal solver to invert the Laplacian psi                         
c
      subroutine tribjj(r,n,uu,m)                                  
      include 'param.f'
      parameter (m1mh=m1m/2+1)
      common/ctrdps/amps(m2),acps(m1mh,m2),apps(m2)                     
      dimension gammm(m1mh,m2),r(m1mh,m2)          
      dimension bett(m1mh),uu(m1mh,m2)                                    
      do 10 i=1,m                                                       
      bett(i)=acps(i,1)                                                     
      uu(i,1)=r(i,1)/bett(i)                                              
   10 continue                                                          
      do 11 j=2,n                                                       
      do 21 i=1,m                                                       
      gammm(i,j)=apps(j-1)/bett(i)                                            
      bett(i)=acps(i,j)-amps(j)*gammm(i,j)  
      uu(i,j)=(r(i,j)-amps(j)*uu(i,j-1))/bett(i)                              
   21 continue                                                          
   11 continue                                                          
      do 12 j=n-1,1,-1                                                  
      do 22 i=1,m                                                       
      uu(i,j)=uu(i,j)-gammm(i,j+1)*uu(i,j+1)                                 
   22 continue                                                          
   12 continue                                                          
      return                                                            
      end                                                               
c
c  **********************  inizialize field in r-theta section  *********    
c                                                                       
      subroutine dipole
      include 'param.f'
      dimension xx(m2)                                        
c
c     DIPOLE in the r-theta plane
c
c
      pi=2.*asin(1.)
      akmo=3.83
      thet0 = 0.
      ramo = 1.
      vsi=1.
      vomax=0.
      besj0=bessj0(akmo)
      do 917 i=1,n1m
      do 917 j=1,n2
      y1d=rc(j) * cos(thetac(i))-yc1mo
      y2d=rc(j) * sin(thetac(i))-yc2mo
      ramod=sqrt(y1d**2+y2d**2)
      radmod=ramod/ramo
      arbesj=radmod*akmo
      if(ramod.lt..1e-06)then
      vor(i,j)=0.
      else
      if(radmod.ge.1.) then
      vor(i,j)=0.
      else
      stho=y2d/ramod
      if(abs(stho).gt.1.) write(6,719) i,j,y2d,ramod
  719 format(3x,2i4,2x,'err in asin',2e12.4)
      astho=(abs(stho)-1.)
      if(astho.gt..0) then
      stho=stho-.1e-06*stho/abs(stho)
      endif
      thet=asin(stho)
      if(y2d.lt.0..and.y1d.gt.0.) then
      thet=pi*2.-abs(thet)
      endif
      if(y2d.lt.0..and.y1d.lt.0.) then
      thet=pi+abs(thet)
      endif
      if(y2d.ge.0..and.y1d.lt.0.) then
      thet=-abs(thet)+pi
      endif
      sth=sin(thet)*cos(thet0)-sin(thet0)*cos(thet)
      besj1=bessj1(arbesj)
c     write(6,*)i,j,vsi,velmo,sth,arbesj,besj1,besj0,akmo,ramo
      vom=2.*vsi*velmo*sth*besj1/besj0*akmo/ramo
      vor(i,j)=vom
      vorma=max(abs(vor(i,j)),vorma)
      endif
      endif
  917 continue

           vorin=0.                                                        
           vormax=0.                                                         
           vormin=0.                                                         
           do 42 j=1,n2                                                      
           do 42 i=1,n1m                                                     
           if(vor(i,j).gt.vormax) then                                
             vormax=vor(i,j)                                            
             jmax=j                                                            
             imax=i                                                            
           endif                               
             if(vor(i,j).lt.vormin) then                                
               vormin=vor(i,j)                                            
               jmin=j  
               imin=i 
             endif                               
        vorin=vorin+vor(i,j)*rc(j)/(dx1*dx2)                                
        psi(i,j)=0.
   42 continue                                                          
c                                                                       
       write(6,701) vormax,imax,jmax                                     
       write(32,701) vormax,imax,jmax                                    
 701  format(3x,'vor.  max=  ',e11.4,2x,'in i,j',2i4)                        
       write(6,702) vormin,imin,jmin                                     
       write(32,702) vormin,imin,jmin                                    
 702  format(3x,'vor.  min=  ',e11.4,2x,'in i,j',2i4)                        
       write(6,713) vorin                                        
       write(32,713)vorin                                       
 713  format(1x,'circolazione iniziale  =  ',e10.3)
c
c   calculation of the streamfunction
c
       call pscalc
c                                                                       
      vormax=0.                                                         
      psimax=0.                                                         
      do 101 j=2,n2m                                                    
      jm=jmv(j)
      do 101 i=1,n1m                                                    
      psimax=max(abs(psi(i,j)),psimax)                                
      im=imv(i)                                                         
      ip=ipv(i)                                                         
      vor(i,j)=(psi(ip,j)-2.*psi(i,j)+psi(im,j))*dx1q/rc(j)**2           
     1       +(  (psi(i,j+1)-psi(i,j))*rm(j)/(rc(j+1)-rc(j))
     1          -(psi(i,j)-psi(i,jm))*rm(jm)/(rc(j)-rc(jm)) )
     1         /(rc(j)*(rm(j)-rm(jm)))
      if(abs(vor(i,j)).gt.vormax) then                                
      vormax=abs(vor(i,j))                                            
      jmax=j                                                            
      imax=i                                                            
                                    endif                               
  101 continue                                                          
       write(6,721) psimax,vormax,imax,jmax                              
       write(32,721) psimax,vormax,imax,jmax                              
  721 format(1x,'psimax=',e11.4,2x,                                     
     1      'vormax=',e11.4,2x,'at i,j',2i4)                            
      return                                                            
      end                                                               
c                                                                       
c                                                                       
c   ********************* subr fftqps                                   
c  this subroutine perform the calculation of trigz for temperton fft   
c  used to invert the Laplacian psi
c                                                                       
      subroutine fftqps                                                 
      include 'param.f'
      parameter (m1mh=m1m/2+1)
      common/fftcom/ifxz(13),trigz(3*m1m/2+1)                           
      common/wavps/an1(m1mh),ank1(m1mh)                                 
      n1mh=n1m/2+1                                                      
      pi=2.*asin(1.)                                                    
c                                                                       
c     modified wave number definition                                   
c                                                                       
      do 16 i=1,n1mh                                                    
   16 an1(i)=(i-1)*2.*pi                                                
      call fftfax(n1m,ifxz,trigz)                                       
      return                                                            
      end                                                               
c                                                                       
c   ********************* subr pscalc                                   
c   this subroutine performs the solution of the poisson equation       
c   related to the stream function calculation for initial conditions   
c   real  fft for x1 direction                                          
c                                                                       
      subroutine pscalc
      include 'param.f'
      parameter (m1mh=m1m/2+1)
      common/fftcom/ifxz(13),trigz(3*m1m/2+1)                           
      dimension qkr(m1mh,m2),qki(m1mh,m2),fe(m1mh,m2)        
      common/ctrdps/amps(m2),acps(m1mh,m2),apps(m2)                     
      n1mh=n1m/2+1                                                      
c                                                                       
c   fft from physical space to wave numbers space                       
c   applied to vorticity                                                
c                                                                       
      call part1(qkr,qki)                                           
c                                                                       
c   solution of poisson equation real part                              
c                                                                       
       qkrmax=0.
      do 61 j=1,n2,n2m                                                      
      do 61 ii=1,n1mh                                                   
      qkr(ii,j)=0.
   61 continue                                                          
      call tribjj(qkr,n2,fe,n1mh)                        
      do 51 j=1,n2                                                      
      do 51 ii=1,n1mh                                                   
      qkr(ii,j)=fe(ii,j)                                                
   51 continue                                                          
c                                                                       
c   solution of poisson equation immag. part                            
c                                                                       
       qkimax=0.
      do 71 j=1,n2,n2m                                                      
      do 71 ii=1,n1mh                                                   
      qki(ii,j)=0.
   71 continue                                                          
      call tribjj(qki,n2,fe,n1mh)                        
      do 52 j=1,n2                                                      
      do 52 ii=1,n1mh                                                   
      qki(ii,j)=fe(ii,j)                                                
   52 continue                                                          
c                                                                       
c   fft-1 form wave numbers space to physical space                     
c   applied to stream function                                          
c                                                                       
      call part3(qkr,qki)                                           
      return                                                            
      end                                                               
c                                                                       
c  ****************************** subrout matps  ********************** 
c                                                                       
c   in this subr the coefficients of the poisson eq. for psi            
c   are calculated this subr. is called only at the beginning           
c   inside phini.                                                       
c                                                                       
      subroutine matps(i)                                               
      include 'param.f'
      parameter (m1mh=m1m/2+1)
      common/wavps/an1(m1mh),ank1(m1mh)                                 
      common/ctrdps/amps(m2),acps(m1mh,m2),apps(m2)                     
c                                                                       
c  ******** modified wave number                                        
c                                                                       
      akk1=2.*(1.-cos(an1(i)/float(n1m)))*dx1q                           
      ank1(i)=akk1                                                       
c                                                                       
c   tridiagonal matrix coefficients at each k                           
c   and cartesian coordinates in x2                                     
c                                                                       
      do 1 jc=2,n2m                                                     
        acps(i,jc)=+dx2q*(rm(jc)/g2rm(jc)
     1                   +rm(jc-1)/g2rm(jc-1))/rc(jc)/g2rc(jc)
     1   +akk1/rc(jc)**2             
    1 continue                                                          
c
c   boundary conditions
c
        acps(i,1)=-1.
        acps(i,n2)=-1.
      return                                                            
      end                                                               
c                                                                       
c   ********************* subr psini                                    
c   coefficients to invert the tridiagonal matrix
c   of Laplacian psi
c                                                                       
      subroutine psini                                                  
      include 'param.f'
      parameter (m1mh=m1m/2+1)
      common/wavps/an1(m1mh),ank1(m1mh)                                 
      common/ctrdps/amps(m2),acps(m1mh,m2),apps(m2)                     
c                                                                       
      n1mh=n1m/2+1                                                      
      call fftqps                                                       
      do 1 jc=2,n2m                                                     
        amps(jc)=-dx2q*rm(jc-1)/g2rm(jc-1)/(rc(jc)*g2rc(jc))
        apps(jc)=-dx2q*rm(jc)/g2rm(jc)/(rc(jc)*g2rc(jc))                                 
    1 continue                                                          
c     do 2 jc=1,n2,n2m                                                  
        amps(1)=0.
        amps(n2)=-1.
        apps(1)=1.                                                    
        apps(n2)=0.                                                    
c   2 continue                                                          
      do 11 i=1,n1mh                                                    
        call matps(i)                                                   
   11 continue                                                          
      return                                                            
      end                                                               
c                                                                       
c     ************** part1 ************************************         
c                                                                       
      subroutine part1(rea,aim)                                     
      include 'param.f'
      parameter (m1mh=m1m/2+1,m1md=m1m+2)
      common/fftcom/ifxz(13),trigz(3*m1m/2+1)                           
      dimension xrr(m1md,m2),workk(m1,m2),                                     
     &aim(m1m/2+1,m2),rea(m1m/2+1,m2)                                   
c
c     vorticity from the physical to wave number in the
c     theta direction
c    
      n1mh=n1m/2+1                                                      
      n1md=n1m+2                                                        
      n1mdu=n1m-1                                                       
      do 1 j=1,n2                                                       
      xrr(1,j)=vor(n1m,j)                                              
      if(j.eq.1.or.j.eq.n2) xrr(1,j)=0.                                  
c                                                                       
c                                                                       
      do 2 i=1,n1m                                                      
      is=i+1                                                            
      xrr(is,j)=vor(i,j)                                               
      if(j.eq.1.or.j.eq.n2) xrr(is,j)=0.                                 
    2 continue                                                          
      xrr(n1md,j)=vor(1,j)                                             
      if(j.eq.1.or.j.eq.n2) xrr(n1md,j)=0.                               
    1 continue                                                          
c
c   real fft from physical to wave number space
c
      call fft99(xrr,workk,trigz,ifxz,1,m1md,n1m,n2,-1)                   
      qkrmax=0.
      qkimax=0.
      do 592 i=1,n1m/2+1                                                
      ip=2*i                                                            
      id=2*i-1                                                          
      do 592 j=1,n2                                                     
      rea(i,j)=xrr(id,j)                                                 
      aim(i,j)=xrr(ip,j)                                                 
   51 continue                                                          
 592  continue                                                          
      return                                                            
      end                                                               
c     ****************** part3 ***********************************      
c     from psi in wave number to physical space                         
c                                                                       
      subroutine part3(rea,aim)                                     
      include 'param.f'
      parameter (m1mh=m1m/2+1,m1md=m1m+2)
      dimension rea(m1m/2+1,m2),                                  
     &aim(m1m/2+1,m2)                                                   
      common/fftcom/ifxz(13),trigz(3*m1m/2+1)                           
      dimension xrr(m1md,m2),workk(m1,m2)                                      
      n1mh=n1m/2+1                                                      
      n1md=n1m+2                                                        
      n1mdu=n1m-1                                                       
c
c     streamfunction in wave number  obtained by Laplacian
c                                                                       
      do 584 i=1,n1m/2+1                                                
      id=2*i-1                                                          
      ip=2*i                                                            
      do 584 j=1,n2                                                     
      xrr(id,j)=rea(i,j)                                                 
      xrr(ip,j)=aim(i,j)                                                 
 584  continue                                                          
c
c    real fft from wave to physical spoace
c
      call fft99(xrr,workk,trigz,ifxz,1,m1md,n1m,n2,1)                    
      do 579 j=1,n2                                                     
      psi(n1m,j)=xrr(1,j)                                                
      do 579 i=1,n1m-1                                                  
      is=i+1                                                            
      psi(i,j)=xrr(is,j)                                                 
 579  continue                                                          
c                                                                       
      return                                                            
      end                                                               
c
c   BESSEL FUNCTION j0   from NUMERICAL RECEPIES
c
      function bessj0(x)                                                
c
      real*8 y,p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,         
     *    s1,s2,s3,s4,s5,s6                                             
      data p1,p2,p3,p4,p5/1.d0,-.1098628627d-2,.2734510407d-4,          
     *    -.2073370639d-5,.2093887211d-6/, q1,q2,q3,q4,q5/-.1562499995d-
     *1,                                                                
     *    .1430488765d-3,-.6911147651d-5,.7621095161d-6,-.934945152d-7/ 
      data r1,r2,r3,r4,r5,r6/57568490574.d0,-13362590354.d0,651619640.7d
     *0,                                                                
     *    -11214424.18d0,77392.33017d0,-184.9052456d0/,                 
     *    s1,s2,s3,s4,s5,s6/57568490411.d0,1029532985.d0,               
     *    9494680.718d0,59272.64853d0,267.8532712d0,1.d0/               
      if(abs(x).lt.8.)then                                              
        y=x**2                                                          
        bessj0=(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))                    
     *      /(s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6)))))                      
      else                                                              
        ax=abs(x)                                                       
        z=8./ax                                                         
        y=z**2                                                          
        xx=ax-.785398164                                                
        bessj0=sqrt(.636619772/ax)*(cos(xx)*(p1+y*(p2+y*(p3+y*(p4+y     
     *      *p5))))-z*sin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))           
      endif                                                             
      return                                                            
      end                                                               
c
c
c   BESSEL FUNCTION j1   from NUMERICAL RECEPIES
c
      function bessj1(x)                                                
      real*8 y,p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,         
     *    s1,s2,s3,s4,s5,s6                                             
      data r1,r2,r3,r4,r5,r6/72362614232.d0,-7895059235.d0,242396853.1d0
     *,                                                                 
     *    -2972611.439d0,15704.48260d0,-30.16036606d0/,                 
     *    s1,s2,s3,s4,s5,s6/144725228442.d0,2300535178.d0,              
     *    18583304.74d0,99447.43394d0,376.9991397d0,1.d0/               
      data p1,p2,p3,p4,p5/1.d0,.183105d-2,-.3516396496d-4,.2457520174d-5
     *,                                                                 
     *    -.240337019d-6/, q1,q2,q3,q4,q5/.04687499995d0,-.2002690873d-3
     *,                                                                 
     *    .8449199096d-5,-.88228987d-6,.105787412d-6/                   
      if(abs(x).lt.8.)then                                              
        y=x**2                                                          
        bessj1=x*(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))                  
     *      /(s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6)))))                      
      else                                                              
        ax=abs(x)                                                       
        z=8./ax                                                         
        y=z**2                                                          
        xx=ax-2.356194491                                               
        bessj1=sqrt(.636619772/ax)*(cos(xx)*(p1+y*(p2+y*(p3+y*(p4+y     
     *      *p5))))-z*sin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))           
     *      *sign(1.,x)                                                 
      endif                                                             
      return                                                            
      end                                                               
