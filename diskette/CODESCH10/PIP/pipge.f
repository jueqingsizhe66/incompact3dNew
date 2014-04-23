c************************************************************************
c                                                                       *
c  ****************************** subrout divgu  **********************  *
c                                                                       *
c     this subroutine calculates divgu(q).                               *
c                                                                       *
c************************************************************************
c
c     Here it is performed the calculation of the divergence of
c     the intermediate velocity field. This term is the source in the
c     Poisson equation for the pressure correction and for the
c     projection of the non-solenoidal velocity field into a solenoidal one
c
      subroutine divgu(al)                                   
      include 'param.f'
      do kc=1,n3m                                                    
      kp=kpv(kc)                                                        
      do jc=1,n2m                                                    
      jp=jc+1                                                           
      usrnu1=dx1
      usrnu2=dx2/g2rm(jc)*rm(jc)
      usrnu3=dx3*rm(jc)*rm(jc)
      do ic=1,n1m                                                   
      ip=ipv(ic)                                                        
      dqcap=+(q1(ip,jc,kc)-q1(ic,jc,kc))*usrnu1
     1      +(q2(ic,jp,kc)-q2(ic,jc,kc))*usrnu2
     1      +(q3(ic,jc,kp)-q3(ic,jc,kc))*usrnu3                          
      qcap(ic,jc,kc)=dqcap/(dt*al)                                      
      enddo
      enddo
      enddo
      return                                                            
      end                                                               
c                                                                       *
c  ****************************** subrout divuck  ***********************
c                                                                       *
c     this subroutine checks divgu(q) to see the local residue.          *
c                                                                       *
c************************************************************************
      subroutine divuck(qmax,qtot)                             
      include 'param.f'
c                                                                       
c     This is just a check on the divergence of the free-divergent
c     velocity field. The calculation is stopped if QMAX > RESID
c                                                                       
      qtot=0.                                                           
      qmax=0.                                                           
      do kc=1,n3m                                                    
      kp=kpv(kc)                                                        
      do jc=1,n2m                                                    
      jp=jc+1                                                           
      usrnu1=dx1
      usrnu2=dx2/g2rm(jc)*rm(jc)
      usrnu3=dx3*rm(jc)*rm(jc)
      uvol=1./(dx1*dx3*usrnu2)
      do ic=1,n1m                                                   
      ip=ipv(ic)                                                        
      dqca1= (q1(ip,jc,kc)-q1(ic,jc,kc))*usrnu1
      dqca2= (q2(ic,jp,kc)-q2(ic,jc,kc))*usrnu2
      dqca3= (q3(ic,jc,kp)-q3(ic,jc,kc))*usrnu3
      dqcap= dqca1+dqca2+dqca3
      qtot=qtot+dqcap*uvol
      qmax=amax1(abs(dqcap),qmax)
      enddo
      enddo
      enddo
c     write(62,162)imxq,jmxq,kmxq,qtot,qmax,dq1m,dq2m,dq3m
  162 format(1x,' DIVGCK',3i4,2x,5e12.4)
  163 format(1x,' STOP IN DIVGCK',3i4,2x,5e12.4)
c
c    divergence to high
c    the calculation stops
c
      if(qmax.gt..1e-01) then
      write(6,163)imxq,jmxq,kmxq,qtot,qmax,dq1m,dq2m,dq3m
      stop
                         endif
      return                                                            
      end                                                               
c                                                                       *
c  ****************************** subrout divgco  ***********************
c                                                                       *
c     this subroutine checks divg(q) for the field obtained
c     in a previous simulation with radial coordinate calculated
c     in cordino
c                                                                       *
c************************************************************************
      subroutine divgco(qmax,qtot,n2o,dx1o,dx3o,n1lm,n3lm)
      include 'param.f'
c                                                                       
c     This is just a check on the divergence 
c                                                                       
      n2mo=n2o-1
      dx2o=float(n2mo)
      qtot=0.                                                           
      qmax=0.                                                           
      uvol=1./(dx1o*dx2o*dx3o)
      do 11 kc=1,n3lm                                                    
      kp=kpv(kc)                                                        
      do 11 jc=1,n2mo                                                   
      jp=jc+1                                                           
      usrnu1=dx1
      usrnu2=dx2o/g2rmo(jc)*rmo(jc)
      usrnu3=dx3o*rmo(jc)*rmo(jc)
      do 11 ic=1,n1lm                                                   
      ip=ipv(ic)                                                        
      dqca1= (q1(ip,jc,kc)-q1(ic,jc,kc))*usrnu1
      dqca2= (q2(ic,jp,kc)-q2(ic,jc,kc))*usrnu2
      dqca3= (q3(ic,jc,kp)-q3(ic,jc,kc))*usrnu3
      dqcap= dqca1+dqca2+dqca3
      qtot=qtot+dqcap*rmo(jc)*g2rmo(jc)/dx2o/dx3o*uvol
      qmax=amax1(abs(dqcap),qmax)                                       
   11 continue
      return                                                            
      end                                                               
c***********************************************************************
c                                                                       *
c  ****************************** subrout indic **********************  *
c                                                                       *
c     in this subroutine the indices ip,im,jp,jm,kp,km are calculated.  *
c                                                                       *
c************************************************************************
      subroutine indic                                                  
      include 'param.f'
c                                                                       
c   azimuthal periodic direction                                               
c                                                                       
      do 1 ic=1,n1m                                                     
      imv(ic)=ic-1                                                      
      if(ic.eq.1) imv(ic)=n1m                                           
      ipv(ic)=ic+1                                                      
      if(ic.eq.n1m) ipv(ic)=1                                           
    1 continue                                                          
c                                                                       
c   vertical periodic direction                                               
c                                                                       
      do 4 kc=1,n3m
      kmv(kc)=kc-1                                                      
      kpv(kc)=kc+1                                                     
      if(kc.eq.1) kmv(kc)=n3m                                           
      if(kc.eq.n3m) kpv(kc)=1                                           
    4 continue                                                          
c                                                                       
c     direction normal to the radial boundary 
c                                                                       
      do 3 jc=1,n2m                                                     
      jmv(jc)=jc-1                                                      
      jpv(jc)=jc+1                                                      
      if(jc.eq.1) jmv(jc)=jc                                            
      if(jc.eq.n2m) jpv(jc)=jc                                          
    3 continue                                                          
c                                                                       
c   indices for the axis of symmetry and the external wall                                                   
c                                                                       
      do 15 jc=1,n2m                                                    
      jpc(jc)=jpv(jc)-jc                                                
      jmc(jc)=jc-jmv(jc)                                                
      jup(jc)=1-jpc(jc)                                                 
      jum(jc)=1-jmc(jc)                                                 
   15 continue                                                          
c
c    azimuthal index to evaluate the velocity at the axis in the
c    non-linear terms
c
      do i=1,n1m
      isym(i) = i + n1m/2
      if(isym(i).gt.n1m) isym(i) = isym(i) - n1m
      enddo
      return                                                            
      end                                                               
c************************************************************************
c                                                                       *
c  ****************************** subrout cflu  **********************   *
c                                                                       *
c************************************************************************
      subroutine cflu(cflm)                                     
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
c                                                                       
      do k=1,n3m                                                      
      kp=kpv(k)                                                      
      do j=1,n2m                                                      
      jp=jpv(j)                                                         
      usrnu=dx2/g2rm(j)
      do i=1,n1m                                                    
      ip=ipv(i)                                                         
      qcf=( abs(q1(i,j,k)+q1(ip,j,k))*0.5*dx1/rm(j)**2
     1     +abs(q2(i,j,k)+q2(i,jp,k))*0.5*usrnu/rm(j)                   
     1     +abs(q3(i,j,k)+q3(i,j,kp))*0.5*dx3  )
      cflm=max(cflm,qcf)                                              
      enddo
      enddo
      enddo
      return                                                            
      end                                                               
