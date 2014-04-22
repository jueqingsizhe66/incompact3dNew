c************************************************************************
c                                                                       *
c ****************************** subrout hdnl1  **********************  *
c    non-linear terms  aand                     
c    first derivatives in the viscous terms   *
c    To speed-up the code the if statements to treat the axis 
c    can be eliminated by doing a further do loop
c
c*************************************************************************
      subroutine hdnl1
      include 'param.f'
c 
c     compute the non-linear terms by centered difference.              
c     h term for the q1 momentum equation at i,j+1/2,k+1/2              
c
      do 10 ic=1,n1m                                                    
      ip=ipv(ic)                                                        
      im=imv(ic)                                                        
      do 10 kc=1,n3m                                                    
      km=kmv(kc)                                                        
      kp=kpv(kc)                                                        
      do 10 jc=1,n2m                                                    
      jmm=jmv(jc)                                                       
      jpp=jpv(jc)                                                       
      jp=jc+1                                                           
c                                                                       
c    q1 q1 term                                                         
c                                                                       
      h11=( (q1(ip,jc,kc)+q1(ic,jc,kc))*(q1(ip,jc,kc)+q1(ic,jc,kc))     
     1     -(q1(im,jc,kc)+q1(ic,jc,kc))*(q1(im,jc,kc)+q1(ic,jc,kc))     
     1    )*dx1/rm(jc)**2*0.25                                             
c                                                                       
c   q1 q2 term  usual                                                   
c                                                                       
       if(jc.ne.1) then
      q1nn=(q1(ic,jpp,kc)/rm(jpp)+q1(ic,jc,kc)/rm(jc))*0.5
      q1ss=(q1(ic,jc,kc)/rm(jc)+q1(ic,jmm,kc)/rm(jmm))*0.5
      h12d=( (q2(ic,jp,kc)+q2(im,jp,kc))*q1nn
     1      -(q2(ic,jc,kc)+q2(im,jc,kc))*q1ss
     1     )*dx2*0.5/g2rm(jc)
        else
      q1nn=(q1(ic,2,kc)/rm(2)+q1(ic,1,kc)/rm(1))*0.5
      h12d=( (q2(ic,2,kc)+q2(im,2,kc))*q1nn
     1     )*dx2*0.5/g2rm(jc)
       endif
c
c                                                                       
c   q1 q2 term  not differentiated                                      
c                                                                       
       if(jc.ne.1) then
      h12n=  q1(ic,jc,kc) * (
     1        q2(ic,jp,kc)/rc(jp)+q2(ic,jc,kc)/rc(jc)+
     1        q2(im,jp,kc)/rc(jp)+q2(im,jc,kc)/rc(jc)  )*0.25 / rm(jc)   
        else
       imsy = isym(imv(ic))
       q2s1 = (q2(ic,2,kc) - q2(isym(ic),2,kc))*0.5/rc(2)
       q2s2 = (q2(im,2,kc) - q2(imsy,2,kc))*0.5/rc(2)
      h12n=  q1(ic,1,kc) * (
     1        ((q2(ic,2,kc)+q2(im,2,kc))/rc(2)+q2s1 + q2s2)  
     1      )*0.25 / rm(1)
       end if
c                                                                       
c   q1 q3 term                                                          
c                                                                       
      h13=((q3(ic,jc,kp)+q3(im,jc,kp))*(q1(ic,jc,kp)+q1(ic,jc,kc))      
     3    -(q3(ic,jc,kc)+q3(im,jc,kc))*(q1(ic,jc,kc)+q1(ic,jc,km))      
     1    )*dx3*0.25                                                    
c
c     Coriolis term  (BE CAREFULL!!!!! THIS ROSSBY IS DEFINED AS THE
c                     INVERSE OF THE CONVENTIONAL ONE)
c
c           Ros  q_r
c
      htco=   +ros *(  q2(ic,jp,kc)+q2(ic,jc,kc)
     1               + q2(im,jp,kc)+q2(im,jc,kc)
     1                        )*.25

      hq1=(h11+h12d+h12n)+h13+htco                                      
      dq(ic,jc,kc)=-hq1                                                  
   10 continue                                                          
      if(ren.lt.1.e17) then
      do 11 ic=1,n1m                                                    
      ip=ipv(ic)                                                        
      im=imv(ic)                                                        
      do 11 kc=1,n3m                                                    
      km=kmv(kc)                                                        
      kp=kpv(kc)                                                        
      do 11 jc=1,n2m                                                    
      jp=jc+1                                                           
c                                                                       
c   first derivative of q2 with respect to x1                           
c                                                                       
      if(jc.eq.1) then
       imsy = isym(imv(ic))
       q2s1 = (q2(ic,2,kc) - q2(isym(ic),2,kc))*0.5/rc(2)
       q2s2 = (q2(im,2,kc) - q2(imsy,2,kc))*0.5/rc(2)
         q2e=(q2(ic,2,kc)/rc(2) + q2s1      )*0.5                             
         q2w=(q2(im,2,kc)/rc(2) + q2s2      )*0.5                             
                   else
      q2e=(q2(ic,jp,kc)/rc(jp)+q2(ic,jc,kc)/rc(jc))*0.5 
      q2w=(q2(im,jp,kc)/rc(jp)+q2(im,jc,kc)/rc(jc))*0.5
                   endif
      d11q1e=2.*(q2e-q2w)*dx1/rm(jc)                                 
      dq(ic,jc,kc)=d11q1e/ren+dq(ic,jc,kc)                              
   11 continue                                                          
      end if
      return                                                            
      end                                                               
c************************************************************************
c                                                                       *
c ****************************** subrout hdnl2  **********************  *
c    the non-linear terms and 
c    first derivatives in the viscous terms   *
c    To speed-up the code the if statements to treat the point near the axis 
c    can be eliminated by doing a further do loop
c                                                                       *
c************************************************************************
      subroutine hdnl2
      include 'param.f'
c                                                                       
c     compute the non-linear terms by centered difference.              
c     h term for the q2 momentum equation at i+1/2,j,k+1/2              
c
      do 20 ic=1,n1m                                                    
      ip=ipv(ic)                                                        
      imm=imv(ic)                                                       
      ipp=ipv(ic)                                                       
      do 20 kc=1,n3m                                                    
      km=kmv(kc)                                                        
      kp=kpv(kc)                                                        
      do 20 jc=2,n2m                                                    
      jm=jc-1                                                           
      jp=jc+1                                                           
c                                                                       
c     q2 q1 term                                                        
c                                                                       
      h21=( (q1(ipp,jc,kc)/rm(jc)+q1(ipp,jm,kc)/rm(jm)) 
     1     *(q2(ipp,jc,kc)+q2(ic,jc,kc))                                
     1     -(q1(ic,jc,kc)/rm(jc)+q1(ic,jm,kc)/rm(jm))  
     1     *(q2(ic,jc,kc)+q2(imm,jc,kc))                                
     1    )*dx1*0.25/rc(jc)
c                                                                       
c     q2 q2 term                                                        
c                                                                       
      if(jc.eq.2) then
       imsy = isym(imv(ic))
       q2s1 = (q2(ic,2,kc) - q2(isym(ic),2,kc))*0.5/rc(2)
      h22d=( (q2(ic,3,kc)/rc(3)+q2(ic,2,kc)/rc(2)) 
     1     *(q2(ic,3,kc)+q2(ic,2,kc))
     1     -(q2(ic,2,kc)/rc(2)+q2s1) 
     1     *(q2(ic,2,kc)+q2(ic,1,kc))
     1    )*dx2*0.25/g2rc(jc)                                                    
       else
      h22d=( (q2(ic,jp,kc)/rc(jp)+q2(ic,jc,kc)/rc(jc))  
     1     *(q2(ic,jp,kc)+q2(ic,jc,kc))
     1     -(q2(ic,jc,kc)/rc(jc)+q2(ic,jm,kc)/rc(jm))  
     1     *(q2(ic,jc,kc)+q2(ic,jm,kc))
     1    )*dx2*0.25/g2rc(jc)
       endif
c                                                                       
c   q1 q1 term  not differentiated                                      
c                                                                       
      h11n=( (q1(ipp,jc,kc)/rm(jc)+q1(ipp,jm,kc)/rm(jm)
     1       +q1(ic,jc,kc)/rm(jc)+q1(ic,jm,kc)/rm(jm))*.25)**2
c                                                                       
c     q2 q3 term                                                        
c                                                                       
      h23=((q3(ic,jc,kp)+q3(ic,jm,kp))*(q2(ic,jc,kp)+q2(ic,jc,kc))      
     3    -(q3(ic,jc,kc)+q3(ic,jm,kc))*(q2(ic,jc,kc)+q2(ic,jc,km))      
     1    )*dx3/4.                                                      
c
c                    (BE CAREFULL!!!!! THIS ROSSBY IS DEFINED AS THE
c                     INVERSE OF THE CONVENTIONAL ONE)
c
c          -Ros  q_t
c
      hrco=-ros
     1     * (q1(ip,jc,kc)+q1(ip,jm,kc)
     1       +q1(ic,jc,kc)+q1(ic,jm,kc))*.25

      hq2=(h21+h22d)+h23-h11n+hrco
      dph(ic,jc,kc)=-hq2                                                  
   20 continue                                                          
      if(ren.lt.1.e17) then
      do 21 ic=1,n1m                                                    
      imm=imv(ic)                                                       
      ipp=ipv(ic)                                                       
      do 21 kc=1,n3m                                                    
      km=kmv(kc)                                                        
      kp=kpv(kc)                                                        
      do 21 jc=2,n2m                                                    
      jm=jc-1                                                           
      jp=jc+1                                                           
      jpp=jpv(jc)                                                       
c                                                                       
c   second derivative of q1 with respect to x2                          
c                                                                       
      q1e=(q1(ipp,jc,kc)/rm(jc)+q1(ipp,jm,kc)/rm(jm))*0.5 
      q1w=(q1(ic,jc,kc)/rm(jc)+q1(ic,jm,kc)/rm(jm))*0.5  
      d11q2e=-2.*(q1e-q1w)*dx1/rc(jc)                                   
      dph(ic,jc,kc)=d11q2e/ren+dph(ic,jc,kc)                              
   21 continue                                                          
      endif 
      return                                                            
      end                                                               
c************************************************************************
c                                                                       *
c ****************************** subrout hdnl3  **********************  *
c    non-linear terms  aand                     
c    first derivatives in the viscous terms   *
c                                                                       *
c************************************************************************
      subroutine hdnl3
      include 'param.f'
c                                                                       
c     compute the non-linear terms by centered difference.              
c                                                                       
c     h term for the q3 momentum equation at i+1/2,j+1/2,k              
c
      do 30 ic=1,n1m                                                    
      imm=imv(ic)                                                       
      ipp=ipv(ic)                                                       
      do 30 kc=1,n3m                                                    
      km=kmv(kc)                                                        
      kp=kpv(kc)                                                        
      do 30 jc=1,n2m                                                    
      jp=jc+1                                                           
      jmm=jmv(jc)                                                       
      jpp=jpv(jc)                                                       
      ugmmm=1./rm(jc)                                                   
c                                                                       
c    q3 q2 term                                                         
c    with q2=0 on the walls                                             
c                                                                       
      h32=(  (  (q2(ic,jp,kc)+q2(ic,jp,km))*0.5                         
     1         *(q3(ic,jpp,kc)+q3(ic,jc,kc))*0.5  )                     
     1      -(  (q2(ic,jc,kc)+q2(ic,jc,km))*0.5                         
     1         *(q3(ic,jc,kc)+q3(ic,jmm,kc))*0.5  )                     
     1    )*dx2/g2rm(jc)
c                                                                       
c    q3 q1 term                                                         
c                                                                       
      h31=(       (q1(ipp,jc,kc)+q1(ipp,jc,km))                         
     1        *   (q3(ipp,jc,kc)+q3(ic,jc,kc))                          
     1           -(q1(ic,jc,kc)+q1(ic,jc,km))                           
     1        *   (q3(ic,jc,kc)+q3(imm,jc,kc))                          
     1    )*dx1*0.25/rm(jc)                                                    
c    q3 q3 term                                                         
c                                                                       
      h33=((q3(ic,jc,kp)+q3(ic,jc,kc))*(q3(ic,jc,kp)+q3(ic,jc,kc))      
     3    -(q3(ic,jc,kc)+q3(ic,jc,km))*(q3(ic,jc,kc)+q3(ic,jc,km))      
     1    )*dx3*0.25                                                    
      hq3=(h31+h32)*ugmmm+h33                                           
      qcap(ic,jc,kc)=-hq3                                                 
   30 continue                                                          
      return                                                            
      end                                                               
