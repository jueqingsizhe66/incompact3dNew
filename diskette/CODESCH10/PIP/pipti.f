c************************************************************************
c                                                                       *
c  ************************ subrout invinv  **********************      *
c                                                                       *
c  this subroutine performs the inversion of the q_i momentum equations *
c  by the third order low storage runge-kutta scheme                    *
c  INVISCID CASE to check energy conservation property of
c  the numerical scheme
c                                                                       *
c************************************************************************
      subroutine invinv(al,ga,ro,ns)                      
      include 'param.f'
      common/timw/timew
c                                                                       
c  ********* compute the q1 component 
c  dq=qhat-q(n)                                                         
c                                                                       
      q1max=0.
      q2max=0.
      q3max=0.
      gp1max=0.
      do kc=1,n3m                                                   
      do jc=1,n2m                                                    
      do ic=1,n1m                                                   
      im=imv(ic)
c                                                                       
c   grad(pr) along 1                                                    
c                                                                       
      dpx11=(pr(ic,jc,kc)-pr(im,jc,kc))*dx1
      gp1max=max(abs(dpx11),gp1max)
      gradp=dpx11*al
      rhsdd=(ga*dq(ic,jc,kc)+ro*ru1(ic,jc,kc)                   
     1          -gradp)*dt
         q1(ic,jc,kc)=q1(ic,jc,kc) + rhsdd
      q1max=max(abs(q1(ic,jc,kc)),q1max)
      ru1(ic,jc,kc)=dq(ic,jc,kc)
      enddo                                                             
      enddo                                                             
      enddo                                                             
c                                                                       
c  ********* compute the q2 component 
c  dq=qhat-q(n)                                                         
c                                                                       
      do kc=1,n3m                                                    
            do jc=2,n2m                                                   
            jm=jc-1                                                           
                do ic=1,n1m                                                   
c                                                                       
c   component of grap(pr) along 2 direction                             
c                                                                       
      dpx22=(pr(ic,jc,kc)-pr(ic,jm,kc))*dx2/g2rc(jc)
      gradp=dpx22*al*rc(jc)
      rhsdd=(ga*dph(ic,jc,kc)+ro*ru2(ic,jc,kc)                   
     1      -gradp)*dt                              
         q2(ic,jc,kc)=q2(ic,jc,kc) + rhsdd
      q2max=max(abs(q2(ic,jc,kc)),q2max)
      ru2(ic,jc,kc)=dph(ic,jc,kc)                                        
                 enddo
             enddo
      enddo
        do ic=1,n1m
            do kc=1,n3m                                                     
      q2(ic,1,kc)=0.
      q2(ic,n2,kc)=0.
            enddo                                                             
        enddo                                                             
c                                                                       
c  ********* compute the q3 component 
c  everything at i,j+1/2,k+1/2                                          
c                                                                       
      if(n3m.gt.1) then
      q3tot=0.
      slorf=0.
      do kc=1,n3m
            do jc=1,n2m
                  do ic=1,n1m
		  q3tot=q3tot+q3(ic,jc,kc)*volz(jc)
		  enddo
            enddo
      enddo
      do kc=1,n3m
            do jc=1,n2m
                  do ic=1,n1m
                  im=imv(ic)
                  ip=ipv(ic)
         iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2m
      rhs(iadd)=(ga*qcap(ic,jc,kc)+ro*ru3(ic,jc,kc))*dt
      slorf=slorf+(ga*qcap(ic,jc,kc)+ro*ru3(ic,jc,kc))*volz(jc)/al
      ru3(ic,jc,kc)=qcap(ic,jc,kc)
                  enddo                                                             
            enddo                                                             
      enddo                                                             
       dplfo=slorf/(pi*alx3d)
       dptot=dplfo
      do kc=1,n3m
      km=kmv(kc)
            do jc=1,n2m
                  do ic=1,n1m
c                                                                       
c  component of grad(pr) along x3 direction                             
c                                                                       
      dpx33=(pr(ic,jc,kc)-pr(ic,jc,km))*dx3
      gradp=(dpx33+dptot)*al
         iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2m
      rhs(iadd)=rhs(iadd)-gradp*dt
         q3(ic,jc,kc)=q3(ic,jc,kc) + rhs(iadd)
      q3max=max(abs(q3(ic,jc,kc)),q3max)
                  enddo  
            enddo                                                             
      enddo                                                             
       q3tot=q3tot/(pi*alx3d)
       if(ns.eq.3) then
       write(75,109) timew,dplfo,q3tot,gp1max,q1max,q2max,q3max
  109 format(3x,'inviscid ',7e12.4)
                   endif
                         endif
      return                                                            
      end                                                               
