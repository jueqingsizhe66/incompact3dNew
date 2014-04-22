c  ****************************** subrout strper  **********************  
c                                                                       
c     this subroutine calculates the  strain rate tensor at the center of
c     the cell. Perodic box
c                                                                       
      subroutine strper(q)                    
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m                                   
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q                            
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)      
      dimension q(ndv,m1,m2,m3)
      common/strain/st(m1,m2,m3,6)
      common/rhsc/rhs(m1,m2,m3)
      common/sij2po/sijji
      vl123=1./float(n1m*n2m*n3m)

c
c   st(l  =Sll  at i+1/2,j+1/2,k+1/2
c   st(2  =S22  at i+1/2,j+1/2,k+1/2
c   st(3  =S33  at i+1/2,j+1/2,k+1/2
c
      sijji=0.
      do 4 kc=1,n3m 
      kp=kpv(kc)                                                    
      do 4 jc=1,n2m                                                     
      jp=jpv(jc)
      do 4 ic=1,n1m                                                     
      ip=ipv(ic)                                                        
      st(ic,jc,kc,1)=(q(1,ip,jc,kc)-q(1,ic,jc,kc))*dx1
      st(ic,jc,kc,2)=(q(2,ic,jp,kc)-q(2,ic,jc,kc))*dx2
      st(ic,jc,kc,3)=(q(3,ic,jc,kp)-q(3,ic,jc,kc))*dx3
      sijji=sijji+ 
     1    st(ic,jc,kc,1)*st(ic,jc,kc,1)+
     1    st(ic,jc,kc,2)*st(ic,jc,kc,2)+
     1    st(ic,jc,kc,3)*st(ic,jc,kc,3)

    4 continue                                                          
c
c   st(4 = S12    at  i,j,k+1/2
c   st(6 = S32    at  i+1/2,j,k
c   st(5 = S13    at  i,j+1/2,k
c
      do 6 kc=1,n3m 
      km=kmv(kc)                                                    
      do 6 jc=1,n2m                                                     
      jm=jmv(jc)
      do 6 ic=1,n1m                                                     
      im=imv(ic)                                                        
c
      st(ic,jc,kc,4)=((q(1,ic,jc,kc)-q(1,ic,jm,kc))*dx2
     1              +(q(2,ic,jc,kc)-q(2,im,jc,kc))*dx1)*0.5
c
      st(ic,jc,kc,5)=((q(1,ic,jc,kc)-q(1,ic,jc,km))*dx3
     1              +(q(3,ic,jc,kc)-q(3,im,jc,kc))*dx1)*0.5
c
      st(ic,jc,kc,6)=((q(3,ic,jc,kc)-q(3,ic,jm,kc))*dx2
     1              +(q(2,ic,jc,kc)-q(2,ic,jc,km))*dx3)*0.5
      sijji=sijji+2.*(
     1    st(ic,jc,kc,4)*st(ic,jc,kc,4)+
     1    st(ic,jc,kc,5)*st(ic,jc,kc,5)+
     1    st(ic,jc,kc,6)*st(ic,jc,kc,6) )
    6 continue                                                          
      sijji=sijji*vl123 
c
c   st(4 = S12    at  i,1,k+1/2
c   st(6 = S32    at  i+1/2,1,k
c
c   evaluate s12 at the cell center
c
      do 20 kc=1,n3m
      do 20 jc=1,n2m   
      jp=jpv(jc)
      do 20 ic=1,n1m             
      ip=ipv(ic)                                                        
      rhs(ic,jc,kc)=(st(ic,jc,kc,4)+st(ic,jp,kc,4)+
     1               st(ip,jc,kc,4)+st(ip,jp,kc,4))*0.25
   20 continue                             
      do 21 kc=1,n3m
      do 21 jc=1,n2m   
      do 21 ic=1,n1m             
      st(ic,jc,kc,4)=rhs(ic,jc,kc)           
   21 continue                             
c
c   evaluate s13 at the cell center
c
      do 30 kc=1,n3m
      kp=kpv(kc)                                                    
      do 30 jc=1,n2m   
      do 30 ic=1,n1m             
      ip=ipv(ic)                                                        
      rhs(ic,jc,kc)=(st(ic,jc,kc,5)+st(ic,jc,kp,5)+
     1               st(ip,jc,kc,5)+st(ip,jc,kp,5))*0.25
   30 continue                             
      do 35 kc=1,n3m
      do 35 jc=1,n2m   
      do 35 ic=1,n1m             
      st(ic,jc,kc,5)=rhs(ic,jc,kc)           
   35 continue                             
c
c   evaluate s23 at the cell center
c
      do 50 kc=1,n3m
      kp=kpv(kc)                                                    
      do 50 jc=1,n2m   
      jp=jpv(jc)
      do 50 ic=1,n1m             
      rhs(ic,jc,kc)=(st(ic,jc,kc,6)+st(ic,jp,kc,6)+
     1               st(ic,jc,kp,6)+st(ic,jp,kp,6))*0.25
   50 continue                             
      do 57 kc=1,n3m
      do 57 jc=1,n2m   
      do 57 ic=1,n1m             
      st(ic,jc,kc,6)=rhs(ic,jc,kc)           
   57 continue                             
      return                                                            
      end                                                               
