c***********************************************
c     Quantities for LES s.g.s
c***********************************************
      subroutine quales
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      common/d13/alx1,alx3
      common/metria/caj(m2),cac(m2)
      common/qles1/deltax1,deltay1,deltaz1,ell1c(m2)
      common/qles2/deltax2,deltay2,deltaz2,ell2c(m2)
c      
c     calculate two length scale
      deltax1=alx1/float(n1m)
      deltaz1=alx3/float(n3m)
      deltax2=2.*deltax1
      deltaz2=2.*deltaz1
      do 5 jc=1,n2m
       ell1c(jc)=(deltax1*deltaz1*caj(jc)/dx2)**.66667
       ell2c(jc)=(deltax2*deltaz2*caj(jc)/dx2)**.66667
   5  continue
      return
      end
c***********************************************
c     dynamic subgrid-scale stress model 
c     Lilley version
c***********************************************
      subroutine dynamic(q,pr,ncount)
      include 'param.f'
      common/rhsc/rhs(m1,m2,m3)
      common/sma/vis(m1,m2,m3)
      dimension q(ndv,m1,m2,m3),pr(m1,m2,m3)
      real*4 mij,lijmij,mijmij
      common/strain/st(m1,m2,m3,6)
      common/mijdy/mij(m1,m2,m3,6)
      common/qdyn/g(m1,m2,m3)
      common/smapra/gold(m1,m2,m3)
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/metria/caj(m2),cac(m2)
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/ledat/csma,cvisc,iles
      common/outd1/cdyn(m2),alij(m2),amij(m2)
      common/outd2/cs(m2),lijmij(m2),mijmij(m2)
      common/outd3/accdyn(m2),acclij(m2),accmij(m2)
      common/qles1/deltax1,deltay1,deltaz1,ell1c(m2)
      common/qles2/deltax2,deltay2,deltaz2,ell2c(m2)
c
c      computation of Sij
c
      call strai(q)
c      
      do 10 k=1,n3m
      do 10 j=1,n2m
      do 10 i=1,n1m
      g(i,j,k)=2.*(st(i,j,k,1)**2+st(i,j,k,2)**2+st(i,j,k,3)**2)+ 
     1         4.*(st(i,j,k,4)**2+st(i,j,k,5)**2+st(i,j,k,6)**2)
      g(i,j,k)=sqrt(g(i,j,k)) 
   10 continue
c     now g contains |S|=sqrt(2*Sij*Sij) 
c     gold stores |S| for the computation of viscosity at the end
c     of the subroutine
      do 11 k=1,n3m
      do 11 j=1,n2m
      do 11 i=1,n1m
      gold(i,j,k)= g(i,j,k)
   11 continue
c     pr=|S|Sij
      do 21 n=1,6
      do 20 k=1,n3m
      do 20 j=1,n2m
      do 20 i=1,n1m
      pr(i,j,k)=st(i,j,k,n)*g(i,j,k)
 20   continue
      call filter(pr,mij(1,1,1,n))
 21   continue
c     now mij contains (|S|Sij)^
      do 31 n=1,6
      do 30 k=1,n3m
      do 30 j=1,n2m
      do 30 i=1,n1m
      pr(i,j,k)=st(i,j,k,n)
 30   continue
      call filter(pr,st(1,1,1,n))
 31   continue
c     now st contains Sij^
      do 32 k=1,n3m
      do 32 j=1,n2m
      do 32 i=1,n1m
      rhs(i,j,k)=2.*(st(i,j,k,1)**2+st(i,j,k,2)**2+st(i,j,k,3)**2)+
     1           4.*(st(i,j,k,4)**2+st(i,j,k,5)**2+st(i,j,k,6)**2)
      rhs(i,j,k)=sqrt(rhs(i,j,k)) 
 32   continue
c     now rhs contains |S|^
      do 33 n=1,6
      do 33 k=1,n3m
      do 33 j=1,n2m
      do 33 i=1,n1m
      amijlo=ell2c(j)*rhs(i,j,k)*st(i,j,k,n)-ell1c(j)*mij(i,j,k,n)
      mij(i,j,k,n)=amijlo/ell1c(j)
 33   continue
c   now mij=l2/l1*|S|^Sij^-(|S|Sij)^
c   we have to contract mij with itself
      do 41 k=1,n3m
      do 41 j=1,n2m
      do 41 i=1,n1m
      pr(i,j,k)=(mij(i,j,k,1)**2+mij(i,j,k,2)**2+mij(i,j,k,3)**2)+
     1       2.*(mij(i,j,k,4)**2+mij(i,j,k,5)**2+mij(i,j,k,6)**2)
 41   continue
c     mijmij(j) average mijmij on plane xz
      do 45 j=1,n2m
       mijmij(j)=0.
      do 46 k=1,n3m
      do 46 i=1,n1m
       mijmij(j)=mijmij(j)+pr(i,j,k)
 46   continue
c
c  accmij is averaging also in time
c
      mijmij(j)=mijmij(j)/float(n1m*n3m)
       accmij(j)=accmij(j)+mijmij(j)
       amij(j)=accmij(j)/ncount
 45   continue
c     we have to contract lij with mij
c    apply the test filter to q
c    we have to define the q at the center
c     st(.,4)<------q1
c     st(.,5)<------q2
c     st(.,6)<------q3
      do 48 k=1,n3m
      kp=kpv(k)
      do 48 j=1,n2m
      do 48 i=1,n1m
      ip=ipv(i)
      st(i,j,k,4)=.5*(q(1,ip,j,k)+q(1,i,j,k))
      st(i,j,k,5)=.5*(q(2,i,j+1,k)+q(2,i,j,k))
      st(i,j,k,6)=.5*(q(3,i,j,kp)+q(3,i,j,k))
 48   continue
      call filter(st(1,1,1,4),st(1,1,1,1))        
      call filter(st(1,1,1,5),st(1,1,1,2))        
      call filter(st(1,1,1,6),st(1,1,1,3))        
c**************************** component 11 ********************** 
      do 50 k=1,n3m
      do 50 j=1,n2m
      do 50 i=1,n1m
      pr(i,j,k)=st(i,j,k,4)*st(i,j,k,4)
 50   continue
      call filter(pr,vis)        
      do 52 k=1,n3m
      do 52 j=1,n2m
      do 52 i=1,n1m
      g(i,j,k)=(vis(i,j,k)-st(i,j,k,1)*st(i,j,k,1))*mij(i,j,k,1)
 52   continue
c**************************** component 22 ********************** 
      do 60 k=1,n3m
      do 60 j=1,n2m
      do 60 i=1,n1m
      pr(i,j,k)=st(i,j,k,5)*st(i,j,k,5)
 60   continue
      call filter(pr,vis)        
      do 62 k=1,n3m
      do 62 j=1,n2m
      do 62 i=1,n1m
      g(i,j,k)=(vis(i,j,k)-st(i,j,k,2)*st(i,j,k,2))*mij(i,j,k,2)
     1          +g(i,j,k)        
 62   continue
c**************************** component 33 ********************** 
      do 70 k=1,n3m
      do 70 j=1,n2m
      do 70 i=1,n1m
      pr(i,j,k)=st(i,j,k,6)*st(i,j,k,6)
 70   continue
      call filter(pr,vis)        
      do 72 k=1,n3m
      do 72 j=1,n2m
      do 72 i=1,n1m
      g(i,j,k)=(vis(i,j,k)-st(i,j,k,3)*st(i,j,k,3))*mij(i,j,k,3)
     1          +g(i,j,k)        
 72   continue
c**************************** component 12 ********************** 
      do 80 k=1,n3m
      do 80 j=1,n2m
      do 80 i=1,n1m
      pr(i,j,k)=st(i,j,k,4)*st(i,j,k,5)
 80   continue
      call filter(pr,vis)        
      do 82 k=1,n3m
      do 82 j=1,n2m
      do 82 i=1,n1m
      g(i,j,k)=2.*(vis(i,j,k)-st(i,j,k,1)*st(i,j,k,2))*mij(i,j,k,4)
     1        +g(i,j,k)       
 82   continue
c**************************** component 13 ********************** 
      do 90 k=1,n3m
      do 90 j=1,n2m
      do 90 i=1,n1m
      pr(i,j,k)=st(i,j,k,4)*st(i,j,k,6)
 90   continue
      call filter(pr,vis)        
      do 92 k=1,n3m
      do 92 j=1,n2m
      do 92 i=1,n1m
      g(i,j,k)=2.*(vis(i,j,k)-st(i,j,k,1)*st(i,j,k,3))*mij(i,j,k,5)
     1         +g(i,j,k)       
 92   continue
c**************************** component 23 ********************** 
      do 100 k=1,n3m
      do 100 j=1,n2m
      do 100 i=1,n1m
      pr(i,j,k)=st(i,j,k,5)*st(i,j,k,6)
100   continue
      call filter(pr,vis)        
      do 102 k=1,n3m
      do 102 j=1,n2m
      do 102 i=1,n1m
      g(i,j,k)=2.*(vis(i,j,k)-st(i,j,k,2)*st(i,j,k,3))*mij(i,j,k,6)
     1         +g(i,j,k)       
102   continue
c     now g contains lij*mij
c      average numerator lijmij on plane xz
      do 110 j=1,n2m
      lijmij(j)=0.
      do 120 k=1,n3m
      do 120 i=1,n1m
      lijmij(j)=lijmij(j)+g(i,j,k)
 120  continue
      lijmij(j)=lijmij(j)/float(n1m*n3m)
c
c  acclij is averaging also in time
c
       acclij(j)=acclij(j)+lijmij(j)
       alij(j)=acclij(j)/ncount
 110  continue
c       calculate smagorinsky constant
       do 200 j=1,n2m
c
c   here cs(j) is c_S*ell1c(j)
c
        cs(j)=-.5*lijmij(j)/mijmij(j)
       accdyn(j)=accdyn(j)+cs(j)
       cdyn(j)=accdyn(j)/ncount
 200   continue
c
c   here the constant has been averaged in time
c
       do 300 k=1,n3m
       do 300 j=1,n2m
       do 300 i=1,n1m
       vis(i,j,k)=cdyn(j)*gold(i,j,k)+cvisc
c
c   with the next the average is only in x1-x3 planes
c      vis(i,j,k)=cs(j)*gold(i,j,k)+cvisc
c
       if(vis(i,j,k).lt.0.) vis(i,j,k)=0.
 300   continue
       return
       end 
c  ****************************** subrout filter  **********************  
      subroutine filter(u,uf)
      include 'param.f'
c     calculates filtered function using a box filter in physical
c     space. No filtering in the direction normal to wall
      dimension u(m1,m2,m3),uf(m1,m2,m3)
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/dim/n1,n1m,n2,n2m,n3,n3m
      do 10 kc=1,n3m
       kp=kpv(kc)
       km=kmv(kc)
      do 10 jc=1,n2m
      do 10 ic=1,n1m
       ip=ipv(ic)
       im=imv(ic)
      uf(ic,jc,kc)=.0625*(u(im,jc,kp)+2.*u(ic,jc,kp)+u(ip,jc,kp)+
     1                 2.*u(im,jc,kc)+4.*u(ic,jc,kc)+2.*u(ip,jc,kc)+
     1                    u(im,jc,km)+2.*u(ic,jc,km)+u(ip,jc,km))
  10  continue
      return
      end
c  ****************************** subrout strain  **********************  
c                                                                       
c     this subroutine calculates the  strain rate tensor at the center of
c     the cell
c                                                                       
      subroutine strai(q)                    
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m                                   
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q                            
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)      
      dimension q(ndv,m1,m2,m3)
      common/strain/st(m1,m2,m3,6)
      common/rhsc/rhs(m1,m2,m3)
      common/metria/caj(m2),cac(m2)
      common/islwal/islv1s,islv1n,islv3s,islv3n
c
c   st(l  =Sll  at i+1/2,j+1/2,k+1/2
c   st(2  =S22  at i+1/2,j+1/2,k+1/2
c   st(3  =S33  at i+1/2,j+1/2,k+1/2
c
      do 4 kc=1,n3m 
      kp=kpv(kc)                                                    
      do 4 jc=1,n2m                                                     
      jp=jc+1
      do 4 ic=1,n1m                                                     
      ip=ipv(ic)                                                        
      st(ic,jc,kc,1)=(q(1,ip,jc,kc)-q(1,ic,jc,kc))*dx1
      st(ic,jc,kc,2)=(q(2,ic,jp,kc)-q(2,ic,jc,kc))*dx2/caj(jc)
      st(ic,jc,kc,3)=(q(3,ic,jc,kp)-q(3,ic,jc,kc))*dx3
    4 continue                                                          
c
c   st(4 = S12    at  i,j,k+1/2
c   st(6 = S32    at  i+1/2,j,k
c   st(5 = S13    at  i,j+1/2,k
c
      do 6 kc=1,n3m 
      km=kmv(kc)                                                    
      do 6 jc=2,n2m                                                     
      jm=jc-1                                                       
      do 6 ic=1,n1m                                                     
      im=imv(ic)                                                        
      st(ic,jc,kc,4)=(q(1,ic,jc,kc)-q(1,ic,jm,kc))*dx2/cac(jc)
     1               +(q(2,ic,jc,kc)-q(2,im,jc,kc))*dx1
      st(ic,jc,kc,6)=(q(3,ic,jc,kc)-q(3,ic,jm,kc))*dx2/cac(jc)
     1               +(q(2,ic,jc,kc)-q(2,ic,jc,km))*dx3
    6 continue                                                          
      do 8 kc=1,n3m 
      km=kmv(kc)                                                    
      do 8 jc=1,n2m                                                     
      do 8 ic=1,n1m                                                     
      im=imv(ic)                                                        
      st(ic,jc,kc,5)=(q(1,ic,jc,kc)-q(1,ic,jc,km))*dx3
     1              +(q(3,ic,jc,kc)-q(3,im,jc,kc))*dx1
    8 continue                                                          
c
c   st(4 = S12    at  i,1,k+1/2
c   st(6 = S32    at  i+1/2,1,k
c   
      jc=1
      do 5 kc=1,n3m
      do 5 ic=1,n1m                                                     
      st(ic,jc,kc,4)=q(1,ic,jc,kc)*dx2/cac(jc)*2.
      st(ic,jc,kc,6)=q(3,ic,jc,kc)*dx2/cac(jc)*2.
    5 continue                             
c
c   st(4 = S12    at  i,n2,k+1/2
c   st(6 = S32    at  i+1/2,n2,k
c   
      jc=n2                                                     
      jm=n2m                                                     
      do 15 kc=1,n3m
      do 15 ic=1,n1m             
      st(ic,jc,kc,4)=-q(1,ic,jm,kc)*dx2/cac(jc)*2.
      st(ic,jc,kc,6)=-q(3,ic,jm,kc)*dx2/cac(jc)*2.
   15 continue                             
c
c   evaluate s12 at the cell center
c
      do 20 kc=1,n3m
      do 20 jc=1,n2m   
      jp=jc+1                                                           
      do 20 ic=1,n1m             
      ip=ipv(ic)                                                        
      rhs(ic,jc,kc)=(st(ic,jc,kc,4)+st(ic,jp,kc,4)+
     1               st(ip,jc,kc,4)+st(ip,jp,kc,4))*0.25
   20 continue                             
      do 21 kc=1,n3m
      do 21 jc=1,n2m   
      do 21 ic=1,n1m             
      st(ic,jc,kc,4)=.5*rhs(ic,jc,kc)           
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
      st(ic,jc,kc,5)=.5*rhs(ic,jc,kc)           
   35 continue                             
c
c   evaluate s23 at the cell center
c
      do 50 kc=1,n3m
      kp=kpv(kc)                                                    
      do 50 jc=1,n2m   
      jp=jc+1                                                           
      do 50 ic=1,n1m             
      rhs(ic,jc,kc)=(st(ic,jc,kc,6)+st(ic,jp,kc,6)+
     1               st(ic,jc,kp,6)+st(ic,jp,kp,6))*0.25
   50 continue                             
      do 57 kc=1,n3m
      do 57 jc=1,n2m   
      do 57 ic=1,n1m             
      st(ic,jc,kc,6)=.5*rhs(ic,jc,kc)           
   57 continue                             
      return                                                            
      end                                                               

