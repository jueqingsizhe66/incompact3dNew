c***************************************
c    HERE THERE ARE A BUNCH OF TRIDIAGONAL SOLVERS
C    WITH INNER DO LOOPS
C    THESE SOLVERS ARE FOR PERIODIC OR NOT-PERIODIC
C    DIRECTIONS
C********************************************
caaaa************************************************************************
c                                                                       *
c  ****************************** subrout trvbkj  ********************** *
c                                                                       *
c************************************************************************
      subroutine trvbkj(n,mi,mf)
c
      include 'param.f'
      common/tvbkj/amj(m3,m2),acj(m3,m2),apj(m3,m2),fj(m3,m2)
      dimension btj(m3),gmkj(m3,m2)
      do 10 i=mi,mf
      btj(i)=acj(i,1)
      fj(i,1)=fj(i,1)/btj(i)
   10 continue
      do 11 j=2,n
      do 21 i=mi,mf
      gmkj(i,j)=apj(i,j-1)/btj(i)
      btj(i)=acj(i,j)-amj(i,j)*gmkj(i,j)
      fj(i,j)=(fj(i,j)-amj(i,j)*fj(i,j-1))/btj(i)
   21 continue
   11 continue
      do 12 j=n-1,1,-1
      do 22 i=mi,mf
      fj(i,j)=fj(i,j)-gmkj(i,j+1)*fj(i,j+1)
   22 continue
   12 continue
      return
      end
c  ****************************** subrout trvpjk  **********************
c
      subroutine trvpjk(j1,j2,mi,mf )
      include 'param.f'
c
      common/tvpjk/amk(m2,m3),ack(m2,m3),apk(m2,m3),fk(m2,m3)
     1         ,qek(m2,m3),qk(m2,m3),sk(m2,m3)
      dimension fn(m2)
      ja = j1 + 1
      jj = j1 + j2
      do 20 k=mi,mf
      qk(k,j1) = -apk(k,j1)/ack(k,j1)
      sk(k,j1) = - amk(k,j1)/ack(k,j1)
      fn(k) = fk(k,j2)
      fk(k,j1) = fk(k,j1)/ack(k,j1)
   20 continue
c
c     forward elimination sweep
c
      do 10 j=ja,j2
      do 21 k=mi,mf
      pp =1./( ack(k,j) + amk(k,j)*qk(k,j-1))
      qk(k,j) = - apk(k,j)*pp
      sk(k,j) = - amk(k,j)*sk(k,j-1)*pp
      fk(k,j) = ( fk(k,j) - amk(k,j)*fk(k,j-1))*pp
   21 continue
   10 continue
c
c     backward pass
c
      do 22 k=mi,mf
      sk(k,j2) = 1.
      qek(k,j2) = 0.
   22 continue
      do 11 i=ja,j2
      j = jj - i
      do 23 k=mi,mf
      sk(k,j) = sk(k,j) + qk(k,j)*sk(k,j+1)
      qek(k,j) = fk(k,j) + qk(k,j)*qek(k,j+1)
   23 continue
   11 continue
      do 24 k=mi,mf
      fk(k,j2)=(fn(k) - apk(k,j2)*qek(k,j1) - amk(k,j2)*qek(k,j2-1))
     &       /(apk(k,j2)*sk(k,j1) + amk(k,j2)*sk(k,j2-1)  +ack(k,j2))
   24 continue
c
c     backward elimination pass
c
      do 12 i=ja,j2
      j = jj -i
      do 25 k=mi,mf
      fk(k,j) = fk(k,j2)*sk(k,j) + qek(k,j)
   25 continue
   12 continue
      return
      end
c************************************************************************
c                                                                       *
c  ****************************** subrout tripvi  ***********************
c                                                                       *
c************************************************************************
      subroutine tripvi(ji,jf,ni,nf)                                    
      include 'param.f'
      common/trici/ami(m2,m1),aci(m2,m1),api(m2,m1),
     1           fi(m2,m1),fei(m2,m1),qq(m2,m1),ss(m2,m1)
      dimension fnn(m2)
c                                                                       
c     vectorized for right hand side and coefficients                   
c                                                                       
      ja = ji + 1                                                       
      jj = ji + jf                                                      
      do 20 k=ni,nf                                                     
      qq(k,ji) = -api(k,ji)/aci(k,ji)                                        
      ss(k,ji) = -ami(k,ji)/aci(k,ji)                                        
      fnn(k) = fi(k,jf)                                                   
      fi(k,ji) = fi(k,ji)/aci(k,ji)                                         
   20 continue                                                          
c                                                                       
c     forward elimination sweep                                         
c                                                                       
      do 10 j=ja,jf                                                     
      do 21 k=ni,nf                                                     
      pp =1./( aci(k,j) + ami(k,j)*qq(k,j-1))                              
      qq(k,j) = - api(k,j)*pp                                            
      ss(k,j) = - ami(k,j)*ss(k,j-1)*pp                                   
      fi(k,j) = ( fi(k,j) - ami(k,j)*fi(k,j-1))*pp                         
   21 continue                                                          
   10 continue                                                          
c                                                                       
c     backward pass                                                     
c                                                                       
      do 22 k=ni,nf                                                     
      ss(k,jf) = 1.                                                      
      fei(k,jf) = 0.                                                     
   22 continue                                                          
      do 11 i=ja,jf                                                     
      j = jj - i                                                        
      do 23 k=ni,nf                                                     
      ss(k,j) = ss(k,j) + qq(k,j)*ss(k,j+1)                                 
      fei(k,j) = fi(k,j) + qq(k,j)*fei(k,j+1)                               
   23 continue                                                          
   11 continue                                                          
      do 24 k=ni,nf                                                     
      fi(k,jf)=(fnn(k) - api(k,jf)*fei(k,ji) - ami(k,jf)*fei(k,jf-1))           
     &       /(api(k,jf)*ss(k,ji) + ami(k,jf)*ss(k,jf-1)  +aci(k,jf))           
   24 continue                                                          
c                                                                       
c     backward elimination pass                                         
c                                                                       
      do 12 i=ja,jf                                                     
      j = jj -i                                                         
      do 25 k=ni,nf                                                     
      fi(k,j) = fi(k,jf)*ss(k,j) + fei(k,j)                                 
   25 continue                                                          
   12 continue                                                          
      return                                                            
      end                                                               
