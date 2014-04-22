c**********************************************************************
C   TRIDIAGONAL SOLVERS FOR PERIODIC CONDITIONS
c   DEPEND IN WHICH DIRECTION ARE DONE THE INNER LOOPS
c  ****************************** subrout trvpjk  **********************
c
c
      subroutine trvpjk(j1,j2,mi,mf )
      include 'param.f'
c
      common/tvpjk/amk(m2,m3),ack(m2,m3),apk(m2,m3),fk(m2,m3)
     1         ,qek(m2,m3),qk(m2,m3),sk(m2,m3)
      common/util1/fnn(m2),p(m2),bet(m1),u(m1,m2),x(m1),fn(m2),gmm(m2),
     1 gm(m1,m2)
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
      common/util/fnn(m2),p(m2),bet(m1),u(m1,m2),x(m1),fn(m2),gmm(m2),
     1 gm(m1,m2)
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
      p(k) =1./( aci(k,j) + ami(k,j)*qq(k,j-1))                              
      qq(k,j) = - api(k,j)*p(k)                                            
      ss(k,j) = - ami(k,j)*ss(k,j-1)*p(k)                                   
      fi(k,j) = ( fi(k,j) - ami(k,j)*fi(k,j-1))*p(k)                         
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
c************************************************************************
c                                                                       *
c  ****************************** subrout trvpkj  ***********************
c                                                                       *
c************************************************************************
      subroutine trvpkj(ji,jf,ni,nf)                                    
      include 'param.f'
      common/tvpkj/amj(m3,m2),acj(m3,m2),apj(m3,m2),
     1           fj(m3,m2),fej(m3,m2),qj(m3,m2),sj(m3,m2)
      dimension fnn(m3)
c                                                                       
c     vectorized for right hand side and coefficients                   
c                                                                       
      ja = ji + 1                                                       
      jj = ji + jf                                                      
      do 20 k=ni,nf                                                     
      qj(k,ji) = -apj(k,ji)/acj(k,ji)                                        
      sj(k,ji) = -amj(k,ji)/acj(k,ji)                                        
      fnn(k) = fj(k,jf)                                                   
      fj(k,ji) = fj(k,ji)/acj(k,ji)                                         
   20 continue                                                          
c                                                                       
c     forward elimination sweep                                         
c                                                                       
      do 10 j=ja,jf                                                     
      do 21 k=ni,nf                                                     
      p =1./( acj(k,j) + amj(k,j)*qj(k,j-1))                              
      qj(k,j) = - apj(k,j)*p                                            
      sj(k,j) = - amj(k,j)*sj(k,j-1)*p                                   
      fj(k,j) = ( fj(k,j) - amj(k,j)*fj(k,j-1))*p                         
   21 continue                                                          
   10 continue                                                          
c                                                                       
c     backward pass                                                     
c                                                                       
      do 22 k=ni,nf                                                     
      sj(k,jf) = 1.                                                      
      fej(k,jf) = 0.                                                     
   22 continue                                                          
      do 11 i=ja,jf                                                     
      j = jj - i                                                        
      do 23 k=ni,nf                                                     
      sj(k,j) = sj(k,j) + qj(k,j)*sj(k,j+1)                                 
      fej(k,j) = fj(k,j) + qj(k,j)*fej(k,j+1)                               
   23 continue                                                          
   11 continue                                                          
      do 24 k=ni,nf                                                     
      fj(k,jf)=(fnn(k) - apj(k,jf)*fej(k,ji) - amj(k,jf)*fej(k,jf-1))           
     &       /(apj(k,jf)*sj(k,ji) + amj(k,jf)*sj(k,jf-1)  +acj(k,jf))           
   24 continue                                                          
c                                                                       
c     backward elimination pass                                         
c                                                                       
      do 12 i=ja,jf                                                     
      j = jj -i                                                         
      do 25 k=ni,nf                                                     
      fj(k,j) = fj(k,jf)*sj(k,j) + fej(k,j)                                 
   25 continue                                                          
   12 continue                                                          
      return                                                            
      end                                                               
c************************************************************************
c                                                                       *
c  ****************************** subrout tripvv  ***********************
c                                                                       *
c************************************************************************
      subroutine tripvv(ami,aci,api,fi,ji,jf,ni,nf)
      include 'param.f'
      dimension ami(m1,m2),aci(m1,m2),api(m1,m2),fi(m1,m2)
      dimension fei(m1,m2),qq(m1,m2),ss(m1,m2)
      dimension fnn(m1),p(m1)
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
      p(k) =1./( aci(k,j) + ami(k,j)*qq(k,j-1))                              
      qq(k,j) = - api(k,j)*p(k)                                            
      ss(k,j) = - ami(k,j)*ss(k,j-1)*p(k)                                   
      fi(k,j) = ( fi(k,j) - ami(k,j)*fi(k,j-1))*p(k)                         
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
