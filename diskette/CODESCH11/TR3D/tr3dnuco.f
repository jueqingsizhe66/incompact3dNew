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
      pi=2.*asin(1.)                                                    
        if(irid.eq.0) then
      dx1=2.*pi/float(n1m)
                      else
      dx1=2.*pi/(float(n1m)*float(lamb))
                      endif
      dx2=1./float(n2m)                                                 
      dx3=1./float(n3m)                                               
      write(6,100)dx1,dx2,dx3                           
      write(32,100)dx1,dx2,dx3                           
 100  format(3x,'mesh size: d_theta= ',e9.4,' d_r= ',e9.4,' dz= ',e9.4)
      write(6,*)'  '
      write(32,*) '  '
      dx1=1./dx1                                                        
      dx2=1./dx2                                                        
      dx3=1./dx3                                                        
      dx1q=dx1*dx1                                                      
      dx2q=dx2*dx2                                                      
      dx3q=dx3*dx3                                                      
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
c     non-uniform grids in r and in z
c
      include 'param.f'
      dimension etaz(m3)
      dimension eta(2*m2)
      dimension xt1(2*m2),xt2(2*m2)
      dimension xt31(m3),xt32(m3)
c
      open(79,file='radstr.out')
      open(78,file='axistr.out')
      pi=2.*asin(1.)                                                    
      if (istr.lt.0) then
      do 111 j=1,n2
       x2=(j-1)/dx2
       eta(j)=x2
       rc(j)=eta(j)*rext
  111 continue
       endif
      if (istr.eq.0) then
      tstr2=tanh(strr)
      do 112 j=1,n2
       x2=(j-1)/dx2
       eta(j)=tanh(strr*x2)/tstr2
       rc(j)=eta(j)*rext
  112 continue
                     endif
      if (istr.eq.1) then
      n2tm=n2t-1
      n2vm=n2v-1
      dx2t=n2tm
      if(n2t.ne.n2v) then
      etdp=float(n2vm)/float(n2tm)
                    endif
      tstr=tanh(strr*etdp)
      do 113 j=1,n2t
       x2=(j-1)/dx2t
      xt1(j)=rmed1/rext*tanh(strr*x2)/tstr
  113 continue
      do 114 j=1,n2t
       x2=(j-1)/dx2t
      xt2(j)=1./xt1(n2t)+(1.-1./xt1(n2t))
     1     *tanh(strb*(x2-1.))/tanh(strb*(etdp-1.))
       eta(j)= xt1(j)*xt2(j)  
       rc(j)=eta(j)*rext
  114 continue
                      endif
      if (istr.eq.2) then
      dx2t=float(n2m)
      tstr=tanh(strr*etdp)
      do j=1,n2
       x2=(j-1)/dx2t
      xt1(j)=rmed1/rext*tanh(strr*x2)/tstr
      enddo
      do j=1,n2
       x2=(j-1)/dx2t
      xt2(j)=1./xt1(n2)+(1.-1./xt1(n2))
     1     *tanh(strb*(x2-1.))/tanh(strb*(etdp-1.))
       eta(j)= xt1(j)*xt2(j)
       rc(j)=eta(j)*rext
      enddo
                      endif
      do 115 j=1,n2
       x2=(j-1)/dx2
       yd=1.-rc(j)
      write(79,201) x2,rc(j)
  201 format(3x,5e12.5)
115   continue
      do 11 j=1,n2                                                      
      ragc(j)=rc(j)                                                        
   11 continue                                                          
      rc(1)=0.                                                          
c$$$$$$$$$$$$$Computation of geometry terms for non-uniform grids$$$$$$$$
      do 12 j=1,n2m                                                     
      ragm(j)=(ragc(j)+ragc(j+1))*0.5 
      rm(j)=ragm(j)                                                          
      g2rm(j)=(rc(j+1)-rc(j))*dx2
   12 continue                                                          
      do 122 j=2,n2m
      g2rc(j)=(rc(j+1)-rc(j-1))*dx2*0.5
122   continue
      g2rc(1)=(rc(2)-rc(1))*dx2
      g2rc(n2)= (rc(n2)-rc(n2m))*dx2
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
c     do jc=1,n2
c     write(78,201) rc(jc),rm(jc),g2rc(jc),g2rm(jc)
c     enddo
      do 13 i=1,n1                                                      
      thetac(i)=-(i-1)/dx1                                              
   13 continue                                                          
      do 14 i=1,n1m                                                     
      thetam(i)=-(i-1+0.5)/dx1                                          
   14 continue                                                          
      if(istr3.eq.0) then
      tstr3=tanh(str3)
      do k=1,n3
      x3=float(k-1)/float(n3m)
      etaz(k)=alx3*(1.+tanh(str3*(x3-1.0))/tstr3)
      zz(k)=etaz(k)                              
      write(78,201) x3,zz(k) 
      enddo
                     endif
      if(istr3.eq.2) then
      tstr3=tanh(str3*0.5)
      do k=1,n3
      x3=float(k-1)/float(n3m)
      etaz(k)=alx3*0.5*(1.+tanh(str3*(x3-0.5))/tstr3)
      zz(k)=etaz(k)                              
      write(78,201) x3,zz(k) 
      enddo
                     endif
      if (istr3.lt.0) then
      do k=1,n3
      x3=float(k-1)/float(n3m)
      etaz(k)=alx3*x3
      zz(k)=etaz(k)
      write(78,201) x3,zz(k)
      enddo
                      endif
      if (istr3.eq.1) then
      dx3t=n3m
      tstr3=tanh(str3*etdp3)
      do k=1,n3
       x3=float(k-1)/dx3t
      xt31(k)=rmed31/alx3*tanh(str3*x3)/tstr3
      end do
      do k=1,n3
       x3=(k-1)/dx3t
      xt32(k)=1./xt31(n3)+(1.-1./xt31(n3))
     1     *tanh(strb3*(x3-1.))/tanh(strb3*(etdp3-1.))
      etaz(k)= xt31(k)*xt32(k)
       zz(k)=etaz(k)*alx3
      write(78,201) x3,zz(k)
       end do
                 end if
      do k=1,n3m
      zzm(k)=(zz(k)+zz(k+1))*0.5
c$$$$$$$$$$$$$Computation of metric for non-uniform grids in x3$$

      g3rm(k)=(zz(k+1)-zz(k))*dx3
      enddo
      do k=2,n3m
      g3rc(k)=(zz(k+1)-zz(k-1))*dx3*0.5
      enddo
      g3rc(1)=(zz(2)-zz(1))*dx3
      g3rc(n3)= (zz(n3)-zz(n3m))*dx3
c
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

      close(79)
      close(78)
      return                                                            
      end                                                               
