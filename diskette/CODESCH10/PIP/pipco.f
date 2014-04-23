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
      dx1=2.*pi/float(n1m)
      dx2=1./float(n2m)                                                 
      dx3=alx3/float(n3m)                                               
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
c
      include 'param.f'
      dimension eta(2*m2)
      dimension xt1(2*m2),xt2(2*m2)
      common/npjet/n2t,n2v
c
      open(77,file='radjet.out')
      open(79,file='radstr.out')
      open(78,file='rg2.out')
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
      dx2t=n2tm
      if(n2t.ne.n2) then
      etdp=float(n2m)/float(n2tm)
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
      write(77,203) j,rc(j)
  203 format(3x,i4,3x,5e12.5)
  114 continue
      endif
      if (istr.eq.2) then
      stra=strr
      etda=0.25
      etdb=1.-0.25
      tstra=tanh(stra*etda)
      tstrb=tanh(strb*etda)
      n2mh=n2m/2+1
      do j=1,n2mh
       x2=(j-1)/dx2
      rc(j)=(1.-rmed1)/2.*(1.-tanh(stra*(etda-x2))/tstra)
      write(77,203) j,rc(j)
      enddo
      do j=n2mh,n2
       x2=(j-1)/dx2
      rc(j)=(1.-rmed1)+rmed1*(1.-tanh(strb*(etdb-x2))/tstrb)/2.
      write(77,203) j,rc(j)
      enddo
                     endif
      do 115 j=1,n2
       x2=(j-1)/dx2
       yd=1.-rc(j)
      write(79,201) x2,yd,rc(j)
  201 format(3x,5e12.5)
115   continue
      do 11 j=1,n2                                                      
      ragc(j)=rc(j)                                                        
   11 continue                                                          
      rc(1)=0.                                                          
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
      do jc=1,n2
      write(78,201) rc(jc),rm(jc),g2rc(jc),g2rm(jc)
      enddo
      do 13 i=1,n1                                                      
      thetac(i)=-(i-1)/dx1                                              
   13 continue                                                          
      do 14 i=1,n1m                                                     
      thetam(i)=-(i-1+0.5)/dx1                                          
   14 continue                                                          
      do 2 k=1,n3                                                       
      zz(k)=(k-1)/dx3                                                   
   2  continue                                                          
      close(79)
      close(78)
      return                                                            
      end                                                               
c************************************************************************
c                                                                       *
c   ************** subroutine cordino                                   *
c                                                                       *
c************************************************************************
      subroutine cordino(n2o)
c
c     Physical coordinates of a previous simulation
c     This coordinate is necessary to interpolate from
c     the old grid to the new grid and to perform
c     the simulation with the new radial distribution
c     of grid points calculated in the cordin routine
c
      include 'param.f'
      dimension eta(2*m2)
      dimension xt1(2*m2),xt2(2*m2)
c
      n2mo=n2o-1
      dx2o=float(n2mo)
      open(77,file='rajeto.out')
      open(79,file='rastro.out')
      if (istro.lt.0) then
      do 111 j=1,n2o
       x2=(j-1)/dx2o
       eta(j)=x2
       rco(j)=eta(j)*rexto
  111 continue
       endif
      if (istro.eq.0) then
      tstr2=tanh(strro)
      do 112 j=1,n2o
       x2=(j-1)/dx2o
       eta(j)=tanh(strro*x2)/tstr2
       rco(j)=eta(j)*rexto
  112 continue
       endif
      if (istro.eq.1) then
      n2tm=n2to-1
      dx2t=n2tm
      if(n2to.ne.n2o) then
      etdp=float(n2mo)/float(n2tm)
                    endif
      tstr=tanh(strro*etdp)
      do 113 j=1,n2to
       x2=(j-1)/dx2t
      xt1(j)=rmed1o/rexto*tanh(strro*x2)/tstr
  113 continue
      do 114 j=1,n2to
       x2=(j-1)/dx2t
      xt2(j)=1./xt1(n2to)+(1.-1./xt1(n2to))
     1     *tanh(strbo*(x2-1.))/tanh(strbo*(etdp-1.))
       eta(j)= xt1(j)*xt2(j)
       rco(j)=eta(j)*rext
      write(77,203) j,rco(j)
  203 format(3x,i4,3x,5e12.5)
  114 continue
      endif
      do 115 j=1,n2o
       x2=(j-1)/dx2o
       yd=1.-rco(j)
      write(79,201) x2,yd,rco(j)
  201 format(3x,5e12.5)
115   continue
      rco(1)=0.
      do 12 j=1,n2mo
      rmo(j)=(rco(j)+rco(j+1))*0.5

c$$$$$$$$$$$$$Computation of geometry terms for non-uniform grids$$$$$$$$

      g2rmo(j)=(rco(j+1)-rco(j))*dx2o

c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

   12 continue

      do 122 j=2,n2mo
c$$$$$$$$$$$$$Computation of geometry terms for non-uniform grids$$$$$$$$

      g2rco(j)=(rco(j+1)-rco(j-1))*dx2o*0.5
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
122   continue
      g2rco(1)=(rco(2)-rco(1))*dx2o
      g2rco(n2o)= (rco(n2o)-rco(n2mo))*dx2o
      close(79)
      close(78)
      return
      end
