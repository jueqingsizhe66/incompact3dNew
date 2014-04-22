      program scdf                                                    
c                                                                       
c     program to evaluate the accuracy of the second
c     derivative when non-uniform grids are used
c     IN this code a function appropriate for a mixing
c     layer is used
c                                                                       
      include 'param.f'
      character *60 namfile
      write(6,*)'enter alfu, al, yprer'
      read(5,*)alfu,al,yprer
      namfile='accuracy.plo'
      open(20,file=namfile)
  191 continue
      write(6,*)' n1'
      read(5,*)n1
      n1m=n1-1
      dx=1./float(n1-1)                                                 
      dx=1./dx                              
      write(*,*) dx
      dxq=dx*dx
      write(*,*) dx,dxq
      write(6,*)'dx,dxq',dx,dxq
c                                                                       
      call coordi 
      write(6,*)' coordi evaluated'
      call metric                               
      write(6,*)' metric evaluated'
      call funct 
c                                                                       
c     evaluation of the discretized second derivative
c                                                                       
      call calc
c                                                                       
c                                                                       
      call outp
      write(6,*)' enter 1 for a further case'
      read(5,*)imor
      if(imor.eq.1) go to 191
      stop                                                              
      end                                                               
c                                                                       
c                                                                       
c     ***************************************************************   
      subroutine outp
c                                                                       
      include 'param.f'
      character *4 npoi
      character *60 namfile
c                                                                       
  133 format(i4.4)
      write(npoi,133) n1m
      an1=n1m
      namfile='errdis'//npoi//'.plo'
      open(16,file=namfile)
c
c   write the files with the results
c
      ! add
      write(16,*) 'csi   er1    er2   er3'
      write(17,*) 'csi   sd1(i)   sd2(i)   sd3(i)'
      write(20,*) 'an1   er1    er2   er3'
      do  i=2,n1m
      csi=(i-1)/float(n1m)
      if(abs(sda(i)).gt.0.) then
      er1=abs((sd1(i)-sda(i))/sda(i))
      er2=abs((sd2(i)-sda(i))/sda(i))
      er3=abs((sd3(i)-sda(i))/sda(i))
                            else
      er1=0.
      er2=0.
      er3=0.
                            endif
      if(y(i).ge.0.) then
      write(16,794)csi,er1,er2,er3                                      
      write(17,794)csi,sd1(i),sd2(i),sda(i)                             
                     endif
      if(y(i).ge.yprer.and.y(i-1).le.yprer) then
      write(20,794) an1,er1,er2,er3
                                            endif
 794  format(4x,6e12.5)                                  
      enddo
      return                                                            
      end                                                               
c                                                                       
c                                                                       
c                                                                       
c     **********************************************************        
      subroutine funct
      include 'param.f'
      character *60 namfile
      namfile='funct.plo'
      open(18,file=namfile)
c                                                                       
c     here the function of which the derivatives are 
c     calculated is given together with its first and
c     second derivatives
c     Here the case of an hiperbolic tangent is
c     assigned typical of a mixing layer
c                                                                       
      n1mh=n1m/2+1
      gam=tanh(alfu)
      write(18,*) 'y(i), u(i)'   !add
      do i=1,n1 
      if(y(i).lt.0.) then  ! y > 0  y 距离中线的距离 
      eta=y(i)+1.  ! y+1
      u(i)=-1.+tanh(eta*alfu)/gam
      dpu=alfu/(gam*cosh(alfu*eta)**2)
      dsu=-alfu*dpu*2.*tanh(eta*alfu)
                      else
      eta=y(i)-1.  ! y-1
      u(i)=1.+tanh(eta*alfu)/gam
      dpu=alfu/(gam*cosh(alfu*eta)**2)
      dsu=-alfu*dpu*2.*tanh(eta*alfu)
                      endif
      if(i.eq.n1mh) then
      u(i)=0.
      dsu=0.
                    endif
      sda(i)=dsu
      if(y(i).ge.0.) then
      write(18,110)y(i),u(i)
                     endif
  110 format(5e12.5)
      enddo
      return                                                            
      end                                                               
c                                                                       
c                                                                       
c     **********************************************************        
      subroutine coordi
      include 'param.f'
      pi=2.*asin(1.)
c                                                                       
c     here the coordinate transformation is given
c     has been chosen the case where clustering is 
c     near the value of Y((n1-1)/2+1)=0  
c                                                                       
      !clustering in the center   here y is similar to x
      !dmay   dcay 只用于sd3可以暂时不用考虑
      gam=tanh(al)                                                      
      ! 这边可以说明！计算域网格是i  物理域网格是y(i),对应
      ! 书本上13页的x 
      alet=2.*al
      write(9,*) 'csi  y(i)'   ! add
      do i=1,n1                                                      
      csi=(i-1)/float(n1m)                                   
      if(csi.le.0.5) then
      eta=alet*csi
      y(i)=-1+tanh(eta)/gam                                   
      dcay(i)=alet/(gam*cosh(eta)**2)
                     else
      eta=alet*(csi-1.)
      y(i)=1.+tanh(eta)/gam                                   
      dcay(i)=alet/(gam*cosh(eta)**2)
                     endif
      write(9,110)csi,y(i)
  110 format(5e12.5)
      enddo
c                                                                       
      do i=1,n1m                                                     
      csi=(i-1+0.5)/float(n1m)                                   
      if(csi.le.0.5) then
      eta=alet*csi
      dmay(i)=alet/(gam*cosh(eta)**2)
                     else
      eta=alet*(csi-1.)
      dmay(i)=alet/(gam*cosh(eta)**2)
                     endif
      enddo
c                                                                       
      return                                                            
      end                                                               
c                                                                       
c     ***************************************************************   
c                                                                       
      subroutine metric
c
c      In this routine the first derivative of the physical coordinates
c      with respect to the computational grid are evaluated
c      at the grid points dcy and at the mid points dmy
c      ddcy is the second derivative of the physical coordinate
c
      include 'param.f'
      do i=2,n1m                                                    
      dcy(i)=(y(i+1)-y(i-1))*0.5*dx              
      ddcy(i)=(y(i+1)-2.*y(i)+y(i-1))*dxq              
      enddo
      do i=1,n1m                                                     
      dmy(i)=(y(i+1)-y(i))*dx                                           
      enddo
      return                                                            
      end                                                               
c                                                                       
c                                                                       
c     ************************************************************      
c                                                                       
      subroutine calc
      include 'param.f'
      write(64,*) '64 where',dx,dxq
      write(63,*) 'sd1(i)   sd2(i)   sda(i)'  !add
      do i=2,n1m                                                    
c                                                                       
c     evaluation of the two expressions of second derivates
c                                                                       
      write(64,*)i,dmy(i),dcy(i),ddcy(i),u(i)
      aa1=dxq/(dcy(i)*dmy(i))                                         
      bb1=-dxq*(1./dmy(i)+1./dmy(i-1))/(dcy(i))                          
      cc1=dxq/(dcy(i)*dmy(i-1))                                       
      sd1(i)=(aa1*u(i+1)+bb1*u(i)+cc1*u(i-1))
      aa2=dxq/dcy(i)**2-dx*ddcy(i)/dcy(i)**3
      bb2=-2.*dxq/dcy(i)**2
      cc2=dxq/dcy(i)**2+dx*ddcy(i)/dcy(i)**3
      sd2(i)=(aa2*u(i+1)+bb2*u(i)+cc2*u(i-1))
      aa3=dxq/(dcay(i)*dmay(i))                                         
      bb3=-dxq*(1./dmay(i)+1./dmay(i-1))/(dcay(i))                          
      cc3=dxq/(dcay(i)*dmay(i-1))                                       
      sd3(i)=(aa3*u(i+1)+bb3*u(i)+cc3*u(i-1))
      write(63,*)i,sd1(i),sd2(i),sda(i)
      enddo
c                                                                       
c     the solutions are written in a file together
c     with the exact value
c                                                                       
      return                                                            
      end                                                               
