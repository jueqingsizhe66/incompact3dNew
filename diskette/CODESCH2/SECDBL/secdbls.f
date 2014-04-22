      program scdf                                                    
c
c     program to evaluate the accuracy of the second
c     derivative when non-uniform grids are used
c     IN this code a function appropriate for a boundary
c     layer is used
c
      include 'param.f'
      character *60 namfile
      write(6,*)'enter  al, yprer'
      read(5,*)al,yprer
      namfile='accuracy.plo'
      open(20,file=namfile)
      namfile='accurwal.plo'
      open(21,file=namfile)
      namfile='accurupp.plo'
      open(22,file=namfile)
  191 continue
      write(6,*)' n1'
      read(5,*)n1
      n1m=n1-1
      dx=1./float(n1-1)                                                 
      dx=1./dx                              
      dxq=dx*dx
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
c
c   write the files with the results
c

      open(16,file=namfile)
      do  i=1,n1m
      csi=(i-1+0.5)/float(n1m)
      if(abs(sda(i)).gt.0.) then
      er1=abs((sd1(i)-sda(i))/sda(i))
      er3=abs((sd3(i)-sda(i))/sda(i))
                            else
      er1=0.
      er3=0.
                            endif
      if(y(i).ge.0.) then
      write(16,794)csi,er1,er3                                      
      write(17,794)ys(i),sd1(i),sd3(i),sda(i)                             
                     endif
      if(y(i).ge.yprer.and.y(i-1).le.yprer) then
      write(20,794) an1,er1,er3
                                            endif
      if(i.eq.1) then
      write(21,794) an1,er1,er3
                                            endif
      if(i.eq.n1m) then
      write(22,794) an1,er1,er3
                                            endif
 794  format(4x,6e12.5)                                  
      enddo
      return                                                            
      end                                                               
c                                                                       
c                                                                       
c     **********************************************************        
      subroutine funct
      include 'param.f'
c                                                                       
c     here the function of which the derivatives are
c     calculated is given together with its first and
c     second derivatives
c     Here the case of a turbulent velocity profile
c     in a boundary layer is assigned
c
      do i=1,n1m 
      eta=ys(i)
      u(i)=eta**(1./7.)                                   
      sda(i)=1./7.*(-1.+1./7.)*eta**(-2.+1./7.)
      write(10,110)y(i),u(i)
  110 format(5e12.5)
      enddo
      uup=1.
      uwall=0.
      return                                                            
      end                                                               
c                                                                       
c                                                                       
c     **********************************************************        
      subroutine coordi
      include 'param.f'
      pi=2.*asin(1.)
c                                                                       
c
c     here the coordinate transformation is given
c     has been chosen the case where clustering is
c     near the value of Y(1)=0. (representative of 
c     a solid wall.
c
      gam=tanh(al)                                                      
      do i=1,n1                                                      
      csi=(i-1)/float(n1m)                                   
      eta=al*(csi-1.)
      y(i)=(1+tanh(eta)/gam)
      dcay(i)=al/(gam*cosh(eta)**2)
      write(9,110)csi,y(i)
  110 format(5e12.5)
      enddo
c                                                                       
      do i=1,n1m                                                     
      csi=(i-1+0.5)/float(n1m)                                   
      eta=al*(csi-1.)
      ys(i)=(1+tanh(eta)/gam)
      dmay(i)=al/(gam*cosh(eta)**2)
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
      dcy(i)=(ys(i)-ys(i-1))*dx              
c     dcy(i)=(y(i+1)-y(i-1))*dx*0.5
      enddo
      dcy(1)=(ys(2)-y(1))*dx*2
      dcy(n1)=(y(n1)-ys(n1m))*dx*2
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
      write(64,*)dx,dxq
      do i=2,n1m-1                                                  
c
c     evaluation of the two expressions of second derivates
c
      aaa=dxq/dmy(i)                                         
      aa1=aaa/dcy(i+1)
      bb1=-aaa*(1./dcy(i+1)+1./dcy(i))
      cc1=aaa/dcy(i)
      sd1(i)=(aa1*u(i+1)+bb1*u(i)+cc1*u(i-1))
      aaa=dxq/dmay(i)                                         
      aa3=aaa/dcay(i+1)
      bb3=-aaa*(1./dcay(i+1)+1./dcay(i))
      cc3=aaa/dcay(i)
      sd3(i)=(aa3*u(i+1)+bb3*u(i)+cc3*u(i-1))
      enddo
      aaa=2.*dxq/(dcy(2)+dcy(1)*0.5)
      aa1=aaa/dcy(2)                                      
      bb1=-aaa*(1./dcy(2)+2./dcy(1))
      cc1=aaa*2./dcy(1)
      sd1(1)=(aa1*u(2)+bb1*u(1)+cc1*uwall)
      aaa=2.*dxq/(dcay(2)+dcay(1)*0.5)
      aa3=aaa/dcay(2)                                      
      bb3=-aaa*(1./dcay(2)+2./dcay(1))
      cc3=aaa*2./dcay(1)
      sd3(1)=(aa3*u(2)+bb3*u(1)+cc3*uwall)
      aaa=2.*dxq/(dcy(n1m)+dcy(n1)*0.5)
      cc1=aaa/dcy(n1m)                                      
      bb1=-aaa*(1./dcy(n1m)+2./dcy(n1))
      aa1=aaa*2./dcy(n1)
      sd1(n1m)=(aa1*uup+bb1*u(n1m)+cc1*u(n1m-1))
      aaa=2.*dxq/(dcay(n1m)+dcay(n1)*0.5)
      cc3=aaa/dcay(n1m)                                      
      bb3=-aaa*(1./dcay(n1m)+2./dcay(n1))
      aa3=aaa*2./dcay(n1)
      sd3(n1m)=(aa3*uup+bb3*u(n1m)+cc3*u(n1m-1))
      return                                                            
      end                                                               
