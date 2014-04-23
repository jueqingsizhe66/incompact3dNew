      program trasc                                                    
c                                                                       
c     program to verify non-uniforms  coordinates                
c                                                                       
      include 'param.f'
      character *60 namfile
  301 format(a4)
      open(15,file='trascor.d')
      read(15,301) dummy
      read(15,*)n1
      n1m=n1-1
      dx=1./float(n1-1)                                                 
      dx=1./dx                              
      dxq=dx*dx
c                                                                       
c                                                                       
      call coordi 
      call metric                               
      end                                                               
      function dgal(al,zitaf)
c
c   analytical derivatives of the function f 
c   to compare the discrete derivatives
c
      alzitf=al*zitaf
      dg1=1./(sinh(alzitf)*cosh(alzitf))
      dg2=-alzitf/sinh(alzitf)**2
      dg3=-alzitf/cosh(alzitf)**2
      dgal=dg1+dg2+dg3
      return
      end
c
c
      function gal(al,zitaf)
c
c   function f for which the derivatives are evaluated
c
      alzitf=al*zitaf
      gal=al/(cosh(alzitf)*sinh(alzitf))
      return
      end

c                                                                       
c                                                                       
c     **********************************************************        
      subroutine coordi
c                                                                       
c     In this routine the grid in the physical space are
c     calculated and these are functions of the grid points
c     in the computational space.
c                                                                       
      include 'param.f'
      dimension xt1(m1),xt2(m1)
      character *60 namfile
      pi=2.*asin(1.)
  110 format(5e12.5)
c                                                                       
c     Transformation 2.4.a                                              
c                                                                       
      namfile='coord1a.plo'
      open(18,file=namfile)
  301 format(a4)
      read(15,301) dummy
      read(15,*) al
      gam=tanh(al*0.5)                                                      
      do i=1,n1                                                      
      csi=(i-1)/float(n1m)                                   
      eta=al*(csi-0.5)
      xt1(i)=0.5*(1.+tanh(eta)/gam)
      y(i)=(-0.5+xt1(i))*2.
      write(18,110)csi,y(i)
      enddo
c                                                                       
c     Transformation 2.4.b                                              
c                                                                       
      namfile='coord1b.plo'
      open(18,file=namfile)
      read(15,301) dummy
      read(15,*)  nc,xc,xf,alin
      csic=(nc-1)/float(n1m)*0.5+0.5
      etac=abs(csic-0.5)
      etaf=0.5
      zitaf=etaf-etac
      al=alin
      beta=xc/etac
      iter=1
      val=beta/(xf-xc)  
c     write(6,*)val,zitaf
   33 continue
      res=gal(al,zitaf)-val
      if(abs(res).lt..1e-04) then
      go to 35
                             else
      tang=dgal(al,zitaf)
      dal=res/tang
      write(6,*)iter,res,al,tang,dal
      al=al-dal
      iter=iter+1
      if(iter.gt.21) go to 35
      go to 33
                             endif
   35 continue
      write(6,*)'convergence',iter,res,al
      gam=tanh(al*zitaf)                                                      
      do i=1,n1                                                      
      csi=(i-1)/float(n1m)                                   
      eta=csi-0.5
      if(abs(eta).le.etac)  then
      y(i)=eta*xc/etac
                            else 
      siet=eta/abs(eta) 
      if(siet.ge.0) then
      zita=(eta-etac)
      azit=al*(zita-zitaf)
      y(i)=(xf+tanh(azit)/gam*(xf-xc))
                            else 
      zita=(eta+etac)
      azit=al*(zita+zitaf)
      y(i)=(-xf+tanh(azit)/gam*(xf-xc))
                            endif                                   
                            endif                                   
      write(18,110)csi,y(i)
      enddo
c                                                                       
c     Transformation 2.4.c                                              
c                                                                       
      namfile='coord1c.plo'
      open(18,file=namfile)
      read(15,301) dummy
      read(15,*)  al      
      gam=tanh(al)                                                      
      alet=2.*al
      do i=1,n1                                                      
      csi=(i-1)/float(n1m)                                   
      if(csi.le.0.5) then
      eta=alet*csi
      y(i)=-1+tanh(eta)/gam                                   
                     else
      eta=alet*(csi-1.)
      y(i)=1.+tanh(eta)/gam                                   
                     endif
      write(18,110)csi,y(i)
      enddo
c                                                                       
c                                                                       
c     Transformation 2.4.d                                              
c                                                                       
      namfile='coord1d.plo'
      open(18,file=namfile)
      read(15,301) dummy
      read(15,*)  al      
      gam=atanh(al*0.25)                                                      
      do i=1,n1                                                      
      csi=(i-1)/float(n1m)                                   
      if(csi.le.0.5) then
      eta=(0.25-csi)
      alet=al*eta
      y(i)=0.5*(-1-atanh(alet)/gam)                                  
                     else
      eta=(csi-0.75)
      alet=al*eta
      y(i)=0.5*(+1+atanh(alet)/gam)                                  
                     endif
      write(18,110)csi,y(i)
      enddo
c                                                                       
c                                                                       
c     Transformation 2.4.e                                              
c                                                                       
      namfile='coord1e.plo'
      open(18,file=namfile)
      read(15,301) dummy
      read(15,*)  al,hbig,hsma
      gam=tanh(al*0.25)                                                      
      do i=1,n1                                                      
      csi=(i-1)/float(n1m)                                   
      if(csi.le.0.5) then
      eta=(0.25-csi)
      alet=al*eta
      y(i)=0.5*(+1.-tanh(alet)/gam)*(hbig-hsma)
                     else
      eta=(csi-0.75)
      alet=al*eta
      y(i)=hbig-hsma+hsma*(1.+tanh(alet)/gam)                                  
                     endif
      write(18,110)csi,y(i)
      enddo
c                                                                       
c     Transformation 2.4.f                                              
c                                                                       
      namfile='coord1f.plo'
      open(18,file=namfile)
      read(15,301) dummy
      read(15,*)  al1,al2,xc,xf,nc
      csic=(nc-1)/float(n1m)
      gam=tanh(al1*csic)                                                      
      do i=1,n1
      csi=(i-1)/float(n1m)
      xt1(i)=xc/xf*tanh(al1*csi)/gam
      enddo
      write(6,*)'2.4.f  ',csic,xt1(n1),xt1(1)
      do i=1,n1
      csi=(i-1)/float(n1m)
      xt2(i)=1./xt1(n1)+(1.-1./xt1(n1))
     1     *tanh(al2*(csi-1.))/tanh(al2*(csic-1.))
      y(i)=xt1(i)*xt2(i)*xf
      write(18,110)csi,y(i)
      enddo
      write(6,*)'2.4.f  ',xt2(1),xt2(n1)
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
