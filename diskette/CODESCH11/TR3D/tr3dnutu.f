c                                                                       
c  ************************* subrout veltur  **********************  
c                                                                       
c     this subroutine calculates the averaged velocity in r-z plane    
c     the averages is in the azimuthal direction
c                                                                       
      subroutine veltur                    
      include 'param.f'
      common/prcost/pcost
      avgn=1./(float(n1m))
      vl123=1./(rc(n2)*n1m*n2m*n3m)
      volto=1./(pi*alx3d)
      pcost=0.
      do jc=1,n2m
      do kc=1,n3m
      do ic=1,n1m
      pcost=pcost+pr(ic,jc,kc)
      enddo
      enddo
      enddo
      pcost=pcost/float(n1m*n2m*n3m)
c
c  velocities at the center of the cell
c
      do jc=1,n2m                                                     
      jm=jmv(jc)                                                        
      jp=jc+1
      do kc=1,n3m 
      do ic=1,n1m                                                     
      if(jc.eq.1) then
      q2akc= (q2(ic,jp,kc) - q2(isym(ic),jp,kc))*0.5/rc(jp)
      dph(ic,jc,kc)=(q2akc+q2(ic,jp,kc)/rc(jp))*.5
       else
      dph(ic,jc,kc)=
     1  (    q2(ic,jc,kc)/rc(jc)+q2(ic,jp,kc)/rc(jp))*0.5
      end if
      dq(ic,jc,kc)=
     1  (    q1(ic,jc,kc)+q1(ipv(ic),jc,kc) )*0.5
      qcap(ic,jc,kc)=(q3(ic,jc,kc)+q3(ic,jc,kc+1) )*0.5 
      enddo
      enddo
      enddo
      enej=0.
      enpsv=0.
c
c   mean velocity distribution in r-z planes
c
      do jc=1,n2m                                                     
      jm=jmv(jc)                                                        
      jp=jc+1
      do kc=1,n3m                                                     
      vmed(1,kc,jc)=0.
      vmed(2,kc,jc)=0.
      vmed(3,kc,jc)=0.
      pmed(kc,jc)=0.
      pscmed(kc,jc)=0.
      do ic=1,n1m 
      pscmed(kc,jc)=pscmed(kc,jc)+psc(ic,jc,kc)
      pmed(kc,jc)=pmed(kc,jc)+pr(ic,jc,kc)-pcost
      vmed(2,kc,jc)=vmed(2,kc,jc)+dph(ic,jc,kc)
      vmed(1,kc,jc)=vmed(1,kc,jc)+dq(ic,jc,kc)
      vmed(3,kc,jc)=vmed(3,kc,jc)+qcap(ic,jc,kc)
      enddo
      vmed(1,kc,jc)=vmed(1,kc,jc)*avgn
      vmed(2,kc,jc)=vmed(2,kc,jc)*avgn
      vmed(3,kc,jc)=vmed(3,kc,jc)*avgn
      pmed(kc,jc)=pmed(kc,jc)*avgn
      pscmed(kc,jc)=pscmed(kc,jc)*avgn
      enddo
      enddo
c
c   now the rms distributions are evaluated
c
      do jc=1,n2m                                                     
      do kc=1,n3m                                                     
      enpsc(kc,jc)=0.
      do l=1,6
      turstr(l,kc,jc)=0.
      enddo
      do ic=1,n1m
      q1p=dq(ic,jc,kc)-vmed(1,kc,jc)
      q2p=dph(ic,jc,kc)-vmed(2,kc,jc)
      q3p=qcap(ic,jc,kc)-vmed(3,kc,jc)
      enejp=(q1p**2+q2p**2+q3p**2)*vl123*rm(jc)*g2rm(jc)
      enej=enej+enejp
      turstr(1,kc,jc)=turstr(1,kc,jc)+q1p**2
      turstr(2,kc,jc)=turstr(2,kc,jc)+q2p**2
      turstr(3,kc,jc)=turstr(3,kc,jc)+q3p**2
      turstr(4,kc,jc)=turstr(4,kc,jc)-q1p*q2p
      turstr(5,kc,jc)=turstr(5,kc,jc)-q3p*q2p
      turstr(6,kc,jc)=turstr(6,kc,jc)-q3p*q1p
      pscp=psc(ic,jc,kc)-pscmed(kc,jc)
      enpsv=(pscp**2)*vl123*rm(jc)*g2rm(jc)+enpsv
      enpsc(kc,jc)=enpsc(kc,jc)+pscp**2
      enddo
      do l=1,3
      turstr(l,kc,jc)=sqrt( turstr(l,kc,jc)*avgn )
      enddo
      turstr(4,kc,jc)= turstr(4,kc,jc)*avgn 
      turstr(5,kc,jc)= turstr(5,kc,jc)*avgn 
      turstr(6,kc,jc)= turstr(6,kc,jc)*avgn 
      enpsc(kc,jc)=enpsc(kc,jc)*avgn
      enddo
      enddo
      enav=enej
      return                                                            
      end                                                               
c
c                                                                       
c  ************************ subrout vorc  **********************  
c                                                                       
c     this subroutine calculates the averaged vorticity               
c     distribution in r-z plane by averaging in theta
c                                                                       
      subroutine vorc
      include 'param.f'
      common/axv/q3ax(m3)
      avgn=1./(float(n1m))
c                                                                       
c  ***********  compute the azimuthal vorticity component               
c               at         i+1/2,j,k 

c                                                                       
      do kc=1,n3
      q3ax(kc)=0.
      do ic=1,n1m
      q3ax(kc)=q3ax(kc)+q3(ic,1,kc)
      enddo
      q3ax(kc)=q3ax(kc)/n1m
      enddo
c
c  inside the field
c
      do ic=1,n1m
      im=imv(ic)
      do jc=2,n2m
      jm=jmv(jc)
      kc=1
      if(inslws.eq.1) then
         dq2x3=q2(ic,jc,1)*dx3*2./rc(jc)/g3rc(1)
      else
         dq2x3=0.
      endif
      dq3x2=(q3(ic,jc,kc)-q3(ic,jm,kc))*dx2/g2rc(jc)
            iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2
      rhs(iadd)=(dq2x3-dq3x2)
      kc=n3
      if(inslwn.eq.1) then
         dq2x3=-q2(ic,jc,n3m)*dx3*2./rc(jc)/g3rc(n3)
      else
         dq2x3=0.
      endif
      dq3x2=(q3(ic,jc,kc)-q3(ic,jm,kc))*dx2/g2rc(jc)
            iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2
      rhs(iadd)=(dq2x3-dq3x2)
      do kc=2,n3m
      km=kmv(kc)
      dq2x3=(q2(ic,jc,kc)-q2(ic,jc,km))*dx3/rc(jc)/g3rc(kc)
      dq3x2=(q3(ic,jc,kc)-q3(ic,jm,kc))*dx2/g2rc(jc)
            iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2
      rhs(iadd)=(dq2x3-dq3x2)
      enddo
      enddo
      enddo
c
c  At  the axis  
c
      jc=1
      jp=jc+1
      do ic=1,n1m
      kc=1
      if(inslws.eq.1) then
      q2akc= (q2(ic,jp,kc) - q2(isym(ic),jp,kc))*0.5/rc(jp)
      q2akm= 0.
      dq2x3=(q2akc-q2akm)*dx3*2./g3rc(kc)
      else
         dq2x3=0.
      endif
      q3axi= (q3(ic,jc,kc) + q3(isym(ic),jc,kc))*0.5
      dq3x2=(q3(ic,jc,kc)-q3axi)*dx2*2./g2rc(jc)
            iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2
      rhs(iadd)=dq2x3-dq3x2
      kc=n3
      if(inslwn.eq.1) then
      q2akc= 0.
      q2akm= (q2(ic,jp,km) - q2(isym(ic),jp,km))*0.5/rc(jp)
      dq2x3=(q2akc-q2akm)*dx3*2./g3rc(kc)
      else
         dq2x3=0.
      endif
      q3axi= (q3(ic,jc,kc) + q3(isym(ic),jc,kc))*0.5
      dq3x2=(q3(ic,jc,kc)-q3axi)*dx2*2./g2rc(jc)
            iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2
      rhs(iadd)=dq2x3-dq3x2
      do kc=2,n3m
      km=kmv(kc)
      q2akc= (q2(ic,jp,kc) - q2(isym(ic),jp,kc))*0.5/rc(jp)
      q2akm= (q2(ic,jp,km) - q2(isym(ic),jp,km))*0.5/rc(jp)
      q3axi= (q3(ic,jc,kc) + q3(isym(ic),jc,kc))*0.5
      dq2x3=(q2akc-q2akm)*dx3/g3rc(kc)
      dq3x2=(q3(ic,jc,kc)-q3axi)*dx2*2./g2rc(jc)
            iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2
      rhs(iadd)=dq2x3-dq3x2
      enddo
      enddo
      jc=n2
      jm=jc-1
      do ic=1,n1m
      do kc=1,n3
c  At the wall  (no-slip)
      dq3x2=(-q3(ic,jm,kc))*dx2*2./g2rc(jc)
            iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2
      rhs(iadd)=-dq3x2
      enddo
      enddo
c
c   azimuthal vorticity at the cell centre 
c
      vomax(1)=-100.
      vomin(1)=+100.
      do jc=1,n2m                                               
      do kc=1,n3m
      vorv(1,kc,jc)=0.
      do ic=1,n1m                                                   
      ip=ipv(ic)                                              
            iada=ic+(jc-1)*n1m+(kc-1)*n1m*n2
            iadc=ic+(jc-1)*n1m+(kc)*n1m*n2
            iadb=ic+(jc)*n1m+(kc-1)*n1m*n2
            iadd=ic+(jc)*n1m+(kc)*n1m*n2
      dq(ic,jc,kc)=(rhs(iada)+rhs(iadb)+
     1              rhs(iadc)+rhs(iadd))*0.25
      vomax(1)=max(vomax(1),dq(ic,jc,kc))
      vomin(1)=min(vomin(1),dq(ic,jc,kc))
      vorv(1,kc,jc)=vorv(1,kc,jc)+dq(ic,jc,kc)
      enddo
      vorv(1,kc,jc)=vorv(1,kc,jc)*avgn
      enddo
      enddo
c                                                                  
c  ***********  compute the radial  vorticity component            
c                           at  i,j+1/2,k
c                                                                  
      do jc=1,n2m
      do ic=1,n1m
      im=imv(ic)
      kc=1
      if(inslws.eq.1) then
         dq1x3=q1(ic,jc,1)*dx3*2.*rm(jc)/g3rc(kc)
      else
         dq1x3=0.
      endif
      dq3x1=(q3(ic,jc,kc)-q3(im,jc,kc))*dx1/rm(jc)
      rhs(iadd)=(dq3x1-dq1x3)
      kc=n3
      if(inslwn.eq.1) then
         dq1x3=-q1(ic,jc,n3m)*dx3*2.*rm(jc)/g3rc(kc)
      else
         dq1x3=0.
      endif
      dq3x1=(q3(ic,jc,kc)-q3(im,jc,kc))*dx1/rm(jc)
            iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2m
      rhs(iadd)=(dq3x1-dq1x3)
      do kc=2,n3m
      km=kmv(kc)
      dq3x1=(q3(ic,jc,kc)-q3(im,jc,kc))*dx1/rm(jc)
      dq1x3=(q1(ic,jc,kc)-q1(ic,jc,km))*dx3/g3rc(kc)
            iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2m
      rhs(iadd)=(dq3x1-dq1x3)
      enddo
      enddo
      enddo

c
c   vorticity at the cell centre 
c
      vomax(2)=-100.
      vomin(2)=+100.
      do jc=1,n2m                                               
      do kc=1,n3m
      vorv(2,kc,jc)=0.
      do ic=1,n1m                                                   
      kp=kpv(kc)
      ip=ipv(ic)                                                  
            iada=ic+(jc-1)*n1m+(kc-1)*n1m*n2m
            iadb=ic+1+(jc-1)*n1m+(kc-1)*n1m*n2m
            iadc=ic+(jc-1)*n1m+(kc)*n1m*n2m
            iadd=ic+1+1+(jc-1)*n1m+(kc)*n1m*n2m
      dph(ic,jc,kc)=(rhs(iada)+rhs(iadb)+
     1               rhs(iadc)+rhs(iadd))*0.25
      vomax(2)=max(vomax(2),dph(ic,jc,kc))
      vomin(2)=min(vomin(2),dph(ic,jc,kc))
      vorv(2,kc,jc)=vorv(2,kc,jc)+dph(ic,jc,kc)
      enddo
      vorv(2,kc,jc)=vorv(2,kc,jc)*avgn
      enddo
      enddo
c                                                                  
c  ***********  compute the vertical vorticity component        
c               at         i,j,k+1/2        
c
      do kc=1,n3m                                                
      do ic=1,n1m                                               
      im=imv(ic)                                                  
      do jc=2,n2m                                               
      jm=jc-1                                                  
      dq1x2=(q1(ic,jc,kc)*rm(jc)-q1(ic,jm,kc)*rm(jm))*dx2/g2rc(jc)
      dq2x1=(q2(ic,jc,kc)-q2(im,jc,kc))*dx1/rc(jc)                 
      vorz=(dq1x2-dq2x1)/rc(jc)
            iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2
      rhs(iadd)=vorz
      enddo
      jc=n2
      jm=n2m
      dq1x2=-1./rc(jc)*q1(ic,jm,kc)*dx2*2./g2rc(jc)*rm(jm)
            iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2
      rhs(iadd)=dq1x2
      enddo
      enddo
      jc=1
      do kc=1,n3m
      vozcm=0.
      do ic=1,n1m
            iadd=ic+(2-1)*n1m+(kc-1)*n1m*n2
      vozcm=vozcm+rhs(iadd)
      enddo
      do ic=1,n1m
            iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2
      rhs(iadd)=vozcm/n1m
      enddo
      enddo

c
c   vorticity at the cell centre 
c
      vomax(3)=-100.
      vomin(3)=+100.
      do jc=1,n2m                                               
      do kc=1,n3m
      vorv(3,kc,jc)=0.
      do ic=1,n1m                                                   
      ip=ipv(ic)                                                  
            iada=ic+(jc-1)*n1m+(kc-1)*n1m*n2
            iadb=ic+(jc)*n1m+(kc-1)*n1m*n2
            iadc=ic+1+(jc-1)*n1m+(kc-1)*n1m*n2
            iadd=ic+1+(jc)*n1m+(kc-1)*n1m*n2
      qcap(ic,jc,kc)=(rhs(iada)+rhs(iadb)+
     1                rhs(iadc)+rhs(iadd))*0.25
      vomax(3)=max(vomax(3),qcap(ic,jc,kc))
      vomin(3)=min(vomin(3),qcap(ic,jc,kc))
      vorv(3,kc,jc)=vorv(3,kc,jc)+qcap(ic,jc,kc)
      enddo
      vorv(3,kc,jc)=vorv(3,kc,jc)*avgn
      enddo
      enddo
      return                                                            
      end                                                               
c
c                                                                       
c  ************************ subrout vortur  **********************  
c                                                                       
c     this subroutine calculates the rms of the vorticity fluctuations 
c                                                                       
      subroutine vortur
      include 'param.f'
      avgn=1./(float(n1m))
      call vorc
      ensj=0.
      do jc=1,n2m                                                     
      jp=jc+1
      volz=rm(jc)*g2rm(jc)/(dx1*dx2*dx3)
      do kc=1,n3m                                                     
      do l=1,6
      vorstr(l,kc,jc)=0.
      enddo
      do ic=1,n1m
      q1p=dq(ic,jc,kc)-vorv(1,kc,jc)
      q2p=dph(ic,jc,kc)-vorv(2,kc,jc)
      q3p=qcap(ic,jc,kc)-vorv(3,kc,jc)
      vorstr(1,kc,jc)=vorstr(1,kc,jc)+q1p**2
      vorstr(2,kc,jc)=vorstr(2,kc,jc)+q2p**2
      vorstr(3,kc,jc)=vorstr(3,kc,jc)+q3p**2
      vorstr(4,kc,jc)=vorstr(4,kc,jc)-q1p*q2p
      vorstr(5,kc,jc)=vorstr(5,kc,jc)-q3p*q2p
      vorstr(6,kc,jc)=vorstr(6,kc,jc)-q3p*q1p
      ensjp=(q1p**2+q2p**2+q3p**2)*volz
      ensj=ensj+ensjp
      enddo
      do l=1,3
      vorstr(l,kc,jc)=sqrt( vorstr(l,kc,jc)*avgn )
      enddo
      vorstr(4,kc,jc)= vorstr(4,kc,jc)*avgn 
      vorstr(5,kc,jc)= vorstr(5,kc,jc)*avgn 
      vorstr(6,kc,jc)= vorstr(6,kc,jc)*avgn 
      enddo
      enddo
      return                                                            
      end                                                               
