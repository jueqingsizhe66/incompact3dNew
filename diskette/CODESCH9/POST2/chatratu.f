c                                                                       
c  ************************* subrout,velc  **********************  
c                                                                       
c     this subroutine calculates the velocity at the centre of the cell
c     and then the fluctuations
c                                                                       
      subroutine velc                    
      include 'param.f'
      avgn=1./(float(n1m*n3m))
      vl123=1./(n1m*n2m*n3m)
      volto=1./(alx1*alx3*2.)
c
c  velocities at the center of the cell
c
      do jc=1,n2m                                                     
      jm=jmv(jc)                                                        
      jp=jc+1
      do l=1,3
      vmed(l,jc)=0.
      vmepo(l,jc)=0.
      enddo
      do kc=1,n3m 
      do ic=1,n1m                                                     
      dphc=
     1  (    q2(ic,jc,kc)+q2(ic,jp,kc))*0.5
      q2akc=q2(ic,jc,kc)
      dqc=
     1  (    q1(ic,jc,kc)+q1(ipv(ic),jc,kc) )*0.5
      qcapc=(q3(ic,jc,kc)+q3(ic,jc,kpv(kc)) )*0.5 
      vmed(2,jc)=vmed(2,jc)+dphc
      vmed(1,jc)=vmed(1,jc)+dqc
      vmed(3,jc)=vmed(3,jc)+qcapc
      vmepo(2,jc)=vmepo(2,jc)+q2akc
      vmepo(1,jc)=vmepo(1,jc)+q1(ic,jc,kc)
      vmepo(3,jc)=vmepo(3,jc)+q3(ic,jc,kc)
      enddo
      enddo
      do l=1,3
      vmed(l,jc)=vmed(l,jc)*avgn
      vmepo(l,jc)=vmepo(l,jc)*avgn
      enddo
      enddo
      pcost=0.
      do jc=1,n2m
      do kc=1,n3m
      do ic=1,n1m
      pcost=pcost+pr(ic,jc,kc)
      enddo
      enddo
      enddo
      pcost=pcost/float(n1m*n2m*n3m)
      vit(1)=0.
      vit(2)=0.
      vit(3)=0.
      vit(4)=0.
      vit(5)=0.
c
c   mean pressure and passive scalar profiles
c
      do jc=1,n2m                                                     
      jm=jmv(jc)                                                        
      jp=jc+1
      pmed(jc)=0.
      do kc=1,n3m 
      do ic=1,n1m                                                     
      prver=pr(ic,jc,kc)
      pmed(jc)=pmed(jc)+prver
      vit(1)=q1(ic,jc,kc)*volz(jc)+vit(1)
      vit(2)=(q2(ic,jp,kc)+q2(ic,jc,kc))*volz(jc)+vit(2)
      vit(3)=q3(ic,jc,kc)*volz(jc)+vit(3)
      vit(4)=pr(ic,jc,kc)*volz(jc)+vit(4)
      vit(5)=pr(ic,jc,kc)*volz(jc)+vit(5)
      enddo
      enddo
      pmed(jc)=pmed(jc)*avgn
      enddo
      vit(1)=vit(1)*volto
      vit(2)=vit(2)*volto
      vit(3)=vit(3)*volto
      vit(4)=vit(4)*volto
      vit(5)=vit(5)*volto
      return
      end
c                                                                       
c  ************************* subrout,veltur  **********************  
c                                                                       
c     this subroutine calculates the turbulent stress                  
c                                                                       
      subroutine veltur(ntime)
      include 'param.f'
      common/phto/phtot(m1,m2,m3)
      common/prcost/pcost
      character*60 filns1,filts1,filns2,filts2,filns3,filts3
      character*2 pnt
      avgn=1./(float(n1m*n3m))
      vl123=1./(n1m*n2m*n3m)
      volto=1./(alx1*alx3*2.)
  133 format(1x,'veltur',2x,5e10.3)
      hare=(rav/ren)**2/pram
      enej=0.
      do jc=1,n2m                                                     
      jp=jc+1
      do l=1,4
      ske(l,jc)=0.
      fla(l,jc)=0.
      pvc(l,jc)=0.
      enddo
      do l=1,5
      vtv(l,jc)=0.
      enddo
      do l=1,6
      turstr(l,jc)=0.
      enddo
      do kc=1,n3m
      do ic=1,n1m                                                     
c
c    calculation velocities at  the cell centre
c
      dphc=
     1  (    q2(ic,jc,kc)+q2(ic,jp,kc))*0.5
      dqc=
     1  (    q1(ic,jc,kc)+q1(ipv(ic),jc,kc) )*0.5
      qcapc=(q3(ic,jc,kc)+q3(ic,jc,kpv(kc)) )*0.5 
c
c   calculation fluctuating velocities
c
      q1p=dqc-vmed(1,jc)
      q2p=dphc-vmed(2,jc)
      q3p=qcapc-vmed(3,jc)
      prver=pr(ic,jc,kc)
      pcp=prver-pmed(jc)
      enejp=(q1p**2+q2p**2+q3p**2)*vl123*caj(jc)
      enej=enej+enejp
      phtot(ic,jc,kc)=pcp+(q1p**2+q2p**2+q3p**2)*0.5
c
c   calculation normal and Reynolds stress
c
      turstr(1,jc)=turstr(1,jc)+q1p**2
      turstr(2,jc)=turstr(2,jc)+q2p**2
      turstr(3,jc)=turstr(3,jc)+q3p**2
      turstr(4,jc)=turstr(4,jc)-q1p*q2p
      turstr(5,jc)=turstr(5,jc)-q3p*q2p
      turstr(6,jc)=turstr(6,jc)-q3p*q1p
c
c   evaluation of pressure velocity correlations
c
      pvc(1,jc)=pvc(1,jc)+q1p*pcp
      pvc(2,jc)=pvc(2,jc)+q2p*pcp
      pvc(3,jc)=pvc(3,jc)+q3p*pcp
      pvc(4,jc)=pvc(4,jc)+pcp*pcp
c
c   evaluation of velocity skewness
c
      ske(1,jc)=ske(1,jc)+q1p**3
      ske(2,jc)=ske(2,jc)+q2p**3
      ske(3,jc)=ske(3,jc)+q3p**3
      ske(4,jc)=ske(4,jc)+pcp**3
c
c   evaluation of velocity flatness
c
      fla(1,jc)=fla(1,jc)+q1p**4
      fla(2,jc)=fla(2,jc)+q2p**4
      fla(3,jc)=fla(3,jc)+q3p**4
      fla(4,jc)=fla(4,jc)+pcp**4
c
c   evaluation of mean normal velocity
c   vmean2 =<V'jV'jV'2>/<V'jV'j>
c
      vtv(1,jc)=-q2p*(q1p**2+q2p**2+q3p**2)+vtv(1,jc)
      vtv(2,jc)=-q1p*q2p**2+vtv(2,jc)
      vtv(3,jc)=-q3p*q2p**2+vtv(3,jc)
      vtv(4,jc)=-q2p*q2p**2+vtv(4,jc)
      vtv(5,jc)=-pcp*q2p**2+vtv(5,jc)
      enddo
      enddo
      do l=1,3
      turstr(l,jc)=sqrt( turstr(l,jc)*avgn )
      enddo
      vtv(1,jc)= vtv(1,jc)*avgn 
      vtv(2,jc)= vtv(2,jc)*avgn 
      vtv(3,jc)= vtv(3,jc)*avgn 
      vtv(4,jc)= vtv(4,jc)*avgn 
      vtv(5,jc)= vtv(5,jc)*avgn 
      pvc(1,jc)= pvc(1,jc)*avgn 
      pvc(2,jc)= pvc(2,jc)*avgn 
      pvc(3,jc)= pvc(3,jc)*avgn 
      pvc(4,jc)= pvc(4,jc)*avgn 
      turstr(4,jc)= turstr(4,jc)*avgn 
      turstr(5,jc)= turstr(5,jc)*avgn 
      turstr(6,jc)= turstr(6,jc)*avgn 
      ske(1,jc)= ske(1,jc)*avgn 
      ske(2,jc)= ske(2,jc)*avgn 
      ske(3,jc)= ske(3,jc)*avgn 
      ske(4,jc)= ske(4,jc)*avgn 
      fla(1,jc)= fla(1,jc)*avgn 
      fla(2,jc)= fla(2,jc)*avgn 
      fla(3,jc)= fla(3,jc)*avgn 
      fla(4,jc)= fla(4,jc)*avgn 
      enddo
      enav=enej
  612 format(1x,2e15.5)
      return                                                            
      end                                                               
c
c                                                                       
c  ************************ subrout vorc  **********************  
c                                                                       
c     this subroutine calculates the azimuthal vorticity               
c                                                                       
      subroutine vorc
      include 'param.f'
      avgn=1./(float(n1m*n3m))
c                                                                       
c  ***********  compute the spanwise vorticity component               
c               at         i+1/2,j,k 

c                                                                       
c  inside the field
      do jc=2,n2m
      jm=jmv(jc)
      vompo(1,jc)=0.
            do ic=1,n1m
      im=imv(ic)
                  do kc=1,n3m
      km=kmv(kc)
      dq3x2=(q3(ic,jc,kc)-q3(ic,jm,kc))/(y2s(jc)-y2s(jm))
      dq2x3=(q2(ic,jc,kc)-q2(ic,jc,km))*dx3
      rhs(ic,jc,kc)=-(dq2x3-dq3x2)
      vompo(1,jc)=vompo(1,jc)+rhs(ic,jc,kc)
                  enddo
            enddo
      vompo(1,jc)=vompo(1,jc)*avgn
      enddo
      jc=1
      vompo(1,jc)=0.
      jp=jc+1
            do ic=1,n1m
                  do kc=1,n3m
      km=kmv(kc)
      dq3x2=q3(ic,jc,kc)/(y2s(1)-y(1))
      dq2x3=(q2(ic,jc,kc)-q2(ic,jc,km))*dx3
      rhs(ic,jc,kc)=-(dq2x3-dq3x2)
      vompo(1,jc)=vompo(1,jc)+rhs(ic,jc,kc)
                  enddo
            enddo
      vompo(1,jc)=vompo(1,jc)*avgn
      jc=n2
      vompo(1,jc)=0.
      jm=jc-1
            do ic=1,n1m
                  do kc=1,n3m
c  At the wall  (no-slip)
      dq3x2=(-q3(ic,jm,kc))/(y(n2)-y2s(n2m))
      dq2x3=(q2(ic,jc,kc)-q2(ic,jc,km))*dx3
      rhs(ic,jc,kc)=-(dq2x3-dq3x2)
      vompo(1,jc)=vompo(1,jc)+rhs(ic,jc,kc)
                  enddo
            enddo
      vompo(1,jc)=vompo(1,jc)*avgn
c
c   spanwise vorticity at the cell centre 
c
      do jc=2,n2m-1                                             
      vorv(1,jc)=0.
            do kc=1,n3m                                                   
      kp=kpv(kc)                                              
                  do ic=1,n1m
      dqc=(rhs(ic,jc,kc)+rhs(ic,jc+1,kc)+
     1     rhs(ic,jc,kp)+rhs(ic,jc+1,kp))*0.25
      vorv(1,jc)=vorv(1,jc)+dqc
                  enddo
            enddo
      vorv(1,jc)=vorv(1,jc)*avgn
      enddo
c
c    wall spanwise. vort
c
      jc=1                                               
      vorv(1,jc)=0.
            do kc=1,n3m                                                   
      kp=kpv(kc)                                              
                  do ic=1,n1m
      dqc=(rhs(ic,jc,kc)+rhs(ic,jc+1,kc)+
     1     rhs(ic,jc,kp)+rhs(ic,jc+1,kp))*0.25
      vorv(1,jc)=vorv(1,jc)+dqc
      vocl(1,ic,kc)=dqc
                  enddo
            enddo
      vorv(1,jc)=vorv(1,jc)*avgn
      do jc=1,n2m                                               
      vorv(1,jc)=0.
            do kc=1,n3m                                                   
      kp=kpv(kc)                                              
                  do ic=1,n1m
      dqc=(rhs(ic,jc,kc)+rhs(ic,jc+1,kc)+
     1     rhs(ic,jc,kp)+rhs(ic,jc+1,kp))*0.25
      vorv(1,jc)=vorv(1,jc)+dqc
      vowa(1,ic,kc)=dqc
                  enddo
            enddo
      vorv(1,jc)=vorv(1,jc)*avgn
      enddo
c                                                                  
c  ***********  compute the normal  vorticity component            
c                           at  i,j+1/2,k
c                                                                  
      do jc=1,n2m
      vompo(2,jc)=0.
            do ic=1,n1m
      im=imv(ic)
                  do kc=1,n3m
      km=kmv(kc)
      dq3x1=(q3(ic,jc,kc)-q3(im,jc,kc))*dx1
      dq1x3=(q1(ic,jc,kc)-q1(ic,jc,km))*dx3
      rhs(ic,jc,kc)=-(dq3x1-dq1x3)
      vompo(2,jc)=vompo(2,jc)+rhs(ic,jc,kc)
                  enddo
            enddo
      vompo(2,jc)=vompo(2,jc)*avgn
      enddo

c
c   vorticity at the cell centre inner field
c
      do jc=2,n2m-1
      vorv(2,jc)=0.
                  do kc=1,n3m                                                   
      kp=kpv(kc)
                  do ic=1,n1m
      ip=ipv(ic)                                                  
      dphc=(rhs(ic,jc,kc)+rhs(ip,jc,kc)+
     1      rhs(ic,jc,kp)+rhs(ip,jc,kp))*0.25
      vorv(2,jc)=vorv(2,jc)+dphc
                  enddo
                  enddo
      vorv(2,jc)=vorv(2,jc)*avgn
      enddo
c
c   vorticity at the cell centre upper wall        
c
      jc=n2m
      vorv(2,jc)=0.
            do kc=1,n3m                                                   
      kp=kpv(kc)
                  do ic=1,n1m
      ip=ipv(ic)                                                  
      dphc=(rhs(ic,jc,kc)+rhs(ip,jc,kc)+
     1      rhs(ic,jc,kp)+rhs(ip,jc,kp))*0.25
      vowa(2,ic,kc)=dphc
      vorv(2,jc)=vorv(2,jc)+dphc
                  enddo
            enddo
      vorv(2,jc)=vorv(2,jc)*avgn
c
c   vorticity at the cell centre lower wall
c
      jc=1
      vorv(2,jc)=0.
                  do kc=1,n3m                                                   
      kp=kpv(kc)
            do ic=1,n1m
      ip=ipv(ic)                                                  
      dphc=(rhs(ic,jc,kc)+rhs(ip,jc,kc)+
     1      rhs(ic,jc,kp)+rhs(ip,jc,kp))*0.25
      vocl(2,ic,kc)=dphc
      vorv(2,jc)=vorv(2,jc)+dphc
            enddo
                  enddo
      vorv(2,jc)=vorv(2,jc)*avgn
c                                                                  
c  ***********  compute the streamwise vorticity component        
c               at         i,j,k+1/2        
c
      do jc=2,n2m                                               
      vompo(3,jc)=0.
      jm=jc-1                                                  
            do kc=1,n3m                                                
                  do ic=1,n1m                                               
      im=imv(ic)                                                  
      dq1x2=(q1(ic,jc,kc)-q1(ic,jm,kc))/(y2s(jc)-y2s(jm))
      dq2x1=(q2(ic,jc,kc)-q2(im,jc,kc))*dx1
      vorz=-(dq1x2-dq2x1)
      rhs(ic,jc,kc)=vorz
      vompo(3,jc)=vompo(3,jc)+rhs(ic,jc,kc)
                  enddo
            enddo
      vompo(3,jc)=vompo(3,jc)*avgn
      enddo
      jc=n2
      vompo(3,jc)=0.
      jm=n2m
            do kc=1,n3m                                                
                  do ic=1,n1m                                               
      im=imv(ic)                                                  
      dq1x2=-q1(ic,jm,kc)/(y(jc)-y2s(jm))
      dq2x1=(q2(ic,jc,kc)-q2(im,jc,kc))*dx1
      vorz=-(dq1x2-dq2x1)
      rhs(ic,jc,kc)=vorz
      vompo(3,jc)=vompo(3,jc)+rhs(ic,jc,kc)
                  enddo
            enddo
      vompo(3,jc)=vompo(3,jc)*avgn
      jc=1
      vompo(3,jc)=0.
            do kc=1,n3m                                                
                  do ic=1,n1m                                               
      im=imv(ic)                                                  
      dq1x2=q1(ic,jc,kc)/(y2s(jc)-y(jc))
      dq2x1=(q2(ic,jc,kc)-q2(im,jc,kc))*dx1
      vorz=-(dq1x2-dq2x1)
      rhs(ic,jc,kc)=vorz
      vompo(3,jc)=vompo(3,jc)+rhs(ic,jc,kc)
                  enddo
            enddo
      vompo(3,jc)=vompo(3,jc)*avgn
c
c   vorticity at the cell centre 
c
      do jc=2,n2m-1
      vorv(3,jc)=0.
            do kc=1,n3m                                                   
                  do ic=1,n1m
      ip=ipv(ic)                                                  
      qcapc=(rhs(ic,jc,kc)+rhs(ic,jc+1,kc)+
     1       rhs(ip,jc,kc)+rhs(ip,jc+1,kc))*0.25
      vorv(3,jc)=vorv(3,jc)+qcapc
                  enddo
            enddo
      vorv(3,jc)=vorv(3,jc)*avgn
      enddo
c
c   vorticity at the cell centre center line
c
      jc=1                                               
      vorv(3,jc)=0.
            do kc=1,n3m                                                   
                  do ic=1,n1m
      ip=ipv(ic)                                                  
      qcapc=(rhs(ic,jc,kc)+rhs(ic,jc+1,kc)+
     1       rhs(ip,jc,kc)+rhs(ip,jc+1,kc))*0.25
      vorv(3,jc)=vorv(3,jc)+qcapc
      vocl(3,ic,kc)=qcapc
                  enddo
            enddo
      vorv(3,jc)=vorv(3,jc)*avgn
c
c   vorticity at the cell centre wall    
c
      jc=n2m                                               
      vorv(3,jc)=0.
            do kc=1,n3m                                                   
                  do ic=1,n1m
      ip=ipv(ic)                                                  
      qcapc=(rhs(ic,jc,kc)+rhs(ic,jc+1,kc)+
     1       rhs(ip,jc,kc)+rhs(ip,jc+1,kc))*0.25
      vowa(3,ic,kc)=qcapc
      vorv(3,jc)=vorv(3,jc)+qcapc
                  enddo
            enddo
      vorv(3,jc)=vorv(3,jc)*avgn
      return                                                            
      end                                                               
c
c                                                                       
c  ************************ subrout quasfa  **********************  
c                                                                       
c     this subroutine calculates the averages of certain quantities 
c                                                                       
      subroutine quasfa(ntime)
      include 'param.f'
      avgn=1./(float(n1m*n3m))
      do jc=n2m,1,-1                                                 
      jp=jc+1
      jm=jc-1
            do l=1,12
      quaav(l,jc)=0.
            enddo
      do kc=1,n3m
      kp=kpv(kc)
      km=kmv(kc)
      do ic=1,n1m                                                     
      ip=ipv(ic)
      im=imv(ic)
      if(jc.ge.2.and.jc.le.n2m-1) then
c
c   spanwise vorticity
c
      dq2x3c=(q2(ic,jc,kc)-q2(ic,jc,km))*dx3
      dq3x2c=(q3(ic,jc,kc)-q3(ic,jm,kc))/(y2s(jc)-y2s(jm))
      dq2x3a=(q2(ic,jp,kc)-q2(ic,jp,km))*dx3
      dq3x2a=(q3(ic,jp,kc)-q3(ic,jc,kc))/(y2s(jp)-y2s(jc))
      dq2x3d=(q2(ic,jc,kp)-q2(ic,jc,kc))*dx3
      dq3x2d=(q3(ic,jc,kp)-q3(ic,jm,kp))/(y2s(jc)-y2s(jm))
      dq2x3b=(q2(ic,jp,kp)-q2(ic,jp,kc))*dx3
      dq3x2b=(q3(ic,jp,kp)-q3(ic,jc,kp))/(y2s(jp)-y2s(jc))
      dv3x2c=(vmepo(3,jc)-vmepo(3,jm))/(y2s(jc)-y2s(jm))
      dv3x2a=(vmepo(3,jp)-vmepo(3,jc))/(y2s(jp)-y2s(jc))
      dv3mx2=(dv3x2a+dv3x2c)*0.5
      dq2x3=(dq2x3c+dq2x3a+dq2x3b+dq2x3d)*0.25
      dq3x2=(dq3x2c+dq3x2a+dq3x2b+dq3x2d)*0.25
      dqc=-(dq2x3-(dq3x2-dv3mx2))
      q1p=dqc
      s23=(dq2x3+(dq3x2-dv3mx2))*0.5
c
c   normal vorticity
c
      dq3x1c=(q3(ic,jc,kc)-q3(im,jc,kc))*dx1
      dq1x3c=(q1(ic,jc,kc)-q1(ic,jc,km))*dx3
      dq3x1a=(q3(ip,jc,kc)-q3(ic,jc,kc))*dx1
      dq1x3a=(q1(ip,jc,kc)-q1(ip,jc,km))*dx3
      dq3x1b=(q3(ip,jc,kp)-q3(ic,jc,kp))*dx1
      dq1x3b=(q1(ip,jc,kp)-q1(ip,jc,kc))*dx3
      dq3x1d=(q3(ic,jc,kp)-q3(im,jc,kp))*dx1
      dq1x3d=(q1(ic,jc,kp)-q1(ic,jc,kc))*dx3
      dq1x3=(dq1x3c+dq1x3a+dq1x3b+dq1x3d)*0.25
      dq3x1=(dq3x1c+dq3x1a+dq3x1b+dq3x1d)*0.25
      dphc=-(dq3x1-dq1x3)
      q2p=dphc
      s13=(dq1x3+dq3x1)*0.5
c
c   streamwise vorticity
c
      dq1x2c=(q1(ic,jc,kc)-q1(ic,jm,kc))/(y2s(jc)-y2s(jm))
      dq2x1c=(q2(ic,jc,kc)-q2(im,jc,kc))*dx1
      dq1x2a=(q1(ic,jp,kc)-q1(ic,jc,kc))/(y2s(jp)-y2s(jc))
      dq2x1a=(q2(ic,jp,kc)-q2(im,jp,kc))*dx1
      dq1x2d=(q1(ip,jc,kc)-q1(ip,jm,kc))/(y2s(jc)-y2s(jm))
      dq2x1d=(q2(ip,jc,kc)-q2(ic,jc,kc))*dx1
      dq1x2b=(q1(ip,jp,kc)-q1(ip,jc,kc))/(y2s(jp)-y2s(jc))
      dq2x1b=(q2(ip,jp,kc)-q2(ic,jp,kc))*dx1
      dv1x2c=(vmepo(1,jc)-vmepo(1,jm))/(y2s(jc)-y2s(jm))
      dv1x2a=(vmepo(1,jp)-vmepo(1,jc))/(y2s(jp)-y2s(jc))
      dv1mx2=(dv1x2a+dv1x2c)*0.5
      dq1x2=(dq1x2c+dq1x2a+dq1x2b+dq1x2d)*0.25
      dq2x1=(dq2x1c+dq2x1a+dq2x1b+dq2x1d)*0.25
      qcapc=-(dq1x2-dv1mx2-dq2x1)
      q3p=qcapc
      s12=(dq1x2-dv1mx2+dq2x1)*0.5
      dq1x1=(q1(ip,jc,kc)-q1(ic,jc,kc))*dx1
      dq3x3=(q3(ic,jc,kp)-q3(ic,jc,kc))*dx3
      dv2mx2=(vmepo(2,jp)-vmepo(2,jc))/(yp2(jp)-yp2(jc))
      dq2x2=(q2(ic,jp,kc)-q2(ic,jc,kc))/(yp2(jp)-yp2(jc))
      s11=dq1x1
      s22=dq2x2-dv2mx2
      s33=dq3x3
      dislo=2.*(s11**2+s22**2+s33**2+
     1      2.*(s12**2+s13**2+s23**2))
                                  else
      if(jc.eq.1) then
c
c   fluctuating vorticities at the cell centre  wall
c
c
c   spanwise vorticity
c
      dq2x3c=(q2(ic,jc,kc)-q2(ic,jc,km))*dx3
      dq3x2c=(q3(ic,jc,kc))/(y2s(jc)-yp2(jc))
      dq2x3a=(q2(ic,jp,kc)-q2(ic,jp,km))*dx3
      dq3x2a=(q3(ic,jp,kc)-q3(ic,jc,kc))/(y2s(jp)-y2s(jc))
      dq2x3d=(q2(ic,jc,kp)-q2(ic,jc,kc))*dx3
      dq3x2d=(q3(ic,jc,kp))/(y2s(jc)-yp2(jc))
      dq2x3b=(q2(ic,jp,kp)-q2(ic,jp,kc))*dx3
      dq3x2b=(q3(ic,jp,kp)-q3(ic,jc,kp))/(y2s(jp)-y2s(jc))
      dq2x3=(dq2x3c+dq2x3a+dq2x3b+dq2x3d)*0.25
      dv3x2c=(vmepo(3,jc))/(y2s(jc)-yp2(jc))
      dv3x2a=(vmepo(3,jp)-vmepo(3,jc))/(y2s(jp)-y2s(jc))
      dv3mx2=(dv3x2a+dv3x2c)*0.5
      dq3x2=(dq3x2c+dq3x2a+dq3x2b+dq3x2d)*0.25
      dqc=-(dq2x3-(dq3x2-dv3mx2))
      q1p=dqc
      s23=(dq2x3+(dq3x2-dv3mx2))*0.5
c
c   normal vorticity
c
      dq3x1c=(q3(ic,jc,kc)-q3(im,jc,kc))*dx1
      dq1x3c=(q1(ic,jc,kc)-q1(ic,jc,km))*dx3
      dq3x1a=(q3(ip,jc,kc)-q3(ic,jc,kc))*dx1
      dq1x3a=(q1(ip,jc,kc)-q1(ip,jc,km))*dx3
      dq3x1b=(q3(ip,jc,kp)-q3(ic,jc,kp))*dx1
      dq1x3b=(q1(ip,jc,kp)-q1(ip,jc,kc))*dx3
      dq3x1d=(q3(ic,jc,kp)-q3(im,jc,kp))*dx1
      dq1x3d=(q1(ic,jc,kp)-q1(ic,jc,kc))*dx3
      dq1x3=(dq1x3c+dq1x3a+dq1x3b+dq1x3d)*0.25
      dq3x1=(dq3x1c+dq3x1a+dq3x1b+dq3x1d)*0.25
      dphc=-(dq3x1-dq1x3)
      q2p=dphc
      s13=(dq1x3+dq3x1)*0.5
c
c   streamwise vorticity
c
      dq1x2c=(q1(ic,jc,kc))/(y2s(jc)-yp2(jc))
      dq2x1c=(q2(ic,jc,kc)-q2(im,jc,kc))*dx1
      dq1x2a=(q1(ic,jp,kc)-q1(ic,jc,kc))/(y2s(jp)-y2s(jc))
      dq2x1a=(q2(ic,jp,kc)-q2(im,jp,kc))*dx1
      dq1x2d=(q1(ip,jc,kc))/(y2s(jc)-yp2(jc))
      dq2x1d=(q2(ip,jc,kc)-q2(ic,jc,kc))*dx1
      dq1x2b=(q1(ip,jp,kc)-q1(ip,jc,kc))/(y2s(jp)-y2s(jc))
      dq2x1b=(q2(ip,jp,kc)-q2(ic,jp,kc))*dx1
      dv1x2c=(vmepo(1,jc))/(y2s(jc)-yp2(jc))
      dv1x2a=(vmepo(1,jp)-vmepo(1,jc))/(y2s(jp)-y2s(jc))
      dv1mx2=(dv1x2a+dv1x2c)*0.5
      dq1x2=(dq1x2c+dq1x2a+dq1x2b+dq1x2d)*0.25
      dq2x1=(dq2x1c+dq2x1a+dq2x1b+dq2x1d)*0.25
      qcapc=-(dq1x2-dv1mx2-dq2x1)
      qcapcp=-(dq1x2-dv1mx2*(1-itot)-dq2x1)
      q3p=qcapc
      s12=(dq1x2-dv1mx2+dq2x1)*0.5
      dq1x1=(q1(ip,jc,kc)-q1(ic,jc,kc))*dx1
      dq3x3=(q3(ic,jc,kp)-q3(ic,jc,kc))*dx3
      dv2mx2=(vmepo(2,jp)-vmepo(2,jc))/(yp2(jp)-yp2(jc))
      dq2x2=(q2(ic,jp,kc)-q2(ic,jc,kc))/(yp2(jp)-yp2(jc))-dv2mx2
      s11=dq1x1
      s22=dq2x2-dv2mx2
      s33=dq3x3
      dislo=2.*(s11**2+s22**2+s33**2+
     1      2.*(s12**2+s13**2+s23**2))
                  endif
      if(jc.eq.n2m) then
c
c   fluctuating vorticities at the cell centre  wall
c
c
c   spanwise vorticity
c
      dq2x3c=(q2(ic,jc,kc)-q2(ic,jc,km))*dx3
      dq3x2c=(q3(ic,jc,kc)-q3(ic,jm,kc))/(y2s(jc)-y2s(jm))
      dq2x3a=(q2(ic,jp,kc)-q2(ic,jp,km))*dx3
      dq3x2a=(-q3(ic,jc,kc))/(yp2(jp)-y2s(jc))
      dq2x3d=(q2(ic,jc,kp)-q2(ic,jc,kc))*dx3
      dq3x2d=(q3(ic,jc,kp)-q3(ic,jm,kp))/(y2s(jc)-y2s(jm))
      dq2x3b=(q2(ic,jp,kp)-q2(ic,jp,kc))*dx3
      dq3x2b=(-q3(ic,jc,kp))/(yp2(jp)-y2s(jc))
      dv3x2c=(vmepo(3,jc)-vmepo(3,jm))/(y2s(jc)-y2s(jm))
      dv3x2a=(-vmepo(3,jc))/(yp2(jp)-y2s(jc))
      dv3mx2=(dv3x2a+dv3x2c)*0.5
      dq2x3=(dq2x3c+dq2x3a+dq2x3b+dq2x3d)*0.25
      dq3x2=(dq3x2c+dq3x2a+dq3x2b+dq3x2d)*0.25
      dqc=-(dq2x3-(dq3x2-dv3mx2))
      q1p=dqc
      s23=(dq2x3+(dq3x2-dv3mx2))*0.5
c
c   normal vorticity
c
      dq3x1c=(q3(ic,jc,kc)-q3(im,jc,kc))*dx1
      dq1x3c=(q1(ic,jc,kc)-q1(ic,jc,km))*dx3
      dq3x1a=(q3(ip,jc,kc)-q3(ic,jc,kc))*dx1
      dq1x3a=(q1(ip,jc,kc)-q1(ip,jc,km))*dx3
      dq3x1b=(q3(ip,jc,kp)-q3(ic,jc,kp))*dx1
      dq1x3b=(q1(ip,jc,kp)-q1(ip,jc,kc))*dx3
      dq3x1d=(q3(ic,jc,kp)-q3(im,jc,kp))*dx1
      dq1x3d=(q1(ic,jc,kp)-q1(ic,jc,kc))*dx3
      dq1x3=(dq1x3c+dq1x3a+dq1x3b+dq1x3d)*0.25
      dq3x1=(dq3x1c+dq3x1a+dq3x1b+dq3x1d)*0.25
      dphc=-(dq3x1-dq1x3)
      q2p=dphc
      s13=(dq1x3+dq3x1)*0.5
c
c   streamwise vorticity
c
      dq1x2c=(q1(ic,jc,kc)-q1(ic,jm,kc))/(y2s(jc)-y2s(jm))
      dq2x1c=(q2(ic,jc,kc)-q2(im,jc,kc))*dx1
      dq1x2a=(-q1(ic,jc,kc))/(yp2(jp)-y2s(jc))
      dq2x1a=(q2(ic,jp,kc)-q2(im,jp,kc))*dx1
      dq1x2d=(q1(ip,jc,kc)-q1(ip,jm,kc))/(y2s(jc)-y2s(jm))
      dq2x1d=(q2(ip,jc,kc)-q2(ic,jc,kc))*dx1
      dq1x2b=(-q1(ip,jc,kc))/(yp2(jp)-y2s(jc))
      dq2x1b=(q2(ip,jp,kc)-q2(ic,jp,kc))*dx1
      dv1x2c=(vmepo(1,jc)-vmepo(1,jm))/(y2s(jc)-y2s(jm))
      dv1x2a=(-vmepo(1,jc))/(yp2(jp)-y2s(jc))
      dv1mx2=(dv1x2a+dv1x2c)*0.5
      dq1x2=(dq1x2c+dq1x2a+dq1x2b+dq1x2d)*0.25
      dq2x1=(dq2x1c+dq2x1a+dq2x1b+dq2x1d)*0.25
      qcapc=-(dq1x2-dv1mx2-dq2x1)
      qcapcp=-(dq1x2-dv1mx2*(1-itot)-dq2x1)
      q3p=qcapc
      s12=(dq1x2-dv1mx2+dq2x1)*0.5
      dq1x1=(q1(ip,jc,kc)-q1(ic,jc,kc))*dx1
      dq3x3=(q3(ic,jc,kp)-q3(ic,jc,kc))*dx3
      dv2mx2=(vmepo(2,jp)-vmepo(2,jc))/(yp2(jp)-yp2(jc))
      dq2x2=(q2(ic,jp,kc)-q2(ic,jc,kc))/(yp2(jp)-yp2(jc))
      s11=dq1x1
      s22=dq2x2-dv2mx2
      s33=dq3x3
      dislo=2.*(s11**2+s22**2+s33**2+
     1      2.*(s12**2+s13**2+s23**2))
                    endif
                                  endif
      voplo=sqrt(q1p**2+q2p**2+q3p**2)
      enslo=voplo**2
      plap=enslo/4.-dislo
c
c   average of vortisity enstr, dissipation rate , strains
c   and Laplacian of pressure
c
      quaav(1,jc)=q1p+quaav(1,jc)
      quaav(2,jc)=q2p+quaav(2,jc)
      quaav(3,jc)=q3p+quaav(3,jc)
      quaav(4,jc)=s11+quaav(4,jc)
      quaav(5,jc)=s22+quaav(5,jc)
      quaav(6,jc)=s33+quaav(6,jc)
      quaav(7,jc)=s12+quaav(7,jc)
      quaav(8,jc)=s23+quaav(8,jc)
      quaav(9,jc)=s13+quaav(9,jc)
      quaav(10,jc)=dislo+quaav(10,jc)
      quaav(11,jc)=enslo+quaav(11,jc)
      quaav(12,jc)=plap+quaav(12,jc)
             enddo
             enddo
          do l=1,12
      quaav(l,jc)=quaav(l,jc)*avgn
          enddo
      enddo
      write(74,*)'  field at ntime=',ntime,' averages'
      do jc=1,n2m
      write(74,174)jc,(quaav(l,jc),l=1,12)
  174 format(i4,2x,12e10.3)
      enddo
      return                                                            
      end                                                               
c
c                                                                       
c  ************************ subrout vortur  **********************  
c                                                                       
c     this subroutine calculates the fluctuating 
c     vorticity components and strain rate
c                                                                       
      subroutine vortur(ntime)
      include 'param.f'
      common/phto/phtot(m1,m2,m3)
      common/ipdf/iquapd
      avgn=1./(float(n1m*n3m))
      asma1=1.
      orma1=1.
   83 format(i3.3)
  144 format(9e14.6)
      hare=(rav/ren)**2/pram
  117 format(3x,i5,3x,4e15.6)
      do n=1,mpq
      enssma(n)=0.
      dissma(n)=0.
      esprma(n)=0.
      espma(n)=-100.
      espmi(n)=+100.
      enddo
      npq=1
      do jc=n2m,1,-1                                                 
      jp=jc+1
      jm=jc-1
      dissj(jc)=0.
      ensy(jc)=0.
      enspr(jc)=0.
            do l=1,3
      vdtomt(l,jc)=0.
      prgrd(l,jc)=0.
            enddo
            do l=1,6
      vorstr(l,jc)=0.
      vcromt(l,jc)=0.
            enddo
            do l=1,12
      quafo(l,jc)=0.
      quacu(l,jc)=0.
      quasq(l,jc)=0.
            enddo
      do kc=1,n3m
      kp=kpv(kc)
      km=kmv(kc)
      do ic=1,n1m                                                     
      ip=ipv(ic)
      im=imv(ic)
      if(jc.ge.2.and.jc.le.n2m-1) then
c
c   spanwise vorticity
c
      dq2x3c=(q2(ic,jc,kc)-q2(ic,jc,km))*dx3
      dq3x2c=(q3(ic,jc,kc)-q3(ic,jm,kc))/(y2s(jc)-y2s(jm))
      dq2x3a=(q2(ic,jp,kc)-q2(ic,jp,km))*dx3
      dq3x2a=(q3(ic,jp,kc)-q3(ic,jc,kc))/(y2s(jp)-y2s(jc))
      dq2x3d=(q2(ic,jc,kp)-q2(ic,jc,kc))*dx3
      dq3x2d=(q3(ic,jc,kp)-q3(ic,jm,kp))/(y2s(jc)-y2s(jm))
      dq2x3b=(q2(ic,jp,kp)-q2(ic,jp,kc))*dx3
      dq3x2b=(q3(ic,jp,kp)-q3(ic,jc,kp))/(y2s(jp)-y2s(jc))
      dv3x2c=(vmepo(3,jc)-vmepo(3,jm))/(y2s(jc)-y2s(jm))
      dv3x2a=(vmepo(3,jp)-vmepo(3,jc))/(y2s(jp)-y2s(jc))
      dv3mx2=(dv3x2a+dv3x2c)*0.5
      dq2x3=(dq2x3c+dq2x3a+dq2x3b+dq2x3d)*0.25
      dq3x2=(dq3x2c+dq3x2a+dq3x2b+dq3x2d)*0.25
      dqc=-(dq2x3-(dq3x2-dv3mx2))
      dqcp=-(dq2x3-(dq3x2-dv3mx2*(1-itot)))
      q1p=dqc
      q1pp=dqcp
      s23=(dq2x3+(dq3x2-dv3mx2))*0.5
      s23p=(dq2x3+(dq3x2-dv3mx2*(1-itot)))*0.5
c
c   normal vorticity
c
      dq3x1c=(q3(ic,jc,kc)-q3(im,jc,kc))*dx1
      dq1x3c=(q1(ic,jc,kc)-q1(ic,jc,km))*dx3
      dq3x1a=(q3(ip,jc,kc)-q3(ic,jc,kc))*dx1
      dq1x3a=(q1(ip,jc,kc)-q1(ip,jc,km))*dx3
      dq3x1b=(q3(ip,jc,kp)-q3(ic,jc,kp))*dx1
      dq1x3b=(q1(ip,jc,kp)-q1(ip,jc,kc))*dx3
      dq3x1d=(q3(ic,jc,kp)-q3(im,jc,kp))*dx1
      dq1x3d=(q1(ic,jc,kp)-q1(ic,jc,kc))*dx3
      dq1x3=(dq1x3c+dq1x3a+dq1x3b+dq1x3d)*0.25
      dq3x1=(dq3x1c+dq3x1a+dq3x1b+dq3x1d)*0.25
      dphc=-(dq3x1-dq1x3)
      q2p=dphc
      s13=(dq1x3+dq3x1)*0.5
c
c   streamwise vorticity
c
      dq1x2c=(q1(ic,jc,kc)-q1(ic,jm,kc))/(y2s(jc)-y2s(jm))
      dq2x1c=(q2(ic,jc,kc)-q2(im,jc,kc))*dx1
      dq1x2a=(q1(ic,jp,kc)-q1(ic,jc,kc))/(y2s(jp)-y2s(jc))
      dq2x1a=(q2(ic,jp,kc)-q2(im,jp,kc))*dx1
      dq1x2d=(q1(ip,jc,kc)-q1(ip,jm,kc))/(y2s(jc)-y2s(jm))
      dq2x1d=(q2(ip,jc,kc)-q2(ic,jc,kc))*dx1
      dq1x2b=(q1(ip,jp,kc)-q1(ip,jc,kc))/(y2s(jp)-y2s(jc))
      dq2x1b=(q2(ip,jp,kc)-q2(ic,jp,kc))*dx1
      dv1x2c=(vmepo(1,jc)-vmepo(1,jm))/(y2s(jc)-y2s(jm))
      dv1x2a=(vmepo(1,jp)-vmepo(1,jc))/(y2s(jp)-y2s(jc))
      dv1mx2=(dv1x2a+dv1x2c)*0.5
      dq1x2=(dq1x2c+dq1x2a+dq1x2b+dq1x2d)*0.25
      dq2x1=(dq2x1c+dq2x1a+dq2x1b+dq2x1d)*0.25
      qcapc=-(dq1x2-dv1mx2-dq2x1)
      qcapcp=-(dq1x2-dv1mx2*(1-itot)-dq2x1)
      q3p=qcapc
      q3pp=qcapcp
      s12=(dq1x2-dv1mx2+dq2x1)*0.5
      s12p=(dq1x2-dv1mx2*(1-itot)+dq2x1)*0.5
      dq1x1=(q1(ip,jc,kc)-q1(ic,jc,kc))*dx1
      dq3x3=(q3(ic,jc,kp)-q3(ic,jc,kc))*dx3
      dv2mx2=(vmepo(2,jp)-vmepo(2,jc))/(yp2(jp)-yp2(jc))
      dq2x2=(q2(ic,jp,kc)-q2(ic,jc,kc))/(yp2(jp)-yp2(jc))
      s11=dq1x1
      s22=dq2x2-dv2mx2
      s22p=dq2x2-dv2mx2*(1-itot)
      s33=dq3x3
      dislo=2.*(s11**2+s22**2+s33**2+
     1      2.*(s12**2+s13**2+s23**2))
      enspro=q1p**2*s11+q2p**2*s22+q3p**2*s33+
     1      2.*(q1p*q2p*s12+q2p*q3p*s23+q1p*q3p*s13)
                                  else
      if(jc.eq.1) then
c
c   fluctuating vorticities at the cell centre  wall
c
c
c   spanwise vorticity
c
      dq2x3c=(q2(ic,jc,kc)-q2(ic,jc,km))*dx3
      dq3x2c=(q3(ic,jc,kc))/(y2s(jc)-yp2(jc))
      dq2x3a=(q2(ic,jp,kc)-q2(ic,jp,km))*dx3
      dq3x2a=(q3(ic,jp,kc)-q3(ic,jc,kc))/(y2s(jp)-y2s(jc))
      dq2x3d=(q2(ic,jc,kp)-q2(ic,jc,kc))*dx3
      dq3x2d=(q3(ic,jc,kp))/(y2s(jc)-yp2(jc))
      dq2x3b=(q2(ic,jp,kp)-q2(ic,jp,kc))*dx3
      dq3x2b=(q3(ic,jp,kp)-q3(ic,jc,kp))/(y2s(jp)-y2s(jc))
      dq2x3=(dq2x3c+dq2x3a+dq2x3b+dq2x3d)*0.25
      dv3x2c=(vmepo(3,jc))/(y2s(jc)-yp2(jc))
      dv3x2a=(vmepo(3,jp)-vmepo(3,jc))/(y2s(jp)-y2s(jc))
      dv3mx2=(dv3x2a+dv3x2c)*0.5
      dq3x2=(dq3x2c+dq3x2a+dq3x2b+dq3x2d)*0.25
      dqc=-(dq2x3-(dq3x2-dv3mx2))
      dqcp=-(dq2x3-(dq3x2-dv3mx2*(1-itot)))
      q1p=dqc
      q1pp=dqcp
      s23=(dq2x3+(dq3x2-dv3mx2))*0.5
      s23p=(dq2x3+(dq3x2-dv3mx2*(1-itot)))*0.5
c
c   normal vorticity
c
      dq3x1c=(q3(ic,jc,kc)-q3(im,jc,kc))*dx1
      dq1x3c=(q1(ic,jc,kc)-q1(ic,jc,km))*dx3
      dq3x1a=(q3(ip,jc,kc)-q3(ic,jc,kc))*dx1
      dq1x3a=(q1(ip,jc,kc)-q1(ip,jc,km))*dx3
      dq3x1b=(q3(ip,jc,kp)-q3(ic,jc,kp))*dx1
      dq1x3b=(q1(ip,jc,kp)-q1(ip,jc,kc))*dx3
      dq3x1d=(q3(ic,jc,kp)-q3(im,jc,kp))*dx1
      dq1x3d=(q1(ic,jc,kp)-q1(ic,jc,kc))*dx3
      dq1x3=(dq1x3c+dq1x3a+dq1x3b+dq1x3d)*0.25
      dq3x1=(dq3x1c+dq3x1a+dq3x1b+dq3x1d)*0.25
      dphc=-(dq3x1-dq1x3)
      q2p=dphc
      s13=(dq1x3+dq3x1)*0.5
c
c   streamwise vorticity
c
      dq1x2c=(q1(ic,jc,kc))/(y2s(jc)-yp2(jc))
      dq2x1c=(q2(ic,jc,kc)-q2(im,jc,kc))*dx1
      dq1x2a=(q1(ic,jp,kc)-q1(ic,jc,kc))/(y2s(jp)-y2s(jc))
      dq2x1a=(q2(ic,jp,kc)-q2(im,jp,kc))*dx1
      dq1x2d=(q1(ip,jc,kc))/(y2s(jc)-yp2(jc))
      dq2x1d=(q2(ip,jc,kc)-q2(ic,jc,kc))*dx1
      dq1x2b=(q1(ip,jp,kc)-q1(ip,jc,kc))/(y2s(jp)-y2s(jc))
      dq2x1b=(q2(ip,jp,kc)-q2(ic,jp,kc))*dx1
      dv1x2c=(vmepo(1,jc))/(y2s(jc)-yp2(jc))
      dv1x2a=(vmepo(1,jp)-vmepo(1,jc))/(y2s(jp)-y2s(jc))
      dv1mx2=(dv1x2a+dv1x2c)*0.5
      dq1x2=(dq1x2c+dq1x2a+dq1x2b+dq1x2d)*0.25
      dq2x1=(dq2x1c+dq2x1a+dq2x1b+dq2x1d)*0.25
      qcapc=-(dq1x2-dv1mx2-dq2x1)
      qcapcp=-(dq1x2-dv1mx2*(1-itot)-dq2x1)
      q3p=qcapc
      q3pp=qcapcp
      s12=(dq1x2-dv1mx2+dq2x1)*0.5
      s12p=(dq1x2-dv1mx2*(1-itot)+dq2x1)*0.5
      dq1x1=(q1(ip,jc,kc)-q1(ic,jc,kc))*dx1
      dq3x3=(q3(ic,jc,kp)-q3(ic,jc,kc))*dx3
      dv2mx2=(vmepo(2,jp)-vmepo(2,jc))/(yp2(jp)-yp2(jc))
      dq2x2=(q2(ic,jp,kc)-q2(ic,jc,kc))/(yp2(jp)-yp2(jc))-dv2mx2
      s11=dq1x1
      s22=dq2x2-dv2mx2
      s22p=dq2x2-dv2mx2*(1-itot)
      s33=dq3x3
      dislo=2.*(s11**2+s22**2+s33**2+
     1      2.*(s12**2+s13**2+s23**2))
      enspro=q1p**2*s11+q2p**2*s22+q3p**2*s33+
     1      2.*(q1p*q2p*s12+q2p*q3p*s23+q1p*q3p*s13)
      fr1wlw=dq1x2
      fr3wlw=dq3x2
      dpx1c=(pr(ic,jc,kc)-pr(imv(ic),jc,kc))*dx1
      dpx3c=(pr(ic,jc,kc)-pr(ic,jc,kmv(kc)))*dx3
      if(iquapd.eq.1) then
      nfile=70
      write(nfile)q1pp,q3pp,s12p,s23p,s13
     1           ,pr(ic,jc,kc),dpx1c,dpx3c
                      endif
                  endif
      if(jc.eq.n2m) then
c
c   fluctuating vorticities at the cell centre  wall
c
c
c   spanwise vorticity
c
      dq2x3c=(q2(ic,jc,kc)-q2(ic,jc,km))*dx3
      dq3x2c=(q3(ic,jc,kc)-q3(ic,jm,kc))/(y2s(jc)-y2s(jm))
      dq2x3a=(q2(ic,jp,kc)-q2(ic,jp,km))*dx3
      dq3x2a=(-q3(ic,jc,kc))/(yp2(jp)-y2s(jc))
      dq2x3d=(q2(ic,jc,kp)-q2(ic,jc,kc))*dx3
      dq3x2d=(q3(ic,jc,kp)-q3(ic,jm,kp))/(y2s(jc)-y2s(jm))
      dq2x3b=(q2(ic,jp,kp)-q2(ic,jp,kc))*dx3
      dq3x2b=(-q3(ic,jc,kp))/(yp2(jp)-y2s(jc))
      dv3x2c=(vmepo(3,jc)-vmepo(3,jm))/(y2s(jc)-y2s(jm))
      dv3x2a=(-vmepo(3,jc))/(yp2(jp)-y2s(jc))
      dv3mx2=(dv3x2a+dv3x2c)*0.5
      dq2x3=(dq2x3c+dq2x3a+dq2x3b+dq2x3d)*0.25
      dq3x2=(dq3x2c+dq3x2a+dq3x2b+dq3x2d)*0.25
      dqc=-(dq2x3-(dq3x2-dv3mx2))
      dqcp=-(dq2x3-(dq3x2-dv3mx2*(1-itot)))
      q1p=dqc
      q1pp=dqcp
      s23=(dq2x3+(dq3x2-dv3mx2))*0.5
      s23p=(dq2x3+(dq3x2-dv3mx2*(1-itot)))*0.5
c
c   normal vorticity
c
      dq3x1c=(q3(ic,jc,kc)-q3(im,jc,kc))*dx1
      dq1x3c=(q1(ic,jc,kc)-q1(ic,jc,km))*dx3
      dq3x1a=(q3(ip,jc,kc)-q3(ic,jc,kc))*dx1
      dq1x3a=(q1(ip,jc,kc)-q1(ip,jc,km))*dx3
      dq3x1b=(q3(ip,jc,kp)-q3(ic,jc,kp))*dx1
      dq1x3b=(q1(ip,jc,kp)-q1(ip,jc,kc))*dx3
      dq3x1d=(q3(ic,jc,kp)-q3(im,jc,kp))*dx1
      dq1x3d=(q1(ic,jc,kp)-q1(ic,jc,kc))*dx3
      dq1x3=(dq1x3c+dq1x3a+dq1x3b+dq1x3d)*0.25
      dq3x1=(dq3x1c+dq3x1a+dq3x1b+dq3x1d)*0.25
      dphc=-(dq3x1-dq1x3)
      q2p=dphc
      s13=(dq1x3+dq3x1)*0.5
c
c   streamwise vorticity
c
      dq1x2c=(q1(ic,jc,kc)-q1(ic,jm,kc))/(y2s(jc)-y2s(jm))
      dq2x1c=(q2(ic,jc,kc)-q2(im,jc,kc))*dx1
      dq1x2a=(-q1(ic,jc,kc))/(yp2(jp)-y2s(jc))
      dq2x1a=(q2(ic,jp,kc)-q2(im,jp,kc))*dx1
      dq1x2d=(q1(ip,jc,kc)-q1(ip,jm,kc))/(y2s(jc)-y2s(jm))
      dq2x1d=(q2(ip,jc,kc)-q2(ic,jc,kc))*dx1
      dq1x2b=(-q1(ip,jc,kc))/(yp2(jp)-y2s(jc))
      dq2x1b=(q2(ip,jp,kc)-q2(ic,jp,kc))*dx1
      dv1x2c=(vmepo(1,jc)-vmepo(1,jm))/(y2s(jc)-y2s(jm))
      dv1x2a=(-vmepo(1,jc))/(yp2(jp)-y2s(jc))
      dv1mx2=(dv1x2a+dv1x2c)*0.5
      dq1x2=(dq1x2c+dq1x2a+dq1x2b+dq1x2d)*0.25
      dq2x1=(dq2x1c+dq2x1a+dq2x1b+dq2x1d)*0.25
      qcapc=-(dq1x2-dv1mx2-dq2x1)
      qcapcp=-(dq1x2-dv1mx2*(1-itot)-dq2x1)
      q3p=qcapc
      q3pp=qcapcp
      s12=(dq1x2-dv1mx2+dq2x1)*0.5
      s12p=(dq1x2-dv1mx2*(1-itot)+dq2x1)*0.5
      dq1x1=(q1(ip,jc,kc)-q1(ic,jc,kc))*dx1
      dq3x3=(q3(ic,jc,kp)-q3(ic,jc,kc))*dx3
      dv2mx2=(vmepo(2,jp)-vmepo(2,jc))/(yp2(jp)-yp2(jc))
      dq2x2=(q2(ic,jp,kc)-q2(ic,jc,kc))/(yp2(jp)-yp2(jc))
      s11=dq1x1
      s22=dq2x2-dv2mx2
      s22p=dq2x2-dv2mx2*(1-itot)
      s33=dq3x3
      dislo=2.*(s11**2+s22**2+s33**2+
     1      2.*(s12**2+s13**2+s23**2))
      enspro=q1p**2*s11+q2p**2*s22+q3p**2*s33+
     1      2.*(q1p*q2p*s12+q2p*q3p*s23+q1p*q3p*s13)
      fr1wup=dq1x2
      fr3wup=dq3x2
      dpx1c=(pr(ic,jc,kc)-pr(imv(ic),jc,kc))*dx1
      dpx3c=(pr(ic,jc,kc)-pr(ic,jc,kmv(kc)))*dx3
      if(iquapd.eq.1) then
      nfile=69
      write(nfile)q1pp,q3pp,s12p,s23p,s13
     1           ,pr(ic,jc,kc),dpx1c,dpx3c
                      endif
                    endif
                                  endif
      rhs(ic,jc,kc)=q3p
c
c   normal vortity correlation
c
      vorstr(1,jc)=vorstr(1,jc)+q1p**2
      vorstr(2,jc)=vorstr(2,jc)+q2p**2
      vorstr(3,jc)=vorstr(3,jc)+q3p**2
c
c   other vortity correlation
c
      vorstr(4,jc)=vorstr(4,jc)-q1p*q2p
      vorstr(5,jc)=vorstr(5,jc)-q3p*q2p
      vorstr(6,jc)=vorstr(6,jc)-q3p*q1p
      ensjp=dislo*volz(jc)
      ensj=ensj+ensjp
      espjp=enspro*volz(jc)
      espj=espj+espjp
c
c   velocity vortity correlation
c
      vrcen=(q2(ic,jc,kc)+q2(ic,jp,kc))*0.5
      vtcen=(q1(ic,jc,kc)+q1(ipv(ic),jc,kc) )*0.5
      vzcen=(q3(ic,jc,kc)+q3(ic,jc,kpv(kc)) )*0.5 
      v1p=vtcen-vmed(1,jc)
      v2p=vrcen-vmed(2,jc)
      v3p=vzcen-vmed(3,jc)
      v1pp=vtcen-vmed(1,jc)*(1-itot)
      v2pp=vrcen-vmed(2,jc)*(1-itot)
      v3pp=vzcen-vmed(3,jc)*(1-itot)
      voplo=sqrt(q1p**2+q2p**2+q3p**2)
      veplo=sqrt(v1p**2+v2p**2+v3p**2)
      denhe=veplo*voplo
c
c   helicity density  
c
      hel1=v1p*q1p
      hel2=v2p*q2p
      hel3=v3p*q3p
      dhel1=v1p*q1p/denhe
      dhel2=v2p*q2p/denhe
      dhel3=v3p*q3p/denhe
      hede=dhel1+dhel2+dhel3
      vdtomt(1,jc)=vdtomt(1,jc)+hel1
      vdtomt(2,jc)=vdtomt(2,jc)+hel2
      vdtomt(3,jc)=vdtomt(3,jc)+hel3
c
c   density   of the Q X OMEGA
c
      vcromt(1,jc)=vcromt(1,jc)+v2p*q1p
      vcromt(2,jc)=vcromt(2,jc)+v1p*q2p
      vcromt(3,jc)=vcromt(3,jc)+v3p*q2p
      vcromt(4,jc)=vcromt(4,jc)+v2p*q3p
      vcromt(5,jc)=vcromt(5,jc)+v1p*q3p
      vcromt(6,jc)=vcromt(6,jc)+v3p*q1p
      vomx=(v2p*q1p-v1p*q2p)/denhe
      vomr=(v1p*q3p-v3p*q1p)/denhe
      vomt=(v3p*q2p-v2p*q3p)/denhe
c
c  gradient  pressure velocity correlationns
c
      dpvthc=((phtot(ip,jc,kc)+phtot(ic,jc,kc))*q1(ip,jc,kc)
     1      -(phtot(im,jc,kc)+phtot(ic,jc,kc))*q1(ic,jc,kc))
     1      *0.5*dx1
      dpvaxc=((phtot(ic,jc,kp)+phtot(ic,jc,kc))*q3(ic,jc,kp)
     1      -(phtot(ic,jc,km)+phtot(ic,jc,kc))*q3(ic,jc,kc))
     1      *0.5*dx3
      if(jc.ge.2.and.jc.le.n2m-1) then
      dpvrac=((phtot(ic,jp,kc)+phtot(ic,jc,kc))*q2(ic,jp,kc)  
     1      -(phtot(ic,jm,kc)+phtot(ic,jc,kc))*q2(ic,jc,kc))
     1      *0.5*dx2/caj(jc)
                                  else
      if(jc.eq.1) then
c
c  gradient  pressure velocity correlationns wall
c
      dpvrac=(phtot(ic,jp,kc)+phtot(ic,jc,kc))*q2(ic,jp,kc)  
     1      *0.5*dx2/caj(jc)
                  endif
      if(jc.eq.n2m) then
c
c  gradient  pressure velocity correlationns wall
c
      dpvrac=-(phtot(ic,jm,kc)+phtot(ic,jc,kc))*q2(ic,jc,kc)
     1      *0.5*dx2/caj(jc)
                   endif
                                  endif
      prgrd(1,jc)=prgrd(1,jc)+dpvthc
      prgrd(2,jc)=prgrd(2,jc)+dpvrac
      prgrd(3,jc)=prgrd(3,jc)+dpvaxc
      enslo=voplo**2
      dissj(jc)=dissj(jc)+dislo*avgn
      ensy(jc)=ensy(jc)+enslo*avgn
      enspr(jc)=enspr(jc)+enspro*avgn
      plap=enslo/4.-dislo
c
c   flatness and skewness of enstr, dissipation rate and strains
c   (quasq rms)  (quacu skewnes)  (quafo flatness)
      quasq(1,jc)=q1p**2+quasq(1,jc)
      quasq(2,jc)=q2p**2+quasq(2,jc)
      quasq(3,jc)=q3p**2+quasq(3,jc)
      quasq(4,jc)=s11**2+quasq(4,jc)
      quasq(5,jc)=s22**2+quasq(5,jc)
      quasq(6,jc)=s33**2+quasq(6,jc)
      quasq(7,jc)=s12**2+quasq(7,jc)
      quasq(8,jc)=s23**2+quasq(8,jc)
      quasq(9,jc)=s13**2+quasq(9,jc)
      quasq(10,jc)=(dislo-quaav(10,jc))**2+quasq(10,jc)
      quasq(11,jc)=(enslo-quaav(11,jc))**2+quasq(11,jc)
      quasq(12,jc)=(plap-quaav(12,jc))**2+quasq(12,jc)
      quacu(1,jc)=q1p**3+quacu(1,jc)
      quacu(2,jc)=q2p**3+quacu(2,jc)
      quacu(3,jc)=q3p**3+quacu(3,jc)
      quacu(4,jc)=s11**3+quacu(4,jc)
      quacu(5,jc)=s22**3+quacu(5,jc)
      quacu(6,jc)=s33**3+quacu(6,jc)
      quacu(7,jc)=s12**3+quacu(7,jc)
      quacu(8,jc)=s23**3+quacu(8,jc)
      quacu(9,jc)=s13**3+quacu(9,jc)
      quacu(10,jc)=(dislo-quaav(10,jc))**3+quacu(10,jc)
      quacu(11,jc)=(enslo-quaav(11,jc))**3+quacu(11,jc)
      quacu(12,jc)=(plap-quaav(12,jc))**3+quacu(12,jc)
      quafo(1,jc)=q1p**4+quafo(1,jc)
      quafo(2,jc)=q2p**4+quafo(2,jc)
      quafo(3,jc)=q3p**4+quafo(3,jc)
      quafo(4,jc)=s11**4+quafo(4,jc)
      quafo(5,jc)=s22**4+quafo(5,jc)
      quafo(6,jc)=s33**4+quafo(6,jc)
      quafo(7,jc)=s12**4+quafo(7,jc)
      quafo(8,jc)=s23**4+quafo(8,jc)
      quafo(9,jc)=s13**4+quafo(9,jc)
      quafo(10,jc)=(dislo-quaav(10,jc))**4+quafo(10,jc)
      quafo(11,jc)=(enslo-quaav(11,jc))**4+quafo(11,jc)
      quafo(12,jc)=(plap-quaav(12,jc))**4+quafo(12,jc)
      if(jc.eq.jprq(npq)) then
c
c    Q X omega 
c
      if(iquapd.eq.1) then
      nfile=71+npq
      write(nfile)v1pp,v2pp,v3pp,q1pp,q2p,q3pp,
     1          s11,s22p,s33,s12p,s23p,s13
     1           ,phtot(ic,jc,kc)
                      endif
      enssma(npq)=max(enslo,enssma(npq))
      dissma(npq)=max(dislo,dissma(npq))
      esprma(npq)=max(abs(enspro),esprma(npq))
      espma(npq)=max(enspro,espma(npq))
      espmi(npq)=min(enspro,espmi(npq))
                          endif
      enddo
      enddo
      if(jc.eq.jprq(npq)) then
      npq=npq+1
                        endif
      do l=1,3
      vorstr(l,jc)=sqrt( vorstr(l,jc)*avgn )
      vdtomt(l,jc)=vdtomt(l,jc)*avgn
      prgrd(l,jc)=prgrd(l,jc)*avgn
      enddo
      vorstr(4,jc)= vorstr(4,jc)*avgn 
      vorstr(5,jc)= vorstr(5,jc)*avgn 
      vorstr(6,jc)= vorstr(6,jc)*avgn 
      do l=1,6
      vcromt(l,jc)=vcromt(l,jc)*avgn
      enddo
          do l=1,12
      quafo(l,jc)=quafo(l,jc)*avgn
      quacu(l,jc)=quacu(l,jc)*avgn
      quasq(l,jc)=quasq(l,jc)*avgn
          enddo
          do l=1,12
      sququ(l,jc)=sqrt(quasq(l,jc))
      flaqu(l,jc)=quafo(l,jc)/sququ(l,jc)**4 
      skequ(l,jc)=quacu(l,jc)/sququ(l,jc)**3
          enddo
      enddo
  174 format(i4,2x,12e10.3)
      write(74,*)'  field at ntime=',ntime,' rms'
      do jc=1,n2m
      write(74,174)jc,(quasq(l,jc),l=1,12)
      enddo
      write(74,*)'  field at ntime=',ntime,' skew'
      do jc=1,n2m
      write(74,174)jc,(quacu(l,jc),l=1,12)
      enddo
      write(74,*)'  field at ntime=',ntime,' flat'
      do jc=1,n2m
      write(74,174)jc,(quafo(l,jc),l=1,12)
      enddo
      return                                                            
      end                                                               
c
c  ******************** subrout taver *********************
c  summation to perform the time averages
c
      subroutine taver
c
      include 'param.f'
c
      dissto=disst+dissto
      enavo=enav+enavo
      do  l=1,3
      do  j=1,n2
      vpmeo(l,j)=(vpmeo(l,j)+vmepo(l,j))
      vpome(l,j)=(vpome(l,j)+vompo(l,j))
      enddo
      enddo
      do  l=1,3
      do  j=1,n2m
      vmeo(l,j)=(vmeo(l,j)+vmed(l,j))
      voro(l,j)=(voro(l,j)+vorv(l,j))
      vdtomo(l,j)=(vdtomo(l,j)+vdtomt(l,j))
      prgro(l,j)=(prgro(l,j)+prgrd(l,j))
      enddo
      enddo
      do  l=1,4
      do  j=1,n2m
      skeo(l,j)=(skeo(l,j)+ske(l,j))
      flao(l,j)=(flao(l,j)+fla(l,j))
      pvmo(l,j)=(pvmo(l,j)+pvc(l,j))
      enddo
      enddo
          do l=1,12
      do  j=1,n2m
      sququo(l,j)=sququo(l,j)+sququ(l,j)
      flaquo(l,j)=flaquo(l,j)+flaqu(l,j)
      skequo(l,j)=skequo(l,j)+skequ(l,j)
      enddo
          enddo
      do  l=1,6
      do  j=1,n2m
      tursto(l,j)=(tursto(l,j)+turstr(l,j))
      vorsto(l,j)=(vorsto(l,j)+vorstr(l,j))
      stmedo(l,j)=(stmedo(l,j)+stmed(l,j))
      vcromo(l,j)=(vcromo(l,j)+vcromt(l,j))
      enddo
      enddo
      do  j=1,n2m
      pmeo(j)=(pmeo(j)+pmed(j))
      dissjo(j)=(dissjo(j)+dissj(j))
      ensyo(j)=(ensyo(j)+ensy(j))
      enspo(j)=(enspo(j)+enspr(j))
      enddo
      return
      end
c
c  ******************** subrout vetrve *********************
c
      subroutine vetrve
c
      include 'param.f'
c
      do  l=1,5
      do  j=1,n2m
      vtvo(l,j)=(vtvo(l,j)+vtv(l,j))
      enddo
      enddo
      do  j=1,n2m
      dvto(1,j)=(dvto(1,j)
     1          +turstr(1,j)**2+turstr(2,j)**2+turstr(3,j)**2)
      dvto(2,j)=(dvto(2,j)+turstr(4,j))
      dvto(3,j)=(dvto(3,j)+turstr(5,j))
      dvto(4,j)=(dvto(4,j)+turstr(2,j)**2)
      dvto(5,j)=(dvto(5,j)+pvc(2,j))
      enddo
      return
      end
