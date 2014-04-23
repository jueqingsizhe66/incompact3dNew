c                                                                       
c  ************************* subrout,velc  **********************  
c                                                                       
c     this subroutine calculates the mean velocity at the centre of the cell
c     and at the j location where are defined
c                                                                       
      subroutine velc                    
      include 'param.f'
      avgn=1./(float(n1m*n3m))
      vl123=1./(n1m*n2m*n3m)
      volto=1./(alx1*alx3*2.)
c
c  velocities at the center of the cell  vmed(
c  velocities at the location   vmepo(
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
c   mean pressure and total velocities
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
c     this subroutine calculates the correlations of
c     the fluctuating velocities
c                                                                       
      subroutine veltur(ntime)
      include 'param.f'
      common/phto/phtot(m1,m2,m3)
      common/prcost/pcost
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
      vtv(l,jc)=0.
      ske(l,jc)=0.
      fla(l,jc)=0.
      pvc(l,jc)=0.
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
c   triple velocity with u'_2
c
      vtv(1,jc)=q2p*q1p**2+vtv(1,jc)
      vtv(2,jc)=q2p*q2p**2+vtv(2,jc)
      vtv(3,jc)=q2p*q3p**2+vtv(3,jc)
      vtv(4,jc)=q3p*q2p**2+vtv(4,jc)
      enddo
      enddo
      do l=1,3
      turstr(l,jc)=sqrt( turstr(l,jc)*avgn )
      enddo
      do l=4,6
      turstr(l,jc)= turstr(l,jc)*avgn 
      enddo
      do l=1,4
      vtv(l,jc)=vtv(l,jc)*avgn
      ske(l,jc)=ske(l,jc)*avgn 
      fla(l,jc)=fla(l,jc)*avgn 
      pvc(l,jc)=pvc(l,jc)*avgn 
      enddo
      enddo
      enav=enej
  612 format(1x,2e15.5)
      return                                                            
      end                                                               
c
c                                                                       
c  ************************ subrout vorc  **********************  
c                                                                       
c     this subroutine calculates the mean vorticity               
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
c  ************************ subrout vortur  **********************  
c                                                                       
c     this subroutine calculates the fluctuating  vorticity               
c     and the relative correlations
c     everything at the cell centre
c                                                                       
      subroutine vortur(ntime)
      include 'param.f'
      avgn=1./(float(n1m*n3m))
   83 format(i3.3)
  144 format(9e14.6)
  117 format(3x,i5,3x,4e15.6)
      do n=1,mpq
      enssma(n)=0.
      dissma(n)=0.
      esprma(n)=0.
      espma(n)=-100.
      espmi(n)=+100.
      enddo
      npq=1
      do jc=1,n2m
      jp=jc+1
      jm=jc-1
      dissj(jc)=0.
      ensy(jc)=0.
      enspr(jc)=0.
            do l=1,9
      stmed(l,jc)=0.
            enddo
            do l=1,3
      vdtomt(l,jc)=0.
            enddo
            do l=1,6
      vorstr(l,jc)=0.
      vcromt(l,jc)=0.
            enddo
      do ne=1,4
      do nt=1,2
      budg(ne,nt,jc)=0.
      enddo
      enddo
      do kc=1,n3m
      kp=kpv(kc)
      km=kmv(kc)
      do ic=1,n1m                                                     
      ip=ipv(ic)
      im=imv(ic)
      v2cen=(q2(ic,jc,kc)+q2(ic,jp,kc))*0.5
      v1cen=(q1(ic,jc,kc)+q1(ip,jc,kc))*0.5
      v3cen=(q3(ic,jc,kc)+q3(ic,jc,kp))*0.5
      v1p=v1cen-vmed(1,jc)
      v2p=v2cen-vmed(2,jc)
      v3p=v3cen-vmed(3,jc)
      prp=pr(ic,jc,kc)-pmed(jc)
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
      vor1p=dqc
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
      vor2p=dphc
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
      vor3p=qcapc
      s12=(dq1x2-dv1mx2+dq2x1)*0.5
      s12p=(dq1x2-dv1mx2*(1-itot)+dq2x1)*0.5
      dq1x1=(q1(ip,jc,kc)-q1(ic,jc,kc))*dx1
      dq3x3=(q3(ic,jc,kp)-q3(ic,jc,kc))*dx3
      dv2mx2=(vmepo(2,jp)-vmepo(2,jc))/(yp2(jp)-yp2(jc))
      dq2x2=(q2(ic,jp,kc)-q2(ic,jc,kc))/(yp2(jp)-yp2(jc))
      s11=dq1x1
      s22=dq2x2-dv2mx2
      s33=dq3x3
      dq3x2p=dq3x2-dv3mx2
      dq1x2p=dq1x2-dv1mx2
      dq2x2p=dq2x2-dv2mx2
                                  else
      if(jc.eq.1) then
c
c   fluctuating vorticities at the cell centre  near
c   the wall
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
      vor1p=dqc
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
      vor2p=dphc
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
      vor3p=qcapc
      s12=(dq1x2-dv1mx2+dq2x1)*0.5
      dq1x1=(q1(ip,jc,kc)-q1(ic,jc,kc))*dx1
      dq3x3=(q3(ic,jc,kp)-q3(ic,jc,kc))*dx3
      dv2mx2=(vmepo(2,jp)-vmepo(2,jc))/(yp2(jp)-yp2(jc))
      dq2x2=(q2(ic,jp,kc)-q2(ic,jc,kc))/(yp2(jp)-yp2(jc))-dv2mx2
      s11=dq1x1
      s22=dq2x2-dv2mx2
      s33=dq3x3
      dq3x2p=dq3x2-dv3mx2
      dq1x2p=dq1x2-dv1mx2
      dq2x2p=dq2x2-dv2mx2
      dislo=2.*(s11**2+s22**2+s33**2+
     1      2.*(s12**2+s13**2+s23**2))
      enspro=vor1p**2*s11+vor2p**2*s22+vor3p**2*s33+
     1      2.*(vor1p*vor2p*s12+vor2p*vor3p*s23+vor1p*vor3p*s13)
      voplo=vor1p**2+vor2p**2+vor3p**2
      fr1wlw=dq1x2
      fr3wlw=dq3x2
      dpx1c=(pr(ic,jc,kc)-pr(imv(ic),jc,kc))*dx1
      dpx3c=(pr(ic,jc,kc)-pr(ic,jc,kmv(kc)))*dx3
                  endif
      if(jc.eq.n2m) then
c
c   fluctuating vorticities at the cell centre  
c    near the wall
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
      vor1p=dqc
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
      vor2p=dphc
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
      vor3p=qcapc
      s12=(dq1x2-dv1mx2+dq2x1)*0.5
      s12p=(dq1x2-dv1mx2*(1-itot)+dq2x1)*0.5
      dq1x1=(q1(ip,jc,kc)-q1(ic,jc,kc))*dx1
      dq3x3=(q3(ic,jc,kp)-q3(ic,jc,kc))*dx3
      dv2mx2=(vmepo(2,jp)-vmepo(2,jc))/(yp2(jp)-yp2(jc))
      dq2x2=(q2(ic,jp,kc)-q2(ic,jc,kc))/(yp2(jp)-yp2(jc))
      s11=dq1x1
      s22=dq2x2-dv2mx2
      s33=dq3x3
      dq3x2p=dq3x2-dv3mx2
      dq1x2p=dq1x2-dv1mx2
      dq2x2p=dq2x2-dv2mx2
      dislo=2.*(s11**2+s22**2+s33**2+
     1      2.*(s12**2+s13**2+s23**2))
      enspro=vor1p**2*s11+vor2p**2*s22+vor3p**2*s33+
     1      2.*(vor1p*vor2p*s12+vor2p*vor3p*s23+vor1p*vor3p*s13)
      dpx1c=(pr(ic,jc,kc)-pr(imv(ic),jc,kc))*dx1
      dpx3c=(pr(ic,jc,kc)-pr(ic,jc,kmv(kc)))*dx3
                    endif
                                  endif
c
c   for the budget of <u'_3u'_3>
c
      pdu3=2.*dq3x3*prp
      eps33=2.*(dq3x1**2+dq3x2p**2+dq3x3**2)
c
c   for the budget of <u'_1u'_1>
c
      pdu1=2.*dq1x1*prp
      eps11=2.*(dq1x1**2+dq1x2p**2+dq1x3**2)
c
c   for the budget of <u'_2u'_2>
c
      pdu2=2.*dq2x2p*prp
      eps22=2.*(dq2x1**2+dq2x2p**2+dq2x3**2)
c
c   for the budget of <u'_2u'_3>
c
      pdu23=dq2x3*prp+dq3x2p*prp
      eps23=2.*(dq3x1*dq2x1+dq3x2p*dq2x2p+dq3x3*dq2x3)

      rhs(ic,jc,kc)=vor3p
c
c   normal vortity correlation
c
      vorstr(1,jc)=vorstr(1,jc)+vor1p**2
      vorstr(2,jc)=vorstr(2,jc)+vor2p**2
      vorstr(3,jc)=vorstr(3,jc)+vor3p**2
c
c   other vortity correlation
c
      vorstr(4,jc)=vorstr(4,jc)-vor1p*vor2p
      vorstr(5,jc)=vorstr(5,jc)-vor3p*vor2p
      vorstr(6,jc)=vorstr(6,jc)-vor3p*vor1p
      ensjp=dislo*volz(jc)
      ensj=ensj+ensjp
      espjp=enspro*volz(jc)
      espj=espj+espjp
      enslo=voplo
      dissj(jc)=dissj(jc)+dislo*avgn
      ensy(jc)=ensy(jc)+enslo*avgn
      enspr(jc)=enspr(jc)+enspro*avgn
c
c   evaluation of the averages budg(neq,nte,jc
c   neq=3 (nte= 1 press vel)  (=2  turb diff.)   (=3    dissip.)
c
      budg(3,1,jc)=budg(3,1,jc)+pdu3*avgn
      budg(3,2,jc)=budg(3,2,jc)-eps33*avgn
      budg(2,1,jc)=budg(2,1,jc)+pdu2*avgn
      budg(2,2,jc)=budg(2,2,jc)-eps22*avgn
      budg(1,1,jc)=budg(1,1,jc)+pdu1*avgn
      budg(1,2,jc)=budg(1,2,jc)-eps11*avgn
      budg(4,1,jc)=budg(4,1,jc)+pdu23*avgn
      budg(4,2,jc)=budg(4,2,jc)-eps23*avgn

c
c   helicity density
c
      hel1=v1p*vor1p
      hel2=v2p*vor2p
      hel3=v3p*vor3p
      vdtomt(1,jc)=vdtomt(1,jc)+hel1
      vdtomt(2,jc)=vdtomt(2,jc)+hel2
      vdtomt(3,jc)=vdtomt(3,jc)+hel3
c
c   density   of the Q X OMEGA
c
      vcromt(1,jc)=vcromt(1,jc)+v2p*vor1p
      vcromt(2,jc)=vcromt(2,jc)+v1p*vor2p
      vcromt(3,jc)=vcromt(3,jc)+v3p*vor2p
      vcromt(4,jc)=vcromt(4,jc)+v2p*vor3p
      vcromt(5,jc)=vcromt(5,jc)+v1p*vor3p
      vcromt(6,jc)=vcromt(6,jc)+v3p*vor1p
c
c   mean velocity gradient  
c
      stmed(1,jc)=stmed(1,jc)+dq1x1
      stmed(2,jc)=stmed(2,jc)+dq1x2
      stmed(3,jc)=stmed(3,jc)+dq1x3
      stmed(4,jc)=stmed(4,jc)+dq2x1
      stmed(5,jc)=stmed(5,jc)+dq2x2
      stmed(6,jc)=stmed(6,jc)+dq2x3
      stmed(7,jc)=stmed(7,jc)+dq3x1
      stmed(8,jc)=stmed(8,jc)+dq3x2
      stmed(9,jc)=stmed(9,jc)+dq3x3
      enddo
      enddo
      do l=1,3
      vorstr(l,jc)=sqrt( vorstr(l,jc)*avgn )
      vdtomt(l,jc)=vdtomt(l,jc)*avgn
      enddo
      do l=1,9
      stmed(l,jc)=stmed(l,jc)*avgn
      enddo
      do l=1,6
      vcromt(l,jc)=vcromt(l,jc)*avgn
      enddo
      vorstr(4,jc)= vorstr(4,jc)*avgn 
      vorstr(5,jc)= vorstr(5,jc)*avgn 
      vorstr(6,jc)= vorstr(6,jc)*avgn 
      enddo
      return                                                            
      end                                                               
c
c  ******************** subrout taver *********************
c  evaluates the summation of the correlations for each
c  field in the outh routine these quantities are
c  divided by the number of fields
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
      enddo
      enddo
      do  l=1,4
      do  j=1,n2m
      skeo(l,j)=(skeo(l,j)+ske(l,j))
      flao(l,j)=(flao(l,j)+fla(l,j))
      pvmo(l,j)=(pvmo(l,j)+pvc(l,j))
      vtvo(l,j)=(vtvo(l,j)+vtv(l,j))
      enddo
      enddo
      do  l=1,6
      do  j=1,n2m
      tursto(l,j)=(tursto(l,j)+turstr(l,j))
      vorsto(l,j)=(vorsto(l,j)+vorstr(l,j))
      vcromo(l,j)=(vcromo(l,j)+vcromt(l,j))
      enddo
      enddo
      do  l=1,9
      do  j=1,n2m
      stmedo(l,j)=(stmedo(l,j)+stmed(l,j))
      enddo
      enddo
      do  j=1,n2m
      pmeo(j)=(pmeo(j)+pmed(j))
      dissjo(j)=(dissjo(j)+dissj(j))
      ensyo(j)=(ensyo(j)+ensy(j))
      enspo(j)=(enspo(j)+enspr(j))
      enddo
      do ne=1,4
      do nt=1,2
      do  j=1,n2m
      budgo(ne,nt,j)=budgo(ne,nt,j)+budg(ne,nt,j)
      enddo
      enddo
      enddo
      return
      end
