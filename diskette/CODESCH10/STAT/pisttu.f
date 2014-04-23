c                                                                       
c  ************************* subrout,velc  **********************  
c                                                                       
c     this routine calculates the velocity at the cell centre vmed       
c     and at the location where are defined  vmepo
c                                                                       
      subroutine velc                    
      include 'param.f'
      avgn=1./(float(n1m*n3m))
      vl123=1./(rc(n2)*n1m*n2m*n3m)
      volto=1./(pi*alx3d)
c
c   axial velocity at the axis
c
      do kc=1,n3m
      q3ax(kc)=0.
      do ic=1,n1m
      q3ax(kc)=q3ax(kc)+q3(ic,1,kc)
      enddo
      q3ax(kc)=q3ax(kc)/n1m
      enddo
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
      if(jc.eq.1) then
      q2akc= (q2(ic,jp,kc) - q2(isym(ic),jp,kc))*0.5/rc(jp)
      dphc=(q2akc+q2(ic,jp,kc)/rc(jp))*.5
       else
      dphc=
     1  (    q2(ic,jc,kc)/rc(jc)+q2(ic,jp,kc)/rc(jp))*0.5
      q2akc=q2(ic,jc,kc)/rc(jc)
      end if
      dqc=
     1  (    q1(ic,jc,kc)+q1(ipv(ic),jc,kc) )*0.5/rm(jc)
      qcapc=(q3(ic,jc,kc)+q3(ic,jc,kpv(kc)) )*0.5 
      vmed(2,jc)=vmed(2,jc)+dphc
      vmed(1,jc)=vmed(1,jc)+dqc
      vmed(3,jc)=vmed(3,jc)+qcapc
      vmepo(2,jc)=vmepo(2,jc)+q2akc
      vmepo(1,jc)=vmepo(1,jc)+q1(ic,jc,kc)/rm(jc)
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
      prver=pr(ic,jc,kc)+0.5*(
     1            -(ros*0.5*rm(jc))**2)
      pmed(jc)=pmed(jc)+prver
      vit(1)=q1(ic,jc,kc)/rm(jc)*volz(jc)+vit(1)
      vit(2)=(q2(ic,jp,kc)+q2(ic,jc,kc))/rm(jc)*volz(jc)+vit(2)
      vit(3)=q3(ic,jc,kc)*volz(jc)+vit(3)
      vit(4)=pr(ic,jc,kc)*volz(jc)+vit(4)
      enddo
      enddo
      pmed(jc)=pmed(jc)*avgn
      enddo
      vit(1)=vit(1)*volto
      vit(2)=vit(2)*volto
      vit(3)=vit(3)*volto
      vit(4)=vit(4)*volto
      return
      end
c                                                                       
c  ************************* subrout,veltur  **********************  
c                                                                       
c     this subroutine calculates the turbulent stress                  
c     and a large number of high order correlations
c                                                                       
      subroutine veltur(ntime)
      include 'param.f'
      common/phto/phtot(m1,m2,m3)
      common/prcost/pcost
      character*60 filns1,filts1,filns2,filts2,filns3,filts3
      character*2 pnt
      write(pnt,77) ntime
   77 format(i2.2)
      filns1='nstr1.'//pnt
      filns2='nstr2.'//pnt
      filns3='nstr3.'//pnt
      filts1='tstr1.'//pnt
      filts2='tstr2.'//pnt
      filts3='tstr3.'//pnt
      avgn=1./(float(n1m*n3m))
      vl123=1./(rc(n2)*n1m*n2m*n3m)
      volto=1./(pi*alx3d)
  133 format(1x,'veltur',2x,5e10.3)
      enej=0.
      enpsv=0.
      do jc=1,n2m                                                     
      jp=jc+1
      do l=1,4
      ske(l,jc)=0.
      fla(l,jc)=0.
      pvc(l,jc)=0.
      enddo
      do l=1,4
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
      if(jc.eq.1) then
      q2akc= (q2(ic,jp,kc) - q2(isym(ic),jp,kc))*0.5/rc(jp)
      dphc=(q2akc+q2(ic,jp,kc)/rc(jp))*.5
       else
      dphc=
     1  (    q2(ic,jc,kc)/rc(jc)+q2(ic,jp,kc)/rc(jp))*0.5
      end if
      dqc=
     1  (    q1(ic,jc,kc)+q1(ipv(ic),jc,kc) )*0.5/rm(jc)
      qcapc=(q3(ic,jc,kc)+q3(ic,jc,kpv(kc)) )*0.5 
c
c   calculation fluctuating velocities
c
      q1p=dqc-vmed(1,jc)
      q2p=dphc-vmed(2,jc)
      q3p=qcapc-vmed(3,jc)
      prver=pr(ic,jc,kc)+0.5*(
     1            -(ros*0.5*rm(jc))**2)
      pcp=prver-pmed(jc)
      enejp=(q1p**2+q2p**2+q3p**2)*vl123*rm(jc)*g2rm(jc)
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
c   evaluation of mean radial velocity
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
      if(mod(ntime,1).eq.0) then
c
c   CALCULATION OF JOINT  PDF at
c        certain radial positions
c
      npq=1
      do jc=n2m,1,-1
      if(jc.eq.jprq(npq)) then
      do kc=1,n3m
      do ic=1,n1m                                                     
c
c    calculation velocities at  the cell centre
c
      if(jc.eq.1) then
      q2akc= (q2(ic,jp,kc) - q2(isym(ic),jp,kc))*0.5/rc(jp)
      dphc=(q2akc+q2(ic,jp,kc)/rc(jp))*.5
       else
      dphc=
     1  (    q2(ic,jc,kc)/rc(jc)+q2(ic,jp,kc)/rc(jp))*0.5
      end if
      dqc=
     1  (    q1(ic,jc,kc)+q1(ipv(ic),jc,kc) )*0.5/rm(jc)
      qcapc=(q3(ic,jc,kc)+q3(ic,jc,kpv(kc)) )*0.5 
c
c   calculation fluctuating velocities
c
      q1p=dqc-vmed(1,jc)
      q2p=dphc-vmed(2,jc)
      q3p=qcapc-vmed(3,jc)
c
c   joint pdf of v'r v'z
c     abscissa v'z/rms v'z  ordinate -v'r/rms v'r
c   
      asc=q3p/turstr(3,jc)
      ord=-q2p/turstr(2,jc)
      if(abs(asc).gt.asma) asc=asma
      if(abs(ord).gt.orma) ord=orma
      ala=liapdh*(asma+asc)/asma+1.
      lla=ala
      alo=liopdh*(orma+ord)/orma+1.
      llo=alo
      ncpdf(npq,1,lla,llo)=ncpdf(npq,1,lla,llo)+1
c
c   joint pdf of v't v'z
c     abscissa v'z/rms v'z  ordinate v't/rms v't
c   
      ord=q1p/turstr(1,jc)
      if(abs(ord).gt.orma) ord=orma
      alo=liopdh*(orma+ord)/orma+1.
      llo=alo
      ncpdf(npq,2,lla,llo)=ncpdf(npq,2,lla,llo)+1
c
c   joint pdf of v't v'r
c     abscissa v'r/rms v'r  ordinate v't/rms v't
c   
      asc=q2p/turstr(2,jc)
      if(abs(asc).gt.asma) asc=asma
      ala=liapdh*(asma+asc)/asma+1.
      lla=ala
      ncpdf(npq,3,lla,llo)=ncpdf(npq,3,lla,llo)+1
      enddo
      enddo
            npq=npq+1
                           endif
      enddo
                     endif
  612 format(1x,2e15.5)
      return                                                            
      end                                                               
c
c                                                                       
c  ************************ subrout vorc  **********************  
c                                                                       
c     this subroutine calculates the vorticity components
c                                                                       
      subroutine vorc
      include 'param.f'
      avgn=1./(float(n1m*n3m))
c                                                                       
c  ***********  compute the azimuthal vorticity component               
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
      dq2x3=(q2(ic,jc,kc)-q2(ic,jc,km))*dx3/rc(jc)
      dq3x2=(q3(ic,jc,kc)-q3(ic,jm,kc))/(rm(jc)-rm(jm))
      rhs(ic,jc,kc)=(dq2x3-dq3x2)
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
c  At  the axis  
      q2akc= (q2(ic,jp,kc) - q2(isym(ic),jp,kc))*0.5/rc(jp)
      q2akm= (q2(ic,jp,km) - q2(isym(ic),jp,km))*0.5/rc(jp)
      q3axi= (q3(ic,jc,kc) + q3(isym(ic),jc,kc))*0.5
      dq2x3=(q2akc-q2akm)*dx3
      dq3x2=(q3(ic,jc,kc)-q3axi)/(rm(jc)-rc(jc))
      rhs(ic,jc,kc)=dq2x3-dq3x2
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
      dq3x2=(-q3(ic,jm,kc))/(rc(n2)-rm(n2m))
      rhs(ic,jc,kc)=-dq3x2
      vompo(1,jc)=vompo(1,jc)+rhs(ic,jc,kc)
                  enddo
            enddo
      vompo(1,jc)=vompo(1,jc)*avgn
c
c   azimuthal vorticity at the cell centre 
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
c    centerline azimuth. vort
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
c  ***********  compute the radial  vorticity component            
c                           at  i,j+1/2,k
c                                                                  
      do jc=1,n2m
      vompo(2,jc)=0.
            do ic=1,n1m
      im=imv(ic)
                  do kc=1,n3m
      km=kmv(kc)
      dq3x1=(q3(ic,jc,kc)-q3(im,jc,kc))*dx1/rm(jc)
      dq1x3=(q1(ic,jc,kc)-q1(ic,jc,km))*dx3/rm(jc)
      rhs(ic,jc,kc)=(dq3x1-dq1x3)
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
c   vorticity at the cell centre wall        
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
c   vorticity at the cell centre center line
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
c  ***********  compute the axial vorticity component        
c               at         i,j,k+1/2        
c
      do jc=2,n2m                                               
      vompo(3,jc)=0.
      jm=jc-1                                                  
            do kc=1,n3m                                                
                  do ic=1,n1m                                               
      im=imv(ic)                                                  
      dq1x2=(q1(ic,jc,kc)-q1(ic,jm,kc))/(rm(jc)-rm(jm))
      dq2x1=(q2(ic,jc,kc)-q2(im,jc,kc))*dx1/rc(jc)                 
      vorz=(dq1x2-dq2x1)/rc(jc)
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
      dq1x2=-q1(ic,jm,kc)/(rc(jc)-rm(jm))
      rhs(ic,jc,kc)=dq1x2
      vompo(3,jc)=vompo(3,jc)+rhs(ic,jc,kc)
                  enddo
            enddo
      vompo(3,jc)=vompo(3,jc)*avgn
      jc=1
      vompo(3,jc)=0.
            do kc=1,n3m
      vozcm=0.
                  do ic=1,n1m
      vozcm=vozcm+rhs(ic,2,kc)
                  enddo
      do ic=1,n1m
      rhs(ic,jc,kc)=vozcm/n1m
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
c     this subroutine calculates the fluctuating vorticity               
c                                                                       
      subroutine vortur(ntime)
      include 'param.f'
      dimension nppcf(mcqu),npdf(mqu)
      common/phto/phtot(m1,m2,m3)
      avgn=1./(float(n1m*n3m))
      asma1=1.
      orma1=1.
   83 format(i3.3)
  144 format(9e14.6)
      ensj=0.
c
c   calculation of maximum enstrophy   
c
      dissma=0.
      do jc=n2m,2,-1                                                 
      jp=jc+1
      jm=jc-1
      do kc=1,n3m
      kp=kpv(kc)
      km=kmv(kc)
      do ic=1,n1m                                                     
      ip=ipv(ic)
      im=imv(ic)
      if(jc.ge.2.and.jc.le.n2m-1) then
c
c   azimuthal vorticity
c
      dq2x3c=(q2(ic,jc,kc)-q2(ic,jc,km))*dx3/rc(jc)
      dq3x2c=(q3(ic,jc,kc)-q3(ic,jm,kc))/(rm(jc)-rm(jm))
      dq2x3a=(q2(ic,jp,kc)-q2(ic,jp,km))*dx3/rc(jp)
      dq3x2a=(q3(ic,jp,kc)-q3(ic,jc,kc))/(rm(jp)-rm(jc))
      dq2x3d=(q2(ic,jc,kp)-q2(ic,jc,kc))*dx3/rc(jc)
      dq3x2d=(q3(ic,jc,kp)-q3(ic,jm,kp))/(rm(jc)-rm(jm))
      dq2x3b=(q2(ic,jp,kp)-q2(ic,jp,kc))*dx3/rc(jp)
      dq3x2b=(q3(ic,jp,kp)-q3(ic,jc,kp))/(rm(jp)-rm(jc))
      dqc=((dq2x3c-dq3x2c)+(dq2x3a-dq3x2a)
     1    +(dq2x3b-dq3x2b)+(dq2x3d-dq3x2d))*0.25
      q1p=dqc-vorv(1,jc)
c
c   radial vorticity
c
      dq3x1c=(q3(ic,jc,kc)-q3(im,jc,kc))*dx1/rm(jc)
      dq1x3c=(q1(ic,jc,kc)-q1(ic,jc,km))*dx3/rm(jc)
      dq3x1a=(q3(ip,jc,kc)-q3(ic,jc,kc))*dx1/rm(jc)
      dq1x3a=(q1(ip,jc,kc)-q1(ip,jc,km))*dx3/rm(jc)
      dq3x1b=(q3(ip,jc,kp)-q3(ic,jc,kp))*dx1/rm(jc)
      dq1x3b=(q1(ip,jc,kp)-q1(ip,jc,kc))*dx3/rm(jc)
      dq3x1d=(q3(ic,jc,kp)-q3(im,jc,kp))*dx1/rm(jc)
      dq1x3d=(q1(ic,jc,kp)-q1(ic,jc,kc))*dx3/rm(jc)
      dphc=((dq3x1c-dq1x3c)+(dq3x1a-dq1x3a)
     1     +(dq3x1b-dq1x3b)+(dq3x1d-dq1x3d))*0.25
      q2p=dphc-vorv(2,jc)
c
c   axial vorticity
c
      dq1x2c=(q1(ic,jc,kc)-q1(ic,jm,kc))/(rm(jc)-rm(jm))
      dq2x1c=(q2(ic,jc,kc)-q2(im,jc,kc))*dx1/rc(jc)                 
      dq1x2a=(q1(ic,jp,kc)-q1(ic,jc,kc))/(rm(jp)-rm(jc))
      dq2x1a=(q2(ic,jp,kc)-q2(im,jp,kc))*dx1/rc(jp)                 
      dq1x2d=(q1(ip,jc,kc)-q1(ip,jm,kc))/(rm(jc)-rm(jm))
      dq2x1d=(q2(ip,jc,kc)-q2(ic,jc,kc))*dx1/rc(jc)                 
      dq1x2b=(q1(ip,jp,kc)-q1(ip,jc,kc))/(rm(jp)-rm(jc))
      dq2x1b=(q2(ip,jp,kc)-q2(ic,jp,kc))*dx1/rc(jp)                 
      qcapc=((dq1x2c-dq2x1c)/rc(jc)+(dq1x2a-dq2x1a)/rc(jp)
     1      +(dq1x2b-dq2x1b)/rc(jp)+(dq1x2d-dq2x1d)/rc(jc))*0.25
      q3p=qcapc-vorv(3,jc)
                                  else
      if(jc.eq.n2m) then
c
c   fluctuating vorticities at the cell centre  wall
c
      q1p=vowa(1,jc,kc)-vorv(1,jc)
      q2p=vowa(2,jc,kc)-vorv(2,jc)
      q3p=vowa(3,jc,kc)-vorv(3,jc)
                    endif
                                  endif
      enslo=(q1p**2+q2p**2+q3p**2)
      dissma=max(enslo,dissma)
      enddo
      enddo
      enddo
      npq=1
      do jc=n2m,1,-1                                                 
      jp=jc+1
      jm=jc-1
      dissj(jc)=0.
            do l=1,3
      vdtomt(l,jc)=0.
      prgrd(l,jc)=0.
            enddo
            do l=1,6
      vorstr(l,jc)=0.
      vcromt(l,jc)=0.
            enddo
      do kc=1,n3m
      kp=kpv(kc)
      km=kmv(kc)
      do ic=1,n1m                                                     
      ip=ipv(ic)
      im=imv(ic)
      if(jc.ge.2.and.jc.le.n2m-1) then
c
c   azimuthal vorticity
c
      dq2x3c=(q2(ic,jc,kc)-q2(ic,jc,km))*dx3/rc(jc)
      dq3x2c=(q3(ic,jc,kc)-q3(ic,jm,kc))/(rm(jc)-rm(jm))
      dq2x3a=(q2(ic,jp,kc)-q2(ic,jp,km))*dx3/rc(jp)
      dq3x2a=(q3(ic,jp,kc)-q3(ic,jc,kc))/(rm(jp)-rm(jc))
      dq2x3d=(q2(ic,jc,kp)-q2(ic,jc,kc))*dx3/rc(jc)
      dq3x2d=(q3(ic,jc,kp)-q3(ic,jm,kp))/(rm(jc)-rm(jm))
      dq2x3b=(q2(ic,jp,kp)-q2(ic,jp,kc))*dx3/rc(jp)
      dq3x2b=(q3(ic,jp,kp)-q3(ic,jc,kp))/(rm(jp)-rm(jc))
      dqc=((dq2x3c-dq3x2c)+(dq2x3a-dq3x2a)
     1    +(dq2x3b-dq3x2b)+(dq2x3d-dq3x2d))*0.25
      q1p=dqc-vorv(1,jc)
c
c   radial vorticity
c
      dq3x1c=(q3(ic,jc,kc)-q3(im,jc,kc))*dx1/rm(jc)
      dq1x3c=(q1(ic,jc,kc)-q1(ic,jc,km))*dx3/rm(jc)
      dq3x1a=(q3(ip,jc,kc)-q3(ic,jc,kc))*dx1/rm(jc)
      dq1x3a=(q1(ip,jc,kc)-q1(ip,jc,km))*dx3/rm(jc)
      dq3x1b=(q3(ip,jc,kp)-q3(ic,jc,kp))*dx1/rm(jc)
      dq1x3b=(q1(ip,jc,kp)-q1(ip,jc,kc))*dx3/rm(jc)
      dq3x1d=(q3(ic,jc,kp)-q3(im,jc,kp))*dx1/rm(jc)
      dq1x3d=(q1(ic,jc,kp)-q1(ic,jc,kc))*dx3/rm(jc)
      dphc=((dq3x1c-dq1x3c)+(dq3x1a-dq1x3a)
     1     +(dq3x1b-dq1x3b)+(dq3x1d-dq1x3d))*0.25
      q2p=dphc-vorv(2,jc)
c
c   axial vorticity
c
      dq1x2c=(q1(ic,jc,kc)-q1(ic,jm,kc))/(rm(jc)-rm(jm))
      dq2x1c=(q2(ic,jc,kc)-q2(im,jc,kc))*dx1/rc(jc)                 
      dq1x2a=(q1(ic,jp,kc)-q1(ic,jc,kc))/(rm(jp)-rm(jc))
      dq2x1a=(q2(ic,jp,kc)-q2(im,jp,kc))*dx1/rc(jp)                 
      dq1x2d=(q1(ip,jc,kc)-q1(ip,jm,kc))/(rm(jc)-rm(jm))
      dq2x1d=(q2(ip,jc,kc)-q2(ic,jc,kc))*dx1/rc(jc)                 
      dq1x2b=(q1(ip,jp,kc)-q1(ip,jc,kc))/(rm(jp)-rm(jc))
      dq2x1b=(q2(ip,jp,kc)-q2(ic,jp,kc))*dx1/rc(jp)                 
      qcapc=((dq1x2c-dq2x1c)/rc(jc)+(dq1x2a-dq2x1a)/rc(jp)
     1      +(dq1x2b-dq2x1b)/rc(jp)+(dq1x2d-dq2x1d)/rc(jc))*0.25
      q3p=qcapc-vorv(3,jc)
                                  else
      if(jc.eq.1) then
c
c   fluctuating vorticities at the cell centre  centerline
c
      q1p=vocl(1,jc,kc)-vorv(1,jc)
      q2p=vocl(2,jc,kc)-vorv(2,jc)
      q3p=vocl(3,jc,kc)-vorv(3,jc)
                  endif
      if(jc.eq.n2m) then
c
c   fluctuating vorticities at the cell centre  wall
c
      q1p=vowa(1,jc,kc)-vorv(1,jc)
      q2p=vowa(2,jc,kc)-vorv(2,jc)
      q3p=vowa(3,jc,kc)-vorv(3,jc)
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
      enslo=(q1p**2+q2p**2+q3p**2)
      dissma=max(enslo,dissma)
      ensjp=enslo*volz(jc)
      ensj=ensj+ensjp
c
c   velocity vortity correlation
c
      if(jc.eq.1) then
      q2akc= (q2(ic,jp,kc) - q2(isym(ic),jp,kc))*0.5
      vrcen=(q2akc+q2(ic,jp,kc))*.5/rc(jp)
       else
      vrcen=(q2(ic,jc,kc)/rc(jc)+q2(ic,jp,kc)/rc(jp))*0.5
      end if
      vtcen=(q1(ic,jc,kc)+q1(ipv(ic),jc,kc) )*0.5/rm(jc)
      vzcen=(q3(ic,jc,kc)+q3(ic,jc,kpv(kc)) )*0.5 
      v1p=vtcen-vmed(1,jc)
      v2p=vrcen-vmed(2,jc)
      v3p=vzcen-vmed(3,jc)
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
      dpvthc=((phtot(ip,jc,kc)+phtot(ic,jc,kc))*q1(ip,jc,kc)/rm(jc)
     1      -(phtot(im,jc,kc)+phtot(ic,jc,kc))*q1(ic,jc,kc)/rm(jc))
     1      *0.5*dx1/rm(jc)
      dpvaxc=((phtot(ic,jc,kp)+phtot(ic,jc,kc))*q3(ic,jc,kp)
     1      -(phtot(ic,jc,km)+phtot(ic,jc,kc))*q3(ic,jc,kc))
     1      *0.5*dx3
      if(jc.ge.2.and.jc.le.n2m-1) then
      dpvrac=((phtot(ic,jp,kc)+phtot(ic,jc,kc))*q2(ic,jp,kc)  
     1      -(phtot(ic,jm,kc)+phtot(ic,jc,kc))*q2(ic,jc,kc))
     1      *0.5*dx2*rm(jc)/g2rm(jc)
                                  else
      if(jc.eq.1) then
c
c  gradient  pressure velocity correlationns centerline
c
      dpvrac=(phtot(ic,jp,kc)+phtot(ic,jc,kc))*q2(ic,jp,kc)  
     1      *0.5*dx2*rm(jc)/g2rm(jc)
                  endif
      if(jc.eq.n2m) then
c
c  gradient  pressure velocity correlationns wall
c
      dpvrac=-(phtot(ic,jm,kc)+phtot(ic,jc,kc))*q2(ic,jc,kc)
     1      *0.5*dx2*rm(jc)/g2rm(jc)
                   endif
                                  endif
      prgrd(1,jc)=prgrd(1,jc)+dpvthc
      prgrd(2,jc)=prgrd(2,jc)+dpvrac
      prgrd(3,jc)=prgrd(3,jc)+dpvaxc

      dissj(jc)=dissj(jc)+enslo*avgn
      if(jc.eq.jprq(npq)) then
      qpdf(1)=hede
c
c    Q X omega 
c
      qpdf(4)=vomx
      qpdf(2)=vomt
      qpdf(3)=vomr
      qpdf(5)=dhel1
      qpdf(6)=dhel2
      qpdf(7)=dhel3
      qpdf(8)=enslo
c
c   calculation of modified pressure gradients
c
      dpth=(phtot(ip,jc,kc)-phtot(im,jc,kc))*0.5*dx1/rm(jc)
      dpax=(phtot(ic,jc,kp)-phtot(ic,jc,km))*0.5*dx3
      if(jc.ge.2.and.jc.le.n2m-1) then
      upd2(jc)=dx2*rc(jc)/g2rc(jc)
      upd2(jp)=dx2*rc(jp)/g2rc(jp)
      dprrp=(phtot(ic,jp,kc)-phtot(ic,jc,kc))
      v2dvgp=dprrp*upd2(jp)
      dprrm=(phtot(ic,jc,kc)-phtot(ic,jm,kc))
      v2dvgm=dprrm*upd2(jc)
      dpra=(v2dvgp+v2dvgm)*0.5
                                  else
      if(jc.eq.1) then
c
c   pressure gradient at the cell centre  centerline
c
      upd2(jp)=dx2*rc(jp)/g2rc(jp)
      dprrp=(phtot(ic,jp,kc)-phtot(ic,jc,kc))
      v2dvgp=dprrp*upd2(jp)
      v2dvgm=0.
      dpra=(v2dvgp+v2dvgm)*0.5
                  endif
      if(jc.eq.n2m) then
c
c   pressure gradient at the cell centre  wall
c
      upd2(jc)=dx2*rc(jc)/g2rc(jc)
      v2dvgp=0.
      dprrm=(phtot(ic,jc,kc)-phtot(ic,jm,kc))
      v2dvgm=dprrm*upd2(jc)
      dpra=(v2dvgp+v2dvgm)*0.5
                    endif
                                  endif

c
c     -grad(p+q^2/2)
c
      qpdf(11)=-dpax/denhe
      qpdf(9)=-dpth/denhe
      qpdf(10)=-dpra/denhe
c
c     QX OMEGA -grad(p+q^2/2)
c
      qpdf(14)=vomx-dpax/denhe
      qpdf(12)=vomt-dpth/denhe
      qpdf(13)=vomr-dpra/denhe
c
c   pdf of helicity and components of Q X omega
c
      ntot(npq)=ntot(npq)+1
         do l=1,7
      if(abs(qpdf(l)).le.1.) then
      all=lipdh*(1.+qpdf(l))+1.
      ll=all
      nllc(npq,l,ll)=nllc(npq,l,ll)+1
       npou(npq,l)=npou(npq,l)+1
                             else
       nou(npq,l)=nou(npq,l)+1
c      write(91,*)npq,ic,kc,l,qpdf(l),ll,nou(npq,l)
                             endif
         enddo
      if(qpdf(8).le..1e-04) qpdf(8)=.1e-04
      all=lipdn*(1.+(alog(qpdf(8))/alog(10.)-3.)/5.)+1.
      ll=all
      nllc(npq,8,ll)=nllc(npq,8,ll)+1
         do l=9,mqu
      if(abs(qpdf(l)).le.1.) then
      all=lipdh*(1.+qpdf(l))+1.
      ll=all
      nllc(npq,l,ll)=nllc(npq,l,ll)+1
       npou(npq,l)=npou(npq,l)+1
                             else
       nou(npq,l)=nou(npq,l)+1
                             endif
         enddo
      if(mod(ntime,1).eq.0) then
c
c   joint pdf of helicity density  and dissipation
c     abscissa helicity density     ordinate dissipation  
c   
      asc=hede
      ord=enslo
      ala=liapdh*(asma1+asc)/asma1+1.
      lla=ala
      if(ord.le..1e-01) ord=.1e-01
      alo=liopdn*(1.+(alog(ord)/alog(10.)-3.)/5.)+1.
c     alo=liopdn*ord+1.
      llo=alo
      ncpdf(npq,5,lla,llo)=ncpdf(npq,5,lla,llo)+1
c
c   joint pdf of Q X Omeg density  and dissipation
c     abscissa helicity density     ordinate dissipation  
c   
      asc=vomx+vomt+vomr
      ala=liapdh*(asma1+asc)/asma1+1.
      lla=ala
      ncpdf(npq,9,lla,llo)=ncpdf(npq,9,lla,llo)+1
c
c   joint pdf of helicity density  and vomx
c     abscissa helicity density     ordinate vomx  
c   
      asc=hede
      ord=vomx
      ala=liapdh*(asma1+asc)/asma1+1.
      lla=ala
      alo=liopdh*(orma1+ord)/orma1+1.
      llo=alo
      ncpdf(npq,6,lla,llo)=ncpdf(npq,6,lla,llo)+1
c
c   joint pdf of helicity density  and vomt
c     abscissa helicity density     ordinate vomt  
c   
      asc=hede
      ord=vomt
      ala=liapdh*(asma1+asc)/asma1+1.
      lla=ala
      alo=liopdh*(orma1+ord)/orma1+1.
      llo=alo
      ncpdf(npq,7,lla,llo)=ncpdf(npq,7,lla,llo)+1
c
c   joint pdf of helicity density  and vomr
c     abscissa helicity density     ordinate vomr  
c   
      asc=hede
      ord=vomr
      ala=liapdh*(asma1+asc)/asma1+1.
      lla=ala
      alo=liopdh*(orma1+ord)/orma1+1.
      llo=alo
      ncpdf(npq,8,lla,llo)=ncpdf(npq,8,lla,llo)+1
                          endif
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
      enddo
      if(mod(ntime,1).eq.0) then
c
c   CALCULATION OF JOINT PDF at
c        certain radial positions
c    of fluctuting pressure and axial vorticity
c
      npq=1
      do jc=n2m,1,-1
      if(jc.eq.jprq(npq)) then
      do kc=1,n3m
      do ic=1,n1m                                                     
c
c   fluctuating vorticities at thge cell centre
c
      q3p=rhs(ic,jc,kc)-vorv(3,jc)
      prver=pr(ic,jc,kc)+0.5*(
     1            -(ros*0.5*rm(jc))**2)
      pcp=prver-pmed(jc)
c
c   joint pdf of p'  and omega'z
c     abscissa omega'z/rms omega'z  ordinate p'/rms p'
c   
      asc=q3p/vorstr(3,jc)
      ord=pcp/sqrt(pvc(4,jc))
      if(abs(asc).gt.asma) asc=asma
      if(abs(ord).gt.orma) ord=orma
      ala=liapdh*(asma+asc)/asma+1.
      lla=ala
      alo=liopdh*(orma+ord)/orma+1.
      llo=alo
      ncpdf(npq,4,lla,llo)=ncpdf(npq,4,lla,llo)+1
      enddo
      enddo
            npq=npq+1
                           endif
      enddo
      write(99,*)' npdf '
      do ny=1,nyf,10
      do npp=1,8
      npdf(npp)=0
      do ll=1,lipdf
      npdf(npp)=npdf(npp)+nllc(ny,npp,ll)
      enddo
      enddo
      write(99,134)ny,(npdf(nl),nl=1,8)
      enddo
      write(99,*)' nppcf '
      do ny=1,nyf,10
      do npp=1,9
      nppcf(npp)=0
      do la=1,liapdf
      do lo=1,liopdf
      nppcf(npp)=nppcf(npp)+ncpdf(ny,npp,la,lo)
      enddo
      enddo
      enddo
      write(99,134)ny,(nppcf(nl),nl=1,9)
  134 format('ny=',i2,3x,3x,9i8)
      enddo
                endif
      return                                                            
      end                                                               
c
c  ******************** subrout taver *********************
c   summation of the statistics to perform time average
c
      subroutine taver
c
      include 'param.f'
c
      disstpo=disstps+disstpo
      enpsvo=enpsv+enpsvo
      dissto=disst+dissto
      enavo=enav+enavo
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
c
c  ****************************** subrout pdfini **********************
c  defines intervals for the pdf
c
      subroutine pdfini   
      include 'param.f'
c************************************************************************
301   format(a4)
      read(15,301) dummy
      read(15,*)lipdf,liapdf,liopdf
      read(15,301) dummy
      read(15,*)asma,orma 
   83 format(i3.3)
      lipdh=(lipdf-1)/2
      lipdn=lipdf-1
      do ll=1,lipdf
c     aleh(ll)=(ll-1)/float(lipdn)
      aleh(ll)=10**(3.-5.*(1.-(ll-1)/float(lipdn)))
      enddo
      do ll=1,lipdf
      alev(ll)=(ll-1-lipdh)/float(lipdh)
      enddo
c     write(6,232)alev(1),alev(lipdf),lipdf,lipdn,lipdh
  232 format('  min max alev',2e12.4,3x,3i5)
      npq=1
      do jc=n2m,1,-1
      if(jc.eq.jprq(npq)) then
      nyf=npq
      ntot(npq)=0
            do l=1,mqu
      nou(npq,l)=0
      npou(npq,l)=0
            enddo
      ny=npq
         do ll=1,lipdf
            do np=1,mqu
      nllc(ny,np,ll)=0
            enddo
         enddo
            do np=1,mqu
      ntoc(ny,np)=0
      ntop(ny,np)=0
      ntom(ny,np)=0
            enddo
            do np=1,mqu
      pdft(ny,np)=0.
            enddo
      npq=npq+1
                           endif

      enddo
      asma1=1.
      orma1=1.
      liopdh=(liopdf-1)/2
      liopdn=liopdf-1
      liapdh=(liapdf-1)/2
      liapdn=liapdf-1
      do ll=1,liapdf
      alas1(ll)=(ll-1-liapdh)/float(liapdh)*asma1
      enddo
      do ll=1,liopdf
      alor0(ll)=(3-5.*(1.-(ll-1)/float(liopdn)))
c     alor0(ll)=(ll-1)/float(liopdn)*orma1
      enddo
      do ll=1,liopdf
      alor1(ll)=(ll-1-liopdh)/float(liopdh)*orma1
      enddo
      do ll=1,liapdf
      alas(ll)=(ll-1-liapdh)/float(liapdh)*asma
      enddo
      do ll=1,liopdf
      alor(ll)=(ll-1-liopdh)/float(liopdh)*orma
      enddo
c     write(6,332)alas(1),alas(liapdf),liapdf,liapdn,liapdh
  332 format('  min max alas',2e12.4,3x,3i5)
c     write(6,432)alor(1),alor(liopdf),liopdf,liopdn,liopdh
  432 format('  min max alas',2e12.4,3x,3i5)
c
c   here the file to visualize the joint pdf by TURB3D are
c   written
c
      filval='covcrpdf.dat'
      open(26,file=filval,form='formatted')
      write(26,*)liapdn,liopdn
      write(26,*) ((alas(la),la=1,liapdn),lo=1,liopdn)
     1         ,((alor(lo),la=1,liapdn),lo=1,liopdn)
      close(26)
      filval='cohelpdf.dat'
      open(26,file=filval,form='formatted')
      write(26,*)liapdn,liopdn
      write(26,*) ((alas1(la),la=1,liapdn),lo=1,liopdn)
     1         ,((alor1(lo),la=1,liapdn),lo=1,liopdn)
      close(26)
      filval='codispdf.dat'
      open(26,file=filval,form='formatted')
      write(26,*)liapdn,liopdn
      write(26,*) ((alas1(la),la=1,liapdn),lo=1,liopdn)
     1         ,((alor0(lo),la=1,liapdn),lo=1,liopdn)
      close(26)
      npq=1
      do jc=n2m,1,-1
      if(jc.eq.jprq(npq)) then
         do la=1,liapdf
         do lo=1,liopdf
            do np=1,mcqu
      ncpdf(npq,np,la,lo)=0
            enddo
         enddo
         enddo
            do np=1,mcqu
      ntcpdf(npq,np)=0
            enddo
            do np=1,mcqu
      cpdft(npq,np)=0.
            enddo
      npq=npq+1
                           endif

      enddo
      return
      end

c
c  ****************************** subrout pdffin **********************
c  here the pdf are calculated and are written
c  the joint pdf can be analized by the TURB3D package
c
      subroutine pdffin(time)
      include 'param.f'
      character*3 njpse
      itim=nint(time)
      write(pntim,67) itim
   67 format(i4.4)
      open(18,file='pdfglob.out')
      write(18,*) ' total pdf'
      npq=1
      do jc=n2m,1,-1
      if(jc.eq.jprq(npq)) then
                      do ll=1,lipdf
            do np=1,mqu
      ntoc(npq,np)=ntoc(npq,np)+nllc(npq,np,ll)
            enddo
                      enddo
                      do ll=1,lipdh
            do np=1,mqu
      ntom(npq,np)=ntom(npq,np)+nllc(npq,np,ll)
            enddo
                      enddo
                      do ll=lipdh+2,lipdf
            do np=1,mqu
      ntop(npq,np)=ntop(npq,np)+nllc(npq,np,ll)
            enddo
                      enddo
            do np=1,mqu
      ntom(npq,np)=ntom(npq,np)+nllc(npq,np,lipdh+1)*0.5
      ntop(npq,np)=ntop(npq,np)+nllc(npq,np,lipdh+1)*0.5
            enddo
      ntoz(npq)=0.
            do np=1,mqu
      ntoz(npq)=ntoz(npq)+ntoc(npq,np)
            enddo
      write(jsp,177) jc
  177 format(i2.2)
c
c    the pdf are written
c
      filval='disspdf'//jsp//'.'//pntim
      open(27,file=filval)
      filval='helccpdf'//jsp//'.'//pntim
      open(24,file=filval)
      filval='helicpdf'//jsp//'.'//pntim
      open(25,file=filval)
      filval='vcrompdf'//jsp//'.'//pntim
      open(26,file=filval)
      filval='prgrdpdf'//jsp//'.'//pntim
      open(28,file=filval)
      filval='trasfpdf'//jsp//'.'//pntim
      open(29,file=filval)
                 do ll=1,lipdn
      alevc=alev(ll)
            do np=1,mqu
      pdfvo(npq,np)=nllc(npq,np,ll)/float(ntoc(npq,np))
      pdft(npq,np)=pdft(npq,np)+pdfvo(npq,np)
            enddo
      write(27,178)aleh(ll),(pdfvo(npq,np),np=8,8)
      write(25,178)alevc,(pdfvo(npq,np),np=1,1)
      write(26,178)alevc,(pdfvo(npq,np),np=2,4)
      write(28,178)alevc,(pdfvo(npq,np),np=9,11)
      write(29,178)alevc,(pdfvo(npq,np),np=12,14)
      write(24,178)alevc,(pdfvo(npq,np),np=5,7)
  178 format(2x,e12.5,3x,6e14.5)
                  enddo
      write(18,198) (pdft(npq,np),np=1,mqu)
  198 format(2x,15e9.3)
      npq=npq+1
                           endif
            enddo
      close(24)
      close(25)
      close(26)
      close(28)
      close(29)
      close(27)
      npq=1
c
c   the number of events which give the pdf are written
c
      do jc=n2m,1,-1
      if(jc.eq.jprq(npq)) then
      write(jsp,177) jc
      filval='disspi'//jsp//'.'//pntim
      open(27,file=filval)
      filval='helccpi'//jsp//'.'//pntim
      open(24,file=filval)
      filval='helicpi'//jsp//'.'//pntim
      open(25,file=filval)
      filval='vcrompi'//jsp//'.'//pntim
      open(26,file=filval)
      filval='prgrdpi'//jsp//'.'//pntim
      open(28,file=filval)
      filval='trasdpi'//jsp//'.'//pntim
      open(29,file=filval)
                 do ll=1,lipdn
      alevc=alev(ll)
      write(27,178)aleh(ll),(float(nllc(npq,np,ll)),np=8,8)
      write(25,178)alevc,(float(nllc(npq,np,ll)),np=1,1)
      write(26,178)alevc,(float(nllc(npq,np,ll)),np=2,4)
      write(28,178)alevc,(float(nllc(npq,np,ll)),np=9,11)
      write(29,178)alevc,(float(nllc(npq,np,ll)),np=12,14)
      write(24,178)alevc,(float(nllc(npq,np,ll)),np=5,7)
                  enddo
      write(18,199) (ntoc(npq,np),np=1,mqu)
c     write(6,199) (ntoc(npq,np),np=1,mqu)
  199 format(2x,15i9)
      npq=npq+1
                           endif
            enddo
      close(24)
      close(25)
      close(26)
      close(28)
      close(29)
      close(27)
      do npq=1,nyf
      jcpr=jprq(npq)
            do np=1,mqu
      pdfpap(npq,np)=ntop(npq,np)/float(ntoc(npq,np))
      pdfpam(npq,np)=ntom(npq,np)/float(ntoc(npq,np))
            enddo
      write(18,194) jcpr,(pdfpam(npq,np),np=1,mqu)
      write(18,196) jcpr,(pdfpap(npq,np),np=1,mqu)
  196 format(1x,i3,1x,'+',15e9.3)
  194 format(1x,i3,1x,'-',15e9.3)
      enddo
c
c   Calculation of the joint pdf
c
c     write(6,*) ' total condit. pdf'
      write(18,*) ' total condit. pdf'
      npq=1
      do jc=n2m,1,-1
      if(jc.eq.jprq(npq)) then
                      do la=1,liapdf
                 do lo=1,liopdf
            do np=1,mcqu
      ntcpdf(npq,np)=ntcpdf(npq,np)+ncpdf(npq,np,la,lo)
            enddo
                 enddo
                      enddo
c 
c   here the joint pdf are  written to be 
c   analized by the TURB3D package
c
      write(jsp,177) jc
      filval='varcpdf'//jsp//'.'//pntim
      open(25,file=filval,form='formatted')
                 do la=1,liapdn
              do lo=1,liopdn
            do np=1,mcqu
      cpdfv(npq,np,la,lo)=ncpdf(npq,np,la,lo)/float(ntcpdf(npq,np))
      cpdft(npq,np)=cpdft(npq,np)+cpdfv(npq,np,la,lo)
            enddo
              enddo
                  enddo
      write(25,*)liapdn,liopdn
      write(25,*)ren,ros,time,ros
      write(25,*) ((cpdfv(npq,1,la,lo),la=1,liapdn),lo=1,liopdn)
     1         ,((cpdfv(npq,2,la,lo),la=1,liapdn),lo=1,liopdn)
     1         ,((cpdfv(npq,3,la,lo),la=1,liapdn),lo=1,liopdn)
     1         ,((cpdfv(npq,4,la,lo),la=1,liapdn),lo=1,liopdn)
      close(25)
      filval='vchelcpdf'//jsp//'.'//pntim
      open(25,file=filval,form='formatted')
      write(25,*)liapdn,liopdn
      write(25,*)ren,ros,time,ros
      write(25,*) 
     1          ((cpdfv(npq,6,la,lo),la=1,liapdn),lo=1,liopdn)
     1         ,((cpdfv(npq,7,la,lo),la=1,liapdn),lo=1,liopdn)
     1         ,((cpdfv(npq,8,la,lo),la=1,liapdn),lo=1,liopdn)
      close(25)
      filval='vchediscpdf'//jsp//'.'//pntim
      open(25,file=filval,form='formatted')
      write(25,*)liapdn,liopdn
      write(25,*)ren,ros,time,ros
      write(25,*) ((cpdfv(npq,5,la,lo),la=1,liapdn),lo=1,liopdn)
     1         ,((cpdfv(npq,9,la,lo),la=1,liapdn),lo=1,liopdn)
      close(25)
      write(18,198) (cpdft(npq,np),np=1,mcqu)
c     write(6,198) (cpdft(npq,np),np=1,mcqu)
      npq=npq+1
                           endif
            enddo
      return
      end




