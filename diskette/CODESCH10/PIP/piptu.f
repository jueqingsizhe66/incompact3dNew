c                                                                       
c  ************************* subrout energy  **********************  
c                                                                       
c     this subroutine calculates the turbulent kinetic energy          
c     and the mean velocity and pressure radial profiles
c                                                                       
      subroutine energy                    
      include 'param.f'
      common/prcost/pcost
      common/enetot/enet
      avgn=1./(float(n1m*n3m))
      vl123=1./(rc(n2)*n1m*n2m*n3m)
      volto=1./(pi*alx3d)
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
     1  (    q1(ic,jc,kc)+q1(ipv(ic),jc,kc) )*0.5/rm(jc)
      qcap(ic,jc,kc)=(q3(ic,jc,kc)+q3(ic,jc,kpv(kc)) )*0.5 
      enddo
      enddo
      enddo
      enej=0.
      enet=0.
      enpsv=0.
c
c   mean velocity radial profiles
c
      do jc=1,n2m                                                     
      jm=jmv(jc)                                                        
      jp=jc+1
      vmed(1,jc)=0.
      vmed(2,jc)=0.
      vmed(3,jc)=0.
      pmed(jc)=0.
      do kc=1,n3m 
      do ic=1,n1m                                                     
      pmed(jc)=pmed(jc)+pr(ic,jc,kc)
      vmed(2,jc)=vmed(2,jc)+dph(ic,jc,kc)
      vmed(1,jc)=vmed(1,jc)+dq(ic,jc,kc)
      vmed(3,jc)=vmed(3,jc)+qcap(ic,jc,kc)
      enddo
      enddo
      vmed(1,jc)=vmed(1,jc)*avgn
      vmed(2,jc)=vmed(2,jc)*avgn
      vmed(3,jc)=vmed(3,jc)*avgn
      pmed(jc)=pmed(jc)*avgn
      enddo
      do jc=1,n2m                                                     
      do kc=1,n3m
      do ic=1,n1m                                                     
      q1p=dq(ic,jc,kc)-vmed(1,jc)
      q2p=dph(ic,jc,kc)-vmed(2,jc)
      q3p=qcap(ic,jc,kc)-vmed(3,jc)
      enejp=(q1p**2+q2p**2+q3p**2)*vl123*rm(jc)*g2rm(jc)
      enej=enej+enejp
      q1t=dq(ic,jc,kc)
      q2t=dph(ic,jc,kc)
      q3t=qcap(ic,jc,kc)
      enetp=(q1t**2+q2t**2+q3t**2)*vl123*rm(jc)*g2rm(jc)
      enet=enet+enetp
      enddo
      enddo
      enddo
      enav=enej
      return                                                            
      end                                                               
c                                                                       
c  ************************* subrout veltur  **********************  
c                                                                       
c     this subroutine calculates the turbulent stresses
c     skewness and flatness are also calculated
c                                                                       
      subroutine veltur                    
      include 'param.f'
      common/prcost/pcost
      common/enetot/enet
      avgn=1./(float(n1m*n3m))
      vl123=1./(rc(n2)*n1m*n2m*n3m)
      volto=1./(pi*alx3d)
      call totvel
      call energy
c
c   mean velocity radial profiles
c
      do jc=1,n2m                                                     
      jm=jmv(jc)                                                        
      jp=jc+1
      vmed(1,jc)=0.
      vmed(2,jc)=0.
      vmed(3,jc)=0.
      pmed(jc)=0.
      do kc=1,n3m 
      do ic=1,n1m                                                     
      pmed(jc)=pmed(jc)+pr(ic,jc,kc)
      vmed(2,jc)=vmed(2,jc)+dph(ic,jc,kc)
      vmed(1,jc)=vmed(1,jc)+dq(ic,jc,kc)
      vmed(3,jc)=vmed(3,jc)+qcap(ic,jc,kc)
      enddo
      enddo
      vmed(1,jc)=vmed(1,jc)*avgn
      vmed(2,jc)=vmed(2,jc)*avgn
      vmed(3,jc)=vmed(3,jc)*avgn
      pmed(jc)=pmed(jc)*avgn
      enddo
      do jc=1,n2m                                                     
      pstr(jc)=0.
      do l=1,4
      ske(l,jc)=0.
      fla(l,jc)=0.
      enddo
      do l=1,6
      turstr(l,jc)=0.
      enddo
      do kc=1,n3m
      do ic=1,n1m                                                     
      q1p=dq(ic,jc,kc)-vmed(1,jc)
      q2p=dph(ic,jc,kc)-vmed(2,jc)
      q3p=qcap(ic,jc,kc)-vmed(3,jc)
      pcp=pr(ic,jc,kc)-pmed(jc)
      pstr(jc)=pstr(jc)+pcp*pcp
      turstr(1,jc)=turstr(1,jc)+q1p**2
      turstr(2,jc)=turstr(2,jc)+q2p**2
      turstr(3,jc)=turstr(3,jc)+q3p**2
      turstr(4,jc)=turstr(4,jc)-q1p*q2p
      turstr(5,jc)=turstr(5,jc)-q3p*q2p
      turstr(6,jc)=turstr(6,jc)-q3p*q1p
      ske(1,jc)=ske(1,jc)+q1p**3
      ske(2,jc)=ske(2,jc)+q2p**3
      ske(3,jc)=ske(3,jc)+q3p**3
      ske(4,jc)=ske(4,jc)+pcp**3
      fla(1,jc)=fla(1,jc)+q1p**4
      fla(2,jc)=fla(2,jc)+q2p**4
      fla(3,jc)=fla(3,jc)+q3p**4
      fla(4,jc)=fla(4,jc)+pcp**4
      enddo
      enddo
      do l=1,3
      turstr(l,jc)=sqrt( turstr(l,jc)*avgn )
      enddo
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
      pstr(jc)= pstr(jc)*avgn 
      enddo
      return                                                            
      end                                                               
c
c  ******************** subrout taver *********************
c   this routine adds the statistics to perform the
c   averages in time
c
      subroutine taver
c
      include 'param.f'
c
      enavo=enav+enavo
      do  l=1,3
      do  j=1,n2m
      vmeo(l,j)=(vmeo(l,j)+vmed(l,j))
      enddo
      enddo
      do  l=1,4
      do  j=1,n2m
      skeo(l,j)=(skeo(l,j)+ske(l,j))
      flao(l,j)=(flao(l,j)+fla(l,j))
      enddo
      enddo
      do  l=1,6
      do  j=1,n2m
      tursto(l,j)=(tursto(l,j)+turstr(l,j))
      enddo
      enddo
      do  j=1,n2m
      psto(j)=(psto(j)+pstr(j))
      pmeo(j)=(pmeo(j)+pmed(j))
      enddo
      return
      end
c
c  ******************** subrout tavwri *********************
c
c   write quantities to evaluate averages in time
c   when the calculation restarts
c
      subroutine tavwri(nav,time)
c
      include 'param.f'
      character*10 tawrfi
      tawrfi='tavwri.out'
      itime=nint(time)
      write(ipfi,82)itime
   82 format(i4.4)
      namfi3='tavwri'//ipfi//'.dat'
      open(70,file=tawrfi,form='unformatted')
      write(6,*) 'the continuation file',namfi3,' was written at ',
     1          'nav=',nav
      write(70)enavo,cfo,nav,dp3mo
      do j=1,n2m
      write(70)
     1(vmeo(l,j),l=1,3),
     1(tursto(l,j),l=1,6),(skeo(l,j),l=1,4),(flao(l,j),l=1,4)
     1,pmeo(j),psto(j)
      enddo
      close(70)
      return
      end
c
c  ******************** subrout tavrea *********************
c
c   read quantities to evaluate averages in time
c   when the calculation restarts
c
      subroutine tavrea(nav,nap)
c
      include 'param.f'
      dimension fo(2*m2),yn(2*m2),fni(2*m2),yu(2*m2)
      dimension y2sol(2*m2)
      common/dp3t/dp3th
      character*10 tawrfi
      n2lm=n2l-1
      tawrfi='tavwri.out'
      open(70,file=tawrfi,form='unformatted')
      write(6,*)'  read past averages from ',tawrfi
      tawrfi='vel.read'
      open(61,file=tawrfi)
      tawrfi='skflv.read'
      open(62,file=tawrfi)
      read(70)enavo,cfo,nav,dp3mo
      write(61,133)enavo,cfo,nav,dp3mo
  133 format(2x,2e11.4,2x,i4,2x,e12.5)
      dp3th=dp3mo*nap/nav
      do j=1,n2lm
      read(70)
     1(vmeo(l,j),l=1,3),
     1(tursto(l,j),l=1,6),(skeo(l,j),l=1,4),(flao(l,j),l=1,4)
     1,pmeo(j),psto(j)
      write(61,134) 
     1(vmeo(l,j),l=1,3),
     1pmeo(j),psto(j),
     1(tursto(l,j),l=1,6)
      write(62,134) 
     1(skeo(l,j),l=1,4),(flao(l,j),l=1,4)
      enddo
  134 format(1x,11e11.4)
      close(61)
      close(62)
  202 format(1x,12e10.3)
      close(70)
      if(n2lm.ne.n2m) then
      n2om=n2lm
c
c  mean velocities
c
      do j=1,n2om
      y2sol(j)=rmo(j)
      enddo
      do l=1,3
      do j=1,n2om
      fo(j)=vmeo(l,j)
      enddo
      pn2=0.
      dp1=0.
      call spline(y2sol,fo,n2om,dp1,dpn2,yu)
c   
      do j=1,n2m
      xx=rm(j)
      call splint(y2sol,fo,yu,n2om,xx,yy)
      fni(j)=yy
      enddo
      do j=1,n2m
      vmeo(l,j)=fni(j)
      enddo
      enddo
c
c  mean skewness
c
      do l=1,4
      do j=1,n2om
      fo(j)=skeo(l,j)
      enddo
      dpn2=0.
      dp1=0.
      call spline(y2sol,fo,n2om,dp1,dpn2,yu)
c
      do j=1,n2m
      xx=rm(j)
      call splint(y2sol,fo,yu,n2om,xx,yy)
      fni(j)=yy
      enddo
      do j=1,n2m
      skeo(l,j)=fni(j)
      enddo
      enddo
c
c  mean flattness
c
      do l=1,4
      do j=1,n2om
      fo(j)=flao(l,j)
      enddo
      dpn2=0.
      dp1=0.
      call spline(y2sol,fo,n2om,dp1,dpn2,yu)
c
      do j=1,n2m
      xx=rm(j)
      call splint(y2sol,fo,yu,n2om,xx,yy)
      fni(j)=yy
      enddo
      do j=1,n2m
      flao(l,j)=fni(j)
      enddo
      enddo
c
c  mean velocity rms
c
      do l=1,6
      do j=1,n2om
      fo(j)=tursto(l,j)
      enddo
      dpn2=0.
      dp1=0.
      call spline(y2sol,fo,n2om,dp1,dpn2,yu)
c
      do j=1,n2m
      xx=rm(j)
      call splint(y2sol,fo,yu,n2om,xx,yy)
      fni(j)=yy
      enddo
      do j=1,n2m
      tursto(l,j)=fni(j)
      enddo
      enddo
c
c  mean pressure
c
      do j=1,n2om
      fo(j)=pmeo(j)
      enddo
      pn2=0.
      dp1=0.
      call spline(y2sol,fo,n2om,dp1,dpn2,yu)
c
      do j=1,n2m
      xx=rm(j)
      call splint(y2sol,fo,yu,n2om,xx,yy)
      fni(j)=yy
      enddo
      do j=1,n2m
      pmeo(j)=fni(j)
      enddo
                  endif
      return
      end
