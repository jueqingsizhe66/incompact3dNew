c**********************************************************************ni*
c                                                                      *
c     this routine gives the initial condition for the 
c     simulation of a 3D tripole stability simulation
c                                                                      * 
c***********************************************************************
      subroutine intrip
      include 'param.f'
      common/vorini/vorin(m1,m2)
      write(6,*)  '  '
      write(32,*) '  '
      write(6,*) ' I N I T I A L   C O N D I T I O N '
      write(32,*) ' I N I T I A L   C O N D I T I O N '
      write(6,*)  '  '
      write(32,*) '  '
c
c    velocity and vorticity distribution
c
      call inivtr
c
c*******************      passive scalar      ************************
c
      pscmax=0.
      do  j=1,n2m
      if(alpha.gt.0.) then
      ram=rc(j)/sqrt(crad)
      xex=ram**alpha*alpha
                      else
      ram=rc(j)
      xex=ram
                      endif
      if(xex.lt.rap) then
      pscas=1.
                     endif
      if(xex.ge.rap.and.xex.le.rapsc) then
      yd=(xex-rap)/(rapsc-rap)
      fy=3.*yd**2-2.*yd**3
      pscas=1.-fy
                              endif
      if(xex.gt.rapsc) then
      pscas=0.
                       endif
      do  k=1,n3m
      do  i=1,n1m
      psc(i,j,k)=pscas
      pscmax=amax1(pscmax,abs(psc(i,j,k)))                                 
      if(i.eq.1.and.k.eq.1) write(29,128) rm(j),psc(i,j,k)
  128 format(3x,5e11.4)
      enddo
      enddo
      enddo
      write(6,*)'   pscmax=',pscmax
      close(29)
      return                                                            
      end                                                               

c***********************************************************************
c                                                                      *
c    the axial velocity distribution is given it can have 
c    a Gaussian distribution if alpha < 0.
c    or the distribution of an isolated vortex that is function
c    of alpha. This permits to investigate the eventual formation
c    of  a tripolar structure
c                                                                      *
c***********************************************************************
      subroutine inivtr
      include 'param.f'
      dimension q1j(m2)
      common/vorini/vorin(m1,m2)
c                                                                       
      pi=2.*asin(1.)                                                    
c     open(29,file='pert.junk')
c     open(28,file='vpro.junk')
      ampl=amp                                                          
      m=n1m                                                              
      do k=1,n3m                                                    
c                                                                       
c     q1 comp.   
c
      if(alpha.gt.0.0) write(6,*)' isolated vortex with al=',alpha
      if(alpha.eq.0.0) stop
      q1max=0.
      q1pt=0.
      iseed=-1
      do i=1,n1m
      if(alpha.gt.0.) then
      rar=rm(1)/sqrt(crad)
      xexp=rar**alpha
      q1j(1)=0.5*rm(1)*exp(-xexp)
      q1(i,1,k)=q1j(1)*rm(1)
                      else
      rar=rm(1)/(crad*sqrt(2.))
      xexp=rar**2
      q1j(1)=1./(pi*rm(1))*(1.-exp(-xexp))
      q1(i,1,k)=q1j(1)*rm(1)
                      endif
      if(i.eq.1.and.k.eq.1) write(27,128) rm(1),q1j(1)
      enddo
      enddo
      do j=2,n2m
      if(alpha.gt.0.) then
      ram=rc(j)/sqrt(crad)
      xem=ram**alpha*alpha
      xex=((xem-rap)/sigr)**2
      if(xex.gt.25.) then
      aexp=0.                                                           
                     else
      aexp=exp(-xex)
                     endif
                      else
      ram=rc(j)/(crad*sqrt(2.))
      xem=ram**2
      xex=((xem-rap)/sigr)**2
      if(xex.gt.25.) then
      aexp=0.                                                           
                     else
      aexp=exp(-xex)
                     endif
                     endif
c     write(29,129) j,rc(j),ram,xex,aexp
  129 format(3x,i4,3x,5e11.4)
c
c    a random or a sin perturbation are added
c    in the axial as in the azimuthal directions
c    Remember that the rand generator depends on 
c    the computer available
c
      do k=1,n3m                                                    
      zm=(zz(k)+zz(k+1))*0.5
      if(iran3.eq.0) then
      verper=1.+amp3*sin(2.*pi*zm/alx3*nwa3)
                    else
      verper=1.+amp3*(0.5-rand())*2.
                    endif
      if(n3m.eq.1) verper=1.
      do i=1,n1m
      if(iran1.eq.0) then
      ragp=rc(j)+amp*sin(thetac(i)*nwa)*verper*aexp
                    else
      ragp=rc(j)+amp*(0.5-rand())*2.*verper*aexp
                    endif
      if(alpha.gt.0.) then
      rar=ragp/sqrt(crad)
      xexp=rar**alpha
                     else
      rar=ragp/(crad*sqrt(2.))
      xexp=rar**2
                      endif
      if(xexp.gt.30.) then
      axexp=0.
                      else
      axexp=exp(-xexp)
                      endif
      if(alpha.gt.0.) then
      vorc=(1.-rar**alpha*0.5*alpha)*axexp
                     else
      vorc=1./(pi*crad**2)*axexp
                     endif
      if(vorc.gt.0.) vorin(i,j)=vorc
      drvel=rc(j)*vorc/dx2*g2rc(j)
      q1j(j)=(q1j(j-1)*rm(j-1)+drvel)/rm(j)   
      q1(i,j,k)=q1j(j)*rm(j)
      q1max=amax1(q1max,abs(q1j(j)))
      if(i.eq.1.and.k.eq.1) write(28,128) rc(j),vorc
      if(i.eq.1.and.k.eq.1) write(27,128) rm(j),q1j(j)
  128 format(3x,5e11.4)
      enddo
      enddo
      enddo
c                                                                       
c  for the no-slip bottom wall q1 goes to zero
c  with a parabolic profile
c
      blwa=alx3d*.1
      write(26,*)blwa
      if(inslws.eq.1) then
      do k=1,n3m                                                    
      zm=(zz(k)+zz(k+1))*0.5
      zds=(zz(k)+zz(k+1))*0.5-zz(1)
      if(zds.le.blwa) then
      diwas=zds/blwa
                     else
      diwas=1.
                     endif
      zdn=zz(n3)-(zz(k)+zz(k+1))*0.5
      if(zdn.le.blwa) then
      diwan=zdn/blwa
                     else
      diwan=1.
                     endif
      dampwa=diwas**2*diwan**2
      write(26,126)k,zm,zds,zdn,diwas,diwan,dampwa
  126 format(2x,i3,3x,7e12.6)
      do j=1,n2m
      do i=1,n1m
      q1(i,j,k)=q1(i,j,k)*dampwa
      enddo
      enddo
      enddo
                      endif
      close(26)
c                                                                       
c   perturbation on the axial velocity
c
      q3max=0.
      do kc=1,n3
      do jc=1,n2m                                                   
      jp=jc+1                                                           
      do ic=1,n1m
      ip=ipv(ic)
      q3(ic,jc,kc)=amp3*sin(2.*pi*nwa3*zz(kc)/alx3)    
      q3max=amax1(q3max,abs(q3(ic,jc,kc)))
      enddo
      enddo
      enddo
      iseed=-1
c
c    radial velocity from continuity
c     q2 comp.                                                          
c
      q2max=0.                                                          
      do kc=1,n3m
      usrnu3=dx3/g3rm(kc)
      kp=kc+1
      do ic=1,n1m
      ip=ipv(ic)  
      q2(ic,1,kc)=0.  
      do jc=1,n2m
      jp=jc+1
      usrnu1=dx1/rm(jc)
      usrnu2=dx2/g2rm(jc)/rm(jc)
      dqca1= (q1(ip,jc,kc)-q1(ic,jc,kc))*usrnu1
      dqca2= (q2(ic,jp,kc)-q2(ic,jc,kc))*usrnu2
      dqca3= (q3(ic,jc,kp)-q3(ic,jc,kc))*usrnu3
      q2(ic,jp,kc)=q2(ic,jc,kc)-(dqca3+dqca1)/usrnu2
      q2max=max(q2max,abs(q2(ic,jc,kc)))                                 
      enddo
      enddo
      enddo
      call divgck(qmax,qtot)                                   
      write(6,705) qmax,qtot
      write(32,705) qmax,qtot
  705 format(3x,'from inqpr  qmax and qtot  =',2e11.4)
      write(6,700) iran1,iran3,q1max,q2max,q3max
      write(32,700) iran1,iran3,q1max,q2max,q3max
  700 format(1x,'iran1=',i3,2x,' iran3=',i3,2x,
     1      'vmx1=',e11.4,2x,'vmx2=',e11.4,2x,'vmx3=',e11.4)
      close(29)
      close(28)
      close(27)
c
c   q1 is divided by rm since in this code 
c                  q1=v_theta    
c
      do i=1,n1m
      do j=1,n2m
      do k=1,n3m                                                    
      q1(i,j,k)=q1(i,j,k)/rm(j)
      enddo
      enddo
      enddo
      return                                                            
      end                                                               
