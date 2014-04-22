c***********************************************************************
c                                                                       *
c  ****************************** subrout indic **********************  *
c                                                                       *
c     in this subroutine the indices ip,im,jp,jm,kp,km are calculated.  *
c                                                                       *
c************************************************************************
      subroutine indic                                                  
      include 'param.f'
c                                                                       
c                                                                       
c   azimuthal periodic direction                                               
c                                                                       
      do 1 ic=1,n1m                                                     
      imv(ic)=ic-1                                                      
      if(ic.eq.1) imv(ic)=n1m                                           
      ipv(ic)=ic+1                                                      
      if(ic.eq.n1m) ipv(ic)=1                                           
    1 continue                                                          
c                                                                       
c   vertical periodic direction                                               
c                                                                       
      do 4 kc=1,n3m
      kmv(kc)=kc-1                                                      
      kpv(kc)=kc+1                                                     
      if(kc.eq.1) kmv(kc)=n3m                                           
      if(kc.eq.n3m) kpv(kc)=1                                           
    4 continue                                                          
                                                                        
c     direction normal to the radial boundary 
                                                                        
      do 3 jc=1,n2m                                                     
      jmv(jc)=jc-1                                                      
      jpv(jc)=jc+1                                                      
      if(jc.eq.1) jmv(jc)=jc                                            
      if(jc.eq.n2m) jpv(jc)=jc                                          
    3 continue                                                          
c                                                                       
c   indices for the axis of symmetry and the external wall                                                   
c                                                                       
      do 15 jc=1,n2m                                                    
      jpc(jc)=jpv(jc)-jc                                                
      jmc(jc)=jc-jmv(jc)                                                
      jup(jc)=1-jpc(jc)                                                 
      jum(jc)=1-jmc(jc)                                                 
   15 continue                                                          
      do i=1,n1m
      isym(i) = i + n1m/2
      if(isym(i).gt.n1m) isym(i) = isym(i) - n1m
      enddo
      return                                                            
      end                                                               
c                                                                       *
c  ****************************** subrout divgck  ***********************
c                                                                       *
c     this subroutine checks divg(q) to see the local residue.          *
c                                                                       *
c************************************************************************
      subroutine divgck(qmax,qtot)                             
      include 'param.f'
c                                                                       
c     This is just a check on the divergence of the free-divergent
c     velocity field. The calculation is stopped if QMAX > RESID
c                                                                       
      qtot=0.                                                           
      qmax=0.                                                           
      uvol=1./(dx1*dx2*dx3)
      do 11 kc=1,n3m                                                    
      kp=kpv(kc)                                                        
      do 11 jc=1,n2m                                                    
      jp=jc+1                                                           
      usrnu1=dx1/rm(jc)**2
      usrnu2=dx2/g2rm(jc)/rm(jc)
      do 11 ic=1,n1m                                                   
      ip=ipv(ic)                                                        
      dqca1= (q1(ip,jc,kc)-q1(ic,jc,kc))*usrnu1
      dqca2= (q2(ic,jp,kc)-q2(ic,jc,kc))*usrnu2
      dqca3= (q3(ic,jc,kp)-q3(ic,jc,kc))*dx3
      dqcap= dqca1+dqca2+dqca3
      qtot=qtot+dqcap*rm(jc)*g2rm(jc)/dx2/dx3*uvol
      qmax=amax1(abs(dqcap),qmax)                                       
   11 continue
c     write(62,*)' DIVGCK',imxq,jmxq,kmxq,qtot,qmax
      return                                                            
      end                                                               
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
      if(irid.eq.1) then
      dx1=2.*pi/(float(n1m)*float(lamb))     
       else
      dx1=2.*pi/float(n1m)
      endif
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
      common/npjet/n2t
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

c$$$$$$$$$$$$$Computation of geometry terms for non-uniform grids$$$$$$$$

      g2rm(j)=(rc(j+1)-rc(j))*dx2

c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

   12 continue                                                          

      do 122 j=2,n2m
c$$$$$$$$$$$$$Computation of geometry terms for non-uniform grids$$$$$$$$

      g2rc(j)=(rc(j+1)-rc(j-1))*dx2*0.5
c$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
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
      subroutine openfi                                                 *
c                                                                       *
c************************************************************************
      include 'param.f'
      open(46,file='nfpipe')
      read(46,'(a)')filcnw
      read(46,'(a)')filcnr
      read(46,'(a)')filth
      read(46,'(a)')filvm
      read(46,'(a)')filpo
      read(46,'(a)')filen
      read(46,'(a)')filet
      read(46,'(a)')filer
      read(46,'(a)')filez
      read(46,'(a)')filed
      read(46,'(a)')filev
      open(32,file=filth)
      open(33,file=filvm)
      open(34,file=filpo)
      open(39,file=filen)
      open(40,file=filet)
      open(41,file=filer)
      open(42,file=filez)
      open(49,file=filed)
      open(50,file='piqm.out')
      open(59,file='piav.out')
      rewind 12
      rewind 33
      rewind 34
      rewind 32
      rewind 39
      rewind 40
      rewind 41
      rewind 42
      rewind 49
      rewind 48
      return
      end   
c************************************************************************
c                                                                       *
c  ****************************** subrout inirea ********************** *
c     read input flow fields                                            *
c                                                                       *
c************************************************************************
      subroutine inirea(ntil,time,ntt,ncount,nap)
      include 'param.f'
      common/avgin/avpscn(m2),vnew(m2)
      character*20 namdir
      character*67 namfil
      itime=nint(time)
      write(ipfi,82)itime
   82 format(i4.4)
      namdir='field'//ipfi//'.dat'
      namfil='../fields/field'//ipfi//'.dat'
      write(6,*)'     legge da ',namdir
      open(13,file=namfil,form='unformatted')
      nfil=13                                                           
      read(nfil) n1l,n2l,n3l                                        
      read(nfil) epsil,lamb,re,timl 
      write(6,*)' legge da ',namfil,'     ntime=',ntil,
     1          '  n1l,n2l,n3l ',n1l,n2l,n3l
c
c    here the large fields with the three velocity components
c    and pressure are read
c    this part should be changed if the restarting file 
c    does not contain the q2 velocity component
c
      read(nfil)  (((q1(i,j,k),i=1,n1l),j=1,n2l),k=1,n3l),           
     1            (((q2(i,j,k),i=1,n1l),j=1,n2l),k=1,n3l),        
     1            (((q3(i,j,k),i=1,n1l),j=1,n2l),k=1,n3l),
     1            (((pr(i,j,k),i=1,n1l),j=1,n2l),k=1,n3l)
      read(nfil) ntii,ntt,npp
      if(ncount.eq.0) then
      enavo=0.
      disstpo=0.
      enpsvo=0.
      dissto=0.
      cfo=0.
      do  l=1,3
      do  j=1,n2m
      vmeo(l,j)=0.
      voro(l,j)=0.
      vdtomo(l,j)=0.
      enddo
      enddo
      do  l=1,4
      do  j=1,n2m
      pvmo(l,j)=0.
      skeo(l,j)=0.
      flao(l,j)=0.
      enddo
      enddo
      do  l=1,6
      do  j=1,n2m
      tursto(l,j)=0.
      vorsto(l,j)=0.
      vcromo(l,j)=0.
      enddo
      enddo
      do  j=1,n2m
      pmeo(j)=0.
      dissjo(j)=0.
      enddo
               endif
      return                                                            
      end                                                               
c************************************************************************
c                                                                       *
c  ****************************** subrout initia ********************** *
c     initial zero conditions in the whole field.                       *
c                                                                       *
c************************************************************************
      subroutine initia
      include 'param.f'
      do 4 j=1,n2                                                       
      do 4 i=1,n1                                                     
      do 4 k=1,n3                                                       
      pr(i,j,k)=0.                                                      
      q1(i,j,k)=0.                                                      
      q2(i,j,k)=0.                                                      
      q3(i,j,k)=0.                                                      
    4 continue                                                          
      do  l=1,5
      do  j=1,n2m
      vtvo(l,j)=0.
      dvto(l,j)=0.
      enddo
      enddo
      return                                                            
      end                                                               
