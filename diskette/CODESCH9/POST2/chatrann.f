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
      usrnu1=dx1
      do 11 kc=1,n3m                                                    
      kp=kpv(kc)                                                        
      do 11 jc=1,n2m                                                    
      jp=jc+1                                                           
      usrnu2=dx2/caj(jc)
      do 11 ic=1,n1m                                                   
      ip=ipv(ic)                                                        
      dqca1= (q1(ip,jc,kc)-q1(ic,jc,kc))*usrnu1
      dqca2= (q2(ic,jp,kc)-q2(ic,jc,kc))*usrnu2
      dqca3= (q3(ic,jc,kp)-q3(ic,jc,kc))*dx3
      dqcap= dqca1+dqca2+dqca3
      qtot=qtot+dqcap*caj(jc)/dx2/dx3*uvol
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
      dx1=alx1/float(n1m)
      dx2=1./float(n2m)                                                 
      dx3=alx3/float(n3m)                                               
      write(6,100)dx1,dx2,dx3                           
      write(32,100)dx1,dx2,dx3                           
 100  format(3x,'mesh size: dx1= ',e9.4,' dx2= ',e9.4,' dx3= ',e9.4)
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
      dimension eta(m2)
      common/strpar/str2
c

      tstr2=tanh(str2*0.5)
      do 63 j=1,n2
      x2=(j-1)/float(n2m)
      eta(j)=0.5*(1.+tanh(str2*(x2-0.5))/tstr2)
   63 continue
      do 65 j=1,n2
      x2=(j-1)/float(n2m)
      yp2(j)=(-0.5+eta(j))*2.
      y(j)=yp2(j)
      write(68,*) x2,y(j)
   65 continue
c     print *,'y al centro',y(n2/2+1)
      do 67 j=1,n2m
      y2s(j)=(yp2(j)+yp2(j+1))*0.5
   67 continue
      if(mod(n2,2).eq.1) then
      print *,'n2 dispari forza la simmetria nella coordinata y' 
      y(n2/2+1)=0.
      do j=1,n2/2
      y(n2-j+1)=-y(j)
      end do
      do j=1,n2/2
      y2s(n2-j)=-y2s(j)
      end do
      end if
      open(47,file='ymesh')
      do j=1,n2m
      dy=y(j+1)-y(j) 
      write(47,*) y2s(j),dy
      end do
      close(47)
      do j=1,n2m/2+1
      ydp=(y(j)-y(1) )*180
      if(ydp.gt.20) then
      jpri=j
      go to 23
                    endif
      end do
   23 continue
      kpri(1)=n3m/4+1
      ipri(1)=n1m/4+1
      kpri(2)=n3m/4+1
      ipri(2)=3*n1m/4+1
      kpri(3)=n3m/2+1
      ipri(3)=n1m/2+1
      kpri(4)=3*n3m/4+1
      ipri(4)=n1m/4+1
      kpri(5)=3*n3m/4+1
      ipri(5)=3*n1m/4+1
      jin=1
      jfi=n2m/2
      do l=1,5
      iin=ipri(l)-3
      ifi=ipri(l)+3
      kin=kpri(l)-3
      kfi=kpri(l)+3
      npoin=0
      do i =iin,ifi
      do j =jin,jfi
      do k =kin,kfi
      npoin=npoin+1
      enddo
      enddo
      enddo
      enddo
      write(6,*)' print npoin=',npoin
c
      do  i=1,n1
          yp1(i) = float(i-1)/float(n1m) * alx1
      enddo
      do  k=1,n3
          yp3(k) = float(k-1)/float(n3m) * alx3
      enddo
      do  i=1,n1m
          y1s(i) = 0.5*(yp1(i) + yp1(i+1) )
      enddo
      do  k=1,n3m
          y3s(k) = 0.5*(yp3(k) + yp3(k+1) )
      enddo
c
c  *********                 i,j+1/2
c
      do 1 j=1,n2m
c
c  ********* derivatives of cartesian corrdinates interior points
c  cn2 deriv respect to x2
c
      cn22=(y(j+1)-y(j))*dx2
   10 continue
      caj(j)=cn22
    1 continue
c
c  *********                 i,j
c
      do 6 j=1,n2
      if(j.gt.1.and.j.lt.n2) go to 31
c
c  ********* derivatives of cartesian corrdinates horizontal boundary
c            three points backw.
      if(j.eq.1) js=1
      if(j.eq.n2) js=-1
      cn22=-js*(y(j)-y(j+js))*dx2
      go to 34
   31 continue
c
c  ********* derivatives of cartesian corrdinates interior points
c
      cn22=(y(j+1)-y(j-1))*dx2*0.5
   60 continue
   34 continue
      cac(j)=cn22
    6 continue
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
      open(62,file='enspeck.out')
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
      open(74,file='flatcheck.out')
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
      common/ichtim/itimcf,itimcf1
      common/qwallo/q1s(m1,m3),q2s(m1,m3),q3s(m1,m3)
      common/qwalup/q1n(m1,m3),q2n(m1,m3),q3n(m1,m3)
      common/walfr/cfuo,cflo
      character*80 namfil
      itime=nint(time)
      write(ipfi,82)itime
   82 format(i4.4)
c     write(6,*)'     in inirea da ',lipdf
c
c   this file was written during the simulation
c   by the code in the directory CODESCH9/CHA
c   This file is a restarting file
c   of large dimensions since the pressure and
c   q2 are given. These are not necessary
c   since q2 can be obtained by the divg if
c   only q1 and q3 re read the pressure can be
c   obtained by performing one time step advancement
c   by the code in CODESCH9/CHA
c
      namfil='../field'//ipfi//'.dat'
      write(6,*)'in inirea field ntil=',ntil,namfil
      open(13,file=namfil,form='unformatted')
      nfil=13                                                           
      read(nfil)n1lm,n2l,n3lm
c     write(6,*)'from ',namfil,'n1lm,n2l,n3lm',n1lm,n2l,n3lm
      read(nfil)ntiml,timl,ene0,dp3nsl,enav,cfn,dtl
      read(nfil) (((q1(i,j,k),i=1,n1lm),j=1,n2l),k=1,n3lm),
     1            (((q2(i,j,k),i=1,n1lm),j=1,n2l),k=1,n3lm),
     1            (((q3(i,j,k),i=1,n1lm),j=1,n2l),k=1,n3lm),
     1            (((pr(i,j,k),i=1,n1lm),j=1,n2l),k=1,n3lm)
c     write(6,*)' read q_i and p'
      read(nfil)ntiml,timl,vit(1),vit(2),vit(3)
     1         ,dp3nl,cfluw,enel,vmax(2),vmax(3),cfllw
      read(nfil) ((q1s(i,k),i=1,n1lm),k=1,n3lm),
     1            ((q2s(i,k),i=1,n1lm),k=1,n3lm),
     1            ((q3s(i,k),i=1,n1lm),k=1,n3lm),
     1            ((q1n(i,k),i=1,n1lm),k=1,n3lm),
     1            ((q2n(i,k),i=1,n1lm),k=1,n3lm),
     1            ((q3n(i,k),i=1,n1lm),k=1,n3lm)
      close(nfil)

  161 format(5e12.4)
      if(ncount.eq.0) then
      enavo=0.
      disstpo=0.
      enpsvo=0.
      dissto=0.
      cfo=0.
      cfuo=0.
      cflo=0.
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
          do l=1,12
      do  j=1,n2m
      sququo(l,j)=0.
      flaquo(l,j)=0.
      skequo(l,j)=0.
      enddo
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
