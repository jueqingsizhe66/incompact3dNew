c************************************************************************
c                                                                       *
c     ********* subrout pricor ******************                       *
c    write the x-y coordinates to visualize the
c    time evolution of a dipole moving in the r-theta plane
c                                                                       *
c************************************************************************
      subroutine pricor
      include 'param.f'
      common/cooxy/xxm(m1,m2),yym(m1,m2)
      do j=1,n2
      do i=1,n1m
      thetal=thetac(i)
      xxm(i,j)=rc(j)*cos(thetal)
      yym(i,j)=rc(j)*sin(thetal)
      enddo
      enddo
      do i=1,n1m
      xxm(i,1)=0.
      yym(i,1)=0.
      enddo
      do j=1,n2
      xxm(n1,j)=xxm(1,j)
      yym(n1,j)=yym(1,j)
      enddo
      namfile='cordxy.dat'
      open(18,file=namfile,form='unformatted')
      write(6,*) 'in prico xy' ,n1,n2m
      write(18) n1,n2
      write(18)
     1   ((xxm(i,j),i=1,n1),j=1,n2),
     1   ((yym(i,j),i=1,n1),j=1,n2)
      close(18)
      return
      end
c************************************************************************
c                                                                       *
c     ********* subrout outhdi ******************                       *
c                                                                       *
c************************************************************************
      subroutine outhdi
      include 'param.f'
      common/cooxy/xxm(m1,m2),yym(m1,m2)
      common/timw/timew
c
c  ***********  compute the axial vorticity component
c               at         i,j,k+1/2
c
      kc=1
      do ic=1,n1m
      im=imv(ic)
      do jc=2,n2m
      jm=jc-1
      dq1x2=(q1(ic,jc,kc)-q1(ic,jm,kc))*dx2/g2rc(jc)
      dq2x1=(q2(ic,jc,kc)-q2(im,jc,kc))*dx1/rc(jc)
      vorz=(dq1x2-dq2x1)/rc(jc)
      vor(ic,jc)=vorz
      enddo
      jc=n2
      jm=n2m
      dq1x2=-1./rc(jc)*q1(ic,jm,kc)*dx2*2./g2rc(jc)
      vor(ic,jc)=dq1x2*(1-islip)
      enddo
      jc=1
      vozcm=0.
      do ic=1,n1m
      vozcm=vozcm+vor(ic,2)
      enddo
      do ic=1,n1m
      vor(ic,jc)=vozcm/n1m
      enddo
      ampsi = -1.
c
c    evaluates the streamfunction 
c
      do ic=1,n1m
         psi(ic,1) = 0.
         do jc=2,n2
            jm = jc -1
            psi(ic,jc) = psi(ic,jm) + q1(ic,jm,1)*g2rm(jm)/dx2/rm(jm)
            ampsi = max(ampsi,psi(ic,jc))
         end do
      end do

      vormax=0.
      psimax=0.
      do 101 j=1,n2
      jm=jmv(j)
      do 101 i=1,n1
      psimax=max(abs(psi(i,j)),psimax)
      if(abs(vor(i,j)).gt.vormax) then
      vormax=abs(vor(i,j))
      jmax=j
      imax=i
                                    endif
  101 continue
       xmax=xxm(imax,jmax)
       ymax=yym(imax,jmax)
c
c   write the maxima value in time
c
       write(32,721) timew,psimax,vormax,xmax,ymax
       write(6,721) timew,psimax,vormax,xmax,ymax
  721 format(1x,'t=',e11.4,2x,'psimax=',e11.4,2x,
     1      'vormax=',e11.4,2x,'at x,y',2e11.4)
      return
      end   
c************************************************************************
c                                                                       *
c     ********* subrout oudip ******************                       *
c                                                                       *
c************************************************************************
      subroutine oudip
      include 'param.f'
      common/timw/timew
      character*80 namfil
      common/vcar/voc1(m1,m2),voc2(m1,m2)
c
c     velocities  at the cell point i,j
c
      vtmax=0.
      vrmax=0.
      kc=1
      do 4 jc=2,n2m
      jm=jmv(jc)
      do 4 ic=1,n1m
      im=imv(ic)
      voc2(ic,jc)=
     1  (    q2(ic,jc,kc)+q2(imv(ic),jc,kc)
     1      )*0.5/rc(jc)
      voc1(ic,jc)=
     1  (    q1(ic,jc,kc)/rm(jc)+q1(ic,jm,kc)/rm(jm)
     1      )*0.5
      vrmax=max(vrmax,voc2(ic,jc))
      vtmax=max(vtmax,voc1(ic,jc))
    4 continue
      do 5 jc=1,n2,n2m
      do 5 ic=1,n1m
      voc2(ic,jc)=0.
      if(jc.eq.1) voc1(ic,jc)=0.
      if(jc.eq.n2) voc1(ic,jc)=voc1(ic,jc-1)
    5 continue
c 
      do 14 jc=1,n2
      voc2(n1,jc)=voc2(1,jc)
      voc1(n1,jc)=voc1(1,jc)
   14 continue
c
c     Cartesiane velocity components
c
      vxmax=0.
      vymax=0.
      do 15 jc=1,n2
      do 15 ic=1,n1
      vx=voc2(ic,jc)*cos(thetac(ic))-voc1(ic,jc)*sin(thetac(ic))
      vy=-voc2(ic,jc)*sin(thetac(ic))-voc1(ic,jc)*cos(thetac(ic))
      voc2(ic,jc)=vy
      voc1(ic,jc)=vx
      vxmax=max(vxmax,vx)
      vymax=max(vymax,vy)
   15 continue

      itime=nint(timew*10.)
      write(ipfi,82)itime
   82 format(i4.4)
c
c   write a 2D field to be visualized by TURB3D
c
      namfil='vorpsi'//ipfi//'.dat'
      open(13,file=namfil,form='unformatted')
      write(6,*) 'the file vorpsi',ipfi
     1          ,'was written at t=',
     1          timew,'vemax  ',vxmax,vymax
     1          ,'vmax  ',vtmax,vrmax
      do j=1,n2
      vor(n1,j)=vor(1,j)
      psi(n1,j)=psi(1,j)
      enddo
      nfil=13
      rewind(nfil)
      write(nfil) n1,n2
      write(nfil) ros,alx3d,ren,timew
      write(nfil) ((voc1(i,j),i=1,n1),j=1,n2),
     1            ((voc2(i,j),i=1,n1),j=1,n2),
     1            ((psi(i,j),i=1,n1),j=1,n2),
     1            ((vor(i,j),i=1,n1),j=1,n2)
      close(13)
      return
      end
c************************************************************************
c                                                                       *
c    *****   subro outiqu   ********************************              *
c    
c                                                                       *
c************************************************************************
      subroutine outiqu(time)
      include 'param.f'
      numf=90
      do jp=1,njprs
      j=npjp(jp)
      do k=kpi,n3m,kpr
      call quapt(numf,j,k,time)
      enddo
      numf=numf+1
      enddo
      return
      end
c
c  ************************* subrout quapt  **********************
c
c     this subroutine calculates the quantities at the
c     cell center to be written at each time steps
c     for postprocessing time statistics
c
      subroutine quapt(numf,jc,kc,time)
      include 'param.f'
      jp=jc+1
      jm=jc-1
      kp=kpv(kc)
      km=kmv(kc)
      do ic=1,n1m,ipr
      ip=ipv(ic)
      im=imv(ic)
      q2m=(q2(ic,jc,kc)/rc(jc)+q2(ic,jp,kc)/rc(jp))*0.5
      q1m=(q1(ic,jc,kc)+q1(ip,jc,kc) )*0.5/rm(jc)
      q3m=(q3(ic,jc,kc)+q3(ic,jc,kp) )*0.5
      dq2x3=(q2(ic,jc,kc)-q2(ic,jc,km))*dx3/rc(jc)
      dq3x2=(q3(ic,jc,kc)-q3(ic,jm,kc))/(rm(jc)-rm(jm))
      aicjckc=(dq2x3-dq3x2)
      dq2x3=(q2(ic,jp,kc)-q2(ic,jp,km))*dx3/rc(jp)
      dq3x2=(q3(ic,jp,kc)-q3(ic,jc,kc))/(rm(jp)-rm(jc))
      aicjpkc=(dq2x3-dq3x2)
      dq2x3=(q2(ic,jp,kp)-q2(ic,jp,kc))*dx3/rc(jp)
      dq3x2=(q3(ic,jp,kp)-q3(ic,jc,kp))/(rm(jp)-rm(jc))
      aicjpkp=(dq2x3-dq3x2)
      dq2x3=(q2(ic,jc,kp)-q2(ic,jc,kc))*dx3/rc(jc)
      dq3x2=(q3(ic,jc,kp)-q3(ic,jm,kp))/(rm(jc)-rm(jm))
      aicjckp=(dq2x3-dq3x2)
      voam=(aicjckc+aicjpkc+aicjckp+aicjpkp)*0.25
      dq3x1=(q3(ic,jc,kc)-q3(im,jc,kc))*dx1/rm(jc)
      dq1x3=(q1(ic,jc,kc)-q1(ic,jc,km))*dx3/rm(jc)
      ricjckc=(dq3x1-dq1x3)
      dq3x1=(q3(ip,jc,kc)-q3(ic,jc,kc))*dx1/rm(jc)
      dq1x3=(q1(ip,jc,kc)-q1(ip,jc,km))*dx3/rm(jc)
      ripjckc=(dq3x1-dq1x3)
      dq3x1=(q3(ic,jc,kp)-q3(im,jc,kp))*dx1/rm(jc)
      dq1x3=(q1(ic,jc,kp)-q1(ic,jc,kc))*dx3/rm(jc)
      ricjckp=(dq3x1-dq1x3)
      dq3x1=(q3(ip,jc,kp)-q3(ic,jc,kp))*dx1/rm(jc)
      dq1x3=(q1(ip,jc,kp)-q1(ip,jc,kc))*dx3/rm(jc)
      ripjckp=(dq3x1-dq1x3)
      vorm=(ricjckc+ripjckc+ricjckp+ripjckp)*0.25
      dq1x2=(q1(ic,jc,kc)-q1(ic,jm,kc))/(rm(jc)-rm(jm))
      dq2x1=(q2(ic,jc,kc)-q2(im,jc,kc))*dx1/rc(jc)
      xicjckc=(dq1x2-dq2x1)/rc(jc)
      dq1x2=(q1(ip,jc,kc)-q1(ip,jm,kc))/(rm(jc)-rm(jm))
      dq2x1=(q2(ip,jc,kc)-q2(ic,jc,kc))*dx1/rc(jc)
      xipjckc=(dq1x2-dq2x1)/rc(jc)
      dq1x2=(q1(ic,jp,kc)-q1(ic,jc,kc))/(rm(jp)-rm(jc))
      dq2x1=(q2(ic,jp,kc)-q2(im,jp,kc))*dx1/rc(jp)
      xicjpkc=(dq1x2-dq2x1)/rc(jp)
      dq1x2=(q1(ip,jp,kc)-q1(ip,jc,kc))/(rm(jp)-rm(jc))
      dq2x1=(q2(ip,jp,kc)-q2(ic,jp,kc))*dx1/rc(jp)
      xipjpkc=(dq1x2-dq2x1)/rc(jp)
      voxm=(xicjckc+xicjpkc+xipjckc+xipjpkc)*0.25
      write(numf)time,kc,ic,q1m,q2m,q3m,pr(ic,jc,kc),
     1           voam,vorm,voxm
      enddo
      return
      end
c************************************************************************
c                                                                       *
c    *****   subro outh   ********************************              *
c   some statistics are written during the simualtions to
c   check the goodness of the simulation
c   These are written and calculate only when iav=1
c                                                                       *
c************************************************************************
      subroutine outh(time,nav,ntime,cflm)        
      include 'param.f'
      character*60 filpvp,filvv,filuvt,filvwt,filvtv,filvtd,filvtn
      character*60 fivevo,filhel,filpro  
      common/utapr/utap
      character*60 filebu1,filebu2,filebu3,filebu4,filebu5
      dimension skp(4),flp(4),vmp(3),tstp(6)
      dimension bskp(3),bflp(3),bmp(3),btstp(6)
      dimension vistzr(m2),visttr(m2)
      common/nbal/nba
      common/enetot/enet
      volto=1./(pi*alx3d)
      n1mh=n1m/2
      if(iav.eq.1) then
       nav=nav+1
      dp3mo=dp3mo+dp3ns
      dp3mav=dp3mo/nav
      cfavp=(abs(dp3mav/2.))
      utap=sqrt(abs(cfavp))
                   endif
      n1vm=n1m
      n2vm=n2m
      n3vm=n3m
c
c   fluctuating velocity components
c
      call veltur
      if(iav.eq.1) then
c         print *,'calcolo taver ', nav
      call taver
                   endif

c
c     write(16,*)' in outh nba,',nba,' nav',nav,'t=',time
      enavp=enavo/nav
      if(nav.gt.1) then
      if(icfl.eq.1) then
      if(amod(time,tpin).lt.dt) go to 306
      go to 305
                     else
      if(mod(ntime,npouth).eq.0) then
           go to 306
                               endif
                     endif
      go to 305
  306 continue
      write(6,*)' in outh nav',nav,
     1    'averages for vertic. transp'
      itim=nint(time+dt)
      write(pntim,77) itim
   77 format(i4.4)

c
c   the profiles are written
c
      filve='oupivm.'//pntim
      filuvp='oupiuvp.'//pntim
      filuvt='oupiuvt.'//pntim
      filvwt='oupivwt.'//pntim
      filnsp='oupinsp.'//pntim
      filvmp='oupivmp.'//pntim
      filns='oupins.'//pntim
      filts='oupits.'//pntim
      filsk='oupisk.'//pntim
      filfl='oupifl.'//pntim
      if(ren.lt..1e07) then
      do j=2,n2m
      vistzr(j)=(vmeo(3,j)-vmeo(3,j-1))*dx2/g2rc(j)/nav/ren
      visttr(j)=(vmeo(1,j)/rm(j)-vmeo(1,j-1)/rm(j-1))
     1          *rc(j)*dx2/g2rc(j)/nav/ren
      enddo
      vistzr(1)=0.
      vistzr(n2)=(-vmeo(3,n2-1))*dx2/g2rc(n2)/nav/ren
      visttr(1)=0.
      visttr(n2)=(-vmeo(1,n2-1)/rm(n2-1))
      open(81,file=filvmp,form='formatted')
      open(89,file=filuvt,form='formatted')
      open(87,file=filvwt,form='formatted')
      open(82,file=filuvp,form='formatted')
      open(83,file=filnsp,form='formatted')
                       endif
      open(51,file=filve,form='formatted')
      open(55,file=filns,form='formatted')
      open(56,file=filts,form='formatted')
      open(57,file=filsk,form='formatted')
      open(58,file=filfl,form='formatted')
      do  j=1,n2m
      pmp=pmeo(j)/nav
      psp=psto(j)/nav
      do  l=1,4
      skp(l)=skeo(l,j)/nav
      flp(l)=flao(l,j)/nav
      enddo
      do  l=1,3
      vmp(l)=vmeo(l,j)/nav
      enddo
      do  l=1,3
      tstp(l)=(tursto(l,j)/nav)
      enddo
      do  l=4,6
      tstp(l)=tursto(l,j)/nav
      enddo
      drl=rc(j+1)-rc(j)
      suu=skp(1)/tstp(1)**3.
      svv=skp(2)/tstp(2)**3.
      sww=skp(3)/tstp(3)**3.
      spp=skp(4)/psp**(1.5) 
      fuu=flp(1)/tstp(1)**4.
      fvv=flp(2)/tstp(2)**4.
      fww=flp(3)/tstp(3)**4.
      fpp=flp(4)/psp**2.
      yd=1.-rm(j) 
      umpl=vmp(3)/utap
      uupl=tstp(1)/utap
      vvpl=tstp(2)/utap
      wwpl=tstp(3)/utap
      uvpl=-tstp(5)/utap/utap
      uwpl=-tstp(4)/utap/utap
      vwpl=-tstp(6)/utap/utap
      ypl=yd*utap*ren
      uplo=1./.41*alog(ypl)+5.5
      toszrv=(vistzr(j)+vistzr(j+1))*0.5/utap/utap
      tostrv=(visttr(j)+visttr(j+1))*0.5/utap/utap
      toszrt=tstp(5)/utap/utap
      tostrt=tstp(4)/utap/utap
      totszr=(toszrv+toszrt) 
      totstr=(tostrv+tostrt) 
      if(ren.lt..1e07) then
      write(81,612)ypl,umpl,uplo
      write(83,612)ypl,uupl,vvpl,wwpl
      write(82,612)yd,uvpl,uwpl,vwpl
      write(87,612)yd,toszrv,toszrt,totszr
      write(89,612)yd,tostrv,tostrt,totstr
                       endif
      write(51,612)yd,(vmp(l),l=1,3),pmp
      write(55,612)yd,(tstp(l),l=1,3)
      write(56,612)yd,(tstp(l),l=4,6)
      write(57,612)yd,suu,svv,sww,spp
      write(58,612)yd,fuu,fvv,fww,fpp
  611 continue
  612 format(1x,e12.4,2x,9(1x,e12.5)) 
      enddo
      if(ren.lt..1e07) then
      close(81)
      close(82)
      close(83)
      close(87)
      close(89)
      close(84)
      close(85)
      close(86)
      close(88)
                       endif
      close(52)
      close(51)
      close(55)
      close(56)
      close(57)
      close(58)
  983 format(3x,e11.4,i5,9(1x,e12.5))
  783 format(3x,10(1x,e12.5))
  785 format(3x,7(1x,e12.5),3x,i5)
  305         continue     
      if(ren.lt..1e07) then
      retap=utap*ren
      uitap=vmeo(3,1)/nav/utap
      reucl=ren*vmeo(3,1)/nav
                       endif
      write(59,785)time,uitap,retap,utap,cfavp,dp3mav,reucl
     1             ,nav
			endif
      utat=sqrt(abs(dp3ns/2.))
c     write(32,783)time,(vit(l),l=1,4),enav,utat
c    1 ,dp3ns,cflm
      write(49,783)time,utat
      write(6,783)time,vit(2),(vit(l),l=1,3,2),enav,utat,utap
     1 ,cflm,dt
      write(50,783)time,vit(2),(vit(l),l=1,3,2),enav,utat
     1 ,dp3ns,cflm
      return                                                            
      end                                                               
c************************************************************************
c                                                                       *
c    *****   subro outhth ********************************              *
c    this write pressure gradient and turb. kinetic energy
c    friction velocity and total mass at certain time intervals
c    to check the right time evolution
c                                                                       *
c************************************************************************
      subroutine outhth(time,nap,ntime,cflm)        
      include 'param.f'
      common/dp3t/dp3th
      common/enetot/enet
      call totvel
      call energy
      volto=1./(pi*alx3d)
      utat=sqrt(abs(dp3ns/2.))
      if(iav.eq.1) then
       nap=nap+1
      cfo=cfnp+cfo
      dp3th=dp3th+dp3ns
      cfav=cfo/nap
      dp3mav=dp3th/nap
      cfavp=(abs(dp3mav/2.))
      utap=sqrt(abs(cfavp))
      write(49,783)time,utat
      if(ren.le.1.e08)  then
      write(6,783)time,(vit(l),l=1,3,2),dp3ns,enav,utap
     1 ,utat,nap,dt
      write(32,783)time,(vit(l),l=1,3,2),dp3ns,enav,utap
     1 ,utat,nap,dt
                       else
      write(6,784)time,(vit(l),l=1,3,2),dp3ns,enet,utat
     1 ,cflm,dt
      write(32,784)time,(vit(l),l=1,3,2),dp3ns,enet,utat
     1 ,cflm,dt
                       endif
                  else
      if(ren.le.1.e08)  then
      write(49,783)time,utat
      write(6,784)time,(vit(l),l=1,3,2),dp3ns,(vmax(l),l=1,3),enav,dt
                       else
      write(6,784)time,vit(3),vmed(3,n2m),(vmax(l),l=1,3)
     1 ,enet,enav,dt
      write(32,784)time,vit(3),vmed(3,n2m),(vmax(l),l=1,3)
     1 ,enet,enav,dt
                       endif
                  endif
  783 format(3x,7(1x,e12.5),3x,i5,e11.4)
  784 format(3x,8(1x,e12.5),5x,e11.4)
      return                                                            
      end                                                               
c************************************************************************
c                                                                       *
c  ****************************** subrout vmaxv **********************  *
c                                                                       *
c************************************************************************
      subroutine vmaxv(n1mv,n2mv,n3mv)
      include 'param.f'
c
c     find the maximum velocities in order to check convergence or
c     to derive stability conditions.
c
      vca=0.
c
      do 311 l=1,ndv
      vmax(l)=0.
      do k=1,n3mv
      do j=2,n2mv
      do i=1,n1mv
      if(l.eq.2) vca=q2(i,j,k)/rc(j)
      if(abs(vca).ge.vmax(l)) then
      vmax(l)=abs(vca)
                              endif
      enddo
      enddo
      enddo
      do k=1,n3mv
      do j=1,n2mv
      do i=1,n1mv
      if(l.eq.1) vca=q1(i,j,k)/rm(j)
      if(l.eq.3) vca=q3(i,j,k)
      if(abs(vca).ge.vmax(l)) then
      vmax(l)=abs(vca)
                              endif
      enddo
      enddo
      enddo
  311 continue
      return
      end
c************************************************************************
c                                                                       *
c  ****************************** subrout contwr ********************** *
c     print flow fields                                                 *
c                                                                       *
c************************************************************************
      subroutine contwr(ntime,time,ntt,nav)             
      include 'param.f'
      common/irety/ireq2,iwrq2
      character*18 pnamfil
      character*80 namfil
      do i=1,n1m
      do j=1,n2
      q1(i,j,n3)=q1(i,j,1)
      q2(i,j,n3)=q2(i,j,1)
      q3(i,j,n3)=q3(i,j,1)
      pr(i,j,n3)=pr(i,j,1)
      enddo
      enddo
      do j=1,n2
      do k=1,n3
      q1(n1,j,k)=q1(1,j,k)
      q2(n1,j,k)=q2(1,j,k)
      q3(n1,j,k)=q3(1,j,k)
      pr(n1,j,k)=pr(1,j,k)
      enddo
      enddo
      prma=-1000.
      prmi=+1000.
      do kc=1,n3l
      do ic=1,n1l
      do jc=n2l,1,-1
      prma=max(prma,pr(ic,jc,kc))
      prmi=min(prmi,pr(ic,jc,kc))
      enddo
      enddo
      enddo
      itime=nint(time)
      write(ipfi,82)itime
   82 format(i4.4)
      if(iav.eq.1) then
      namfil='field'//ipfi//'.dat'
      pnamfil='field'//ipfi//'.dat'
      open(13,file=namfil,form='unformatted')
                   else
      pnamfil=filcnw
      open(13,file=filcnw,form='unformatted')
                   endif
      write(6,*) pnamfil,'written at t=',
     1          time, ' prma mi=',prma,prmi
      nfil=13
      rewind(nfil)
      write(nfil) n1,n2,n3
      write(nfil) ros,alx3d,ren,time
      if(iwrq2.eq.1) then
c
c    large memory occupancy
c
      write(nfil) (((q1(i,j,k),i=1,n1),j=1,n2),k=1,n3),           
     1            (((q2(i,j,k),i=1,n1),j=1,n2),k=1,n3),        
     1            (((q3(i,j,k),i=1,n1),j=1,n2),k=1,n3),
     1            (((pr(i,j,k),i=1,n1),j=1,n2),k=1,n3)
                     else
c
c    reduced memory occupancy
c    the pressure is not necessary for restarting files
c    even for post processing can be saved but then
c    the advancement of a time step should be performed
c
      write(nfil) (((q1(i,j,k),i=1,n1),j=1,n2),k=1,n3),           
     1            (((q3(i,j,k),i=1,n1),j=1,n2),k=1,n3),
     1            (((pr(i,j,k),i=1,n1),j=1,n2),k=1,n3)
                     endif
      write(nfil) ntime,ntt,nav
      close(nfil)                                                       
      return                                                            
      end                                                               
c************************************************************************
c                                                                       *
c  ****************************** subrout inirea ********************** *
c     read input flow fields                                            *
c                                                                       *
c************************************************************************
      subroutine inirea(ntii,time,ntt,ncount,nav)
      include 'param.f'
      dimension q1ppr(m2),q2ppr(m2),q3ppr(m2)
      dimension q1rmo(m2),q2rmo(m2),q3rmo(m2)
      dimension q1rmn(m2),q2rmn(m2),q3rmn(m2)
      common/avgin/avpscn(m2),vnew(m2)
      common/irety/ireq2,iwrq2
      open(13,file=filcnr,form='unformatted')
      nfil=13                                                           
      read(nfil) n1l,n2l,n3l                                        
      n1lm=n1l-1
      n2lm=n2l-1
      n3lm=n3l-1
      n1vm=n1lm
      n2vm=n2lm
      n3vm=n3lm
      read(nfil) a1,a1,a1,time
      write(6,*)'ireq2=',ireq2,
     1      ' legge da ',filcnr
      write(6,*)'  n1l,n2l,n3l,time ',n1l,n2l,n3l,time
c
c   large memory occupancy
c
      if(ireq2.eq.1) then
      read(nfil)  (((q1(i,j,k),i=1,n1l),j=1,n2l),k=1,n3l),           
     1            (((q2(i,j,k),i=1,n1l),j=1,n2l),k=1,n3l),        
     1            (((q3(i,j,k),i=1,n1l),j=1,n2l),k=1,n3l),
     1            (((pr(i,j,k),i=1,n1l),j=1,n2l),k=1,n3l)
      q1mal=0.
      q2mal=0.
      q3mal=0.
      do j=1,n2l
      do i=1,n1l
      do k=1,n3l
      q1mal=max(q1mal,abs(q1(i,j,k)))
      q2mal=max(q2mal,abs(q2(i,j,k)))
      q3mal=max(q3mal,abs(q3(i,j,k)))
      enddo
      enddo
      enddo
      write(6,*)' maxima ',q1mal,q2mal,q3mal

                     else
c
c   reduced memory occupancy
c
      read(nfil)  (((q1(i,j,k),i=1,n1l),j=1,n2l),k=1,n3l),           
     1            (((q3(i,j,k),i=1,n1l),j=1,n2l),k=1,n3l),
     1            (((pr(i,j,k),i=1,n1l),j=1,n2l),k=1,n3l)
      prma=-1000.
      prmi=+1000.
      q2ma=-1000.
      q2mi=+1000.
      do kc=1,n3lm
      do ic=1,n1lm
      do jc=n2l,1,-1
      prma=max(prma,pr(ic,jc,kc))
      prmi=min(prmi,pr(ic,jc,kc))
      q2ma=max(q2ma,q2(ic,jc,kc))
      q2mi=min(q2mi,q2(ic,jc,kc))
      enddo
      enddo
      enddo
c
c   q2 from divergency q2 is supposed 0 at the wall
c
      do kc=1,n3lm
      do ic=1,n1lm
      q2(ic,n2,kc)=0.
      enddo
      enddo
      do jc=n2lm,1,-1
      jp=jc+1
      usrnu1=dx1
      usrnu2=dx2/g2rm(jc)*rm(jc)
      usrnu3=dx3*rm(jc)*rm(jc)
      do kc=1,n3lm
      kp=kpv(kc)
      do ic=1,n1lm
      ip=ipv(ic)
      dqcap=+(q1(ip,jc,kc)-q1(ic,jc,kc))*usrnu1
     1      +(q3(ic,jc,kp)-q3(ic,jc,kc))*usrnu3
      q2(ic,jc,kc)=q2(ic,jp,kc)+dqcap/usrnu2
      enddo
      enddo
      enddo
      q2madi=-1000.
      q2midi=+1000.
      do kc=1,n3lm
      do ic=1,n1lm
      do jc=n2l,1,-1
      q2madi=max(q2madi,q2(ic,jc,kc))
      q2midi=min(q2midi,q2(ic,jc,kc))
      enddo
      enddo
      enddo
      write(6,*)'prma q2ma e prmi q2mi',prma,q2madi,prmi,q2midi
                         endif
      if(ireset.eq.0) then

c
c  a field previously generated in a pipe of length L
c  is used to restart the simulation in a pipe of 
c  length 2*L
c
      read(nfil) ntii,ntt
      nav=0
      time=0.
      write(6,*)' ntii,in inirea',ntii,ntt,nav,time
      ntii=0
      ntt=0
      if(n3m.eq.2*n3lm) then
c
c   copies the velocity field in the first section in the
c   added section
c
      do j=1,n2l
      q1rmo(j)=0.
      q2rmo(j)=0.
      q3rmo(j)=0.
      do i=1,n1l
      do k=1,n3l
      kk=k+n3lm
      q1(i,j,kk)=q1(i,j,k)
      q2(i,j,kk)=q2(i,j,k)
      q3(i,j,kk)=q3(i,j,k)
      q1rmo(j)=max(q1rmo(j),abs(q1(i,j,k)))
      q2rmo(j)=max(q2rmo(j),abs(q2(i,j,k)))
      q3rmo(j)=max(q3rmo(j),abs(q3(i,j,k)))
      enddo
      enddo
      enddo
c
c   add a random disturbance to the field to distroy  the
c   periodicity and evaluates the averaged random disturbance
c
      do j=1,n2l
      q1ppr(j)=0.
      q3ppr(j)=0.
      do i=1,n1l
      do k=1,n3
c     q1per=vper*(-1.+2.*ranf())
      q1per=vper*(-1.+2.*rand())
      q1ppr(j)=q1ppr(j)+q1per
c     q3per=vper*(-1.+2.*ranf())
      q3per=vper*(-1.+2.*rand())
      q3ppr(j)=q3ppr(j)+q3per
      q1(i,j,k)=q1(i,j,k)*(1.+q1per)
      q3(i,j,k)=q3(i,j,k)*(1.+q3per)
      enddo
      enddo
      q1ppr(j)=q1ppr(j)/float(n1l*n3)
      q3ppr(j)=q3ppr(j)/float(n1l*n3)
      enddo
c
c   subctracts  the averaged random disturbance
c   to have zero mass of the disturbance
c
      do j=1,n2l
      do i=1,n1l
      do k=1,n3
      q1(i,j,k)=q1(i,j,k)-q1ppr(j)
      q3(i,j,k)=q3(i,j,k)-q3ppr(j)
      enddo
      enddo
      enddo
      do j=1,n2l
      q1rmn(j)=0.
      q2rmn(j)=0.
      q3rmn(j)=0.
      do i=1,n1m
      do k=1,n3m
      q1rmn(j)=max(q1rmn(j),abs(q1(i,j,k)))
      q3rmn(j)=max(q3rmn(j),abs(q3(i,j,k)))
      enddo
      enddo
      enddo
      write(78,*)' old and new maxima of velocity'
      do j=1,n2l
      write(78,121)j,q1rmo(j),q1rmn(j),q3rmo(j),q3rmn(j)
  121 format(3x,i4,3x,6e12.5)
      enddo
                   endif
                      else
      read(nfil) ntii,ntt,nav
      write(6,*)' ntii,in inirea',ntii,ntt,nav,time
                      endif
      open(78,file='inico.out')
      write(78,*)' in inirea',ntii,ntt,time
      avgnl=1./float(n1lm*n3lm)
c     call vmaxv(n1vm,n2vm,n3vm)
c     write(78,783)time,(vmax(l),l=1,3),(vmaxo(l),l=1,3)
  783 format(3x,'ini ',7(1x,e12.4))
      if(n2m.eq.n2lm) then
      dx1o=2.*pi/float(n1lm)
      dx3o=alx3/float(n3lm)
      dx1o=1./dx1o
      dx3o=1./dx3o
      volto=1./(pi*alx3d)
      vit3ol=0.
      vit2ol=0.
      vit1ol=0.
      do jc=1,n2m
      jp=jc+1
      voz=rm(jc)*g2rm(jc)/(dx1o*dx2*dx3o)
      do kc=1,n3lm
      do ic=1,n1lm
      vit3ol=q3(ic,jc,kc)*voz+vit3ol
      vit1ol=q1(ic,jc,kc)*voz+vit1ol
      vit2ol=(q2(ic,jp,kc)+q2(ic,jc,kc))*voz+vit2ol
      enddo
      enddo
      enddo
      vit3ol=vit3ol*volto
      vit2ol=vit2ol*volto
      vit1ol=vit1ol*volto
                     endif
c
c    here the different options to perform the simulations
c    from field generated by coarse or fine simulations
c    interpoaltions in theta and x3 are bilinear and only
c    grids with ratio 2 are considered
c    the routine with names starting with in fine to coarse
c    the routine with names starting with ex coarse to fine 
c    in each routine only q1 and q3 evaluated
c    q2 obtained by divergence
C
      if(n1lm.eq.2*n1m.and.n3lm.eq.2*n3m) then
      write(6,*)' from fine grid n1lm=',n1lm,' n3lm=',n3lm,
     1          ' to coarse grid n1m=',n1m,' n3m=',n3m
      call inthze(n1lm,n3lm)
						    endif
      if(n1lm.eq.2*n1m.and.n3lm.eq.n3m) then
      write(6,*)' from fine grid n1lm=',n1lm,' n3lm=',n3lm,
     1          ' to coarse grid n1m=',n1m,' n3m=',n3m
      call inthet(n1lm)
          				    endif
      if(n1lm.eq.n1m.and.n3lm.eq.2*n3m) then
      write(6,*)' from fine grid n1lm=',n1lm,' n3lm=',n3lm,
     1          ' to coarse grid n1m=',n1m,' n3m=',n3m
      call inzeta(n3lm)
					    endif
      if(n1m.eq.2*n1lm.and.n3m.eq.n3lm) then
      write(6,*)' from coarse grid n1lm=',n1lm,' n3lm=',n3lm,
     1          ' to fine grid n1m=',n1m,' n3m=',n3m
      call exthet(n1lm,n3lm)
						    endif
      if(n1m.eq.n1lm.and.n3m.eq.2*n3lm) then
      write(6,*)' from coarse grid n1lm=',n1lm,' n3lm=',n3lm,
     1          ' to fine grid n1m=',n1m,' n3m=',n3m
      call exzeta(n1lm,n3lm)
						    endif
      if(n1m.eq.2*n1lm.and.n3m.eq.2*n3lm) then
      write(6,*)' from coarse grid n1lm=',n1lm,' n3lm=',n3lm,
     1          ' to fine grid n1m=',n1m,' n3m=',n3m
      call exthze(n1lm,n3lm)
						    endif
c
c   the distribution of radial coordinates is changed
c
      if(n2m.ne.n2lm)  then
      n2o=n2l
      call cordino(n2o)
      dx2o=float(n2lm)
      dx1o=2.*pi/float(n1lm)
      dx3o=alx3/float(n3lm)
      dx1o=1./dx1o
      dx3o=1./dx3o
      volto=1./(pi*alx3d)
      vit3ol=0.
      vit2ol=0.
      vit1ol=0.
      do jc=1,n2lm
      jp=jc+1
      voz=rmo(jc)*g2rmo(jc)/(dx1o*dx2o*dx3o)
      do kc=1,n3lm
      do ic=1,n1lm
      vit3ol=q3(ic,jc,kc)*voz+vit3ol
      vit1ol=q1(ic,jc,kc)*voz+vit1ol
      vit2ol=(q2(ic,jp,kc)+q2(ic,jc,kc))*voz+vit2ol
      enddo
      enddo
      enddo
      vit3ol=vit3ol*volto
      vit2ol=vit2ol*volto
      vit1ol=vit1ol*volto
      call divgco(qmaxo,qtoto,n2o,dx1o,dx3o,n1lm,n3lm)                             
c
c   spline interpolation
c
      call intrpr(n2o)
c
c   divergence check
c
      call divuck(qmax,qtot)                             
          write(6,900)qmaxo,qtoto,qmax,qtot
  900 format(3x,' old dmax = ',e11.4,' old dtot = ',e11.4
     1     /,3x,' new dmax = ',e11.4,' new dtot = ',e11.4)
                       endif
      vit3ne=0.
      vit2ne=0.
      vit1ne=0.
      do jc=1,n2m
      jp=jc+1
      voz=rm(jc)*g2rm(jc)/(dx1*dx2*dx3)
      do kc=1,n3m
      do ic=1,n1m
      vit1ne=q1(ic,jc,kc)*voz+vit1ne
      vit3ne=q3(ic,jc,kc)*voz+vit3ne
      vit2ne=(q2(ic,jp,kc)+q2(ic,jc,kc))*voz+vit2ne
      enddo
      enddo
      enddo
      vit3ne=vit3ne*volto
      vit2ne=vit2ne*volto
      vit1ne=vit1ne*volto
      write(78,*) 'in inirea q3 total old = ',vit3ol,
     1           '  new=',vit3ne
      write(78,*) 'in inirea q2 total old = ',vit2ol,
     1           '  new=',vit2ne
      write(78,*) 'in inirea q1 total old = ',vit1ol,
     1           '  new=',vit1ne
      write(6,*) 'in inirea q3 total old = ',vit3ol,
     1           '  new=',vit3ne
      write(6,*) 'in inirea q2 total old = ',vit2ol,
     1           '  new=',vit2ne
      write(6,*) 'in inirea q1 total old = ',vit1ol,
     1           '  new=',vit1ne
      avgn=1./(float(n1m*n3m))
      if (ireset.eq.1) then                                             
      ntii=0                                                            
      endif                                                             
      write(78,*) 'in inirea denmax e min = ',ddd,ddm
      if(iav.eq.0) then
      call dp3ini
      write(6,*) ' dp3ns from the field=',dp3ns
      nav=0
           call outh(time,0,ntime,0.)
                     endif
      if(iav.eq.0) write(78,*) 'end inico.out'
      if(iav.eq.0) close(78)
      return                                                            
      end                                                               
c************************************************************************
c                                                                       *
c ****************************** subrout dp3ini  ********************** *
c  this subroutine performs the calculation of the pressure gradients 
c  of the restarting file
c                                                                       *
c************************************************************************
      subroutine dp3ini
      include 'param.f'
      common/cor3j/ap3j(m2),ac3j(m2),am3j(m2)
      write(78,*) 'in dp3ini ',ren 
      alre=1./ren                                                       
c     c
c   add second derivatie in r
c
      s3tot=0.
      do kc=1,n3m
            do jc=2,n2m-1
            jm=jmv(jc)                                                       
            jp=jpv(jc)                                                       
                  do ic=1,n1m
c
c   22 second derivatives of q3
c
      dq32=(q3(ic,jp,kc)*ap3j(jc)
     1     -q3(ic,jc,kc)*(ap3j(jc)+am3j(jc))
     1     +q3(ic,jm,kc)*am3j(jc))*alre
      s3tot=s3tot+dq32*volz(jc)
                  enddo  
            enddo                                                             
      enddo                                                             
c
c   22 second derivatives of q3  at r=n2m
c
      jc=n2m
      jm=jc-1                                                      
      do kc=1,n3m
                  do ic=1,n1m
c
c   22 second derivatives of q3
c
      dq32=(-q3(ic,jc,kc)*(ap3j(jc)+ac3j(jc))
     1      +q3(ic,jm,kc)*am3j(jc))*alre
      s3tot=s3tot+dq32*volz(jc)
                  enddo                                                             
      enddo                                                             
c     c
c   add second derivatie in r at r=0
c
      jc=1
      jp=jc+1                                                      
      do kc=1,n3m
                  do ic=1,n1m
c
c   22 second derivatives of q3
c
      dq32=(q3(ic,jp,kc)*ap3j(jc)
     1     -q3(ic,jc,kc)*ap3j(jc))*alre
      s3tot=s3tot+dq32*volz(jc)
                  enddo                                                             
      enddo                                                             
       dp3ns=s3tot/(pi*alx3d)
      write(78,*) 'in dp3ini ',dp3ns 
      return                                                            
      end                                                               
c************************************************************************
c                                                                       *
c  ****************************** subrout inthze ********************** *
c     interpolates in coarser grid both in Theta and Zeta direction     *
c                                                                       *
c************************************************************************
      subroutine inthze(n1lm,n3lm) 
      include 'param.f'
      do j=1,n2m
      do i=1,n1m
      inm=2*(i-1)
      if(i.eq.1) inm=2*n1m
      inc=2*(i-1)+1
      inp=2*i
      do k=1,n3m
      knm=2*(k-1)+1
      knp=2*k
      dq(i,j,k)=(q1(inp,j,knp)+q1(inp,j,knm)
     1           +q1(inm,j,knp)+q1(inm,j,knm)
     1       +2.*(q1(inc,j,knp)+q1(inc,j,knm)) )/8.
      pr(i,j,k)=0.
      enddo
      enddo
      enddo
      do j=1,n2m
      do k=1,n3m
      knm=2*(k-1)
      if(k.eq.1) knm=2*n3m
      knc=2*(k-1)+1
      knp=2*k
      do i=1,n1m
      inm=2*(i-1)+1
      inp=2*i
      qcap(i,j,k)=(q3(inp,j,knp)+q3(inp,j,knm)
     1           +q3(inm,j,knp)+q3(inm,j,knm)
     1       +2.*(q3(inm,j,knc)+q3(inp,j,knc)) )/8.
      enddo
      enddo
      enddo
      call q2div
      return                                                            
      end                                                               
c************************************************************************
c                                                                       *
c  ****************************** subrout inthet ********************** *
c     interpolates in coarser grid  in Theta  direction                 *
c                                                                       *
c************************************************************************
      subroutine inthet(n1lm) 
      include 'param.f'
      do j=1,n2m
      do i=1,n1m
      inm=2*(i-1)
      if(i.eq.1) inm=2*n1m
      inc=2*(i-1)+1
      inp=2*i
      do k=1,n3m
      dq(i,j,k)=(q1(inp,j,k)+q1(inm,j,k)+2.*q1(inc,j,k) )/4.
      pr(i,j,k)=0.
      enddo
      enddo
      enddo
      do j=1,n2m
      do k=1,n3m
      do i=1,n1m
      inm=2*(i-1)+1
      inp=2*i
      qcap(i,j,k)=(q3(inp,j,k)+q3(inp,j,k))/2.
      enddo
      enddo
      enddo
      call q2div
      return                                                            
      end                                                               
c************************************************************************
c                                                                       *
c  ****************************** subrout q2div ********************** *
c    evaluates q2 from div q=0                                          *
c                                                                       *
c************************************************************************
      subroutine q2div        
      include 'param.f'
      do kc=1,n3m                                                    
      kp=kpv(kc)                                                        
      do ic=1,n1m                                                   
      ip=ipv(ic)                                                        
      dph(ic,1,kc)=0.
      do jc=1,n2m                                                    
      jp=jc+1                                                           
      usrnu1=dx1/rm(jc)**2
      usrnu2=dx2/g2rm(jc)/rm(jc)
      dqcap=(dq(ip,jc,kc)-dq(ic,jc,kc))*usrnu1
     1     +(qcap(ic,jc,kp)-qcap(ic,jc,kc))*dx3                          
      dph(ic,jp,kc)=dph(ic,jc,kc)-dqcap/usrnu2
      enddo
      enddo
      enddo
      do kc=1,n3m                                                    
      do ic=1,n1m                                                   
      do jc=1,n2                                                     
      q2(ic,jc,kc)=dph(ic,jc,kc)
      enddo
      enddo
      enddo
      do kc=1,n3m                                                    
      do ic=1,n1m                                                   
      do jc=1,n2m
      q1(ic,jc,kc)=dq(ic,jc,kc)
      q3(ic,jc,kc)=qcap(ic,jc,kc)
      enddo
      enddo
      enddo
      return                                                            
      end                                                               
c************************************************************************
c                                                                       *
c  ****************************** subrout inzeta ********************** *
c     interpolates in coarser grid in  Zeta direction                   *
c                                                                       *
c************************************************************************
      subroutine inzeta(n3lm) 
      include 'param.f'
      do j=1,n2m
      do i=1,n1m
      do k=1,n3m
      knm=2*(k-1)+1
      knp=2*k
      dq(i,j,k)=(q1(i,j,knp)+q1(i,j,knm))/2.
      pr(i,j,k)=0.
      enddo
      enddo
      enddo
      do j=1,n2m
      do k=1,n3m
      knm=2*(k-1)
      if(k.eq.1) knm=2*n3m
      knc=2*(k-1)+1
      knp=2*k
      do i=1,n1m
      qcap(i,j,k)=(q3(i,j,knp)+q3(i,j,knm)+2.*q3(i,j,knc))/4.
      enddo
      enddo
      enddo
      call q2div
      return                                                            
      end                                                               
c************************************************************************
c                                                                       *
c  ****************************** subrout exthze ********************** *
c     extrapolates fro coarse grid in Theta and Zeta direction          *
c     values of q1 and q3 components in a fine grid in Theta and Zeta   *
c                                                                       *
c************************************************************************
      subroutine exthze(n1lm,n3lm) 
      include 'param.f'
      do j=1,n2m
      do in=1,n1lm
      ifn=2*in-1
      ifp=2*in
      ip=in+1
      if(in.eq.n1lm) ip=1
      do kn=1,n3lm
      kp=kn+1
      kfn=2*kn
      kfp=2*kn+1
      if(kn.eq.n3lm) then
      kfp=1
      kp=1
                    endif
      dq(ifn,j,kfn)=(3.*q1(in,j,kn)+1.*q1(in,j,kp))/4.
      dq(ifn,j,kfp)=(1.*q1(in,j,kn)+3.*q1(in,j,kp))/4.
      dq(ifp,j,kfn)=(3.*q1(in,j,kn)+1.*q1(in,j,kp))/8.
     1              +(3.*q1(ip,j,kn)+1.*q1(ip,j,kp))/8.
      dq(ifp,j,kfp)=(1.*q1(in,j,kn)+3.*q1(in,j,kp))/8.
     1              +(1.*q1(ip,j,kn)+3.*q1(ip,j,kp))/8.
      pr(ifn,j,kfn)=0.
      pr(ifn,j,kfp)=0.
      pr(ifp,j,kfn)=0.
      pr(ifp,j,kfp)=0.
      enddo
      enddo
      enddo
      do j=1,n2m
      do kn=1,n3lm
      kfn=2*kn-1
      kfp=2*kn
      kp=kn+1
      if(kn.eq.n3lm) kp=1
      do in=1,n1lm
      ip=in+1
      ifn=2*in
      ifp=2*in+1
      if(in.eq.n1lm) then
      ifp=1
      ip=1
                    endif
      qcap(ifn,j,kfn)=(3.*q3(in,j,kn)+1.*q3(ip,j,kn))/4.
      qcap(ifp,j,kfn)=(1.*q3(in,j,kn)+3.*q3(ip,j,kn))/4.
      qcap(ifn,j,kfp)=(3.*q3(in,j,kn)+1.*q3(ip,j,kn))/8.
     1              +(3.*q3(in,j,kp)+1.*q3(ip,j,kp))/8.
      qcap(ifp,j,kfp)=(1.*q3(in,j,kn)+3.*q3(ip,j,kn))/8.
     1              +(1.*q3(in,j,kp)+3.*q3(ip,j,kp))/8.
      enddo
      enddo
      enddo
      call q2div
      return                                                            
      end                                                               
c************************************************************************
c                                                                       *
c  ****************************** subrout exthet ********************** *
c     extrapolates fro coarse grid in Theta and Zeta direction          *
c     values of q1 and q3 components in a fine grid in Theta    *
c                                                                       *
c************************************************************************
      subroutine exthet(n1lm,n3lm) 
      include 'param.f'
      do j=1,n2m
      do in=1,n1lm
      ifn=2*in-1
      ifp=2*in
      ip=in+1
      if(in.eq.n1lm) ip=1
      do k=1,n3m
      dq(ifn,j,k)=q1(in,j,k)
      dq(ifp,j,k)=(q1(in,j,k)+q1(ip,j,k))/2.
      pr(ifn,j,k)=0.
      pr(ifp,j,k)=0.
      enddo
      enddo
      enddo
      do j=1,n2m
      do k=1,n3m
      do in=1,n1lm
      ip=in+1
      ifn=2*in
      ifp=2*in+1
      if(in.eq.n1lm) then
      ifp=1
      ip=1
                    endif
      qcap(ifn,j,k)=(3.*q3(in,j,k)+1.*q3(ip,j,k))/4.
      qcap(ifp,j,k)=(1.*q3(in,j,k)+3.*q3(ip,j,k))/4.
      enddo
      enddo
      enddo
      call q2div
      return                                                            
      end                                                               
c************************************************************************
c                                                                       *
c  ****************************** subrout exzeta ********************** *
c     extrapolates fro coarse grid in Theta and Zeta direction          *
c     values of q1 and q3 components in a fine grid in  Zeta   *
c                                                                       *
c************************************************************************
      subroutine exzeta(n1lm,n3lm) 
      include 'param.f'
      do j=1,n2m
      do i=1,n1m
      do kn=1,n3lm
      kp=kn+1
      kfn=2*kn
      kfp=2*kn+1
      if(kn.eq.n3lm) then
      kfp=1
      kp=1
                    endif
      dq(i,j,kfn)=(3.*q1(i,j,kn)+1.*q1(i,j,kp))/4.
      dq(i,j,kfp)=(1.*q1(i,j,kn)+3.*q1(i,j,kp))/4.
      pr(i,j,kfn)=0.
      pr(i,j,kfp)=0.
      enddo
      enddo
      enddo
      do j=1,n2m
      do kn=1,n3lm
      kfn=2*kn-1
      kfp=2*kn
      kp=kn+1
      if(kn.eq.n3lm) kp=1
      do i=1,n1m
      qcap(i,j,kfn)=q3(i,j,kn)
      qcap(i,j,kfp)=(q3(i,j,kn)+q3(i,j,kp))/2.
      enddo
      enddo
      enddo
      call q2div
      return                                                            
      end                                                               
c
c
c  *************  subroutine intrpr  ************
c
      subroutine intrpr(n2o)
      include'param.f'
      dimension fo(2*m2),yn(2*m2),fni(2*m2),yu(2*m2)
      dimension y2sol(2*m2)
c
c in this subr evaluate the velocity, read from a continuation
c file, by an cubic spline interpolation in the new grid
c
c
c
      do j=1,n2o
      y2sol(j)=rmo(j)
      enddo
      n2om=n2o-1
      write(6,*)'n2o=',n2o
      do 1 k=1,n3m
      do 1 i=1,n1m
      do 10 j=1,n2om
      fo(j)=q1(i,j,k)
   10 continue
      dpn2=0.
      dp1=0.
c
c
      call spline(y2sol,fo,n2om,dp1,dpn2,yu)
c
      do 11 j=1,n2m
      xx=rm(j)
      call splint(y2sol,fo,yu,n2om,xx,yy)
      fni(j)=yy
   11 continue
      do 12 j=1,n2m
      dq(i,j,k)=fni(j)
      pr(i,j,k)=0.
   12 continue
    1 continue
c
c
      do 2 k=1,n3m
      do 2 i=1,n1m
      do 20 j=1,n2om
      fo(j)=q3(i,j,k)
   20 continue
      dpn2=0.
      dp1=0.
c
c
      call spline(y2sol,fo,n2om,dp1,dpn2,yu)
c
      do 21 j=1,n2m
      xx=rm(j)
      call splint(y2sol,fo,yu,n2om,xx,yy)
      fni(j)=yy
   21 continue
      do 22 j=1,n2m
      qcap(i,j,k)=fni(j)
   22 continue
    2 continue
c
c
c
      call q2div
c
c
      return
      end
c
c    this routines was taken by NUMERICAL RECEPIES
c
      subroutine spline(xin,y,n,rp1,rpn,y2)
      include'param.f'
      dimension xin(2*m2),y(2*m2),y2(2*m2),u(2*m2)
      if (rp1.gt..99e30) then
      y2(1)=0.
      u(1)=0.
      else
      y2(1)=-0.5
      u(1)=(3./(xin(2)-xin(1)))*((y(2)-y(1))
     *     /(xin(2)-xin(1))-rp1)
      endif
      do 11 i=2,n-1
      sig=(xin(i)-xin(i-1))/(xin(i+1)-xin(i-1))
      pnn=sig*y2(i-1)+2.
      y2(i)=(sig-1.)/pnn
      u(i)=(6.*((y(i+1)-y(i))
     *   /(xin(i+1)-xin(i))-(y(i)-y(i-1))
     *   /(xin(i)-xin(i-1)))
     *   /(xin(i+1)-xin(i-1))-sig*u(i-1))/pnn
11    continue
      if (rpn.gt..99e30) then
      qn=0.
      un=0.
      else
      qn=0.5
      un=(3./(xin(n)-xin(n-1)))*(rpn-(y(n)-y(n-1))
     *   /(xin(n)-xin(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do 12 k=n-1,1,-1
      y2(k)=y2(k)*y2(k+1)+u(k)
12    continue
      return
      end
c
c    this routines was taken by NUMERICAL RECEPIES
c
      subroutine splint(xain,ya,y2a,n,x,y)
      dimension xain(1),ya(1),y2a(1)
      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xain(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xain(khi)-xain(klo)
      if (h.eq.0.) pause 'bad xa input.'
      a=(xain(khi)-x)/h
      b=(x-xain(klo))/h
      y=a*ya(klo)+b*ya(khi)+
     *      ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
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
      ru1(i,j,k)=0.                                                     
      ru2(i,j,k)=0.                                                     
      ru3(i,j,k)=0.                                                     
    4 continue                                                          
      return                                                            
      end                                                               
