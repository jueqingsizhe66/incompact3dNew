c************************************************************************
c                                                                       *
c    *****   subro outh   ********************************              *
c                                                                       *
c************************************************************************
      subroutine outh(time,nav,ntime,cflm,nvv,navbu)        
      include 'param.f'
      character*60 filpvp,filvv,filuvt,filvwt,filvtv,filvtd,filvtn
      character*60 fivevo,filhel,filpro,filprt,filvmt,filvcock
      character*60 filpgr
      common/utapr/utap
      character*60 filebu1,filebu2,filebu3,filebu4,filebu5
      dimension skp(4),flp(4),vmp(3),tstp(6),vstp(6),vop(3)
      dimension pvp(4),pco(4),bud(27),hel(3),vvo(6),pgrp(3)
      dimension vistzr(m2),visttr(m2),gpr(m2),drvr2(m2),tvittr(m2)
      dimension dvzvrs(m2),dvtm(m2),vstc(m2),vrvxc(m2),visdif(m2)
      dimension prx(m2),prt(m2),prxa(m2),prta(m2)
      dimension vxc(m2),dstr(3,m2)
      dimension quax(3,m2),quat(3,m2),quar(3,m2)
      common/nbal/nba,ibudg
      character*60 filns1,filts1,filns2,filts2,filns3,filts3 
      n1mh=n1m/2
      n3mh=n3m/2
       nav=nav+1
       navbu=navbu+1
       nvv=nvv+1
      cfo=cfnp+cfo
      dp3mo=dp3mo+dp3ns
      qsouro=qsouro+qsourc
      cfav=cfo/nav
      dp3mav=dp3mo/nav
      cfavp=(abs(dp3mav/2.))
      utal=sqrt(abs(cfnp))
      utam=sqrt(abs(cfav))
      utap=sqrt(abs(cfavp))
      n1vm=n1m
      n2vm=n2m
      n3vm=n3m
c
c   statistics of velocity fluctuations 
c   and relative pdf
c
      call veltur(ntime)
c
c   statistics of vorticity fluctuations 
c   and relative pdf
c
      call vortur(ntime)
cmax
      if(iav.eq.1) then
c         print *,'calcolo taver ', nav
      call taver
      call vetrve
                   endif

c
c
c     write(16,*)' in outh nba,',nba,' nav',nav,'t=',time
      enavp=enavo/nav
      disstp=dissto/nav
      if(nav.ge.1) then
c     write(6,*)' in outh nav=',nav,
c    1   ' number samples  ','   nvv=',nvv,
c    1    'for vertic. transp'
      do jc=2,n2m-1
      jm=jc-1
      jp=jc+1
      enekm=(tursto(1,jm)**2+tursto(2,jm)**2+tursto(3,jm)**2)*0.5/nav**2
      enekc=(tursto(1,jc)**2+tursto(2,jc)**2+tursto(3,jc)**2)*0.5/nav**2
      enekp=(tursto(1,jp)**2+tursto(2,jp)**2+tursto(3,jp)**2)*0.5/nav**2
      visdif(jc)=((enekp-enekc)/(rm(jp)-rm(jc))
     1           -(enekc-enekm)/(rm(jc)-rm(jm)))/(rc(jp)-rc(jc))
      enddo
      jc=1
      jp=jc+1
      enekc=(tursto(1,jc)**2+tursto(2,jc)**2+tursto(3,jc)**2)*0.5/nav**2
      enekp=(tursto(1,jp)**2+tursto(2,jp)**2+tursto(3,jp)**2)*0.5/nav**2
      visdif(jc)=((enekp-enekc)/(rm(jp)-rm(jc))
     1           -enekc/(rm(jc)-rc(jc)))/(rc(jp)-rc(jc))
      jc=n2m
      jm=jc-1
      enekm=(tursto(1,jm)**2+tursto(2,jm)**2+tursto(3,jm)**2)*0.5/nav**2
      enekc=(tursto(1,jc)**2+tursto(2,jc)**2+tursto(3,jc)**2)*0.5/nav**2
      visdif(jc)=-(enekc-enekm)/(rm(jc)-rm(jm))/(rc(jp)-rc(jc))
      do j=2,n2m
      dvzvrs(j)=(tursto(5,j)-tursto(5,j-1))/(rm(j)-rm(j-1))/nav
      vistzr(j)=(vmeo(3,j)-vmeo(3,j-1))/(rm(j)-rm(j-1))/nav/ren
      dvtm(j)=(vmeo(1,j)*rm(j)-vmeo(1,j-1)*rm(j-1))
     1          /rc(j)/(rm(j)-rm(j-1))/nav/ren
      vstc(j)=-(tursto(4,j)*rm(j)+tursto(4,j-1)*rm(j-1))/nav*0.5
      visttr(j)=(vmeo(1,j)/rm(j)-vmeo(1,j-1)/rm(j-1))
     1          *rc(j)/(rm(j)-rm(j-1))/nav/ren
      tvittr(j)=(visttr(j)+tursto(4,j)/nav)*rc(j)**2
      gpr(j)=(pmeo(j)-pmeo(j-1))*rc(j)/(rm(j)-rm(j-1))/nav
      tstpc=(tursto(2,j)/nav)
      tstpm=(tursto(2,j-1)/nav)
      drvr2(j)=(tstpc**2*rm(j)-tstpm**2*rm(j-1))
     1        /(rm(j)-rm(j-1))/nav/ren
      quat(1,j)=-(tursto(4,j)*rm(j)*vmeo(1,j)
     1           +tursto(4,j-1)*rm(j-1)*vmeo(1,j-1))/nav**2*0.5
      quax(1,j)=-(tursto(5,j)*rm(j)*vmeo(3,j)
     1           +tursto(5,j-1)*rm(j-1)*vmeo(3,j-1))/nav**2*0.5
      quat(2,j)=-(tursto(4,j)*rm(j)**2
     1           +tursto(4,j-1)*rm(j-1)**2)/nav*0.5
      quax(2,j)=-(tursto(5,j)*rm(j)
     1           +tursto(5,j-1)*rm(j-1))/nav*0.5
      quax(3,j)=(vmeo(3,j)+vmeo(3,j-1))*0.5/nav
      quat(3,j)=(vmeo(1,j)+vmeo(1,j-1))*0.5/nav
      quar(2,j)=(tursto(2,j)**2-tursto(1,j)**2-tursto(3,j)**2
     1          +tursto(2,j-1)**2-tursto(1,j-1)**2-tursto(3,j-1)**2)
     1          /nav**2*0.25
      enddo
      dvzvrs(1)=tursto(5,1)/(rm(1)-rc(j-1))/nav
      dvzvrs(n2)=-tursto(5,n2-1)/(rc(n2)-rm(n2-1))/nav
      do l=1,3
      quax(l,1)=0.
      quat(l,1)=0.
      quar(l,1)=0.
      quax(l,n2)=0.
      quat(l,n2)=0.
      enddo
      quar(2,n2)=quar(2,n2m)
      tvittr(1)=0.
      tvittr(n2)=0.
      drvr2(1)=0.
      drvr2(n2)=0.
      gpr(1)=0.
      gpr(n2)=0.
      vstc(1)=0.
      vstc(n2)=0.
      dvtm(1)=0.
      dvtm(n2)=-vmeo(1,n2-1)*rm(n2-1)/(rc(n2)-rm(n2m))/nav/ren
      vistzr(1)=0.
      visttr(1)=0.
      vistzr(n2)=-vmeo(3,n2-1)/(rc(n2)-rm(n2m))/nav/ren
      visttr(n2)=-vmeo(1,n2-1)/(rc(n2)-rm(n2m))/nav/ren
c
c  Reynolds stress differentiation
c
      do j=1,n2m
      dstr(1,j)=(quax(2,j+1)-quax(2,j))/(rc(j+1)-rc(j))/rm(j)
      dstr(2,j)=(quat(2,j+1)-quat(2,j))/(rc(j+1)-rc(j))/rm(j)**2
      dstr(3,j)=(quar(2,j+1)-quar(2,j))/(rc(j+1)-rc(j))
     1         +(tursto(2,j)**2-tursto(1,j)**2)/nav**2/(2.*rm(j))
      prx(j)=(quax(1,j+1)-quax(1,j))/(rc(j+1)-rc(j))/rm(j)
      prt(j)=(quat(1,j+1)-quat(1,j))/(rc(j+1)-rc(j))/rm(j)
      prxa(j)=tursto(5,j)*(quax(3,j+1)-quax(3,j))/(rc(j+1)-rc(j))/nav
      prta(j)=-tursto(4,j)*vmeo(1,j)/rm(j)/nav**2
     1        +tursto(4,j)*(quat(3,j+1)-quat(3,j))/(rc(j+1)-rc(j))/nav
      enddo
      itim=nint(time)
      write(pntim,77) itim
   77 format(i4.4)
      if(mod(ntime,npouth).eq.0) then
c
c   the outputs are writtenc
c
      filve='oupivm.'//pntim
      filuvp='oupiuvp.'//pntim
      filuvt='oupiuvt.'//pntim
      filvwt='oupivwt.'//pntim
      filvmt='oupivmt.'//pntim
      filnsp='oupinsp.'//pntim
      filvtb='oupivtb.'//pntim
      filpvp='oupipvp.'//pntim
      filvsp='oupivsp.'//pntim
      filvmp='oupivmp.'//pntim
      filvtn='oupivtn.'//pntim
      filvtv='oupivtv.'//pntim
      filvtd='oupivtd.'//pntim
      filvcock='oupivcock.'//pntim
      filnu='oupinu.'//pntim
      filvo='oupivo.'//pntim
      filvn='oupivn.'//pntim
      filvt='oupivt.'//pntim
      filvv='oupivv.'//pntim
      filns='oupins.'//pntim
      filts='oupits.'//pntim
      filsk='oupisf.'//pntim
      filfl='oupifl.'//pntim
      filpro='oupipr1.'//pntim
      filprt='oupipr2.'//pntim
      filpgr='oupipgr.'//pntim
      filhel='oupihel.'//pntim
      fivevo='oupivev.'//pntim
      open(61,file=filpgr)
      open(91,file=fivevo)
      open(92,file=filpro)
      open(93,file=filhel)
      open(94,file=filprt)
      open(95,file=filvcock)
      open(81,file=filvmp)
      open(82,file=filuvp)
      open(83,file=filnsp)
      open(84,file=filvsp)
      open(85,file=filvtb)
      open(86,file=filpvp)
      open(89,file=filuvt)
      open(80,file=filvmt)
      open(87,file=filvwt)
      open(88,file=filvv)
      open(48,file=filvtv)
      open(47,file=filvtd)
      open(45,file=filvtn)
      open(51,file=filve)
      open(52,file=filvo)
      open(53,file=filvn)
      open(54,file=filvt)
      open(57,file=filsk)
      open(58,file=filfl)
      filns='oupins.'//pntim
      open(55,file=filns)
      filts='oupits.'//pntim
      open(56,file=filts)
                    		endif
      if(mod(ntime,npouth).eq.0) then
      filns1='oupins1.'//pntim
      open(33,file=filns1)
      filts1='oupits1.'//pntim
      open(35,file=filts1)
      filns2='oupins2.'//pntim
      open(37,file=filns2)
      filts2='oupits2.'//pntim
      open(39,file=filts2)
      filns3='oupins3.'//pntim
      open(38,file=filns3)
      filts3='oupits3.'//pntim
      open(36,file=filts3)
                                                    endif
      rewind(42)
      rewind(62)
      do  j=1,n2m
      pmp=pmeo(j)/nav
      do  l=1,4
      skp(l)=skeo(l,j)/nav
      flp(l)=flao(l,j)/nav
      pvp(l)=pvmo(l,j)/nav
      enddo
      do  l=1,3
      vmp(l)=vmeo(l,j)/nav
      vop(l)=voro(l,j)/nav
      enddo
      umpl=vmp(3)/utap
      do  l=1,3
      tstp(l)=(tursto(l,j)/nav)
      vstp(l)=(vorsto(l,j)/nav)
      hel(l)=(vdtomo(l,j)/nav)
      pgrp(l)=(prgro(l,j)/nav)
      enddo
      do  l=1,6
      vvo(l)=(vcromo(l,j)/nav)
      enddo
      do  l=4,6
      tstp(l)=tursto(l,j)/nav
      vstp(l)=vorsto(l,j)/nav
      enddo
      uupl=tstp(1)/utap
      vvpl=tstp(2)/utap
      wwpl=tstp(3)/utap
      dissp=dissjo(j)/nav*2./ren
      drl=rc(j+1)-rc(j)
      akl=(1./(dissp*ren**3))**(1./4.)
      otpl=(vstp(1))/(utap**2*ren)
      orpl=(vstp(2))/(utap**2*ren)
      ozpl=(vstp(3))/(utap**2*ren)
      uvpl=-tstp(5)/utap/utap
      uwpl=-tstp(4)/utap/utap
      vwpl=-tstp(6)/utap/utap
      suu=skp(1)/tstp(1)**3.
      svv=skp(2)/tstp(2)**3.
      sww=skp(3)/tstp(3)**3.
      spp=skp(4)/pvp(4)**(1.5)
      fuu=flp(1)/tstp(1)**4.
      fvv=flp(2)/tstp(2)**4.
      fww=flp(3)/tstp(3)**4.
      fpp=flp(4)/pvp(4)**2.
      yd=1.-rm(j)
      ypl=yd*utap*ren
      uplo=1./.41*alog(ypl)+5.5
      toszrv=(vistzr(j)+vistzr(j+1))*0.5/utap/utap
      tostrv=(visttr(j)+visttr(j+1))*0.5/utap/utap
      toszrt=tstp(5)/utap/utap
      tostrt=tstp(4)/utap/utap
      totszr=(toszrv+toszrt)
      totstr=(tostrv+tostrt)
      dqrp=(drvr2(j)+drvr2(j+1))*0.5
      qt2p=-tstp(1)**2
      rovt=-ros*vmp(1)*rm(j)
      vtqu=-vmp(1)**2
      gprp=(gpr(j)+gpr(j+1))*0.5
      bavr=gprp+vtqu+rovt+qt2p+dqrp
      bavt=(tvittr(j+1)-tvittr(j))/(rc(j+1)-rc(j))
      vtm1=(dvtm(j+1)-dvtm(j))/(rc(j+1)-rc(j))
      vtm2=-(vstc(j+1)-vstc(j))/(rc(j+1)-rc(j))
      vtm3=-tursto(4,j)/nav/rm(j)
      bavtm=vtm1+vtm2+vtm3
      pvtd1=dvto(1,j)/nvv
      pvtd2=dvto(2,j)/nvv
      pvtd3=dvto(3,j)/nvv
      pvtd4=dvto(4,j)/nvv
      pvtd5=dvto(5,j)/nvv
      pvtn1=vtvo(1,j)/nvv
      pvtn2=vtvo(2,j)/nvv
      pvtn3=vtvo(3,j)/nvv
      pvtn4=vtvo(4,j)/nvv
      pvtn5=vtvo(5,j)/nvv
      pvt1=vtvo(1,j)/dvto(1,j)
      pvt2=vtvo(2,j)/dvto(2,j)
      pvt3=vtvo(3,j)/dvto(3,j)
      pvt4=vtvo(4,j)/dvto(4,j)
      pvt5=vtvo(5,j)/dvto(5,j)
      if(mod(ntime,npouth).eq.0) then
      write(81,612)ypl,umpl,uplo
      write(83,612)ypl,uupl,vvpl,wwpl
      write(82,612)yd,uvpl,uwpl,vwpl
      write(84,612)ypl,otpl,orpl,ozpl
      write(85,612)ypl,dissp,akl,drl
      write(86,612)yd,(pvp(l),l=1,4)
      write(87,612)yd,toszrv,toszrt,totszr
      write(89,612)yd,tostrv,tostrt,totstr,bavt
      write(80,612)yd,vtm1,vtm2,vtm3,bavtm
      write(88,612)yd,dqrp,qt2p,rovt,vtqu,gprp,bavr
      write(51,612)yd,(vmp(l),l=1,3),pmp
      write(52,612)yd,(vop(l),l=1,3)
      write(53,612)yd,(vstp(l),l=1,3)
      write(54,612)yd,(vstp(l),l=4,6)
      write(55,612)yd,(tstp(l),l=1,3)
      write(56,612)yd,(tstp(l),l=4,6)
      write(57,612)yd,suu,svv,sww,spp
      write(58,612)yd,fuu,fvv,fww,fpp
      write(48,612)yd,pvt1,pvt2,pvt3,pvt4,pvt5
      write(47,612)yd,pvtd1,pvtd2,pvtd3,pvtd4,pvtd5
      write(45,612)yd,pvtn1,pvtn2,pvtn3,pvtn4,pvtn5
      write(91,612)yd,(vvo(l),l=1,6)
      write(93,612)yd,(hel(l),l=1,3)
      prod11=prx(j)
      prod1a=vmp(3)*dstr(1,j)
      prod1b=vmp(3)*(vvo(2)-vvo(1))
      prod1c=prxa(j)     
      prod22=prt(j)
      prod2c=prta(j)     
      prod2a=vmp(1)*dstr(2,j)
      prod2b=vmp(1)*(vvo(4)-vvo(3))
      prod1=prod1b-prod11
      prod2=prod2b-prod22
      vcomx=(vvo(2)-vvo(1))
      vcomt=(vvo(4)-vvo(3))
      vcomr=(vvo(6)-vvo(5))
      write(92,612)yd,prod11,prod1b,prod1,prod1a,prod1c
      write(94,612)yd,prod22,prod2b,prod2,prod2a,prod2c
      write(95,612)yd,dstr(1,j),vcomx,dstr(2,j),vcomt,dstr(3,j),vcomr
      write(61,612)yd,prod1b,prod11,pgrp(2),visdif(j),prod2b,prod22
                                endif
  612 format(1x,e12.4,2x,9(1x,e12.5))
      if(mod(ntime,npouth).eq.0) then
      write(33,612)yd,tstp(1)
      write(35,612)yd,tstp(4)
      write(37,612)yd,tstp(2)
      write(39,612)yd,tstp(5)
      write(38,612)yd,tstp(3)
      write(36,612)yd,tstp(6)
                             endif
      enddo
      close(61)
      close(91)
      close(92)
      close(93)
      close(94)
      close(95)
      close(45)
      close(47)
      close(48)
      close(71)
      close(72)
      close(73)
      close(81)
      close(82)
      close(83)
      close(84)
      close(85)
      close(86)
      close(87)
      close(88)
      close(89)
      close(80)
      close(51)
      close(52)
      close(53)
      close(54)
      close(57)
      close(58)
      close(55)
      close(56)
      if(mod(ntime,nprde).eq.0.or.mod(ntime,npouth).eq.0) then
      close(33)
      close(35)
      close(37)
      close(39)
      close(38)
      close(36)
                                                          endif
      if(mod(ntime,npouth).eq.0) then
      if(icorspe.eq.1) then
c     write(6,*)' in outh nba',nba,
c    1          ' number samples to average spectra'
      namfile='kolmog.out'
      open(58,file=namfile)
      do k=1,k1max
      ake=k
      enko=.1e-02*ake**(-5./3.)
      write(58,133)ake,enko
      enddo
      npq=1
c
c   the one dimensional azimuthal spectra are wriiten
c
      do j=n2m,1,-1
      if(j.eq.jprq(npq)) then
      write(jsp,177) j
  177 format(i2.2)
      namfile='envot'//jsp//'.'//pntim
      open(55,file=namfile)
      namfile='enpt'//jsp//'.'//pntim
      open(57,file=namfile)
      namfile='encot'//jsp//'.'//pntim
      open(56,file=namfile)
      do k=2,k1max
      ake=(k-1)
      penep=enepto(j,k)/nba
      penv1=env11t(j,k)/nba
      penv2=env22t(j,k)/nba
      penv3=env33t(j,k)/nba
      pene1=ene1to(j,k)/nba
      pene2=ene2to(j,k)/nba
      pene3=ene3to(j,k)/nba
      write(55,133)ake,penv1,penv2,penv3
      write(57,133)ake,penep
      write(56,133)ake,pene1,pene2,pene3
  133 format(9e13.5)
      enddo
      close(58)
      close(57)
      close(56)
      close(55)
      namfile='entze'//jsp//'.'//pntim
      open(56,file=namfile)
      k=1
      ake=(k-1)
      penep=enepto(j,k)/nba
      penv1=env11t(j,k)/nba
      penv2=env22t(j,k)/nba
      penv3=env33t(j,k)/nba
      pene1=ene1to(j,k)/nba
      pene2=ene2to(j,k)/nba
      pene3=ene3to(j,k)/nba
      write(56,133)ake,pene1,pene2,pene3,penv1,penv2,penv3,penep
      npq=npq+1
                        endif
      enddo
      npq=1
c
c   the one dimensional axial spectra are wriiten
c
      do j=n2m,1,-1
      if(j.eq.jprq(npq)) then
      write(jsp,177) j
      namfile='envoz'//jsp//'.'//pntim
      open(55,file=namfile)
      namfile='enpz'//jsp//'.'//pntim
      open(57,file=namfile)
      namfile='encoz'//jsp//'.'//pntim
      open(56,file=namfile)
      pi=2.*asin(1.)
      do k=2,k3max
      ake=(k-1)*2.*pi/alx3d
      penet=(ene1zo(j,k)+ene2zo(j,k)+ene3zo(j,k))/nba
      penep=enepzo(j,k)/nba
      penv1=env11z(j,k)/nba
      penv2=env22z(j,k)/nba
      penv3=env33z(j,k)/nba
      pene1=ene1zo(j,k)/nba
      pene2=ene2zo(j,k)/nba
      pene3=ene3zo(j,k)/nba
      write(55,133)ake,penv1,penv2,penv3
      write(57,133)ake,penep
      write(56,133)ake,pene1,pene2,pene3
      enddo
      close(58)
      close(57)
      close(56)
      close(55)
      namfile='enzze'//jsp//'.'//pntim
      open(56,file=namfile)
      k=1
      ake=(k-1)*2.*pi/alx3d
      penep=enepzo(j,k)/nba
      penv1=env11z(j,k)/nba
      penv2=env22z(j,k)/nba
      penv3=env33z(j,k)/nba
      pene1=ene1zo(j,k)/nba
      pene2=ene2zo(j,k)/nba
      pene3=ene3zo(j,k)/nba
      write(56,133)ake,pene1,pene2,pene3,penv1,penv2,penv3,penep
      npq=npq+1
                        endif
      enddo
      namfile='ten1.'//pntim
      open(55,file=namfile)
      namfile='ten2.'//pntim
      open(56,file=namfile)
      namfile='ten3.'//pntim
      open(57,file=namfile)
      do j=1,n2m
      yd=1.-rm(j) 
      pet1=e1to(j)/nba
      pet2=e2to(j)/nba
      pet3=e3to(j)/nba
      pest1=e1stto(j)/nba
      pest2=e2stto(j)/nba
      pest3=e3stto(j)/nba
      pesz1=e1stzo(j)/nba
      pesz2=e2stzo(j)/nba
      pesz3=e3stzo(j)/nba
      write(55,133)yd,pet1,pest1,pesz1
      write(57,133)yd,pet3,pest3,pesz3
      write(56,133)yd,pet2,pest2,pesz2
	 enddo
      close(57)
      close(56)
      close(55)
      npq=1
c
c    velocity and vorticity azimuthal correlations are written
c
      do j=n2m,1,-1
      if(j.eq.jprq(npq)) then
      jrpl=j
      write(jrp,157) jrpl
  157 format(i3.3)
      namfile='viith'//jrp//'.'//pntim
      open(83,file=namfile)
      namfile='riith'//jrp//'.'//pntim
      open(81,file=namfile)
      do i=1,n1mh
      rtplus=-rm(j)*thetac(i)*utap*ren
      pr11=r11the(j,i)/nba
      pr22=r22the(j,i)/nba
      pr33=r33the(j,i)/nba
      pppc=pcothe(j,i)/nba
      write(81,133)rtplus,pr11,pr22,pr33,pppc
      pv11=v11the(j,i)/nba
      pv22=v22the(j,i)/nba
      pv33=v33the(j,i)/nba
      write(83,133)rtplus,pv11,pv22,pv33
      enddo
      close(83)
      close(81)
      npq=npq+1
                        endif
      enddo
c
c    velocity and vorticity axial correlations are written
c
      npq=1
      do j=n2m,1,-1
      if(j.eq.jprq(npq)) then
      jrpl=j
      write(jrp,157) jrpl
      namfile='viiax'//jrp//'.'//pntim
      open(83,file=namfile)
      namfile='riiax'//jrp//'.'//pntim
      open(81,file=namfile)
      do k=1,n3mh
      zzmm=(zz(k)+zz(k+1))*0.5
c     rtplus=zzmm*utap*ren
      rtplus=zzmm
      pr11=r11axi(j,k)/nba
      pr22=r22axi(j,k)/nba
      pr33=r33axi(j,k)/nba
      pppc=pcoaxi(j,k)/nba
      write(81,133)rtplus,pr11,pr22,pr33,pppc
      pv11=v11axi(j,k)/nba
      pv22=v22axi(j,k)/nba
      pv33=v33axi(j,k)/nba
      write(83,133)rtplus,pv11,pv22,pv33
      enddo
      close(83)
      close(81)
      npq=npq+1
                        endif
      enddo
		endif
                               endif
  983 format(3x,e11.4,i5,9(1x,e12.5))
  783 format(1x,11(e11.4))
  785 format(3x,7(1x,e12.5),3x,i5)
      retau=utam*ren
      retap=utap*ren
      uitau=vmeo(3,1)/nav/utam
      uitap=vmeo(3,1)/nav/utap
      reucl=ren*vmeo(3,1)/nav
      write(59,785)time,uitap,retap,utap,cfav,dp3mav,reucl
     1             ,nav
			endif
      vstzrn=-vmeo(3,n2-1)/(rc(n2)-rm(n2m))/nav/ren
      vsttrn=-vmeo(1,n2-1)/(rc(n2)-rm(n2m))/nav/ren
      tauwzr=sqrt(abs(vstzrn))
      tauwtr=sqrt(abs(vsttrn))
      cfzr=tauwzr**2*2./vit(3)**2
      cftr=tauwtr**2*2./vit(3)**2
      utat=sqrt(abs(dp3ns/2.))
      ttat=qsourc/(2.*utat)
      write(50,783)time,vit(5),(vit(l),l=1,3,2),enavp,utat
     1 ,dp3ns,cflm
      write(32,783)time,(vit(l),l=1,5),enavp,utat
     1 ,dp3ns
      write(6,783)time,cfzr,cftr,uitap,enavp,utat,utap
     1 ,dissma,retap
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
      do j=1,n2mv
      do i=1,n1mv
      if(l.eq.1) vca=q1(i,j,k)
      if(l.eq.2) vca=q2(i,j,k)
      if(l.eq.3) vca=q3(i,j,k)
      if(abs(vca).ge.vmax(l)) then
      vmax(l)=abs(vca)
                              endif
      enddo
      enddo
      enddo
      vmaxo(l)=0.
      do k=1,n3mv
      do j=2,n2mv
      do i=1,n1mv
      if(l.eq.2) vca=q2(i,j,k)/rc(j)
      if(abs(vca).ge.vmaxo(l)) then
      vmaxo(l)=abs(vca)
                              endif
      enddo
      enddo
      enddo
      do k=1,n3mv
      do j=1,n2mv
      do i=1,n1mv
      if(l.eq.1) vca=q1(i,j,k)/rm(j)
      if(l.eq.3) vca=q3(i,j,k)
      if(abs(vca).ge.vmaxo(l)) then
      vmaxo(l)=abs(vca)
                              endif
      enddo
      enddo
      enddo
  311 continue
      return
      end

