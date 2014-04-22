c******************************* subrout wstre **********************
c
      subroutine wstre
c
      include 'param.f'
      common/wallst/cflw,cfuw
c
c     lower and upper walls shear
c
      dql=0.
      dqu=0.
      dyl=(y2s(1)-y(1))
      dyu=-(y2s(n2m)-y(n2))
      do 450 kc=1,n3m
      do 450 ic=1,n1m
      dql=(q3(ic,1,kc)+q3(ic,1,kpv(kc)))*.5/dyl+dql
      dqu=-(q3(ic,n2m,kc)+q3(ic,n2m,kpv(kc)))*.5/dyu+dqu
  450 continue
      cfl1=dql/(re*n1m*n3m)
      cfu1=dqu/(re*n1m*n3m)
      cfuw=cfu1
      cflw=cfl1
      return
      end
c************************************************************************
c                                                                       *
c    *****   subro outh   ********************************              *
c                                                                       *
c************************************************************************
      subroutine outh(time,nav,ntime,cflm,nvv,navbu)        
      include 'param.f'
      character*60 filpvp,filvv,filuvt,filvwt,filvtv,filvtd,filvtn
      character*60 fivevo,filhel,filpro,filqst,filvmt,filvcock
      character*60 filpgr,filqno
      common/utapr/utap
      character*60 filebu1,filebu2,filebu3,filebu4,filebu5
      dimension skp(4),flp(4),vmp(3),tstp(6),vstp(6),vop(3)
      dimension sqpq(12),skpq(12),flpq(12)
      dimension pvp(4),pco(4),bud(27),hel(3),vvo(6),pgrp(3)
      dimension vistzr(m2),visttr(m2),gpr(m2),drvr2(m2),tvittr(m2)
      dimension dvzvrs(m2),dvtm(m2),vstc(m2),vrvxc(m2),visdif(m2)
      dimension prx(m2),prt(m2),prxa(m2),prta(m2)
      dimension vxc(m2),dstr(3,m2)
      dimension quax(3,m2),quat(3,m2),quar(3,m2)
      common/nbal/nba,ibudg
      common/walfr/cfuo,cflo
      common/wallst/cflw,cfuw
      character*60 filns1,filts1,filns2,filts2,filns3,filts3 
      character*3 njpse
      n1mh=n1m/2
      n3mh=n3m/2
       nav=nav+1
       navbu=navbu+1
       nvv=nvv+1
      call wstre
      cfo=cfnp+cfo
      cfuo=cfuw+cfuo
      cflo=cflw+cflo
      dp3mo=dp3mo+dp3ns
      cfav=cfo/nav
      cfuav=cfuo/nav
      cflav=cflo/nav
      dp3mav=dp3mo/nav
      cfavp=(abs(dp3mav))
      utaul=sqrt(abs(cflav))
      utauu=sqrt(abs(cfuav))
      utap=sqrt(abs(cfavp))
      n1vm=n1m
      n2vm=n2m
      n3vm=n3m
c
c   velocity fluctuations statistics
c
      call veltur(ntime)
c
c   evaluation of the averages in planes of certain quantities
c
      call quasfa(ntime)
c
c   evaluation of the statistics of vorticity and related quantities
c
      call vortur(ntime)
      if(iav.eq.1) then
c         calcolo taver 
      call taver
      call vetrve
                   endif

c
c
      enavp=enavo/nav
      disstp=dissto/nav
      if(nav.ge.1) then
c     write(6,*)' in outh nav=',nav,
c    1   ' number samples  ','   nvv=',nvv,
c    1    'for vertic. transp'
      do j=2,n2m
      vistzr(j)=(vmeo(3,j)-vmeo(3,j-1))/(y2s(j)-y2s(j-1))/nav/ren
      quat(1,j)=-(tursto(4,j)*vmeo(1,j)
     1           +tursto(4,j-1)*vmeo(1,j-1))/nav**2*0.5
      quax(1,j)=-(tursto(5,j)*vmeo(3,j)
     1           +tursto(5,j-1)*vmeo(3,j-1))/nav**2*0.5
      quat(2,j)=-(tursto(4,j)
     1           +tursto(4,j-1))/nav*0.5
      quax(2,j)=-(tursto(5,j)
     1           +tursto(5,j-1))/nav*0.5
      quax(3,j)=(vmeo(3,j)+vmeo(3,j-1))*0.5/nav
      quat(3,j)=(vmeo(1,j)+vmeo(1,j-1))*0.5/nav
      quar(2,j)=(tursto(2,j)**2-tursto(1,j)**2-tursto(3,j)**2
     1          +tursto(2,j-1)**2-tursto(1,j-1)**2-tursto(3,j-1)**2)
     1          /nav**2*0.25
      enddo
      do l=1,3
      quax(l,1)=0.
      quat(l,1)=0.
      quar(l,1)=0.
      quax(l,n2)=0.
      quat(l,n2)=0.
      enddo
      quar(2,n2)=quar(2,n2m)
      vistzr(1)=0.
      vistzr(1)=vmeo(3,2)/(y2s(2)-y(1))/nav/ren
      vistzr(n2)=-vmeo(3,n2-1)/(y(n2)-y2s(n2m))/nav/ren
c
c  Reynolds stress differentiation
c
      do j=1,n2m
      dstr(1,j)=(quax(2,j+1)-quax(2,j))/(y(j+1)-y(j))
      dstr(2,j)=(quat(2,j+1)-quat(2,j))/(y(j+1)-y(j))
      dstr(3,j)=(quar(2,j+1)-quar(2,j))/(y(j+1)-y(j))
     1         +(tursto(2,j)**2-tursto(1,j)**2)/nav**2/(2.)
      prx(j)=(quax(1,j+1)-quax(1,j))/(y(j+1)-y(j))
      prt(j)=(quat(1,j+1)-quat(1,j))/(y(j+1)-y(j))
      prxa(j)=tursto(5,j)*(quax(3,j+1)-quax(3,j))/(y(j+1)-y(j))/nav
      prta(j)=-tursto(4,j)*vmeo(1,j)/nav**2
     1        +tursto(4,j)*(quat(3,j+1)-quat(3,j))/(y(j+1)-y(j))/nav
      enddo
      itim=nint(time)
      write(pntim,77) itim
   77 format(i4.4)
      if(mod(ntime,npouth).eq.0) then
      filve='oupivm.'//pntim
      filuvp='oupiuvp.'//pntim
      filuvt='oupiuvt.'//pntim
      filvwt='oupivwt.'//pntim
      filvmt='oupiens.'//pntim
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
      filpro='oupipr1.'//pntim
      filqst='oupiqst.'//pntim
      filqno='oupiqno.'//pntim
      filpgr='oupipgr.'//pntim
      filhel='oupihel.'//pntim
      fivevo='oupivev.'//pntim
      filsk='oupisqvo.'//pntim
      open(81,file=filsk)
      filsk='oupisqst.'//pntim
      open(82,file=filsk)
      filsk='oupisqee.'//pntim
      open(83,file=filsk)
      filsk='oupisfvo.'//pntim
      open(61,file=filsk)
      filsk='oupisfst.'//pntim
      open(62,file=filsk)
      filsk='oupisfee.'//pntim
      open(63,file=filsk)
      filsk='oupiflvo.'//pntim
      open(64,file=filsk)
      filsk='oupiflst.'//pntim
      open(65,file=filsk)
      filsk='oupiflee.'//pntim
      open(66,file=filsk)
      open(91,file=fivevo)
      open(92,file=filpro)
      open(93,file=filhel)
      open(94,file=filqst)
      open(96,file=filqno)
      open(95,file=filvcock)
      open(115,file=filvmp)
      open(116,file=filuvp)
      open(117,file=filnsp)
      open(118,file=filvsp)
      open(119,file=filvtb)
      open(129,file=filvmt)
      open(120,file=filpvp)
      open(114,file=filvwt)
      open(48,file=filvtv)
      open(47,file=filvtd)
      open(45,file=filvtn)
      open(51,file=filve)
      open(52,file=filvo)
      open(53,file=filvn)
      open(54,file=filvt)
      filsk='oupisf.'//pntim
      filfl='oupifl.'//pntim
      open(57,file=filsk)
      open(58,file=filfl)
      filns='oupins.'//pntim
      open(55,file=filns)
      filts='oupits.'//pntim
      open(56,file=filts)
                    		endif
      if(mod(ntime,nprde).eq.0) then
      filns1='oupins1.'//pntim
      open(121,file=filns1)
      filts1='oupits1.'//pntim
      open(111,file=filts1)
      filns2='oupins2.'//pntim
      open(122,file=filns2)
      filts2='oupits2.'//pntim
      open(112,file=filts2)
      filns3='oupins3.'//pntim
      open(123,file=filns3)
      filts3='oupits3.'//pntim
      open(113,file=filts3)
                                 endif
      if(mod(ntime,npouth).eq.0) then
      do  j=1,n2m
      if(y(j).le.0.) then
      yd=y2s(j)-y(1)
                     else
      yd=y(n2)-y2s(j)
                     endif
      vmp(1)=vpmeo(1,j)/nav
      vmp(3)=vpmeo(3,j)/nav
      vop(2)=vpome(2,j)/nav
      write(94,612)yd,vmp(1),vmp(3),vop(2)
      enddo
      do  j=1,n2
      if(y(j).le.0.) then
      yd=y(j)-y(1)
                     else
      yd=y(n2)-y(j)
                     endif
      vmp(2)=vpmeo(2,j)/nav
      vop(1)=vpome(1,j)/nav
      vop(3)=vpome(3,j)/nav
      write(96,612)yd,vmp(2),vop(1),vop(3)
      enddo
                                 endif
      do  j=1,n2m
      pmp=pmeo(j)/nav
      do  l=1,12
      skpq(l)=skequo(l,j)/nav
      flpq(l)=flaquo(l,j)/nav
      sqpq(l)=sququo(l,j)/nav
      enddo
      do  l=1,4
      skp(l)=skeo(l,j)/nav
      flp(l)=flao(l,j)/nav
      pvp(l)=pvmo(l,j)/nav
      enddo
      do  l=1,3
      vmp(l)=vmeo(l,j)/nav
      vop(l)=voro(l,j)/nav
      enddo
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
      umpl=vmp(3)/utap
      uupl=tstp(1)/utap
      vvpl=tstp(2)/utap
      wwpl=tstp(3)/utap
      dissp=dissjo(j)/nav*2./ren
      sijsij=dissjo(j)/nav
      essp=ensyo(j)/nav
      essprp=enspo(j)/nav
      drl=y(j+1)-y(j)
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
      if(j.le.n2m/2) then
      yd=y2s(j)-y(1)
      ypl=yd*utap*ren
                     else
      yd=-y2s(j)+y(n2)
      ypl=yd*utap*ren
                     endif
      uplo=1./.41*alog(ypl)+5.5
      toszrv=(vistzr(j)+vistzr(j+1))*0.5/utap/utap
      toszrt=tstp(5)/utap/utap
      totszr=(toszrv+toszrt)
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
      write(115,612)ypl,umpl,uplo
      write(117,612)ypl,uupl,vvpl,wwpl
      write(116,612)yd,uvpl,uwpl,vwpl
      write(118,612)ypl,otpl,orpl,ozpl
      write(119,612)ypl,dissp,akl,drl
      write(129,612)yd,sijsij,essp,essprp
      write(120,612)yd,(pvp(l),l=1,4)
      write(114,612)yd,toszrv,toszrt,totszr
      write(81,612)yd,(sqpq(l),l=1,3)
      write(82,612)yd,(sqpq(l),l=4,9)
      write(83,612)yd,(sqpq(l),l=10,12)
      write(61,612)yd,(skpq(l),l=1,3)
      write(62,612)yd,(skpq(l),l=4,9)
      write(63,612)yd,(skpq(l),l=10,12)
      write(64,612)yd,(flpq(l),l=1,3)
      write(65,612)yd,(flpq(l),l=4,9)
      write(66,612)yd,(flpq(l),l=10,12)
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
      prod1=prod1b-prod11
      vcomx=(vvo(2)-vvo(1))
      vcomt=(vvo(4)-vvo(3))
      vcomr=(vvo(6)-vvo(5))
      write(92,612)yd,prod11,prod1b,prod1,prod1a,prod1c
      write(95,612)yd,dstr(1,j),vcomx,dstr(2,j),vcomt,dstr(3,j),vcomr
                                endif
  612 format(1x,e12.4,2x,9(1x,e12.5))
      if(mod(ntime,nprde).eq.0) then
      write(121,612)yd,tstp(1)
      write(111,612)yd,tstp(4)
      write(122,612)yd,tstp(2)
      write(112,612)yd,tstp(5)
      write(123,612)yd,tstp(3)
      write(113,612)yd,tstp(6)
                             endif
      enddo
      close(81)
      close(82)
      close(83)
      close(91)
      close(92)
      close(93)
      close(94)
      close(95)
      close(45)
      close(47)
      close(48)
      close(115)
      close(116)
      close(117)
      close(118)
      close(119)
      close(129)
      close(120)
      close(114)
      close(64)
      close(65)
      close(66)
      close(61)
      close(62)
      close(63)
      close(51)
      close(52)
      close(53)
      close(54)
      close(57)
      close(58)
      close(55)
      close(56)
      if(mod(ntime,nprde).eq.0) then
      close(121)
      close(111)
      close(122)
      close(112)
      close(123)
      close(113)
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
      write(6,*)'k1max,k3max  ',k1max,k3max
      npq=1
      do j=n2m,1,-1
c
c   write spectra and cospectra in k1 wave number
c
      if(j.eq.jprq(npq)) then
      write(njpse,177) j
  177 format(i3.3)
      namfile='speprk1'//njpse//'.'//pntim
      open(69,file=namfile)
      namfile='spevvk1'//njpse//'.'//pntim
      open(55,file=namfile)
c     write(6,*)j,namfile
      namfile='speook1'//njpse//'.'//pntim
      open(56,file=namfile)
      namfile='cospvvk1'//njpse//'.'//pntim
      open(57,file=namfile)
      namfile='cospook1'//njpse//'.'//pntim
      open(54,file=namfile)
      namfile='cospvo1k1'//njpse//'.'//pntim
      open(51,file=namfile)
      namfile='cospvo2k1'//njpse//'.'//pntim
      open(52,file=namfile)
      namfile='cospvo3k1'//njpse//'.'//pntim
      open(53,file=namfile)
      do k=2,k1max
      ake=(k-1)*2./alx1d
      psoo11=env11t(j,k)/nba
      psoo22=env22t(j,k)/nba
      psoo33=env33t(j,k)/nba
      psoo12=env12t(j,k)/nba
      psoo13=env13t(j,k)/nba
      psoo23=env23t(j,k)/nba
      pspre1=enepto(j,k)/nba
      psvv11=ene1to(j,k)/nba
      psvv22=ene2to(j,k)/nba
      psvv33=ene3to(j,k)/nba
      psvv12=en12to(j,k)/nba
      psvv13=en13to(j,k)/nba
      psvv23=en23to(j,k)/nba
      psvo11=evo11t(j,k)/nba
      psvo12=evo12t(j,k)/nba
      psvo13=evo13t(j,k)/nba
      psvo21=evo21t(j,k)/nba
      psvo22=evo22t(j,k)/nba
      psvo23=evo23t(j,k)/nba
      psvo31=evo31t(j,k)/nba
      psvo32=evo32t(j,k)/nba
      psvo33=evo33t(j,k)/nba
      write(69,133)ake,pspre1
      write(55,133)ake,psvv11,psvv22,psvv33
      write(56,133)ake,psoo11,psoo22,psoo33
      write(57,133)ake,psvv12,psvv13,psvv23
      write(54,133)ake,psoo12,psoo13,psoo23
      write(51,133)ake,psvo11,psvo12,psvo13
      write(52,133)ake,psvo21,psvo22,psvo23
      write(53,133)ake,psvo31,psvo32,psvo33
  133 format(9e13.5)
      enddo
      close(69)
      close(58)
      close(56)
      close(55)
      close(57)
      close(54)
      close(53)
      close(52)
      close(51)
      namfile='enk1ze'//njpse//'.'//pntim
      open(56,file=namfile)
      k=1
      ake=(k-1)
      penv1=env11t(j,k)/nba
      penv2=env22t(j,k)/nba
      penv3=env33t(j,k)/nba
      pene1=ene1to(j,k)/nba
      pene2=ene2to(j,k)/nba
      pene3=ene3to(j,k)/nba
      write(56,133)ake,pene1,pene2,pene3,penv1,penv2,penv3
      close(56)
      npq=npq+1
                        endif
      enddo
      npq=1
      do j=n2m,1,-1
      if(j.eq.jprq(npq)) then
      write(njpse,177) j
      namfile='speprk3'//njpse//'.'//pntim
      open(69,file=namfile)
      namfile='spevvk3'//njpse//'.'//pntim
      open(55,file=namfile)
      namfile='speook3'//njpse//'.'//pntim
      open(56,file=namfile)
      namfile='cospvvk3'//njpse//'.'//pntim
      open(57,file=namfile)
      namfile='cospook3'//njpse//'.'//pntim
      open(54,file=namfile)
      namfile='cospvo1k3'//njpse//'.'//pntim
      open(51,file=namfile)
      namfile='cospvo2k3'//njpse//'.'//pntim
      open(52,file=namfile)
      namfile='cospvo3k3'//njpse//'.'//pntim
      open(53,file=namfile)
      do k=2,k3max
      ake=(k-1)*2./alx3d
      psoo11=env11z(j,k)/nba
      psoo22=env22z(j,k)/nba
      psoo33=env33z(j,k)/nba
      psoo12=env12z(j,k)/nba
      psoo13=env13z(j,k)/nba
      psoo23=env23z(j,k)/nba
      pspre3=enepzo(j,k)/nba
      psvv11=ene1zo(j,k)/nba
      psvv22=ene2zo(j,k)/nba
      psvv33=ene3zo(j,k)/nba
      psvv12=en12zo(j,k)/nba
      psvv13=en13zo(j,k)/nba
      psvv23=en23zo(j,k)/nba
      psvo11=evo11z(j,k)/nba
      psvo12=evo12z(j,k)/nba
      psvo13=evo13z(j,k)/nba
      psvo21=evo21z(j,k)/nba
      psvo22=evo22z(j,k)/nba
      psvo23=evo23z(j,k)/nba
      psvo31=evo31z(j,k)/nba
      psvo32=evo32z(j,k)/nba
      psvo33=evo33z(j,k)/nba
      write(69,133)ake,pspre3
      write(55,133)ake,psvv11,psvv22,psvv33
      write(56,133)ake,psoo11,psoo22,psoo33
      write(57,133)ake,psvv12,psvv13,psvv23
      write(54,133)ake,psoo12,psoo13,psoo23
      write(51,133)ake,psvo11,psvo12,psvo13
      write(52,133)ake,psvo21,psvo22,psvo23
      write(53,133)ake,psvo31,psvo32,psvo33
      enddo
      close(69)
      close(58)
      close(56)
      close(55)
      close(57)
      close(54)
      close(53)
      close(52)
      close(51)
      namfile='enk3ze'//njpse//'.'//pntim
      open(56,file=namfile)
      k=1
      ake=(k-1)*2.*pi/alx3d
      penv1=env11z(j,k)/nba
      penv2=env22z(j,k)/nba
      penv3=env33z(j,k)/nba
      pene1=ene1zo(j,k)/nba
      pene2=ene2zo(j,k)/nba
      pene3=ene3zo(j,k)/nba
      write(56,133)ake,pene1,pene2,pene3,penv1,penv2,penv3
      close(56)
      npq=npq+1
                        endif
      enddo
      npq=1
      do j=n2m,1,-1
c
c   write correlation in x1 
c
      if(j.eq.jprq(npq)) then
      write(njpse,177) j
      namfile='corrprx1'//njpse//'.'//pntim
      open(69,file=namfile)
      namfile='corrvvx1'//njpse//'.'//pntim
      open(55,file=namfile)
      namfile='corroox1'//njpse//'.'//pntim
      open(56,file=namfile)
      namfile='crcovvx1'//njpse//'.'//pntim
      open(57,file=namfile)
      namfile='crcooox1'//njpse//'.'//pntim
      open(54,file=namfile)
      namfile='crcovo1x1'//njpse//'.'//pntim
      open(51,file=namfile)
      namfile='crcovo2x1'//njpse//'.'//pntim
      open(52,file=namfile)
      namfile='crcovo3x1'//njpse//'.'//pntim
      open(53,file=namfile)
      namfile='corr0x1'//njpse//'.'//pntim
      open(58,file=namfile)
      pspre1=pcox1(j,1)
      psvv11=r11x1(j,1)
      psvv22=r22x1(j,1)
      psvv33=r33x1(j,1)
      psoo11=v11x1(j,1)
      psoo22=v22x1(j,1)
      psoo33=v33x1(j,1)
      psoo12=v12x1(j,1)
      psoo13=v31x1(j,1)
      psoo23=v23x1(j,1)
      psvv12=r12x1(j,1)
      psvv13=r31x1(j,1)
      psvv23=r23x1(j,1)
      psvo11=vo11x1(j,1)
      psvo12=vo12x1(j,1)
      psvo13=vo13x1(j,1)
      psvo21=vo21x1(j,1)
      psvo22=vo22x1(j,1)
      psvo23=vo23x1(j,1)
      psvo31=vo31x1(j,1)
      psvo32=vo32x1(j,1)
      psvo33=vo33x1(j,1)
      dpl=0.
      write(69,133)dpl,pspre1
      write(55,133)dpl,psvv11,psvv22,psvv33
      write(56,133)dpl,psoo11,psoo22,psoo33
      write(57,133)dpl,psvv12,psvv13,psvv23
      write(54,133)dpl,psoo12,psoo13,psoo23
      write(51,133)dpl,psvo11,psvo12,psvo13
      write(52,133)dpl,psvo21,psvo22,psvo23
      write(53,133)dpl,psvo31,psvo32,psvo33
      do k=1,k1max
      dpl=(k-1)/dx1*utap*ren
      pspre1=pcox1(j,k)/nba
      psvv11=r11x1(j,k)/nba
      psvv22=r22x1(j,k)/nba
      psvv33=r33x1(j,k)/nba
      psoo11=v11x1(j,k)/nba
      psoo22=v22x1(j,k)/nba
      psoo33=v33x1(j,k)/nba
      psoo12=v12x1(j,k)/nba
      psoo13=v31x1(j,k)/nba
      psoo23=v23x1(j,k)/nba
      psvv12=r12x1(j,k)/nba
      psvv13=r31x1(j,k)/nba
      psvv23=r23x1(j,k)/nba
      psvo11=vo11x1(j,k)/nba
      psvo12=vo12x1(j,k)/nba
      psvo13=vo13x1(j,k)/nba
      psvo21=vo21x1(j,k)/nba
      psvo22=vo22x1(j,k)/nba
      psvo23=vo23x1(j,k)/nba
      psvo31=vo31x1(j,k)/nba
      psvo32=vo32x1(j,k)/nba
      psvo33=vo33x1(j,k)/nba
      write(69,133)dpl,pspre1
      write(55,133)dpl,psvv11,psvv22,psvv33
      write(56,133)dpl,psoo11,psoo22,psoo33
      write(57,133)dpl,psvv12,psvv13,psvv23
      write(54,133)dpl,psoo12,psoo13,psoo23
      write(51,133)dpl,psvo11,psvo12,psvo13
      write(52,133)dpl,psvo21,psvo22,psvo23
      write(53,133)dpl,psvo31,psvo32,psvo33
      enddo
      close(69)
      close(58)
      close(56)
      close(55)
      close(57)
      close(54)
      close(53)
      close(52)
      close(51)
      npq=npq+1
                        endif
      enddo
      npq=1
      do j=n2m,1,-1
c
c   write correlation in x3 
c
      if(j.eq.jprq(npq)) then
      write(njpse,177) j
      namfile='corrprx3'//njpse//'.'//pntim
      open(69,file=namfile)
      namfile='corrvvx3'//njpse//'.'//pntim
      open(55,file=namfile)
      namfile='corroox3'//njpse//'.'//pntim
      open(56,file=namfile)
      namfile='crcovvx3'//njpse//'.'//pntim
      open(57,file=namfile)
      namfile='crcooox3'//njpse//'.'//pntim
      open(54,file=namfile)
      namfile='crcovo1x3'//njpse//'.'//pntim
      open(51,file=namfile)
      namfile='crcovo2x3'//njpse//'.'//pntim
      open(52,file=namfile)
      namfile='crcovo3x3'//njpse//'.'//pntim
      open(53,file=namfile)
      namfile='corr0x3'//njpse//'.'//pntim
      open(58,file=namfile)
      pspre3=pcox3(j,1)
      psvv11=r11x3(j,1)
      psvv22=r22x3(j,1)
      psvv33=r33x3(j,1)
      psoo11=v11x3(j,1)
      psoo22=v22x3(j,1)
      psoo33=v33x3(j,1)
      psoo12=v12x3(j,1)
      psoo13=v31x3(j,1)
      psoo23=v23x3(j,1)
      psvv12=r12x3(j,1)
      psvv13=r31x3(j,1)
      psvv23=r23x3(j,1)
      psvo11=vo11x3(j,1)
      psvo12=vo12x3(j,1)
      psvo13=vo13x3(j,1)
      psvo21=vo21x3(j,1)
      psvo22=vo22x3(j,1)
      psvo23=vo23x3(j,1)
      psvo31=vo31x3(j,1)
      psvo32=vo32x3(j,1)
      psvo33=vo33x3(j,1)
      dpl=0.
      write(69,133)dpl,pspre3
      write(55,133)dpl,psvv11,psvv22,psvv33
      write(56,133)dpl,psoo11,psoo22,psoo33
      write(57,133)dpl,psvv12,psvv13,psvv23
      write(54,133)dpl,psoo12,psoo13,psoo23
      write(51,133)dpl,psvo11,psvo12,psvo13
      write(52,133)dpl,psvo21,psvo22,psvo23
      write(53,133)dpl,psvo31,psvo32,psvo33
      do k=1,k3max
      dpl=(k-1)/dx3*utap*ren
      pspre3=pcox3(j,k)/nba
      psvv11=r11x3(j,k)/nba
      psvv22=r22x3(j,k)/nba
      psvv33=r33x3(j,k)/nba
      psoo11=v11x3(j,k)/nba
      psoo22=v22x3(j,k)/nba
      psoo33=v33x3(j,k)/nba
      psoo12=v12x3(j,k)/nba
      psoo13=v31x3(j,k)/nba
      psoo23=v23x3(j,k)/nba
      psvv12=r12x3(j,k)/nba
      psvv13=r31x3(j,k)/nba
      psvv23=r23x3(j,k)/nba
      psvo11=vo11x3(j,k)/nba
      psvo12=vo12x3(j,k)/nba
      psvo13=vo13x3(j,k)/nba
      psvo21=vo21x3(j,k)/nba
      psvo22=vo22x3(j,k)/nba
      psvo23=vo23x3(j,k)/nba
      psvo31=vo31x3(j,k)/nba
      psvo32=vo32x3(j,k)/nba
      psvo33=vo33x3(j,k)/nba
      write(69,133)dpl,pspre3
      write(55,133)dpl,psvv11,psvv22,psvv33
      write(56,133)dpl,psoo11,psoo22,psoo33
      write(57,133)dpl,psvv12,psvv13,psvv23
      write(54,133)dpl,psoo12,psoo13,psoo23
      write(51,133)dpl,psvo11,psvo12,psvo13
      write(52,133)dpl,psvo21,psvo22,psvo23
      write(53,133)dpl,psvo31,psvo32,psvo33
      enddo
      close(69)
      close(58)
      close(56)
      close(55)
      close(57)
      close(54)
      close(53)
      close(52)
      close(51)
      npq=npq+1
                        endif
      enddo
		endif
                               endif
  983 format(3x,e11.4,i5,9(1x,e12.5))
  783 format(1x,11(e11.4))
  785 format(3x,7(1x,e12.5),3x,i5)
      write(59,785)time,utaul,utauu,utap,cfou,cfol,dp3mav
     1             ,nav
			endif
      vstzrn=-vmeo(3,n2-1)/(y(n2)-y2s(n2m))/nav/ren
      vsttrn=-vmeo(1,n2-1)/(y(n2)-y2s(n2m))/nav/ren
      tauwzr=sqrt(abs(vstzrn))
      tauwtr=sqrt(abs(vsttrn))
      cfzr=tauwzr**2*2./vit(3)**2
      cftr=tauwtr**2*2./vit(3)**2
      utat=sqrt(abs(dp3ns))
      retat=re*utap
      write(50,783)time,vit(5),(vit(l),l=1,3,2),enavp,utat
     1 ,dp3ns,cflm
      write(32,783)time,(vit(l),l=1,5),enavp,utat
     1 ,dp3ns
      write(6,783)time,cfzr,cftr,utat,utap
     1 ,enavp,retat
      do n=1,npqf
      write(6,787)time,jprq(n),dissma(n),enssma(n),esprma(n),espmi(n),espma(n)
      enddo
  787 format(e11.4,2x,i5,2x,6e12.4)
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
      if(l.eq.2) vca=q2(i,j,k)
      if(abs(vca).ge.vmaxo(l)) then
      vmaxo(l)=abs(vca)
                              endif
      enddo
      enddo
      enddo
      do k=1,n3mv
      do j=1,n2mv
      do i=1,n1mv
      if(l.eq.1) vca=q1(i,j,k)
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

