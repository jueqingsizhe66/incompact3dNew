c******************************* subrout wstre **********************
c
      subroutine wstre
c
      include 'param.f'
      common/wallst/cflw,cfuw
c
c     lower and upper walls shear
c     averaged in the homogeneous directions
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
      common/utapr/utap
      dimension skp(4),flp(4),vmp(3),tstp(6),vstp(6),vop(3)
      dimension pvp(4),pco(4),hel(3),vvo(6),pgrp(3)
      dimension dustr(6,m2),bud(5,6),stmedp(9),dtrco(4,m2),vtvp(4)
      common/nbal/nba,ibudg
      common/walfr/cfuo,cflo
      common/wallst/cflw,cfuw
      character*60 filebu1,filebu2,filebu3,filebu4,filebu5 
      character*60 filebp1,filebp2,filebp3,filebp4,filebp5 
      character*60 filestm,filevtv,filepvp,filedust,filetrco
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
c   evaluates the velocity correlations
c
      call veltur(ntime)
c
c   evaluates the vorticity correlations
c
      call vortur(ntime)
c
c     evaluates the time averages of the statistics
c
      call taver
c
c
      enavp=enavo/nav
      disstp=dissto/nav
c     write(6,*)' in outh nav=',nav,
c    1   ' number samples  ','   nvv=',nvv,
c    1    'for vertic. transp'
c
c    evaluates the derivatives of the statistics for the
c    budgets
c
      isign=1. 
      do j=2,n2m
      dustr(1,j)=((tursto(1,j)/nav)**2-(tursto(1,j-1)/nav)**2)
     1            /(y2s(j)-y2s(j-1))*isign
      dustr(2,j)=((tursto(2,j)/nav)**2-(tursto(2,j-1)/nav)**2)
     1            /(y2s(j)-y2s(j-1))*isign
      dustr(3,j)=((tursto(3,j)/nav)**2-(tursto(3,j-1)/nav)**2)
     1            /(y2s(j)-y2s(j-1))*isign
      dustr(4,j)=-(tursto(5,j)-tursto(5,j-1))/nav
     1            /(y2s(j)-y2s(j-1))*isign
      dustr(5,j)=(pvmo(2,j)-pvmo(2,j-1))/(y2s(j)-y2s(j-1))/nav*isign
      dustr(6,j)=(pvmo(3,j)-pvmo(3,j-1))/(y2s(j)-y2s(j-1))/nav*isign
      dtrco(1,j)=(vtvo(1,j)-vtvo(1,j-1))/(y2s(j)-y2s(j-1))/nav*isign
      dtrco(2,j)=(vtvo(2,j)-vtvo(2,j-1))/(y2s(j)-y2s(j-1))/nav*isign
      dtrco(3,j)=(vtvo(3,j)-vtvo(3,j-1))/(y2s(j)-y2s(j-1))/nav*isign
      dtrco(4,j)=(vtvo(4,j)-vtvo(4,j-1))/(y2s(j)-y2s(j-1))/nav*isign
      enddo
      dustr(1,1)=(tursto(1,1)/nav)**2/(y2s(1)-y(1))*isign
      dustr(2,1)=(tursto(2,1)/nav)**2/(y2s(1)-y(1))*isign
      dustr(3,1)=(tursto(3,1)/nav)**2/(y2s(1)-y(1))*isign
      dustr(4,1)=-tursto(5,1)/(y2s(1)-y(1))/nav*isign
      dustr(5,1)=pvmo(2,1)/(y2s(1)-y(1))/nav*isign
      dustr(6,1)=pvmo(3,1)/(y2s(1)-y(1))/nav*isign
      dtrco(1,1)=vtvo(1,1)/(y2s(1)-y(1))/nav*isign
      dtrco(2,1)=vtvo(2,1)/(y2s(1)-y(1))/nav*isign
      dtrco(3,1)=vtvo(3,1)/(y2s(1)-y(1))/nav*isign
      dtrco(4,1)=vtvo(4,1)/(y2s(1)-y(1))/nav*isign
      isign=1. 
      dustr(1,n2)=-(tursto(1,n2-1)/nav)**2/(y(n2)-y2s(n2m))*isign
      dustr(2,n2)=-(tursto(2,n2-1)/nav)**2/(y(n2)-y2s(n2m))*isign
      dustr(3,n2)=-(tursto(3,n2-1)/nav)**2/(y(n2)-y2s(n2m))*isign
      dustr(4,n2)=tursto(5,n2-1)/(y(n2)-y2s(n2m))/nav*isign
      dustr(5,n2)=-pvmo(2,n2-1)/(y(n2)-y2s(n2m))/nav*isign
      dustr(6,n2)=-pvmo(3,n2-1)/(y(n2)-y2s(n2m))/nav*isign
      dtrco(1,n2)=-vtvo(1,n2-1)/(y(n2)-y2s(n2m))/nav*isign
      dtrco(2,n2)=-vtvo(2,n2-1)/(y(n2)-y2s(n2m))/nav*isign
      dtrco(3,n2)=-vtvo(3,n2-1)/(y(n2)-y2s(n2m))/nav*isign
      dtrco(4,n2)=-vtvo(4,n2-1)/(y(n2)-y2s(n2m))/nav*isign
      itim=nint(time)
      write(pntim,77) itim
   77 format(i4.4)
c
c   print the budgets and some statistics
c
      if(mod(ntime,npouth).eq.0) then
      filebu1='u1-bud.'//pntim
      filebu2='u2-bud.'//pntim
      filebu3='u3-bud.'//pntim
      filebu4='u23-bud.'//pntim
      filebu5='k-bud.'//pntim
      filebp1='u1-bud.plo'
      filebp2='u2-bud.plo'
      filebp3='u3-bud.plo'
      filebp4='u23-bud.plo'
      filebp5='k-bud.plo'
      filestm='stmed.'//pntim
      filevtv='vtvmed.'//pntim
      filepvp='pvpmed.'//pntim
      filedust='dustmed.'//pntim
      filetrco='dtrcmed.'//pntim
      filns='turstr.'//pntim
      open(85,file=filebp5)
      open(81,file=filebp1)
      open(82,file=filebp2)
      open(83,file=filebp3)
      open(84,file=filebp4)
      open(95,file=filebu5)
      open(91,file=filebu1)
      open(92,file=filebu2)
      open(93,file=filebu3)
      open(94,file=filebu4)
      open(96,file=filestm)
      open(98,file=filevtv)
      open(99,file=filepvp)
      open(88,file=filedust)
      open(89,file=filetrco)
      open(97,file=filns)
      scalbu=1./(re*utap**4)
      write(66,*)'scalbu= ',scalbu
      do jc=1,n2
      write(88,100)yp2(jc),(dustr(l,jc),l=1,6) 
      write(89,100)yp2(jc),(dtrco(l,jc),l=1,4) 
       enddo
      do jc=1,n2m
      if(y2s(jc).le.0.) then
      yd=y2s(jc)-y(1)
      isign=1. 
                     else
      yd=y(n2)-y2s(jc)
      isign=1. 
                     endif
      budgo(3,3,jc)=-(dtrco(3,jc+1)+dtrco(3,jc))*0.5
      budgo(3,4,jc)=2.*tursto(5,jc)/nav*stmedo(8,jc)/nav
      budgo(3,5,jc)=(dustr(3,jc+1)-dustr(3,jc))/(yp2(jc+1)-yp2(jc))/re
      budgo(3,6,jc)=0.
      budgo(2,3,jc)=-(dtrco(2,jc+1)+dtrco(2,jc))*0.5
      budgo(2,4,jc)=0.
      budgo(2,6,jc)=-2.*(dustr(5,jc+1)+dustr(5,jc))*0.5
      budgo(2,5,jc)=(dustr(2,jc+1)-dustr(2,jc))/(yp2(jc+1)-yp2(jc))/re
      budgo(1,3,jc)=-(dtrco(1,jc+1)+dtrco(1,jc))*0.5
      budgo(1,6,jc)=0.
      budgo(1,4,jc)=0.
      dust11=(dustr(1,jc+1)-dustr(1,jc))/(yp2(jc+1)-yp2(jc))
      budgo(1,5,jc)=dust11/re
      budgo(4,3,jc)=-(dtrco(4,jc+1)+dtrco(4,jc))*0.5
      budgo(4,4,jc)=-(tursto(2,jc)/nav)**2*stmedo(8,jc)/nav*isign
      budgo(4,5,jc)=(dustr(4,jc+1)-dustr(4,jc))/(yp2(jc+1)-yp2(jc))/re
      budgo(4,6,jc)=(dustr(6,jc+1)+dustr(6,jc))*0.5
c
c   the budgets are scaled as in Moser et al (1987)
c   ne = 1 spanwise               u1u1
c   ne = 2 normal                 u2u2
c   ne = 3 streamwise             u3u3
c   ne = 4 u2u3 correlations
c   ne = 5 kinetic energy         k
c
      do ne=1,4
      bud(ne,2)=budgo(ne,2,jc)/nav/re*scalbu
      bud(ne,1)=budgo(ne,1,jc)/nav*scalbu
      do nt=3,6
      bud(ne,nt)=budgo(ne,nt,jc)*scalbu
      enddo
      enddo
      do ns=1,9
      stmedp(ns)=stmedo(ns,jc)/nav*isign
      enddo
      do nt=1,6
      bud(5,nt)=(bud(1,nt)+bud(2,nt)+bud(3,nt))*0.5
      enddo
      do  l=1,4
      vtvp(l)=(vtvo(l,jc)/nav)
      pvp(l)=(pvmo(l,jc)/nav)
      enddo
      do  l=1,3
      tstp(l)=(tursto(l,jc)/nav)
      enddo
      do  l=4,6
      tstp(l)=tursto(l,jc)/nav
      enddo
c     prod3=2.*tstp(5)*stmedp(8)
c     prod3n=prod3*scalbu
c     write(66,*)jc,prod3,budgo(3,4,jc),prod3n,bud(3,4)
      write(95,100)yd,(bud(5,nt),nt=1,6)
      write(94,100)yd,(bud(4,nt),nt=1,6)
      write(93,100)yd,(bud(3,nt),nt=1,6)
      write(92,100)yd,(bud(2,nt),nt=1,6)
      write(91,100)yd,(bud(1,nt),nt=1,6)
      write(96,100)yd,(stmedp(ns),ns=1,9)
      write(99,100)yd,(pvp(ns),ns=1,4)
      write(98,100)yd,(vtvp(ns),ns=1,4)
      write(97,100)yd,(tstp(ns),ns=1,6)
      if(y2s(jc).le.0.) then
      yplu=yd*re*utap
      write(85,100)yplu,(bud(5,nt),nt=1,6)
      write(84,100)yplu,(bud(4,nt),nt=1,6)
      write(83,100)yplu,(bud(3,nt),nt=1,6)
      write(82,100)yplu,(bud(2,nt),nt=1,6)
      write(81,100)yplu,(bud(1,nt),nt=1,6)
                        endif
      enddo
100   format(3x,10(1x,e12.5))
101   format(3x,8(1x,e12.5))
      close(85)
      close(81)
      close(82)
      close(83)
      close(84)
      close(95)
      close(91)
      close(92)
      close(93)
      close(94)
      close(96)
      close(97)
      close(98)
      close(89)
      close(88)
      close(99)
                                 endif
c
c   evaluates some global quantities for each field read
c
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
      write(32,783)time,(vit(l),l=1,5),enavp,utat,utap
     1 ,dp3ns
      write(6,783)time,cfzr,cftr,utat,utap
     1 ,enavp,retat
      do n=1,npqf
      write(6,787)time,jprq(n),dissma(n),enssma(n),esprma(n),espmi(n),espma(n)
      enddo
  787 format(e11.4,2x,i5,2x,6e12.4)
  983 format(3x,e11.4,i5,9(1x,e12.5))
  783 format(1x,11(e11.4))
  785 format(3x,7(1x,e12.5),3x,i5)
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
c     find the maximum velocities 
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

