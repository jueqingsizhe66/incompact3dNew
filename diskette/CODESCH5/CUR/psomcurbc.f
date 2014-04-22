c
c  ****************************** subrout bcnsvi  **********************
c
c  this subroutine calculate the b.c for psi and vor for j=1 and j=n2
c
      subroutine bcnsvi(gam,ro,vor,psi,y,dvor,psc,ft)
      include 'param.f'
      dimension psc(m1,m2)
      dimension vor(m1,m2),psi(m1,m2),y(ndd,m1,m2),dvor(m1,m2)
      dimension q1s(m1,3),q2s(m1,3),q1ss(m1),q2ss(m1,3)
      dimension q1n(m1,3),q2n(m1,3),q1ns(m1),q2ns(m1,3)
      common/vorwal/vorlo(m1),vorup(m1)
      common/dvowal/dvorbs(m1),dvorbn(m1)
      common/dim/n1,n1m,n2,n2m
      common/mesh/dx1,dx1q,dx2,dx2q
      common/metrbw/gwa112(ndd,m1)
      common/psbou/psbi(2,m1),psbj(2,m2)
      common/sliwal/insls,insln
      common/metrip/alfi(ndd,ndd,m1,m2),alfj(ndd,ndd,m1,m2)
      common/metrst/gccc(m1,m2)
      common/vopsin/psinf(m2),vobj(m2)
      common/uiwai/uinfi(m1),vinfi(m1)
      common/uiwal/uinfs(m1),vinfs(m1)
      common/slewa/xls(m1)
      common/lewai/xl1,xl2
      common/pscwal/psclo(m1),pscup(m1)
      common/dpswal/dpscbs(m1),dpscbn(m1)
      common/psiwal/psins(m1)
      common/bnool/dvnoro(m1),dpnoro(m1),dpsnoo(m1)
      common/bdqnor/dpnor(m1)
      common/convel/cout
      common/velcve/cven(3,m1)
      common/tstep/dt
      dimension cn1(2),cn2(2)
      if(insls.eq.0) then
c
c   free-slip condition south  boundary
c
      do i=2,n1m
      vorjb=-gwa112(1,i)*(psi(i,2)-psi(i,1))*dx2
      dvorbs(i)=vorjb-vorlo(i)
      vorlo(i)=vorjb
      enddo
      do i=1,n1
      psbi(1,i)=0.
      psi(i,1)=psbi(1,i)
      enddo
                     endif
      if(insls.eq.1) then
c
c  eventual localised transpiration
c
      do i=1,n1
      vinfs(i)=vinfi(i)*ft
      uinfs(i)=uinfi(i)*ft
      enddo
      j=1
      js=1
      psbi(1,1)=0.
      do i =1,n1m
      do nd=1,ndd
      cn2(nd)=-js*((y(nd,i,j)-y(nd,i,j+js))
     1            +(y(nd,i+1,j)-y(nd,i+1,j+js)))*dx2*0.5
      cn1(nd)=(y(nd,i+1,j)-y(nd,i,j))*dx1
      enddo
      q2=cn2(1)*uinfs(i)+cn2(2)*vinfs(i)
      psbi(1,i+1)=psbi(1,i)-q2/dx2
      enddo
      do i=1,n1
      psi(i,1)=psbi(1,i)
      enddo
c
c   no-slip condition south  boundary
c
      do j=1,2
      do i=1,n1
      q1s(i,j)=(psi(i,j+1)-psi(i,j))*dx2
      enddo
      enddo
      do j=1,3
      do i=1,n1m
      q2s(i,j)=-(psi(i+1,j)-psi(i,j))*dx1
      enddo
      enddo
      do i=1,n1m
      q1ss(i)=(q1s(i+1,2)+q1s(i,2)
     1        +q1s(i,1)+q1s(i+1,1))*0.25
      enddo
      do j=1,2
      do i=2,n1m
      q2ss(i,j)=(q2s(i-1,j+1)+q2s(i,j+1)
     1          +q2s(i-1,j)+q2s(i,j))*0.25
      enddo
      enddo
      do i=2,n1m
      ual2=-1./gccc(i,2)
      ual1=-1./gccc(i,1)
      dq21=-(ual2*(alfi(1,1,i,2)*q2s(i,2)
     1           -alfi(1,1,i-1,2)*q2s(i-1,2))
     1      -ual1*(alfi(1,1,i,1)*q2s(i,1)
     1           -alfi(1,1,i-1,1)*q2s(i-1,1)))*dx1
      dq11=+ual2*(alfi(1,2,i,2)*q1ss(i)
     1           -alfi(1,2,i-1,2)*q1ss(i-1))*dx1
      dq22=-(ual2*(alfj(2,1,i,2)*q2ss(i,2)
     1            -alfj(2,1,i,1)*q2ss(i,1))
     1    -2.*ual1*alfj(2,1,i,1)*q2ss(i,1)) *dx2
      dq12=-(ual2*(alfj(2,2,i,2)*q1s(i,2)
     1            -alfj(2,2,i,1)*q1s(i,1))
     1    -2.*ual1*alfj(2,2,i,1)*q1s(i,1))  *dx2
      vorjb=dq21+dq11+dq22+dq12
      dvorbs(i)=(vorjb-vorlo(i))
      vorlo(i)=vorjb
      enddo
                endif
c
c    vorjb     is vor(i,2 -vor(i,1
c    dvorbs     is dvor(i,1 -dvor(i,r21
c
c
c   passive scalar inlet on the wall south bound.
c
      dpscma=0.
      do i =1,n1
      psins(i)=0.
      dpscbs(i)=0.
      if(xls(i).ge.xl1.and.xls(i).le.xl2) then
      xd=(xls(i)-xl1)/(xl2-xl1)         
      psins(i)=xd*(1.-xd)/0.25*ft
      dpscbs(i)=psins(i)-psc(i,1)
      dpscma=max(dpscma,dpscbs(i))
                                endif
      enddo
c     write(68,*)'ft,dpscma',ft,dpscma
      if(insln.eq.0) then
c
c   free-slip condition north boundary
c
      do i=2,n1m
      vorjb=-gwa112(2,i)*(psi(i,n2)-psi(i,n2m))*dx2
      dvorbn(i)=vorjb-vorup(i)
      vorup(i)=vorjb
      enddo
      do i=1,n1
      psbi(2,i)=psbj(1,n2)
      psi(i,n2)=psbi(2,i)
      enddo
                endif
      if(insln.eq.1) then
c
c   no-slip condition north boundary
c
      do jn=1,2
      j=n2-jn
      q2n(1,jn)=0.
      q2n(n1,jn)=0.
      do i=1,n1
      q1n(i,jn)=(psi(i,j+1)-psi(i,j))*dx2
      enddo
      enddo
      do jn=1,3
      j=n2-jn+1
      do i=1,n1m
      q2n(i,jn)=-(psi(i+1,j)-psi(i,j))*dx1
      enddo
      enddo
      do i=1,n1m
      q1ns(i)=(q1n(i,2)+q1n(i+1,2)
     1        +q1n(i,1)+q1n(i+1,1))*0.25
      enddo
      do j=1,2
      do i=2,n1m
      q2ns(i,j)=(q2n(i-1,j+1)+q2n(i,j+1)
     1          +q2n(i-1,j)+q2n(i,j))*0.25
      enddo
      enddo
      do i=2,n1m
      ual2=-1./gccc(i,n2m)
      ual1=-1./gccc(i,n2)
      dq21=-ual2*(alfi(1,1,i,n2m)*q2n(i,2)
     1           -alfi(1,1,i-1,n2m)*q2n(i-1,2))*dx1
      dq11=+ual2*(alfi(1,2,i,n2m)*q1ns(i)
     1           -alfi(1,2,i-1,n2m)*q1ns(i-1))*dx1
      dq22=-(ual2*(alfj(2,1,i,n2m-1)*q2ns(i,2)
     1            -alfj(2,1,i,n2m)*q2ns(i,1))
     1    -2.*ual1*alfj(2,1,i,n2m)*q2ns(i,1)) *dx2
      dq12=-(ual2*(alfj(2,2,i,n2m-1)*q1n(i,2)
     1            -alfj(2,2,i,n2m)*q1n(i,1))
     1    -2.*ual1*alfj(2,2,i,n2m)*q1n(i,1))  *dx2
      vorjb=dq21+dq11+dq22+dq12
      dvorbn(i)=-(vorjb-vorup(i))
      vorup(i)=vorjb
      enddo
      do i=1,n1
      psbi(2,i)=psbj(1,n2)
      psi(i,n2)=psbi(2,i)
      enddo
c
c    vorjb     is vor(i,n2m -vor(i,n2
c    dvorbn     is dvor(i,n2 -dvor(i,n2m
c
                     endif
      if(insln.eq.2) then
      do i=1,n1
c
c    radiative b.c. for dvor
c
        dvoox2 = (vor(i,n2)-vor(i,n2m))*dx2
        dvorbn(i) = -dt*(gam*dvoox2+ro*dvnoro(i))*cven(1,i)
        dvnoro(i)=dvoox2
c
c    radiative b.c. for dpsi
c 
        dpsox2 = (psi(i,n2)-psi(i,n2m))*dx2
        dpnor(i) = -dt*(gam*dpsox2+ro*dpnoro(i))*cven(2,i)
        dpnoro(i)=dpsox2
        psbi(2,i)=psi(i,n2)+dpnor(i)
c
c    radiative b.c. for dpsc
c
        dpscx2 = (psc(i,n2)-psc(i,n2m))*dx2
        dpscbn(i) = -dt*(gam*dpscx2+ro*dpsnoo(i))*cven(3,i)
        dpsnoo(i)=dpscx2
      enddo
                     endif
      return
      end
c
c  ****************************** subrout bcwest  **********************
c
c  this subroutine calculate the b.c for psi and vor for i=1 
c  west boundary
c
      subroutine bcwest(gam,ro,vor,psi,y,dvor,psc,ft)
      include 'param.f'
      dimension psc(m1,m2)
      dimension vor(m1,m2),psi(m1,m2),y(ndd,m1,m2),dvor(m1,m2)
      common/dim/n1,n1m,n2,n2m
      common/mesh/dx1,dx1q,dx2,dx2q
      common/metrip/alfi(ndd,ndd,m1,m2),alfj(ndd,ndd,m1,m2)
      common/metrst/gccc(m1,m2)
      common/vopsin/psinf(m2),vobj(m2)
      common/pscinfl/pscinf(m2)
      common/uinfl/uinfl(m2),vinfl(m2)
      common/inwes/inwest 
      common/tstep/dt
      common/bdqinl/dvwes(m2),dpwes(m2),dpswes(m2)
      common/psbou/psbi(2,m1),psbj(2,m2)
      common/bqinl/vwes(m2),pwes(m2),pswes(m2)
      common/binol/dvweso(m2),dpweso(m2),dpsweo(m2)
      common/convel/cout
      common/velcho/choe(3,m2),chow(3,m2)
      dimension cn1(2),cn2(2)
      if(inwest.eq.0) then
c
c  inlet conditions
c
      do j=1,n2
      psi(1,j)=psbj(1,j)
      dvor(1,j)=0.
      vor(1,j)=vobj(j)
      enddo
      do j=1,n2
      psc(1,j)=0.
      enddo
                      endif
      if(inwest.eq.1) then
      do j=1,n2
c
c    radiative b.c. for dvor
c
        dvoox1 = (vor(1,j)-vor(2,j))*dx1
        dvwes(j) = -dt*(gam*dvoox1+ro*dvweso(j))*chow(1,j)
        dvweso(j)=dvoox1
c
c    radiative b.c. for dpsi
c 
        dpsox1 = (psi(1,j)-psi(2,j))*dx1
        dpwes(j) = -dt*(gam*dpsox1+ro*dpweso(j))*chow(2,j)
        dpweso(j)=dpsox1
        psbj(1,j)=psi(1,j)+dpwes(j)
c
c    radiative b.c. for dpsc
c
        dpscx1 = (psc(1,j)-psc(2,j))*dx1
        dpswes(j) = -dt*(gam*dpscx1+ro*dpsweo(j))*chow(3,j)
        dpsweo(j)=dpscx1
c     if(j.eq.5) write(21,133)j,dvwes(j),dpsox1,dpwes(j),psbj(1,j)
c    1           ,dpswes(j)
  133 format(3x,'wes',i4,3x,5e12.4)
      enddo
                      endif
      return
      end
c
c  **********************  set boundary conditions  *********
c   est boundary for i=n1
c
      subroutine bcest(gam,ro,vor,psi,psc)
      include 'param.f'
      dimension vor(m1,m2),psi(m1,m2),psc(m1,m2)
      common/dim/n1,n1m,n2,n2m
      common/tstep/dt
      common/mesh/dx1,dx1q,dx2,dx2q
      common/bdqout/dvest(m2),dpest(m2),dpsest(m2)
      common/psbou/psbi(2,m1),psbj(2,m2)
      common/bqout/vest(m2),pest(m2),psest(m2)
      common/bouol/dvesto(m2),dpesto(m2),dpseso(m2)
      common/convel/cout
      common/velcho/choe(3,m2),chow(3,m2)
      do j=1,n2
c
c    radiative b.c. for dvor
c
        dvoox1 = (vor(n1,j)-vor(n1m,j))*dx1
        dvest(j) = -dt*(gam*dvoox1+ro*dvesto(j))*choe(1,j)
        dvesto(j)=dvoox1
c
c    radiative b.c. for dpsi
c 
        dpsox1 = (psi(n1,j)-psi(n1m,j))*dx1
        dpest(j) = -dt*(gam*dpsox1+ro*dpesto(j))*choe(2,j)
        dpesto(j)=dpsox1
        psbj(2,j)=psi(n1,j)+dpest(j)
c
c    radiative b.c. for dpsc
c
        dpscx1 = (psc(n1,j)-psc(n1m,j))*dx1
        dpsest(j) = -dt*(gam*dpscx1+ro*dpseso(j))*choe(3,j)
        dpseso(j)=dpscx1
c     if(j.eq.5) write(21,133)j,dvest(j),
c    1     dpsox1,dpest(j),psbj(2,j),dpsest(j)
  133 format(3x,'est',i4,3x,5e12.4)
      enddo
      return
      end
c
c  **********************  set values at the old time step *********
c   first derivatives along x2 for i=1-n1 
c    north boundary
c
      subroutine cvnool(vor,psi,psc)
      include 'param.f'
      dimension vor(m1,m2),psi(m1,m2),psc(m1,m2)
      common/dim/n1,n1m,n2,n2m
      common/mesh/dx1,dx1q,dx2,dx2q
      common/firdve/dvox2n(m1),dpsx2n(m1),dssx2n(m1)
      common/qveold/voven(m1),psven(m1),ssven(m1)
      common/metrst/gccc(m1,m2)
      common/sliwal/insls,insln
c
c    quantities at the old time step
c    at the bound. on the upper wall
c
      if(insln.eq.2) then
      do i=1,n1
      udnn2=dx2*0.5/gccc(i,n2m)
      dvox2n(i)=(vor(i,n2)-vor(i,n2-2))*udnn2
      voven(i)=vor(i,n2m)
      dpsx2n(i)=(psi(i,n2)-psi(i,n2-2))*udnn2
      psven(i)=psi(i,n2m)
      dssx2n(i)=(psc(i,n2)-psc(i,n2-2))*udnn2
      ssven(i)=psc(i,n2m)
      enddo
                     endif
      return
      end
c
c  **********************  set values at the old time step *********
c   first derivatives along x1 for j=1-n2 
c    west boundary
c
      subroutine cvweol(vor,psi,psc)
      include 'param.f'
      dimension vor(m1,m2),psi(m1,m2),psc(m1,m2)
      common/dim/n1,n1m,n2,n2m
      common/mesh/dx1,dx1q,dx2,dx2q
      common/firdho/dvox1e(m2),dpsx1e(m2),dvox1w(m2),dpsx1w(m2)
     1             ,dssx1e(m2),dssx1w(m2)
      common/qhoold/vohoe(m2),pshoe(m2),vohow(m2),pshow(m2)
     1             ,sshoe(m2),sshow(m2)
      common/metrst/gccc(m1,m2)
      common/inwes/inwest 
c
c    quantities at the old time step
c    at the bound. on the right
c
      do j=1,n2
      udnn1=dx1*0.5/gccc(n1m,j)
      dvox1e(j)=(vor(n1,j)-vor(n1-2,j))*udnn1
      vohoe(j)=vor(n1m,j)
      dpsx1e(j)=(psi(n1,j)-psi(n1-2,j))*udnn1
      pshoe(j)=psi(n1m,j)
      dssx1e(j)=(psc(n1,j)-psc(n1-2,j))*udnn1
      sshoe(j)=psc(n1m,j)
      enddo
c
c    quantities at the old time step
c    at the bound. on the left
c
      if(inwest.eq.1) then
      do j=1,n2
      udnn1=dx1*0.5/gccc(2,j)
      dvox1w(j)=-(vor(3,j)-vor(1,j))*udnn1
      vohow(j)=vor(2,j)
      dpsx1w(j)=-(psi(3,j)-psi(1,j))*udnn1
      pshow(j)=psi(2,j)
      dssx1w(j)=-(psc(3,j)-psc(1,j))*udnn1
      sshow(j)=psc(2,j)
      enddo
                      endif
      return
      end
c
c  ******* computes radiation velocity  *********
c  ******* at the vertical boundaries
c
      subroutine cvelno(vor,psi,psc,y)
c
c  calculation of b.c radiative velocity Orlansky 
c  computational velocity  to be used in the routine
c   bcnsvi
      include 'param.f'
      dimension y(ndd,m1,m2)
      dimension vor(m1,m2),psi(m1,m2),psc(m1,m2)
      common/dim/n1,n1m,n2,n2m
      common/tstep/dt
      common/mesh/dx1,dx1q,dx2,dx2q
      common/firdve/dvox2n(m1),dpsx2n(m1),dssx2n(m1)
      common/qveold/voven(m1),psven(m1),ssven(m1)
      common/velcve/cven(3,m1)
      common/metrst/gccc(m1,m2)
      common/sliwal/insls,insln
      common/convel/cout
      dimension cvnoma(3)
c
c
c  calculation of radiative velocity for omega
c  at the boundary on the upper boundary north
 
c
      compv=1./(dt*dx2)
      if(insln.eq.2) then
      do l=1,3
      cvnoma(l)=0.
      enddo
      do i=1,n1
      udnn2=dx2*0.5/gccc(i,n2m)
c     compv=(y(2,i,n2)-y(2,i,n2m))/dt
      dqx2n=(vor(i,n2)-vor(i,n2-2))*udnn2
      dqt=(vor(i,n2m)-voven(i))/dt
      dqx2=(dqx2n+dvox2n(i))*0.5
       if(abs(dqx2).le..1e-05) dqx2=.1e-05
         covel1=-dqt/dqx2
         if(covel1.gt.0.)then
           if(covel1.gt.compv) covel1=compv
           cven(1,i)=covel1
                           else
           cven(1,i)=0.
                          endif
c      cven(1,i)=cout
      cvnoma(1)=max(cvnoma(1),cven(1,i))
      enddo
c
c  calculation of radiative velocity for psi
c        upper boundary north
      do i=1,n1
      udnn2=dx2*0.5/gccc(i,n2m)
c     compv=(y(2,i,n2)-y(2,i,n2m))/dt
      dqx2n=(psi(i,n2)-psi(i,n2-2))*udnn2
      dqt=(psi(i,n2m)-psven(i))/dt
      dqx2=(dqx2n+dpsx2n(i))*0.5
       if(abs(dqx2).le..1e-05) dqx2=.1e-05
         covel2=-dqt/dqx2
         if(covel2.gt.0.)then
           if(covel2.gt.compv) covel2=compv
           cven(2,i)=covel2
                           else
           cven(2,i)=0.
                          endif
c      cven(2,i)=cout
      cvnoma(2)=max(cvnoma(2),cven(2,i))
      enddo
c
c  calculation of radiative velocity for passive scalar
c  at the boundary on the upper boundary north
 
c
      do i=1,n1
      udnn2=dx2*0.5/gccc(i,n2m)
c     compv=(y(2,i,n2)-y(2,i,n2m))/dt
      dqx2n=(psc(i,n2)-psc(i,n2-2))*udnn2
      dqt=(psc(i,n2m)-ssven(i))/dt
      dqx2=(dqx2n+dssx2n(i))*0.5
       if(abs(dqx2).le..1e-05) dqx2=.1e-05
         covel1=-dqt/dqx2
         if(covel1.gt.0.)then
           if(covel1.gt.compv) covel1=compv
           cven(3,i)=covel1
                           else
           cven(3,i)=0.
                          endif
c      cven(3,i)=cout
      cvnoma(3)=max(cvnoma(3),cven(3,i))
      enddo
                          endif
c     write(69,103)(cvnoma(l),l=1,3),compv
  103 format(3x,'no',4e12.4)
      return
      end
c
c  **********************  computes radiation velocity  *********
c  ******* at the west boundary to be used in bcwest 
c
      subroutine cvelwe(vor,psi,psc,y)
      include 'param.f'
      dimension y(ndd,m1,m2)
      dimension vor(m1,m2),psi(m1,m2),psc(m1,m2)
      common/dim/n1,n1m,n2,n2m
      common/tstep/dt
      common/mesh/dx1,dx1q,dx2,dx2q
      common/firdho/dvox1e(m2),dpsx1e(m2),dvox1w(m2),dpsx1w(m2)
     1             ,dssx1e(m2),dssx1w(m2)
      common/qhoold/vohoe(m2),pshoe(m2),vohow(m2),pshow(m2)
     1             ,sshoe(m2),sshow(m2)
      common/metrst/gccc(m1,m2)
      common/inwes/inwest 
      common/velcho/choe(3,m2),chow(3,m2)
      dimension cvwema(3),cvesma(3)
      common/convel/cout
c
c
c  calculation of b.c radiative velocity Orlansky JFM pp.254
c
c
c  calculation of radiative velocity for omega
c  at the boundary on the right
 
c
      compv=1./(dt*dx1)
      do l=1,3
      cvesma(l)=0.
      enddo
      do j=1,n2
      udnn1=dx1*0.5/gccc(n1m,j)
c     compv=(y(1,n1,j)-y(1,n1m,j))/dt
      dq1x1e=(vor(n1,j)-vor(n1-2,j))*udnn1
      dq1t=(vor(n1m,j)-vohoe(j))/dt
      dq1x1=(dq1x1e+dvox1e(j))*0.5
       if(abs(dq1x1).le..1e-05) dq1x1=.1e-05
         covel1=-dq1t/dq1x1
         if(covel1.gt.0.)then
           if(covel1.gt.compv) covel1=compv
           choe(1,j)=covel1
                           else
           choe(1,j)=0.
                          endif
           choe(1,j)=cout
      cvesma(1)=max(cvesma(1),choe(1,j))
      enddo
c
c  calculation of radiative velocity for psi
c
      do j=1,n2
      udnn1=dx1*0.5/gccc(n1m,j)
c     compv=(y(1,n1,j)-y(1,n1m,j))/dt
      dq2x1e=(psi(n1,j)-psi(n1-2,j))*udnn1
      dq2t=(psi(n1m,j)-pshoe(j))/dt
      dq2x1=(dq2x1e+dpsx1e(j))*0.5
       if(abs(dq2x1).le..1e-05) dq2x1=.1e-05
         covel2=-dq2t/dq2x1
         if(covel2.gt.0.)then
           if(covel2.gt.compv) covel2=compv
           choe(2,j)=covel2
                           else
           choe(2,j)=0.
                           endif       
           choe(2,j)=cout
      cvesma(2)=max(cvesma(2),choe(2,j))
      enddo
c
c  calculation of radiative velocity for psc
c
      do j=1,n2
      udnn1=dx1*0.5/gccc(n1m,j)
c     compv=(y(1,n1,j)-y(1,n1m,j))/dt
      dq2x1e=(psc(n1,j)-psc(n1-2,j))*udnn1
      dq2t=(psc(n1m,j)-sshoe(j))/dt
      dq2x1=(dq2x1e+dssx1e(j))*0.5
       if(abs(dq2x1).le..1e-05) dq2x1=.1e-05
         covel2=-dq2t/dq2x1
         if(covel2.gt.0.)then
           if(covel2.gt.compv) covel2=compv
           choe(3,j)=covel2
                           else
           choe(3,j)=0.
                           endif       
           choe(3,j)=cout
      cvesma(3)=max(cvesma(3),choe(3,j))
      enddo
      write(69,103)(cvesma(l),l=1,3),compv
  103 format(3x,'es',4e12.4)
c
c  calculation of radiative velocity for omega
c  at the boundary on the left
c
 
c
      if(inwest.eq.1) then
      compv=1./(dt*dx1)
      do l=1,3
      cvwema(l)=0.
      enddo
      do j=1,n2
      udnn1=dx1*0.5/gccc(2,j)
c     compv=(y(1,2,j)-y(1,1,j))/dt
      dq1x1e=-(vor(3,j)-vor(1,j))*udnn1
      dq1t=(vor(2,j)-vohow(j))/dt
      dq1x1=-(dq1x1e+dvox1w(j))*0.5
       if(abs(dq1x1).le..1e-05) dq1x1=.1e-05
         covel1=-dq1t/dq1x1
         if(covel1.lt.0.)then
           if(abs(covel1).gt.compv) covel1=compv
           chow(1,j)=covel1
                           else
           chow(1,j)=0.
                          endif
           chow(1,j)=cout
      cvwema(1)=max(cvwema(1),chow(1,j))
      enddo
c
c  calculation of radiative velocity for psi
c
      do j=1,n2
      udnn1=dx1*0.5/gccc(2,j)
c     compv=(y(1,2,j)-y(1,1,j))/dt
      dq2x1e=-(psi(3,j)-psi(1,j))*udnn1
      dq2t=(psi(2,j)-pshow(j))/dt
      dq2x1=-(dq2x1e+dpsx1w(j))*0.5
       if(abs(dq2x1).le..1e-05) dq2x1=.1e-05
         covel2=-dq2t/dq2x1
         if(covel2.lt.0.)then
           if(abs(covel2).gt.compv) covel2=compv
           chow(2,j)=covel2
                           else
           chow(2,j)=0.
                          endif
           chow(2,j)=cout
      cvwema(2)=max(cvwema(2),chow(2,j))
      enddo
c
c  calculation of radiative velocity for psc
c
      do j=1,n2
      udnn1=dx1*0.5/gccc(2,j)
c     compv=(y(1,2,j)-y(1,1,j))/dt
      dq2x1e=-(psc(3,j)-psc(1,j))*udnn1
      dq2t=(psc(2,j)-sshow(j))/dt
      dq2x1=-(dq2x1e+dssx1w(j))*0.5
       if(abs(dq2x1).le..1e-05) dq2x1=.1e-05
         covel2=-dq2t/dq2x1
         if(covel2.lt.0.)then
           if(abs(covel2).gt.compv) covel2=compv
           chow(3,j)=covel2
                           else
           chow(3,j)=0.
                          endif
           chow(3,j)=cout
      cvwema(3)=max(cvwema(3),chow(3,j))
      enddo
                     endif
      write(69,105)(cvwema(l),l=1,3),compv
  105 format(3x,'we',4e12.4)
      return
      end
c
c  ******* initial convection velocities
c   are used in the first step since the
c   derivatives at the old step are unknown
c
      subroutine  incvel
      include 'param.f'
      common/dim/n1,n1m,n2,n2m
      common/velcho/choe(3,m2),chow(3,m2)
      common/velcve/cven(3,m1)
      common/convel/cout
      do l=1,3
      do i=1,n1
      cven(l,i)=cout
      enddo
      enddo
      do l=1,3
      do j=1,n2
      choe(l,j)=cout
      chow(l,j)=cout
      enddo
      enddo
      return
      end
