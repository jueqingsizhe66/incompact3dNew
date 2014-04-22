      program main
      include 'param.f'
      common/d2/nstop,nprint,ntst,npin,nmolp,ntra,mnpi
      common/dim/n1,n1m,n2,n2m
      common/wrre/nwrit,nread,nfilr
      common/tscoe/ga(3),ro(3),nsst
      common/parmod/yc1mo,yc2mo,akmo,ramo,velmo,vsi
      common/chnlc/chal,chbe,chga
      common/tstep/dt
      common/angmd/thet0
      common/reyn/beta,re
      common/mlv/mlev,iwmg,mlw,maxcmg,epsm
      common/maxit/nmaxi
      common/omsor/omeg
      common/corbum/bum1,sigb1,bum2,sigb2,xp1bu,xp2bu
      common/cordat/alx1i,alx1f,alx2i,alx2f
      common/sliwal/insls,insln
      common/inwes/inwest
      common/convel/cout
      common/lewai/xl1,xl2
      common/vistr/vampu,vampv
      common/icase/imod,istrea
      common/timinf/tau1
      common/schm/schmid
      common/strx1/str1,xcr1,etr1,istr1
      common/strx2/str2,xcr2,etr2,istr2
      common/icflp/icflm,cflma,dtl
      common/tfinst/tprin,tpin,tfin
      common/itypco/itypc
      open(15,file='psomcur.d')
      read(15,*)n1,n2,nsst
      read(15,*)nwrit,nread,nfilr,nplot,mnpi
      read(15,*)dt,ntst   ,nprint,npin,nmolp,ntra
      read(15,*)icflm,cflma,dtl
      read(15,*)tprin,tpin,tfin
      read(15,*)nstop
      read(15,*)chal,chbe,chga
      read(15,*)insls,insln,inwest
      read(15,*)re,schmid
      read(15,*)imod,istrea,vampu,vampv,xl1,xl2,tau1,cout      
      read(15,*)yc1mo,yc2mo,ramo,velmo,vsi
      read(15,*)thet0
      read(15,*)alx1i,alx1f,alx2i,alx2f
      read(15,*)str1,xcr1,etr1,istr1
      read(15,*)str2,xcr2,etr2,istr2
      read(15,*)itypc
      read(15,*)bum1,sigb1,xp1bu
      read(15,*)bum2,sigb2,xp2bu
      read(15,*)mlev,iwmg,mlw,maxcmg,nmaxi,omeg,epsm
      akmo=3.83/ramo
      pi=2.*asin(1.)
      write(6,201)thet0
      thet0=thet0*pi/180.
  201 format(10x,' dipole  vort stream funct. tht0=',e10.4)
      n1m=n1-1
      n2m=n2-1
      if(nsst.eq.1) go to 10
      ga(1)=8./15.
      ga(2)=5./12.
      ga(3)=3./4.
      ro(1)=0.
      ro(2)=-17./60.
      ro(3)=-5./12.
      write(6,100) (ga(n),n=1,nsst),(ro(n),n=1,nsst)
  100 format(10x,'third order runge gutta=',4x,
     1       'gam=',3f8.3,4x,'ro',3f8.3)
      go to 12
   10 continue
      ga(1)=1.5
      ga(2)=0.
      ga(3)=0.
      ro(1)=-.5
      ro(2)=0.
      ro(3)=0.
      write(6,110) ga(1),ro(1)
  110 format(10x,'adams bashfort=',4x,
     1       'gam=',f8.3,4x,'ro',f8.3)
   12 continue
      beta=dt/re*0.5                                                    
      call solve
      stop
      end
      subroutine solve
c
c     code for computation of two-dimensional incompressible flows
c     in general   coordinates.  with vorticity stream function.
c
c     this code solves flow fields with different possibilities
c     of boundary conditions in x1 and x2.
c     the time discretization is the third order runge kutta.
c     this code permits also to check the accuracy of the 
c     arakawa scheme to solve inviscid flows.
c
      include 'param.f'
      common/wrre/nwrit,nread,nfilr
      common/mesh/dx1,dx1q,dx2,dx2q
      dimension y(ndd,m1,m2)
      dimension vor(m1,m2),psi(m1,m2),dvor(m1,m2)
      dimension psc(m1,m2)
      common/dim/n1,n1m,n2,n2m
      common/d2/nstop,nprint,ntst,npin,nmolp,ntra,mnpi
      common/tstep/dt
      common/reyn/beta,re
      common/nonlt/rt(m1,m2),ru(m1,m2)
      common/tscoe/ga(3),ro(3),nsst
      common/camo/y1ma,y2ma,y1mi,y2mi,vori,ensti
      common/circul/vorip,vorim
      common/chnlc/chal,chbe,chga
      common/sliwal/insls,insln
      common/mgoui/ncm
      common/mgout/reml1,reml2,rrm
      common/ntpr/ntime,ns
      common/vpmima/vorma,vormi,pscma,pscmi
      common/icflp/icflm,cflma,dtl
      common/tfinst/tprin,tpin,tfin
      data istop,enen/0,0./
      pi=2.*asin(1.)
      npfile=10*npin
c
c     step  and messh sizes calculations
c
      call meshes
c
c    assigned the coordinate transformation by
c    algebraic expressions. The coordinate in general
c    is non-orthogonal
c
      call coordi(y)
      call metric(y)
c     call pmetri(y)
  758 format(1x,2i3)
      write(6,754)n1,n2
  754 format(10x,'n1=',i4,2x,'n2=',i4,3x)
      write(6,755) dx1,dx2,dt,ntst
  755 format(3x,'dx1=',e10.3,3x,'dx2=',e10.3,3x
     1      ,'dt=',e10.3,3x,'ntst=',i4)
  170 continue
      open(35,file='multgr.out')
      open(20,file='thisti.out')
      nti=0
      time=0.
      ntime=0
      ntii=0
      n1mh=n1m/2+1
c
c     initial conditions
c
      call phini
      if(nread.eq.0) call initia(vor,psi,y,dvor,psc)
      if(nread.eq.1) call inirea(vor,psi,time,ntii,psc)
      call enerca(vor,psi,enen0,y,psc)
      if(enen0.lt..1e-10) enen0=.1e-10
      call carmod(vor,psi,y)
c     write(20,158)ntime,time,enen0,vori,ensti,y1ma,y2ma,y1mi,y2mi
      write(6,158)ntime,time,enen0,vori,ensti,y1ma,y2ma,y1mi,y2mi
c
c    print the initial field
c
      call outpf(vor,psi,time,psc)
      nstea=nti
c
c  ********* starts the time dependent calculation ***
c
      if(icflm.eq.1) then
      nstop=100000   
      ntstf=100000   
      ntii=1
                     else
      nstop=nstop+ntii
      ntstf=ntii+ntst
      ntii=ntii+1
                     endif
      write(6,711)nprint,ntii,ntstf,nstop,nti,dt,chal,chbe,chga
  711 format(1x,'in cond',5i8,2x,e11.3,/3x,'chal=',e8.2,3x,'chbe=',e8.2,
     1       3x,'chga=',e8.2)
c     write(20,717)n1,n2,dt,chal,chbe,chga
  717 format(1x,'modops',3x,'n1xn2=',i3,'x',i3,3x,'dt=',e10.3
     1   ,3x,'chal=',e8.2,3x,'chbe=',e8.2,3x,'chga=',e8.2)
  181 continue
c
c    the time dipendent calculation starts
c
      call  incvel 
      do 350 ntime=ntii,ntstf
      if(ntime.gt.nstop) istop=1
c
c  ********* calculation of the vorticity and stream function
c
      tin=time
c
c   evaluation of quantities for the radiative b.c.
c
      call  cvnool(vor,psi,psc)
      call  cvweol(vor,psi,psc)
      do 2000 ns=1,nsst
      alp=(ga(ns)+ro(ns))
      gam=ga(ns)
      rho=ro(ns)
      tino=tin
      tin=tin+alp*dt
      call tschem(vor,psi,alp,gam,rho,y,dvor,tin,psc)
 2000 continue
c
c    radiative velocity
c
      call  cvelno(vor,psi,psc,y)
      call  cvelwe(vor,psi,psc,y)
      time=time+dt
      vomack=0.
      do i=1,n1
      do j=1,n2
      vomack=max(vomack,abs(vor(i,j)))
      enddo
      enddo
      if(vomack.gt..1e+04) go to 178
c
c     print time history of some global quantity 
c     
      if(icflm.eq.1) then

c   Courant number calculation 

      call cfl(psi,cflmm)
c     print*,'88888******   dt=',dt,'  cflmm',cflmm

      dt=cflma/cflmm
      cflm=cflma
      if(dt.gt.dtl) dt=dtl
      if(dt.lt..1e-05)      then
      write(6,368)dt
  368 format(3x,'end calcul for dt=',e12.4)
      stop
                            endif
                      endif
      if(icflm.eq.1) then
      if(amod(time,tpin).lt.dt) go to 306
      go to 305
                     else
      ntim=ntime-1
      if(mod(ntime,npin).eq.0) then
      call cfl(psi,cflmm)
      cflm=cflmm*dt
      go to 306       
                                endif
      go to 305
                     endif
  306 continue
      timepr=time
c
c
c     total energy and other global quantities
c
      call enerca(vor,psi,enen,y,psc)
c
c   evaluates some global quantities for the case of adipoleS
c
      call carmod(vor,psi,y)
  938 continue
c
c   write on file thist
c
      if(re.gt..1e06) then
      write(20,157)timepr,vorip,vorim,enen,ensti
     1             ,ncm,reml1,rrm,cflm
      write(6,158)ntime,timepr,vorip,vorim,enen,ensti
     1            ,vorma,pscma,vormi,pscmi,cflm
  158 format(2x,i4,1x,e10.3,2x,4(1x,e11.5),1x,5e10.3)
  157 format(2x,e11.4,1x,4(1x,e11.5),1x,i5,3x,2e11.5,e10.3)
                      else
      write(20,257)timepr,vori,enen,ensti
     1             ,ncm,reml1,rrm,cflm
      write(6,258)ntime,timepr,vori,enen,ensti
     1            ,vorma,pscma,vormi,pscmi,cflm
  258 format(2x,i4,1x,e10.3,2x,3(1x,e11.5),1x,5e10.3)
  257 format(2x,e11.4,1x,3(1x,e11.5),1x,i5,3x,2e11.5,e10.3)
                      endif
  305 continue
      if(istop.eq.1) go to 156
c
c    write the flow field
c
      if(icflm.eq.1) then
      if(amod(time,tprin).lt.dt) go to 301
            go to 300
                     else
      if(mod(ntime,nprint).eq.0) go to 301
            go to 300
                     endif
  301 continue
      call outpf(vor,psi,time,psc)
  300 continue
            if(time.gt.tfin.and.icflm.eq.1) then
      call outpf(vor,psi,time,psc)
      write(6,367)time
  367 format(3x,'end calcul for t=',e12.4)
      stop
                                            endif
  350 continue
  156 continue
      write(6,159)nstea
  159 format(/20x,'time convergency',3x,'nstea=',i4)
 1159 format(2x,i5,3e12.4)
      write(6,909)
  909 format(//20x,'steady solution')
      if(iexi.eq.1) stop
  180 continue
      iexi=1
      go to 167
  178 continue
      write(6,*)' program stops for too high vorticity'
      call outpf(vor,psi,time,psc)
      stop
  167 continue
      return
      end
c
c  **************  subrout tschem
c in this routine the vorticity stream function equations are
c solved by evaluating before the passive scalar and then
c  the vorticity field.
c  The streamfunction is calculated by a multigrid
c  method
c
      subroutine tschem(vor,psi,al,ga,ro,y,dvor,time,psc)
      include 'param.f'
      dimension dvor(m1,m2),vor(m1,m2),psi(m1,m2),y(ndd,m1,m2)
      dimension psc(m1,m2)
      common/dim/n1,n1m,n2,n2m
      common/reyn/beta,re
      common/schm/schmid
      common/mgoui/ncm
      common/mgout/reml1,reml2,rrm
      common/ntpr/ntime,ns
      common/dvowal/dvorbs(m1),dvorbn(m1)
      common/psbou/psbi(2,m1),psbj(2,m2)
      common/dpswal/dpscbs(m1),dpscbn(m1)
      common/bdqout/dvest(m2),dpest(m2),dpsest(m2)
      common/bdqinl/dvwes(m2),dpwes(m2),dpswes(m2)
      common/vopsin/psinf(m2),vobj(m2)
c
c   solve the passive scalar convected by the velocity field
c   non-linear terms for psc
c
      call hdnl(psc,psi,dvor)
c  *****  solve the dvor=vor(n+1-vor(n) potential vorticity
      visc=1./(schmid*re)
c
c   inlet stream  passing from zero to one smootly
c
      call extflo(time,ft)
c
c    boundary conditions
c
      call bcwest(ga,ro,vor,psi,y,dvor,psc,ft)
      call bcest(ga,ro,vor,psi,psc)
      call bcnsvi(ga,ro,vor,psi,y,dvor,psc,ft)
      if(re.lt.1.e+05) then
c
c   viscous case
c
      call invtps(dvor,al,ga,ro,psc,visc)
                       else
c
c   inviscid case
c
      call invtpi(dvor,al,ga,ro,psc)
                       endif
c
c   time integration  adams baschfort or third ord runge kutta nonlinear
c    vorticity Jacobian
c
      call hdnl(vor,psi,dvor)
c  *****  solve the dvor=vor(n+1-vor(n) potential vorticity
      if(re.lt.1.e+05) then
      visc=1./re
c
c    the terms of the cross derivatives are added to the 
c    non-linera termas by Arakawa
c    if the coorsdinates are orthogonal teh cross-derivatives are
c                   zero
c
      call hdnlvi(vor,dvor,visc)
      call invtrv(dvor,al,ga,ro,vor,visc)
                       else
c
c   vorticity in the inviscid case
c
      call invtri(dvor,al,ga,ro,vor)
                       endif
c
c  ********* calculation of the stream function
c
      call phcalc(vor,psi)
      return
      end
c
c  **************  subrout extflo
c    time variation from t=0 up to t=tau
c
      subroutine extflo(tin,ft)
      include 'param.f'
      common/timinf/tau1
      if(tin.le.tau1) then
      tina=tin/tau1
      ft=3.*tina**2-2.*tina**3
                      else
      ft=1.
                      endif
      return
      end

