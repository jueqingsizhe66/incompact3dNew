      subroutine solve
c     code for computation of two-dimensional incompressible flows
c     in cartesian coordinates.  with vorticity stream function.
c
c     this code solves flow fields with different boundary cond. in x1 
c     and in x2 as described in Chapter 4.
c     this code uses Arakawa schemes and
c     the time discretization is the third order runge kutta.
c     or Adams Bashfort
c
c    the  initial  condition are in inicha    ****************
c
c    the  beta effect  is considered
c  ***********************************************************
c
c    Poisson eq is solved by the phcal subr. that uses FFT routines
c
c  *********************************************************+
      include 'param.f'
      dimension vor(m1,m2),psi(m1,m2),vorq(m1,m2)
      dimension dvor(m1,m2)
      dimension ru(m1,m2),rt(mpsc,m1,m2),psc(mpsc,m1,m2)
      common/wrre/nread
      common/mesh/dx1,dx1q,dx2,dx2q
      common/metri/g2m(m1),g2c(m1)
      common/dim/n1,n1m,n2,n2m
      common/d1/alx1i,alx1f,alx2i,alx2f
      common/d2/nfield,ntst,nhist,nscrn
      common/tstep/dt
      common/pi/pi
      common/waves/an(m1)
      common/tscoe/ga(3),ro(3),nsst
      common/parmod/yc1mo,yc2mo,akmo,ramo,velmo,vsi
      common/camoc/vocmax,y1cma,y2cma,vocmin,y1cmi,y2cmi
      common/psmoc/pscmax,y1pma,y2pma,pscmin,y1pmi,y2pmi
      common/sfmoc/sfcmax,y1sma,y2sma,sfcmin,y1smi,y2smi
      common/chnlc/chal,chbe,chga
      common/visct/re,ekmn
      common/qhoold/vohoe(m2),pshoe(m2),vohow(m2),pshow(m2)
     1             ,sshoe(mpsc,m2),sshow(mpsc,m2)
      common/qveold/voven(m1),psven(m1),voves(m1),psves(m1)
     1             ,ssven(mpsc,m1),ssves(mpsc,m1)
      common/inmodo/inmod
      common/camo/cflm
      common/icflp/icflm,cflma,dtl
      common/tfinst/tprin,tpin,tfin
      common/tpsta/tpstar
      common/iperi/ib2per
      common/indpe2/n2i,n2f
      common/ftooo/fto
c
      pi=acos(-1.)
c
c     step  and messh sizes calculations
c
      call meshes
  758 format(1x,2i3)
      write(6,754)n1,n2,str1
c     write(6,754)n1,n2
  754 format(10x,'n1=',i3,2x,'n2=',i3,3x,'y1 coordinate with str1='
     1       ,e8.3)
      write(6,755) dx1,dx2,dt,ntst,re
  755 format(3x,'dx1=',e10.3,3x,'dx2=',e10.3,3x
     1      ,'dt=',e10.3,3x,'ntst=',i4,3x,'re=',e10.3)
      call indic
      call coordi
      call maketop
      time=0.
      ntime=0
      n1mh=n1m/2+1
c
c   routine to evaluate quantities necessary for FFT
c
      call phini
c
c     sets initial outflo velocity = 0
c
      call incveve
      call inchove
c
c
c     initial conditions  generated(initia   nread=0) 
c           or read from file(inirea   nread=1)
c
      if(nread.eq.0) then
      call initia(vor,psi,ru,rt,psc)
      x1mao=0.
      x2mao=0.
      call cfield(vor,psi,psc)
      call outpf(vor,psi,time,vorq,psc)
                     else
      timei=nscrn
      call inirea(vor,psi,timei,vorq,psc,ru,rt)
      call cfield(vor,psi,psc)
      time=timei
                     endif

      x1mao=0.
      x2mao=0.
c     write(6,815) vorip,vma,ensti,enerv
  815 format(3x,'in cir,vma,enst,ene',2x,4e12.4)
c
c
c  ********* starts the time dependent calculation ***
c
      write(6,711)nfield,ntst,dt,chal,chbe,chga
  711 format(1x,'in cond',2i8,2x,e11.3,/3x,'chal=',e8.2,3x,'chbe=',e8.2,
     1       3x,'chga=',e8.2)
  181 continue
      cflm=0.
      tin=0.
  717 format(1x,2i4,3x,3e15.4)
      call outth(ntime,time,vor,psi)
      if(icflm.eq.1) then
      ntst=1000000
      write(6,*)' calculation at dt variable by fixing cfl=',cflma
                     endif
      fto=0.
      call potvor(vor,vorq)
c
c    the dependent calculation starts
c
      do 350 ntime=1,ntst
c
c******    stores old values 
c
      call choold(vorq,psi,psc) 
      call cveold(vorq,psi,psc) 
c
c  ********* calculation of the vorticity and stream function
c
      do 2000 ns=1,nsst
      alp=(ga(ns)+ro(ns))
      gam=ga(ns)
      rho=ro(ns)
      tin=tin+alp*dt
      call tschem(vor,psi,ru,alp,gam,rho,dvor,psc,rt,tin,vorq)
 2000 continue
      call potvor(vor,vorq)
      call cvevel(vorq,psi,psc)
      call chovel(vorq,psi,psc)
      time=time+dt
c
c     print some quantity max(divg),total energy,max(veloc),
c     and their gradient in time.
c
      if(icflm.eq.1) then
C
c   courant number calculation 
C
      call cfl(psi,cflmm)
c 
C
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
      if(mod(ntime,nhist).eq.0) then
      call cfl(psi,cflmm)
      cflm=cflmm*dt
      go to 306       
                                endif
      go to 305
                     endif
  306 continue
c
c     characteristics of the flowfield
c
      call cfield(vor,psi,psc)
      if(vocmax.gt..1e04) go to 451    
      call outth(ntime,time,vor,psi)
  158 format(1x,i4,1x,e10.3,1x,10(1x,e10.4))
  305 continue
c
c    write the flow field 
c
      if(time.ge.tpstar) then
      if(icflm.eq.1) then
      if(amod(time,tprin).lt.dt) go to 301
            go to 300
                     else
      if(mod(ntime,nfield).eq.0) go to 301
            go to 300
                     endif
  301 continue
      call outpf(vor,psi,time,vorq,psc)
  300 continue
                     endif
            if(time.gt.tfin.and.icflm.eq.1) then
      call outpf(vor,psi,time,vorq,psc)
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
      go to 452
  451 continue
      write(6,109) vma
  109 format(//20x,'vmax = ',e12.4)
      call outpf(vor,psi,time,vorq,psc)
  452 continue
      return
      end
c
c    time perturbations at the inlet given by the following 
c    routines
c
c  **************  subrout extflo
c
      subroutine extflo(tin,ft)
      include 'param.f'
      common/timinf/tau1,ampltp,timepe,dtau,time0
      timefi=timepe+2.*tau1+3.*dtau
      time1=timepe+dtau
      time2=timepe+dtau+tau1
      time3=timepe+2.*dtau+tau1
      time4=timepe+2.*(dtau+tau1)
      if(time0.gt.0.) then
        if(tin.le.time0) then
          tina=tin/time0
          fti=3.*tina**2-2.*tina**3
        else
          fti=1.
        endif
      else
        fti=0.
      endif
      if(tin.gt.time1.and.tin.le.time2) then
      ftp=1.
                                        endif
      if(tin.gt.time3.and.tin.le.time4) then
      ftp=-1.
                                        endif
      if(tin.le.timepe.or.tin.ge.timefi) then
      ftp=0.
                                        endif
      if(tin.ge.timepe.and.tin.le.time1) then
      tina=(tin-timepe)/dtau
      ftp=3.*tina**2-2.*tina**3
                       endif
      if(tin.gt.time2.and.tin.le.time3) then
      tina=(tin-time2)/dtau
      ftp=1.-2.*(3.*tina**2-2.*tina**3)
                                        endif
      if(tin.gt.time4.and.tin.le.timefi) then
      tina=(tin-time4)/dtau
      ftp=-1.+(3.*tina**2-2.*tina**3)
                                          endif
      ft=fti+ftp*ampltp
      if(tin.le.timefi) then
c     write(6,*) ' TIN = ',tin,' FT = ',ft
c     write(91,131)tin,ft
  131 format(3x,'time perturb',2e12.4)
                                          endif
      return
      end
c 
c             subroutine extflo1
c
c
      subroutine extflo1(tin,ft)
      include 'param.f'
      common/timinf/tau1,ampltp,timepe,dtau,time0
      common/timinf1/om,eps
      time1=timepe+dtau
      if(tin.ge.timepe.and.tin.le.time1) then
      tina=(tin-timepe)/dtau
      ftp=1+(eps*sin(om*tin))*(3.*tina**2-2.*tina**3)
                      endif
      if(tin.gt.time1) then
      ftp=1+(eps*sin(om*tin))
                      endif
  131 format(3x,'time perturb',2e12.4)
      return
      end
c
c
c             subroutine extflo2
c
      subroutine extflo2(tin,ft)
      include 'param.f'
      common/pi/pi
      common/timinf/tau1,ampltp,timepe,dtau,time0
      common/timinf1/om,eps
      pi=2.*asin(1.)
      ftp=1+(eps*sin(om*2*pi*tin))
      ft=ftp
  131 format(3x,'time perturb',2e12.4)
      return
      end
c 
c
c  **************  subrout tschem
c
      subroutine tschem(vor,psi,ru,al,ga,ro,dvor,psc,rt,time,vorq)
      include 'param.f'
      dimension dvor(m1,m2),vor(m1,m2),psi(m1,m2),vorq(m1,m2)
      common/dsebns/d2psud(m1),d2pnor(m1)
      common/bchodv/dvsud(m1),dvnor(m1)
      common/bchodp/dpsud(m1),dpnor(m1)
      common/bchods/dssud(mpsc,m1),dsnor(mpsc,m1)
      common/bcvedv/dvest(m2),dvwes(m2)
      common/bcvedp/dpest(m2),dpwes(m2)
      common/bcveds/dsest(mpsc,m2),dswes(mpsc,m2)
      common/dim/n1,n1m,n2,n2m
      dimension ru(m1,m2),dpsi(m1,m2),rt(mpsc,m1,m2),psc(mpsc,m1,m2)
      dimension pscap(m1,m2)
      common/visct/re,ekmn
      common/pscqu/sch(mpsc),scla,npscf
      common/ydif/ydi(m2)
      common/quainf/qinp(6,md2)
      common/ipert/ipert
      common/iperi/ib2per
      common/indpe2/n2i,n2f
      dimension psio(m2)
      common/indwa/jmv(m2),jpv(m2)
      common/cwalv/inbcvs,inbcvn,inbcvw,inbcve
      common/inflow/voinf(m2),psinf(m2),ssinf(mpsc,m2)
      common/ftooo/fto
      common/inmodo/inmod
      if(ipert.eq.1) then
      call extflo(time,ft)
      endif
      if(ipert.eq.2) then
      call extflo1(time,ft)
      endif
      if(ipert.eq.3) then
      call extflo2(time,ft)
      endif
      dft=ft-fto
      fto=ft
c
c   ******  vor, psi ,psc values on the walls j=1 j=n2  *************
c
      call potvor(vor,vorq)
C     write(6,*) ' INMOD ',inmod
C     write(6,*) ' TIME ',time,' TIMEO ',timeo
      if(inmod.eq.-5) then
      call bcspaml(time,timeo)
      timeo = time
      end if
      call bdqpsn(ga,ro,vorq,psi,psc,dft)
c
c   ******  vor, psi ,psc values on inf outf  i=1 i=n1  *************
c
      call bdqpew(ga,ro,vorq,psi,psc,dft)           
      if(inbcvw.eq.-2) then
      do j=1,n2
        qinp(1,j)=dvwes(j)
        qinp(2,j)=dswes(1,j)
        qinp(3,j)=dpwes(j)
      enddo
                       endif
      if(inbcvs.eq.-2) then
      do i=1,n1
        qinp(1,i)=dvsud(i)
        qinp(2,i)=dssud(1,i)
        qinp(3,i)=dpsud(i)
      enddo
                       endif
c
      n1mh=n1m/2+1
c     do i=1,n1m
c       write(90,*)i,vor(i,1),dvsud(i),psi(i,1),dpsud(i)
c     end do
c     close(90)
      do npsc=1,npscf !{ the do for npsc starts
        visc=1./(re*sch(npsc))
        do i=1,n1
          do j=1,n2f
            pscap(i,j)=psc(npsc,i,j)
          enddo
        enddo
c
c   ****** passive scalar calculation   *********
c
c   time integration  adams bashfort or third ord runge kutta nonlinear
c
      call hdnl(pscap,psi,dvor)

c
c  moves q to vor to compute the viscous terms and add to the non-linear
c  terms in subroutine invtps
c
c
c  *****  solve the dqpsc=psc(n+1-psc(n) passive scalar equation
c
c
c 
      if(ib2per.eq.1) then
        call invpps(dvor,al,ga,ro,pscap,rt,visc,npsc)    
      else
        call invtps(dvor,al,ga,ro,pscap,rt,visc,npsc)    
      endif
c
c  calculation of the vorticity
c
        do i=1,n1
          do j=1,n2f
            psc(npsc,i,j)=dvor(i,j)+pscap(i,j)
          enddo
        enddo
      enddo    !} end of do for npsc
      do i=1,n1
        do j=1,n2f
          dvor(i,j)=0.
        enddo
      enddo
      visc=1./re
c
c   ****** vorticity calculation   *********
c
c   time integration  adams baschfort or third ord runge kutta nonlinear
c
      call hdnl(vorq,psi,dvor)
c
c
c  moves q to vor to compute the viscous terms and add to the non-linear
c  terms in subroutine invtv
c
c
c  *****  solve the dqvor=vor(n+1-vor(n) momentum equation
c
c
c 
      if(ib2per.eq.1) then
      call invpv(dvor,al,ga,ro,vor,ru,visc,ekmn)    
                      else
      call invtv(dvor,al,ga,ro,vor,ru,visc,ekmn)    
                      endif
c
c  calculation of the vorticity
c
      do i=1,n1
        do j=1,n2f
          vor(i,j)=dvor(i,j)+vor(i,j)
        enddo
      enddo
c     write(6,130)time,ft,fto,dft,voinf(2),dvwes(2)
c    1           ,dvor(1,2),vor(1,2)
  130 format(4e10.2,4x,3e12.4,2x,3e10.3)
c
c  add right hand side for streamfunction calv.
c
      if(ib2per.eq.0) then
        call dsepns 
        do j=2,n2m
          do i=2,n1m
            dvor(i,j)=dvor(i,j)+
     %                ydi(j)*(d2pnor(i)-d2psud(i))+d2psud(i) 
          end do
        end do
      endif
c
c  boundary conditions for quantity dphi which has
c  identically zero b.c. at the upper and lower level
c
c
c  ********* calculation of the stream function 
c            tridiag in vertical x2
c
      if(ib2per.eq.0) then
      call bdphew
      call phcalc(dvor,dpsi)
                      else
      call bdppew(psi)
      call phcalp(vor,psi)
                      endif
c
c  moves vor to q to advance in time
c
      if(ib2per.eq.0) then
  181 format(11e11.3)
c       write(6,181)time,dpsud(1),dpnor(1),psi(1,1),psi(1,n2),
c    1      ydi(1),ydi(n2),dpsi(1,1),dpsi(1,n2)
        do i=1,n1
          do j=1,n2
            psi(i,j)=psi(i,j)+dpsi(i,j)+
     %               ydi(j)*(dpnor(i)-dpsud(i))+dpsud(i)
          enddo
        enddo
      endif
c
c******   update vor(boundary)=vor(boundary)+dvor(boundary)  
c
      call boucqt
      if(inbcvs.eq.-2) then
      do i=1,n1
      qinp(4,i)=vor(i,1)
      qinp(5,i)=psc(1,i,1)
      qinp(6,i)=psi(i,1)
      enddo
                        endif
      if(inbcvw.eq.-2) then
      do j=1,n2f
      qinp(4,j)=vor(1,j)
      qinp(5,j)=psc(1,1,j)
      qinp(6,j)=psi(1,j)
      enddo
                        endif
      vomal=0.
      psmal=0.
      do i=1,n1
      do j=1,n2
      psmal=max(abs(psi(i,j)),psmal) 
      vomal=max(abs(vor(i,j)),vomal)
      enddo
      enddo
      write(99,*)' max quant',psmal,vomal
      return
      end
c
c  **********************  set initial outflo velocity = 0  *********
c    for radiative boundary conditions
c
      subroutine incveve
      include 'param.f'
      parameter(m1m=m1-1,m2m=m2-1)
      common/dim/n1,n1m,n2,n2m
      common/velcve/cven(2+mpsc,m1),cves(2+mpsc,m1)
      common/pscqu/sch(mpsc),scla,npscf
c
c   set initial radiative b.c. =0
c
      lfi=2+npscf
      do i=1,n1
        do l=1,lfi
       cven(l,i)=0.
       cves(l,i)=0.
        enddo
      enddo
      return
      end                         
c
c  **********************  set initial outflo velocity = 0  *********
c    for radiative boundary conditions
c
      subroutine inchove
      include 'param.f'
      parameter(m1m=m1-1,m2m=m2-1)
      common/dim/n1,n1m,n2,n2m
      common/velcho/choe(2+mpsc,m2),chow(2+mpsc,m2)
      common/pscqu/sch(mpsc),scla,npscf
c
c   set initial radiative b.c. =0
c
      lfi=2+npscf
      do j=1,n2
        do l=1,lfi
       choe(l,j)=0.
       chow(l,j)=0.
        enddo
      enddo
      return
      end                         
