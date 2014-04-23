      program main
      include 'param.f'
      common/d1/alx1i,alx1f,alx2i,alx2f
      common/strpar/str1,xcra,etra,istr
      common/d2/nfield,ntst,nhist,nscrn
      common/dim/n1,n1m,n2,n2m
      common/wrre/nread
      common/tscoe/ga(3),ro(3),nsst
      common/chnlc/chal,chbe,chga
      common/tstep/dt
      common/visct/re,ekmn
      common/parmod/yc1mo,yc2mo,akmo,ramo,velmo,vsi
      common/inmodo/inmod
      common/angmd/thet0
      common/cwalv/inbcvs,inbcvn,inbcvw,inbcve
      common/cwalp/inbcps,inbcpn,inbcpw,inbcpe
      common/cwals/inbcss,inbcsn,inbcsw,inbcse
      common/pscqu/sch(mpsc),scla,npscf
      common/anmod/ar1,ar2
      common/betef/beta
      common/pi/pi
      common/nmolt/nmolp
      common/topoda/dim1i,dim1f,dim2i,dim2f,hmax,dhmax
      common/itop/itopog
      common/parinf/u0,blve,cpp 
      common/ipert/ipert
      common/timinf/tau1,ampltp,timepe,dtau,time0
      common/timinf1/om,eps
      common/icflp/icflm,cflma,dtl
      common/tfinst/tprin,tpin,tfin
      common/tpsta/tpstar
      common/iperi/ib2per
      common/indpe2/n2i,n2f
      common/wavnu/akx,nvor
      common/forci/forc(m1,m2),fff
      common/pamlsp/ampl1,ampl2,omeg1,omeg2
      open(20,file='thist.out')
      open(33,file='traj1.out')
      open(34,file='traj2.out')
      open(35,file='quavt.out')
      open(36,file='quavco.out')
      open(37,file='quapsc.out')
      open(38,file='enesps.out')
c
c
      open(15,file='psomgbc.d')
      read(15,*) n1,n2,nsst
      read(15,*) nread
      read(15,*) nmolp
      read(15,*) dt,ntst,nfield,nhist,nscrn
      read(15,*) icflm,cflma,dtl
      read(15,*) tprin,tpin,tfin,tpstar
      read(15,*) alx1i,alx1f,str1,xcra,etra,istr
      read(15,*) alx2i,alx2f
      read(15,*) re,beta,ekmn
      read(15,*) chal,chbe,chga
      read(15,*) inbcvs,inbcvn,inbcvw,inbcve
      read(15,*) inbcps,inbcpn,inbcpw,inbcpe
      read(15,*) inbcss,inbcsn,inbcsw,inbcse
      read(15,*) ib2per
      read(15,*) npscf,scla
      read(15,*) (sch(npsc),npsc=1,npscf)
      read(15,*) inmod,yc1mo,yc2mo,ramo,velmo,vsi
      read(15,*) thet0
      read(15,*) ar1,ar2
      read(15,*) ipert
      read(15,*) itopog,dim1i,dim1f,dim2i,dim2f,hmax,dhmax
      read(15,*) u0,blve,cpp
      read(15,*) tau1,ampltp,timepe,dtau,time0
      read(15,*) om,eps
      read(15,*) nvor,akx,alx2d
c
      if(inmod.eq.10) then
      fff=thet0
      write(6,*)'forcing with fff=',fff
                      endif
      akmo=3.83711/ramo
      pi=2.*asin(1.)
      if(inmod.eq.-4) then
      alx2=2.*pi*alx2d/akx*nvor
      alx2i=-alx2*0.5
      alx2f=+alx2*0.5
                      endif
      if(inmod.eq.-5) then
c
c   spatial mixing layer 
c   amplitude and frequency of inlet perturbations
c
      ampl1=yc1mo
      ampl2=yc2mo
      omeg1=ramo
      omeg2=velmo
                      endif
      thet0=thet0*pi/180.
      if(inmod.eq.1) write(6,201)thet0
  201 format(3x,'modone vort stream funct. tht0=',e10.4)
      n1m=n1-1
      n2m=n2-1
c
c    ib2per=1 allows periodicity in x2
c
      if(ib2per.eq.1) then
      n2i=1
      n2f=n2m
                      else
      n2i=2
      n2f=n2
                       endif
      write(6,*)'     ib2per,n2i,n2f',ib2per,n2i,n2f
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
      if(inmod.eq.-3) then
      open(51,file='circp.out')
      open(52,file='circm.out')
      open(53,file='troma.out')
      open(54,file='tromi.out')
      open(55,file='arep1.out')
      open(56,file='arep2.out')
      open(57,file='trpma.out')
      open(58,file='floin.out')
      open(59,file='floou.out')

                      endif
      if(inmod.eq.-4) then
      open(51,file='deltm.out')
      open(52,file='omdbr.out')
      open(53,file='omdce.out')
      open(54,file='smdbr.out')
      open(55,file='smdce.out')
                      endif
      call solve
c
      close(15)
      close(20)
      close(12)
      close(13)
c     close(30)
      stop
      end
c
c  ****************************** subrout initia **********************
c
c   initial conditions  of vortity
c
      subroutine initia(vor,psi,ru,rt,psc)
      include 'param.f'
      dimension vor(m1,m2),psi(m1,m2),psc(mpsc,m1,m2)
      dimension ru(m1,m2),rt(mpsc,m1,m2)
      common/dim/n1,n1m,n2,n2m
      common/indbo/imv(m1),ipv(m1)
      common/mesh/dx1,dx1q,dx2,dx2q
      common/parmod/yc1mo,yc2mo,akmo,ramo,velmo,vsi
      common/inmodo/inmod
      common/pi/pi
      common/coor/yp1(m1),yp2(m2)
      common/metri/g2m(m1),g2c(m1)
      common/pscqu/sch(mpsc),scla,npscf
      common/inflow/voinf(m2),psinf(m2),ssinf(mpsc,m2)
      common/infmxl/voins(m1),psins(m1),ssins(mpsc,m1)
      common/infmxp/voinp(m1),psifr(m1),
     %              psifi(m1),psiar(m1),psiai(m1),psinp(m1)
      dimension uinfl(m2),uins(m1),uinpe(m1),uinc(m1)
      common/timinf/tau1,ampltp,timepe,dtau,time0
      common/parinf/u0,blve,cpp
      common/iperi/ib2per
      common/indpe2/n2i,n2f
      common/forci/forc(m1,m2),fff
      dimension u(m1)
c
c if  inmod=1 Lamb dipole
c
c
      if(inmod.eq.1) then
      write(6,*)'  Lamb dipole  y1c,y2c,v,r,k,vs'
      write(6,106)yc1mo,yc2mo,velmo,ramo,akmo,vsi
  106 format(2x,6e10.3)
      call dipol(vor,ramo,vsi,yc1mo,yc2mo,akmo,velmo,psc)
c
                     endif
c
c
c if  inmod=0 Stern modon
c
      if(inmod.eq.0) then
      write(6,*)'  Stern modon   y1c,y2c,v,r,vs'
      write(6,106)yc1mo,yc2mo,velmo,ramo,vsi
      call modon(vor,ramo,vsi,yc1mo,yc2mo,velmo,psc)
                     endif
c
c if  inmod=-1 tripole
c
c
      if(inmod.eq.-1) then
      write(6,*)'  TRipol  y1c,y2c,crad,alpha'
      call tripol(vor,yc1mo,yc2mo,ramo,velmo,vsi,ru,psc)
                     endif
c
c if  inmod=-2 monopole
c
      if(inmod.eq.-2) then
      write(6,*)' MONOPOL    y1c,y2c,sig,amp '
      write(6,106)yc1mo,yc2mo,ramo,velmo
      call monpol(vor,yc1mo,yc2mo,ramo,velmo,psc,vsi)
                     endif
c
c   gets psi from vorticity inverting Poisson eq.
c   in  this case with b.c. zero
c
c
      if(inmod.gt.-3) then
      if(ib2per.eq.0) then
      call phcalc(vor,psi)
                      endif
      if(ib2per.eq.1) then
      call phcalp(vor,psi)                                               
                      endif
      write(6,*) ' psi evaluated in initia'
                     endif
c
c   inlet conditions at the west boundary
c
c
c if  inmod=-3 inlet conditions
c
      if(inmod.eq.-3) then
      write(6,*)' INLET VELOCITY    u0,blve ',u0,blve
      write(6,*)' time var tau1,ampltp,timepe,dtau,time0',
     1          tau1,ampltp,timepe,dtau,time0
      call inltop(u0,blve,vor,psi,psc)
                      endif
      if(inmod.eq.-7) then
      write(6,*)' INLET VELOCITY    u0,blve ',u0,blve
      write(6,*)' time var tau1,ampltp,timepe,dtau,time0',
     1          tau1,ampltp,timepe,dtau,time0
      call inljet(u0,blve,cpp,vor,psi,psc)
                      endif
      if(inmod.eq.-8) then
      write(6,*)' INLET VELOCITY    u0,blve ',u0,blve
      write(6,*)' time var tau1,ampltp,timepe,dtau,time0',
     1          tau1,ampltp,timepe,dtau,time0
      call inljetsud(u0,blve,cpp,vor,psi,psc)
                      endif
      if(inmod.eq.-4) then
      if(ib2per.eq.0) then 
      write(6,*)' mixlay time developing calculation'
      write(6,*)' ampl fourier mode yc1mo ,yc2mo',yc1mo,yc2mo 
      write(6,*)' stops ib2per must be =1 '
      stop
                      endif
      write(6,*)' mixlay time developing     '
      write(6,*)' mixlay time developing vor calc enter psi    '
      call inimlt(vor,psi,psc,yc1mo,yc2mo)
                      endif
      if(inmod.eq.-5) then
      if(ib2per.eq.1) then 
      write(6,*)' mixlay space developing calculation'
      write(6,*)' stops ib2per must be =0 '
      stop
                      endif
      write(6,*)' mixlay space developing '
      call inlmls(vor,psi,psc,u0,blve)
c
c   these commented instructions permit to perform
c   simulations of space developing mixing layers
c   with the vorticity at the inlet in the
c   whole domain
c
c     do j=1,n2
c       do i=1,n1
c         vor(i,j)= voins(i)
c         psi(i,j)= psins(i)
c         do npsc=1,npscf
c           psc(npsc,i,j)= ssinf(npsc,i)
c         enddo
c       enddo
c     enddo
                      endif
      if(inmod.eq.-6) then
c
c   this flow has not been tested carefully
c
      if(ib2per.eq.1) then 
      write(6,*)' planar jet space developing calculation'
      write(6,*)' stops ib2per must be =0 '
      stop
                      endif
      write(6,*)' planar jet space developing '
      call injets(vor,psi,psc,u0,blve)
                      endif
      do j=n2i,n2f
      do i=1,n1
      forc(i,j)=0.
      enddo
      enddo
      if(inmod.eq.10) then
c
c   this flow has not been tested carefully
c
      pi=2.*asin(1.)
      write(6,*)'in initia fff=',fff,yp2(1),yp2(n2)
      write(59,*)'in initia fff=',fff,yp2(1),yp2(n2)
      forma=0.
      do j=1,n2
      write(59,*)j,yp2(j),forcy
      do i=1,n1
      forcy=-fff*sin(pi*yp2(j))*sin(pi*yp1(i))
      vor(i,j)=0.
      psi(i,j)=0.
      forc(i,j)=forcy
      forma=max(abs(forc(i,j)),forma)
      do npsc=1,npscf
      psc(npsc,i,j)=0.5
      enddo
      enddo
      enddo
      write(6,*)' geostrophic with forcing forma=',forma
                     endif
c
c non-linear terms set equal to zero at t=0
c
      do j=n2i,n2f
      do i=1,n1
      ru(i,j)=0.
      enddo
      enddo
      do npsc=1,npscf
      do j=n2i,n2f
      do i=1,n1
      rt(npsc,i,j)=0.
      enddo
      enddo
      enddo
      return
      end
c
c  ****************************** subrout inltop  **********************
c
c   initial vort psi and psc for the inlet prof. over a topography
c   flow described in Sect. 4.12.2e
c
      subroutine inltop(u0,blve,vor,psi,psc)
      include 'param.f'
      dimension vor(m1,m2),psi(m1,m2),psc(mpsc,m1,m2)
      common/dim/n1,n1m,n2,n2m
      common/mesh/dx1,dx1q,dx2,dx2q
      common/coor/yp1(m1),yp2(m2)
      common/pi/pi
      common/indpe2/n2i,n2f
      common/pscqu/sch(mpsc),scla,npscf
      common/inflow/voinf(m2),psinf(m2),ssinf(mpsc,m2)
      dimension uinfl(m2)
      common/timinf/tau1,ampltp,timepe,dtau,time0
      common/cwalv/inbcvs,inbcvn,inbcvw,inbcve
      common/cwalp/inbcps,inbcpn,inbcpw,inbcpe
      common/cwals/inbcss,inbcsn,inbcsw,inbcse
      common/topoda/dim1i,dim1f,dim2i,dim2f,hmax,dhmax
      common/pscste/psstep
      if(inbcvn.eq.0.or.inbcvn.eq.2) then
      do j=1,n2m
      ym=yp2(n2)-(yp2(j+1)+yp2(j))*0.5
      uinfl(j)=u0*exp(-ym/blve)
      enddo
      do j=1,n2f
      ym=yp2(n2)-(yp2(j+1)+yp2(j))*0.5
      ssinf(1,j)=exp(-ym/blve)
      enddo
                     endif
      if(inbcvs.eq.0.or.inbcvs.eq.2) then
      do j=1,n2m
      ym=(yp2(j+1)+yp2(j))*0.5
      uinfl(j)=u0*exp(-ym/blve)
      enddo
      do j=1,n2f
      ssinf(1,j)=exp(-yp2(j)/blve)
      enddo
                     endif
      if(inbcvs.eq.-1) then
      do j=1,n2m
      ym=(yp2(j+1)+yp2(j))*0.5
      ydim=ym/blve
      if(ydim.ge.1.) ydim=1.
      uinfl(j)=u0*ydim**(1./7.)
      enddo
      do j=n2i,n2f
      ssinf(1,j)=uinfl(j)
      enddo
                     endif
      if(inbcvs.eq.0.or.inbcvs.eq.2) then
      psinf(1)=0.
      do j=1,n2m
      psinf(j+1)=psinf(j)+uinfl(j)/dx2
      enddo
      ybou=yp2(1)
                     endif
      if(inbcvs.eq.-1.or.inbcvs.eq.-2) then
      psinf(1)=0.
      do j=1,n2m
      psinf(j+1)=psinf(j)+uinfl(j)/dx2
      enddo
      ybou=yp2(1)
                     endif
      if(inbcvn.eq.0.or.inbcvn.eq.2) then
      psinf(n2)=0.
      do j=n2m,1,-1
      psinf(j)=psinf(j+1)-uinfl(j)/dx2
      enddo
      ybou=yp2(n2)
                     endif
      do j=2,n2m
      voinf(j)=-(psinf(j+1)-2.*psinf(j)+psinf(j-1))*dx2q
      enddo
      voinf(1)=voinf(2)
      voinf(n2)=voinf(n2m)
      ssima=0.
      psima=0.
      velma=0.
      vorma=0.
      open(47,file='infprof.out')
      write(47,*)'inbcvs,inbcvn  ',inbcvs,inbcvn
      do j=1,n2f
c     uinfa=u0*exp(-(yp2(j)-ybou)/blve)
      write(47,133) yp2(j),uinfl(j),psinf(j),voinf(j),
     1             (ssinf(n,j),n=1,npscf)
      vorma=max(vorma,abs(voinf(j)))
      velma=max(velma,abs(uinfl(j)))
      psima=max(psima,abs(psinf(j)))
      ssima=max(ssima,abs(ssinf(1,j)))
      enddo
      close(47)
  133 format(3x,8e12.4)
      write(6,*)' maximum inlet qua   ',velma,psima,vorma,ssima,mpsc,
     1                  npscf
      do i=1,n1
      do j=1,n2f
      psi(i,j)=0.
      vor(i,j)=0.
      enddo
      enddo
      do i=1,n1
      do j=1,n2f
      psc(1,i,j)=0.
      enddo
      enddo
      psstep=(dim1i-yp1(1))/(yp1(n1)-yp1(1))
      do j=1,n2f
      do i=1,n1
      x1nd=(yp1(i)-yp1(1))/(yp1(n1)-yp1(1))
      psc(2,i,j)=x1nd
      enddo
      enddo
      do j=n2i,n2f
      ssinf(2,j)=psc(2,1,j)
      enddo
      write(6,*)' values at the step of second scalar=',psstep
      return
      end
c
c  ****************************** subrout inljet  **********************
c
c   initial vort psi and psc for the inlet parabolic prof. over a topography
c
c
      subroutine inljet(u0,blve,cpp,vor,psi,psc)
      include 'param.f'
      dimension vor(m1,m2),psi(m1,m2),psc(mpsc,m1,m2)
      common/dim/n1,n1m,n2,n2m
      common/mesh/dx1,dx1q,dx2,dx2q
      common/coor/yp1(m1),yp2(m2)
      common/pi/pi
      common/indpe2/n2i,n2f
      common/pscqu/sch(mpsc),scla,npscf
      common/inflow/voinf(m2),psinf(m2),ssinf(mpsc,m2)
      dimension uinfl(m2)
      common/timinf/tau1,ampltp,timepe,dtau,time0
      common/cwalv/inbcvs,inbcvn,inbcvw,inbcve
      common/cwalp/inbcps,inbcpn,inbcpw,inbcpe
      common/cwals/inbcss,inbcsn,inbcsw,inbcse
      common/topoda/dim1i,dim1f,dim2i,dim2f,hmax,dhmax
      common/pscste/psstep
      if(inbcvn.eq.0.or.inbcvn.eq.2) then
          do j=1,n2m
           ym=yp2(n2)-(yp2(j+1)+yp2(j))*0.5
           uinfl(j)=u0*exp(-ym/blve)
          enddo
          do j=1,n2f
           ym=yp2(n2)-(yp2(j+1)+yp2(j))*0.5
           ssinf(1,j)=exp(-ym/blve)
          enddo
      endif
c 
c     PARABOLIC PROFILE: in the input file ccp is the
c     parabolic profile centre, u0 is the maximum of
c     velocity and blve is the width of the profile.
c 
      if(inbcvs.eq.0.or.inbcvs.eq.2.or.inbcvs.eq.1) then
      write(115,*)u0,blve,cpp
          ycpar=cpp 
          yipar=ycpar-blve*0.5
          yfpar=ycpar+blve*0.5
          do j=1,n2m
              uinfl(j)=0
              ym=(yp2(j+1)+yp2(j))*0.5
               if(ym.ge.yipar.and.ym.le.yfpar) then
                 uinfl(j)=-u0*((ym-yfpar)*
     1                    (ym-yipar))/((blve/2.)**2)
               endif
          enddo
          do j=1,n2f
             ssinf(1,j)=0
             if(yp2(j).ge.yipar.and.yp2(j).le.yfpar) then
               ssinf(1,j)=-u0*((yp2(j)-yfpar)*
     1                    (yp2(j)-yipar))/((blve/2.)**2)
             endif
             write(111,*)yp2(j),ssinf(1,j)
         enddo
       endif
c ..............................
      if(inbcvs.eq.-1) then
      do j=1,n2m
      ym=(yp2(j+1)+yp2(j))*0.5
      ydim=ym/blve
      if(ydim.ge.1.) ydim=1.
      uinfl(j)=u0*ydim**(1./7.)
      enddo
      do j=n2i,n2f
      ssinf(1,j)=uinfl(j)
      enddo
                     endif
      if(inbcvs.eq.0.or.inbcvs.eq.2.or.inbcvs.eq.1) then
      psinf(1)=0.
      do j=1,n2m
      psinf(j+1)=psinf(j)+uinfl(j)/dx2
      enddo
      ybou=yp2(1)
                     endif
      if(inbcvs.eq.-1.or.inbcvs.eq.-2) then
      psinf(1)=0.
      do j=1,n2m
      psinf(j+1)=psinf(j)+uinfl(j)/dx2
      enddo
      ybou=yp2(1)
                     endif
      do j=1,n2
      psinf(j)=psinf(j)-psinf(n2)*0.5
      enddo
      do j=2,n2m
      voinf(j)=-(psinf(j+1)-2.*psinf(j)+psinf(j-1))*dx2q
      enddo
      voinf(1)=voinf(2)
      voinf(n2)=voinf(n2m)
      ssima=0.
      psima=0.
      velma=0.
      vorma=0.
      open(47,file='infprof.out')
      write(47,*)'inbcvs,inbcvn  ',inbcvs,inbcvn
      do j=1,n2f
      write(47,133) yp2(j),uinfl(j),psinf(j),voinf(j),
     1             (ssinf(n,j),n=1,npscf)
      vorma=max(vorma,abs(voinf(j)))
      velma=max(velma,abs(uinfl(j)))
      psima=max(psima,abs(psinf(j)))
      ssima=max(ssima,abs(ssinf(1,j)))
      enddo
      close(47)
  133 format(3x,8e12.4)
      write(6,*)' maximum inlet qua   ',velma,psima,vorma,ssima,mpsc,
     1                  npscf
      do i=1,n1
      do j=1,n2f
      psi(i,j)=0.
      vor(i,j)=0.
      enddo
      enddo
      do i=1,n1
      do j=1,n2f
      psc(1,i,j)=0.
      enddo
      enddo
      psstep=(dim1i-yp1(1))/(yp1(n1)-yp1(1))
      do j=1,n2f
      do i=1,n1
      x1nd=(yp1(i)-yp1(1))/(yp1(n1)-yp1(1))
      psc(2,i,j)=x1nd
      enddo
      enddo
      do j=n2i,n2f
      ssinf(2,j)=psc(2,1,j)
      enddo
      write(6,*)' values at the step of second scalar=',psstep
      return
      end
c  
c     In this subroutine we use a jet with a
c     PARABOLIC PROFILE: in the input file ccp is the
c     parabolic profile centre, u0 is the maximum of
c     velocity and blve is the width of the profile.
c 
      subroutine inljetsud(u0,blve,cpp,vor,psi,psc)
      include 'param.f'
      dimension vor(m1,m2),psi(m1,m2),psc(mpsc,m1,m2)
      common/dim/n1,n1m,n2,n2m
      common/indbo/imv(m1),ipv(m1)
      common/mesh/dx1,dx1q,dx2,dx2q
      common/coor/yp1(m1),yp2(m2)
      common/metri/g2m(m1),g2c(m1)
      common/pscqu/sch(mpsc),scla,npscf
      common/infmxl/voins(m1),psins(m1),ssins(mpsc,m1)
c     common/infmxp/voinp(m1),psifr(m1),
c    %              psifi(m1),psiar(m1),psiai(m1),psinp(m1)
      dimension uins(m1),uinpe(m1),uinc(m1)
      common/timinf/tau1,ampltp,timepe,dtau,time0
      common/iperi/ib2per
      common/indpe2/n2i,n2f
      common/pscste/psstep
      dimension u(m1)
      write(6,*)' plane jet space developing     '
         xcpar=cpp
         xipar=xcpar-blve*0.5
         xfpar=xcpar+blve*0.5
         do i=1,n1m
            uins(i)=0
            xm=(yp1(i+1)+yp1(i))*0.5
            if(xm.ge.xipar.and.xm.le.xfpar) then
             uins(i)=-u0*((xm-xfpar)*
     1                    (xm-xipar))/((blve/2.)**2)
            endif
         enddo
          do i=1,n1
              ssins(1,i)=uins(i)
         enddo
      psins(1)=0.
      do i=1,n1m
      psins(i+1)=psins(i)-uins(i)/dx1*g2m(i)
      enddo
      do i=1,n1
      psins(i)=psins(i)-psins(n1)*0.5
      enddo
c ......................
      vomas=0.
      velma=0.
      do i=2,n1m
      ip=i+1
      im=i-1
      dpsx1p=(psins(ip)-psins(i))/g2m(i)*dx1q/g2c(i)
      dpsx1m=(psins(i)-psins(im))/g2m(im)*dx1q/g2c(i)
      dpsx1=dpsx1p-dpsx1m
      voins(i)=-( dpsx1)
      vomas=max(abs(voins(i)),vomas)
      velma=max(abs(u(i)),velma)
      uinc(i)=(uins(i)+uins(i-1))*0.5
        write(44,*) yp1(i),voins(i)
      enddo
      uinc(1)=uins(2)
      uinc(n1)=uins(n1m)
      voins(1)=voins(2)
      voins(n1)=voins(n1m)
      do npsc=1,npscf
      do i=1,n1
      ssins(npsc,i)=abs(uinc(i)-uinc(1))
      enddo
      enddo
      do i=1,n1
      do j=n2i,n2f
      vor(i,j)=0.
      psi(i,j)=0.
      enddo
      enddo
      do npsc=1,npscf
      do i=1,n1
      do j=n2i,n2f
      psc(npsc,i,j)=0.
      enddo
      enddo
      enddo
      ssima=0.
      psima=0.
      vorma=0.
      open(47,file='infprofr.out')
      do j=1,n1
      write(47,135) yp1(j),uins(j),psins(j),voins(j),
     1             (ssins(n,j),n=1,npscf)
      vorma=max(vorma,abs(voins(j)))
      psima=max(psima,abs(psins(j)))
      ssima=max(ssima,abs(ssins(1,j)))
      enddo
      close(47)
  135 format(3x,8e12.4)
      write(6,*)' maximum inlet mxl   ',velma,psima,vorma,ssima,mpsc,
     1                  npscf
c     write(6,*)' values at the step of second scalar=',psstep
      return
      end
c
c
c  ****************************** subrout quainl **********************
c
c  this subroutine calculate the integral  quantities
c  vorticity and enstrophy and energy
c
c
c     controllare il modo di calcolare l'energia
c
c
      subroutine quainl(time,vor,psi,psc)
      include 'param.f'
      common/d1/alx1i,alx1f,alx2i,alx2f
      dimension vor(m1,m2),psi(m1,m2),psc(mpsc,m1,m2)
      common/vely/v2(m1,m2)
      common/dim/n1,n1m,n2,n2m
      common/mesh/dx1,dx1q,dx2,dx2q
      common/metri/g2m(m1),g2c(m1)
      common/indbo/imv(m1),ipv(m1)
      common/coor/yp1(m1),yp2(m2)
      common/coorh/ym1(m1)
      common/indpe2/n2i,n2f
      common/parinf/u0,blve,cpp
      common/topoda/dim1i,dim1f,dim2i,dim2f,hmax,dhmax
      common/pscste/psstep
      n1m=n1-1
      n1mt=n1m
      vortp=0.
      vortm=0.
      pma=-.1e+04
      vma=-.1e+04
      vmi=+.1e+04
      ares1p=0.
      ares2p=0.
      udx=1./(dx2)
      uste=2.*blve
      do j=1,n2
      do i=1,n1
      undx=udx*(ym1(i)-ym1(i-1))
      if(vor(i,j).lt.vmi) then
      vmi=vor(i,j)
      yv1mi=yp1(i)
      yv2mi=yp2(j)
                          endif
      if(yp2(j).gt.uste.and.yp1(i).lt.dim1i) then
      if(psc(1,i,j).gt.vma) then
      pma=psc(1,i,j)
      yp1ma=yp1(i)
      yp2ma=yp2(j)
                          endif
      if(vor(i,j).gt.vma) then
      vma=vor(i,j)
      yv1ma=yp1(i)
      yv2ma=yp2(j)
                          endif
      if(vor(i,j).gt.0.) then
      vortp=vortp+vor(i,j)*undx
                         endif
      if(psc(1,i,j).gt.0.1) then
      ares1p=ares1p+undx
                           endif
                                              endif
      if(vor(i,j).le.0.) vortm=vortm+vor(i,j)*undx
      if(yp1(i).lt.dim1i) then
      if(psc(2,i,j).gt.psstep) then
      ares2p=ares2p+undx
                               endif
                           endif
      enddo
      enddo
      psin=psi(1,n2m)-psi(1,1)
      psou=psi(n1,n2m)-psi(n1,1)
  139 format(4e15.6)
      write(51,139)time,vortp
      write(52,139)time,vortm
      write(53,139)yv1ma,yv2ma 
      write(54,139)yv1mi,yv2mi 
      write(55,139)time,ares1p
      write(56,139)time,ares2p
      write(57,139)yp1ma,yp2ma 
      write(58,139)time,psin 
      write(59,139)time,psou
      return
      end
c
c
c  ****************************** subrout inimlt  **********************
c
c   initial vort  psi and psc for the time developing mixing layer
c
      subroutine inimlt(vor,psi,psc,a1mpl,a2mpl)
      include 'param.f'
      dimension vor(m1,m2),psi(m1,m2),psc(mpsc,m1,m2)
      common/dim/n1,n1m,n2,n2m
      common/mesh/dx1,dx1q,dx2,dx2q
      common/coor/yp1(m1),yp2(m2)
      common/pi/pi
      common/indpe2/n2i,n2f
      dimension u(m1)
      common/timinf/tau1,ampltp,timepe,dtau,time0
      common/cwalv/inbcvs,inbcvn,inbcvw,inbcve
      common/cwalp/inbcps,inbcpn,inbcpw,inbcpe
      common/cwals/inbcss,inbcsn,inbcsw,inbcse
      common/inmodo/inmod
      call uini(u)
      call vorini(vor,u,psi,a1mpl,a2mpl)
      avomax=0.
      vomax=-100.
      vomin=100.
      sfmax=-100.
      sfmin=100.
      scmax=-100.
      scmin=100.
      do i=1,n1
      do j=n2i,n2f
      vomax=max(vor(i,j),vomax)
      vomin=min(vor(i,j),vomin)
      sfmax=max(psi(i,j),sfmax)
      sfmin=min(psi(i,j),sfmin)
      avomax=max(abs(vor(i,j)),avomax)
      enddo
      enddo
      pi=2.*asin(1.)
      do i=1,n1
      do j=n2i,n2f
      arg=yp1(i)*sqrt(pi)
      psc(1,i,j)=0.5*(1+erf(arg))
      scmax=max(psc(1,i,j),scmax)
      scmin=min(psc(1,i,j),scmin)
      enddo
      enddo
      write(6,*)' inmod=',inmod,'  vomax,sfmax,scmax ',vomax,sfmax,scmax
      write(6,*)'  vomin,sfmin,scmin ',vomin,sfmin,scmin
      return
      end
c                                                                       
c
c  ****************************** subrout inlmls **********************
c
c   initial conditions  for the space developing mixing layer
c
      subroutine inlmls(vor,psi,psc,u0,blve)
      include 'param.f'
      dimension vor(m1,m2),psi(m1,m2),psc(mpsc,m1,m2)
      complex u1(m1),u2(m1)
      common/dim/n1,n1m,n2,n2m
      common/indbo/imv(m1),ipv(m1)
      common/mesh/dx1,dx1q,dx2,dx2q
      common/inmodo/inmod
      common/pi/pi
      common/coor/yp1(m1),yp2(m2)
      common/metri/g2m(m1),g2c(m1)
      common/pscqu/sch(mpsc),scla,npscf
      common/infmxl/voins(m1),psins(m1),ssins(mpsc,m1)
      common/infmxp/voinp(m1),psifr(m1),
     %              psifi(m1),psiar(m1),psiai(m1),psinp(m1)
      common/infvop/voinfr(m1),voinfi(m1),
     %              voinar(m1),voinai(m1)
      dimension uins(m1),uinpe(m1),uinc(m1)
      common/timinf/tau1,ampltp,timepe,dtau,time0
      common/iperi/ib2per
      common/indpe2/n2i,n2f
      dimension u(m1)
      write(6,*)' mixlay space developing     '
      u2i=blve/(1.-blve)
      u1i=u0+u2i
      write(6,*)' u1=',u1i,'   u2=',u2i,' u1/u2=',blve
      write(6,*)' time var tau1,ampltp,timepe,dtau,time0',
     1          tau1,ampltp,timepe,dtau,time0
      do i=1,n1m
      ym=(yp1(i)+yp1(i+1))*0.5
      uins(i)=((1+blve)/(1-blve)+tanh(2.*ym))*0.5
      ydb=ym/blve
c     uinpe(i)=exp(-(ydb**2))
      uinpe(i)=0.
      enddo
      call uinper(u1,u2)
      psins(1)=0.
      psifr(1)=0.
      psifi(1)=0.
      psiar(1)=0.
      psiai(1)=0.
      do i=2,n1
      psins(i)=psins(i-1)-uins(i-1)/dx1*g2m(i-1)
      psifr(i)=psifr(i-1)-real(u1(i-1))/dx1*g2m(i-1)
      psifi(i)=psifi(i-1)-imag(u1(i-1))/dx1*g2m(i-1)
      psiar(i)=psiar(i-1)-real(u2(i-1))/dx1*g2m(i-1)
      psiai(i)=psiai(i-1)-imag(u2(i-1))/dx1*g2m(i-1)
      enddo
      do i=1,n1
        write(42,*) yp1(i),psins(i),psifr(i),psifi(i),
     %              psiar(i),psiai(i)
      end do
      close(42)
      vomas=0.
      velma=0.
      do i=2,n1m                                                    
      ip=i+1                                                            
      im=i-1                                                            
      dpsx1p=(psins(ip)-psins(i))/g2m(i)*dx1q/g2c(i)              
      dpsx1m=(psins(i)-psins(im))/g2m(im)*dx1q/g2c(i)
      dpsx1=dpsx1p-dpsx1m
      voins(i)=-( dpsx1)                 
      dpsx1p=(psifr(ip)-psifr(i))/g2m(i)*dx1q/g2c(i)              
      dpsx1m=(psifr(i)-psifr(im))/g2m(im)*dx1q/g2c(i)
      dpsx1=dpsx1p-dpsx1m
      voinfr(i)=-( dpsx1)                 
      dpsx1p=(psifi(ip)-psifi(i))/g2m(i)*dx1q/g2c(i)              
      dpsx1m=(psifi(i)-psifi(im))/g2m(im)*dx1q/g2c(i)
      dpsx1=dpsx1p-dpsx1m
      voinfi(i)=-( dpsx1)                 
      dpsx1p=(psiar(ip)-psiar(i))/g2m(i)*dx1q/g2c(i)              
      dpsx1m=(psiar(i)-psiar(im))/g2m(im)*dx1q/g2c(i)
      dpsx1=dpsx1p-dpsx1m
      voinar(i)=-( dpsx1)                 
      dpsx1p=(psiai(ip)-psiai(i))/g2m(i)*dx1q/g2c(i)              
      dpsx1m=(psiai(i)-psiai(im))/g2m(im)*dx1q/g2c(i)
      dpsx1=dpsx1p-dpsx1m
      voinai(i)=-( dpsx1)                 
      vomas=max(abs(voins(i)),vomas)
      velma=max(abs(uins(i)),velma)
      uinc(i)=(uins(i)+uins(i-1))*0.5
        write(43,155) yp1(i),voins(i),voinfi(i),
     %              voinar(i),voinai(i),uinc(i)
  155 format(3x,7e12.5)
      enddo
      uinc(1)=uins(2)
      uinc(n1)=uins(n1m)
      voins(1)=voins(2)
      voins(n1)=voins(n1m)
      voinfr(1)=voinfr(2)
      voinfr(n1)=voinfr(n1m)
      voinfi(1)=voinfi(2)
      voinfi(n1)=voinfi(n1m)
      voinar(1)=voinar(2)
      voinar(n1)=voinar(n1m)
      voinai(1)=voinai(2)
      voinai(n1)=voinai(n1m)
      do i=1,n1,n1m
        write(43,155) yp1(i),voins(i),voinfr(i),voinfi(i),
     %              voinar(i),voinai(i)
      end do
      do npsc=1,npscf
      do i=1,n1
      ssins(npsc,i)=abs(uinc(i)-uinc(1))
      enddo
      enddo
      write(43,*)'  ssins done'
      do i=1,n1
      do j=n2i,n2f
c     vor(i,j)=0.
c     psi(i,j)=0.
      vor(i,j)=voins(i)
      psi(i,j)=psins(i)
      do npsc=1,npscf
c     psc(npsc,i,j)=0.
      psc(npsc,i,j)=ssins(npsc,i)
      enddo
      enddo
      enddo
      write(43,*)'  psi and vor =0 in the field '
      close(43)
      ssima=0.
      psima=0.
      vorma=0.
      open(47,file='in5prof.out')
      do i=1,n1
      write(47,135) yp1(i),uinc(i),psins(i),voins(i),
     1             (ssins(n,i),n=1,npscf),uinpe(i),psinp(i),voinp(i)
      vorma=max(vorma,abs(voins(i)))
      psima=max(psima,abs(psins(i)))
      ssima=max(ssima,abs(ssins(1,i)))
      enddo
      close(47)
  135 format(3x,8e12.4)
      write(6,*)' maximum inlet mxl   ',velma,psima,vorma,ssima,mpsc,
     1                  npscf
      return
      end
c                                                                       
c    ******************* uinper *******************************           
c   perturbation for the  space developing mixing layer
c                                                                       
      subroutine uinper(u1,u2)                                                
      include 'param.f'
      parameter (nmax=500)
      common/dim/n1,n1m,n2,n2m                                          
      common/mesh/dx1,dx1q,dx2,dx2q                                     
      common/d1/alx1i,alx1f,alx2i,alx2f
      common/coor/yp1(m1),yp2(m2)
      common/metri/g2m(m1),g2c(m1)
      dimension u(m1),v(m1)                                            
      complex u1(m1),u2(m1)
      dimension y(nmax),f2(nmax),f3(nmax),f4(nmax),f5(nmax)   
      dimension y2(nmax)                                                  
      common/per/fzc(m1),fzs(m1)                                        
c                                                                       
c                                                                       
      nn2=78                                                            
c                                                                       
      open(69,file='modfon.data')
      do 829 j=1,nn2                                                    
      read(69,*,end=1000) y(j),f2(j),f3(j),f4(j),f5(j)                         
  829 continue                                                          
 1000 close(69)
      npun = j-1
      yma = -10.
      ymi = +10.
      do i=1,npun
        yma = max(yma,y(i))
        ymi = min(ymi,y(i))
      end do
      do i=1,n1m
        if(yp1(i).lt.ymi) then
          iymi = i
        end if
      end do
      do i=1,n1m
        if(yp1(i).gt.yma) then
          iyma = i
          go to 2000
        end if
      end do
 2000 continue
c                                                                       
      call spline(y,f2,nn2,0.,0.,y2)                                      
c                                                                       
      do i=1,n1m                                                        
      xx=(yp1(i+1)+yp1(i))*.5
      call splint(y,f2,y2,nn2,xx,yy)                                    
      u(i)=yy                                                           
      end do                                                            
c                                                                       
c                                                                       
      call spline(y,f3,nn2,0.,0.,y2)                                      
c                                                                       
      do i=1,n1m                                                        
      xx=(yp1(i+1)+yp1(i))*.5
      call splint(y,f3,y2,nn2,xx,yy)                                    
      fzc(i)=yy                                                         
      end do                                                            
      do i=1,n1m                                                        
        are=u(i)*cos(fzc(i))
        aim=u(i)*sin(fzc(i))
        u1(i) = cmplx(are,aim)
      end do                                                            
c                                                                       
c                                                                       
      call spline(y,f4,nn2,0.,0.,y2)                                      
c                                                                       
      do i=1,n1m                                                        
      xx=(yp1(i+1)+yp1(i))*.5
      call splint(y,f4,y2,nn2,xx,yy)                                    
      v(i)=yy                                                         
      end do                                                            
      enfon = 0.
      do i=iymi,iyma                                                        
        enfon = enfon + 0.5*( u(i)**2 + v(i)**2 )/dx1*g2m(i)
      end do                                                            
      write(6,*) ' ENFON = ',enfon
c                                                                       
      open(69,file='modarm.data')
      do 830 j=1,nn2                                                    
      read(69,*) y(j),f2(j),f3(j),f4(j),f5(j)                         
  830 continue                                                          
      close(69)
c                                                                       
      call spline(y,f2,nn2,0.,0.,y2)                                      
c                                                                       
      do i=1,n1m                                                        
      xx=(yp1(i+1)+yp1(i))*.5
      call splint(y,f2,y2,nn2,xx,yy)                                    
      u(i)=yy                                                           
      end do                                                            
c                                                                       
c                                                                       
      call spline(y,f3,nn2,0.,0.,y2)                                      
c                                                                       
      do i=1,n1m                                                        
      xx=(yp1(i+1)+yp1(i))*.5
      call splint(y,f3,y2,nn2,xx,yy)                                    
      fzc(i)=yy                                                         
      end do                                                            
      do i=1,n1m                                                        
        are=u(i)*cos(fzc(i))
        aim=u(i)*sin(fzc(i))
        u2(i) = cmplx(are,aim)
      end do                                                            
c                                                                       
      call spline(y,f4,nn2,0.,0.,y2)                                      
c                                                                       
      do i=1,n1m                                                        
      xx=(yp1(i+1)+yp1(i))*.5
      call splint(y,f4,y2,nn2,xx,yy)                                    
      v(i)=yy                                                         
      end do                                                            
      enarm = 0.
      do i=iymi,iyma                                                        
        enarm = enarm + 0.5*( u(i)**2 + v(i)**2 )/dx1*g2m(i)
      end do                                                            
      write(6,*) ' ENARM = ',enarm
CC    pause
c                                                                       
      do i=iyma,n1m                                                        
         u(i) = u(iyma-1)
         u1(i) = u1(iyma-1)
         u2(i) = u2(iyma-1)
      end do
      do i=1,iymi
         u(i) = u(iymi+1)
         u1(i) = u1(iymi+1)
         u2(i) = u2(iymi+1)
      end do
c     do i=1,n1
c       write(433,*) yp1(i),u(i),real(u1(i)),imag(u1(i)),
c    %              real(u2(i)),imag(u2(i))
c     end do
c     close(433)

      return                                                            
      end                                                               
c                                                                       
c
c  ****************************** subrout injets **********************
c
c   initial conditions  planar jet space developing
c
      subroutine injets(vor,psi,psc,u0,blve)
      include 'param.f'
      dimension vor(m1,m2),psi(m1,m2),psc(mpsc,m1,m2)
      common/dim/n1,n1m,n2,n2m
      common/indbo/imv(m1),ipv(m1)
      common/mesh/dx1,dx1q,dx2,dx2q
      common/inmodo/inmod
      common/pi/pi
      common/coor/yp1(m1),yp2(m2)
      common/metri/g2m(m1),g2c(m1)
      common/pscqu/sch(mpsc),scla,npscf
      common/infmxl/voins(m1),psins(m1),ssins(mpsc,m1)
      common/infmxp/voinp(m1),psifr(m1),
     %              psifi(m1),psiar(m1),psiai(m1),psinp(m1)
      dimension uins(m1),uinpe(m1),uinc(m1)
      common/timinf/tau1,ampltp,timepe,dtau,time0
      common/iperi/ib2per
      common/indpe2/n2i,n2f
      dimension u(m1)
      write(6,*)' plane jet space developing     '
         ycpar=u0
         yipar=ycpar-blve*0.5
         yfpar=ycpar+blve*0.5
         do i=1,n1m
            uins(i)=0
            uinpe(i)=0.
           ym=(yp1(i+1)+yp1(i))*0.5
           if(ym.ge.yipar.and.ym.le.yfpar) then
           uins(i)=(ym-yipar)*(yfpar-ym)*4.
           endif
         enddo
      write(6,*)' ycpar=',ycpar,'   blve=',blve,'  icc=',icc
      n1mh=n1m/2+1
      psins(1)=0.
      do i=1,n1m
      psins(i+1)=psins(i)+uins(i)/dx1*g2m(i)
      enddo
c     do i=1,n1
c     psins(i)=psins(i)-psins(n1)*0.5
c     enddo
      write(6,*)' time var tau1,ampltp,timepe,dtau,time0',
     1          tau1,ampltp,timepe,dtau,time0
      do i=1,n1
c     psinp(i)=psins(i)
      psinp(i)= 0.
      enddo
      vomas=0.
      velma=0.
      do i=2,n1m                                                    
      ip=i+1                                                            
      im=i-1                                                            
      dpsx1p=(psins(ip)-psins(i))/g2m(i)*dx1q/g2c(i)              
      dpsx1m=(psins(i)-psins(im))/g2m(im)*dx1q/g2c(i)
      dpsx1=dpsx1p-dpsx1m
      voins(i)=-( dpsx1)                 
      dpsx1p=(psinp(ip)-psinp(i))/g2m(i)*dx1q/g2c(i)              
      dpsx1m=(psinp(i)-psinp(im))/g2m(im)*dx1q/g2c(i)
      dpsx1=dpsx1p-dpsx1m
      voinp(i)=-( dpsx1)                 
      vomas=max(abs(voins(i)),vomas)
      velma=max(abs(u(i)),velma)
      uinc(i)=(uins(i)+uins(i-1))*0.5
      enddo
      uinc(1)=uins(2)
      uinc(n1)=uins(n1m)
      voins(1)=voins(2)
      voins(n1)=voins(n1m)
      do npsc=1,npscf
      do i=1,n1
      ssins(npsc,i)=abs(uinc(i)-uinc(1))
      enddo
      enddo
      do i=1,n1
      do j=n2i,n2f
c     vor(i,j)=voins(i)
c     psi(i,j)=psins(i)
      vor(i,j)=0.
      psi(i,j)=0.
      enddo
      enddo
      do npsc=1,npscf
      do i=1,n1
      do j=n2i,n2f
      psc(npsc,i,j)=ssins(npsc,i)
c     psc(npsc,i,j)=0.
      enddo
      enddo
      enddo
      ssima=0.
      psima=0.
      vorma=0.
      open(47,file='in6prof.out')
      do j=1,n1
      write(47,135) yp1(j),uins(j),psins(j),voins(j),
     1             (ssins(n,j),n=1,npscf),uinpe(j),psinp(j),voinp(j)
      vorma=max(vorma,abs(voins(j)))
      psima=max(psima,abs(psins(j)))
      ssima=max(ssima,abs(ssins(1,j)))
      enddo
      close(47)
  135 format(3x,8e12.4)
      write(6,*)' maximum inlet mxl   ',velma,psima,vorma,ssima,mpsc,
     1                  npscf
      return
      end
c
c ***************************************************************+
c    evaluates the topography for the flow in Sect.4.12.2e
c
      subroutine maketop
      include 'param.f'
      common/dim/n1,n1m,n2,n2m
      common/topog/htop(m1,m2)
      dimension hbou(m1),jth(m1)
      common/topoda/dim1i,dim1f,dim2i,dim2f,hmax,dhmax
      common/itop/itopog
      common/coor/yp1(m1),yp2(m2)
      if(itopog.eq.0) then
      do j=1,n2
      do i=1,n1
      htop(i,j)=0.
      topoma=0.
      enddo
      enddo
                      else
      write(6,*)'in maketo  dim1i,dim1f,hmax,dhmax',
     1                      dim1i,dim1f,hmax,dhmax
      topoma=-100.
      topomi=100.
      do j=1,n2
      do i=1,n1
      yd1=yp1(i)-yp1(1)
      if(yd1.ge.dim1i)  then
      dd1=(yp1(i)-dim1i)/dim1f
      if(j.eq.1) then
      write(67,*)'**************',yp1(i),dim1i,dim1f,dd1
                 endif
      if(dd1.ge.1.) then
      fdd1=1.
                    else
      fdd1=3.*dd1**2-2.*dd1**3
                    endif
      htop(i,j)=hmax+dhmax*fdd1
                       endif
      if(yd1.le.dim1i)  then
      dd1=yd1
      htop(i,j)=hmax
                        endif
      if(j.eq.1) then
      write(67,*)yd1,dd1,fdd1,htop(i,j)
                 endif
      topoma=max(htop(i,j),topoma)
      topomi=min(htop(i,j),topomi)
      enddo
      enddo
                      endif
      write(20,*)'  max htop  ',topoma,'   min  htop',topomi
      write(6,*)'  max htop  ',topoma,'   min  htop',topomi
      return
      end
c
c
c  ****************************** subrout dipol  **********************
c
c   initial vort. and psi for Lamb dipole of Batchlor book pg 535
c
      subroutine dipol(vor,ramo,vsi,yc1mo,yc2mo,akmo,velmo,psc)
      include 'param.f'
      dimension vor(m1,m2),psc(mpsc,m1,m2)
      common/dim/n1,n1m,n2,n2m
      common/coor/yp1(m1),yp2(m2)
      common/angmd/thet0
      common/anmod/ar1,ar2
      common/pi/pi
      common/indpe2/n2i,n2f
      common/pscqu/sch(mpsc),scla,npscf
      common/inmodo/inmod
c
      ji=1
      jf=n2
      ii=1
      if=n1
      vorma=-100.
      vormi=+100.
      do 917 i=ii,if
      do 917 j=n2i,n2f
         y1d=yp1(i)-yc1mo
      y2d=yp2(j)-yc2mo
      ramod=sqrt(y1d**2+y2d**2)
      if(ramod.lt..1e-06)then
      vor(i,j)=0.
      else
      stho=y2d/ramod
      if(abs(stho).gt.1.) write(6,719) i,j,y2d,ramod
  719 format(3x,2i4,2x,'err in asin',2e12.4)
      astho=(abs(stho)-1.)
c
c     definizione theta
c
      if(astho.gt..0) stho=stho-.1e-06*stho/abs(stho)
      thet=asin(stho)
      if(y2d.lt.0..and.y1d.gt.0.) thet=pi*2.-abs(thet)
      if(y2d.lt.0..and.y1d.lt.0.) thet=pi+abs(thet)
      if(y2d.ge.0..and.y1d.lt.0.) thet=-abs(thet)+pi
c
      sth=sin(thet)*cos(thet0)-sin(thet0)*cos(thet)
      ramo=sqrt(ar1**2*ar2**2/(ar2**2+(ar1**2-ar2**2)*sth**2))
      radmod=ramod/ramo
      if(radmod.ge.1.) then
      vor(i,j)=0.
      else
c
c   
c   Bessel functions **********************
c   
      arbesj=radmod*akmo
      vor(i,j)=2.*vsi*velmo*sth*bessj1(arbesj)/bessj0(akmo)*akmo
      vorma=max(vor(i,j),vorma)
      vormi=min(vor(i,j),vormi)
      endif
      endif
  917 continue
      write(6,791) vorma,vormi
  791 format(3x,'***  in dipol vormax =',e12.4,'  vormi = ',e12.4)
      vomax=max(abs(vorma),abs(vormi))
      scmax=0.
      do i=1,n1
      do j=n2i,n2f
      psc(1,i,j)=abs(vor(i,j))/vomax
      scmax=max(abs(psc(1,i,j)),scmax)
      enddo
      enddo
      if(npscf.eq.2) then
      do j=n2i,n2f
      do i=1,n1
      psc(2,i,j)=0.
      yd=(yp1(n1)-yp1(i))/scla
      if(yd.le.1.) then
      fy=3.*yd**2-2.*yd**3
      psc(2,i,j)=1.-fy
                     endif
      enddo
      enddo
                     endif
      write(6,*)'inmod=',inmod,'npscf=',npscf,
     1          ' voma,scma',vomax,scmax
      return
      end
c
c  ****************************** subrout modon  **********************
c
c   initial vort. and psi for Stern modon 
c
      subroutine modon(vor,ramo,vsi,yc1mo,yc2mo,velmo,psc)
      include 'param.f'
      dimension vor(m1,m2),psc(mpsc,m1,m2)
      common/pscqu/sch(mpsc),scla,npscf
      common/dim/n1,n1m,n2,n2m
      common/coor/yp1(m1),yp2(m2)
      common/angmd/thet0
      common/anmod/ar1,ar2
      common/betef/beta
      common/inmodo/inmod
      common/indpe2/n2i,n2f
      ramo=sqrt(ar1**2*ar2**2)
      call disp(velmo,ramo,akmo)
      pi=2.*asin(1.)
      besj0=bessj1(akmo)
      uokc=beta/sqrt(velmo)
      besk0=bessk1(uokc)
      write(6,*)' for modon    ',uokc,velmo,akmo
      vorma=-100.
      vormi=+100.
      do 917 i=1,n1m   
      do 917 j=2,n2m
      y1d=yp1(i)-yc1mo
      y2d=yp2(j)-yc2mo
      ramod=sqrt(y1d**2+y2d**2)
      if(ramod.lt..1e-06)then
      vor(i,j)=0.
                         else
      stho=y2d/ramod
      astho=(abs(stho)-1.)
      if(astho.gt..0) then
      stho=stho-.1e-06*stho/abs(stho)
                      endif
      thet=asin(stho)
      if(y2d.lt.0..and.y1d.gt.0.) then
      thet=pi*2.-abs(thet)
      endif
      if(y2d.lt.0..and.y1d.lt.0.) then
      thet=pi+abs(thet)
      endif
      if(y2d.ge.0..and.y1d.lt.0.) then
      thet=-abs(thet)+pi
      endif
      sth=sin(thet)*cos(thet0)-sin(thet0)*cos(thet)
      radmod=ramod/ramo
      if(radmod.le.1.) then
      arbesj=radmod*akmo
      besj1=bessj1(arbesj)
      vor(i,j)=vsi*beta*ramo*sth*besj1/besj0
                       else
      ark1=uokc*radmod
      besk1=bessk1(ark1)
      vor(i,j)=vsi*beta*ramo*sth*besk1/besk0
                       endif
      vorma=max(vor(i,j),vorma)
      vormi=min(vor(i,j),vormi)
                       endif
  917 continue
      write(6,791) vorma,vormi
  791 format(3x,' modon vorma=',e12.4,' vormi=',e12.4)
      vomax=max(abs(vorma),abs(vormi))
      scmax=0.
      do i=1,n1
      do j=n2i,n2f
      psc(1,i,j)=abs(vor(i,j))/vomax
      scmax=max(abs(psc(1,i,j)),scmax)
      enddo
      enddo
      if(npscf.eq.2) then
      do j=n2i,n2f
      do i=1,n1
      psc(2,i,j)=0.
      yd=(yp1(n1)-yp1(i))/scla
      if(yd.le.1.) then
      fy=3.*yd**2-2.*yd**3
      psc(2,i,j)=1.-fy
                     endif
      enddo
      enddo
                     endif
      write(6,*)'inmod=',inmod,'npscf=',npscf,
     1          ' voma,scma',vomax,scmax
      return
      end
      subroutine disp(velmo,ramo,akmo)
c     find the value of flambda for modon given beta,
c     a (radius) and c (speed)
      common/c1/c,rho,a
      common/betef/beta
      external func
      a=ramo
      c=velmo 
      tol=.000001
      rho=sqrt(beta/c)
      fj11=3.83171
      fj21=5.13   
      fj12=7.01559
      x1=fj11/a
      x2=fj21/a
100   continue
      y1=func(x1)
      y2=func(x2)
      akmo=ZBRENT(FUNC,X1,X2,TOL) 
      return
      end
c
c c     ********************* monpol *******************************
c    vorticity distribution for  a monopol as in Sect.4.11b.d
c
      subroutine monpol(vor,y1c,y2c,crad,amp,psc,vsi)
      include 'param.f'
      dimension vor(m1,m2),psc(mpsc,m1,m2)
      common/dim/n1,n1m,n2,n2m
      common/coor/yp1(m1),yp2(m2)
      common/angmd/thet0
      common/anmod/ar1,ar2
      common/mesh/dx1,dx1q,dx2,dx2q
      common/indpe2/n2i,n2f
      common/inmodo/inmod
      vorpm=0.
      vorma=-1.e04
      vormi=+1.e04
      m=n2
      ampl=amp
      do 594 i=1,n1
      vor(i,1)=0.0
      vor(i,n2)=0.0
      do 596 j=2,n2m
      y1r=yp1(i)-y1c
      y2r=yp2(j)-y2c
      ra2=(y1r**2+y2r**2)
      rad=sqrt(ra2/(ar1**2*ar2**2)) 
      rar=rad/(sqrt(2.)*crad)
      if(rad.gt.20.) then
      aexp=0.
      ampl=0.
      else
      aexp=exp(-rar**2)
      ampl=amp
      endif
      if(vsi.gt.0.) then
      vor(i,j)=ampl*aexp
                    else
      vor(i,j)=(1.-rar**2)*ampl*aexp
                    endif
      vormi=min(vormi,vor(i,j))
      vorma=max(vorma,vor(i,j))
 596  continue
 594  continue
      vomax=max(abs(vorma),abs(vormi))
      write(6,788)vorma,vormi
 788  format(
     1     2x,'vma=',e11.4,2x,'vmi=',e11.4)
      svor=0.
      do 196 j=1,n2
      do 196 i=1,n1m
      svor=svor+vor(i,j)
 196  continue
      svor=svor/float(n2*n1m)
      write(6,798)svor
      write(16,798)svor
 798  format(1x,'svor=',e11.4)
      scmax=0.
      do i=1,n1
      do j=n2i,n2f
      psc(1,i,j)=abs(vor(i,j))/vomax
      scmax=max(abs(psc(1,i,j)),scmax)
      enddo
      enddo
      write(6,*)'inmod=',inmod,' voma,scma',vomax,scmax
      return
      end
c
c c     ********************* tripol *******************************
c   vorticity distribution for an isolated vortex as in Sect.4.11b.c
c
      subroutine tripol(vor,y1c,y2c,crad,alpha,amp,vorp,psc)
      include 'param.f'
      dimension vor(m1,m2),vorp(m1,m2),psc(mpsc,m1,m2)
      common/dim/n1,n1m,n2,n2m
      common/coor/yp1(m1),yp2(m2)
      common/angmd/thet0
      common/anmod/ar1,ar2
      common/mesh/dx1,dx1q,dx2,dx2q
      common/indpe2/n2i,n2f
      common/inmodo/inmod
      dimension x(m2)
      vorpm=0.
      vorma=-1.e04
      vormi=+1.e04
      write(6,106)y1c,y2c,crad,alpha,amp
  106 format(3x,6e12.4)
      idum=-9 
      npa=m2
      do i=1,n1
      idum=idum+2*i
      call gerand(idum,npa,x)
      do j=1,n2
      y1r=yp1(i)-y1c
      y2r=yp2(j)-y2c
      ra2=(y1r**2+y2r**2)
      rad=sqrt(ra2)/sqrt(ar1**2*ar2**2)
      rar=rad/sqrt(crad)
      if(rad.gt.4.) then
      aexp=0.
      ampl=0.
      else
      aexp=exp(-rar**alpha)
      ampl=amp*exp(-(rar**alpha*alpha-2.)**2)
      endif
      vor(i,j)=(1.-rar**alpha*0.5*alpha)*aexp
      vormi=min(vormi,vor(i,j))
      vorma=max(vorma,vor(i,j))
c     xr=rand()
      xr=x(j)
c     write(98,198)i,j,rad,rar,aexp,ampl,xr
  198 format(2i5,3x,8e12.4)
      xra=xr*ampl
      vorp(i,j)=xra
      vorpm=vorpm+vorp(i,j)
      enddo
      enddo
      vorpm=vorpm/float(n1*n2)
      write(6,788)vorma,vormi,vorpm
 788  format(
     1       2x,'vma=',e11.4,2x,'vmi=',e11.4,2x,'vorpm=',e11.4)
      vomax=max(abs(vorma),abs(vormi))
      scmax=0.
      do i=1,n1
      do j=n2i,n2f
      psc(1,i,j)=abs(vor(i,j))/vomax
      scmax=max(abs(psc(1,i,j)),scmax)
      enddo
      enddo
      svor=0.
      vorpma=-1.e04
      vorpmi=+1.e04
      do j=1,n2
      do i=1,n1
      vorpe=vorp(i,j)-vorpm
      vorpmi=min(vorpmi,vorpe)
      vorpma=max(vorpma,vorpe)
      vor(i,j)=vor(i,j)+vorp(i,j)-vorpm
      svor=svor+vor(i,j)
      vorp(i,j)=vor(i,j)
      enddo
      enddo
      svor=svor/float(n2*n1)
      write(6,798)svor,vorpmi,vorpma
      write(16,798)svor,vorpmi,vorpma
 798  format(1x,'svor=',e11.4,
     1 2x,'vorpma=',e11.4,2x,'vorpmi=',e11.4)
      write(6,*)'inmod=',inmod,' voma,scma',vomax,scmax
      return
      end
c                                                                       
c                                                                       
c    ******************* uini *******************************           
c   velocity distribution for time developing mixing layers
c   as given by Rogers and Moser
c                                                                       
      subroutine uini(u)                                                
      include 'param.f'
      parameter (nmax=500)
      common/dim/n1,n1m,n2,n2m                                          
      common/mesh/dx1,dx1q,dx2,dx2q                                     
      common/d1/alx1i,alx1f,alx2i,alx2f
      common/coor/yp1(m1),yp2(m2)
      dimension u(m1)                                                  
      dimension y(nmax),f2(nmax),f3(nmax),f4(nmax),f5(nmax)                       
      dimension y2(nmax)                                                  
      common/per/fzc(m1),fzs(m1)                                        
c                                                                       
c                                                                       
      nn2=97                                                            
c                                                                       
      open(69,file='inmosit.d')
      do 829 j=1,nn2                                                    
      read(69,100) y(j),f2(j),f3(j),f4(j),f5(j)                         
  829 continue                                                          
  100 format(2x,5e14.7)                                                 
c                                                                       
      call spline(y,f2,nn2,0.,0.,y2)                                      
c                                                                       
      do i=1,n1m                                                        
      xx=(yp1(i+1)+yp1(i))*.5
      call splint(y,f2,y2,nn2,xx,yy)                                    
      u(i)=yy                                                           
      end do                                                            
c                                                                       
c                                                                       
      call spline(y,f3,nn2,0.,0.,y2)                                      
c                                                                       
      do j=1,n1                                                         
      xx=yp1(j)                                                          
      call splint(y,f3,y2,nn2,xx,yy)                                    
      fzc(j)=yy                                                         
      end do                                                            
c                                                                       
      call spline(y,f4,nn2,0.,0.,y2)                                      
c                                                                       
      do j=1,n1                                                         
      xx=yp1(j)                                                          
      call splint(y,f4,y2,nn2,xx,yy)                                    
      fzs(j)=yy                                                         
      end do                                                            
c                                                                       
c                                                                       
      return                                                            
      end                                                               
c                                                                       
c     ********************* vorini *******************************      
c    vorticity calculation from velocity profiles given in the uini routine
c                                                                       
      subroutine vorini(vor,u,psi,a1mpl,a2mpl)
      include 'param.f'
      common/dim/n1,n1m,n2,n2m                                          
      common/mesh/dx1,dx1q,dx2,dx2q                                     
      dimension vor(m1,m2),v2(m1),u(m1),psi(m1,m2),pso(m1),voo(m1)
      common/per/fzc(m1),fzs(m1)                                        
      common/wavnu/akx,nvor
      common/coor/yp1(m1),yp2(m2)
      common/coorh/ym1(m1)
      common/metri/g2m(m1),g2c(m1)
      common/indwa/jmv(m2),jpv(m2)
      common/bvedph/dpest(m2),dpwes(m2)
      pi=2.*asin(1.)
      dudxm=0.                                                          
c
c  attention in our convection the streamwise velocity is the 
c                 v2   component
c  to have the same sign of omega z v2 changes sign
c  as well as fzc and fzs
c
      do i=1,n1                                                  
      v2(i)=u(i)
      enddo
      do i=1,n1m-1                                                  
      dudx=  ( v2(i+1)-v2(i) )*dx1/g2m(i)
      dudxm=amax1(dudxm,abs(dudx))                                      
      enddo
      delom=2./dudxm                                                    
      ipr=n1m/2+1
      dxin=g2m(ipr)/dx1
      write(6,*) 'dudxm, delom  dxin ',dudxm,delom,dxin
      do j=1,n2m                                                    
      psi(1,j)=0.0                                                      
      do i=2,n1m/2+1
      psi(i,j)=psi(i-1,j)-v2(i-1)*(yp1(i)-yp1(i-1))
      enddo
      psi(n1,j)=0.0                                                     
      do i=n1m,n1m/2+1,-1
      psi(i,j)=psi(i+1,j)+v2(i)*(yp1(i+1)-yp1(i))
      enddo
      enddo
      pig=2.*asin(1.)                                                   
      psfma=0.                                                          
      do 196 j=1,n2                                                     
      do 196 i=1,n1m                                                    
      psfma=amax1(abs(psi(i,j)),psfma)                                  
 196  continue                                                          
      vomax=0.                                                          
      circ=0.                                                           
      idum=-9 
      npa=m1
      jpr=n2m/2+1
      do 597 j=1,n2m                                                    
      jp=jpv(j)                                                         
      jm=jmv(j)                                                         
      vor(1,j)=0.0                                                      
      vor(n1,j)=0.0                                                     
      do 597 i=2,n1m                                                    
      ip=i+1                                                            
      im=i-1                                                            
      dpsx1p=(psi(ip,j)-psi(i,j))/(yp1(ip)-yp1(i))                
      dpsx1m=(psi(i,j)-psi(im,j))/(yp1(i)-yp1(im))  
      dpsx1=(dpsx1p-dpsx1m)/(ym1(i)-ym1(i-1))
      dpsx2=(psi(i,jm)-2.*psi(i,j)+psi(i,jp))*dx2q                
      if(nvor.eq.1) then
      xx=(yp1(i)-yp1(n1m/2+1))
      vor(i,j)=-( dpsx2+dpsx1)                 
     1           -fzc(i)*cos(akx*(j-1)/dx2)                             
     1           -fzs(i)*sin(akx*(j-1)/dx2)                             
     1           -rand( )*.10*exp(-2.*xx*xx)
                     else
      vor(i,j)=-( dpsx2+dpsx1)                 
     1       -fzc(i)*cos(2.*pi*nvor*(j-1)/float(n2m))*a2mpl
     1       -fzs(i)*sin(2.*pi*nvor*(j-1)/float(n2m))*a2mpl
     1       -fzc(i)*cos(2.*pi*(j-1)/float(n2m))*a1mpl                        
     1       -fzs(i)*sin(2.*pi*(j-1)/float(n2m))*a1mpl                       
                     endif
      circ=circ+vor(i,j)                                                
      vomax=amax1(abs(vor(i,j)),vomax)                                  
 597  continue                                                          
      do j=1,n2m
      dpest(j)=psi(n1,j)
      dpwes(j)=psi(1,j)
      enddo
      call phcalp(vor,psi)                                               
      circn=0.                                                          
      do 599 j=1,n2m                                                    
      jp=jpv(j)                                                         
      jm=jmv(j)                                                         
      vor(1,j)=0.0                                                      
      vor(n1,j)=0.0                                                     
      do 599 i=2,n1m                                                    
      ip=i+1                                                            
      im=i-1                                                            
      dpsx2=(psi(i,jm)-2.*psi(i,j)+psi(i,jp))*dx2q                
      dpsx1p=(psi(ip,j)-psi(i,j))/(yp1(ip)-yp1(i))                
      dpsx1m=(psi(i,j)-psi(im,j))/(yp1(i)-yp1(im))  
      dpsx1=(dpsx1p-dpsx1m)/(ym1(i)-ym1(i-1))
      vor(i,j)=-( dpsx2+dpsx1)                 
      circn=circn+vor(i,j)                                              
 599  continue                                                          
      write(6,722)circ,circn                                            
 722  format(3x,'  circ  circn   ',2x,2e12.4)                           
      write(6,711)psfma,vomax                                           
 711  format(3x,'  psimax  vomax   ',2x,2e12.4)                         
      return                                                            
      end                                                               
c
c  ****************************** subrout quamlt **********************
c
c  this subroutine calculate the output quantities
c  for the time developing mixing layer
c
c
c
      subroutine quamlt(time,vor,psi,psc)
      include 'param.f'
      parameter(m1m=m1-1,m2m=m2-1)
      common/d1/alx1i,alx1f,alx2i,alx2f
      dimension vor(m1,m2),psi(m1,m2),psc(mpsc,m1,m2),v2me(m1)
      dimension en1(m2m/2+1),en2(m2m/2+1)
      common/indbo/imv(m1),ipv(m1)
      common/vely/v2(m1,m2)
      common/dim/n1,n1m,n2,n2m
      common/mesh/dx1,dx1q,dx2,dx2q
      common/metri/g2m(m1),g2c(m1)
      common/indwa/jmv(m2),jpv(m2)
      common/coor/yp1(m1),yp2(m2)
      common/coorh/ym1(m1)
      common/indpe2/n2i,n2f
      do i=1,n1m
      v2me(i)=0.
      do j=1,n2m
      v2l=-(psi(i+1,j)-psi(i,j))/(yp1(i+1)-yp1(i))
      v2me(i)=v2me(i)+v2l
      enddo
      v2me(i)=v2me(i)/n2m
      enddo
      deltm=0.
      do i=1,n1m
      deltm=deltm+(1.-v2me(i)**2)*(yp1(i+1)-yp1(i))/4.
      enddo
      ommdbr=vor(n1m/2+1,1)
      ommdce=vor(n1m/2+1,n2m/2+1)
      i=n1m/2+1
      ip=i+1
      im=i-1
      j=1
      jp=jpv(j)
      jm=jmv(j)
      dpsx2=(psi(i,jm)-2.*psi(i,j)+psi(i,jp))*dx2q                
      dpsx1p=(psi(ip,j)-psi(i,j))/g2m(i)*dx1q/g2c(i)              
      dpsx1m=(psi(i,j)-psi(im,j))/g2m(im)*dx1q/g2c(i)
      dpsx1=dpsx1p-dpsx1m
      stmdbr=dpsx2-dpsx1
      j=n2m/2+1
      jp=jpv(j)
      jm=jmv(j)
      dpsx2=(psi(i,jm)-2.*psi(i,j)+psi(i,jp))*dx2q                
      dpsx1p=(psi(ip,j)-psi(i,j))/g2m(i)*dx1q/g2c(i)              
      dpsx1m=(psi(i,j)-psi(im,j))/g2m(im)*dx1q/g2c(i)
      dpsx1=dpsx1p-dpsx1m
      stmdce=dpsx2-dpsx1
      do j=1,n2m                                                    
      jp = jpv(j)
      do i=1,n1m                                                    
      v2(i,j)=(psi(i,jp)-psi(i,j))*dx2                                  
      end do
      end do
      call spect(v2,en1)                                                 
      a10u = 2.* sqrt( en1(2) )
      a20u = 2.* sqrt( en1(3) )
      aaa = en1(2)
      bbb = en1(3)
      do j=1,n2m                                                    
      do i=1,n1m                                                    
      ip=ipv(i)
      v2(i,j)=-(psi(ip,j)-psi(i,j))*dx1/g2m(i)
      end do
      end do
      call spect(v2,en2)                                                 
      a10v = 2.* sqrt( en2(2) )
      a20v = 2.* sqrt( en2(3) )
      a10 = 2.* sqrt( (aaa + en2(2))/2. )
      a20 = 2.* sqrt( (bbb + en2(3))/2. )
  139 format(8e15.6)
      write(51,139)time,deltm,a10, a20,a10u,a20u,a10v,a20v
      write(52,139)time,ommdbr
      write(53,139)time,ommdce 
      write(54,139)time,stmdbr 
      write(55,139)time,stmdce
      return
      end

c                                                                       
c                                                                       
c    **************** spect*********************************            
c    necessary for the outputs of the time developing mixing layers     
c                                                                       
      subroutine spect(uu,en)                                           
      include'param.f'
      parameter(m1m=m1-1,m2m=m2-1)
      parameter (m2mh=m2m/2+1,m2md=m2m+2)  
      dimension uu(m1,m2),en(m2m/2+1)
      complex uucom(m2m/2+1,m1)
      common/mesh/dx1,dx1q,dx2,dx2q                                     
      common/metri/g2m(m1),g2c(m1)
      common/fftcom/ifxz(13),trigz(3*m2m/2+1)                             
      common/dim/n1,n1m,n2,n2m                                          
      real xr(m2md,m1),work(m2,m1)
      common/d1/alx1i,alx1f,alx2i,alx2f
      call fftfax(n2m,ifxz,trigz)
      alx2=alx2f-alx2i                                                  
      n2mh=n2m/2+1
      n2md=n2m+2
      n2mdu=n2md+2
      n2mu=n2m-1
      do 2 i=1,n1
      xr(1,i)=uu(i,n2m)
      do 1 j=1,n2m
      js=j+1
      xr(js,i)=uu(i,j)
   1  continue
      xr(n2md,i)=uu(i,1)
      do 3 j=n2mdu,m2md
      xr(j,i)=0.      
   3  continue
   2  continue
      call fft99(xr,work,trigz,ifxz,1,m2md,n2m,n1,-1)
      do 4 kk=1,n2mh
      kd=2*kk-1
      kp=2*kk
      do 4 i=1,n1
      uucom(kk,i)=cmplx(xr(kd,i),xr(kp,i))
   4  continue
      do 240 j=1,n2m/2+1                                                
      en(j)=0.                                                          
      do 241 i=1,n1m                                                    
      en(j)=en(j)+uucom(j,i)*conjg(uucom(j,i))/dx1*g2m(i)               
 241  continue                                                          
 240  continue                                                          
      return                                                            
      end                                                               
