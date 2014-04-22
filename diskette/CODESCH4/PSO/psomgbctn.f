c
c     *********************** subrout invtps **********************
c
c      this subroutine performs the computation of the
c      passive scalar it is similar to invtv which solves
c      dvor=vor(n+1-vor(n
c
c     *************************************************************
      subroutine invtps(dvor,al,ga,ro,vor,ru,visc,npsc)
c    
      include 'param.f'
      dimension dvor(m1,m2),vor(m1,m2)
      dimension ru(mpsc,m1,m2)
      dimension ami(md2),aci(md2),api(md2)
      dimension amj(md2),acj(md2),apj(md2),fj(md2,md2)
      common/dim/n1,n1m,n2,n2m
      common/tstep/dt
      common/pscqu/sch(mpsc),scla,npscf
      common/mesh/dx1,dx1q,dx2,dx2q
      common/metri/g2m(m1),g2c(m1)
      common/bchods/dssud(mpsc,m1),dsnor(mpsc,m1)
      common/bcveds/dsest(mpsc,m2),dswes(mpsc,m2)
      common/cwals/inbcss,inbcsn,inbcsw,inbcse
      m2v=m2
      m1v=m1
c
c     *********************************************
c       add to the non linear term the diffusive
c       derivatives 11 and 22
c     *********************************************
      do 1 jc=2,n2m
      jm=jc-1
      jp=jc+1
      do 1 ic=2,n1m
      ip=ic+1
      im=ic-1
      visder=dx1q*((vor(ip,jc)-vor(ic,jc))/g2m(ic)
     1            -(vor(ic,jc)-vor(im,jc))/g2m(ic-1))/g2c(ic)
     1      +dx2q*(vor(ic,jm)-2.*vor(ic,jc)+vor(ic,jp))
      rhs=(ga*dvor(ic,jc)+ro*ru(npsc,ic,jc))*dt +
     1     (visder*visc)*dt*al 
      ru(npsc,ic,jc)=dvor(ic,jc)
      dvor(ic,jc)=rhs
    1 continue
c
c     *****************************************************
c        compute dpsc sweeping along the x2 direction
c               wall boundaries
c     *****************************************************
      betady=dx2q*0.5*visc*dt*al
      accy=1.+betady*2.
      betadx=dx1q*0.5*visc*al*dt
      do 20 jc=2,n2m
      jm=jc-1
      acj(jc)=accy
      amj(jc)=-betady
      apj(jc)=-betady
   20 continue
c
c   here the different possibilities for the north south
c   boundary conditions are available
c
      if(inbcsn.eq.1) then
      acj(n2)=1.
      apj(n2)=0.
      amj(n2)=0.
                      endif
      if(inbcsn.eq.0) then
      acj(n2)=1.
      apj(n2)=0.
      amj(n2)=-1.
                      endif
      if(inbcss.eq.1.or.inbcss.eq.-2) then
      acj(1)=1.
      apj(1)=0.
      amj(1)=0.
                      endif
      if(inbcss.eq.0) then
      acj(1)=1.
      apj(1)=-1.
      amj(1)=0.
                      endif
c
c    here the second derivatives are computed for the 
c    changes in b.c. due to the factorization
c
      do ic=2,n1m
      dvor(ic,n2)=dsnor(npsc,ic)-betadx*(
     1           (dsnor(npsc,ic+1)-dsnor(npsc,ic))/g2m(ic-1)
     1          -(dsnor(npsc,ic)-dsnor(npsc,ic-1))/g2m(ic))
      dvor(ic,1)=dssud(npsc,ic)-betadx*(
     1           (dssud(npsc,ic+1)-dssud(npsc,ic))/g2m(ic-1)
     1          -(dssud(npsc,ic)-dssud(npsc,ic-1))/g2m(ic))
      enddo
      do ic=2,n1m
      do jc=1,n2
      fj(ic,jc)=dvor(ic,jc)
      enddo
      enddo
c!
      call tribk(amj,acj,apj,fj,n2,2,n1m)
c!
      do ic=2,n1m
      do jc=1,n2
      dvor(ic,jc)=fj(ic,jc)
      enddo
      enddo
c
c     **************************************************
c        compute dpsc  sweeping in the x1 direction
c     **************************************************      
      do 3 ic=2,n1m
      aci(ic)=1.+betadx*(1./g2m(ic-1)+1./g2m(ic))/g2c(ic)
      ami(ic)=-betadx/g2m(ic-1)/g2c(ic)
      api(ic)=-betadx/g2m(ic)/g2c(ic)
    3 continue
c
c   here the different possibilities for the east west 
c   boundary conditions are available
c
      ic=1
      if(inbcsw.eq.0) then
      aci(ic)=1.
      ami(ic)=0.
      api(ic)=-1.
                      endif
      if(inbcsw.eq.1.or.inbcsw.eq.-2) then
      aci(ic)=1.
      ami(ic)=0.
      api(ic)=0.
                      endif
      ic=n1
      if(inbcse.eq.1) then
      aci(ic)=1.
      ami(ic)=0.
      api(ic)=0.
                      endif
      if(inbcse.eq.0) then
      aci(ic)=1.
      ami(ic)=-1.
      api(ic)=0.
                      endif
      do 21 jc=2,n2m
        fj(jc,1)=dswes(npsc,jc)
        fj(jc,n1)=dsest(npsc,jc)
      do 21 ic=2,n1m
        fj(jc,ic)=dvor(ic,jc)
   21 continue     
c!
      call tribk(ami,aci,api,fj,n1,2,n2m)
c!
      do ic=2,n1m
      do jc=2,n2m
         dvor(ic,jc)=fj(jc,ic)
      enddo
      enddo
c
c  set the boundary conditions
c  
      if(inbcsn.eq.1) then
      do ic=2,n1m
      dvor(ic,n2)=dsnor(npsc,ic)
      enddo
                      endif
      if(inbcsn.eq.0) then
      do ic=2,n1m
      dvor(ic,n2)=dvor(ic,n2m)+dsnor(npsc,ic)
      enddo
                      endif
      if(inbcss.eq.1.or.inbcss.eq.-2) then
      do ic=2,n1m
      dvor(ic,1)=dssud(npsc,ic)
      enddo
                      endif
      if(inbcss.eq.0) then
      do ic=2,n1m
      dvor(ic,1)=dvor(ic,2)+dssud(npsc,ic)
      enddo
                      endif
      do jc=1,n2
        dvor(n1,jc)=dsest(npsc,jc)
        dvor(1,jc)=dswes(npsc,jc)
      enddo
      return
      end
c
c
c     ************************* subrout invtv **********************
c
c       this subroutine performs the computation of the
c       dvor=vor(n+1-vor(n
c
c     ***************************************************************
      subroutine invtv(dvor,al,ga,ro,vor,ru,visc,ekmn)
c    
      include 'param.f'
      dimension dvor(m1,m2),vor(m1,m2)
      dimension ru(m1,m2)
      dimension ami(md2),aci(md2),api(md2)
      dimension amj(md2),acj(md2),apj(md2),fj(md2,md2)
      common/dim/n1,n1m,n2,n2m
      common/tstep/dt
      common/indbo/imv(m1),ipv(m1)
      common/mesh/dx1,dx1q,dx2,dx2q
      common/metri/g2m(m1),g2c(m1)
      common/bchodv/dvsud(m1),dvnor(m1)
      common/bcvedv/dvest(m2),dvwes(m2)
      common/cwalv/inbcvs,inbcvn,inbcvw,inbcve
      m2v=m2
      m1v=m1
c
c     *********************************************
c       add to the non linear term the diffusive
c       derivatives 11 and 22
c     *********************************************
c
      do 1 ic=2,n1m
      ip=ic+1
      im=ic-1
      do 1 jc=2,n2m
      jm=jc-1
      jp=jc+1
      visder=dx1q*((vor(ip,jc)-vor(ic,jc))/g2m(ic)
     1            -(vor(ic,jc)-vor(im,jc))/g2m(ic-1))/g2c(ic)
     1      +dx2q*(vor(ic,jm)-2.*vor(ic,jc)+vor(ic,jp))
      rhs=(ga*dvor(ic,jc)+ro*ru(ic,jc))*dt +
     1     (visder*visc)*dt*al - ekmn*vor(ic,jc)*dt*al
      ru(ic,jc)=dvor(ic,jc)
      dvor(ic,jc)=rhs
    1 continue
c
c     *******************************************************
c          compute dvor sweeping along the x2 direction
c               wall boundaries
c     *******************************************************
      betady=dx2q*0.5*visc*dt*al
      ecco=dt*al*ekmn*0.25
      accy=1.+betady*2.+ecco
      betadx=dx1q*0.5*visc*al*dt
      do 20 jc=2,n2m
      jm=jc-1
      acj(jc)=accy
      amj(jc)=-betady
      apj(jc)=-betady
   20 continue
      if(inbcvn.eq.1) then
c
c   radiative
c
      acj(n2)=1.
      apj(n2)=0.
      amj(n2)=0.
      do ic=2,n1m
      dvor(ic,n2)=dvnor(ic)*(1.+ecco)-betadx*(
     1           (dvnor(ic+1)-dvnor(ic))/g2m(ic-1)
     1          -(dvnor(ic)-dvnor(ic-1))/g2m(ic))
      enddo
                      endif
      if(inbcvn.eq.0) then
c
c   free slip   omega=0
c
      acj(n2)=1.
      apj(n2)=0.
      amj(n2)=0.
      do ic=2,n1m
      dvor(ic,n2)=dvnor(ic)
      enddo
                      endif
      if(inbcvn.eq.-1.or.inbcvn.eq.2) then
c
c   free slip   d omega/dx2=0 or no-slip
c
      acj(n2)=1.
      apj(n2)=-1.
      amj(n2)=0.
      do ic=2,n1m
      dvor(ic,n2)=dvnor(ic)
      enddo
                      endif
      if(inbcvs.eq.1) then
c
c   radiative b.c.   sud
c
      acj(1)=1.
      apj(1)=0.
      amj(1)=0.
      do ic=2,n1m
      dvor(ic,1)=dvsud(ic)*(1.+ecco)-betadx*(
     1           (dvsud(ic+1)-dvsud(ic))/g2m(ic-1)
     1          -(dvsud(ic)-dvsud(ic-1))/g2m(ic))
      enddo
                      endif
      if(inbcvs.eq.-1.or.inbcvs.eq.2) then
c
c  no-slip b.c.   and zero flux
c
      acj(1)=1.
      apj(1)=-1.
      amj(1)=0.
      do ic=2,n1m
      dvor(ic,1)=dvsud(ic)
      enddo
                      endif
      if(inbcvs.eq.0.or.inbcvs.eq.-2) then
c
c  free-slip b.c.   sud  or inlet
c
      acj(1)=1.
      apj(1)=0.
      amj(1)=0.
      do ic=2,n1m
      dvor(ic,1)=dvsud(ic)
      enddo
                      endif
   22 continue
      dvoma=0.
      do ic=2,n1m
      do jc=1,n2
      fj(ic,jc)=dvor(ic,jc)
      dvoma=max(abs(dvor(ic,jc)),dvoma)
      enddo
      enddo
c!
      call tribk(amj,acj,apj,fj,n2,2,n1m)
c!
      dvoman=0.
      do ic=2,n1m
      do jc=1,n2
      dvor(ic,jc)=fj(ic,jc)
      dvoman=max(abs(dvor(ic,jc)),dvoman)
      enddo
      enddo
c
c     ***************************************************
c         compute dvor  sweeping in the x1 direction
c     ***************************************************      
      do ic=2,n1m
      aci(ic)=1.+ecco+betadx*(1./g2m(ic-1)+1./g2m(ic))/g2c(ic)
      ami(ic)=-betadx/g2m(ic-1)/g2c(ic)
      api(ic)=-betadx/g2m(ic)/g2c(ic)
      enddo
      ic=1
      if(inbcvw.eq.-1) then
c
c   no-slip
c
      aci(ic)=-1.
      api(ic)=1
      ami(ic)=0.
                       endif
      if(inbcvw.eq.0.or.inbcvw.eq.-2.or.inbcvw.eq.1) then
c
c   free-slip   or inlet or radiative
c
      aci(ic)=1.
      ami(ic)=0.
      api(ic)=0.
                       endif
      ic=n1
      if(inbcve.eq.-1) then
c
c   no-slip
c
      aci(ic)=1.
      ami(ic)=-1
      api(ic)=0.
                       endif
      if(inbcve.eq.1.or.inbcve.eq.0) then
c
c  free-slip b.c.   est  or  radiative
c
      aci(ic)=1.
      ami(ic)=0.
      api(ic)=0.
                       endif
      if(inbcve.eq.2) then
c
c   zero flux
c
      aci(ic)=1.
      ami(ic)=-1
      api(ic)=0.
                       endif
      do 21 jc=2,n2m
        fj(jc,1)=dvwes(jc)
        fj(jc,n1)=dvest(jc)
      do 21 ic=2,n1m
        fj(jc,ic)=dvor(ic,jc)
   21 continue     
c! 
      call tribk(ami,aci,api,fj,n1,2,n2m)
c!
      dvomal=0.
      do ic=1,n1
      do jc=2,n2m
         dvor(ic,jc)=fj(jc,ic)
      dvomal=max(abs(dvor(ic,jc)),dvomal)
      enddo
      enddo
      if(inbcvs.ne.-1) then
      do ic=1,n1
      dvor(ic,1)=dvsud(ic)
      enddo
                      endif
      if(inbcvs.eq.-1.or.inbcvs.eq.2) then
      do ic=1,n1
      dvor(ic,1)=dvor(ic,2)+dvsud(ic)
      enddo
                      endif
      if(inbcvn.ne.-1) then
      do ic=1,n1
      dvor(ic,n2)=dvnor(ic)
      enddo
                      endif
      if(inbcvn.eq.-1.or.inbcvn.eq.2) then
      do ic=1,n1
      dvor(ic,n2)=dvor(ic,n2m)+dvnor(ic)
      dvomal=max(abs(dvor(ic,jc)),dvomal)
      enddo
                      endif
      write(98,*)'dvor', dvoma,dvoman,dvomal
      return
      end
c
c
c  ****************************** subrout potvor **********************
c
c  in this subroutine the potential vorticity is calculated
c
      subroutine potvor(vor,vorq)
      include 'param.f'
      common/dim/n1,n1m,n2,n2m
      dimension vorq(m1,m2)
      dimension vor(m1,m2)
      common/betef/beta
      common/topog/htop(m1,m2)
      common/coor/yp1(m1),yp2(m2)
      common/indpe2/n2i,n2f
c
      vorma=0.
      voqma=0.
      do 2 ic=1,n1
      do 2 jc=1,n2f
      vorq(ic,jc)=vor(ic,jc)+beta*yp2(jc)+htop(ic,jc)
      vorma=max(abs(vor(ic,jc)),vorma)
      voqma=max(abs(vorq(ic,jc)),voqma)
    2 continue
c     write(6,*)' potvor',n2f,vorma,voqma
      return
      end
c
c
c  ****************************** subrout hdnlh **********************
c
c  in this subroutine are calculated the non-linear terms,
c  following the arakawa scheme.
c  The topography and the beta term are added to the vorticity
c  to evaluate the potential vorticity Jacobian in hdnl.
c
      subroutine hdnlh(vor,psi,h)
      include 'param.f'
      common/dim/n1,n1m,n2,n2m
      dimension h(m1,m2),vorq(m1,m2)
      dimension vor(m1,m2),psi(m1,m2)
      common/betef/beta
      common/topog/htop(m1,m2)
      common/coor/yp1(m1),yp2(m2)
      common/indpe2/n2i,n2f
c
      do 2 ic=1,n1 
      do 2 jc=1,n2f
      vorq(ic,jc)=vor(ic,jc)+beta*yp2(jc)+htop(ic,jc)
    2 continue  
      call hdnl(vorq,psi,h)
      return
      end
c
c
c  ****************************** subrout hdnl  **********************
c
c  in this subroutine are calculated the non-linear terms,
c  following the arakawa scheme.
c
      subroutine hdnl(vor,psi,h)
      include 'param.f'
      common/dim/n1,n1m,n2,n2m
      common/mesh/dx1,dx1q,dx2,dx2q
      common/metri/g2m(m1),g2c(m1)
      dimension h(m1,m2)
      dimension vor(m1,m2),psi(m1,m2)
      common/indbo/imv(m1),ipv(m1)
      common/chnlc/chal,chbe,chga
      common/indpe2/n2i,n2f
      common/indwa/jmv(m2),jpv(m2)
      common/betef/beta
c
c  **********  compute the non-linear terms by different
c   schemes
c
      udx12=dx2*dx1*0.25
      do 10 ic=2,n1m
      ip=ic+1
      im=ic-1
      udn12=udx12/g2c(ic)
      do 10 jc=n2i,n2m
      jm=jmv(jc)
      jp=jpv(jc)
c
c     terms as Arakawa sect 4.6
c
      h22a=((vor(ip,jc)-vor(im,jc))*(psi(ic,jp)-psi(ic,jm))-
     1      (vor(ic,jp)-vor(ic,jm))*(psi(ip,jc)-psi(im,jc)))
     1     *udn12
      h22b=(vor(ip,jc)*(psi(ip,jp)-psi(ip,jm))-
     1      vor(im,jc)*(psi(im,jp)-psi(im,jm))-
     1      vor(ic,jp)*(psi(ip,jp)-psi(im,jp))+
     1      vor(ic,jm)*(psi(ip,jm)-psi(im,jm)))
     1     *udn12
      h22c=(vor(ip,jp)*(psi(ic,jp)-psi(ip,jc))-
     1      vor(im,jm)*(psi(im,jc)-psi(ic,jm))-
     1      vor(im,jp)*(psi(ic,jp)-psi(im,jc))+
     1      vor(ip,jm)*(psi(ip,jc)-psi(ic,jm)))
     1     *udn12
      hq2=chal*h22a+chbe*h22b+chga*h22c
      h(ic,jc)=-hq2
   10 continue
      return
      end
c
c  ****************************** subrout invpps **********************
c
c   this subroutine performs the computation of the 
c   passive scalar dpsc=psc(n+1-psc(n
c   in the case of ib2per=1 that accounts for periodicity
c   along x2
c
c
      subroutine invpps(dvor,al,ga,ro,vor,ru,visc,npsc)    
      include 'param.f'
      dimension dvor(m1,m2),vor(m1,m2)
      dimension ru(mpsc,m1,m2)
      dimension ami(md2),aci(md2),api(md2)
      dimension amj(md2),acj(md2),apj(md2),fj(md2,md2)
      common/dim/n1,n1m,n2,n2m
      common/tstep/dt
      common/pscqu/sch(mpsc),scla,npscf
      common/mesh/dx1,dx1q,dx2,dx2q
      common/metri/g2m(m1),g2c(m1)
      common/bcveds/dsest(mpsc,m2),dswes(mpsc,m2)
      common/cwals/inbcss,inbcsn,inbcsw,inbcse
      common/indwa/jmv(m2),jpv(m2)
      m2v=m2
      m1v=m1
c
c   add to the non linear term the diffusive derivatives
c   11 1nd 22
c
      do 1 jc=1,n2m
      jm=jmv(jc)
      jp=jpv(jc)
      do 1 ic=2,n1m
      ip=ic+1
      im=ic-1
      visder=dx1q*((vor(ip,jc)-vor(ic,jc))/g2m(ic)
     1            -(vor(ic,jc)-vor(im,jc))/g2m(ic-1))/g2c(ic)
     1      +dx2q*(vor(ic,jm)-2.*vor(ic,jc)+vor(ic,jp))
      rhs=(ga*dvor(ic,jc)+ro*ru(npsc,ic,jc))*dt +
     1     (visder*visc)*dt*al 
      ru(npsc,ic,jc)=dvor(ic,jc)
      dvor(ic,jc)=rhs
    1 continue
c
c  ************ compute dvor sweeping along the x2 direction
c               periodicity     
c
      betady=dx2q*0.5*visc*dt*al
      accy=1.+betady*2.
      betadx=dx1q*0.5*visc*al*dt
      do 20 jc=1,n2m
      acj(jc)=accy
      amj(jc)=-betady
      apj(jc)=-betady
   20 continue
      do ic=2,n1m
      do jc=1,n2m
      fj(ic,jc)=dvor(ic,jc)
      enddo
      enddo
      call tripv(amj,acj,apj,fj,n2m,2,n1m)
      do ic=2,n1m
      do jc=1,n2m
      dvor(ic,jc)=fj(ic,jc)
      enddo
      enddo
c
c  ********* compute dvor  sweeping in the x1 direction
c            more general boundary conditions
c           
      do 3 ic=2,n1m
      aci(ic)=1.+betadx*(1./g2m(ic-1)+1./g2m(ic))/g2c(ic)
      ami(ic)=-betadx/g2m(ic-1)/g2c(ic)
      api(ic)=-betadx/g2m(ic)/g2c(ic)
    3 continue
      ic=1
      if(inbcsw.eq.0) then
      aci(ic)=1.
      ami(ic)=0.
      api(ic)=-1.
                      else
      aci(ic)=1.
      ami(ic)=0.
      api(ic)=0.
                      endif
      ic=n1
      if(inbcse.ge.1) then
      aci(ic)=1.
      ami(ic)=0.
      api(ic)=0.
                      else
      aci(ic)=1.
      ami(ic)=-1.
      api(ic)=0.
                      endif
      do 21 jc=1,n2m
        fj(jc,1)=dswes(npsc,jc)
        fj(jc,n1)=dsest(npsc,jc)
      do 21 ic=2,n1m
        fj(jc,ic)=dvor(ic,jc)
   21 continue     
c
      call tribk(ami,aci,api,fj,n1,1,n2m)
c
      do ic=2,n1m
      do jc=1,n2m
         dvor(ic,jc)=fj(jc,ic)
      enddo
      enddo
      if(inbcse.ge.1) then
      do jc=1,n2m
        dvor(n1,jc)=dsest(npsc,jc)
      enddo
                      else
      do jc=1,n2m
        dvor(n1,jc)=dvor(n1m,jc)+dsest(npsc,jc)
      enddo
                      endif
      if(inbcsw.eq.0) then
      do jc=1,n2m
        dvor(1,jc)=dvor(2,jc)+dswes(npsc,jc)
      enddo
                      else
      do jc=1,n2m
        dvor(1,jc)=dswes(npsc,jc)
      enddo
                      endif
      return
      end
c
c  ****************************** subrout invpv **********************
c
c   this subroutine performs the computation of the dvor=vor(n+1-vor(n
c    this routine is called when ib2per=1 accounting for periodicity
c   in x2
c
c
      subroutine invpv(dvor,al,ga,ro,vor,ru,visc,ekmn)    
      include 'param.f'
      dimension dvor(m1,m2),vor(m1,m2)
      dimension ru(m1,m2)
      dimension ami(md2),aci(md2),api(md2)
      dimension amj(md2),acj(md2),apj(md2),fj(md2,md2)
      common/dim/n1,n1m,n2,n2m
      common/tstep/dt
      common/indbo/imv(m1),ipv(m1)
      common/mesh/dx1,dx1q,dx2,dx2q
      common/metri/g2m(m1),g2c(m1)
      common/bcvedv/dvest(m2),dvwes(m2)
      common/cwalv/inbcvs,inbcvn,inbcvw,inbcve
      common/indwa/jmv(m2),jpv(m2)
      common/indpe2/n2i,n2f
      m2v=m2
      m1v=m1
c
c   add to the non linear term the diffusive derivatives
c   11 1nd 22
c
c     write(58,134)dt,al,visc,ekmn
  134 format(3x,10e12.5)
c     ipr=n1m/2+1
c     jpr=n2m
      do 1 ic=2,n1m
      ip=ic+1
      im=ic-1
      do 1 jc=1,n2m
      jm=jmv(jc)
      jp=jpv(jc)
      visder=dx1q*((vor(ip,jc)-vor(ic,jc))/g2m(ic)
     1            -(vor(ic,jc)-vor(im,jc))/g2m(ic-1))/g2c(ic)
     1      +dx2q*(vor(ic,jm)-2.*vor(ic,jc)+vor(ic,jp))
      rhs=(ga*dvor(ic,jc)+ro*ru(ic,jc))*dt +
     1     (visder*visc)*dt*al - ekmn*vor(ic,jc)*dt*al
      ru(ic,jc)=dvor(ic,jc)
      dvor(ic,jc)=rhs
    1 continue

c
c  ************ compute dvor sweeping along the x2 direction
c               periodicity    
c
      betady=dx2q*0.5*visc*dt*al
      ecco=dt*al*ekmn*0.25
      accy=1.+betady*2.+ecco
      betadx=dx1q*0.5*visc*al*dt
      do 20 jc=1,n2m
      acj(jc)=accy
      amj(jc)=-betady
      apj(jc)=-betady
   20 continue
      do ic=2,n1m
      do jc=1,n2m
      fj(ic,jc)=dvor(ic,jc)
      enddo
      enddo
      call tripv(amj,acj,apj,fj,n2m,2,n1m)
      do ic=2,n1m
      do jc=1,n2m
      dvor(ic,jc)=fj(ic,jc)
      enddo
      enddo
c     dvojpr=dvor(ipr,jpr)
c
c  ********* compute dvor  sweeping in the x1 direction
c           
      do ic=2,n1m
      aci(ic)=1.+ecco+betadx*(1./g2m(ic-1)+1./g2m(ic))/g2c(ic)
      ami(ic)=-betadx/g2m(ic-1)/g2c(ic)
      api(ic)=-betadx/g2m(ic)/g2c(ic)
      enddo
      ic=1
      aci(ic)=1.
      ami(ic)=0.
      api(ic)=0.
      ic=n1
      if(inbcve.eq.-1) then
      aci(ic)=1.
      ami(ic)=-1
      api(ic)=0.
                       else
      aci(ic)=1.
      ami(ic)=0.
      api(ic)=0.
                       endif
      do 21 jc=n2i,n2m
        fj(jc,1)=dvwes(jc)
        fj(jc,n1)=dvest(jc)
      do 21 ic=2,n1m
        fj(jc,ic)=dvor(ic,jc)
   21 continue     
c
      call tribk(ami,aci,api,fj,n1,1,n2m)
c
      do ic=1,n1
      do jc=1,n2m
         dvor(ic,jc)=fj(jc,ic)
      enddo
      enddo
c     dvoipr=dvor(ipr,jpr)
c     write(57,151),rhopr,dvojpr,dvoipr
c 151 format(3x,' in ivtpv',5x,5e12.5)
      return
      end
