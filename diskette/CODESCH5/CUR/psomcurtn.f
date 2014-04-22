c
c  ****************************** subrout hdnlvi **********************
c
c   in this subroutine  the cross derivatives of the viscous terms
c   are added  to the nonlinear terms
c
c
      subroutine hdnlvi(vor,dvor,visc)
      include 'param.f'
      dimension dvor(m1,m2),vor(m1,m2)
      common/dim/n1,n1m,n2,n2m
      common/tstep/dt
      common/cnlev/n1g(0:nlev),n2g(0:nlev),n12g(0:nlev)
      common/comg/co(10,nij)
      igjg(i,j,l)=i+(j-1)*n1g(l)+n12g(l-1)
      do  ic=2,n1m
      ip=ic+1
      im=ic-1
      do  jc=2,n2m
      jm=jc-1
      jp=jc+1
      ij=igjg(ic,jc,1)
      visdcr= +co(1,ij)*vor(ic+1,jc+1)+
     1         co(3,ij)*vor(ic+1,jc-1)+
     1         co(5,ij)*vor(ic-1,jc-1)+
     1         co(7,ij)*vor(ic-1,jc+1)
      dvor(ic,jc)=dvor(ic,jc)+visdcr*visc*dt  
      enddo
      enddo
      return
      end
c
c  ****************************** subrout invtri **********************
c
c   this subroutine performs the computation of the dvor=vor(n+1-vor(n
c   in the inviscid case
c
      subroutine invtri(dvor,al,ga,ro,vor)
      include 'param.f'
      dimension dvor(m1,m2),vor(m1,m2)
      common/sliwal/insls,insln
      common/dvowal/dvorbs(m1),dvorbn(m1)
      common/bdqout/dvest(m2),dpest(m2),dpsest(m2)
      common/bdqinl/dvwes(m2),dpwes(m2),dpswes(m2)
      common/dim/n1,n1m,n2,n2m
      common/nonlt/rt(m1,m2),ru(m1,m2)
      common/tstep/dt
      common/reyn/beta,re
      common/cnlev/n1g(0:nlev),n2g(0:nlev),n12g(0:nlev)
      do ic=2,n1m
      ip=ic+1
      im=ic-1
      do jc=2,n2m
      jm=jc-1
      jp=jc+1
c
c  only the 11 and 22 derivatives
c
      rhs=(ga*dvor(ic,jc)+ro*ru(ic,jc))*dt
      ru(ic,jc)=dvor(ic,jc)
      dvor(ic,jc)=rhs
      enddo
      enddo
      do  jc=2,n2m
          dvor(jc,1)=dvwes(jc)
          dvor(jc,n1)=dvest(jc)
      enddo
      do ic=2,n1m
          dvor(ic,1)=dvorbs(ic)
      enddo
      jc=n2
      do ic=2,n1m
          dvor(ic,n2)=dvorbn(ic)
      enddo
      do j=1,n2
      do i=2,n1m
      vor(i,j)=dvor(i,j)+vor(i,j)
      enddo
      enddo
      do j=1,n2
      vor(1,j)=dvor(1,j)+vor(1,j)
      vor(n1,j)=dvor(n1,j)+vor(n1,j)
      enddo
      return
      end
c
c  ****************************** subrout invtrv **********************
c
c   this subroutine performs the computation of the dvor=vor(n+1-vor(n
c   in the viscous case
c
      subroutine invtrv(dvor,al,ga,ro,vor,visc)
      include 'param.f'
      common/coetrb/am(md1,md1),ac(md1,md1),ap(md1,md1)
     1       ,fi(md1,md1)
      dimension dvor(m1,m2),vor(m1,m2)
      common/sliwal/insls,insln
      common/dvowal/dvorbs(m1),dvorbn(m1)
      common/bdqout/dvest(m2),dpest(m2),dpsest(m2)
      common/bdqinl/dvwes(m2),dpwes(m2),dpswes(m2)
      common/dim/n1,n1m,n2,n2m
      common/nonlt/rt(m1,m2),ru(m1,m2)
      common/tstep/dt
      common/reyn/beta,re
      common/cnlev/n1g(0:nlev),n2g(0:nlev),n12g(0:nlev)
      common/comg/co(10,nij)
      igjg(i,j,l)=i+(j-1)*n1g(l)+n12g(l-1)
      do ic=2,n1m
      ip=ic+1
      im=ic-1
      do jc=2,n2m
      jm=jc-1
      jp=jc+1
      ij=igjg(ic,jc,1)
c
c  only the 11 and 22 derivatives are inserted in the RHS
c
c  the cross derivatives 12 and 21 in hdnlv
c
      visder=+ co(2,ij)*vor(ic+1,jc)+
     1         co(4,ij)*vor(ic,jc-1)+
     1         co(6,ij)*vor(ic-1,jc)+
     1         co(8,ij)*vor(ic,jc+1)+
     1         co(9,ij)*vor(ic,jc)
      rhs=(ga*dvor(ic,jc)+ro*ru(ic,jc))*dt
     1    +visder*dt*al*visc
      ru(ic,jc)=dvor(ic,jc)
      dvor(ic,jc)=rhs
      enddo
      enddo
      betadx=al*visc*dt*0.5
c
c    operator in the 1 direct. for vorticity  equations
c
      do ic=2,n1m
      do jc=2,n2m
      ij=igjg(ic,jc,1)
          ac(jc,ic)=1.-co(9,ij)*betadx
          am(jc,ic)=-co(6,ij)*betadx
          ap(jc,ic)=-co(2,ij)*betadx
      enddo
      enddo
c
c    boundary conditions  fre-slip or radiative
c
      do jc=2,n2m
      ic=1
          am(jc,ic)=0.
          ac(jc,ic)=1.
          ap(jc,ic)=0.
      ic=n1
          ap(jc,ic)=0.
          ac(jc,ic)=1.
          am(jc,ic)=0.
      enddo
      do ic=2,n1m
      do jc=2,n2m
           fi(jc,ic)=dvor(ic,jc)
      enddo
      enddo
      do  jc=2,n2m
          fi(jc,1)=dvwes(jc)
          fi(jc,n1)=dvest(jc)
      enddo
      n=n1
      mi=2
      mf=n2m
      call tribk(n,mi,mf)
      do jc=2,n2m
      do ic=1,n1
      dvor(ic,jc)=fi(jc,ic)
      enddo
      enddo

c
c    operator in the 2 direct. for vorticity  equations
c
      do ic=2,n1m
      do jc=2,n2m
      ij=igjg(ic,jc,1)
          ac(ic,jc)=1.-co(9,ij)*betadx
          am(ic,jc)=-co(4,ij)*betadx
          ap(ic,jc)=-co(8,ij)*betadx
      enddo
      enddo
c
c    boundary conditions no-slip or free-slip for
c    the south boundary
c
      jc=1
      if(insls.eq.0) then
      do ic=2,n1m
          am(ic,jc)=0.
          ac(ic,jc)=1.
          ap(ic,jc)=0.
          fi(ic,1)=dvorbs(ic)
      enddo
                     else
      do ic=2,n1m
          am(ic,jc)=0.
          ac(ic,jc)=1.
          ap(ic,jc)=-1.
          fi(ic,1)=dvorbs(ic)
      enddo
                     endif
      jc=n2
      if(insln.eq.0) then
c
c  free-slip
c
      do ic=2,n1m
          ap(ic,jc)=0.
          ac(ic,jc)=1.
          am(ic,jc)=0.
          fi(ic,n2)=dvorbn(ic)
      enddo
                     endif
      if(insln.eq.1) then
c
c  no-slip
c
      do ic=2,n1m
          ap(ic,jc)=0.
          ac(ic,jc)=1.
          am(ic,jc)=-1.
          fi(ic,n2)=dvorbn(ic)
      enddo
                     endif
      if(insln.eq.2) then
c
c  radiative only on th enorth boundary
c
      do ic=2,n1m
          ap(ic,jc)=0.
          ac(ic,jc)=1.
          am(ic,jc)=0.
          fi(ic,n2)=dvorbn(ic)
      enddo
                     endif
      do ic=2,n1m
      do jc=2,n2m
           fi(ic,jc)=dvor(ic,jc)
      enddo
      enddo
      n=n2
      mi=2 
      mf=n1m
      call tribk(n,mi,mf)
      do j=1,n2
      do i=2,n1m
      vor(i,j)=fi(i,j)+vor(i,j)
      enddo
      enddo
      do j=1,n2
      vor(1,j)=dvor(1,j)+vor(1,j)
      vor(n1,j)=dvor(n1,j)+vor(n1,j)
      enddo
      return
      end
c
c  ****************************** subrout invtpi **********************
c
c   this subroutine performs the computation of the dpsc=psc(n+1-psc(n
c    in the inviscid case
c
      subroutine invtpi(dvor,al,ga,ro,psc)
      include 'param.f'
      dimension dvor(m1,m2),psc(m1,m2)
      common/dpswal/dpscbs(m1),dpscbn(m1)
      common/sliwal/insls,insln
      common/lewai/xl1,xl2
      common/slewa/xls(m1)
      common/bdqout/dvest(m2),dpest(m2),dpsest(m2)
      common/bdqinl/dvwes(m2),dpwes(m2),dpswes(m2)
      common/dim/n1,n1m,n2,n2m
      common/nonlt/rt(m1,m2),ru(m1,m2)
      common/tstep/dt
      do ic=2,n1m
      ip=ic+1
      im=ic-1
      do jc=2,n2m
      jm=jc-1
      jp=jc+1
      rhs=(ga*dvor(ic,jc)+ro*rt(ic,jc))*dt
      rt(ic,jc)=dvor(ic,jc)
      dvor(ic,jc)=rhs
      enddo
      enddo
      do  jc=2,n2m
          dvor(1,jc)=dpswes(jc)
          dvor(n1,jc)=dpsest(jc)
      enddo
      jc=1
      do ic=2,n1m
          dvor(ic,1)=dpscbs(ic)
      enddo
      jc=n2
      do ic=2,n1m
          dvor(ic,n2)=dpscbn(ic)
      enddo
      do i=2,n1m
      do j=1,n2
      psc(i,j)=dvor(i,j)+psc(i,j)
      enddo
      enddo
      do j=1,n2
      psc(1,j)=psc(1,j)+dvor(1,j)
      psc(n1,j)=dvor(n1,j)+psc(n1,j)
      enddo
      return
      end
c
c  ****************************** subrout invtps **********************
c
c   this subroutine performs the computation of the dpsc=psc(n+1-psc(n
c   in the viscous case
c
      subroutine invtps(dvor,al,ga,ro,psc,visc)
      include 'param.f'
      common/coetrb/am(md1,md1),ac(md1,md1),ap(md1,md1)
     1       ,fi(md1,md1)
      dimension dvor(m1,m2),psc(m1,m2)
      common/dpswal/dpscbs(m1),dpscbn(m1)
      common/sliwal/insls,insln
      common/lewai/xl1,xl2
      common/slewa/xls(m1)
      common/bdqout/dvest(m2),dpest(m2),dpsest(m2)
      common/bdqinl/dvwes(m2),dpwes(m2),dpswes(m2)
      common/dim/n1,n1m,n2,n2m
      common/nonlt/rt(m1,m2),ru(m1,m2)
      common/tstep/dt
      common/cnlev/n1g(0:nlev),n2g(0:nlev),n12g(0:nlev)
      common/comg/co(10,nij)
      igjg(i,j,l)=i+(j-1)*n1g(l)+n12g(l-1)
      do ic=2,n1m
      ip=ic+1
      im=ic-1
      do jc=2,n2m
      jm=jc-1
      jp=jc+1
      ij=igjg(ic,jc,1)
c
c  only the 11 and 22 derivatives
c
c  the cross derivatives 12 and 21 in hdnlv
c
      visder=+ co(2,ij)*psc(ic+1,jc)+
     1         co(4,ij)*psc(ic,jc-1)+
     1         co(6,ij)*psc(ic-1,jc)+
     1         co(8,ij)*psc(ic,jc+1)+
     1         co(9,ij)*psc(ic,jc)
      rhs=(ga*dvor(ic,jc)+ro*rt(ic,jc))*dt
     1    +visder*dt*al*visc
      rt(ic,jc)=dvor(ic,jc)
      dvor(ic,jc)=rhs
      enddo
      enddo
      betadx=al*visc*dt*0.5
c
c    operator in the 1 direct. for vorticity  equations
c
      do ic=2,n1m
      do jc=2,n2m
      ij=igjg(ic,jc,1)
          ac(jc,ic)=1.-co(9,ij)*betadx
          am(jc,ic)=-co(6,ij)*betadx
          ap(jc,ic)=-co(2,ij)*betadx
      enddo
      enddo
c
c    boundary conditions on the west east boundaries
c
      do jc=2,n2m
      ic=1
          am(jc,ic)=0.
          ac(jc,ic)=1.
          ap(jc,ic)=0.
      ic=n1
          ap(jc,ic)=0.
          ac(jc,ic)=1.
          am(jc,ic)=0.
      enddo
      do ic=2,n1m
      do jc=2,n2m
           fi(jc,ic)=dvor(ic,jc)
      enddo
      enddo
      do  jc=2,n2m
          fi(jc,1)=dpswes(jc)
          fi(jc,n1)=dpsest(jc)
      enddo
      n=n1
      mi=2
      mf=n2m
      call tribk(n,mi,mf)
      do jc=2,n2m
      do ic=1,n1
      dvor(ic,jc)=fi(jc,ic)
      enddo
      enddo

c
c    operator in the 2 direct. for vorticity  equations
c
      do ic=1,n1
      do jc=2,n2m
      ij=igjg(ic,jc,1)
          ac(ic,jc)=1.-co(9,ij)*betadx
          am(ic,jc)=-co(4,ij)*betadx
          ap(ic,jc)=-co(8,ij)*betadx
      enddo
      enddo
c
c   conditions on the south boundary 
c   the possibilities of a thermal discharge is considered
c
      jc=1
      do ic=2,n1m
      if(xls(ic).ge.xl1.and.xls(ic).le.xl2) then
          am(ic,jc)=0.
          ac(ic,jc)=1.
          ap(ic,jc)=0.
                     else
          am(ic,jc)=0.
          ac(ic,jc)=1.
          ap(ic,jc)=-1.
                     endif
          fi(ic,1)=dpscbs(ic)
      enddo
      jc=n2
      if(insln.lt.2) then
c
c  zero flux
c
      do ic=2,n1m
          ap(ic,jc)=0.
          ac(ic,jc)=1.
          am(ic,jc)=-1.
          fi(ic,n2)=dpscbn(ic)
      enddo
                     endif
      if(insln.lt.2) then
c
c  radiative 
c
      do ic=2,n1m
          ap(ic,jc)=0.
          ac(ic,jc)=1.
          am(ic,jc)=0.
          fi(ic,n2)=dpscbn(ic)
      enddo
                     endif
      do ic=2,n1m
      do jc=2,n2m
           fi(ic,jc)=dvor(ic,jc)
      enddo
      enddo
      n=n2
      mi=2
      mf=n1m
      call tribk(n,mi,mf)
      do i=2,n1m
      do j=1,n2
      psc(i,j)=fi(i,j)+psc(i,j)
      enddo
      enddo
      do j=1,n2
      psc(1,j)=psc(1,j)+dvor(1,j)
      psc(n1,j)=dvor(n1,j)+psc(n1,j)
      enddo
      return
      end
c
c  ****************************** subrout tribk  **********************
c
      subroutine trib(n,mi,mf,am,ac,ap,f)
      include 'param.f'
      dimension am(md1,md1),ac(md1,md1),ap(md1,md1)
     1       ,f(md1,md1)
c
c  ******** reduction of trid. matrix to an upper rigth matrix
c
      do 1 i=2,n
      do 2 k=mi,mf
      ac(k,i)=ac(k,i)*ac(k,i-1)-ap(k,i-1)*am(k,i)
      ap(k,i)=ap(k,i)*ac(k,i-1)
      f(k,i)=f(k,i)*ac(k,i-1)-f(k,i-1)*am(k,i)
    2 continue
    1 continue
c  ******** calculation of the unknown by backward elimination
c
      do 3 k=mi,mf
      f(k,n)=f(k,n)/ac(k,n)
    3 continue
      nm=n-1
      do 10 ii=1,nm
      i=n-ii
      do 11 k=mi,mf
      f(k,i)=(f(k,i)-ap(k,i)*f(k,i+1))/ac(k,i)
   11 continue
   10 continue
      return
      end
c                                                                       *
c************************************************************************
c
c     ROUTINE FOR THE INVERSION OF BOUNDED TRIDIAGONAL MATRICES
c
      subroutine tribk(n,mi,mf)
      include 'param.f'
      common/coetrb/am(md1,md1),ac(md1,md1),ap(md1,md1)
     1       ,fj(md1,md1)
      dimension gmj(md1,md1),bj(md1)
      do 10 i=mi,mf
      bj(i)=ac(i,1)
      fj(i,1)=fj(i,1)/bj(i)
   10 continue
      do 11 j=2,n
      do 21 i=mi,mf
      gmj(i,j)=ap(i,j-1)/bj(i)
      bj(i)=ac(i,j)-am(i,j)*gmj(i,j)
      fj(i,j)=(fj(i,j)-am(i,j)*fj(i,j-1))/bj(i)
   21 continue
   11 continue
      do 12 j=n-1,1,-1
      do 22 i=mi,mf
      fj(i,j)=fj(i,j)-gmj(i,j+1)*fj(i,j+1)
   22 continue
   12 continue
      return
      end
c
c
c  ****************************** subrout hdnl  **********************
c
c  in this subroutine are calculated the non-linear terms,
c  following the Arakawa scheme that is valid even for 
c   general curvilinear coordinates
c
      subroutine hdnl(vor,psi,h)
      include 'param.f'
      common/dim/n1,n1m,n2,n2m
      common/mesh/dx1,dx1q,dx2,dx2q
      dimension vor(m1,m2),psi(m1,m2),h(m1,m2)
      common/chnlc/chal,chbe,chga
      common/metrst/gccc(m1,m2)
      common/metrih/cac(ndd,ndd,m1,m2)
c
c  **********  compute the non-linear terms by different
c   schemes
c
      do 10 ic=2,n1m
      ip=ic+1
      im=ic-1
      do 10 jc=2,n2m
      jm=jc-1
      jp=jc+1
c
c     termsas arakawa pg 129
c
      h22a=((vor(ip,jc)-vor(im,jc))*(psi(ic,jp)-psi(ic,jm))-
     1      (vor(ic,jp)-vor(ic,jm))*(psi(ip,jc)-psi(im,jc)))
     1     *dx2/gccc(ic,jc)*dx1*0.25
      h22b=(vor(ip,jc)*(psi(ip,jp)-psi(ip,jm))-
     1      vor(im,jc)*(psi(im,jp)-psi(im,jm))-
     1      vor(ic,jp)*(psi(ip,jp)-psi(im,jp))+
     1      vor(ic,jm)*(psi(ip,jm)-psi(im,jm)))
     1     *dx2/gccc(ic,jc)*dx1*0.25
      h22c=(vor(ip,jp)*(psi(ic,jp)-psi(ip,jc))-
     1      vor(im,jm)*(psi(im,jc)-psi(ic,jm))-
     1      vor(im,jp)*(psi(ic,jp)-psi(im,jc))+
     1      vor(ip,jm)*(psi(ip,jc)-psi(ic,jm)))
     1     *dx2/gccc(ic,jc)*dx1*0.25
      hq2=chal*h22a+chbe*h22b+chga*h22c
      h(ic,jc)=-hq2
   10 continue
      return
      end
