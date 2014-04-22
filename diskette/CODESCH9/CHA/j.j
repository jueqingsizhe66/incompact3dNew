c
c  ****************************** subrout invtr3  **********************
c   this subroutine performs the inversion of the q3 momentum equation
c   by a factored implicit scheme, only the derivatives 11,22,33 of q3
c   are treated implicitly
c
      subroutine invtr3(q,dq,al,ga,ro,pr,ru)
c

      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      common/amcpj3/amj3(m2),acj3(m2),apj3(m2)
      common/trici/ami(m2,m1),aci(m2,m1),api(m2,m1),
     1           fi(m2,m1),fei(m2,m1),qq(m2,m1),ss(m2,m1)
      common/tvbkj/amj(m3,m2),acj(m3,m2),apj(m3,m2),fj(m3,m2)
      common/tvpjk/amk(m2,m3),ack(m2,m3),apk(m2,m3),fk(m2,m3)
     1         ,qek(m2,m3),qk(m2,m3),sk(m2,m3)
      dimension ru(ndv,m1,m2,m3)
      common/rhsc/rhs(m1,m2,m3)
      dimension q(ndv,m1,m2,m3),pr(m1,m2,m3)
      dimension dq(m1,m2,m3)
      common/metria/caj(m2),cac(m2)
      common/tstep/dt,beta,ren
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/rhs3p/dp3ns
      common/neno/dpnew
      common/s3t/s3tot
c
      common/qwallo/q1s(m1,m3),q2s(m1,m3),q3s(m1,m3)
      common/qwalup/q1n(m1,m3),q2n(m1,m3),q3n(m1,m3)
      common/dqwalo/dq1s(m1,m3),dq2s(m1,m3),dq3s(m1,m3)
      common/dqwalu/dq1n(m1,m3),dq2n(m1,m3),dq3n(m1,m3)
      common/islwal/islv1s,islv1n,islv3s,islv3n
c
c
c  ********* compute the rhs of the factored equation
c  the rhs include h(n), h(n-1)=ru, grad(pr) component
c  and the 11, 22,3 second derivat. of q3(n)
c  everything at i,j+1/2,k=1/2
c
      alre=al/ren
      s3tot=0.
      n2mm=n2m-1
      do kc=1,n3m
      km=kmv(kc)
      kp=kpv(kc)
      do jc=1,n2m
      do ic=1,n1m
      im=imv(ic)
      ip=ipv(ic)
c
c   11 second derivatives of q3
c
      dq31=(q(3,ip,jc,kc)-2.*q(3,ic,jc,kc)+q(3,im,jc,kc))*dx1q
c
c   33 second derivatives of q3
c
      dq33=(q(3,ic,jc,kp)-2.*q(3,ic,jc,kc)+q(3,ic,jc,km))*dx3q
      dcq3=dq31+dq33
      rhs(ic,jc,kc)=(ga*dq(ic,jc,kc)+ro*ru(3,ic,jc,kc)+
     1               alre*dcq3)*dt
      s3tot=s3tot+dcq3*caj(jc)/ren
      ru(3,ic,jc,kc)=dq(ic,jc,kc) 
      enddo
      enddo
      enddo
c
c  add 22 second derivatives of q3 inner field
c
      do kc=1,n3m
      do jc=2,n2mm
      jmm=jmv(jc)
      jpp=jpv(jc)
      do ic=1,n1m
c
c   22 second derivatives of q3
c
      dq32=(apj3(jc)*q(3,ic,jpp,kc)
     1     -acj3(jc)*q(3,ic,jc,kc)
     1     +amj3(jc)*q(3,ic,jmm,kc))*dx2q
      rhs(ic,jc,kc)=rhs(ic,jc,kc)+alre*dq32*dt
      s3tot=s3tot+dq32*caj(jc)/ren
      enddo
      enddo
      enddo
c
c  add 22 second derivatives of q3 lower wall   
c
      jc=1
      jp=jc+1
      do kc=1,n3m
      do ic=1,n1m
c
c   22 second deriv. of q3(n) O(dx2^2)
c
      am=4./(2.*cac(2)+caj(1))*2./cac(jc)
      dq32=(apj3(jc)*q(3,ic,jp,kc)
     1     -acj3(jc)*q(3,ic,jc,kc)
     1     +am*q3s(ic,kc)*(1-islv3s) )*dx2q
      s3tot=s3tot+dq32*caj(jc)/ren
      rhs(ic,jc,kc)=rhs(ic,jc,kc)+alre*dq32*dt
      enddo
      enddo
c
c   add the second derivative d22q3 aqt the upper wall
c
      jc=n2m
      jm=jc-1
      do kc=1,n3m
      do ic=1,n1m
c
c   22 second deriv. of q3(n) O(dx2^2)
c
      ap=4./(2.*cac(n2m)+caj(n2m))*2./cac(n2)
      dq32=(amj3(jc)*q(3,ic,jm,kc)
     1     -acj3(jc)*q(3,ic,jc,kc)
     1     +ap*q3n(ic,kc)*(1-islv3n) )*dx2q
      s3tot=s3tot+dq32*caj(jc)/ren
      rhs(ic,jc,kc)=rhs(ic,jc,kc)+alre*dq32*dt
      enddo
      enddo
      dp3ns=s3tot/(2.*n1m*n2m*n3m)
c     print *,'s3tot=',s3tot
c     print *,'dp3ns=',dp3ns
c
c     print *,'correzione dp3ns =',dp3ns
c     print *,'correzione dp3nsa=',dp3nsa
c     print *,'       differenza=',dp3nsa-dp3ns
c
      do 38 kc=1,n3m
      km=kmv(kc)
      do 38 jc=1,n2m
      do 38 ic=1,n1m
c
c  component of grad(pr) along x3 direction
c
      dpx33=(pr(ic,jc,kc)-pr(ic,jc,km))*dx3
      gradp=(dpx33+dp3ns)*al
      rhs(ic,jc,kc)=rhs(ic,jc,kc)-gradp*dt
   38 continue
c
c  ********* compute the dq3* sweeping in the x3 direction
c
      betadz=dx3q*al*beta
      do 1 ic=1,n1m
      do 15 kc=1,n3m
      do 15 jc=1,n2m
      apk(jc,kc)=-betadz
      ack(jc,kc)=1.+betadz*2.
      amk(jc,kc)=-betadz
      fk(jc,kc)=rhs(ic,jc,kc)
   15 continue
c
      call trvpjk(1,n3m,1,n2m)
c
      do 3 kc=1,n3m
      do 3 jc=1,n2m
      rhs(ic,jc,kc)=fk(jc,kc)
    3 continue
    1 continue
c
c  ************ compute dq3** sweeping along the x1 direction
c
      betadx=dx1q*al*beta
      do 4 kc=1,n3m
      do 14 jc=1,n2m
      do 14 ic=1,n1m
      api(jc,ic)=-betadx
      aci(jc,ic)=1.+betadx*2.
      ami(jc,ic)=-betadx
      fi(jc,ic)=rhs(ic,jc,kc)
   14 continue
c
      call tripvi(1,n1m,1,n2m)
c
      do 6 jc=1,n2m
      do 6 ic=1,n1m
      rhs(ic,jc,kc)=fi(jc,ic)
    6 continue
    4 continue
c
c  ************ compute dq3 sweeping along the x2 direction
c
      aldt=dx2q*al*beta
      ucaj=4./(2.*cac(2)+caj(1))
      am=ucaj/cac(1)*2.
      ucaj=4./(2.*cac(n2m)+caj(n2m))
      ap=ucaj/cac(n2)*2.
      do 10 ic=1,n1m
      do 11 kc=1,n3m
      do 11 jc=1,n2m
      apj(kc,jc)=-aldt*apj3(jc)
      acj(kc,jc)=1.+aldt*acj3(jc)
      amj(kc,jc)=-aldt*amj3(jc)
   11 continue
      do 21 kc=1,n3m
      fj(kc,1)=rhs(ic,1,kc)+am*aldt*dq3s(ic,kc)*(1-islv3s)
      fj(kc,n2m)=rhs(ic,n2m,kc)+aldt*ap*dq3n(ic,kc)*(1-islv3n)
      do 22 jc=2,n2m
      fj(kc,jc)=rhs(ic,jc,kc)
   22 continue
   21 continue
c
      call trvbkj(n2m,1,n3m)
c
      do 30 kc=1,n3m
      do 30 jc=1,n2m
      dq(ic,jc,kc)=fj(kc,jc)
   30 continue
   10 continue
      return
      end
c****************************subrout prcalc ***************************
      subroutine prcalc(pre,pr,al)
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      dimension pre(m1,m2,m3),pr(m1,m2,m3)
      common/tstep/dt,beta,ren
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/metria/caj(m2),cac(m2)
c
c    the pressure is evaluated at the center of the box
c    at the near boundary cells newman b.c. for dph are assumed
c    pre=pressure pr=dph .
c
      be=al*beta
      do 1 kc=1,n3m
      kp=kpv(kc)
      km=kmv(kc)
      do 1 jc=1,n2m
      jm=jmv(jc)
      jp=jpv(jc)
      ac2=1./caj(jc)*dx2q
      co8=float(jp-jc)/cac(jp)
      co4=float(jc-jm)/cac(jc)
      do 1 ic=1,n1m
      ip=ipv(ic)
      im=imv(ic)
      pre(ic,jc,kc)=pre(ic,jc,kc)+pr(ic,jc,kc)-be*(
     1  dx1q*(pr(ip,jc,kc)-2.*pr(ic,jc,kc)+pr(im,jc,kc))+
     1  ac2*(co8*pr(ic,jp,kc)-(co8+co4)*pr(ic,jc,kc)+co4*pr(ic,jm,kc))+
     1  dx3q*(pr(ic,jc,kp)-2.*pr(ic,jc,kc)+pr(ic,jc,km))
     1                                            )
    1 continue
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine boucdq(time,tloc)
      include 'param.f'
      common/qwaloo/q1sold(m1,m3),q2sold(m1,m3),q3sold(m1,m3)
      common/qwalou/q1nold(m1,m3),q2nold(m1,m3),q3nold(m1,m3)
      common/qwallo/q1s(m1,m3),q2s(m1,m3),q3s(m1,m3)
      common/qwalup/q1n(m1,m3),q2n(m1,m3),q3n(m1,m3)
      common/dqwalo/dq1s(m1,m3),dq2s(m1,m3),dq3s(m1,m3)
      common/dqwalu/dq1n(m1,m3),dq2n(m1,m3),dq3n(m1,m3)
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/movwal/tosc,uosc
      common/slotfl/flowq2,tau2
      common/slotii/tim0sl
      common/vslot/v2inf(m1,m3)
      pi=2.*asin(1.)
      omeg1=2.*pi/tosc
      omeg2=2.*pi/tau2
      tpert=tim0sl+tau2
      dq1sma=0.
      dq2sma=0.
      dq3sma=0.
c
c   lower wall
c
      do kc=1,n3m
       do ic=1,n1m
       q1sold(ic,kc)=q1s(ic,kc)
       q2sold(ic,kc)=q2s(ic,kc)
       q3sold(ic,kc)=q3s(ic,kc)
       q1s(ic,kc)=0.
c      q2s(ic,kc)=v2inf(ic,kc)*sin(omeg2*tloc)
       q2s(ic,kc)=v2inf(ic,kc)*tloc
       q3s(ic,kc)=uosc*tloc
       dq1s(ic,kc)=q1s(ic,kc)-q1sold(ic,kc)
       dq2s(ic,kc)=q2s(ic,kc)-q2sold(ic,kc)
       dq3s(ic,kc)=q3s(ic,kc)-q3sold(ic,kc)
       dq1sma=max(abs(dq1s(ic,kc)),dq1sma)
       dq2sma=max(abs(dq2s(ic,kc)),dq2sma)
       dq3sma=max(abs(dq3s(ic,kc)),dq3sma)
       end do
      end do
c
c   upper wall at the moment set as no-slip wall
c   slip conditions for q1 and q3 are directly assigned in invtr1 and invtr3
c
      dq1nma=0.
      dq2nma=0.
      dq3nma=0.
      do kc=1,n3m
       do ic=1,n1m
       q1nold(ic,kc)=q1n(ic,kc)
       q2nold(ic,kc)=q2n(ic,kc)
       q3nold(ic,kc)=q3n(ic,kc)
       q1n(ic,kc)=0.
       q2n(ic,kc)=0.
       q3n(ic,kc)=uosc*tloc
       dq1n(ic,kc)=0.
       dq2n(ic,kc)=0.
       dq3n(ic,kc)=q3n(ic,kc)-q3nold(ic,kc)
       dq1nma=max(abs(dq1n(ic,kc)),dq1nma)
       dq2nma=max(abs(dq2n(ic,kc)),dq2nma)
       dq3nma=max(abs(dq3n(ic,kc)),dq3nma)
       end do
      end do
      if(uosc.gt.0.and.time.lt.tpert) then
      write(95,*)'t,dqisma  ',time,tloc,dq1sma,dq2sma,dq3sma
     1           ,dq1nma,dq2nma,dq3nma
c     print *,'q1 parete inferiore',tloc,sin(omega*tloc)
c     print *,'dq1 parete inferiore',dq1s(8,8)
      endif
      return
      end
