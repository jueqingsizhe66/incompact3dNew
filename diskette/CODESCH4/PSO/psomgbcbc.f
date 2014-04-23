c  **********************  set values at the old time step *********
c
      subroutine cveold(vor,psi,psc)
      include 'param.f'
      parameter(m1m=m1-1,m2m=m2-1)
      dimension vor(m1,m2),psi(m1,m2),psc(mpsc,m1,m2)
      common/dim/n1,n1m,n2,n2m
      common/mesh/dx1,dx1q,dx2,dx2q
      common/metri/g2m(m1),g2c(m1)
      common/firdve/dvox2n(m1),dpsx2n(m1),dvox2s(m1),dpsx2s(m1)
     1             ,dssx2n(mpsc,m1),dssx2s(mpsc,m1)
      common/qveold/voven(m1),psven(m1),voves(m1),psves(m1)
     1             ,ssven(mpsc,m1),ssves(mpsc,m1)
      common/pscqu/sch(mpsc),scla,npscf
c
c    quantities at the old time step
c    at the bound. on the upper wall
c    indicated by NORTH  in Fig.4.1 
c
      do i=1,n1
      dvox2n(i)=(vor(i,n2)-vor(i,n2-2))*dx2*0.5
      voven(i)=vor(i,n2m)
      dpsx2n(i)=(psi(i,n2)-psi(i,n2-2))*dx2*0.5
      psven(i)=psi(i,n2m)
      do npsc=1,npscf
      dssx2n(npsc,i)=(psc(npsc,i,n2)-psc(npsc,i,n2-2))*dx2*0.5
      ssven(npsc,i)=psc(npsc,i,n2m)
      enddo
      enddo
c
c    quantities at the old time step
c    at the bound. on the lower wall
c    indicated by SOUTH  in Fig.4.1 
c
      do i=1,n1
      dvox2s(i)=(vor(i,3)-vor(i,1))*dx2*0.5
      voves(i)=vor(i,2)
      dpsx2s(i)=(psi(i,3)-psi(i,1))*dx2*0.5
      psves(i)=psi(i,2)
      do npsc=1,npscf
      dssx2s(npsc,i)=(psc(npsc,i,3)-psc(npsc,i,1))*dx2*0.5
      ssves(npsc,i)=psc(npsc,i,2)
      enddo
      enddo
      return
      end
c
c  **********************  set values at the old time step *********
c
      subroutine choold(vor,psi,psc)
      include 'param.f'
      parameter(m1m=m1-1,m2m=m2-1)
      dimension vor(m1,m2),psi(m1,m2),psc(mpsc,m1,m2)
      common/dim/n1,n1m,n2,n2m
      common/mesh/dx1,dx1q,dx2,dx2q
      common/metri/g2m(m1),g2c(m1)
      common/firdho/dvox1e(m2),dpsx1e(m2),dvox1w(m2),dpsx1w(m2)
     1             ,dssx1e(mpsc,m2),dssx1w(mpsc,m2)
      common/qhoold/vohoe(m2),pshoe(m2),vohow(m2),pshow(m2)
     1             ,sshoe(mpsc,m2),sshow(mpsc,m2)
      common/pscqu/sch(mpsc),scla,npscf
c
c    quantities at the old time step
c    at the bound. on the right
c    indicated by EST  in Fig.4.1 
c
      udnn1=dx1*0.5/g2c(n1m)
      do j=1,n2
      dvox1e(j)=(vor(n1,j)-vor(n1-2,j))*udnn1
      vohoe(j)=vor(n1m,j)
      dpsx1e(j)=(psi(n1,j)-psi(n1-2,j))*udnn1
      pshoe(j)=psi(n1m,j)
      do npsc=1,npscf
      dssx1e(npsc,j)=(psc(npsc,n1,j)-psc(npsc,n1-2,j))*udnn1
      sshoe(npsc,j)=psc(npsc,n1m,j)
      enddo
      enddo
c
c    quantities at the old time step
c    at the bound. on the left
c    indicated by WEST  in Fig.4.1 
c
      udn1=dx1*0.5/g2c(2)
      do j=1,n2
      dvox1w(j)=(vor(3,j)-vor(1,j))*udn1
      vohow(j)=vor(2,j)
      dpsx1w(j)=(psi(3,j)-psi(1,j))*udn1
      pshow(j)=psi(2,j)
      do npsc=1,npscf
      dssx1w(npsc,j)=(psc(npsc,3,j)-psc(npsc,1,j))*udn1
      sshow(npsc,j)=psc(npsc,2,j)
      enddo
      enddo
      return
      end
c
c  ******* computes radiation velocity  *********
c  ******* at the NORTH SOUTH  boundaries
c
      subroutine cvevel(vor,psi,psc)
      include 'param.f'
      parameter(m1m=m1-1,m2m=m2-1)
      dimension vor(m1,m2),psi(m1,m2),psc(mpsc,m1,m2)
      common/dim/n1,n1m,n2,n2m
      common/coor/yp1(m1),yp2(m2)
      common/tstep/dt
      common/mesh/dx1,dx1q,dx2,dx2q
      common/metri/g2m(m1),g2c(m1)
      common/firdve/dvox2n(m1),dpsx2n(m1),dvox2s(m1),dpsx2s(m1)
     1             ,dssx2n(mpsc,m1),dssx2s(mpsc,m1)
      common/qveold/voven(m1),psven(m1),voves(m1),psves(m1)
     1             ,ssven(mpsc,m1),ssves(mpsc,m1)
      common/velcve/cven(2+mpsc,m1),cves(2+mpsc,m1)
      common/cwalv/inbcvs,inbcvn,inbcvw,inbcve
      common/cwals/inbcss,inbcsn,inbcsw,inbcse
      common/cwalp/inbcps,inbcpn,inbcpw,inbcpe
      common/pscqu/sch(mpsc),scla,npscf
c
c
c  calculation of b.c radiative velocity Orlansky JCP 
c
c  computational velocity compv
c
      compv=(yp2(n2)-yp2(n2m))/dt
      t=t+dt
c
c  calculation of radiative velocity for omega
c  at the boundary on the upper boundary north
 
c
      do i=1,n1
      dqx2n=(vor(i,n2)-vor(i,n2-2))*dx2*0.5
      dqt=(vor(i,n2m)-voven(i))/dt
      dqx2=(dqx2n+dvox2n(i))*0.5
       if(abs(dqx2).le..1e-05) dqx2=.1e-05
         covel1=-dqt/dqx2
         if(covel1.gt.0.)then
           if(covel1.gt.compv) covel1=compv*inbcvn 
           cven(1,i)=covel1*inbcvn 
c
c   this is an option to impose constant radiation velocity
c
c          cven(1,i)=0.6
                           else
           cven(1,i)=0.
                          endif
      enddo
c
c  calculation of radiative velocity for psi
c        upper boundary north
      do i=1,n1
      dqx2n=(psi(i,n2)-psi(i,n2-2))*dx2*0.5
      dqt=(psi(i,n2m)-psven(i))/dt
      dqx2=(dqx2n+dpsx2n(i))*0.5
       if(abs(dqx2).le..1e-05) dqx2=.1e-05
         covel2=-dqt/dqx2
         if(covel2.gt.0.)then
           if(covel2.gt.compv) covel2=compv*inbcpn 
           cven(2,i)=covel2*inbcpn
c
c   this is an option to impose constant radiation velocity
c
c          cven(2,i)=0.6
                           else
           cven(2,i)=0.
                          endif
      enddo
c
c  calculation of radiative velocity for passive scalar
c  at the boundary on the upper boundary north
 
c
      do npsc=1,npscf
      do i=1,n1
      dqx2n=(psc(npsc,i,n2)-psc(npsc,i,n2-2))*dx2*0.5
      dqt=(psc(npsc,i,n2m)-ssven(npsc,i))/dt
      dqx2=(dqx2n+dssx2n(npsc,i))*0.5
       if(abs(dqx2).le..1e-05) dqx2=.1e-05
         covel1=-dqt/dqx2
         if(covel1.gt.0.)then
           if(covel1.gt.compv) covel1=compv*inbcsn 
           cven(2+npsc,i)=covel1*inbcsn 
c
c   this is an option to impose constant radiation velocity
c
c          cven(2+npsc,i)=0.6
                           else
           cven(2+npsc,i)=0.
                          endif
      enddo
      enddo
c
c  calculation of radiative velocity for omega
c  at the boundary on the lower boundary south
c
c
      if(inbcvs.eq.1) then
      do i=1,n1
      dqx2s=(vor(i,3)-vor(i,1))*dx2*0.5
      dqt=(vor(i,2)-voves(i))/dt
      dqx2=(dqx2s+dvox2s(i))*0.5
       if(abs(dqx2).le..1e-05) dqx2=.1e-05
         covel1=-dqt/dqx2
         if(covel1.lt.0.)then
           if(abs(covel1).gt.compv) covel1=-compv*inbcvs
           cves(1,i)=covel1*inbcvs
                           else
           cves(1,i)=0.
                          endif
      enddo
                      endif
      if(inbcvs.eq.0.or.inbcvs.eq.2) then
      do i=1,n1
           cves(1,i)=0.
      enddo
                                     endif
c
c  calculation of radiative velocity for psi
c         lower boundary south
c
      if(inbcps.ge.0) then
      do i=1,n1
      dqx2s=(psi(i,3)-psi(i,1))*dx2*0.5
      dqt=(psi(i,2)-psves(i))/dt
      dqx2=(dqx2s+dpsx2s(i))*0.5
       if(abs(dqx2).le..1e-05) dqx2=.1e-05
         covel2=-dqt/dqx2
         if(covel2.lt.0.)then
           if(abs(covel2).gt.compv) covel2=-compv*inbcps
           cves(2,i)=covel2*inbcps
                           else
           cves(2,i)=0.
                          endif
      enddo
                      endif
c
c  calculation of radiative velocity for psc
c         lower boundary south
c
      do npsc=1,npscf
      if(inbcss.ge.0) then
      do i=1,n1
      dqx2s=(psc(npsc,i,3)-psc(npsc,i,1))*dx2*0.5
      dqt=(psc(npsc,i,2)-ssves(npsc,i))/dt
      dqx2=(dqx2s+dssx2s(npsc,i))*0.5
       if(abs(dqx2).le..1e-05) dqx2=.1e-05
         covel2=-dqt/dqx2
         if(covel2.lt.0.)then
           if(abs(covel2).gt.compv) covel2=-compv*inbcss
           cves(2+npsc,i)=covel2*inbcss
                           else
           cves(2+npsc,i)=0.
                          endif
      enddo
                      endif
      enddo
      return
      end
c
c  **********************  computes radiation velocity  *********
c  ******* at the EST WEST boundaries
c
      subroutine chovel(vor,psi,psc)
      include 'param.f'
      parameter(m1m=m1-1,m2m=m2-1)
      dimension vor(m1,m2),psi(m1,m2),psc(mpsc,m1,m2)
      common/dim/n1,n1m,n2,n2m
      common/coor/yp1(m1),yp2(m2)
      common/tstep/dt
      common/mesh/dx1,dx1q,dx2,dx2q
      common/metri/g2m(m1),g2c(m1)
      common/velcho/choe(2+mpsc,m2),chow(2+mpsc,m2)
      common/firdho/dvox1e(m2),dpsx1e(m2),dvox1w(m2),dpsx1w(m2)
     1             ,dssx1e(mpsc,m2),dssx1w(mpsc,m2)
      common/qhoold/vohoe(m2),pshoe(m2),vohow(m2),pshow(m2)
     1             ,sshoe(mpsc,m2),sshow(mpsc,m2)
      common/cwalv/inbcvs,inbcvn,inbcvw,inbcve
      common/cwalp/inbcps,inbcpn,inbcpw,inbcpe
      common/cwals/inbcss,inbcsn,inbcsw,inbcse
      common/pscqu/sch(mpsc),scla,npscf
c
c
c  calculation of b.c radiative velocity Orlansky JFM pp.254
c
c  computational velocity compv
c
      compv=(yp1(n1)-yp1(n1m))/dt
      t=t+dt
c
c  calculation of radiative velocity for omega
c  at the boundary on the right
 
c
      udnn1=dx1*0.5/g2c(n1m)
      udn1=dx1*0.5/g2c(2)
      if(inbcve.ge.0) then
      do j=1,n2
      dq1x1e=(vor(n1,j)-vor(n1-2,j))*udnn1
      dq1t=(vor(n1m,j)-vohoe(j))/dt
      dq1x1=(dq1x1e+dvox1e(j))*0.5
       if(abs(dq1x1).le..1e-05) dq1x1=.1e-05
         covel1=-dq1t/dq1x1
         if(covel1.gt.0.)then
           if(covel1.gt.compv) covel1=compv*inbcve
           choe(1,j)=covel1*inbcve
                           else
           choe(1,j)=0.
                          endif
      enddo
                       endif
c
c  calculation of radiative velocity for psi
c
      if(inbcpe.ge.0) then
      do j=1,n2
      dq2x1e=(psi(n1,j)-psi(n1-2,j))*udnn1
      dq2t=(psi(n1m,j)-pshoe(j))/dt
      dq2x1=(dq2x1e+dpsx1e(j))*0.5
       if(abs(dq2x1).le..1e-05) dq2x1=.1e-05
         covel2=-dq2t/dq2x1
         if(covel2.gt.0.)then
           if(covel2.gt.compv) covel2=compv*inbcpe
           choe(2,j)=covel2*inbcpe
                           else
           choe(2,j)=0.
                           endif       
      enddo
                       endif
c
c  calculation of radiative velocity for psc
c
      if(inbcse.gt.0) then
      do npsc=1,npscf
      do j=1,n2
      dq2x1e=(psc(npsc,n1,j)-psc(npsc,n1-2,j))*udnn1
      dq2t=(psc(npsc,n1m,j)-sshoe(npsc,j))/dt
      dq2x1=(dq2x1e+dssx1e(npsc,j))*0.5
       if(abs(dq2x1).le..1e-05) dq2x1=.1e-05
         covel2=-dq2t/dq2x1
         if(covel2.gt.0.)then
           if(covel2.gt.compv) covel2=compv*inbcse
           choe(2+npsc,j)=covel2*inbcse
                           else
           choe(2+npsc,j)=0.
                           endif       
      enddo
      enddo
                      endif
      if(inbcse.eq.0) then
      do npsc=1,npscf
      do j=1,n2
           choe(2+npsc,j)=0.
      enddo
      enddo
                      endif
c
c  calculation of radiative velocity for omega
c  at the boundary on the left
c
 
c
      if(inbcvw.ge.0) then
      do j=1,n2
      dq1x1e=(vor(3,j)-vor(1,j))*udn1
      dq1t=(vor(2,j)-vohow(j))/dt
      dq1x1=(dq1x1e+dvox1w(j))*0.5
       if(abs(dq1x1).le..1e-05) dq1x1=.1e-05
         covel1=-dq1t/dq1x1
         if(covel1.lt.0.)then
           if(abs(covel1).gt.compv) covel1=-compv*inbcvw
           chow(1,j)=covel1*inbcvw
                           else
           chow(1,j)=0.
                          endif
      enddo
                     endif
c
c  calculation of radiative velocity for psi
c
      if(inbcpw.ge.0) then
      do j=1,n2
      dq2x1e=(psi(3,j)-psi(1,j))*udn1
      dq2t=(psi(2,j)-pshow(j))/dt
      dq2x1=(dq2x1e+dpsx1w(j))*0.5
       if(abs(dq2x1).le..1e-05) dq2x1=.1e-05
         covel2=-dq2t/dq2x1
         if(covel2.lt.0.)then
           if(abs(covel2).gt.compv) covel2=-compv*inbcpw
           chow(2,j)=covel2*inbcpw
                           else
           chow(2,j)=0.
                          endif
      enddo
                     endif
c
c  calculation of radiative velocity for psc
c
      if(inbcsw.ge.0) then
      do npsc=1,npscf
      do j=1,n2
      dq2x1e=(psc(npsc,3,j)-psc(npsc,1,j))*udn1
      dq2t=(psc(npsc,2,j)-sshow(npsc,j))/dt
      dq2x1=(dq2x1e+dssx1w(npsc,j))*0.5
       if(abs(dq2x1).le..1e-05) dq2x1=.1e-05
         covel2=-dq2t/dq2x1
         if(covel2.lt.0.)then
           if(abs(covel2).gt.compv) covel2=-compv*inbcsw
           chow(2+npsc,j)=covel2*inbcsw
                           else
           chow(2+npsc,j)=0.
                          endif
      enddo
      enddo
                     endif
      return
      end
c
c  ****************************** subrout dsepns  **********************
c
c  this subroutine calculate the second derivative of psi  
c to evaluate the boundary conditions described in Sect.4.4
c            NORTH SOUTH Boundaries
c
      subroutine dsepns
      include 'param.f'
      common/dim/n1,n1m,n2,n2m
      common/bchodp/dpsud(m1),dpnor(m1)
      common/dsebns/d2psud(m1),d2pnor(m1)
      common/mesh/dx1,dx1q,dx2,dx2q
      common/metri/g2m(m1),g2c(m1)
c
c    second derivatives of dpsi at the upper and lower bound.
c
      do i=2,n1m
        d2psud(i)=((dpsud(i+1)-dpsud(i))/g2m(i)
     1            -(dpsud(i)-dpsud(i-1))/g2m(i-1))*dx1q/g2c(i)
        d2pnor(i)=((dpnor(i+1)-dpnor(i))/g2m(i)
     1            -(dpnor(i)-dpnor(i-1))/g2m(i-1))*dx1q/g2c(i)
      enddo
      return
      end
c
c  ****************************** subrout bdqpsn  **********************
c
c  this subroutine calculate the b.c for vor , psi   and psc
c            NORTH SOUTH Boundaries
c
      subroutine bdqpsn(gam,ro,vor,psi,psc,dft)
      include 'param.f'
      dimension psi(m1,m2),vor(m1,m2),v1(2),psc(mpsc,m1,m2)
      common/dim/n1,n1m,n2,n2m
      common/bchov/vsud(m1),vnor(m1)
      common/bchop/psud(m1),pnor(m1)
      common/bchodv/dvsud(m1),dvnor(m1)
      common/bchodp/dpsud(m1),dpnor(m1)
      common/bchods/dssud(mpsc,m1),dsnor(mpsc,m1)
      common/infmxl/voins(m1),psins(m1),ssins(mpsc,m1)
      common/infmxp/voinp(m1),psifr(m1),
     %              psifi(m1),psiar(m1),psiai(m1),psinp(m1)
      common/mesh/dx1,dx1q,dx2,dx2q
      common/metri/g2m(m1),g2c(m1)
      common/tstep/dt
      common/velcve/cven(2+mpsc,m1),cves(2+mpsc,m1)
      common/dvpsvl/dvnoro(m1),dpnoro(m1),dvsudo(m1),dpsudo(m1)
     1         ,dsnoro(mpsc,m1),dssudo(mpsc,m1)
      common/cwalv/inbcvs,inbcvn,inbcvw,inbcve
      common/cwalp/inbcps,inbcpn,inbcpw,inbcpe
      common/cwals/inbcss,inbcsn,inbcsw,inbcse
      common/pscqu/sch(mpsc),scla,npscf
      common/dinfmxp/dvoinp(m1),dpsinp(m1) 
      common/bcvedv/dvest(m2),dvwes(m2)
      common/bcvedp/dpest(m2),dpwes(m2)
c
c    radiative b.c. for d(vorticity)
c    or free-slip  north
c
      if(inbcvn.eq.1.or.inbcvn.eq.0) then
      do i=1,n1
        dvox2n= (vor(i,n2)-vor(i,n2m))*dx2
        dvnor(i)= -dt*(gam*dvox2n+ro*dvnoro(i))*cven(1,i)
        dvnoro(i)=dvox2n
      enddo
                                  else
      write(6,*)'program stops no right b.c. at the north vor'
      stop
                                  endif
      if(inbcpn.eq.1) then
c
c    radiative b.c. for d(stream function)
c    free-slip or radiative   north
c
C     do i=1,n1
C       dpnor(i) = psins(i)*dft
C     enddo
      do i=1,n1
        dpsx2n = (psi(i,n2)-psi(i,n2m))*dx2
        dpnor(i)= -dt*(gam*dpsx2n+ro*dpnoro(i))*cven(2,i)
        dpnoro(i)=dpsx2n
      enddo
                                  endif
      if(inbcpn.eq.0) then
c
      do i=1,n1
        dpnor(i)= dpwes(n2)
      enddo
                                  endif
      if(inbcsn.eq.1.or.inbcsn.eq.0) then
c
c    radiative b.c. for d(passive scalar )
c    free-slip or radiative   north
c
      do npsc=1,npscf
      do i=1,n1
        dssx2n = (psc(npsc,i,n2)-psc(npsc,i,n2m))*dx2
        dsnor(npsc,i)= -dt*(gam*dssx2n+ro*dsnoro(npsc,i))*cven(2+npsc,i)
        dsnoro(npsc,i)=dssx2n
      enddo
      enddo
                                  else
      write(6,*)'program stops no right b.c. at the north psc'
      stop
                                  endif
      if(inbcvs.ge.1.or.inbcvs.eq.0) then
c
c    radiative b.c. for d(vorticity)
c    or free-slip  south
c
      do i=1,n1
        dvox2s= (vor(i,2)-vor(i,1))*dx2
        dvsud(i)= -dt*(gam*dvox2s+ro*dvsudo(i))*cves(1,i)
        dvsudo(i)=dvox2s
      enddo
                                     endif
      if(inbcvs.eq.-1) then
      ugmm2=4./3.
      do ic=2,n1m
      ip=ic+1
      im=ic-1
      do jc=1,2
      v1(jc)=(psi(ic,jc+1)-psi(ic,jc))*dx2
      enddo
      jb=1
      jc=jb
      js=+1
      jp=jc+1
      dq12p=(v1(jp)-v1(jc))
      dq12mw=+2.*v1(jc)
      dq12=(dq12p-dq12mw)*dx2*ugmm2
      dq21=-((psi(ip,jb+js)-psi(ic,jb+js))/g2m(ic)
     1      -(psi(ic,jb+js)-psi(im,jb+js))/g2m(ic-1))*dx1q/g2c(ic)
      vorjb=-js*(-dq12+js*dq21)
      dvsud(ic)=vorjb-dvsudo(ic)
      dvsudo(ic)=vorjb
      enddo
      do ic=1,n1,n1m
      do jc=1,2
      v1(jc)=(psi(ic,jc+1)-psi(ic,jc))*dx2
      enddo
      jb=1
      jc=jb
      js=+1
      jp=jc+1
      dq12p=(v1(jp)-v1(jc))
      dq12mw=+2.*v1(jc)
      dq12=(dq12p-dq12mw)*dx2*ugmm2
      vorjb=-js*(-dq12)
      dvsud(ic)=vorjb-dvsudo(ic)
      dvsudo(ic)=vorjb
      enddo
                      endif
      if(inbcvs.eq.-2) then
c
c    inlet b.c. for d(vorticity) at sud boundary
c
      do i=1,n1
C       dvsud(i) = (voins(i)+voinp(i))*dft
c       dvsud(i) = voins(i)*dft + dvoinp(i)
        dvsud(i) = voins(i)*dft
      enddo
                      endif
      if(inbcps.eq.1) then
      do i=1,n1
        dpsx2s = (psi(i,2)-psi(i,1))*dx2
        dpsud(i)= -dt*(gam*dpsx2s+ro*dpsudo(i))*cves(2,i)
        dpsudo(i)=dpsx2s
      enddo
                                     endif
      if(inbcps.eq.0) then
      do i=1,n1
        dpsud(i)= dpwes(1) 
      enddo
                                     endif
c
c    no-slip b.c. for d(streamfunct) at sud boundary
c
      if(inbcps.eq.-1) then
      do i=1,n1
        dpsud(i)= 0.
      enddo
                      endif
      if(inbcps.eq.-2) then
c
c    inlet b.c. for d(streamfunct) at sud boundary
c
      do i=1,n1
c       dpsud(i) = psins(i)*dft + dpsinp(i)
C       dpsud(i) =( psins(i)+psinp(i))*dft
        dpsud(i) = psins(i)*dft
      enddo
                      endif
      if(inbcss.eq.1.or.inbcss.eq.0) then
      do npsc=1,npscf
      do i=1,n1
        dssx2s = (psc(npsc,i,2)-psc(npsc,i,1))*dx2
        dssud(npsc,i)= -dt*(gam*dssx2s+ro*dssudo(npsc,i))*cves(2+npsc,i)
        dssudo(npsc,i)=dssx2s
      enddo
      enddo
                                     endif
      if(inbcss.eq.-1) then
      do npsc=1,npscf
      do i=1,n1
        dssud(npsc,i)= 0.
      enddo
      enddo
                      endif
      if(inbcss.eq.-2) then
      do npsc=1,npscf
      do i=1,n1
c       dssud(npsc,i)= 0. 
       dssud(npsc,i)= ssins(npsc,i)*dft
      enddo
      enddo
                      endif
      return
      end
c
c  **********************  set boundary conditions  *********
c
c  this subroutine calculate the b.c for vor , psi   and psc
c            EST WEST Boundaries
c
c
      subroutine bdqpew(gam,ro,vor,psi,psc,dft)
      include 'param.f'
      parameter(m1m=m1-1,m2m=m2-1)
      common/pscqu/sch(mpsc),scla,npscf
      dimension vor(m1,m2),psi(m1,m2),psc(mpsc,m1,m2),v2(2)
      common/dim/n1,n1m,n2,n2m
      common/tstep/dt
      common/mesh/dx1,dx1q,dx2,dx2q
      common/metri/g2m(m1),g2c(m1)
      common/bcvedv/dvest(m2),dvwes(m2)
      common/bcvedp/dpest(m2),dpwes(m2)
      common/bchodp/dpsud(m1),dpnor(m1)
      common/bcveds/dsest(mpsc,m2),dswes(mpsc,m2)
      common/bcvev/vwes(m2),vest(m2)
      common/bcvep/pwes(m2),pest(m2)
      common/velcho/choe(2+mpsc,m2),chow(2+mpsc,m2)
      common/dvpsol/dvesto(m2),dpesto(m2),dvweso(m2),dpweso(m2)
     1         ,dsesto(mpsc,m2),dsweso(mpsc,m2)
      common/inflow/voinf(m2),psinf(m2),ssinf(mpsc,m2)
      common/cwalv/inbcvs,inbcvn,inbcvw,inbcve
      common/cwals/inbcss,inbcsn,inbcsw,inbcse
      common/cwalp/inbcps,inbcpn,inbcpw,inbcpe
c
c    radiative b.c. for d(vorticity) at left boundary
c    for inbcvw=0 gives symmetry 0 vorticity
c
      udn1=dx1/g2m(1)
      udnn1=dx1/g2m(n1m)
      if(inbcvw.eq.1.or.inbcvw.eq.0) then
      do j=1,n2
        dvoix1 = (vor(2,j)-vor(1,j))*udn1
        dvwes(j) = -dt*(gam*dvoix1+ro*dvweso(j))*chow(1,j)
        dvweso(j)=dvoix1
      enddo
                                     endif
c
c     no-slip boundary for vorticity
c
      if(inbcvw.eq.-1) then
      ugmm2=4./(g2m(2)+2.*g2c(1))
      do jc=2,n2m
      jp=jc+1
      jm=jc-1
      v21=-(psi(2,jc)-psi(1,jc))*dx1/g2m(1)
      v22=-(psi(3,jc)-psi(2,jc))*dx1/g2m(2)
      dq21m=(v22-v21)*dx1/g2c(2)
      dq21pw=2.*v21*dx1/g2c(1)
      dq21=-(dq21pw-dq21m)*ugmm2*g2m(1)
      dq12=-((psi(2,jp)-psi(2,jc))
     1      -(psi(2,jc)-psi(2,jm)))*dx2q
      vorib=(dq21-dq12)
      dvwes(jc)=vorib-dvesto(jc)
      dvweso(jc)=vorib
      enddo
                       endif
      if(inbcvw.eq.-2) then
c
c    inlet b.c. for d(vorticity) at left boundary
c    
      do j=1,n2
        dvwes(j) = voinf(j)*dft
      enddo
                      endif
c
c    radiative b.c. for d(stream function) at left boundary
c    for inbcpw=0 gives symmetry 0 streamfunction
c
C     if(inbcpw.eq.1.or.inbcpw.eq.0) then
C     do j=1,n2
C       dpsix1 = (psi(2,j)-psi(1,j))*udn1
C       dpwes(j) = -dt*(gam*dpsix1+ro*dpweso(j))*chow(2,j)
C       dpweso(j)=dpsix1
C     enddo
C                                    endif
      if(inbcpw.eq.1) then
        do j=1,n2
          dpsix1 = (psi(2,j)-psi(1,j))*udn1
          dpwes(j) = -dt*(gam*dpsix1+ro*dpweso(j))*chow(2,j)
          dpweso(j)=dpsix1
        enddo
      endif
      if(inbcpw.eq.0) then
        do j=1,n2
          dpwes(j) = dpsud(1)
        enddo
      endif
c
c     no-slip boundary for streamfunct
c
      if(inbcpw.eq.-1) then
      do j=1,n2
        dpwes(j) = 0.
      enddo
                      endif
      if(inbcpw.eq.-2) then
c
c    inlet b.c. for d(streamfunc) at left boundary
c    
      do j=1,n2
        dpwes(j) = psinf(j)*dft
      enddo
                      endif
c
c    radiative b.c. for d(passive scalar ) at left boundary
c
      if(inbcsw.eq.1.or.inbcsw.eq.0) then
      do npsc=1,npscf
      do j=1,n2
       dssix1 = (psc(npsc,2,j)-psc(npsc,1,j))*udn1
       dswes(npsc,j) = -dt*(gam*dssix1+ro*dsweso(npsc,j))*chow(2+npsc,j)
       dsweso(npsc,j)=dssix1
      enddo
      enddo
                                      endif
      if(inbcsw.eq.-1) then
c
c    b.c. for d(passive scala) at left boundary zero flux no-slip
c    
      do npsc=1,npscf
      do j=1,n2
        dswes(npsc,j) = 0.
      enddo
      enddo
                      endif
      if(inbcsw.eq.-2) then
c
c    inlet b.c. for d(passive scala) at left boundary
c    
      do npsc=1,npscf
      do j=1,n2
        dswes(npsc,j) = ssinf(npsc,j)*dft
      enddo
      enddo
                      endif
      if(inbcve.eq.1) then
c
c    radiative b.c. for d(vorticity) at rigth boundary
c    for inbcve=0 gives symmetry 0 vorticity
c
      do j=1,n2
        dvoox1 = (vor(n1,j)-vor(n1m,j))*udnn1
        dvest(j) = -dt*(gam*dvoox1+ro*dvesto(j))*choe(1,j)
        dvesto(j)=dvoox1
      enddo
                                     endif
      if(inbcve.eq.2.or.inbcve.eq.0) then
c
c    radiative b.c. for d(vorticity) at rigth boundary
c    for inbcve=0 gives symmetry 0 vorticity
c
      do j=1,n2
        dvest(j) = 0.
      enddo
                                     endif
      if(inbcve.eq.-1) then
      ugmm2=4./(g2m(n1m)+2.*g2c(n1))
      do jc=2,n2m
      jp=jc+1
      jm=jc-1
      v2n1m=-(psi(n1,jc)-psi(n1m,jc))*dx1/g2m(n1m)
      v2n1d=-(psi(n1m,jc)-psi(n1m-1,jc))*dx1/g2m(n1m-1)
      dq21m=(v2n1m-v2n1d)*dx1/g2c(n1m)
      dq21pw=-2.*v2n1m*dx1/g2c(n1)
      dq21=(dq21pw-dq21m)*ugmm2*g2m(n1m)
      dq12=-((psi(n1m,jp)-psi(n1m,jc))
     1      -(psi(n1m,jc)-psi(n1m,jm)))*dx2q
      vorib=(dq21-dq12)
      dvest(jc)=vorib-dvesto(jc)
      dvesto(jc)=vorib
      enddo
                                     endif
c
c    radiative b.c. for d(stream function) at left boundary
c    for inbcpw=0 gives symmetry 0 streamfunction
c
C     if(inbcpe.eq.1.or.inbcpe.eq.0) then
C     do j=1,n2
C       dpsox1 = (psi(n1,j)-psi(n1m,j))*udnn1
C       dpest(j) = -dt*(gam*dpsox1+ro*dpesto(j))*choe(2,j)
C       dpesto(j)=dpsox1
C     enddo
C                                    endif
c
c    radiative b.c. for d(stream function) at left boundary
c    for inbcpw=0 gives symmetry 0 streamfunction
c
      if(inbcpe.eq.1) then
        do j=1,n2
          dpsox1 = (psi(n1,j)-psi(n1m,j))*udnn1
          dpest(j) = -dt*(gam*dpsox1+ro*dpesto(j))*choe(2,j)
          dpesto(j)=dpsox1
        enddo
      endif
      if(inbcpe.eq.0) then
        do j=1,n2
          dpest(j) = dpsud(n1)
        enddo
      endif
c
c     no-slip boundary for streamfunct
c
      if(inbcpe.eq.-1) then
      do j=1,n2
        dpest(j) = 0.
      enddo
                       endif
      if(inbcse.eq.1.or.inbcse.eq.0) then
      do npsc=1,npscf
      do j=1,n2
       dssox1 = (psc(npsc,n1,j)-psc(npsc,n1m,j))*udnn1
       dsest(npsc,j) = -dt*(gam*dssox1+ro*dsesto(npsc,j))*choe(2+npsc,j)
       dsesto(npsc,j)=dssox1
      enddo
      enddo
                                     endif
      if(inbcse.eq.-1) then
      do npsc=1,npscf
      do j=1,n2
        dsest(npsc,j) = 0.
      enddo
      enddo
                       endif
      return
      end
c
c  **********************  set boundary conditions  *********
c
c   evaluates the boundary conditions for phi at the
c      EST WEST Boundaries as described in Sect.4.4
c
      subroutine bdppew(psi)
      include 'param.f'
      parameter(m1m=m1-1,m2m=m2-1)
      common/dim/n1,n1m,n2,n2m
      common/bcvedp/dpest(m2),dpwes(m2)
      common/bvedph/dphest(m2),dphwes(m2)
      common/ydif/ydi(m2)
      dimension psi(m1,m2)

c
c    b.c. for dphi from physical b.c. for dpsi
c
      do 101 j=1,n2m
        dphest(j)=dpest(j)+psi(n1,j)
        dphwes(j)=dpwes(j)+psi(1,j)
  101 continue
      return
      end
c
c  **********************  set boundary conditions  *********
c
c   evaluates the boundary conditions for phi at the
c      EST WEST Boundaries as described in Sect.4.4
c
      subroutine bdphew
      include 'param.f'
      parameter(m1m=m1-1,m2m=m2-1)
      common/dim/n1,n1m,n2,n2m
      common/bcvedp/dpest(m2),dpwes(m2)
      common/bvedph/dphest(m2),dphwes(m2)
      common/bchodp/dpsud(m1),dpnor(m1)
      common/ydif/ydi(m2)

c
c    b.c. for dphi from physical b.c. for dpsi
c
      do 101 j=1,n2
        dphest(j)=dpest(j)
     1           -(ydi(j)*(dpnor(n1)-dpsud(n1))+dpsud(n1))
        dphwes(j)=dpwes(j)
     1           -(ydi(j)*(dpnor(1)-dpsud(1))+dpsud(1))
  101 continue
      return
      end
c
c  **********************  updates omega psi and psc (boundary) ********
c
      subroutine boucqt
      include 'param.f'
      parameter(m1m=m1-1,m2m=m2-1)
      common/dim/n1,n1m,n2,n2m
      common/bcvedv/dvest(m2),dvwes(m2)
      common/bcvedp/dpest(m2),dpwes(m2)
      common/bcveds/dsest(mpsc,m2),dswes(mpsc,m2)
      common/bcvev/vwes(m2),vest(m2)
      common/bcvep/pwes(m2),pest(m2)
      common/bcves/swes(mpsc,m2),sest(mpsc,m2)
      common/bchov/vsud(m1),vnor(m1)
      common/bchop/psud(m1),pnor(m1)
      common/bchos/ssud(mpsc,m1),snor(mpsc,m1)
      common/bchodv/dvsud(m1),dvnor(m1)
      common/bchodp/dpsud(m1),dpnor(m1)
      common/bchods/dssud(mpsc,m1),dsnor(mpsc,m1)
      common/pscqu/sch(mpsc),scla,npscf
c
      do j=1,n2
         vwes(j) = vwes(j)+dvwes(j)
         vest(j) = vest(j)+dvest(j)
         pwes(j) = pwes(j)+dpwes(j)
         pest(j) = pest(j)+dpest(j)
         do npsc=1,npscf
         swes(npsc,j) = swes(npsc,j)+dswes(npsc,j)
         sest(npsc,j) = sest(npsc,j)+dsest(npsc,j)
         enddo
      enddo
      do i=1,n1
         vnor(i) = vnor(i)+dvnor(i)
         vsud(i) = vsud(i)+dvsud(i)
         pnor(i) = pnor(i)+dpnor(i)
         psud(i) = psud(i)+dpsud(i)
         do npsc=1,npscf
         snor(npsc,i) = snor(npsc,i)+dsnor(npsc,i)
         ssud(npsc,i) = ssud(npsc,i)+dssud(npsc,i)
         enddo
      enddo
      return
      end
c
c  **********************  set b.c. for spa. dev. m.l.  *********
c
C 
      subroutine bcspaml(time,timeo)
      include 'param.f'
      parameter(m1m=m1-1,m2m=m2-1)
      common/dim/n1,n1m,n2,n2m
      common/inflow/voinf(m2),psinf(m2),ssinf(mpsc,m2)
      common/dinfmxp/dvoinp(m1),dpsinp(m1) 
      common/infmxp/voinp(m1),psifr(m1),
     %              psifi(m1),psiar(m1),psiai(m1),psinp(m1)
      common/infvop/voinfr(m1),voinfi(m1),
     %              voinar(m1),voinai(m1)
C
      omeg1 = 0.65
      omeg2 = 0.325
      dft = 1.
      ampl1 = 0.005
      ampl2 = 0.000
      dvom = -10.
      dvm = -10.
      dpm = -10.
c
c    inlet b.c. for d(vorticity) at south boundary
c    
      do i=1,n1
       dvoinp(i) = ampl1*(voinfr(i)*cos(omeg1*time)+
     %                    voinfi(i)*sin(omeg1*time))
     %           + ampl2*(voinar(i)*cos(omeg2*time)+
     %                    voinai(i)*sin(omeg2*time))
     %           - ampl1*(voinfr(i)*cos(omeg1*timeo)+
     %                    voinfi(i)*sin(omeg1*timeo))
     %           - ampl2*(voinar(i)*cos(omeg2*timeo)+
     %                    voinai(i)*sin(omeg2*timeo))
      dvom = max(dvom,abs(voinfr(i)))
      dvm = max(dvm,abs(dvoinp(i)))
      enddo
c
c    inlet b.c. for d(streamfunc) at south boundary
c    
      do i=1,n1
       dpsinp(i) = ampl1*(psifr(i)*cos(omeg1*time)+
     %                    psifi(i)*sin(omeg1*time))
     %           + ampl2*(psiar(i)*cos(omeg2*time)+
     %                    psiai(i)*sin(omeg2*time))
     %           - ampl1*(psifr(i)*cos(omeg1*timeo)+
     %                    psifi(i)*sin(omeg1*timeo))
     %           - ampl2*(psiar(i)*cos(omeg2*timeo)+
     %                    psiai(i)*sin(omeg2*timeo))
      dpm = max(dpm,abs(dpsinp(i)))
c       write(442,*) psifr(i),psifi(i),
c    %              psiar(i),psiai(i)
      enddo
c     close(442)
c     write(6,*) ' DVM = ',dvm,' DPM = ',dpm,' DVOM = ',dvom
c
      return
      end
