c************************************************************************
c
c     this subroutine calculates the solenoidal vel field.
c       q(n+1)=qhat-grad(dph)*dt ,  pr=dph
c    third order runge-kutta is used.
c
c************************************************************************
      subroutine updvp(al)
      include'param.f'
c
c  ***********  compute the q1 velocity component
c               v1dgf=component 1 of grad(dph)
c   q1=v_theta
      do 1 kc=1,n3m
        do 1 jc=1,n2m
          usurm = al*dt*dx1/rm(jc)
          do 1 ic=1,n1m
            im=imv(ic)
            q1(ic,jc,kc)=q1(ic,jc,kc)-
     %      (dph(ic,jc,kc)-dph(im,jc,kc))*usurm
    1 continue
c
c
c  ***********  compute the q2 velocity component
c               v2dgf=component 2 of grad(dph)
      do 2 kc=1,n3m
        do 2 jc=2,n2m
          jm=jc-1
          usurm = al*dt*dx2*rc(jc)/g2rc(jc)
          do 2 ic=1,n1m
            q2(ic,jc,kc)=q2(ic,jc,kc)-
     %      (dph(ic,jc,kc)-dph(ic,jm,kc))*usurm
    2 continue
      do 21 kc=1,n3m
        do 21 ic=1,n1m
          q2(ic,1,kc)=0.
          q2(ic,n2,kc)=0.
   21 continue
c
c  ***********  compute the q3 velocity component
c               q3 is the cartesian component
c               v3dgf=component 3 of grad(dph)
      do 5 kc=2,n3m
      usurm = al*dt*dx3/g3rc(kc)
        km=kc-1
        do 5 jc=1,n2m
          do 5 ic=1,n1m
            q3(ic,jc,kc)=q3(ic,jc,kc)-
     %      (dph(ic,jc,kc)-dph(ic,jc,km))*usurm
    5 continue
      do 210 jc=1,n2m
        do 210 ic=1,n1m
          q3(ic,jc,1)=0.
          q3(ic,jc,n3)=0.
  210 continue
      return
      end
************************************************************************
c   this subroutine performs the calculation of the pressure.
c
************************************************************************
      subroutine prcalc(al)
      include'param.f'
c
c    the pressure is evaluated at the center of the box.
c
c     p^{n+1} = p^{n} + phi^{n+1} - b * Nabla^2 phi^{n+1}
c
      be=al*beta
        do 1 jc=1,n2m
          jm=jmv(jc)
          jp=jpv(jc)
          cthe=dx1q/rm(jc)**2
      do 1 kc=1,n3m
        kp=kpv(kc)
        km=kmv(kc)
            do 1 ic=1,n1m
              ip=ipv(ic)
              im=imv(ic)
              pr(ic,jc,kc)=pr(ic,jc,kc)+dph(ic,jc,kc)-be*(
     %      (dph(ip,jc,kc)-2.*dph(ic,jc,kc)+dph(im,jc,kc))*cthe + 
     %      (dph(ic,jc,kp)*apphk(kc)
     %      +dph(ic,jc,kc)*acphk(kc)
     %      +dph(ic,jc,km)*amphk(kc)) + 
     %      (dph(ic,jp,kc)*apphj(jc)
     %      +dph(ic,jc,kc)*acphj(1,jc)
     %      +dph(ic,jm,kc)*amphj(jc))                    )
    1 continue
      return
      end
c************************************************************************
c
c   this subroutine performs the inversion of the q1 momentum equation
c   by a factored implicit scheme, only the derivatives 11,22,33 of q1
c   are treated implicitly
c       direction x3
c************************************************************************
      subroutine solq1k(al)
      include'param.f'
c
c   compute  from dq** sweeping along the x3 periodic direction
c   
      n3mm=n3m-1
      ugkk=4./3.
      betadx=beta*al
      do 4 kc=1,n3m
        amk(kc)=-am3sk(kc)*betadx
        ack(kc)=1.-ac3sk(kc)*betadx
        apk(kc)=-ap3sk(kc)*betadx
    4 continue
c
      call brtrk(amk,ack,apk,n1m,n2m,n3m)
c
      do kc=1,n3m
       do jc=1,n2m
        do ic=1,n1m
         iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2m
         q1(ic,jc,kc)=q1(ic,jc,kc) + rhs(iadd)
         rhs(iadd) = 0.
        end do
       end do
      end do
c
c
      return
      end
c************************************************************************
c
c   this subroutine performs the inversion of the q1 momentum equation
c   by a factored implicit scheme, only the derivatives 11,22,33 of q1
c   are treated implicitly
c         direction x1
c
c************************************************************************
      subroutine solq1i(al)
      include'param.f'
c
c  ********* compute dq1*  sweeping in the x1 direction
c            
      betadx=beta*dx1q*al
        do jc=1,n2m
            ugmmm=betadx/(rm(jc)**2)
            api(jc)=-ugmmm
            ami(jc)=-ugmmm
            aci(jc)=1.+ugmmm*2.
         end do
c
        call tripvmy(1,n1m,1,n2m,1,n3m,n1m,n2m)
c
      return
      end
c************************************************************************
c
c   this subroutine performs the inversion of the q1 momentum equation
c   by a factored implicit scheme, only the derivatives 11,22,33 of q1
c   are treated implicitly
c
c************************************************************************
      subroutine solq1j(al)
      include'param.f'
c
c  ************ compute dq1  sweeping along the x2 direction
c
      betadx=beta*al
      do 21 jc=1,n2m
      apj(jc)= -ap1j(jc)*betadx
      acj(jc)=1.-ac1j(jc)*betadx
      amj(jc)= -am1j(jc)*betadx
   21 continue
c
      call brtrj(amj,acj,apj,n1m,n2m,n3m)
c
      return
      end
c************************************************************************
c   this subroutine performs the inversion of the q2 momentum equation
c   by a factored implicit scheme, only the derivatives 11,22,33 of q2
c   are treated implicitly
c   in the first part the rhs is calculated
c        direction x3
c
c************************************************************************
      subroutine solq2k(al)
      include'param.f'
c
c  ************ compute dq2 sweeping along the x3 direction
c               
      n3mm=n3m-1
      ugkk=4./3.
      betadx=beta*al
      do 6 kc=1,n3m
        amk(kc)=-am3sk(kc)*betadx
        ack(kc)=1.-ac3sk(kc)*betadx
        apk(kc)=-ap3sk(kc)*betadx
    6 continue
c
      call brtrk(amk,ack,apk,n1m,n2,n3m)
c
      do kc=1,n3m
       do jc=2,n2m
        do ic=1,n1m
         iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2
         q2(ic,jc,kc)=q2(ic,jc,kc) + rhs(iadd)
         rhs(iadd) = 0.
        end do
       end do
      end do
c
      return
      end
c ************************************************************************
c   this subroutine performs the inversion of the q2 momentum equation
c   by a factored implicit scheme, only the derivatives 11,22,33 of q2
c   are treated implicitly
c   in the first part the rhs is calculated
c        direction x1
c
c ************************************************************************
      subroutine solq2i(al)
      include'param.f'
c  ************ compute dq2** sweeping along the x1 direction
c
      betadx=beta*dx1q*al
c
c   coeff. tridiag periodic in x1
c
        do jc=2,n2m
            ugmmm=betadx/(rc(jc)**2)
            api(jc)=-ugmmm
            ami(jc)=-ugmmm
            aci(jc)=1.+ugmmm*2.
         end do
c
        call tripvmy(1,n1m,2,n2m,1,n3m,n1m,n2)
c
      return
      end
c ************************************************************************
c   this subroutine performs the inversion of the q2 momentum equation
c   by a factored implicit scheme, only the derivatives 11,22,33 of q2
c   are treated implicitly
c   in the first part the rhs is calculated
c         direction x2
c
c ************************************************************************
      subroutine solq2j(al)
      include'param.f'
c  ********* compute the dq2* sweeping in the x2 direction
c            wall boundaries direction
      betadx=beta*al
      amj(1)=0.
      apj(1)=0.
      acj(1)=1.
      amj(n2)=0.
      apj(n2)=0.
      acj(n2)=1.
c
      do 2 jc=2,n2m
         apj(jc)=-ap2j(jc)*betadx
         acj(jc)=1.-ac2j(jc)*betadx
         amj(jc)=-am2j(jc)*betadx
    2 continue
c
      call brtrj(amj,acj,apj,n1m,n2,n3m)
c
      return
      end
c ************************************************************************
c   this subroutine performs the inversion of the q3 momentum equation
c   by a factored implicit scheme, only the derivatives 11,22,33 of q3
c   are treated implicitly
c       direction x3
c
c ************************************************************************
      subroutine solq3k(al)
      include'param.f'
c  ********* compute the dq3* sweeping in the x3 direction
c
      betadx=beta*al
c
      do 2 kc=2,n3m
        amk(kc)=-am3ck(kc)*betadx
        ack(kc)=1.-ac3ck(kc)*betadx
        apk(kc)=-ap3ck(kc)*betadx
    2 continue
c
      amk(1)=0.
      apk(1)=0.
      ack(1)=1.
      amk(n3)=0.
      apk(n3)=0.
      ack(n3)=1.
c
      call brtrk(amk,ack,apk,n1m,n2m,n3)
c
      do kc=1,n3
       do jc=1,n2m
        do ic=1,n1m
         iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2m
         q3(ic,jc,kc)=q3(ic,jc,kc) + rhs(iadd)
         rhs(iadd) = 0.
        end do
       end do
      end do
c
      return
      end
c ************************************************************************
c   this subroutine performs the inversion of the q3 momentum equation
c   by a factored implicit scheme, only the derivatives 11,22,33 of q3
c   are treated implicitly
c       direction x1
c
c ************************************************************************
      subroutine solq3i(al)
      include'param.f'
c
c  ************ compute dq3** sweeping along the x1 direction
c             
      betadx=beta*dx1q*al
        do jc=1,n2m
            ugmmm=betadx/(rm(jc)**2)
            api(jc)=-ugmmm
            ami(jc)=-ugmmm
            aci(jc)=1.+ugmmm*2.
         end do
c
        call tripvmy(1,n1m,1,n2m,1,n3m,n1m,n2m)
c
      return
      end
c ************************************************************************
c   this subroutine performs the inversion of the q3 momentum equation
c   by a factored implicit scheme, only the derivatives 11,22,33 of q3
c   are treated implicitly
c       direction x2
c
c ************************************************************************
      subroutine solq3j(al)
      include'param.f'
c
c  ************ compute dq3 sweeping along the x2 direction
c
      betadx=beta*al
      do 21 jc=1,n2m
      apj(jc)=-ap3j(jc)*betadx
      amj(jc)=-am3j(jc)*betadx
      acj(jc)=1.-ac3j(jc)*betadx
   21 continue
c
      call brtrj(amj,acj,apj,n1m,n2m,n3)
c
      return
      end
c***********************************************************************
c                                                                      *
c                       CONDIZIONI AL CONTORNO                         *
c                                                                      *
c***********************************************************************
      subroutine pscbo
      include'param.f'
c
c   in this subrout the fluxes of psc on the upper and lower wall
c    are set equal to zero
c
         do 810 j=1,n2
            do 810 i=1,n1
              pscbn(i,j)=0.
              pscbs(i,j)=0.
  810 continue
      return
      end
c ************************************************************************
c                       SUBROUTINE INVTRps
c   This subroutine performs the computation of he scalar field.
c   For details see the introduction of INVTR1
c
c ************************************************************************
      subroutine invtrps(al,ga,ro,ns)
      include'param.f'
c
      alpec=al/(ren*pran)
c
c  compute the rhs of the factored equation
c  everything at i+1/2,j+1/2,k+1/2
c
      do 18 kc=1,n3m
        km=kmv(kc)
        kp=kpv(kc)
        do 18 jc=1,n2m
          jmm=jmv(jc)
          jpp=jpv(jc)
          jp=jc+1
          aap=apscj(jc)
          aam=amscj(jc)
          aac=acscj(jc)
          udx1q = dx1q/rm(jc)**2
          do 18  ic=1,n1m
            imm=imv(ic)
            ipp=ipv(ic)
            ip=ipp
c
c   11 second derivatives of psc
c
            dq31=(psc(ipp,jc,kc)-2.*psc(ic,jc,kc)+psc(imm,jc,kc))
     %           *udx1q
c
c   22 second derivatives of psc
c
            dq32=psc(ic,jpp,kc)*aap+psc(ic,jc,kc)*aac
     %          +psc(ic,jmm,kc)*aam
c
c   33 second derivatives of psc
c
               dq33=psc(ic,jc,kp)*ap3ssk(kc)
     %              +psc(ic,jc,kc)*ac3ssk(kc)
     %              +psc(ic,jc,km)*am3ssk(kc)
            dcq3=dq31+dq32+dq33
c
c    right hand side of the pscity equation
c
c
            iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2m
c
            rhs(iadd)=(ga*qcap(ic,jc,kc)+ro*rupsc(ic,jc,kc)
     %              +alpec*dcq3)*dt
c
c    updating of the non-linear terms
c
            rupsc(ic,jc,kc)=qcap(ic,jc,kc)
   18 continue
c
      call solpsi(al)
      call solpsj(al)
      call solpsk(al)
      return
      end
c ************************************************************************
c   this subroutine performs the inversion of the q3 momentum equation
c   by a factored implicit scheme, only the derivatives 11,22,33 of q3
c   are treated implicitly
c       direction x2
c
c ************************************************************************
      subroutine solpsj(al)
      include'param.f'
c
c  ************ compute dq3 sweeping along the x2 direction
c
      betadx=0.5*al*dt/(ren*pran)
      do 21 jc=1,n2m
      apj(jc)=-apscj(jc)*betadx
      amj(jc)=-amscj(jc)*betadx
      acj(jc)=1.-acscj(jc)*betadx
   21 continue
c
c    dpsc/dr=0 has been assumed at the external boundary
c
c
      call brtrj(amj,acj,apj,n1m,n2m,n3m)
c
      return
      end
c ************************************************************************
c   this subroutine performs the inversion of the q3 momentum equation
c   by a factored implicit scheme, only the derivatives 11,22,33 of q3
c   are treated implicitly
c       direction x1
c
c ************************************************************************
      subroutine solpsi(al)
      include'param.f'
c
c  ************ compute dq3** sweeping along the x1 direction
c           
      betadx=0.5*dx1q*al*dt/(ren*pran)
        do jc=1,n2m
            ugmmm=betadx/(rm(jc)**2)
            api(jc)=-ugmmm
            ami(jc)=-ugmmm
            aci(jc)=1.+ugmmm*2.
         end do
c
        call tripvmy(1,n1m,1,n2m,1,n3m,n1m,n2m)
c
      return
      end
c ************************************************************************
c   this subroutine performs the inversion of the q3 momentum equation
c   by a factored implicit scheme, only the derivatives 11,22,33 of q3
c   are treated implicitly
c       direction x3
c
c ************************************************************************
      subroutine solpsk(al)
      include'param.f'
c  ********* compute the dq3* sweeping in the x3 direction
c
      m2v=m2
      betadx=0.5*al*dt/(ren*pran)
      ugkk=4./3.
      do 2 kc=1,n3m
        amk(kc)=-am3ssk(kc)*betadx
        ack(kc)=1.-ac3ssk(kc)*betadx
        apk(kc)=-ap3ssk(kc)*betadx
    2 continue
c
      call brtrk(amk,ack,apk,n1m,n2m,n3m)
c
      do kc=1,n3m
       do jc=1,n2m
        do ic=1,n1m
         iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2m
         psc(ic,jc,kc)=psc(ic,jc,kc) + rhs(iadd)
         rhs(iadd) = 0.
        end do
       end do
      end do
      return
      end
c ************************************************************************
c                       SUBROUTINE INVTR1
c   This subroutine performs the computation of Q~~ for the q1 momentum 
c   equation (azimuthal direction) by a factored implicit scheme.
c   Viscous terms are treated implicitly, nonlinear terms explicitly.
c   
c        ~~     n
c   dQ = Q  -  Q    
c
c         alp dt   
c   bet = -------
c          2 Re
c
c   The first equation in the introduction of THSCHEM then becomes
c   
c                
c              2      [            n         n       n-1   alp       2  n ]
c  (1-bet*nabla ) dQ =[-alp*grad (P ) + gam*H + rho*H   + ----- nabla (Q )]*dt
c                     [                                    Re             ]
c  
c                                                                      3 
c   The left hand side of this equation is then factored at an order dt 
c   as follows
c
c           2           2           2
c  (1 -bet*d_x)*(1-bet*d_r)*(1-bet*d_th) dQ = RHS
c
c           2     2       2
c  Where   d_x,  d_r and d_th  are the discrete differential operators of
c  the viscous terms in the axial, radial nad azimuthal directions
c  respectively. RHS is the righ hand side of the momentum equation as
c  is written above.            
c
      subroutine invtr1(al,ga,ro,ns)
      include'param.f'
      alre=al/ren
c
c  compute the rhs of the factored equation everything at i,j+1/2,k+1/2
c
      ugkk=4./3.
      n3mm=n3m-1
        do 25 kc=1,n3m
          km=kmv(kc)
          kp=kpv(kc)
          do 25 jc=1,n2m
             udx1=dx1/rm(jc)*al
             udx1q = dx1q/rm(jc)**2
             jp=jpv(jc)
             jm=jmv(jc)
             do 25 ic=1,n1m
                ip=ipv(ic)
                im=imv(ic)
c
c   11 second deriv. of q1(n)
c
                d11q1=(q1(ip,jc,kc)-q1(ic,jc,kc)*2.+q1(im,jc,kc))*
     %            udx1q
c
c   22 second deriv. of q1(n)
c
               d22q1=q1(ic,jp,kc)*ap1j(jc)
     %              +q1(ic,jc,kc)*ac1j(jc)
     %              +q1(ic,jm,kc)*am1j(jc)
c
c   33 second deriv. of q1(n)
c
               d33q1=q1(ic,jc,kp)*ap3sk(kc)
     %              +q1(ic,jc,kc)*ac3sk(kc)
     %              +q1(ic,jc,km)*am3sk(kc)
c 
c   all viscous terms
c
                dcq1 = d11q1 + d22q1 + d33q1
c
c   grad(pr) along 1
c
                dpx11=(pr(ic,jc,kc)-pr(im,jc,kc))*udx1
c
            iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2m
c
c   evaluation of RHS
c
            rhs(iadd)=(ga*dq(ic,jc,kc)+ro*ru1(ic,jc,kc)
     %                 +alre*dcq1-dpx11)*dt
c
c   updatind the nonlinear terms
c
                ru1(ic,jc,kc)=dq(ic,jc,kc)
   25 continue
c
c
      call solq1i(al)
      call solq1j(al)
      call solq1k(al)
      return
      end
c************************************************************************
c                       SUBROUTINE INVTR2
c   This subroutine performs the computation of Q~~ for the q2 momentum 
c   equation (radial direction) by a factored implicit scheme.
c   For details see the introduction of INVTR1
c   
c
      subroutine invtr2(al,ga,ro,ns)
      include'param.f'
      alre=al/ren
c
c
c  compute the rhs of the factored equation
c  everything at i+1/2,j,k+1/2
c
c    points inside the flowfield
c
      n3mm=n3m-1
      ugkk=4./3.
        do 26 kc=1,n3m
          km=kmv(kc)
          kp=kpv(kc)
          do 26 jc=2,n2m
            jm=jc-1
            jp=jc+1
            udx1q =  dx1q/rc(jc)**2
            udx2 = rc(jc)*al*dx2/g2rc(jc)
             do 26 ic=1,n1m
               im=imv(ic)
               ip=ipv(ic)
c
c   11 second derivative of q2
c
            d11q2=(q2(ip,jc,kc)-2.*q2(ic,jc,kc)+q2(im,jc,kc))*
     %            udx1q
c
c   33 second derivative of q2
c
               d33q2=q2(ic,jc,kp)*ap3sk(kc)
     %              +q2(ic,jc,kc)*ac3sk(kc)
     %              +q2(ic,jc,km)*am3sk(kc)
c
            dcq2 = d11q2 + d33q2
c
c   component of grad(pr) along 2 direction
c
            dpx22=(pr(ic,jc,kc)-pr(ic,jm,kc))*udx2
c
            iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2
c
            rhs(iadd)=(ga*dph(ic,jc,kc)+ro*ru2(ic,jc,kc)
     %                    +alre*dcq2-dpx22)*dt
c
            ru2(ic,jc,kc)=dph(ic,jc,kc)
   26   continue
c
c   22 second derivative of q2
c
         do kc=1,n3m
         km=kmv(kc)
         kp=kpv(kc)
         do ic=1,n1m
            insy = isym(ic)
C
C     jc = 2 (axis of symmetry)
C
            jc = 2
            jp = jc + 1
            jm = jc - 1
            q2s1 = (q2(ic,jc,kc) - q2(insy,jc,kc))*0.5/rc(jc)
            d22q2=(
     %       (q2(ic,jp,kc)/rc(jp)-q2(ic,jc,kc)/rc(jc))*rm(jc)/g2rm(jc)
     %      -(q2(ic,jc,kc)/rc(jc)-q2s1             )*rm(jm)/g2rm(jm)
     %                 )/g2rc(jc) *dx2q
     %       -q2(ic,jc,kc)/rc(jc)**2
c
               iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2
c
               rhs(iadd)=rhs(iadd) +alre*d22q2*dt
         end do
C
C     inner points
C
         do jc=3,n2m
            jm=jc-1
            jp=jc+1
            aap=ap2je(jc)
            aac=ac2je(jc)
            aam=am2je(jc)
               d22q2=q2(ic,jp,kc)*aap+q2(ic,jc,kc)*aac+q2(ic,jm,kc)*aam
c
               iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2
c
               rhs(iadd)=rhs(iadd) +alre*d22q2*dt
         end do
      end do

c
      call solq2i(al)
      call solq2j(al)
      call solq2k(al)
      return
      end
************************************************************************
c                       SUBROUTINE INVTR3
c   This subroutine performs the computation of Q~~ for the q3 momentum 
c   equation (axial direction) by a factored implicit scheme.
c   Viscous terms are treated implicitly, nonlinear terms explicitly.
c   For details see the introduction of INVTR1 
c
      subroutine invtr3(al,ga,ro,ns)
      include'param.f'
      alre=al/ren
c
c  compute the rhs of the factored equation
c  everything at i+1/2,j+1/2,k
c
      do 18 kc=2,n3m
        km=kmv(kc)
        kp=kc+1
        do 18 jc=1,n2m
          jmm=jmv(jc)
          jpp=jpv(jc)
          jp=jc+1
          aap=ap3j(jc)
          aam=am3j(jc)
          aac=ac3j(jc)
          udx1q = dx1q/rm(jc)**2
          do 18  ic=1,n1m
            imm=imv(ic)
            ipp=ipv(ic)
            ip=ipp
c
c   11 second derivatives of q3
c
            dq31=(q3(ipp,jc,kc)-2.*q3(ic,jc,kc)+q3(imm,jc,kc))*udx1q
c
c   22 second derivatives of q3
c
            dq32=aam*q3(ic,jmm,kc)+aac*q3(ic,jc,kc)+aap*q3(ic,jpp,kc)
c
c   33 second derivatives of q3
c
            dq33=q3(ic,jc,kp)*ap3ck(kc)
     %          +q3(ic,jc,kc)*ac3ck(kc)
     %          +q3(ic,jc,km)*am3ck(kc)
c
c   viscous terms
c
            dcq3=dq31+dq32+dq33
c
c  component of grad(pr) along x3 direction
c
            dpx33=(pr(ic,jc,kc)-pr(ic,jc,km))*dx3*al/g3rc(kc)
c
c  right hand side of the momentum equation
c
c
            iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2m
c
            rhs(iadd)=(ga*qcap(ic,jc,kc)+ro*ru3(ic,jc,kc)
     %                    +alre*dcq3-dpx33)*dt 
c
c  updating of the non-linear terms
c
            ru3(ic,jc,kc)=qcap(ic,jc,kc)
   18 continue
      do jc=1,n2m
        do ic=1,n1m
         iaddi=ic+(jc-1)*n1m
         iaddf=ic+(jc-1)*n1m+(n3-1)*n1m*n2m
         rhs(iaddi)=0.
         rhs(iaddf)=0.
        end do
      end do
      call solq3i(al)
      call solq3j(al)
      call solq3k(al)
      return
      end
c************************************************************************
c                                                                       *
c ****************************** subrout coetar  ********************** *
c                                                                       *
c    this subroutine calculates the coefficients for the                *
c    integration in the radial direction with non-uniform coor. trasf.  *
c                                                                       *
c************************************************************************
      subroutine coetar
      include 'param.f'
c
c  ***********  coefficients for q1   inner points
c
c
      do 151 jc=2,n2m-1
      jp=jc+1
      jm=jc-1
      a22=dx2q/g2rm(jc)/rm(jc)
      a22p= +a22*rc(jp)/g2rc(jp)
      a22m= +a22*rc(jc)/g2rc(jc)
      ap1j(jc)=a22p
      am1j(jc)=a22m
      ac1j(jc)=-(a22p+a22m)-1./rm(jc)**2
  151 continue
c
c    external bound. conditions
c
      jc=n2m
      jp=jc+1
      jm=jc-1
      a22=dx2q/g2rm(jc)/rm(jc)**2
      apjjc=a22*rc(jp)**3/rm(jc)/g2rc(jp)*2.
      acjjc=a22*rc(jc)**3/rm(jc)/g2rc(jc)
      amjjc=a22*rc(jc)**3/rm(jm)/g2rc(jc)
       if(inslwr.eq.1) then
        ap1j(jc)=-apjjc
       else
        ap1j(jc)=0.
       end if
      ac1j(jc)=-acjjc
      am1j(jc)=amjjc
c
c    internal boundary conditions
c
      jc=1
      jp=jc+1
      a22=dx2q/g2rm(jc)/rm(jc)
      a22p=a22*rc(jp)/g2rc(jp)
      ap1j(jc)=a22p
      ac1j(jc)=-a22p-1./rm(jc)**2
      am1j(jc)=0.
c
c  ***********  coefficients for q2   
c
      am2j(1)=0.
      ap2j(1)=0.
      ac2j(1)=1.
      am2j(n2)=0.
      ap2j(n2)=0.
      ac2j(n2)=1.
      do 2 jc=2,n2m
      jm=jc-1
      jp=jc+1
      a22=rc(jc)*dx2q/g2rc(jc)
      a22p=1./(rm(jc)*g2rm(jc))
      a22m=1./(rm(jm)*g2rm(jm))
      ap2j(jc)=a22*a22p
      am2j(jc)=a22*a22m
      ac2j(jc)=-(a22*a22p+a22*a22m)
    2 continue
c
c  ***********  coefficients for q2  (explicit part)
c
      do jc=3,n2m
         jm=jc-1
         jp=jc+1
         ugmm = dx2q/g2rc(jc)
         a22p=rm(jc)/(rc(jp)*g2rm(jc))
         a22m=rm(jm)/(rc(jm)*g2rm(jm))
         a22= (rm(jc)/g2rm(jc) + rm(jm)/g2rm(jm))/rc(jc)
         ap2je(jc)=a22p*ugmm
         am2je(jc)=a22m*ugmm
         ac2je(jc)=-a22*ugmm  - 1./rc(jc)**2
      end do
c
c
c  ***********  coefficients for q3   inner points
c  *********** are equal to those for psc
c
      do 51 jc=2,n2m-1
      jp=jc+1
      jm=jc-1
      a22=dx2q/g2rm(jc)/rm(jc)
      a22p= +a22*rc(jp)/g2rc(jp)
      a22m= +a22*rc(jc)/g2rc(jc)
      ap3j(jc)=a22p
      am3j(jc)=a22m
      ac3j(jc)=-(a22p+a22m)
      apscj(jc)=a22p
      amscj(jc)=a22m
      acscj(jc)=-(a22p+a22m)
   51 continue
c     
c    external bound. conditions  q3
c     
      jc=n2m
      jp=jc+1
      jm=jc-1
       ugmm2=dx2q/g2rm(jc)/rm(jc)
       am3j(jc)=rc(jc)/g2rc(jc)*ugmm2
       ap3j(jc)=0.
      if(inslwr.eq.1) then
c                 no-slip wall
       ac3j(jc)= - (ugmm2*rc(jc)/g2rc(jc) +
     %              ugmm2*rc(jp)/g2rc(jp)*2. )
      else
c                 dq3/dr=0 has been assumed at the external boundary
       ac3j(jc)= - ugmm2*rc(jc)/g2rc(jc)
      endif
c     
c    external bound. conditions  psc
c     
      a22=dx2q/g2rm(jc)/rm(jc)
      a22m= +a22*rc(jc)/g2rc(jc)
      apscj(jc)=0.
      amscj(jc)=a22m
      acscj(jc)=-a22m
c     
c    internal boundary conditions   q3
c     
      jc=1
      jp=jc+1
      a22=dx2q/g2rm(jc)/rm(jc)
      a22p= +a22*rc(jp)/g2rc(jp)
      ap3j(jc)=a22p
      am3j(jc)=0.
      ac3j(jc)=-a22p
      apscj(jc)=a22p
      amscj(jc)=0.
      acscj(jc)=-a22p
c
c
c  ***********  coefficients for q3   for x3 differentation
c  c means centered that is at k location
c
      am3ck(1)=0.
      ap3ck(1)=0.
      ac3ck(1)=1.
      am3ck(n3)=0.
      ap3ck(n3)=0.
      ac3ck(n3)=1.
      do kc=2,n3m
      km=kc-1
      kp=kc+1
      a33=dx3q/g3rc(kc)
      a33p=1./g3rm(kc)
      a33m=1./g3rm(km)
      ap3ck(kc)=a33*a33p
      am3ck(kc)=a33*a33m
      ac3ck(kc)=-(ap3ck(kc)+am3ck(kc))
      enddo
c
c
c  **coefficients for q1, q2 ap3sk,am3sk,ac3sk
c   ap3ssk,am3ssk,ac3ssk, psc   for x3 differentation
c  s means staggered that is at k+1/2 location
c
      do kc=2,n3m-1
      kp=kc+1
      km=kc-1
      a33=dx3q/g3rm(kc)
      a33p= +a33/g3rc(kp)
      a33m= +a33/g3rc(kc)
      ap3sk(kc)=a33p
      am3sk(kc)=a33m
      ac3sk(kc)=-(ap3sk(kc)+am3sk(kc))
      ap3ssk(kc)=ap3sk(kc)
      am3ssk(kc)=am3sk(kc)
      ac3ssk(kc)=ac3sk(kc)
      enddo
c     
c    lower wall  bound. conditions  indicated by inslws
c    differemtiation of sec. derivative at 1+1/2
c     
      kc=1
      kp=kc+1
      a33=dx3q/g3rm(kc)
      a33p= +a33/g3rc(kp)
      a33m= +a33/g3rc(kc)
      ap3sk(kc)=a33p
      am3sk(kc)=0.
      ac3sk(kc)=-(a33p+inslws*a33m*2.)
      ap3ssk(kc)=ap3sk(kc)
      am3ssk(kc)=am3sk(kc)
      ac3ssk(kc)=-(a33p)
c     
c    upper wall  bound. conditions  indicated by inslws
c    differemtiation of sec. derivative at n3-1/2
c     
      kc=n3m
      kp=kc+1
      a33=dx3q/g3rm(kc)
      a33p= +a33/g3rc(kp)
      a33m= +a33/g3rc(kc)
      am3sk(kc)=a33m
      ap3sk(kc)=0.
      ac3sk(kc)=-(a33m+inslwn*a33p*2.)
      ap3ssk(kc)=ap3sk(kc)
      am3ssk(kc)=am3sk(kc)
      ac3ssk(kc)=-(a33m)
      return
      end
