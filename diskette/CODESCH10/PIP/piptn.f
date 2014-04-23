c************************************************************************
c                                                                       *
c ****************************** subrout coetar  **********************  *
c                                                                       *
c    this subroutine calculates the coefficients for the              *
c    integration in the radial direction with non-uniform coor. trasf.  *
c   NEW VERSION FOR q\theta=v\theta*r
c    Bound. Cond. 1st oreder accurate
c                                                                       *
c************************************************************************
      subroutine coetar
      include 'param.f'
      common/cor1j/ap1j(m2),ac1j(m2),am1j(m2)
      common/cor2j/ap2j(m2),ac2j(m2),am2j(m2)
      common/cor3j/ap3j(m2),ac3j(m2),am3j(m2)
c   
c  *******   coefficients in several eq. funct of j
c
      do jc=1,n2m
      upd1(jc)=dx1/rm(jc)
      udx1q(jc)=dx1q/rm(jc)**2
      udxi(jc)=dx2*rc(jc)/g2rc(jc)
      udx1(jc)=dx1/rm(jc)*0.25
      udx1n(jc)=dx1/rm(jc)**2*0.25
      udx2(jc)=dx2/rm(jc)/rm(jc)/g2rm(jc)
      vd1d2(jc)=dx1/rm(jc)**2/ren*dx2/g2rm(jc)
      vd1d3(jc)=dx1/rm(jc)/ren*dx3
      vdxi(jc)=dx1q/rm(jc)**2/ren
      vdx1(jc)=2.*dx1/rm(jc)/ren
      ugmv(jc)=1./rm(jc)
      udh2(jc)=dx2/g2rm(jc)*0.25
      vh23(jc)=dx2*dx3/rm(jc)/ren/g2rm(jc)
      vh13(jc)=dx1*dx3/rm(jc)/ren
      volz(jc)=rm(jc)*g2rm(jc)/(dx1*dx2*dx3)
      enddo
      do jc=2,n2m
      a11(jc)=dx1q/rc(jc)**2
      uvdx1(jc)=2.*dx1/rc(jc)/ren
      vdxh2(jc)=dx1/rc(jc)/ren
      upd2(jc)=dx2*rc(jc)/g2rc(jc)
      udx1c(jc)=dx1/rc(jc)*0.25
      udx2c(jc)=1./g2rc(jc)*dx2*0.25
      vdyh2(jc)=rc(jc)*dx2*dx1/ren/g2rc(jc)
      vdyh3(jc)=rc(jc)*dx2*dx3/ren/g2rc(jc)
      volz2(jc)=rc(jc)*g2rc(jc)/(dx1*dx2*dx3)
      enddo
c
c  ***********  coefficients for q1=qt   inner points
c
c
      do 151 jc=2,n2m-1
      jp=jc+1
      jm=jc-1
      a22=dx2q/g2rm(jc)
      a22p= +a22*rc(jp)/g2rc(jp)
      a22m= +a22*rc(jc)/g2rc(jc)
      ap1j(jc)=a22p/rm(jp)
      am1j(jc)=a22m/rm(jm)
      ac1j(jc)=-(a22p+a22m)/rm(jc)-1./rm(jc)**2
  151 continue
c
c   at  r=0    rd/dr=0
c
      jc=1
      jp=jc+1
      a22=dx2q/g2rm(jc)
      a22p=a22*rc(jp)/g2rc(jp)
      ap1j(jc)=a22p/rm(jp)
      ac1j(jc)=-a22p/rm(jc)-1./rm(jc)**2
      am1j(jc)=0.
      jc=n2m
      jm=jc-1
      jp=jc+1
      if(islip.eq.0) then

c
c   at external boundary r d/dr(1/r) qt
c
      a22=dx2q/g2rm(jc)*rm(jc)
      ap1j(jc)=-a22/rc(jp)/g2rc(jp)*2.
      ac1j(jc)=-a22/rc(jc)/g2rc(jc)
      am1j(jc)=a22/rc(jc)/g2rc(jc)
                     else
c
c   at external boundary 1/r d/dr(r^3 d(qt/r^2)/dr)  FREE-SLIP
c   with    rd(qt/r^2)/dr=0 at the wall
c
      jc=n2m
      jm=jc-1
      a22=dx2q/(g2rm(jc)*g2rc(jc)*rm(jc))*rc(jc)**3
      ap1j(jc)=0.
      ac1j(jc)=-a22/rm(jc)**2
      am1j(jc)=a22/rc(jc)**2
                      endif
c
c  ***********  coefficients for q2   inner points
c   d^2/dr^2 (q2) -1/r d/dr(q2)
c

      am2j(1)=0.                                                         
      ap2j(1)=0.                                                         
      ac2j(1)=0.                                                         
      am2j(n2)=0.                                                        
      ap2j(n2)=0.                                                        
      ac2j(n2)=0.                                                        
      do  jc=2,n2m                                                     
      jm=jc-1
      a22=dx2q/g2rc(jc)
      a22dp=dx2/(2.*rc(jc)*g2rc(jc))
      ap2j(jc)=a22/g2rm(jc)-a22dp
      am2j(jc)=a22/g2rm(jm)+a22dp
      ac2j(jc)=-(a22/g2rm(jc)+a22/g2rm(jm))
       enddo                                                            
c
c  ***********  coefficients for q3   inner points
c
      do jc=2,n2m-1
      jp=jc+1                                                           
      a22=dx2q/g2rm(jc)/rm(jc)
      a22p= +a22*rc(jp)/g2rc(jp)
      a22m= +a22*rc(jc)/g2rc(jc)
      ap3j(jc)=a22p
      am3j(jc)=a22m
      ac3j(jc)=-(a22p+a22m)
      enddo                                                             
c                                                                       
c    r=0 gives the following b.c. at axis  equiv. to dq3/dr=0      
c                                                                       
      jc=1
      jp=jc+1
      ugmm2=dx2q/g2rm(jc)/rm(jc)
      am3j(jc)=0.
      ac3j(jc)=ugmm2*rc(jp)/g2rc(jp)
      ap3j(jc)=(rc(jp)/g2rc(jp))*ugmm2
      jc=n2m
      jp=jc+1
      if(islip.eq.0) then
c                                                                       
c    q3=0 has been assumed at the wall boundary                 
c                                                                       
      ugmm2=dx2q/g2rm(jc)/rm(jc)
      am3j(jc)=rc(jc)/g2rc(jc)*ugmm2
      ac3j(jc)=ugmm2*rc(jc)/g2rc(jc)
      ap3j(jc)=ugmm2*rc(jp)/g2rc(jp)*2.
                     else
c
c    dq3/dr=0  at the wall boundary  FREE-SLIP
c
      jc=n2m
      jm=jc-1
      a22=dx2q/(g2rm(jc)*g2rc(jc)*rm(jc))*rc(jc)
      ap3j(jc)=0.
      ac3j(jc)=-a22
      am3j(jc)=a22
                      endif
      return                                                            
      end                                                               
c************************************************************************
c                                                                       *
c ****************************** subrout updvp  **********************  *
c                                                                       *
c    this subroutine calculates the solenoidal vel field.               *
c    q(n+1)=qhat-grad(dph)*dt ,  pr=dph                                 *
c    third order runge-kutta is used.                                   *
c                                                                       *
c************************************************************************
      subroutine updvp(al)                         
      include 'param.f'
c                                                                       
c  ***********  compute the q1 velocity component                       
c               v1dgf=component 1 of grad(dph)                          
      do kc=1,n3m                                                     
         do jc=1,n2m                                                     
            do ic=1,n1m                                                   
               im=imv(ic)                                                        
      dfx11=(dph(ic,jc,kc)-dph(im,jc,kc))
      v1dvgf=dfx11*dx1*al                                            
      q1(ic,jc,kc)=q1(ic,jc,kc)-v1dvgf*dt                              
            enddo
         enddo
      enddo
c                                                                       
c                                                                       
c  ***********  compute the q2 velocity component                       
c               v2dgf=component 2 of grad(dph)                          
      do kc=1,n3m                                                     
         do jc=2,n2m                                                     
         jm=jc-1                                                           
            do ic=1,n1m                                                   
      dfx22=(dph(ic,jc,kc)-dph(ic,jm,kc))
      v2dvgf=(dfx22*upd2(jc))*al                                          
      q2(ic,jc,kc)=q2(ic,jc,kc)-v2dvgf*dt                              
            enddo
         enddo
      enddo
      do  kc=1,n3m                                                    
          do  ic=1,n1m                                                    
      q2(ic,1,kc)=0.                                                    
      q2(ic,n2,kc)=0.                                                   
         enddo
      enddo
c                                                                       
c  ***********  compute the q3 velocity component                       
c               q3 is the cartesian component                           
c               v3dgf=component 3 of grad(dph)                          
c    
      if(n3m.gt.1) then
      do kc= 1,n3m                                                
         km=kmv(kc)                                                           
         do jc=1,n2m                                                     
            do ic=1,n1m                                                   
      dfx33=(dph(ic,jc,kc)-dph(ic,jc,km))*dx3                           
      v3dvgf=dfx33*al                                                   
      q3(ic,jc,kc)=q3(ic,jc,kc)-v3dvgf*dt
            enddo
         enddo
      enddo
                    endif
      return                                                            
      end                                                               
c************************************************************************
c                                                                       *
c  ************************ subrout invtu1  **********************      *
c                                                                       *
c  this subroutine performs the inversion of the q1 momentum equation   *
c  by a factored implicit scheme, only the derivatives 11,22,33 of q1   *
c  are treated implicitly                                               *
c                                                                       *
c************************************************************************
      subroutine invtu1(al,ga,ro,ns)                      
      include 'param.f'
      common/cor1j/ap1j(m2),ac1j(m2),am1j(m2)
      alre=al/ren                                                       
c                                                                       
c  ********* compute the rhs of the factored equation                   
c  the rhs includes h(n), h(n-1)=ru, grad(pr) component                 
c  and the 11, 22,3 second derivat. of q1(n)                            
c  everything at i,j+1/2,k+1/2                                          
c  dq=qhat-q(n)                                                         
c                                                                       
      do kc=1,n3m                                                   
      do jc=1,n2m                                                    
      do ic=1,n1m                                                   
      ip=ipv(ic)                                                        
      im=imv(ic)                                                        
c                                                                       
c   11 second deriv. of q1(n)                                           
c                                                                       
      d11q1=(q1(ip,jc,kc)
     1      -q1(ic,jc,kc)*2.
     1      +q1(im,jc,kc))*udx1q(jc)
c                                                                       
c   grad(pr) along 1                                                    
c                                                                       
      dpx11=(pr(ic,jc,kc)-pr(im,jc,kc))*dx1
      gradp=dpx11*al
         iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2m
      rhs(iadd)=(ga*dq(ic,jc,kc)+ro*ru1(ic,jc,kc)                   
     1             +alre*d11q1-gradp)*dt
      ru1(ic,jc,kc)=dq(ic,jc,kc)                                        
      enddo                                                             
      enddo                                                             
      enddo                                                             
c                                                                       
c   add the 33 derivat. 
c                                                                       
      if(n3m.gt.1) then
      do kc=1,n3m
      km=kmv(kc)                                                        
      kp=kpv(kc)                                                        
      do jc=1,n2m
      do ic=1,n1m                                                   
      d33q1=(q1(ic,jc,kp)
     1      -q1(ic,jc,kc)*2.
     1      +q1(ic,jc,km))*dx3q
         iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2m
      rhs(iadd)=rhs(iadd)+alre*d33q1*dt
      enddo                                                             
      enddo                                                             
      enddo                                                             
                   endif
c
c   add the 22 derivat. 
c
      do kc=1,n3m
      do jc=2,n2m-1
      jp=jc+1
      jm=jc-1
      do ic=1,n1m                                                   
      d22q1=q1(ic,jp,kc)*ap1j(jc)
     1     +q1(ic,jc,kc)*ac1j(jc)
     1     +q1(ic,jm,kc)*am1j(jc)
         iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2m
      rhs(iadd)=rhs(iadd)+alre*d22q1*dt
      enddo                                                             
      enddo                                                             
      enddo                                                             
c
c    outer wall no-slip
c

      jc=n2m
      jm=jc-1
      do kc=1,n3m
      do ic=1,n1m                                                   
      d22q1=+q1(ic,jc,kc)*(ac1j(jc)+ap1j(jc))
     1      +q1(ic,jm,kc)*am1j(jc)
         iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2m
      rhs(iadd)=rhs(iadd)+alre*d22q1*dt
      enddo                                                             
      enddo                                                             
c

c    axis r=0
c

      jc=1
      jp=jc+1
      do kc=1,n3m
      do ic=1,n1m                                                   
      d22q1=q1(ic,jp,kc)*ap1j(jc)
     1     +q1(ic,jc,kc)*ac1j(jc)
         iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2m
      rhs(iadd)=rhs(iadd)+alre*d22q1*dt
      enddo                                                             
      enddo                                                             
c
c
c
      call solq1i(ns) 
c
      if(n3m.gt.1) then
      call solq1k(al)                                                   
                   endif
c
c  the numeration of rhs is changed to speed-up the tridiagonal
c  solver in j  btrjik with a single do instead of a do i and a do k
c
      do kc=1,n3m
      do jc=1,n2m
      do ic=1,n1m                                                   
         iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2m
      dum(ic,jc,kc)=rhs(iadd)
      enddo                                                             
      enddo                                                             
      enddo                                                             
      do kc=1,n3m
      do jc=1,n2m
      do ic=1,n1m                                                   
         iadd=ic+(kc-1)*n1m+(jc-1)*n1m*n3m
      rhs(iadd)=dum(ic,jc,kc)
      enddo                                                             
      enddo                                                             
      enddo                                                             
      call solq1j(al)                                            
      return                                                            
      end                                                               
c************************************************************************
c                                                                       *
c ****************************** subrout solq1k  ********************** *
c                                                                       *
c  this subroutine performs the inversion of the q1 momentum equation   *
c  by a factored implicit scheme, only the derivatives 11,22,33 of q1   *
c  are treated implicitly                                               *
c      direction x3                                                     *
c                                                                       *
c************************************************************************
      subroutine solq1k(al)                                             
      include 'param.f'
c
c    Solve for the equation:
c
c      (       al * dt    d     d  )
c      ( 1 -  --------   ---    --- ) ( D q_t)  =  RHS
c      (       2 * Re   (d z)   dz )
c
c                                                                       
c  ************ compute  from dq** sweeping along the x3 direction      
c     periodic                                                          
      betadx=beta*dx3q*al
      do 4 kc=1,n3m
      amk(kc)=-betadx
      ack(kc)=1.+betadx*2.
      apk(kc)=-betadx
    4 continue
      call tripkji(1,n3m,1,n2m,1,n1m,n1m,n2m)
      return                                                            
      end                                                               
c************************************************************************
c                                                                       *
c ****************************** subrout solq1i  ********************** *
c                                                                       *
c  this subroutine performs the inversion of the q1 momentum equation   *
c  by a factored implicit scheme, only the derivatives 11,22,33 of q1   *
c  are treated implicitly                                               *
c        direction x1                                                   *
c                                                                       *
c************************************************************************
      subroutine solq1i(ns)
      include 'param.f'
c
c    Solve for the equation:
c
c
c      (       al * dt [ d   d ] )
c      ( 1 -  -------- [ ------] ) ( q_t) =  RHS
c      (       2 * Re  [ dt dt ] )
c
c  ********* compute dq1*  sweeping in the x1 direction
c 
      do kc=1,n3m
           do ic=1,n1m
               do jc=1,n2m
         iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2m
      fi(jc,ic)=rhs(iadd)
               enddo
            enddo
      call trpv1ij(1,n1m,1,n2m,ns)
               do ic=1,n1m
            do jc=1,n2m
         iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2m
      rhs(iadd)=fi(jc,ic)
            enddo
                enddo
      enddo
      return
      end
c************************************************************************
c                                                                       *
c ****************************** subrout solq1j  ********************** *
c                                                                       *
c  this subroutine performs the inversion of the q1 momentum equation   *
c  by a factored implicit scheme, only the derivatives 11,22,33 of q1   *
c  are treated implicitly                                               *
c                                                                       *
c************************************************************************
      subroutine solq1j(al)                                      
      include 'param.f'
      common/cor1j/ap1j(m2),ac1j(m2),am1j(m2)
      dimension apj(m2),acj(m2),amj(m2)
c
c    Solve for the equation:
c
c
c      (       al * dt [ d     d     1  ] )
c      ( 1 -  -------- [ --- r---- - ---] ) ( q_t/r) =  RHS
c      (       2 * Re  [ dr    d r   r  ] )
c                                                                       
c  ************ compute dq1  sweeping along the x2 direction            
c               wall boundaries interior points                         
c                                                                       
c
      betadx=beta*al
      do 21 jc=1,n2m-1
      apj(jc)= -ap1j(jc)*betadx
      acj(jc)=1.-ac1j(jc)*betadx
      amj(jc)= -am1j(jc)*betadx
   21 continue
      jc=n2m
      apj(jc)= 0.
      acj(jc)=1.-(ac1j(jc)+ap1j(jc))*betadx
      amj(jc)= -am1j(jc)*betadx
c
      call btrjik(amj,acj,apj,n1m,n2m,n3m)
c
      do kc=1,n3m
       do jc=1,n2m
        do ic=1,n1m
         iadd=ic+(kc-1)*n1m+(jc-1)*n1m*n3m
         q1(ic,jc,kc)=q1(ic,jc,kc) + rhs(iadd)
         rhs(iadd) = 0.
        end do
       end do
      end do
      return                                                            
      end                                                               
c************************************************************************
c                                                                       *
c ****************************** subrout invtu2  ********************** *
c  this subroutine performs the inversion of the q2 momentum equation   *
c  by a factored implicit scheme, only the derivatives 11,22,33 of q2   *
c  are treated implicitly                                               *
c  in the first part the rhs is calculated                              *
c                                                                       *
c************************************************************************
      subroutine invtu2(al,ga,ro,ns)                      
      include 'param.f'
      common/cor2j/ap2j(m2),ac2j(m2),am2j(m2)
      alre=al/ren                                                       
c                                                                       
c                                                                       
c  ********* compute the rhs of the factored equation                   
c  the rhs includes h(n), h(n-1)=ru, grad(pr) component                 
c  and the 11, 22,3 second derivat. of q2(n)                            
c  everything at i,j+1/2,k+1/2                                          
c                                                                       
c    points inside the flowfield                                        
c                                                                       
      do kc=1,n3m                                                    
            do jc=2,n2m                                                   
            jm=jc-1                                                           
            jp=jc+1                                                           
                do ic=1,n1m                                                   
                im=imv(ic)                                                        
                ip=ipv(ic)                                                        
c
c   11 second derivative of q2
c
      d11q2=(q2(ip,jc,kc)
     1      -q2(ic,jc,kc)*2.
     1      +q2(im,jc,kc))*a11(jc)
c                                                                       
c   22 second derivative of q2                                          
c                                                                       
      d22q2=q2(ic,jp,kc)*ap2j(jc)
     1     +q2(ic,jm,kc)*am2j(jc)
     1     +q2(ic,jc,kc)*ac2j(jc)
c                                                                       
c   add 33 second derivative of q2                                      
c                                                                       
       dcq2=(d11q2+d22q2)                                                
c                                                                       
c   component of grap(pr) along 2 direction                             
c                                                                       
      dpx22=(pr(ic,jc,kc)-pr(ic,jm,kc))*dx2/g2rc(jc)
      gradp=dpx22*al*rc(jc)
         iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2
      rhs(iadd)=(ga*dph(ic,jc,kc)+ro*ru2(ic,jc,kc)                   
     1             +alre*dcq2)*dt-gradp*dt                              
      ru2(ic,jc,kc)=dph(ic,jc,kc)                                        
                 enddo
             enddo
      enddo
c                                                                       
c
c                                                                       
c   add 33 second derivative of q2                                      
c                                                                       
      if(n3m.gt.1) then
      do kc=1,n3m
      km=kmv(kc)                                                        
      kp=kpv(kc)                                                        
           do jc=2,n2m
               do ic=1,n1m
      d33q2=(q2(ic,jc,kp)
     1      -q2(ic,jc,kc)*2.    
     1      +q2(ic,jc,km))*dx3q
         iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2
      rhs(iadd)=rhs(iadd)+alre*d33q2*dt
                enddo
             enddo
      enddo
                   endif
c
      call solq2i(ns) 
c
      if(n3m.gt.1) then
      call solq2k(al)                                                   
                   endif
c
c  the numeration of rhs is changed to speed-up the tridiagonal
c  solver in j  btrjik with a single do instead of a do i and a do k
c
      do kc=1,n3m
      do jc=1,n2
      do ic=1,n1m                                                   
         iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2
      dum(ic,jc,kc)=rhs(iadd)
      enddo                                                             
      enddo                                                             
      enddo                                                             
      do kc=1,n3m
      do jc=1,n2
      do ic=1,n1m                                                   
         iadd=ic+(kc-1)*n1m+(jc-1)*n1m*n3m
      rhs(iadd)=dum(ic,jc,kc)
      enddo                                                             
      enddo                                                             
      enddo                                                             
      call solq2j(al)                                            
c
      return                                                            
      end                                                               
c************************************************************************
c                                                                       *
c ****************************** subrout solq2k  ********************** *
c  this subroutine performs the inversion of the q2 momentum equation   *
c  by a factored implicit scheme, only the derivatives 11,22,33 of q2   *
c  are treated implicitly                                               *
c  in the first part the rhs is calculated                              *
c       direction x3                                                    *
c                                                                       *
c************************************************************************
      subroutine solq2k(al)                                             
      include 'param.f'
c
c    Solve for the equation:
c
c      (       al * dt    d  d   )
c      ( 1 -  --------   --- --- ) ( D q_r)  =  RHS
c      (       2 * Re   (d z dz  )
c
c                                                                       
c  ************ compute dq2 sweeping along the x3 direction             
c               periodic                                                
      betadx=beta*dx3q*al
      do 4 kc=1,n3m
      amk(kc)=-betadx
      ack(kc)=1.+betadx*2.
      apk(kc)=-betadx
    4 continue
      call tripkji(1,n3m,1,n2m,1,n1m,n1m,n2)
      return                                                            
      end                                                               
c************************************************************************
c                                                                       *
c ****************************** subrout solq2i  ********************** *
c  this subroutine performs the inversion of the q2 momentum equation   *
c  by a factored implicit scheme, only the derivatives 11,22,33 of q2   *
c  are treated implicitly                                               *
c  in the first part the rhs is calculated                              *
c       direction x1                                                    *
c                                                                       *
c************************************************************************
      subroutine solq2i(ns)
      include 'param.f'
c
c    Solve for the equation:
c
c      (       al * dt  1     d     d  )
c      ( 1 -  -------- ---   ----  --- ) (D q_r)  =  RHS
c      (       2 * Re   r^2  (d t   dt )
c
c  ************ compute dq2** sweeping along the x1 direction
c
c     
      do kc=1,n3m
c
c   coeff. tridiag periodic in x1
c
            do ic=1,n1m
                  do jc=2,n2m
         iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2
      fi(jc,ic)=rhs(iadd)
      enddo
      enddo
      call trpv2ij(1,n1m,2,n2m,ns)
            do jc=2,n2m
                  do ic=1,n1m
         iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2
      rhs(iadd)=fi(jc,ic)
                  enddo
            enddo
      enddo
      return
      end
c************************************************************************
c                                                                       *
c ****************************** subrout solq2j  ********************** *
c  this subroutine performs the inversion of the q2 momentum equation   *
c  by a factored implicit scheme, only the derivatives 11,22,33 of q2   *
c  are treated implicitly                                               *
c  in the first part the rhs is calculated                              *
c        direction x2                                                   *
c                                                                       *
c************************************************************************
      subroutine solq2j(al)                                      
      include 'param.f'
      common/cor2j/ap2j(m2),ac2j(m2),am2j(m2)
      dimension apj(m2),acj(m2),amj(m2)
c
c    Solve for the equation:
c
c      (       al * dt [ d^2     1  d   ] )
c      ( 1 -  -------- [ ---   - -  --  ] ) ( q_r)  =  RHS
c      (       2 * Re  [ dr^2    r  d r ] )
c
c  ********* compute the dq2* sweeping in the x2 direction              
c            wall boundaries direction                                  
      betadx=beta*al
      do jc=1,n2,n2m
      apj(jc)=0.
      acj(jc)=1.
      amj(jc)=0.
      enddo                                                             
            do jc=2,n2m
      apj(jc)=-ap2j(jc)*betadx
      acj(jc)=1.-betadx*ac2j(jc)
      amj(jc)=-am2j(jc)*betadx
             enddo                                                             
      do kc=1,n3m
       do jc=1,n2,n2m
        do ic=1,n1m
         iadd=ic+(kc-1)*n1m+(jc-1)*n1m*n3m
         rhs(iadd) = 0.
        end do
       end do
      end do
c
      call btrjik(amj,acj,apj,n1m,n2,n3m)
c
      do kc=1,n3m
       do jc=1,n2
        do ic=1,n1m
         iadd=ic+(kc-1)*n1m+(jc-1)*n1m*n3m
         q2(ic,jc,kc)=q2(ic,jc,kc) + rhs(iadd)
         rhs(iadd) = 0.
        end do
       end do
      end do
        do ic=1,n1m
            do kc=1,n3m                                                     
      q2(ic,1,kc)=0.
      q2(ic,n2,kc)=0.
            enddo                                                             
        enddo                                                             
      return                                                            
      end                                                               
c************************************************************************
c                                                                       *
c ****************************** subrout invtu3  ********************** *
c  this subroutine performs the inversion of the q3 momentum equation   *
c  by a factored implicit scheme, only the derivatives 11,22,33 of q3   *
c  are treated implicitly                                               *
c                                                                       *
c************************************************************************
      subroutine invtu3(al,ga,ro,ns)                      
      include 'param.f'
      common/cor3j/ap3j(m2),ac3j(m2),am3j(m2)
      common/timw/timew
      alre=al/ren                                                       
c                                                                       
c  ********* compute the rhs of the factored equation                   
c  the rhs includes h(n), h(n-1)=ru, grad(pr) component                 
c  and the 11, 22,3 second derivat. of q3(n)                            
c  everything at i,j+1/2,k+1/2                                          
c                                                                       
      s3tot=0.
      q3tot=0.
      slorf=0.
      do kc=1,n3m
            do jc=1,n2m
                  do ic=1,n1m
		  q3tot=q3tot+q3(ic,jc,kc)*volz(jc)
		  enddo
            enddo
      enddo
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
      dq31=(q3(ip,jc,kc)
     1     -q3(ic,jc,kc)*2.
     1     +q3(im,jc,kc))*udx1q(jc)
c
c   33 second derivatives of q3
c
      dq33=(q3(ic,jc,kp)
     1     -q3(ic,jc,kc)*2.
     1     +q3(ic,jc,km))*dx3q
      dcq3=dq31+dq33                                               
         iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2m
      slorf=slorf+(ga*qcap(ic,jc,kc)+ro*ru3(ic,jc,kc))*volz(jc)/al
      rhs(iadd)=(ga*qcap(ic,jc,kc)+ro*ru3(ic,jc,kc)
     1              +alre*dcq3)*dt
      ru3(ic,jc,kc)=qcap(ic,jc,kc)
                  enddo                                                             
            enddo                                                             
      enddo                                                             
c     c
c   add second derivatie in r
c
      do kc=1,n3m
            do jc=2,n2m-1
            jm=jmv(jc)                                                       
            jp=jpv(jc)                                                       
                  do ic=1,n1m
c
c   22 second derivatives of q3
c
      dq32=(q3(ic,jp,kc)*ap3j(jc)
     1     -q3(ic,jc,kc)*(ap3j(jc)+am3j(jc))
     1     +q3(ic,jm,kc)*am3j(jc))*alre
         iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2m
      rhs(iadd)=rhs(iadd)+dq32*dt
      s3tot=s3tot+dq32*volz(jc)/al
                  enddo  
            enddo                                                             
      enddo                                                             
c
c   22 second derivatives of q3  at r=n2m
c
      jc=n2m
      jm=jc-1                                                      
      do kc=1,n3m
                  do ic=1,n1m
c
c   22 second derivatives of q3
c
      dq32=(-q3(ic,jc,kc)*(ap3j(jc)+ac3j(jc))
     1      +q3(ic,jm,kc)*am3j(jc))*alre
         iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2m
      rhs(iadd)=rhs(iadd)+dq32*dt
      s3tot=s3tot+dq32*volz(jc)/al
                  enddo                                                             
      enddo                                                             
c     c
c   add second derivatie in r at r=0
c
      jc=1
      jp=jc+1                                                      
      do kc=1,n3m
                  do ic=1,n1m
c
c   22 second derivatives of q3
c
      dq32=(q3(ic,jp,kc)*ap3j(jc)
     1     -q3(ic,jc,kc)*ap3j(jc))*alre
         iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2m
      rhs(iadd)=rhs(iadd)+dq32*dt
      s3tot=s3tot+dq32*volz(jc)/al
                  enddo                                                             
      enddo                                                             
       dp3ns=s3tot/(pi*alx3d)
       dplfo=slorf/(pi*alx3d)
       dptot=dplfo+dp3ns
      do kc=1,n3m
      km=kmv(kc)
            do jc=1,n2m
                  do ic=1,n1m
c                                                                       
c  component of grad(pr) along x3 direction                             
c                                                                       
      dpx33=(pr(ic,jc,kc)-pr(ic,jc,km))*dx3
      gradp=(dpx33+dptot)*al
         iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2m
      rhs(iadd)=rhs(iadd)-gradp*dt
                  enddo  
            enddo                                                             
      enddo                                                             
       q3tot=q3tot/(pi*alx3d)
       if(ns.eq.3) then
       write(95,109) timew,dp3ns,dplfo,dptot,q3tot
  109 format(3x,5e12.4)
                   endif

c
c
      call solq3i(ns) 
c
      call solq3k(al)                                                   
c
c  the numeration of rhs is changed to speed-up the tridiagonal
c  solver in j  btrjik with a single do instead of a do i and a do k
c
      do kc=1,n3m
      do jc=1,n2m
      do ic=1,n1m                                                   
         iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2m
      dum(ic,jc,kc)=rhs(iadd)
      enddo                                                             
      enddo                                                             
      enddo                                                             
      do kc=1,n3m
      do jc=1,n2m
      do ic=1,n1m                                                   
         iadd=ic+(kc-1)*n1m+(jc-1)*n1m*n3m
      rhs(iadd)=dum(ic,jc,kc)
      enddo                                                             
      enddo                                                             
      enddo                                                             
      call solq3j(al)                                            
      return                                                            
      end                                                               
c************************************************************************
c                                                                       *
c ****************************** subrout solq3k  ********************** *
c  this subroutine performs the inversion of the q3 momentum equation   *
c  by a factored implicit scheme, only the derivatives 11,22,33 of q3   *
c  are treated implicitly                                               *
c      direction x3                                                     *
c                                                                       *
c************************************************************************
      subroutine solq3k(al)                                             
      include 'param.f'
c
c    Solve for the equation:
c
c      (       al * dt    d   d    )
c      ( 1 -  --------   --  ----- ) ( D q_z)  =  RHS
c      (       2 * Re   (d z  dz   )
c
c  ********* compute the dq3* sweeping in the x3 direction              
c                                                                       
      betadx=beta*dx3q*al
      do 4 kc=1,n3m
      amk(kc)=-betadx
      ack(kc)=1.+betadx*2.
      apk(kc)=-betadx
    4 continue
      call tripkji(1,n3m,1,n2m,1,n1m,n1m,n2m)
      return                                                            
      end                                                               
c************************************************************************
c                                                                       *
c ****************************** subrout solq3i  ********************** *
c  this subroutine performs the inversion of the q3 momentum equation   *
c  by a factored implicit scheme, only the derivatives 11,22,33 of q3   *
c  are treated implicitly                                               *
c      direction x1                                                     *
c                                                                       *
c*************************************************************************
      subroutine solq3i(ns)
      include 'param.f'
c
c    Solve for the equation:
c
c      (       al * dt  1     d    d    )
c      ( 1 -  -------- ---   ----  ---  ) (D q_z)  =  RHS
c      (       2 * Re   r^2  (d t  dt   )
c
c  ************ compute dq3** sweeping along the x1 direction
c               inflow outflow
c

      do kc=1,n3m
            do jc=1,n2m
                  do ic=1,n1m
         iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2m
      fi(jc,ic)=rhs(iadd)
                  enddo
            enddo
c
c   coeff tridiag periodic in x1
c
      call trpv1ij(1,n1m,1,n2m,ns)
            do jc=1,n2m
                  do ic=1,n1m
         iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2m
      rhs(iadd)=fi(jc,ic)
                  enddo
            enddo
      enddo
c
      return   
      end
c************************************************************************
c                                                                       *
c ****************************** subrout solq3j  ********************** *
c  this subroutine performs the inversion of the q3 momentum equation   *
c  by a factored implicit scheme, only the derivatives 11,22,33 of q3   *
c  are treated implicitly                                               *
c      direction x2                                                     *
c                                                                       *
c************************************************************************
      subroutine solq3j(al)                                      
      include 'param.f'
      common/cor3j/ap3j(m2),ac3j(m2),am3j(m2)
      dimension apj(m2),acj(m2),amj(m2)
c
c    Solve for the equation:
c
c      (       al * dt  1  [ d       d   ] )
c      ( 1 -  -------- --- [ --- r ----  ] ) ( q_z)  =  RHS
c      (       2 * Re   r  [ dr      d r ] )
c
c  ************ compute dq3 sweeping along the x2 direction             
c                                                                       
      betadx=beta*al
            do jc=2,n2m-1
      apj(jc)=-ap3j(jc)*betadx
      amj(jc)=-am3j(jc)*betadx
      acj(jc)=1.-(apj(jc)+amj(jc)) 
            enddo                                                             
      jc=1
      apj(jc)=-ap3j(jc)*betadx
      amj(jc)=0.
      acj(jc)=1.-apj(jc)
      jc=n2m                                                    
      apj(jc)=0.
      amj(jc)=-am3j(jc)*betadx
      acj(jc)=1.+(ac3j(jc)+ap3j(jc))*betadx
      call btrjik(amj,acj,apj,n1m,n2m,n3m)
c
      do kc=1,n3m
       do jc=1,n2m
        do ic=1,n1m
         iadd=ic+(kc-1)*n1m+(jc-1)*n1m*n3m
         q3(ic,jc,kc)=q3(ic,jc,kc) + rhs(iadd)
         rhs(iadd) = 0.
        end do
       end do
      end do
      return                                                            
      end                                                               
c************************************************************************
                                                                        *
       subroutine prcalc(al)                                            *
c************************************************************************ 
      include 'param.f'
c                                                                       
c    the pressure is evaluated at the center of the box.                
c    pr=pressure.                                                       
c    remember that pressure correction has only 11,22,33 derivatives bec
c    of the treatment of implicit method in q1,q2,q3.                   
c                                                                       
      be=al*beta                                                        
      do 1 kc=1,n3m                                                     
      kp=kpv(kc)                                                        
      km=kmv(kc)                                                        
      do 1 jc=1,n2m                                                     
      jm=jmv(jc)                                                        
      jp=jpv(jc)                                                        
      ugmmm=1./(g2rm(jc)*rm(jc))*dx2q
      do 1 ic=1,n1m                                                   
      ip=ipv(ic)                                                        
      im=imv(ic)                                                        
      pr(ic,jc,kc)=pr(ic,jc,kc)+dph(ic,jc,kc)-be*                       
     1  ( (dph(ip,jc,kc)-2.*dph(ic,jc,kc)+dph(im,jc,kc))*dx1q/rm(jc)**2
     1   +( (dph(ic,jp,kc)-dph(ic,jc,kc))*rc(jp)/g2rc(jp)               
     1     -(dph(ic,jc,kc)-dph(ic,jm,kc))*rc(jc)/g2rc(jc) )*ugmmm
     1   +(dph(ic,jc,kp)-2.*dph(ic,jc,kc)+dph(ic,jc,km))*dx3q)    
    1 continue                                                          
      return                                                            
      end                                                               
c************************************************************************
c                                                                       *
c ****************************** subrout csoq13i  ********************** *
c                                                                       *
c  this subroutine performs the valculation of the coefficients for     *
c  solq1i and solq13 to spped-up the tridiagonal inversion              *
c                                                                       *
c************************************************************************
      subroutine csoq13i
      include 'param.f'
c
c    Solve for the equation:
c
c      (       al * dt  1     d   d    )
c      ( 1 -  -------- ---   ---  ---  ) (D q_t)  =  RHS
c      (       2 * Re   r^2  (dt  dt   )
c
c
c  ********* compute dq1*  sweeping in the x1 direction
c            inflow outflow
      do ns=1,nsst
      al=alm(ns)
      betadx=beta*al
           do ic=1,n1m
               do jc=1,n2m
      acc=betadx*udx1q(jc)
      apif(jc,ic,ns)=-acc
      amif(jc,ic,ns)=-acc
      acif(jc,ic,ns)=1.-(apif(jc,ic,ns)+amif(jc,ic,ns))
               enddo
            enddo
      call ctpv1ij(1,n1m,1,n2m,ns)
      enddo
      return
      end
c************************************************************************
c                                                                       *
c ****************************** subrout csoq2i  ********************** *
c                                                                       *
c  this subroutine performs the valculation of the coefficients for     *
c  solq2i  to spped-up the tridiagonal inversion              *
c                                                                       *
c************************************************************************
      subroutine csoq2i
      include 'param.f'
c
c    Solve for the equation:
c
c      (       al * dt  1     d   d    )
c      ( 1 -  -------- ---   ---  ---  ) (D q_t)  =  RHS
c      (       2 * Re   r^2  (dt  dt   )
c
c
c  ********* compute dq1*  sweeping in the x1 direction
c            inflow outflow
      do ns=1,nsst
      al=alm(ns)
      betadx=beta*al
           do ic=1,n1m
               do jc=2,n2m
      acc=betadx*a11(jc)
      apiv(jc,ic,ns)=-acc
      amiv(jc,ic,ns)=-acc
      aciv(jc,ic,ns)=1.-(apiv(jc,ic,ns)+amiv(jc,ic,ns))
               enddo
            enddo
      call ctpv2ij(1,n1m,2,n2m,ns)
      enddo
      return
      end


