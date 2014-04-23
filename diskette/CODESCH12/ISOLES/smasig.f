c***********************************************
c     Smagorinsky subgrid-scale stress model 
c***********************************************
      subroutine smarho(q,stst,rho,ncount,time)
      include 'param.f'
      dimension q(ndv,m1,m2,m3),stst(m1,m2,m3)
      dimension rho(m1,m2,m3)
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      common/d123/alx1,alx2,alx3
      common/sma/vis(m1,m2,m3)
      common/smapra/gold(m1,m2,m3)
      common/strain/st(m1,m2,m3,6)
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/ledat/ics0,cvisc,ifiltr,ibox
      common/strat2/igrad,bvais
      common/strat1/istrat,rho0,g,schm
      common/d2/nstop,nprint,ntst,npin,npstf
      common/cosma/csma,pratu
      common/qles1/deltax1,deltay1,deltaz1,ell1c
      common/qles2/deltax2,deltay2,deltaz2,ell2c
c
c      computation of Sij
c
      call strper(q)
c                                                                       
c     this calculates the  strain rate tensor at the center of
c     the cell. Periodic box. The strains are defined as
c     Sij=(dui/dxj + duj/dxi)*0.5

c To evaluate the eddy viscosity now we use the expression
c usually used by the Stanford group see  Piomelli
c This is done to be consistent with the expressions in
c the dynamic model
c 
c     
      do 4 k=1,n3m 
      do 4 j=1,n2m                                                     
      do 4 i=1,n1m                                                     
c
c calculate 2*SijSij at center of cell
c
      stst(i,j,k)=2.*(st(i,j,k,1)**2+st(i,j,k,2)**2+st(i,j,k,3)**2)+
     1           4.*(st(i,j,k,4)**2+st(i,j,k,5)**2+st(i,j,k,6)**2)
    4 continue                                                          
      if(igrad.eq.0) then
      do kc=1,n3m 
      do jc=1,n2m                                                     
      do ic=1,n1m                                                     
      galt=stst(ic,jc,kc)
        gtmp=sqrt(galt) 
c     calculate smagorinsky eddy viscosity
      vis(ic,jc,kc)=((csma*csma)*ell1c)
     &               *gtmp+cvisc
      enddo
      enddo
      enddo
                     else
c      
       visav=0.
       difav=0.
       stcon=0.
       pecon=0.
       npne=0
      do kc=1,n3m 
      do jc=1,n2m                                                     
      do ic=1,n1m                                                     
c
c calculate drho/dz at center of cell.
c recall that we define rho on lower face of cell
c
      drdz=(rho(ic,jc,kpv(kc))-rho(ic,jc,kc))*dx3
      galt=stst(ic,jc,kc)-drdz/pratu
       if(galt.gt.0.) then
        gtmp=sqrt(galt)
      stcon=stcon+stst(ic,jc,kc)
      pecon=pecon-drdz/pratu
         else
        gtmp=0.0
        npne=npne+1
       endif
c     calculate smagorinsky eddy viscosity
      vistu=((csma*csma)*ell1c)*gtmp
           visav=visav+vistu
            vis(ic,jc,kc)=vistu+cvisc
      enddo
      enddo
      enddo
                     endif
      do kc=1,n3m 
      do jc=1,n2m                                                     
      do ic=1,n1m                                                     
            diftu=(vis(ic,jc,kc)-cvisc)/pratu
           difav=difav+diftu
            gold(ic,jc,kc)=diftu+cvisc/schm
      enddo
      enddo
      enddo
       stcon=stcon/float(n1m*n2m*n3m)
       pecon=pecon/float(n1m*n2m*n3m)
       visav=visav/float(n1m*n2m*n3m)
       difav=difav/float(n1m*n2m*n3m)
       write(61,161)time,stcon,pecon,npne,visav,difav
  161 format(2x,'visc smarho',3e12.4,3x,i8,2x,2e12.4)
       return
       end 
c***********************************************
c     Smagorinsky subgrid-scale stress model 
c***********************************************
      subroutine cosmrh(q,stst,rho,ncount,time)
      include 'param.f'
      dimension q(ndv,m1,m2,m3),stst(m1,m2,m3)
      dimension rho(m1,m2,m3)
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      common/d123/alx1,alx2,alx3
      common/sma/vis(m1,m2,m3)
      common/smapra/gold(m1,m2,m3)
      common/strain/st(m1,m2,m3,6)
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/ledat/ics0,cvisc,ifiltr,ibox
      common/strat2/igrad,bvais
      common/strat1/istrat,rho0,g,schm
      common/d2/nstop,nprint,ntst,npin,npstf
      common/cosma/csma,pratu
      common/qles1/deltax1,deltay1,deltaz1,ell1c
      common/qles2/deltax2,deltay2,deltaz2,ell2c
c
c      computation of Sij
c
      call strper(q)
c                                                                       
c     this calculates the  strain rate tensor at the center of
c     the cell. Periodic box. The strains are defined as
c     Sij=(dui/dxj + duj/dxi)*0.5

c To evaluate the eddy viscosity now we use the expression
c usually used by the Stanford group see  Piomelli
c This is done to be consistent with the expressions in
c the dynamic model
c 
c     
      do 4 k=1,n3m 
      do 4 j=1,n2m                                                     
      do 4 i=1,n1m                                                     
c
c calculate 2*SijSij at center of cell
c
      stst(i,j,k)=2.*(st(i,j,k,1)**2+st(i,j,k,2)**2+st(i,j,k,3)**2)+
     1           4.*(st(i,j,k,4)**2+st(i,j,k,5)**2+st(i,j,k,6)**2)
    4 continue                                                          
      if(igrad.eq.0) then
      do kc=1,n3m 
      do jc=1,n2m                                                     
      do ic=1,n1m                                                     
      galt=stst(ic,jc,kc)
        gtmp=sqrt(galt) 
c     calculate smagorinsky eddy viscosity
      vis(ic,jc,kc)=((csma*csma)*ell1c)
     &               *gtmp+cvisc
      enddo
      enddo
      enddo
                     else
c      
       visav=0.
       difav=0.
       stcon=0.
       pecon=0.
       npne=0
      do kc=1,n3m 
      do jc=1,n2m                                                     
      do ic=1,n1m                                                     
c calculate drho/dz at center of cell.
c recall that we define rho on lower face of cell
      drdz=(rho(ic,jc,kpv(kc))-rho(ic,jc,kc))*dx3
      galt=stst(ic,jc,kc)-drdz/pratu
       if(galt.gt.0.) then
        gtmp=sqrt(galt)
      stcon=stcon+stst(ic,jc,kc)
      pecon=pecon-drdz/pratu
         else
        gtmp=0.0
        npne=npne+1
       endif
c     calculate smagorinsky eddy viscosity
      vistu=((csma*csma)*ell1c)*gtmp
           visav=visav+vistu
      enddo
      enddo
      enddo
                     endif
       visav=visav/float(n1m*n2m*n3m)
      do kc=1,n3m
      do jc=1,n2m
      do ic=1,n1m
      vis(ic,jc,kc)=visav+cvisc  
      enddo
      enddo
      enddo
      do kc=1,n3m 
      do jc=1,n2m                                                     
      do ic=1,n1m                                                     
            diftu=(vis(ic,jc,kc)-cvisc)/pratu
           difav=difav+diftu
      enddo
      enddo
      enddo
       difav=difav/float(n1m*n2m*n3m)
      do kc=1,n3m
      do jc=1,n2m
      do ic=1,n1m
      gold(ic,jc,kc)=difav+cvisc/schm
      enddo
      enddo
      enddo
       return
       end 
