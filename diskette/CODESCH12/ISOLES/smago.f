c***********************************************
c     Quantities for LES s.g.s
c     here the Smagorinsky eddy viscosity is calculated
c***********************************************
      subroutine quales 
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      common/d123/alx1,alx2,alx3
      common/qles1/deltax1,deltay1,deltaz1,ell1c
      common/qles2/deltax2,deltay2,deltaz2,ell2c
c      
c     calculate two length scale
c
      deltax1=alx1/float(n1m)
      deltay1=alx2/float(n2m)
      deltaz1=alx3/float(n3m)
       ell1c=(deltax1*deltaz1*deltay1)**.66667
c      
c     calculate tsecond length scale for dynamic model
c
      deltax2=2.*deltax1
      deltay2=2.*deltay1
      deltaz2=2.*deltaz1
       ell2c=(deltax2*deltaz2*deltay2)**.66667
      return
      end
c***********************************************
c    Before the Smagorinsky subgrid model is calculated
c    Then it is averaged to get a constant viscosity
c    This model works well at low R-lambda
c     Constant Smagorinsky 
c***********************************************
      subroutine consma(q,stst,ncount,time)
      include 'param.f'
      dimension q(ndv,m1,m2,m3),stst(m1,m2,m3)
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      common/d123/alx1,alx2,alx3
      common/sma/vis(m1,m2,m3)
      common/strain/st(m1,m2,m3,6)
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/ledat/ics0,cvisc,ifiltr,ibox
      common/d2/nstop,nprint,ntst,npin,npstf
      common/cosma/csma,pratu
      common/qles1/deltax1,deltay1,deltaz1,ell1c
      common/qles2/deltax2,deltay2,deltaz2,ell2c
      common/smapra/gold(m1,m2,m3)
      common/strat2/igrad,bvais
      common/strat1/istrat,rho0,g,schm
c
c      computation of Sij
c
      call strper(q)
c                                                                       
c     this calculates the  strain rate tensor at the center of
c     the cell. Periodic box. The strains are defined as
c     Sij=(dui/dxj + duj/dxi)*0.5
c
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
      do kc=1,n3m 
      do jc=1,n2m                                                     
      do ic=1,n1m                                                     
      galt=stst(ic,jc,kc)
        gtmp=sqrt(galt) 
c     calculate smagorinsky constant
      vis(ic,jc,kc)=((csma*csma)*ell1c)
     &               *gtmp+cvisc
      enddo
      enddo
      enddo
c      
      vismed=0.
      do kc=1,n3m 
      do jc=1,n2m                                                     
      do ic=1,n1m                                                     
      vismed=vismed+vis(ic,jc,kc)
      enddo
      enddo
      enddo
      vl123=1./float(n1m*n2m*n3m)
      do kc=1,n3m 
      do jc=1,n2m                                                     
      do ic=1,n1m                                                     
      vis(ic,jc,kc)=vismed*vl123
      enddo
      enddo
      enddo
c
c   here the eddy diffusivity is calculated when
c   The turbulent Prandtl number is assumed
c
      if(istrat.eq.1.and.igrad.eq.0) then
      do kc=1,n3m 
      do jc=1,n2m                                                     
      do ic=1,n1m                                                     
      gold(ic,jc,kc)=vis(ic,jc,kc)/pratu
      enddo
      enddo
      enddo
                                     endif
       return
       end 
c***********************************************
c     Smagorinsky subgrid-scale stress model 
c***********************************************
      subroutine smago(q,stst,ncount,time)
      include 'param.f'
      dimension q(ndv,m1,m2,m3),stst(m1,m2,m3)
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      common/d123/alx1,alx2,alx3
      common/sma/vis(m1,m2,m3)
      common/strain/st(m1,m2,m3,6)
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/ledat/ics0,cvisc,ifiltr,ibox
      common/d2/nstop,nprint,ntst,npin,npstf
      common/cosma/csma,pratu
      common/qles1/deltax1,deltay1,deltaz1,ell1c
      common/qles2/deltax2,deltay2,deltaz2,ell2c
      common/smapra/gold(m1,m2,m3)
      common/strat2/igrad,bvais
      common/strat1/istrat,rho0,g,schm
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
c
c   here the eddy diffusivity is calculated
c   when the turbulent Prandtl number is assigned
c
      if(istrat.eq.1.and.igrad.eq.0) then
      do kc=1,n3m 
      do jc=1,n2m                                                     
      do ic=1,n1m                                                     
            diftu=(vis(ic,jc,kc)-cvisc)/pratu
            gold(ic,jc,kc)=diftu+cvisc/schm
      enddo
      enddo
      enddo
                                     endif
c      
       return
       end 
