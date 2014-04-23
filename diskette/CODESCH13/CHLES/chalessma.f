c***********************************************
c     Constant Smagorinsky 
c***********************************************
      subroutine consma(q,stst)
      include 'param.f'
      dimension q(ndv,m1,m2,m3),stst(m1,m2,m3)
      common/strain/st(m1,m2,m3,6)
      common/qdyn/g(m1,m2,m3)
      common/smapra/gold(m1,m2,m3)
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/metria/caj(m2),cac(m2)
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/ledat/csma,cvisc,iles
      common/outd1/cdyn(m2),alij(m2),amij(m2)
      common/outd2/cs(m2),lijmij(m2),mijmij(m2)
      common/outd3/accdyn(m2),acclij(m2),accmij(m2)
      common/qles1/deltax1,deltay1,deltaz1,ell1c(m2)
      common/qles2/deltax2,deltay2,deltaz2,ell2c(m2)
      common/sma/vis(m1,m2,m3)
      dimension vismed(m2)
      common/y2sta/y2s(m2)
      common/corpri/yp1(m1),yp2(m2),yp3(m3)
c
c     calculate smagorinsky eddy viscosity
c
      call smago(q,stst)
c      
c   average in x1-x3 planes
c
      do jc=1,n2m                                                     
      vismed(jc)=0.
      do kc=1,n3m 
      do ic=1,n1m                                                     
      vismed(jc)=vismed(jc)+vis(ic,jc,kc)
      enddo
      enddo
      enddo
      vl13=1./float(n1m*n3m)
c
c  set eddy viscosity equal to the averaged
c
      do jc=1,n2m                                                     
      do kc=1,n3m 
      do ic=1,n1m                                                     
      vis(ic,jc,kc)=vismed(jc)*vl13
      enddo
      enddo
      enddo
       return
       end 
c***********************************************
c     Smagorinsky subgrid-scale stress model 
c***********************************************
      subroutine smago(q,stst)
      include 'param.f'
      dimension q(ndv,m1,m2,m3),stst(m1,m2,m3)
      common/strain/st(m1,m2,m3,6)
      common/qdyn/g(m1,m2,m3)
      common/smapra/gold(m1,m2,m3)
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/metria/caj(m2),cac(m2)
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/ledat/csma,cvisc,iles
      common/outd1/cdyn(m2),alij(m2),amij(m2)
      common/outd2/cs(m2),lijmij(m2),mijmij(m2)
      common/outd3/accdyn(m2),acclij(m2),accmij(m2)
      common/qles1/deltax1,deltay1,deltaz1,ell1c(m2)
      common/qles2/deltax2,deltay2,deltaz2,ell2c(m2)
      common/sma/vis(m1,m2,m3)
      dimension vismed(m2)
      common/y2sta/y2s(m2)
      common/corpri/yp1(m1),yp2(m2),yp3(m3)
c
c      computation of Sij
c
      call strai(q)
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
c
c   van-driest damping
c
      dyl=(y2s(1)-yp2(1))
      dyu=-(y2s(n2m)-yp2(n2))
      dql=0.
      dqu=0.
      do kc=1,n3m
      do ic=1,n1m
      dql=q(3,ic,1,kc)/dyl+dql
      dqu=-q(3,ic,n2m,kc)/dyu+dqu
      enddo
      enddo
      cfllk=abs(dql*cvisc/(n1m*n2m))
      cfulk=abs(dqu*cvisc/(n1m*n2m))
      cfavp=(cfllk+cfulk)*0.5
      utap=sqrt(abs(cfavp))
      yplf=utap/cvisc
      apl=0.1*yplf 
      do jc=1,n2m                                                     
      if(jc.le.n2m/2) then
      yd=y2s(jc)-yp2(1)
      ypl=yd*utap/cvisc
                     else
      yd=-y2s(jc)+yp2(n2)
      ypl=yd*utap/cvisc
                     endif
      yex=exp(-ypl/apl)
      damp=(1.-yex**2)**2

      do kc=1,n3m 
      do ic=1,n1m                                                     
      galt=stst(ic,jc,kc)
        gtmp=sqrt(galt) 
c
c     calculate smagorinsky constant
c
      vis(ic,jc,kc)=((csma*csma)*damp*ell1c(jc))
     &               *gtmp+cvisc
      enddo
      enddo
      enddo
       return
       end 
