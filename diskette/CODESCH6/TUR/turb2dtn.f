c  *****************************************************
c  solution of the vorticity equation in the viscous case
c
      subroutine invtrv(dvor,al,ga,ro,vor,ru,visc,dt)
      include 'param.f'
      dimension dvor(m1,m2),vor(m1,m2)
      dimension ru(m1,m2)
      dimension ami(m1),aci(m1),api(m1),fi(m2,m1),qi(m1),si(m1)
     1         ,qei(m2,m1)
      dimension amj(m2),acj(m2),apj(m2),qj(m2),sj(m2)
     1         ,qej(m1,m2)
      common/dim/n1,n1m,n2,n2m
      common/indx1/imv(m1),ipv(m1)
      common/indx2/jmv(m2),jpv(m2)
      common/mesh/dx1,dx1q,dx2,dx2q
      m2v=m2
      m1v=m1
c
c   add to the non linear term the diffusive derivatives
c   11 1nd 22
c
      do 1 jc=1,n2m
      do 1 ic=1,n1m
      ip=ipv(ic)
      im=imv(ic)
      jm=jmv(jc)
      jp=jpv(jc)
      visder=dx1q*(vor(ip,jc)+vor(im,jc))+dx2q*(vor(ic,jm)+vor(ic,jp))
     1      -2.*(dx1q+dx2q)*vor(ic,jc)
      rhs=(ga*dvor(ic,jc)+ro*ru(ic,jc))*dt+visder*dt*al*visc
      ru(ic,jc)=dvor(ic,jc)
      dvor(ic,jc)=rhs
    1 continue
c
c  ********* compute dvor  sweeping in the x1 direction
c            periodic
c
      betadx=al*dx1q*0.5*visc*dt
      do 3 ic=1,n1m
      aci(ic)=1.+betadx*2.
      ami(ic)=-betadx
      api(ic)=-betadx
    3 continue
      do 21 ic=1,n1m
      do 21 jc=1,n2m
      fi(jc,ic)=dvor(ic,jc)
   21 continue
      call tripv(ami,aci,api,fi,qi,si,qei,1,n1m,1,n2m,m2v)
      do 30 ic=1,n1m
      do 30 jc=1,n2m
      dvor(ic,jc)=fi(jc,ic)
   30 continue
c
c  ************ compute dvor sweeping along the x2 direction
c               periodic
c
      betadx=al*dx2q*0.5*visc*dt
      do 20 jc=1,n2m
      acj(jc)=1.+betadx*2.
      amj(jc)=-betadx
      apj(jc)=-betadx
   20 continue
      call tripv(amj,acj,apj,dvor,qj,sj,qej,1,n2m,1,n1m,m1v)
      return
      end
c  *****************************************************
c  solution of the vorticity equation in the inviscid  case
c  to check conservation properties of the Arakawa scheme
c
      subroutine invtri(dvor,al,ga,ro,vor,ru,dt)
      include 'param.f'
      dimension dvor(m1,m2),vor(m1,m2)
      dimension ru(m1,m2)
      common/dim/n1,n1m,n2,n2m
      common/indx1/imv(m1),ipv(m1)
      common/indx2/jmv(m2),jpv(m2)
      common/mesh/dx1,dx1q,dx2,dx2q
      m2v=m2
      m1v=m1
c
c   add to the non linear term the diffusive derivatives
c   11 1nd 22
c
      do 1 jc=1,n2m
      do 1 ic=1,n1m
      rhs=(ga*dvor(ic,jc)+ro*ru(ic,jc))*dt
      ru(ic,jc)=dvor(ic,jc)
      dvor(ic,jc)=rhs
    1 continue
      return
      end
c
c  **************  subrout tschem
c
      subroutine tschem(vor,psi,ru,al,ga,ro)
      include 'param.f'
      dimension dvor(m1,m2),vor(m1,m2),psi(m1,m2)
      common/dim/n1,n1m,n2,n2m
      dimension ru(m1,m2)
      common/visct/re
      common/tstep/dt
c
c   ****** non linear terms calculation   *********
c
      call hdnl(vor,psi,dvor)
c
c  *****  solve the dqvor=vor(n+1-vor(n) momentum equation
c
      visc=1./re
      if(re.gt..1e08) then
      call invtri(dvor,al,ga,ro,vor,ru,dt)
                      else
      call invtrv(dvor,al,ga,ro,vor,ru,visc,dt)
                      endif
      do 115 i=1,n1m
      do 115 j=1,n2m
      vor(i,j)=dvor(i,j)+vor(i,j)
  115 continue
c
c  ********* calculation of the stream functio periodic in x1
c            and tridiag in vertical x2
      call phcal(vor,psi)
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
      dimension h(m1,m2)
      dimension vor(m1,m2),psi(m1,m2)
      common/indx1/imv(m1),ipv(m1)
      common/indx2/jmv(m2),jpv(m2)
      common/in4x1/inv(m1),itv(m1)
      common/in4x2/jnv(m2),jtv(m2)
      common/chnlc/chal,chbe,chga
      common/iark4/iforth
c
c  **********  compute the non-linear terms
c
                            if(iforth.eq.0) then
c
c    second order accurate
c
      ij1=1
      ij2=0

      ad1d2=dx2*dx1*0.25
      do 10 ic=1,n1m
      ip=ipv(ic)
      im=imv(ic)
      do 10 jc=1,n2m
      jm=jmv(jc)
      jp=jpv(jc)
c
c     terms as arakawa pg 129
c
      h22a=((vor(ip,jc)-vor(im,jc))*(psi(ic,jp)-psi(ic,jm))-
     1      (vor(ic,jp)-vor(ic,jm))*(psi(ip,jc)-psi(im,jc)))
     1     *ad1d2
      h22b=(vor(ip,jc)*(psi(ip,jp)-psi(ip,jm))-
     1      vor(im,jc)*(psi(im,jp)-psi(im,jm))-
     1      vor(ic,jp)*(psi(ip,jp)-psi(im,jp))+
     1      vor(ic,jm)*(psi(ip,jm)-psi(im,jm)))
     1     *ad1d2
      h22c=(vor(ip,jp)*(psi(ic,jp)-psi(ip,jc))-
     1      vor(im,jm)*(psi(im,jc)-psi(ic,jm))-
     1      vor(im,jp)*(psi(ic,jp)-psi(im,jc))+
     1      vor(ip,jm)*(psi(ip,jc)-psi(ic,jm)))
     1     *ad1d2
      hq2=chal*h22a+chbe*h22b+chga*h22c
      h(ic,jc)=hq2
   10 continue
                                     endif
      hd1d2=dx2*dx1*0.125
                     if(iforth.eq.1) then
c
c    fourth order accurate
c
      ij1=1
      ij2=1
      do 11 ic=1,n1m
      ip=ipv(ic)
      im=imv(ic)
      it=itv(ic)
      in=inv(ic)
      do 11 jc=1,n2m
      jm=jmv(jc)
      jp=jpv(jc)
      jt=jtv(jc)
      jn=jnv(jc)
c
c     terms as arakawa 
c
      h22a=((vor(ip,jc)-vor(im,jc))*(psi(ic,jp)-psi(ic,jm))-
     1      (vor(ic,jp)-vor(ic,jm))*(psi(ip,jc)-psi(im,jc)))
     1     *ad1d2
      h22b=(vor(ip,jc)*(psi(ip,jp)-psi(ip,jm))-
     1      vor(im,jc)*(psi(im,jp)-psi(im,jm))-
     1      vor(ic,jp)*(psi(ip,jp)-psi(im,jp))+
     1      vor(ic,jm)*(psi(ip,jm)-psi(im,jm)))
     1     *ad1d2
      h22c=(vor(ip,jp)*(psi(ic,jp)-psi(ip,jc))-
     1      vor(im,jm)*(psi(im,jc)-psi(ic,jm))-
     1      vor(im,jp)*(psi(ic,jp)-psi(im,jc))+
     1      vor(ip,jm)*(psi(ip,jc)-psi(ic,jm)))
     1     *ad1d2
      hq2=chal*h22a+chbe*h22b+chga*h22c
      h22d=((vor(ip,jp)-vor(im,jm))*(psi(im,jp)-psi(ip,jm))-
     1      (vor(im,jp)-vor(ip,jm))*(psi(ip,jp)-psi(im,jm)))
     1     *hd1d2
      h22e=(vor(ip,jp)*(psi(ic,jt)-psi(it,jc))-
     1      vor(im,jm)*(psi(in,jc)-psi(ic,jn))-
     1      vor(im,jp)*(psi(ic,jt)-psi(in,jc))+
     1      vor(ip,jm)*(psi(it,jc)-psi(ic,jn)))
     1     *hd1d2
      h22f=(vor(it,jc)*(psi(ip,jp)-psi(ip,jm))-
     1      vor(in,jc)*(psi(im,jp)-psi(im,jm))-
     1      vor(ic,jt)*(psi(ip,jp)-psi(im,jp))+
     1      vor(ic,jn)*(psi(ip,jm)-psi(im,jm)))
     1     *hd1d2
      hq4=chal*h22d+chbe*h22e+chga*h22f
      h(ic,jc)=2.*hq2-hq4
   11 continue 
                                     endif
                     if(iforth.eq.2) then
c
c    fourth order accurate
c
      ij1=0
      ij2=1
      do 12 ic=1,n1m
      ip=ipv(ic)
      im=imv(ic)
      it=itv(ic)
      in=inv(ic)
      do 12 jc=1,n2m
      jm=jmv(jc)
      jp=jpv(jc)
      jt=jtv(jc)
      jn=jnv(jc)
      h22d=((vor(ip,jp)-vor(im,jm))*(psi(im,jp)-psi(ip,jm))-
     1      (vor(im,jp)-vor(ip,jm))*(psi(ip,jp)-psi(im,jm)))
     1     *hd1d2
      h22e=(vor(ip,jp)*(psi(ic,jt)-psi(it,jc))-
     1      vor(im,jm)*(psi(in,jc)-psi(ic,jn))-
     1      vor(im,jp)*(psi(ic,jt)-psi(in,jc))+
     1      vor(ip,jm)*(psi(it,jc)-psi(ic,jn)))
     1     *hd1d2
      h22f=(vor(it,jc)*(psi(ip,jp)-psi(ip,jm))-
     1      vor(in,jc)*(psi(im,jp)-psi(im,jm))-
     1      vor(ic,jt)*(psi(ip,jp)-psi(im,jp))+
     1      vor(ic,jn)*(psi(ip,jm)-psi(im,jm)))
     1     *hd1d2
      hq4=chal*h22d+chbe*h22e+chga*h22f
      h(ic,jc)=hq4
   12 continue 
                  endif
      do 13 ic=1,n1m
      do 13 jc=1,n2m
      h(ic,jc)=-h(ic,jc)
   13 continue
      return
      end
c
c  ****************************** subrout tripv  **********************
c
      subroutine tripv ( a,b,c,f,q,s,qe,j1,j2,m1,m2,mv )
c
      dimension a(1),b(1),c(1),f(mv,1),q(1),s(1),qe(mv,1),fn(1025)
      ja = j1 + 1
      jj = j1 + j2
      q(j1) = -c(j1)/b(j1)
      s(j1) = - a(j1)/b(j1)
      do 20 k=m1,m2
      fn(k) = f(k,j2)
      f(k,j1) = f(k,j1)/b(j1)
   20 continue
c
c     forward elimination sweep
c
      do 10 j=ja,j2
      p =1./( b(j) + a(j)*q(j-1))
      q(j) = - c(j)*p
      s(j) = - a(j)*s(j-1)*p
      do 21 k=m1,m2
      f(k,j) = ( f(k,j) - a(j)*f(k,j-1))*p
   21 continue
   10 continue
c
c     backward pass
c
      s(j2) = 1.
      do 22 k=m1,m2
      qe(k,j2) = 0.
   22 continue
      do 11 i=ja,j2
      j = jj - i
      s(j) = s(j) + q(j)*s(j+1)
      do 23 k=m1,m2
      qe(k,j) = f(k,j) + q(j)*qe(k,j+1)
   23 continue
   11 continue
      do 24 k=m1,m2
      f(k,j2)=(fn(k) - c(j2)*qe(k,j1) - a(j2)*qe(k,j2-1))
     &       /(c(j2)*s(j1) + a(j2)*s(j2-1)  +b(j2))
   24 continue
c
c     backward elimination pass
c
      do 12 i=ja,j2
      j = jj -i
      do 25 k=m1,m2
      f(k,j) = f(k,j2)*s(j) + qe(k,j)
   25 continue
   12 continue
      return
      end
