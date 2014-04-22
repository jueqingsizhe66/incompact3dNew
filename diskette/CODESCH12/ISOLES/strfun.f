c***********************************************
c     Structure function subgrid-scale stress model 
c     by Metais and Lesieur
c     To be used whithout density
c    
c***********************************************
      subroutine strfun(q,g,ncount,time)
      include 'param.f'
      dimension q(ndv,m1,m2,m3)
      dimension g(m1,m2,m3)
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      common/d123/alx1,alx2,alx3
      common/sma/vis(m1,m2,m3)
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/ledat/ics0,cvisc,ifiltr,ibox
      common/cosma/csma
      common/instf/icsf(6),jcsf(6),kcsf(6)
      common/icstf/ics(0:m1),jcs(0:m2),kcs(0:m3)
c
c      computation of Longitudinal structure function
c
      deltax1=alx1/float(n1m)
      deltay1=alx2/float(n2m)
      deltaz1=alx3/float(n3m)
       delta=(deltax1*deltaz1*deltay1)**.33333
      ells=0.1*1.4**(-3./2.)*delta
          do kcc=1,n3m
              do jcc=1,n2m
                  do icc=1,n1m
      strfto=0.
      do l=1,6
      ic=ics(icc+icsf(l))
      jc=jcs(jcc+jcsf(l))
      kc=kcs(kcc+kcsf(l))
      jp=jpv(jc)
      kp=kpv(kc)
      ip=ipv(ic)
      stfu1=(q(1,ip,jc,kc)-q(1,ic,jc,kc))**2
      stfu2=(q(2,ic,jp,kc)-q(2,ic,jc,kc))**2
      stfu3=(q(3,ic,jc,kp)-q(3,ic,jc,kc))**2
      strfto=strfto+stfu1+stfu2+stfu3
      enddo
      g(icc,jcc,kcc)=sqrt(strfto/6.)
                  enddo
              enddo
          enddo
c
c       calculate Structure function eddy viscosity 
c     visma=-1000.
c     vismi=+1000.
       do k=1,n3m
       do j=1,n2m
       do i=1,n1m
       vistu=ells*g(i,j,k)
       vis(i,j,k)=vistu+cvisc
                  enddo
              enddo
          enddo
c     write(6,*)' struct fun ',visma,vismi,cvisc,ells
       return
       end 
c
c
      subroutine indstf
      include 'param.f'
      common/instf/icsf(6),jcsf(6),kcsf(6)
      common/icstf/ics(0:m1),jcs(0:m2),kcs(0:m3)
      common/dim/n1,n1m,n2,n2m,n3,n3m
      icsf(1)=-1
      jcsf(1)=0
      kcsf(1)=0
      icsf(2)=+1
      jcsf(2)=0
      kcsf(2)=0
      jcsf(3)=-1
      icsf(3)=0
      kcsf(3)=0
      jcsf(4)=+1
      icsf(4)=0
      kcsf(4)=0
      kcsf(5)=-1
      icsf(5)=0
      jcsf(5)=0
      kcsf(6)=+1
      icsf(6)=0
      jcsf(6)=0
      do ic=0,n1
      ics(ic)=ic
      if(ic.eq.0) then
      ics(ic)=n1m
                   endif
      if(ic.eq.n1) then
      ics(ic)=1
                   endif
      enddo
      do jc=0,n2
      jcs(jc)=jc
      if(jc.eq.0) then
      jcs(jc)=n2m
                   endif
      if(jc.eq.n2) then
      jcs(jc)=1
                   endif
      enddo
      do kc=0,n3
      kcs(kc)=kc
      if(kc.eq.0) then
      kcs(kc)=n3m
                   endif
      if(kc.eq.n3) then
      kcs(kc)=1
                   endif
      enddo

      return
      end
