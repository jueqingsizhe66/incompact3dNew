c
c  ****************************** subrout initia **********************
c
c   initial  conditions in the whole field.
c   pr(i,j,k)=0 
c   set up  dp3ns. for Poiseuille
c   evaluates the transpired velocity
c
      subroutine initia(q,y,pr)
c
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      dimension q(ndv,m1,m2,m3),pr(m1,m2,m3)
      dimension y(m2)
      common/tstep/dt,beta,ren
      common/rhs3p/dp3ns
c
      call inqpr(q,y)
      call slotin
c
      do 4 k=1,n3
      do 4 j=1,n2
      do 4 i=1,n1
      pr(i,j,k)=0.
    4 continue
      dp3ns=-2./ren
      return
      end

c  ****************************** subrout inqpr **********************
c
c   Initial laminar Poiseuille
c   with superimposed a random disturbance
c
      subroutine inqpr(q,y)
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/metria/caj(m2),cac(m2)
      dimension q(ndv,m1,m2,m3)
      dimension y(m2)
      common/tstep/dt,beta,re
      common/vperin/vper
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/vini0/v30(m2)
      common/velmax/vmax(ndv),vmaxo(ndv)
      dimension v1m(m2),v2m(m2),v3m(m2)
      common/y2sta/y2s(m2)
      dimension q1m(m2),q2m(m2)
c
      pi=2.*asin(1.)
      vol=n1m*n2m*n3m*2.
      sv3ii=0.
      sv1iip=0.
      sv2iip=0.
      sv3iip=0.
       isd=17
       call srand(isd)
      do 919 j=1,n2m
      v1mj=0.
      v2mj=0.
      v3mj=0.
      v30(j)=(1.-y2s(j)**2)
      ydw=1.-abs(y2s(j))
      vperq=vper
ccccccccccccccccccccccccccccccccccccc
      if(ydw.lt.0.025) vperq=vper/5.
ccccccccccccccccccccccccccccccccccccc
      do 917 k=1,n3m
      do 917 i=1,n1m
c     q(1,i,j,k)=vperq*rnd()
c     q(2,i,j,k)=vperq*rnd()
c     q(3,i,j,k)=vperq*rnd()
      q(1,i,j,k)=vperq*rand()
      q(2,i,j,k)=vperq*rand()
      q(3,i,j,k)=vperq*rand()
      v1mj=v1mj+q(1,i,j,k)
      v2mj=v2mj+q(2,i,j,k)
      v3mj=v3mj+q(3,i,j,k)
  917 continue
      v1m(j)=v1mj/(n1m*n3m)
      v2m(j)=v2mj/(n1m*n3m)
      v3m(j)=v3mj/(n1m*n3m)
  919 continue
      do 915 k=1,n3m
      do 915 i=1,n1m
      q(2,i,1,k)=0.
      q(2,i,n2,k)=0.
  915 continue
      do 916 k=1,n3m
      do 916 j=1,n2m
      do 916 i=1,n1m
      sv3ii=sv3ii+v30(j)*caj(j)/vol
      q(1,i,j,k)=q(1,i,j,k)-v1m(j)
      sv1iip=sv1iip+q(1,i,j,k)*caj(j)/vol
      q(2,i,j,k)=q(2,i,j,k)-v2m(j)
      sv2iip=sv2iip+(q(2,i,j,k)*cac(j)
     1              +q(2,i,j+1,k)*cac(j+1))*0.5/vol
      q(3,i,j,k)=q(3,i,j,k)-v3m(j)
      sv3iip=sv3iip+q(3,i,j,k)*caj(j)/vol
      q(3,i,j,k)=(v30(j)+q(3,i,j,k))
  916 continue
      do 920 k=1,n3m
      do 920 i=1,n1m
      q(2,i,1,k)=0.
      q(2,i,n2,k)=0.
  920 continue
      do 44 k=1,n3m
      do 44 j=1,n2m
      q(1,n1,j,k)=q(1,1,j,k)
   44 continue
      do 45 j=1,n2m
      do 45 i=1,n1m
      q(3,i,j,n3)=q(3,i,j,1)
   45 continue
      do 410 j=1,n2m
      vm1m=0.
      vm2m=0.
      do 411 i=1,n1m
      ip=ipv(i)
      do 411 k=1,n3m
      kp=kpv(k)
      vm1m=q(1,i,j,k)+vm1m
      vm2m=q(2,i,j,k)+vm2m
  411 continue
      q1m(j)=vm1m/(n1m*n3m)
      q2m(j)=vm2m/(n1m*n3m)
  410 continue
      do 414 k=1,n3m
      do 414 j=1,n2m
      do 414 i=1,n1m
      q(1,i,j,k)=q(1,i,j,k)-q1m(j)
  414 continue
      do 415 k=1,n3m
      do 415 j=2,n2m
      do 415 i=1,n1m
      q(2,i,j,k)=q(2,i,j,k)-q2m(j)
  415 continue
c
      call vmaxv(q)
c
      write(6,115) sv3ii,sv3iip,sv2iip,sv1iip
      write(61,115) sv3ii,sv3iip,sv2iip,sv1iip
  115 format(3x,' v3i=',e10.3,2x,'v3ip=',e10.3,2x,
     1      'v2ip=',e10.3,2x,'v1ip=',e10.3)
      write(6,116) (vmax(l),l=1,3)
      write(61,116) (vmax(l),l=1,3)
  116 format(3x,'vm(1)=',e11.4,2x,'vm(2)=',e11.4,2x,'vm(3)=',e11.4)
      do 120 j=1,n2m
c     write(6,615)j,q1m(j),q2m(j),v3m(j)
  615 format(3x,'init. prof v3,q1,q2',i3,2x,3e12.4)
  120 continue
      return
      end
c  ****************************** subrout slotin **********************
c
c   In this routine the location for the transpiration is 
c   evaluated
c
      subroutine slotin
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/vslot/v2inf(m1,m3)
      common/y13sta/y1s(m1),y3s(m3)
      common/slotfl/flowq2,tau2
      common/slotpa/y1gsl,y1ssl,y3gsl,y3ssl
      common/slotdi/y1disl,y3disl
      common/d13/alx1,alx3
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      alx3i=0.
      alx3f=alx3i+alx3
      alx1i=0.
      alx1f=alx1i+alx1
      open(68,file='vinf.out')
      d1slo=y1gsl-y1ssl
      d3slo=y3gsl-y3ssl
      p1slo=(y1gsl+y1ssl)*0.5
      p3slo=(y3gsl+y3ssl)*0.5
      riin=y1ssl/alx1*n1m
      rifi=y1gsl/alx1*n1m
      rkin=(y3ssl-alx3i)/(alx3f-alx3i)*n3m
      rkfi=(y3gsl-alx3i)/(alx3f-alx3i)*n3m
      iin=riin+1
      ifi=rifi+1
      kin=rkin+1
      kfi=rkfi+1
       write(68,*)'iin,ifi,kin,kfi  ',iin,ifi,kin,kfi
       write(6,*)'first slot iin,ifi,kin,kfi  ',iin,ifi,kin,kfi
      v2tos1=0.
      do kc=1,n3m
      do ic=1,n1m
      v2inf(ic,kc)=0.
      if(y1s(ic).le.y1gsl.and.y1s(ic).ge.y1ssl) then
      if(y3s(kc).le.y3gsl.and.y3s(kc).ge.y3ssl) then
      if(y1s(ic).le.p1slo) then
      y1d=(y1s(ic)-y1ssl)/(d1slo*0.125)
      if(y1ssl.eq.0.) y1d=1.
                           else
      y1d=-(y1s(ic)-y1gsl)/(d1slo*0.125)
      if(y1gsl.eq.alx1) y1d=1.
                           endif
      if(y3s(kc).le.p3slo) then
      y3d=(y3s(kc)-y3ssl)/(d1slo*0.125)
                           else
      y3d=-(y3s(kc)-y3gsl)/(d1slo*0.125)
                           endif
      v2inf(ic,kc)=-flowq2*ftir(y1d)*ftir(y3d)
                                                endif
                                                endif
      v2tos1=v2tos1+v2inf(ic,kc)/(dx1*dx3)
      enddo
      enddo
      rkdi=(y3disl-alx3i)/(alx3f-alx3i)*n3m
      ridi=(y1disl-alx1i)/(alx1f-alx1i)*n3m
      idisl=ridi
      kdisl=rkdi
       write(6,*)'secon slot idisl,kdisl ',idisl,kdisl        
      is2i=idisl+iin
      is2f=idisl+ifi
      ks2i=kdisl+kin
      ks2f=kdisl+kfi
       write(6,*)'secon slot iin,ifi,kin,kfi  ',is2i,is2f,ks2i,ks2f
      do is2=is2i,is2f
      do ks2=ks2i,ks2f
      is1=is2-idisl
      ks1=ks2-kdisl
      v2inf(is2,ks2)=-v2inf(is1,ks1)  
      enddo
      enddo
      write(68,*)' first slot'
      do kc=kin,kfi
      write(68,133) (v2inf(ic,kc),ic=iin,ifi)
      enddo
      write(68,*)' second slot'
      do kc=ks2i,ks2f
      write(68,133) (v2inf(ic,kc),ic=is2i,is2f)
      enddo
  133 format(11e11.4)
      totv2=0.
      do ic=1,n1m
      do kc=1,n3m
      totv2=totv2+v2inf(ic,kc)
      enddo
      enddo
      write(6,*)'from first slot',v2tos1,'  total v=',totv2
      write(68,*)'from first slot',v2tos1,'  total v=',totv2
      close(68)
      return
      end
c  
c  **************  function ftir
c   a smoothing function
c
      function  ftir(tina)
      if(tina.ge.1.) then
      ftir=1.
                     else
      ftir=3.*tina**2-2.*tina**3
c     ftir=tina
                     endif
      return
      end
