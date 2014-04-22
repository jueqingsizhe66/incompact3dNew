c
c
c  ****************************** subrout outrho **********************
c
c   calculation of total kinetic energy
c           enkin=sum rho0*(v1**2+v2**2+v3**2)*0.5 dxdydz
c   calculation of total potential energy
c           enpot=sum -g/(drho/dz)(rho**2)*0.5 dxdydz
c   calculation of total energy = entot = enkin + enpot
c   vi cartesian components
c
      subroutine outrho(ntime,time,q,rho)
c
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      dimension q(ndv,m1,m2,m3),rho(m1,m2,m3)
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      common/meshu/d1x1,d1x2,d1x3
      common/d123/alx1,alx2,alx3
      common/vmean/vm(ndv),vrms(6)
      common/strat1/istrat,rho0,g,schm
      common/strat2/igrad,bvais
      common/strat3/grbar(m3)
c
      vm1m=0.
      vm2m=0.
      vm3m=0.
      vl123=1./float(n1m*n2m*n3m)
      dx=1./dx1
      dy=1./dx2
      dz=1./dx3
      do 411 k=1,n3m
      do 411 j=1,n2m
      do 411 i=1,n1m
       vm1m=q(1,i,j,k)+vm1m
       vm2m=q(2,i,j,k)+vm2m
       vm3m=q(3,i,j,k)+vm3m
  411 continue
      vm(1)=vm1m*vl123
      vm(2)=vm2m*vl123
      vm(3)=vm3m*vl123
  410 continue
      enej1=0.
      enej2=0.
      enej3=0.
      do 414 k=1,n3m
      do 414 j=1,n2m
      do 414 i=1,n1m
       v1m=q(1,i,j,k)-vm(1)
       v2m=q(2,i,j,k)-vm(2)
       v3m=q(3,i,j,k)-vm(3)
      enej1=enej1+v1m**2
      enej2=enej2+v2m**2
      enej3=enej3+v3m**2
  414 continue
      enki1=rho0*enej1*vl123
      enki2=rho0*enej2*vl123
      enki3=rho0*enej3*vl123
      enkin=(enki1+enki2+enki3)
      entot=enkin
      sum=0.0
      if(igrad.eq.0) then
      do k=1,n3m
      do j=1,n2m
      do i=1,n1m
      sum=sum+rho(i,j,k)**2
      enddo
      enddo
      enddo
                     else
      do k=1,n3m
      scal=1.
      if(bvais.ne.0.) scal=1./grbar(k)
      do j=1,n2m
      do i=1,n1m
      sum=sum+rho(i,j,k)**2*scal
      enddo
      enddo
      enddo
                     endif
      open(43,file='totene.out')
      enpot=-0.5*g*sum*vl123
      entot=enkin+enpot
      write(43,*)time,enki1,enki2,enki3,enpot,entot
      close(43)
      rhom=0.
      vl123=1./float(n1m*n2m*n3m)
      dx=1./dx1
      dy=1./dx2
      dz=1./dx3
      do k=1,n3m
      do j=1,n2m
      do  i=1,n1m
       rhom=rho(i,j,k)+rhom
      enddo
      enddo
      enddo
      rhm=rhom*vl123
      rhorms=0.
      do k=1,n3m
      do j=1,n2m
      do i=1,n1m
       rhom=rho(i,j,k)-rhm
      rhorms=rhorms+rhom**2
      enddo
      enddo
      enddo
      rhorms=rhorms*vl123
      rhoske=0.
      rhofla=0.
      do k=1,n3m
      do j=1,n2m
      do i=1,n1m
       rhom=rho(i,j,k)-rhm
      rhoske=rhoske+(rhom/sqrt(rhorms))**3
      rhofla=rhofla+(rhom/sqrt(rhorms))**4
      enddo
      enddo
      enddo
      rhoske=rhoske*vl123
      rhofla=rhofla*vl123
      open(43,file='rhoflu.out')
      write(43,*)time,rhorms,rhoske,rhofla 
      close(43)
      return
      end
c
c  ****************************** subrout ouprho  **********************
c
      subroutine ouprho(time,nav,rho,qcap)
c
c    write the spectrum of the density
c
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      dimension rho(m1,m2,m3),qcap(m1,m2,m3)
      common/sprh/srho(0:m3)
      common/numond/kx(m1),ky(m2),kz(m3),kkmax
      common/kpmasp/kmax
      common/mtime/multim
      character*17 filcos
      character*4 pntim
      common/strat2/igrad,bvais
      common/strat1/istrat,rho0,g,schm
      common/spemw/akkpp,qq,sig
c
      itim=(time*multim+0.3)
      write(pntim,77) itim
   77 format(i4.4)
      filcos='sperho.'//pntim
      open(58,file=filcos)
      rewind(58)
      call rhphwa(rho,qcap)
      call sperho(qcap)
       akrap=1.
       if(igrad.eq.1.and.istrat.eq.1) akrap=akkpp
       do k=1,kmax
       ak=k/akrap
       write (58,*)ak,srho(k)
       end do
      close(58)
      return
      end
c*************************************************************************
      subroutine rhphwa(rho,qtil)
c
c     density fluct spectrum 
c     in this routine the density from physical to wave number
c
      include 'param.f'
      dimension rho(m1,m2,m3)
      complex qtil(m1,(m2-1)/2+1,m3)
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/fft001/ifxx1(13),trigxx1(2*(m1-1))
      common/fft002/ifxx2(13),trigxx2(3*(m2-1)/2+1)
      common/fft003/ifxx3(13),trigxx3(2*(m3-1))
      common/numond/kx(m1),ky(m2),kz(m3),kkmax
      complex xa(m3-1,m1-1),wor(m3-1,m1-1)
      complex xa2(m1-1,m3-1),wor2(m1-1,m3-1)
      real xr(m2+1,m1-1),work(m2,m1-1)
c
      n2mh=n2m/2+1
c
      do  10 k=1,n3m
c
        do i=1,n1m
         xr(1,i)=rho(i,n2m,k)
         xr(n2m+2,i)=rho(i,1,k)
         do j=1,n2m
          js=j+1
          xr(js,i)=rho(i,j,k)
         enddo
        enddo
c
c   real fft applied from  physical space to wave number 
c   in y direction
c
        call fft99(xr,work,trigxx2,ifxx2,1,m2+1,n2m,n1m,-1)
c
        do j=1,n2mh
         jp=2*j   
         jd=2*j-1
         do i=1,n1m
          qtil(i,j,k)=cmplx(xr(jd,i),xr(jp,i))
         enddo
        enddo
c
 10   continue
c
c   2-d  cfft applied (twice) from
c   physical space to wave number
c
      do 20 j=1,n2mh
c
        do k=1,n3m
         do i=1,n1m
          xa(k,i)=qtil(i,j,k)
         enddo
        enddo
c
        call cfft99(xa,wor,trigxx3,ifxx3,1,m3-1,n3m,n1m,-1)
c
        do k=1,n3m
         do i=1,n1m
          xa2(i,k)=xa(k,i)/float(n3m)
         enddo
        enddo
c
        call cfft99(xa2,wor2,trigxx1,ifxx1,1,m1-1,n1m,n3m,-1)
c
        do i=1,n1m
         do k=1,n3m
          qtil(i,j,k)=xa2(i,k)/float(n1m)
         enddo
        enddo
c
  20   continue
       return
       end 
c*************************************************************************
      subroutine sperho(qtil)
c
c    calculation of density fluct spectrum 
c
      include 'param.f'
      common/sprh/srho(0:m3)
      complex qtil(m1,(m2-1)/2+1,m3)
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/numond/kx(m1),ky(m2),kz(m3),kkmax
c
      n2mh=n2m/2+1
       do kk=1,kkmax
       srho(kk)=0.0
       enddo
       do k=1,n3m
        do j=2,n2mh
         do i=1,n1m
         akk=sqrt(float(kx(i)*kx(i)+ky(j)*ky(j)+kz(k)*kz(k)))
         kk=nint(akk)
         urea=real(qtil(i,j,k))
         uimm=aimag(qtil(i,j,k))
         srho(kk)=srho(kk)+(urea*urea+uimm*uimm)
         end do
        end do
       end do
       j=1
       do k=1,n3m
         do i=1,n1m
         akk=sqrt(float(kx(i)*kx(i)+ky(j)*ky(j)+kz(k)*kz(k)))
         kk=nint(akk)
         urea=real(qtil(i,j,k))
         uimm=aimag(qtil(i,j,k))
         srho(kk)=srho(kk)+(urea*urea+uimm*uimm)*0.5
         end do
       end do
       return
       end 
