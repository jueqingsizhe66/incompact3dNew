c
c  ****************************** subrout cfl  **********************
c
c  in this subroutine is calculated the cfl
c
      subroutine cfl(q,cflm)
c
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      common/tstep/dt,beta,ren
      common/metria/caj(m2),cac(m2)
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/newdat/icost,timeav
      dimension q(ndv,m1,m2,m3)
c
      cflm=0.
      cflma=0.
      do 7 k=1,n3m
      kp=kpv(k)
      do 7 j=1,n2m
      jp=j+1
      sucaj=1./caj(j)
      do 7 i=1,n1m
      ip=ipv(i)
      qcf=(abs((q(1,i,j,k)+q(1,ip,j,k))*dx1)+
     1     abs((q(2,i,j,k)+q(2,i,jp,k))*dx2*sucaj)+
     1     abs((q(3,i,j,k)+q(3,i,j,kp))*dx3))*.5
    7 cflma=max(cflma,qcf)
      cflm=cflma
      return
      end
c
c  ****************************** subrout coordi  **********************
c    this subroutine assign the physical coordinates as function
c    of the computational ones x2
c    the coordinate is clustered near the walls
c
      subroutine coordi(y)
c
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      dimension x(m1),eta(m2)
      common/strpar/str2
      dimension y(m2)
      common/d13/alx1,alx3
      common/y2sta/y2s(m2)
      common/y13sta/y1s(m1),y3s(m3)
      common/corpri/yp1(m1),yp2(m2),yp3(m3)
      common/jprts/ipri(5),kpri(5)
c
      tstr2=tanh(str2*0.5)
      do 63 j=1,n2
      x2=(j-1)/float(n2m)
      eta(j)=0.5*(1.+tanh(str2*(x2-0.5))/tstr2)
   63 continue
      do 65 j=1,n2
      yp2(j)=(-0.5+eta(j))*2.
      y(j)=yp2(j)
   65 continue
c     print *,'y al centro',y(n2/2+1)
      do 67 j=1,n2m
      y2s(j)=(yp2(j)+yp2(j+1))*0.5
   67 continue
c
c     simmetry of the  mesh enforced
c     when n2 odd
c
      if(mod(n2,2).eq.1) then
      y(n2/2+1)=0.
      do j=1,n2/2
      y(n2-j+1)=-y(j)
      end do
      do j=1,n2/2
      y2s(n2-j)=-y2s(j)
      end do
      end if
      open(47,file='ymesh')
      do j=1,n2m
      dy=y(j+1)-y(j) 
      write(47,*) y2s(j),dy
      end do
      close(47)
      do j=1,n2m/2+1
      ydp=(y(j)-y(1) )*180
      if(ydp.gt.20) then
      jpri=j
      go to 23
                    endif
      end do
   23 continue
c
c   this are poits where certain quantities are written at
c   each time step. This is to analize the time signals as
c   done in the experiments
c
      kpri(1)=n3m/4+1
      ipri(1)=n1m/4+1
      kpri(2)=n3m/4+1
      ipri(2)=3*n1m/4+1
      kpri(3)=n3m/2+1
      ipri(3)=n1m/2+1
      kpri(4)=3*n3m/4+1
      ipri(4)=n1m/4+1
      kpri(5)=3*n3m/4+1
      ipri(5)=3*n1m/4+1
      jin=1
      jfi=n2m/2
      do l=1,5
      iin=ipri(l)-3
      ifi=ipri(l)+3
      kin=kpri(l)-3
      kfi=kpri(l)+3
      npoin=0
      do i =iin,ifi
      do j =jin,jfi
      do k =kin,kfi
      npoin=npoin+1
      enddo
      enddo
      enddo
      enddo
      write(6,*)' print npoin=',npoin
c
      do  i=1,n1
          yp1(i) = float(i-1)/float(n1m) * alx1
      enddo
      do  k=1,n3
          yp3(k) = float(k-1)/float(n3m) * alx3
      enddo
      do  i=1,n1m
          y1s(i) = 0.5*(yp1(i) + yp1(i+1) )
      enddo
      do  k=1,n3m
          y3s(k) = 0.5*(yp3(k) + yp3(k+1) )
      enddo
      return
      end

c
c  ****************************** subrout divg  **********************
c
c  this subroutine perform the calculation of divg(q)
c
      subroutine divg(qcap,vq,al)
c
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      dimension qcap(m1,m2,m3),vq(ndv,m1,m2,m3)
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/tstep/dt,beta,ren
      common/metria/caj(m2),cac(m2)
c
c  ***** compute the divg(u) 
      sudtal=1./(dt*al)
      do 11 kc=1,n3m
      kp=kpv(kc)
      do 11 jc=1,n2m
      jp=jc+1
      sucaj=1./caj(jc)
      do 11 ic=1,n1m
      ip=ipv(ic)
      dqcap=(vq(1,ip,jc,kc)-vq(1,ic,jc,kc))*dx1
     1     +(vq(2,ic,jp,kc)-vq(2,ic,jc,kc))*dx2*sucaj
     1     +(vq(3,ic,jc,kp)-vq(3,ic,jc,kc))*dx3
      qcap(ic,jc,kc)=dqcap*sudtal
   11 continue
      return
      end
c  ****************************** subrout divgck  **********************
c
c  this subroutine perform the calculation of divg(q) and checks
c  whether the local mass is conserved
c
      subroutine divgck(vq,qmax)
c
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      dimension vq(ndv,m1,m2,m3)
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/metria/caj(m2),cac(m2)
c
c  ***** compute the divg(q)
      qmax=0.
      do 11 kc=1,n3m
      kp=kpv(kc)
      do 11 jc=1,n2m
      jp=jc+1
      sucaj=1./caj(jc)
      do 11 ic=1,n1m
      ip=ipv(ic)
      dqcap=(vq(1,ip,jc,kc)-vq(1,ic,jc,kc))*dx1
     1     +(vq(2,ic,jp,kc)-vq(2,ic,jc,kc))*dx2*sucaj
     1     +(vq(3,ic,jc,kp)-vq(3,ic,jc,kc))*dx3
      qmax=max(abs(dqcap),qmax)
   11 continue
      return
      end
c
c  ****************************** subrout indic **********************
c
c  in this subroutine the indices ip,im,jp,jm,kp,km are calculated
c  these are necessary when the equation are solved near the
c  walls and for the periodic conditions.
c
      subroutine indic
c
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/injmp/jpppv(m2),jmmmv(m2)
      common/inimp/ipppv(m1),immmv(m1)
c
c
c   periodic directions
c
      n1mm=n1m-1
      do 11 ic=1,n1m
      ipv(ic)=ic+1
      if(ic.eq.n1m) ipv(ic)=1
      imv(ic)=ic-1
      if(ic.eq.1) imv(ic)=n1m
      ip=ipv(ic)
      ipppv(ic)=ip+1
      if(ic.eq.n1mm) ipppv(ic)=1
      im=imv(ic)
      immmv(ic)=im-1
      if(ic.eq.2) immmv(ic)=n1m
   11 continue
      do 2 kc=1,n3m
      kmv(kc)=kc-1
      kpv(kc)=kc+1
      if(kc.eq.1) kmv(kc)=n3m
      if(kc.eq.n3m) kpv(kc)=1
    2 continue
c
c   direction normal to non-slip walls
c
      do 3 jc=1,n2m
      jp=jc+1
      jpppv(jc)=jp+1
      if(jc.eq.n2m) jpppv(jc)=n2
      jmv(jc)=jc-1
      jpv(jc)=jc+1
      if(jc.eq.1) jmv(jc)=jc
      if(jc.eq.n2m) jpv(jc)=jc
    3 continue
      do 5 jc=2,n2m
      jmmmv(jc)=jc-2
      if(jc.eq.2) jmmmv(jc)=1
    5 continue
      return
      end
c
c
c  ****************************** subrout meshes **********************
c
c  generates the mesh inverse of the spatial steps dx1, dx2, dx3,
c  and the squares of the dx1 and dx2.
c
      subroutine meshes
c
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      common/d13/alx1,alx3
      common/meshu/d1x1,d1x2,d1x3
c
      d1x2=1./float(n2m)
      d1x1=alx1/float(n1m)
      d1x3=alx3/float(n3m)
      dx1=1./d1x1
      dx2=1./d1x2
      dx3=1./d1x3
      dx1q=dx1*dx1
      dx2q=dx2*dx2
      dx3q=dx3*dx3
      return
      end
c
c
c  ****************************** subrout metric  **********************
c
c  this subroutine performs the calculation of the metric quantities
c  at j ; j+1/2 
c  ** the quantities with  j(  evaluated at .,j+1/2,.     positiion
c  ** the quantities with  c(  evaluated at .,j,.         positiion
c  ** the metric are evaluated by centered differences.
c
      subroutine metric(y)
c
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      common/metria/caj(m2),cac(m2)
      dimension y(m2)
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/d13/alx1,alx3
c
c  *********                 i,j+1/2
c
      do 1 j=1,n2m
c
c  ********* derivatives of physical coordinate interior points
c   cn2 deriv respect to x2
c
      cn22=(y(j+1)-y(j))*dx2
   10 continue
      caj(j)=cn22
    1 continue
c
c  *********                 i,j
c
      do 6 j=1,n2
      if(j.gt.1.and.j.lt.n2) go to 31
c
c  ********* derivatives of physical corrdinate horizontal boundary
c            three points backw.
      if(j.eq.1) js=1
      if(j.eq.n2) js=-1
      cn22=-js*(y(j)-y(j+js))*dx2
      go to 34
   31 continue
c
c  ********* derivatives of physical corrdinate interior points
c
      cn22=(y(j+1)-y(j-1))*dx2*0.5
   60 continue
   34 continue
      cac(j)=cn22
    6 continue
      return
      end
