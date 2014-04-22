c
      program main
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/dimfi/n1fi,n1mfi,n2fi,n2mfi,n3fi,n3mfi
      common/d123/alx1,alx2,alx3
      common/nfilt/nr1,nr2,nr3
      common/pias/pi
      common/ispec/imic
      common/vpeini/vper,omtres
      common/spemw/akkpp,qq,sig
      common/itypfi/ifilcu,ifilt
      common/reyn/re
      common/prinkm/kpmax,kfpmax
      dimension qmaxr(3)
      character*60 namfile

c
      open(15,file='filtrobook.d')
      write(6,*)'read from filles.d'
      read(15,*) n1,n2,n3,alx1d,alx2d,alx3d
c     write(6,*) n1,n2,n3,alx1d,alx2d,alx3d
      read(15,*) nr1,nr2,nr3,ifilcu,ifilt
c     write(6,*) nr1,nr2,nr3
      read(15,*) imic
c     write(6,*) imic
      if(imic.ge.0) then
      read(15,*) akkpp,qq,sig
c     write(6,*) akkpp,qq,sig
                    endif
      pi=2.*asin(1.)
      n1m=n1-1
      n2m=n2-1
      n3m=n3-1
      n1mfi=n1m/nr1
      n2mfi=n2m/nr2
      n3mfi=n3m/nr3
      n1fi=n1mfi+1
      n2fi=n2mfi+1
      n3fi=n3mfi+1
      kpmax=n1m/2
      kfpmax=n1mfi/2
      alx1=alx1d*pi
      alx2=alx2d*pi
      alx3=alx3d*pi
      write(6,112)alx1d,alx2d,alx3d
  112 format(10x,'box dimension',3x,'lx=',f4.2,'*pi'
     1        ,3x,'lz=',f4.2,'*pi',3x,'lz=',f4.2,'*pi')
      write(6,*) '**************************************'
      write(6,200)
  200 format(10x,'3d isotropic turbulence filtering check')
      write(6,*) '**************************************'
      write(6,*) 'ifilcu=',ifilcu 
      call gcurv
      stop
      end
c
c  ************************************************************
c  ************subrout gcurv **********************************
c  ************************************************************
c
      subroutine gcurv
c
c
      include 'param.f'
      parameter (m1m=m1-1)
c
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      dimension qfi(ndv,m1,m2,m3),rhsfi(m1,m2,m3)
      common/rhsc/rhs(m1,m2,m3)
      dimension q(ndv,m1,m2,m3),qcap(m1,m2,m3)
     1         ,dph(m1,m2,m3),pr(m1,m2,m3)
      common/dimfi/n1fi,n1mfi,n2fi,n2mfi,n3fi,n3mfi
      common/velmax/vmax(ndv),vmaxo(ndv)
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/d123/alx1,alx2,alx3
      common/speene/e(ndv,0:2*m3)
      dimension eq(0:2*m3),eqfi(0:2*m3),eqre(0:2*m3)
      dimension eqq(ndv,0:2*m3),eqqre(ndv,0:2*m3),eqqfi(ndv,0:2*m3)
      common/numond/kx(m1),ky(m2),kz(m3),kkmax
      common/itypfi/ifilcu,ifilt
      common/ispec/imic
      common/pias/pi
      common/reyn/re
      common/prinkm/kpmax,kfpmax
      istop=0
      nat=0
c
c
      npfile=npin*npstf
c
c     step  and mesh sizes calculations
c
c
      call meshes
      call indic
      call coordi
c
      write(6,754)n1m,n2m,n3m
  754 format(10x,'original grid      ',2x
     1      ,5x,'n1m=',i5,2x,'n2m=',i5,2x,'n3m=',i5)
      write(6,756)n1mfi,n2mfi,n3mfi
  756 format(10x,'filtered at        ',2x
     1      ,5x,'n1m=',i5,2x,'n2m=',i5,2x,'n3m=',i5)
      write(6,755) dx1,dx2,dx3
  755 format(3x,'dx1=',e10.3,3x,'dx2=',e10.3,3x
     1      ,3x,'dx3=',e10.3,3x)
      if(ifilt.eq.0) write(6,*)'  test physical filter'
      if(ifilt.eq.1) write(6,*)'  test spectral filter'
      n1mh=n1m/2+1
c
c     initial conditions
c it looks like speini and phini both initialize fft, perhaps speini is redundant 
c here. Check.
c
      call speini
      call phini
      do 22 l=1,ndv
      vmax(l)=0.1e-10
   22 continue
c
c
      if(imic.gt.0) then
      call initur(q,pr,qcap,dph)
                    else
      call inirea(q,pr)
                    endif
      call spectre(q,pr,qcap,dph)
       open(77,file='spec.0000')
       open(78,file='spet.0000')
       do k=1,kpmax
       write (77,*)k,e(1,k),e(2,k),e(3,k)
       write (78,*)k,e(1,k)+e(2,k)+e(3,k)
       end do
       close(77)
       close(78)
c
c    velocity component 1    
c
      do i=1,n1m
           do j=1,n2m
                do k=1,n3m
        rhs(i,j,k)=q(1,i,j,k)
                enddo
           enddo
      enddo
      if(ifilt.eq.1)  then
      call filspe(rhsfi,pr,eq,eqfi,eqre)
       do k=1,kkmax
       eqq(1,k)=eq(k)
       eqqre(1,k)=eqre(k)
       eqqfi(1,k)=eqfi(k)
       end do
                      else
      call filphy(rhsfi)
                      endif
      do i=1,n1m
           do j=1,n2m
                do k=1,n3m
        qfi(1,i,j,k)=rhsfi(i,j,k)
                enddo
           enddo
      enddo
c
c    velocity component 2    
c
      do i=1,n1m
           do j=1,n2m
                do k=1,n3m
        rhs(i,j,k)=q(2,i,j,k)
                enddo
           enddo
      enddo
      if(ifilt.eq.1)  then
      call filspe(rhsfi,pr,eq,eqfi,eqre)
       do k=1,kkmax
       eqq(2,k)=eq(k)
       eqqre(2,k)=eqre(k)
       eqqfi(2,k)=eqfi(k)
       end do
                      else
      call filphy(rhsfi)
                      endif
      do i=1,n1m
           do j=1,n2m
                do k=1,n3m
        qfi(2,i,j,k)=rhsfi(i,j,k)
                enddo
           enddo
      enddo
c
c    velocity component 3    
c
      do i=1,n1m
           do j=1,n2m
                do k=1,n3m
        rhs(i,j,k)=q(3,i,j,k)
                enddo
           enddo
      enddo
      if(ifilt.eq.1)  then
      call filspe(rhsfi,pr,eq,eqfi,eqre)
       do k=1,kkmax
       eqq(3,k)=eq(k)
       eqqre(3,k)=eqre(k)
       eqqfi(3,k)=eqfi(k)
       end do
                      else
      call filphy(rhsfi)
                      endif
      do i=1,n1m
           do j=1,n2m
                do k=1,n3m
        qfi(3,i,j,k)=rhsfi(i,j,k)
                enddo
           enddo
      enddo
      call divc(qfi,qmax)
      print *,'divergenza massima=',qmax
       if(ifilt.eq.1) then
       open(77,file='specsp.0000')
       open(78,file='spetsp.0000')
       do k=1,kpmax
       write (77,*)k,eqq(1,k),eqq(2,k),eqq(3,k)
       write (78,*)k,eqq(1,k)+eqq(2,k)+eqq(3,k)
       end do
       close(77)
       close(78)
       open(77,file='specrefi.0000')
       open(78,file='spetrefi.0000')
       do k=1,kpmax
       write (77,*)k,eqqre(1,k),eqqre(2,k),eqqre(3,k)
       end do
       close(77)
       close(78)
       open(77,file='specnofi.0000')
       open(78,file='spetnofi.0000')
       do k=1,kpmax
       write (77,*)k,eqqfi(1,k),eqqfi(2,k),eqqfi(3,k)
       write (78,*)k,eqqfi(1,k)+eqqfi(2,k)+eqqfi(3,k)
       end do
       close(77)
       close(78)
                      endif
      call spectre(qfi,pr,qcap,dph)
       if(ifilt.eq.1) then
       open(77,file='specspfi.0000')
       open(78,file='spetspfi.0000')
                      else
       open(77,file='specphfi.0000')
       open(78,file='spetphfi.0000')
                      endif
       do k=1,kfpmax
       write (77,*)k,e(1,k),e(2,k),e(3,k)
       write (78,*)k,e(1,k)+e(2,k)+e(3,k)
       end do
       close(77)
       close(78)
      return
      end
c
c  ****************************** subrout inirea **********************
c
c   reading initial  conditions from file
c   when calculation are continued from a previous
c   calculation.
c
      subroutine inirea(q,pr)
c
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      dimension pr(m1,m2,m3)
      dimension q(ndv,m1,m2,m3)
      common/reyn/re
      character*60 namfile
c
      nfil=23
      read(15,100)namfile
      read(15,*)re
  100 format(a60)
      open(23,file=namfile,form='unformatted')
      rewind(nfil)
      read(nfil) n1,n2,n3
      read(nfil) time,ren,time,dum
      read(nfil) (((q(1,i,j,k),i=1,n1),j=1,n2),k=1,n3),
     1            (((q(2,i,j,k),i=1,n1),j=1,n2),k=1,n3),
     1            (((q(3,i,j,k),i=1,n1),j=1,n2),k=1,n3)
     1           ,(((pr(i,j,k),i=1,n1),j=1,n2),k=1,n3)
      close(nfil)
      call divc(q,qmax)
      print *,'divergenza massima=',qmax
      return
      end

c***********************************************************
      subroutine speini
      include 'param.f'
      common/fft001/ifxx1(13),trigxx1(2*(m1-1))
      common/fft002/ifxx2(13),trigxx2(3*(m2-1)/2+1)
      common/fft003/ifxx3(13),trigxx3(2*(m3-1))
      common/dim/n1,n1m,n2,n2m,n3,n3m
      nx3fft=n3m
      call cftfax(nx3fft,ifxx3,trigxx3)
      nx1fft=n1m
      call cftfax(nx1fft,ifxx1,trigxx1)
      nx2fft=n2m
      call fftfax(nx2fft,ifxx2,trigxx2)
      return 
      end
c
c  ****************************** subrout coordi  **********************
c    this subroutine assign the physical coordinates as function
c    of the computational ones x1,x2
c
      subroutine coordi
c
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/d123/alx1,alx2,alx3
      common/corpri/yp1(m1),yp2(m2),yp3(m3)
      common/corpst/ys1(m1),ys2(m2),ys3(m3)
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      common/ispec/imic
      common/pias/pi
c
      do i=1,n1
      x1=float(i-1)/float(n1m)
      yp1(i)=-alx1*0.5+x1*alx1
      enddo
c
      do j=1,n2
      x2=float(j-1)/float(n2m)
      yp2(j)=-alx2*0.5+x2*alx2
      enddo
c
      do k=1,n3
      x3=float(k-1)/float(n3m)
      yp3(k)=-alx3*0.5+x3*alx3
      enddo
c
c
c    at the centre
c
      do i=1,n1m
      x1=(float(i-1)+0.5)/float(n1m)
      ys1(i)=-alx1*0.5+x1*alx1
      enddo
c
      do j=1,n2m
      x2=(float(j-1)+0.5)/float(n2m)
      ys2(j)=-alx2*0.5+x2*alx2
      enddo
c
      do k=1,n3m
      x3=(float(k-1)+0.5)/float(n3m)
      ys3(k)=-alx3*0.5+x3*alx3
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
      dt=.1
c
c  ***** compute the divg(u)
      do 11 kc=1,n3m
      kp=kpv(kc)
      do 11 jc=1,n2m
      jp=jpv(jc)
      do 11 ic=1,n1m
      ip=ipv(ic)
      dqcap=(vq(1,ip,jc,kc)-vq(1,ic,jc,kc))*dx1
     1     +(vq(2,ic,jp,kc)-vq(2,ic,jc,kc))*dx2
     1     +(vq(3,ic,jc,kp)-vq(3,ic,jc,kc))*dx3
      qcap(ic,jc,kc)=dqcap/(dt*al)
   11 continue
      return
      end
c  ****************************** subrout updvp  **********************
c
c  this subroutine calculate the solenoidal vel field
c       q(n+1)=q-grad(dph)*dt ,  pr=dph
c  third order runge-kutta is used.
c
      subroutine updvp(dph,q,al)
c
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      dimension dph(m1,m2,m3),q(ndv,m1,m2,m3)
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      dt=.1
c
c  ***********  compute the q1 velocity component
c               v1dgf=component 1 of grad(dph)
      do 1 kc=1,n3m
      do 1 jc=1,n2m
      do 1 ic=1,n1m
      im=imv(ic)
      dfx11=(dph(ic,jc,kc)-dph(im,jc,kc))*dx1
      q(1,ic,jc,kc)=q(1,ic,jc,kc)-dfx11*dt*al
    1 continue
c
c  ***********  compute the q2 velocity component
c               v2dgf=component 2 of grad(dph)
      do 2 kc=1,n3m
      do 2 jc=1,n2m
      do 2 ic=1,n1m
      jm=jmv(jc)
      dfx22=(dph(ic,jc,kc)-dph(ic,jm,kc))*dx2
      q(2,ic,jc,kc)=q(2,ic,jc,kc)-dfx22*dt*al
    2 continue
c
c  ***********  compute the q3 velocity component
c               q3 is the cartesian component
c               v3dgf=component 3 of grad(dph)
      do 5 kc=1,n3m
      km=kmv(kc)
      do 5 jc=1,n2m
      do 5 ic=1,n1m
      dfx33=(dph(ic,jc,kc)-dph(ic,jc,km))*dx3
      q(3,ic,jc,kc)=q(3,ic,jc,kc)-dfx33*al*dt
    5 continue
      return
      end
c
c  ****************************** subrout divc  **********************
c
c  this subroutine perform the calculation of divg(q)
c
      subroutine divc(vq,qmax)
c
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      dimension vq(ndv,m1,m2,m3)
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
c
c  ***** compute the divg(q)
      qmax=0.
      do 11 kc=1,n3m
      kp=kpv(kc)
      do 11 jc=1,n2m
      jp=jpv(jc)
      do 11 ic=1,n1m
      ip=ipv(ic)
      dqcap=(vq(1,ip,jc,kc)-vq(1,ic,jc,kc))*dx1
     1     +(vq(2,ic,jp,kc)-vq(2,ic,jc,kc))*dx2
     1     +(vq(3,ic,jc,kp)-vq(3,ic,jc,kc))*dx3
      qmax=amax1(abs(dqcap),qmax)
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
c
c
c   periodic directions
c
      do 1 ic=1,n1m
      ipv(ic)=ic+1
      if(ic.eq.n1m) ipv(ic)=1
      imv(ic)=ic-1
      if(ic.eq.1) imv(ic)=n1m
    1 continue
      do 2 kc=1,n3m
      kmv(kc)=kc-1
      kpv(kc)=kc+1
      if(kc.eq.1) kmv(kc)=n3m
      if(kc.eq.n3m) kpv(kc)=1
    2 continue
      do 3 jc=1,n2m
      jpv(jc)=jc+1
      if(jc.eq.n2m) jpv(jc)=1
      jmv(jc)=jc-1
      if(jc.eq.1) jmv(jc)=n2m
    3 continue
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
      common/d123/alx1,alx2,alx3
      common/meshu/d1x1,d1x2,d1x3
c
      d1x1=alx1/float(n1m)
      d1x2=alx2/float(n2m)
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
c  ****************************** subrout vmaxv **********************
c
      subroutine vmaxv(q)
c
      include 'param.f'
      dimension q(ndv,m1,m2,m3)
      common/velmax/vmax(ndv),vmaxo(ndv)
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/metria/caj(m2),cac(m2)
c
c  ***** calculation of the maximum velocities
c  in order to check convergency or for stability calculations
c  to derive stability conditions.
c
      do 311 l=1,ndv
      vmaxo(l)=vmax(l)
      vmax(l)=0.
      do 310 k=1,n3m
      do 310 j=1,n2m
      do 310 i=1,n1m
      if(l.eq.1) vca=q(1,i,j,k)
      if(l.eq.2) vca=q(2,i,j,k)
      if(l.eq.3) vca=q(3,i,j,k)
      vfm=abs(vca)
      vm=vmax(l)
      vmax(l)=amax1(vm,vfm)
  310 continue
  311 continue
      return
      end
c
c  ****************************** subrout initur **********************
c
c   initial zero conditions in the whole field.
c   and at the inner and outer walls for dph(i,j,k),
c   pr(i,j,k), ru(l,i,j,k) and h(l,i,j,k)
c   set up  dp3ns.
c
      subroutine initur(q,pr,qcap,dph)
c
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      dimension q(ndv,m1,m2,m3),pr(m1,m2,m3)
      dimension qcap(m1,m2,m3),dph(m1,m2,m3)
      common/corpri/yp1(m1),yp2(m2),yp3(m3)
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/speene/e(ndv,0:2*m3)
      common/numond/kx(m1),ky(m2),kz(m3),kkmax
      common/pias/pi
c call mic to make initial conditions. Use qcap,dph and pr as work files
      call turbis(qcap,dph,pr,q)
        al=1.
      call divg(qcap,q,al)
c qcap now contains the divergence of the velocity field q
c
c  ********* calculation of the pressure dph by fft in two
c            directions and tridiag in vertical
      call phcalc(qcap,dph)
c
c  ********* calculation of solenoidal vel field
c
      call updvp(dph,q,1.)
      vm1m=0.
      vm2m=0.
      vm3m=0.
      vl123=1./float(n1m*n2m*n3m)
      do 411 k=1,n3m
      do 411 j=1,n2m
      do 411 i=1,n1m
       vm1m=q(1,i,j,k)+vm1m
       vm2m=q(2,i,j,k)+vm2m
       vm3m=q(3,i,j,k)+vm3m
  411 continue
      vm1=vm1m*vl123
      vm2=vm2m*vl123
      vm3=vm3m*vl123
      enej=0.
      print*,'vm1,vm2,vm3'
      print*,vm1,vm2,vm3
      do 414 k=1,n3m
      do 414 j=1,n2m
      do 414 i=1,n1m
       q(1,i,j,k) =q(1,i,j,k)-vm1
       q(2,i,j,k) =q(2,i,j,k)-vm2
       q(3,i,j,k) =q(3,i,j,k)-vm3
414   continue
      vm1m=0.
      vm2m=0.
      vm3m=0.
      vl123=1./float(n1m*n2m*n3m)
      do 511 k=1,n3m
      do 511 j=1,n2m
      do 511 i=1,n1m
       vm1m=q(1,i,j,k)+vm1m
       vm2m=q(2,i,j,k)+vm2m
       vm3m=q(3,i,j,k)+vm3m
  511 continue
      vm1=vm1m*vl123
      vm2=vm2m*vl123
      vm3=vm3m*vl123
      print*,'vm1,vm2,vm3'
      print*,vm1,vm2,vm3
      call divc(q,qmax)
      print *,'divergenza massima=',qmax
      return
      end
*************************************************************************
      subroutine turbis(qtil1,qtil2,qtil3,q)
c make intitial conditions
      include 'param.f'
      dimension q(ndv,m1,m2,m3)
      complex qtil1(m1,(m2-1)/2+1,m3)
      complex qtil2(m1,(m2-1)/2+1,m3)
      complex qtil3(m1,(m2-1)/2+1,m3)
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/numond/kx(m1),ky(m2),kz(m3),kkmax
      common/rot/f0
      common/strat1/istrat,rho0,g,schm
      common/strat2/igrad,bvais
      common/d123/alx1,alx2,alx3
      complex alp,bet,cexp1,cexp2,cexp3,cexp4
      dimension ekk(m1),e(3,m1)
      common/ispec/imic
      common/spemw/akkpp,qq,sig
      common/pias/pi
c
      n2mh=n2m/2+1
       if(imic.eq.0)  then
       open(21,file='speread.dat')
       do kk=1,500
       read(21,*,end=211)k,e(1,k),e(2,k),e(3,k)
       ekk(k)=e(1,k)+e(2,k)+e(3,k)
       enddo
  211 continue
                      endif
       if(imic.eq.2)  then
       
c  here we try  the initial spectrum of Mansour and Wray
c               DNS
       aa=0.
       do kk=1,kkmax
       ak=kk
c      rk=ak/akkpp
       rk=ak
       aint=rk**sig*exp(-sig*rk**2/2.)
       aa=aa+aint
       enddo
       aa1=0.
       eto=0.
       do kk=1,kkmax
       ak=kk
       rk=ak/akkpp
       ek0kl=qq/(2.*aa)/akkpp
       ek1kl=rk**sig*exp(-sig*rk**2/2.)
       aa1=aa1+ek1kl
       ekk(kk)=ek0kl*ek1kl
       eto=eto+ekk(kk)
       enddo
       write(6,*)'initial DNS Mans Wray spectrum'
       write(6,*)'     akkpp,qq,sig',akkpp,qq,sig
       write(6,*)'     aa,aa1,eto',aa,aa1,eto
                       endif
       isd=17
       call srand(isd)
c      call ranset(isd)
       do k=1,n3m
            akz=kz(k)
        do j=2,n2mh
            aky=ky(j)
         do i=1,n1m
            akx=kx(i)
            fk2=float(kx(i)*kx(i)+ky(j)*ky(j)+kz(k)*kz(k))
            fk=sqrt(fk2)
            kk=fk
            fkhor2=float(kx(i)*kx(i)+ky(j)*ky(j))
            fkhor=sqrt(fk2)
            ampk=sqrt(ekk(kk)/(4.*pi*fk2))
c           x1=rnd()
c           x2=rnd()
c           x3=rnd()
            x1=rand()
            x2=rand()
            x3=rand()
c           x1=rand()
c           x2=rand()
c           x3=rand()
            phi=2.*pi*x1
            the1=2.*pi*x2
            the2=2.*pi*x3
            cexp1=cmplx(cos(the1),sin(the1))
            cexp2=cmplx(cos(the2),sin(the2))
            alp=ampk*cexp1*cos(phi)
            bet=ampk*cexp2*sin(phi)
            qtil1(i,j,k)=(alp*fk*aky+bet*akx*aky)/(fk*fkhor)
            qtil2(i,j,k)=(-alp*fk*akx+bet*akz*aky)/(fk*fkhor)
            qtil3(i,j,k)=-(bet*fkhor)/fk
         end do
        end do
       end do
       write(6,*)' first do'
       j=1
            aky=ky(j)
       do k=1,n3m
            akz=kz(k)
         do i=2,n1m
            akx=kx(i)
            fk2=float(kx(i)*kx(i)+ky(j)*ky(j)+kz(k)*kz(k))
            fk=sqrt(fk2)
            kk=fk
            fkhor2=float(kx(i)*kx(i)+ky(j)*ky(j))
            fkhor=sqrt(fk2)
            ampk=sqrt(ekk(kk)/(4.*pi*fk2))
c           x1=rnd()
c           x2=rnd()
c           x3=rnd()
            x1=rand()
            x2=rand()
            x3=rand()
c           x1=rand()
c           x2=rand()
c           x3=rand()
            phi=2.*pi*x1
            the1=2.*pi*x2
            the2=2.*pi*x3
            cexp1=cmplx(cos(the1),sin(the1))
            cexp2=cmplx(cos(the2),sin(the2))
            alp=ampk*cexp1*cos(phi)
            bet=ampk*cexp2*sin(phi)
            qtil1(i,j,k)=(alp*fk*aky+bet*akx*aky)/(fk*fkhor)
            qtil2(i,j,k)=(-alp*fk*akx+bet*akz*aky)/(fk*fkhor)
            qtil3(i,j,k)=-(bet*fkhor)/fk
         end do
       end do
       write(6,*)' second do'

       
212    continue
      call fft3d(q,1,qtil1,+1)
      call fft3d(q,2,qtil2,+1)
      call fft3d(q,3,qtil3,+1)
      return
      end
c***********************************************************   
      subroutine fft3d(a,ll,b,ifft)
c     energy spectrum 
      include 'param.f'
      dimension a(3,m1,m2,m3)
      complex b(m1,(m2-1)/2+1,m3)
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/fft001/ifxx1(13),trigxx1(2*(m1-1))
      common/fft002/ifxx2(13),trigxx2(3*(m2-1)/2+1)
      common/fft003/ifxx3(13),trigxx3(2*(m3-1))
      common/numond/kx(m1),ky(m2),kz(m3),kkmax
      common/speene/e(ndv,0:2*m3)
      complex xa(m3-1,m1-1),wor(m3-1,m1-1)
      complex xa2(m1-1,m3-1),wor2(m1-1,m3-1)
      real xr(m2+1,m1-1),work(m2,m1-1)
c
c     print *,'calcolo dello spettro'
      n2mh=n2m/2+1
c
      if(ifft.eq.-1) then
c     do 1 n=1,3
      do  10 k=1,n3m
c
        do i=1,n1m
         xr(1,i)=a(ll,i,n2m,k)
         xr(n2m+2,i)=a(ll,i,1,k)
         do j=1,n2m
          js=j+1
          xr(js,i)=a(ll,i,j,k)
         enddo
        enddo
c
c   real fft applied from  physical space to wave number 
c   in y direction
c
        call fft99(xr,work,trigxx2,ifxx2,1,m2+1,n2m,n1m,ifft)
c
        do j=1,n2mh
         jp=2*j   
         jd=2*j-1
         do i=1,n1m
          b(i,j,k)=cmplx(xr(jd,i),xr(jp,i))
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
          xa(k,i)=b(i,j,k)
         enddo
        enddo
c
        call cfft99(xa,wor,trigxx3,ifxx3,1,m3-1,n3m,n1m,ifft)
c
        do k=1,n3m
         do i=1,n1m
          xa2(i,k)=xa(k,i)/float(n3m)
         enddo
        enddo
c
        call cfft99(xa2,wor2,trigxx1,ifxx1,1,m1-1,n1m,n3m,ifft)
c
        do i=1,n1m
         do k=1,n3m
          b(i,j,k)=xa2(i,k)/float(n1m)
         enddo
        enddo
c
  20   continue
c
c     b is the velocity component in Fourier space.
c
      
     

      else

      do 41 j=1,n2mh
c
        do i=1,n1m
         do k=1,n3m
          xa2(i,k)=b(i,j,k)
         enddo
        enddo
        call cfft99(xa2,wor2,trigxx1,ifxx1,1,m1-1,n1m,n3m,ifft)
        do k=1,n3m
         do i=1,n1m
          xa(k,i)=xa2(i,k)
         enddo
        enddo
c
        call cfft99(xa,wor,trigxx3,ifxx3,1,m3-1,n3m,n1m,ifft)
c
c
        do i=1,n1m
         do k=1,n3m
          b(i,j,k)=xa(k,i)
         enddo
        enddo
c
c
  41   continue
      do  40 k=1,n3m
c
        do j=1,n2mh
         jp=2*j   
         jd=2*j-1
         do i=1,n1m
          xr(jd,i)=real(b(i,j,k)) 
          xr(jp,i)=aimag(b(i,j,k)) 
         enddo
        enddo
c
c   real fft applied from  physical space to wave number 
c   in y direction
c
        call fft99(xr,work,trigxx2,ifxx2,1,m2+1,n2m,n1m,ifft)
        do i=1,n1m
         do j=1,n2m
          js=j+1
          a(ll,i,j,k)=xr(js,i)
         enddo
        enddo
c
c
 40   continue
       endif
       return
       end 
c
c    ********************  subr spectre
c
      subroutine spectre(q,qtil1,qtil2,qtil3)
c     energy spectrum 
      include 'param.f'
      dimension q(ndv,m1,m2,m3)
      complex qtil1(m1,(m2-1)/2+1,m3)
      complex qtil2(m1,(m2-1)/2+1,m3)
      complex qtil3(m1,(m2-1)/2+1,m3)
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/fft001/ifxx1(13),trigxx1(2*(m1-1))
      common/fft002/ifxx2(13),trigxx2(3*(m2-1)/2+1)
      common/fft003/ifxx3(13),trigxx3(2*(m3-1))
      common/numond/kx(m1),ky(m2),kz(m3),kkmax
      common/speene/e(ndv,0:2*m3)
      complex xa(m3-1,m1-1),wor(m3-1,m1-1)
      complex xa2(m1-1,m3-1),wor2(m1-1,m3-1)
      real xr(m2+1,m1-1),work(m2,m1-1)
c
      
      n2mh=n2m/2+1
c
      do 1 n=1,3
      do  10 k=1,n3m
c
        do i=1,n1m
         xr(1,i)=q(n,i,n2m,k)
         xr(n2m+2,i)=q(n,i,1,k)
         do j=1,n2m
          js=j+1
          xr(js,i)=q(n,i,j,k)
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
          qtil3(i,j,k)=cmplx(xr(jd,i),xr(jp,i))
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
          xa(k,i)=qtil3(i,j,k)
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
          qtil3(i,j,k)=xa2(i,k)/float(n1m)
         enddo
        enddo
c
  20   continue
       do kk=1,kkmax
       e(n,kk)=0.0
       enddo
       do k=1,n3m
        do j=2,n2mh
         do i=1,n1m
         akk=sqrt(float(kx(i)*kx(i)+ky(j)*ky(j)+kz(k)*kz(k)))
         kk=nint(akk)
         urea=real(qtil3(i,j,k))
         uimm=aimag(qtil3(i,j,k))
         e(n,kk)=e(n,kk)+2.*(urea*urea+uimm*uimm)
         end do
        end do
       end do
       j=1
       do k=1,n3m
         do i=1,n1m
         akk=sqrt(float(kx(i)*kx(i)+ky(j)*ky(j)+kz(k)*kz(k)))
         kk=nint(akk)
         urea=real(qtil3(i,j,k))
         uimm=aimag(qtil3(i,j,k))
         e(n,kk)=e(n,kk)+(urea*urea+uimm*uimm)
         end do
       end do
       if(n.eq.1) then
       do j=1,n2mh
           do k=1,n3m
                do i=1,n1m
           qtil1(i,j,k)=qtil3(i,j,k)
               end do
           end do
       end do
                  endif
       if(n.eq.2) then
       do j=1,n2mh
           do k=1,n3m
                do i=1,n1m
           qtil2(i,j,k)=qtil3(i,j,k)
               end do
           end do
       end do
                  endif
  1    continue
       return
       end 
c
c******************   subroutine filspe
c
      subroutine filspe(uf,qtil,eq,eqfi,eqre)
c     3d sharp Fourier cutoff filter
c     u is the resolved field
c     uf is the filtered field
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/rhsc/u(m1,m2,m3)
      dimension uf(m1,m2,m3)
      dimension eq(0:2*m3),eqfi(0:2*m3),eqre(0:2*m3)
      complex qtil(m1,(m2-1)/2+1,m3)
      common/fft001/ifxx1(13),trigxx1(2*(m1-1))
      common/fft002/ifxx2(13),trigxx2(3*(m2-1)/2+1)
      common/fft003/ifxx3(13),trigxx3(2*(m3-1))
      common/qspe/xa(m3-1,m1-1),wor(m3-1,m1-1)
     1             ,xa2(m1-1,m3-1),wor2(m1-1,m3-1)
     1             ,xr(m2+1,m1-1),work(m2,m1-1)
      complex xa,wor,xa2,wor2
      real xr,work
      common/numond/kx(m1),ky(m2),kz(m3),kkmax
      common/dimfi/n1fi,n1mfi,n2fi,n2mfi,n3fi,n3mfi
      common/itypfi/ifilcu,ifilt
c
      n2mh=n2m/2+1
      nxcut=n1mfi
      nycut=n2mfi
      nzcut=n3mfi
      n2mhfi=n2mfi/2+1
      
      qmax=0.
      qfmax=0.

c
      do  10 k=1,n3m
c
        do i=1,n1m
         xr(1,i)=u(i,n2m,k)
         xr(n2m+2,i)=u(i,1,k)
         do j=1,n2m
          js=j+1
          xr(js,i)=u(i,j,k)
          qmax=max(abs(u(i,j,k)),qmax)
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
       emax=0.
       do kk=1,kkmax
       eq(kk)=0.0
       enddo
       do k=1,n3m
        do j=2,n2mh
         do i=1,n1m
         akk=sqrt(float(kx(i)*kx(i)+ky(j)*ky(j)+kz(k)*kz(k)))
         kk=nint(akk)
         urea=real(qtil(i,j,k))
         uimm=aimag(qtil(i,j,k))
         eq(kk)=eq(kk)+2.*(urea*urea+uimm*uimm)
         emax=max(eq(kk),emax)
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
         eq(kk)=eq(kk)+(urea*urea+uimm*uimm)
         emax=max(eq(kk),emax)
         end do
       end do
c
c     qtil is the velocity component in Fourier space.
c     Now apply a  sharp-Fourier cutoff
c     if kx>kxcut q=0.
c     if ky>kycut q=0.
c     if kz>kzcut q=0.
c
c     filter in z direction
c
       do kk=1,kkmax
       eqre(kk)=0.0
       enddo
      do k=nzcut/2+2,n3m-nzcut/2
       do j=1,n2mh
        do i=1,n1m
         akk=sqrt(float(kx(i)*kx(i)+ky(j)*ky(j)+kz(k)*kz(k)))
         kk=nint(akk)
         urea=real(qtil(i,j,k))
         uimm=aimag(qtil(i,j,k))
         if(j.gt.1) then
         eqre(kk)=eqre(kk)+2.*(urea*urea+uimm*uimm)
                    else
         eqre(kk)=eqre(kk)+(urea*urea+uimm*uimm)
                    endif
        qtil(i,j,k)=(0.,0.)
        end do
       end do
      end do
c
c     filter in x direction
c
      do k=1,n3m
       do j=1,n2mh
        do i=nxcut/2+2,n1m-nxcut/2
         akk=sqrt(float(kx(i)*kx(i)+ky(j)*ky(j)+kz(k)*kz(k)))
         kk=nint(akk)
         urea=real(qtil(i,j,k))
         uimm=aimag(qtil(i,j,k))
         if(j.gt.1) then
         eqre(kk)=eqre(kk)+2.*(urea*urea+uimm*uimm)
                    else
         eqre(kk)=eqre(kk)+(urea*urea+uimm*uimm)
                    endif
        qtil(i,j,k)=(0.,0.)
        end do
       end do
      end do
c
c     filter in y direction
c
      do k=1,n3m
       do j=nycut/2+2,n2mh
        do i=1,n1m
         akk=sqrt(float(kx(i)*kx(i)+ky(j)*ky(j)+kz(k)*kz(k)))
         kk=nint(akk)
         urea=real(qtil(i,j,k))
         uimm=aimag(qtil(i,j,k))
         eqre(kk)=eqre(kk)+2.*(urea*urea+uimm*uimm)
        qtil(i,j,k)=(0.,0.)
        end do
       end do
      end do
       do kk=1,kkmax
       eqfi(kk)=0.0
       enddo
       do k=1,n3m
        do j=2,n2mh
         do i=1,n1m
         akk=sqrt(float(kx(i)*kx(i)+ky(j)*ky(j)+kz(k)*kz(k)))
         kk=nint(akk)
         urea=real(qtil(i,j,k))
         uimm=aimag(qtil(i,j,k))
         eqfi(kk)=eqfi(kk)+2.*(urea*urea+uimm*uimm)
         emax=max(eqfi(kk),emax)
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
         eqfi(kk)=eqfi(kk)+(urea*urea+uimm*uimm)
         emax=max(eqfi(kk),emax)
         end do
       end do
c
c   back to the physical space
c
c
c   2-d  cfft applied (twice) from
c   wave number to physical space
c
       do j=1,n2mh
c
            do i=1,n1m
                  do k=1,n3m
         xa(k,i)=qtil(i,j,k)
                  enddo
            enddo
c
       call cfft99(xa,wor,trigxx3,ifxx3,1,m3-1,n3m,n1m,+1)
c
            do i=1,n1m
                  do k=1,n3m
         qtil(i,j,k)=xa(k,i)
                  enddo
            enddo
        enddo
c
       do j=1,n2mh
            do k=1,n3m
                  do i=1,n1m
         xa2(i,k)=qtil(i,j,k)
                  enddo
            enddo
       call cfft99(xa2,wor2,trigxx1,ifxx1,1,m1-1,n1m,n3m,+1)
c
            do i=1,n1m
                  do k=1,n3m
         qtil(i,j,k)=xa2(i,k)
                  enddo
            enddo
c
         enddo
c
c   2-d  fft applied from
c   wave number to physical space
c
      do k=1,n3m
c
            do j=1,n2mh
       jp=2*j
       jd=2*j-1
                  do i=1,n1m
         xr(jd,i)=real(qtil(i,j,k))
         xr(jp,i)=aimag(qtil(i,j,k))
                  enddo
            enddo
c
       call fft99(xr,work,trigxx2,ifxx2,1,m2+1,n2m,n1m,+1)
c
            do i=1,n1m
                  do j=1,n2m
         js=j+1
         uf(i,j,k)=xr(js,i)
          qfmax=max(abs(uf(i,j,k)),qfmax)
                  enddo
            enddo
c
        enddo
c
       return
       end 

c
c****************   subroutine filspa   
c   to be corrected
c
      subroutine filphy(uf)
      include 'param.f'
c     calculates filtered function using a box filter in physical
c     space. Periodic box
      common/rhsc/u(m1,m2,m3)
      dimension uf(m1,m2,m3)
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/nfilt/nr1,nr2,nr3
      nlr1=2*(nr1-1)+1
      nlr2=2*(nr2-1)+1
      nlr3=2*(nr3-1)+1
      suv3=1./float(3*nr3-2)
      suv1=1./float(3*nr1-2)
      suv2=1./float(3*nr2-2)
      write(6,*)'nr_i ',nr1,nr2,nr3,suv1,suv2,suv3
      ufc=0.
      ufk=0.
      ufj=0.
      ufi=0.
      do kc=1,n3m
      do jc=1,n2m
      do ic=1,n1m
      ufc=max(abs(u(ic,jc,kc)),ufc)
        ufsum=0.
        do lrk= 1,nlr3
        lk=lrk-1-(nlr3-1)/2
               kk=kc+lk
               if(kk.gt.n3m) kk=kk-n3m
               if(kk.lt.1) kk=kk+n3m
               weig=2.**(nr3-abs(lk))
        ufsum=ufsum+u(ic,jc,kk)*weig
        if(ic.eq.n1m/2.and.jc.eq.n2m/2.and.kc.eq.n3m/2) then
        write(6,*)u(ic,jc,kk),weig
                                                        endif
        enddo
      uf(ic,jc,kc)=suv3*ufsum
        if(ic.eq.n1m/2.and.jc.eq.n2m/2.and.kc.eq.n3m/2) then
        write(6,*)'kc direct',uf(ic,jc,kc)
                                                        endif
      ufk=max(abs(uf(ic,jc,kc)),ufk)
        enddo
        enddo
        enddo
      do kc=1,n3m
      do jc=1,n2m
      do ic=1,n1m
      u(ic,jc,kc)=uf(ic,jc,kc)
        enddo
        enddo
        enddo
      do kc=1,n3m
      do jc=1,n2m
      do ic=1,n1m
        ufsum=0.
        do lrj= 1,nlr2
        lj=lrj-1-(nlr2-1)/2
               jj=jc+lj
               if(jj.gt.n2m) jj=jj-n2m
               if(jj.lt.1) jj=jj+n2m
               weig=2.**(nr2-abs(lj))
        ufsum=ufsum+u(ic,jj,kc)*weig
        if(ic.eq.n1m/2.and.jc.eq.n2m/2.and.kc.eq.n3m/2) then
        write(6,*)u(ic,jj,kc),weig
                                                        endif
        enddo
      uf(ic,jc,kc)=suv2*ufsum
        if(ic.eq.n1m/2.and.jc.eq.n2m/2.and.kc.eq.n3m/2) then
        write(6,*)'jc direct',uf(ic,jc,kc)
                                                        endif
      ufj=max(abs(uf(ic,jc,kc)),ufj)
        enddo
        enddo
        enddo
      do kc=1,n3m
      do jc=1,n2m
      do ic=1,n1m
      u(ic,jc,kc)=uf(ic,jc,kc)
        enddo
        enddo
        enddo
      do kc=1,n3m
      do jc=1,n2m
      do ic=1,n1m
        ufsum=0.
        do lri= 1,nlr1
        li=lri-1-(nlr1-1)/2
               ii=ic+li
               if(ii.gt.n1m) ii=ii-n1m
               if(ii.lt.1) ii=ii+n1m
               weig=2.**(nr1-abs(li))
        ufsum=ufsum+u(ii,jc,kc)*weig
        if(ic.eq.n1m/2.and.jc.eq.n2m/2.and.kc.eq.n3m/2) then
        write(6,*)u(ii,jc,kc),weig
                                                        endif
        enddo
      uf(ic,jc,kc)=suv1*ufsum
        if(ic.eq.n1m/2.and.jc.eq.n2m/2.and.kc.eq.n3m/2) then
        write(6,*)'ic direct',uf(ic,jc,kc)
                                                        endif
      ufi=max(abs(uf(ic,jc,kc)),ufi)
        enddo
        enddo
        enddo
      write(6,*)'ufc,ufk,ufj,ufi ',ufc,ufk,ufj,ufi
      return
      end





      


