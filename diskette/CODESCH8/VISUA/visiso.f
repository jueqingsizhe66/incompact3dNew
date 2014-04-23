
c
      program main
      common/d1/re
      common/d2/nstop,nprint,ntst,npin,npstf
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/wrre/nwrit,nread
      common/tscoe/ga(3),ro(3),nsst
      common/d123/alx1,alx2,alx3
      common/averou/iav
      common/inior/indrea,icosma,icont
      common/ledat/ics0,cvisc,ifiltr,ibox
      common/cosma/csma
      common/cflco/icfl,cflc,tpin,tprin,tfin
      common/rot/f0
      common/ispec/imic
      common/vpeini/vper
      common/spemw/akkpp,qq,sig
      common/tstep/dt,beta,ren
      common/pardip/thet0,vsi,yc1mo,yc2mo,akmo,velmo
      common/anmod/ar1,ar2
      common/i3dpr/njumk,njumj
      common/j3dpr/mpq,jprq(15)
      pi=2.*asin(1.)
c
      open(15,file='iso.d')
      read(15,*) n1,n2,n3,nsst,alx1d,alx2d,alx3d
      read(15,*) nwrit,nread,iav
      read(15,*) re,vper,dt,ntst,nprint,npin ,npstf
      read(15,*) nstop
      read(15,*) icfl,cflc,tpin,tprin,tfin
      read(15,*) ics0,ifiltr,csma,ibox
      read(15,*) f0
c add rotation to the momentum equation
c coriolis parameter f0, rotation about z-axis i.e.  axis 3
c try writing in inverse sec.
      read(15,*) imic
      if(imic.ge.0) then
      read(15,*) akkpp,qq,sig
                    else
      read(15,*) yc1mo,yc2mo,ramo,velmo,vsi
      read(15,*) thet0
      read(15,*) ar1,ar2
      akmo=3.83711/ramo
      pi=2.*asin(1.)
      thet0=thet0*pi/180.
      yc1mo=yc1mo*pi
      yc2mo=yc2mo*pi
      write(6,201)thet0
  201 format(3x,'modone vort stream funct. tht0=',e10.4)
                    endif
      read(15,*) njumk,njumj
      read(15,*) npq,(jprq(n),n=1,npq)      
      mpq=npq
      cvisc=1./re
      n1m=n1-1
      n2m=n2-1
      n3m=n3-1
      alx1=alx1d*pi
      alx2=alx2d*pi
      alx3=alx3d*pi
      call openfi
      call visua
      stop
      end
c
      subroutine openfi
      common/names/filcnw,filcnr,filth,filou,filuu,filsf
      character*27 filcnw,filcnr,filth,filou,filuu,filsf
      open(46,file='nfisodyn')
c     continuation file
      read(46,'(a)')filcnw
c     restart file
      read(46,'(a)')filcnr
c     time history
      read(46,'(a)')filth
c     rms file
      read(46,'(a)')filou
c     skewness
      read(46,'(a)')filuu
c     flatness
      read(46,'(a)')filsf
      open(23,file=filcnr,form='unformatted')
      open(32,file=filth)
      rewind 23
      rewind 32
      return
      end
c
c  ************************************************************
c  ************subrout solve **********************************
c  ************************************************************
c
      subroutine visua
c
c     code for visualization of vorticity and velocity fields
c     for flows simulated in three periodic box
c
c     this code uses routines in the code iso
c
c     the code has a structure similar to the code
c     for postprocessing
c
      include 'param.f'
      parameter (m1m=m1-1)
c
      common/cflco/icfl,cflc,tpin,tprin,tfin
      common/wrre/nwrit,nread
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      dimension q(ndv,m1,m2,m3),pr(m1,m2,m3)
      dimension qcap(m1,m2,m3)
      common/velmax/vmax(ndv),vmaxo(ndv)
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/d1/re
      common/d123/alx1,alx2,alx3
      common/d2/nstop,nprint,ntst,npin,npstf
      common/tstep/dt,beta,ren
      common/tscoe/ga(3),ro(3),nsst
      common/inener/ene0
      common/averou/iav
      common/sc/sc,sc1
      common/names/filcnw,filcnr,filth,filou,filuu,filsf
      common/jrrr/jri,jrf,djr,irejr,iruuca
      common/ledat/ics0,cvisc,ifiltr,ibox
      common/newdat/icost,timeav
      common/cosma/csma
      common/rot/f0
      common/speene/e(ndv,0:m3)
      common/numond/kx(m1),ky(m2),kz(m3),kkmax
      common/ispec/imic
      dimension vor(m1,m2)
      character*27 filcnw,filcnr,filth,filou,filuu,filsf
      character*4 pntim
c
c
      pi=2.*asin(1.)
c
c     step  and mesh sizes calculations
c
c
      call meshes
      call indic
      call coordi
c
      write(6,754)n1,n2,n3
  754 format(10x,'centered velocities',2x
     1      ,5x,'n1=',i2,2x,'n2=',i3,2x,'n3=',i2)
      write(6,755) re,dx1,dx2,dx3,dt,ntst
  755 format(3x,'re=',e10.3,3x,'dx1=',e10.3,3x,'dx2=',e10.3,3x
     1      ,3x,'dx3=',e10.3,3x,
     1      'dt=',e10.3,3x,'ntst=',i5)
      write(6,756)  f0
  756 format(3x,'f0=',e10.3)
      ren=re
      call inirea(ntii,time,q,pr)
      print *,'leggo la soluzione'
      call divgck(q,qmax)
      print *,'divergenza della soluzione letta=',qmax
      call vmaxv(q)
      write(6,*)'  vmax',(vmax(l),l=1,3)    
      write(57,*)'  vmax',(vmax(l),l=1,3)    
      call outpf(time,enen,nav,q,pr,qcap)
      return
      end
c
c  ****************************** subrout coordi  **********************
c    this subroutine calculates the coordinates and print
c    the coordinates in 2D planes for contour plots visualizations
c    by the TURB3D package
c
      subroutine coordi
c
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/d123/alx1,alx2,alx3
      common/corpri/yp1(m1),yp2(m2),yp3(m3)
      common/cordst/y1s(m1),y3s(m3),y2s(m2)
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      common/ispec/imic
      character*60 namfile
      pi=2.*asin(1.)
c
      do 61 i=1,n1
      x1=float(i-1)/float(n1m)
      yp1(i)=-alx1*0.5+x1*alx1
   61 continue
c
      do 62 j=1,n2
      x2=float(j-1)/float(n2m)
      yp2(j)=-alx2*0.5+x2*alx2
   62 continue
c
      do 63 k=1,n3
      x3=float(k-1)/float(n3m)
      yp3(k)=-alx3*0.5+x3*alx3
   63 continue
      do  i=1,n1m
          y1s(i) = 0.5*(yp1(i) + yp1(i+1) )
      enddo
      do  k=1,n3m
          y3s(k) = 0.5*(yp3(k) + yp3(k+1) )
      enddo
      do j=1,n2m
          y2s(j)=(yp2(j)+yp2(j+1))*0.5
      enddo

      namfile='cord13.dat'
      open(18,file=namfile,form='unformatted')
      aaa=1.
c     write(6,*) 'in prico 13' ,n3m,n1m
      write(18) n1,n3,1
      write(18)
     1            ((yp1(i),i=1,n1),k=1,n3),
     1            ((yp3(k),i=1,n1),k=1,n3),
     1            ((aaa,i=1,n1),k=1,n3)
      close(18)
      namfile='cord12.dat'
      open(18,file=namfile,form='unformatted')
c     write(6,*) 'in prico 12' ,n1,n2m
      write(18) n1,n2,1
      write(18)
     1   ((yp1(i),i=1,n1),j=1,n2),
     1   ((yp2(j),i=1,n1),j=1,n2),
     1   ((aaa,i=1,n1),j=1,n2)
      close(18)
      namfile='cord32.dat'
      open(18,file=namfile,form='unformatted')
c     write(6,*) 'in prico 12' ,n1,n2m
      write(18) n3,n2,1
      write(18)
     1   ((yp3(k),k=1,n3),j=1,n2),
     1   ((yp2(j),k=1,n3),j=1,n2),
     1   ((aaa,k=1,n3),j=1,n2)
      close(18)
      return
      end
c  ****************************** subrout divgck  **********************
c
c  this subroutine perform the calculation of divg(q)
c
      subroutine divgck(vq,qmax)
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
  311 continue
      do 310 k=1,n3m
      do 310 j=1,n2m
      do 310 i=1,n1m
      vca1=q(1,i,j,k)
      vmax(1)=max(vmax(1),abs(vca1))
      vca2=q(2,i,j,k)
      vmax(2)=max(vmax(2),abs(vca2))
      vca3=q(3,i,j,k)
      vmax(3)=max(vmax(3),abs(vca3))
  310 continue
      return
      end
c
c  ****************************** subrout inirea **********************
c
c   reading initial  conditions from file
c   this file is different from the restarting file.
c   It must be produced by running only one time step
c   of the ISO code with a very small Delta t.
c   The user should insert in the code ISO the routine
c   that print the field in the following format.
c   ************   IMPORTANT ******************************
c
      subroutine inirea(ntii,time,q,pr)
c
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      dimension pr(m1,m2,m3)
      dimension q(ndv,m1,m2,m3)
      common/velmax/vmax(ndv),vmaxo(ndv)
      common/cordst/y1s(m1),y3s(m3),y2s(m2)
      common/tstep/dt,beta,ren
      common/wallst/cflw,cfuw
      common/veltot/vit(ndv)
      common/omegas/omr,omi
      common/timref/tref
      common/inener/ene0
      common/eneav/enav,enavo
      common/inior/indrea,icosma,icont
      common/d123/alx1,alx2,alx3
      common/rot/f0
c
      nfil=23
      rewind(nfil)
      read(nfil) n1,n2,n3
      read(nfil) time,ren,time,dum
      read(nfil) (((q(1,i,j,k),i=1,n1),j=1,n2),k=1,n3),
     1            (((q(2,i,j,k),i=1,n1),j=1,n2),k=1,n3),
     1            (((q(3,i,j,k),i=1,n1),j=1,n2),k=1,n3)
     1           ,(((pr(i,j,k),i=1,n1),j=1,n2),k=1,n3)
      close(nfil)
      pmean=0.
      do i=1,n1m
      do j=1,n2m
      do k=1,n3m
      pmean=pmean+pr(i,j,k)
      enddo
      enddo
      enddo
      pmean=pmean/float(n1m*n2m*n3m)
      do i=1,n1m
      do j=1,n2m
      do k=1,n3m
      pr(i,j,k)=pr(i,j,k)-pmean
      enddo
      enddo
      enddo
      write(6,*)' in inirea pmean=',pmean
      return
      end
c
c  ****************************** subrout outpf  **********************
c
      subroutine outpf(time,enen,nav,q,pr,qcap)
c
      include 'param.f'
      parameter (m1m=m1-1,m2m=m2-1)
      common/dim/n1,n1m,n2,n2m,n3,n3m
      dimension q1sh(m1,m2),q2sh(m1,m2),q3sh(m1,m2)
      dimension vo1sh(m1,m2),vo2sh(m1,m2),vo3sh(m1,m2)
      dimension prsh(m1,m2)
      dimension q1se(m1,m2),q2se(m1,m2),q3se(m1,m2)
      dimension vo1se(m1,m2),vo2se(m1,m2),vo3se(m1,m2)
      dimension prse(m1,m2)
      dimension pr(m1,m2,m3)
      dimension qcap(m1,m2,m3)
      dimension q(ndv,m1,m2,m3)
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      common/tstep/dt,beta,ren
      dimension vprms(4),vmp(3),skp(4),flp(4)
      common/filep/ifilp
      common/vmean/vm(ndv),vrms(4)
      common/eneav/enav,enavo
      common/qmean/qm(ndv),qrms(4)
      common/skfl/ske(4),fla(4)
      common/prrm/prm,prms
      common/skflo/skeo(4),flao(4)
      common/prrmo/prmo,prmso
      common/vmeao/vmo(ndv),vrmso(4)
      common/names/filcnw,filcnr,filth,filou,filuu,filsf
      common/outd1/cdyn(m2),ell1(m2),alij(m2),amij(m2)
      common/outd2/cs(m2),flmij(m2),fmmij(m2)
      common/corpri/yp1(m1),yp2(m2),yp3(m3)
      common/speene/e(ndv,0:m3)
      common/numond/kx(m1),ky(m2),kz(m3),kkmax
      common/ledat/ics0,cvisc,ifiltr,ibox
      common/kolqua/diss,eta
      common/ispec/imic
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/i3dpr/njumk,njumj
      common/j3dpr/mpq,jprq(15)
c      
      character*27 filcnw,filcnr,filth,filou,filuu,filsf
      character*17 filcos
      character*4 pntim
      character*6 slp
      character*1 is1s,is3s
      character*1 is1n,is3n
      character*60 namfi3
      character*3 nkpse,njpse
      write(6,*)'enter 1 to shift the xy plot'
      read(5,*) ishift
      itim=time+0.3
      write(pntim,77) itim
   77 format(i4.4)
      nav=1
      if(imic.gt.0) then
      call enerca(q,pr,enen,time)
      call pdiss(q,enen,diss,flamb,eta)
c     print*,'in outpf'
c
      enavp=enavo/nav
      filcos='spec.'//pntim
      open(57,file=filcos)
      rewind(57)
      filcos='speko.'//pntim
      open(59,file=filcos)
      rewind(59)
      do 615 l=1,3
      skp(l)=skeo(l)/nav
      flp(l)=flao(l)/nav
      vmp(l)=vmo(l)/nav
      qrms(l)=sqrt(qrms(l))
  615 vprms(l)=sqrt(vrmso(l)/nav)
      skp(4)=skeo(4)/nav
      flp(4)=flao(4)/nav
      vprms(4)=vrmso(4)/nav
      ppmp=prmso
  611 continue
  612 format(1x,i3,2x,e14.7,2x,9(1x,e11.4)) 
  613 format(e14.7,1x,9(1x,e11.4)) 
      close(42)
      close(62)
c
c    spectra calculation
c
      call spectre(q,qcap)
      speska=(cvisc**5./diss)**(1./4.)
       ee1t0=0.
       ee2t0=0.
       ee3t0=0.
       do k=1,kkmax
       write (57,*)k,e(1,k),e(2,k),e(3,k)
       ee1to=ee1to+e(1,k)
       ee2to=ee2to+e(2,k)
       ee3to=ee3to+e(3,k)
       akol=k*eta
       se1=e(1,k)/speska
       se2=e(2,k)/speska
       se3=e(3,k)/speska
       write (59,*)akol,se1,se2,se3          
       end do
      close(57)
      close(59)
      write(6,133)time,ee1to,ee2to,ee3to,eeto,enen
  133 format('energy',e12.4,3x,3e12.4,3x,'etot sp,phy',2e12.5)
             endif
c
c  THE QUANTITIES ARE EVALUATED AT 
C        THE POINT I,J,K (
c
c  1-2  section at each k position 
c
   83 format(i3.3)
      ki=1+njumk/2
      do k=ki,n3m,njumk
      write(nkpse,83)k
      namfi3='se12pl'//nkpse//'.dat'
c     write(6,*)' written sections aty3=',y3s(k)
      kp=kpv(k)
      km=kmv(k)
      do i=1,n1m
      ip=ipv(i)
      im=imv(i)
      do j=1,n2m
      jp=jpv(j)
      jm=jmv(j)
      q2c=0.25*
     1  (    q(2,i,j,k)+q(2,i,j,km)
     1     + q(2,im,j,k)+q(2,im,j,km))
      q1c=0.25*
     1  (    q(1,i,j,k)+q(1,i,jm,k) 
     1     + q(1,i,j,km)+q(1,i,jm,km))
      q3c=0.25*
     1  (    q(3,i,j,k)+q(3,im,j,k) 
     1     + q(3,i,jm,k)+q(3,im,jm,k))
c
c      OM_1 COMPONENT
c
      om1=(+(q(3,i,j,k)-q(3,i,jm,k))*dx2
     1    -(q(2,i,j,k)-q(2,i,j,km))*dx3
     1    +(q(3,im,j,k)-q(3,im,jm,k))*dx2
     1    -(q(2,im,j,k)-q(2,im,j,km))*dx3)*0.5
c
c      OM_y COMPONENT
c
      om2=(+(q(1,i,j,k)-q(1,i,j,km))*dx3
     1    -(q(3,i,j,k)-q(3,im,j,k))*dx1
     1    +(q(1,i,jm,k)-q(1,i,jm,km))*dx3
     1    -(q(3,i,jm,k)-q(3,im,jm,k))*dx1)*0.5
c
c      OM_z COMPONENT
c
      om3=(-(q(1,i,j,k)-q(1,i,jm,k))*dx2
     1    +(q(2,i,j,k)-q(2,im,j,k))*dx1
     1    -(q(1,i,j,km)-q(1,i,jm,km))*dx2
     1    +(q(2,i,j,km)-q(2,im,j,km))*dx1)*0.5
      q1se(i,j)=q1c
      q2se(i,j)=q2c
      q3se(i,j)=q3c
      prse(i,j)=(pr(i,j,k)+pr(im,j,k)+pr(im,j,km)+pr(i,j,km)
     1          +pr(i,jm,k)+pr(im,jm,k)+pr(im,jm,km)+pr(im,j,km))*0.125
      vo1se(i,j)=om1
      vo2se(i,j)=om2
      vo3se(i,j)=om3
      enddo
      enddo
      i=n1
      do j=1,n2m
      q1se(i,j)=q1se(1,j)
      q2se(i,j)=q2se(1,j)
      q3se(i,j)=q3se(1,j)
      prse(i,j)=prse(1,j)
      vo1se(i,j)=vo1se(1,j)
      vo2se(i,j)=vo2se(1,j)
      vo3se(i,j)=vo3se(1,j)
      enddo
      j=n2
      do i=1,n1
      q1se(i,j)=q1se(i,1)
      q2se(i,j)=q2se(i,1)
      q3se(i,j)=q3se(i,1)
      prse(i,j)=prse(i,1)
      vo1se(i,j)=vo1se(i,1)
      vo2se(i,j)=vo2se(i,1)
      vo3se(i,j)=vo3se(i,1)
      enddo
      if(ishift.eq.0) then
      open(59,file=namfi3,form='unformatted')
      rewind 59
      write(59) n1,n2,1
      write(59) re,re,re,time
      write(59)
     1   ((q1se(i,j),i=1,n1),j=1,n2),
     1   ((q2se(i,j),i=1,n1),j=1,n2),
     1   ((q3se(i,j),i=1,n1),j=1,n2),
     1   ((prse(i,j),i=1,n1),j=1,n2),
     1   ((vo1se(i,j),i=1,n1),j=1,n2),
     1   ((vo2se(i,j),i=1,n1),j=1,n2),
     1   ((vo3se(i,j),i=1,n1),j=1,n2)
                      else
      n1mh=n1m/2+1
      do j=1,n2
      do i=1,n1m
      if(i.le.n1mh) then
      ii=i+n1m/2
                     else
      ii=i-n1m/2
                     endif
      q1sh(ii,j)=q1se(i,j)
      q2sh(ii,j)=q2se(i,j)
      q3sh(ii,j)=q3se(i,j)
      prsh(ii,j)=prse(i,j)
      vo1sh(ii,j)=vo1se(i,j)
      vo2sh(ii,j)=vo2se(i,j)
      vo3sh(ii,j)=vo3se(i,j)
      enddo
      enddo
      open(59,file=namfi3,form='unformatted')
      rewind 59
      write(59) n1,n2,1
      write(59) re,re,re,time
      write(59)
     1   ((q1sh(i,j),i=1,n1),j=1,n2),
     1   ((q2sh(i,j),i=1,n1),j=1,n2),
     1   ((q3sh(i,j),i=1,n1),j=1,n2),
     1   ((prsh(i,j),i=1,n1),j=1,n2),
     1   ((vo1sh(i,j),i=1,n1),j=1,n2),
     1   ((vo2sh(i,j),i=1,n1),j=1,n2),
     1   ((vo3sh(i,j),i=1,n1),j=1,n2)
                      endif
      close(59)
      enddo
c
c  2-3  section
c
      do i=1,n1m,njumk
      write(nkpse,83)i
      namfi3='se32pl'//nkpse//'.dat'
c     write(6,*)' written sections aty3=',y3s(k)
      ip=ipv(i)
      im=imv(i)
      do k=1,n3m
      kp=kpv(k)
      km=kmv(k)
      do j=1,n2m
      jp=jpv(j)
      jm=jmv(j)
      q2c=0.25*
     1  (    q(2,i,j,k)+q(2,i,j,km)
     1     + q(2,im,j,k)+q(2,im,j,km))
      q1c=0.25*
     1  (    q(1,i,j,k)+q(1,i,jm,k) 
     1     + q(1,i,j,km)+q(1,i,jm,km))
      q3c=0.25*
     1  (    q(3,i,j,k)+q(3,im,j,k) 
     1     + q(3,i,jm,k)+q(3,im,jm,k))
c
c      OM_1 COMPONENT
c
      om1=(+(q(3,i,j,k)-q(3,i,jm,k))*dx2
     1    -(q(2,i,j,k)-q(2,i,j,km))*dx3
     1    +(q(3,im,j,k)-q(3,im,jm,k))*dx2
     1    -(q(2,im,j,k)-q(2,im,j,km))*dx3)*0.5
c
c      OM_y COMPONENT
c
      om2=(+(q(1,i,j,k)-q(1,i,j,km))*dx3
     1    -(q(3,i,j,k)-q(3,im,j,k))*dx1
     1    +(q(1,i,jm,k)-q(1,i,jm,km))*dx3
     1    -(q(3,i,jm,k)-q(3,im,jm,k))*dx1)*0.5
c
c      OM_z COMPONENT
c
      om3=(-(q(1,i,j,k)-q(1,i,jm,k))*dx2
     1    +(q(2,i,j,k)-q(2,im,j,k))*dx1
     1    -(q(1,i,j,km)-q(1,i,jm,km))*dx2
     1    +(q(2,i,j,km)-q(2,im,j,km))*dx1)*0.5
      q1se(k,j)=q1c
      q2se(k,j)=q2c
      q3se(k,j)=q3c
      prse(k,j)=(pr(i,j,k)+pr(im,j,k)+pr(im,j,km)+pr(i,j,km)
     1          +pr(i,jm,k)+pr(im,jm,k)+pr(im,jm,km)+pr(im,j,km))*0.125
      vo1se(k,j)=om1
      vo2se(k,j)=om2
      vo3se(k,j)=om3
      enddo
      enddo
      k=n3
      do j=1,n2m
      q1se(k,j)=q1se(1,j)
      q2se(k,j)=q2se(1,j)
      q3se(k,j)=q3se(1,j)
      prse(k,j)=prse(1,j)
      vo1se(k,j)=vo1se(1,j)
      vo2se(k,j)=vo2se(1,j)
      vo3se(k,j)=vo3se(1,j)
      enddo
      j=n2
      do k=1,n3
      q1se(k,j)=q1se(k,1)
      q2se(k,j)=q2se(k,1)
      q3se(k,j)=q3se(k,1)
      prse(k,j)=prse(k,1)
      vo1se(k,j)=vo1se(k,1)
      vo2se(k,j)=vo2se(k,1)
      vo3se(k,j)=vo3se(k,1)
      enddo
      open(59,file=namfi3,form='unformatted')
      rewind 59
      write(59) n3,n2,1
      write(59) re,re,re,time
      write(59)
     1   ((q1se(k,j),k=1,n3),j=1,n2),
     1   ((q2se(k,j),k=1,n3),j=1,n2),
     1   ((q3se(k,j),k=1,n3),j=1,n2),
     1   ((prse(k,j),k=1,n3),j=1,n2),
     1   ((vo1se(k,j),k=1,n3),j=1,n2),
     1   ((vo2se(k,j),k=1,n3),j=1,n2),
     1   ((vo3se(k,j),k=1,n3),j=1,n2)
      close(59)
      enddo
c
c  1-3  section
c
      npq=1
      if(mpq.eq.1) then
      ji=1+njumj/2
      jf=n2m
      nju=njumj
      do j=ji,jf,nju    
      write(nkpse,83)j
      namfi3='se13pl'//nkpse//'.dat'
c     write(6,*)' written sections aty3=',y3s(k)
      jp=jpv(j)
      jm=jmv(j)
      do i=1,n1m
      ip=ipv(i)
      im=imv(i)
      do k=1,n3m
      kp=kpv(k)
      km=kmv(k)
      q2c=0.25*
     1  (    q(2,i,j,k)+q(2,i,j,km)
     1     + q(2,im,j,k)+q(2,im,j,km))
      q1c=0.25*
     1  (    q(1,i,j,k)+q(1,i,jm,k) 
     1     + q(1,i,j,km)+q(1,i,jm,km))
      q3c=0.25*
     1  (    q(3,i,j,k)+q(3,im,j,k) 
     1     + q(3,i,jm,k)+q(3,im,jm,k))
c
c      OM_1 COMPONENT
c
      om1=(+(q(3,i,j,k)-q(3,i,jm,k))*dx2
     1    -(q(2,i,j,k)-q(2,i,j,km))*dx3
     1    +(q(3,im,j,k)-q(3,im,jm,k))*dx2
     1    -(q(2,im,j,k)-q(2,im,j,km))*dx3)*0.5
c
c      OM_y COMPONENT
c
      om2=(+(q(1,i,j,k)-q(1,i,j,km))*dx3
     1    -(q(3,i,j,k)-q(3,im,j,k))*dx1
     1    +(q(1,i,jm,k)-q(1,i,jm,km))*dx3
     1    -(q(3,i,jm,k)-q(3,im,jm,k))*dx1)*0.5
c
c      OM_z COMPONENT
c
      om3=(-(q(1,i,j,k)-q(1,i,jm,k))*dx2
     1    +(q(2,i,j,k)-q(2,im,j,k))*dx1
     1    -(q(1,i,j,km)-q(1,i,jm,km))*dx2
     1    +(q(2,i,j,km)-q(2,im,j,km))*dx1)*0.5
      q1se(i,k)=q1c
      q2se(i,k)=q2c
      q3se(i,k)=q3c
      prse(i,k)=(pr(i,j,k)+pr(im,j,k)+pr(im,j,km)+pr(i,j,km)
     1          +pr(i,jm,k)+pr(im,jm,k)+pr(im,jm,km)+pr(im,j,km))*0.125
      vo1se(i,k)=om1
      vo2se(i,k)=om2
      vo3se(i,k)=om3
      enddo
      enddo
      i=n1
      do k=1,n3m
      q1se(i,k)=q1se(1,k)
      q2se(i,k)=q2se(1,k)
      q3se(i,k)=q3se(1,k)
      prse(i,k)=prse(1,k)
      vo1se(i,k)=vo1se(1,k)
      vo2se(i,k)=vo2se(1,k)
      vo3se(i,k)=vo3se(1,k)
      enddo
      k=n3
      do i=1,n1
      q1se(i,k)=q1se(i,1)
      q2se(i,k)=q2se(i,1)
      q3se(i,k)=q3se(i,1)
      prse(i,k)=prse(i,1)
      vo1se(i,k)=vo1se(i,1)
      vo2se(i,k)=vo2se(i,1)
      vo3se(i,k)=vo3se(i,1)
      enddo
      open(59,file=namfi3,form='unformatted')
      rewind 59
      write(59) n1,n3,1
      write(59) re,re,re,time
      write(59)
     1   ((q1se(i,k),i=1,n1),k=1,n3),
     1   ((q2se(i,k),i=1,n1),k=1,n3),
     1   ((q3se(i,k),i=1,n1),k=1,n3),
     1   ((prse(i,k),i=1,n1),k=1,n3),
     1   ((vo1se(i,k),i=1,n1),k=1,n3),
     1   ((vo2se(i,k),i=1,n1),k=1,n3),
     1   ((vo3se(i,k),i=1,n1),k=1,n3)
      enddo
                   else
      ji=n2m
      jf=1
      nju=-1
      do j=ji,jf,nju    
      if(j.eq.jprq(npq)) then
      write(nkpse,83)j
      namfi3='se13pl'//nkpse//'.dat'
c     write(6,*)' written sections aty3=',y3s(k)
      jp=jpv(j)
      jm=jmv(j)
      do i=1,n1m
      ip=ipv(i)
      im=imv(i)
      do k=1,n3m
      kp=kpv(k)
      km=kmv(k)
      q2c=0.25*
     1  (    q(2,i,j,k)+q(2,i,j,km)
     1     + q(2,im,j,k)+q(2,im,j,km))
      q1c=0.25*
     1  (    q(1,i,j,k)+q(1,i,jm,k) 
     1     + q(1,i,j,km)+q(1,i,jm,km))
      q3c=0.25*
     1  (    q(3,i,j,k)+q(3,im,j,k) 
     1     + q(3,i,jm,k)+q(3,im,jm,k))
c
c      OM_1 COMPONENT
c
      om1=(+(q(3,i,j,k)-q(3,i,jm,k))*dx2
     1    -(q(2,i,j,k)-q(2,i,j,km))*dx3
     1    +(q(3,im,j,k)-q(3,im,jm,k))*dx2
     1    -(q(2,im,j,k)-q(2,im,j,km))*dx3)*0.5
c
c      OM_y COMPONENT
c
      om2=(+(q(1,i,j,k)-q(1,i,j,km))*dx3
     1    -(q(3,i,j,k)-q(3,im,j,k))*dx1
     1    +(q(1,i,jm,k)-q(1,i,jm,km))*dx3
     1    -(q(3,i,jm,k)-q(3,im,jm,k))*dx1)*0.5
c
c      OM_z COMPONENT
c
      om3=(-(q(1,i,j,k)-q(1,i,jm,k))*dx2
     1    +(q(2,i,j,k)-q(2,im,j,k))*dx1
     1    -(q(1,i,j,km)-q(1,i,jm,km))*dx2
     1    +(q(2,i,j,km)-q(2,im,j,km))*dx1)*0.5
      q1se(i,k)=q1c
      q2se(i,k)=q2c
      q3se(i,k)=q3c
      prse(i,k)=(pr(i,j,k)+pr(im,j,k)+pr(im,j,km)+pr(i,j,km)
     1          +pr(i,jm,k)+pr(im,jm,k)+pr(im,jm,km)+pr(im,j,km))*0.125
      vo1se(i,k)=om1
      vo2se(i,k)=om2
      vo3se(i,k)=om3
      enddo
      enddo
      i=n1
      do k=1,n3m
      q1se(i,k)=q1se(1,k)
      q2se(i,k)=q2se(1,k)
      q3se(i,k)=q3se(1,k)
      prse(i,k)=prse(1,k)
      vo1se(i,k)=vo1se(1,k)
      vo2se(i,k)=vo2se(1,k)
      vo3se(i,k)=vo3se(1,k)
      enddo
      k=n3
      do i=1,n1
      q1se(i,k)=q1se(i,1)
      q2se(i,k)=q2se(i,1)
      q3se(i,k)=q3se(i,1)
      prse(i,k)=prse(i,1)
      vo1se(i,k)=vo1se(i,1)
      vo2se(i,k)=vo2se(i,1)
      vo3se(i,k)=vo3se(i,1)
      enddo
      open(59,file=namfi3,form='unformatted')
      rewind 59
      write(59) n1,n3,1
      write(59) re,re,re,time
      write(59)
     1   ((q1se(i,k),i=1,n1),k=1,n3),
     1   ((q2se(i,k),i=1,n1),k=1,n3),
     1   ((q3se(i,k),i=1,n1),k=1,n3),
     1   ((prse(i,k),i=1,n1),k=1,n3),
     1   ((vo1se(i,k),i=1,n1),k=1,n3),
     1   ((vo2se(i,k),i=1,n1),k=1,n3),
     1   ((vo3se(i,k),i=1,n1),k=1,n3)
      close(59)
      npq=npq+1
                   endif
      enddo
                   endif
      return
      end
c
c  ****************************** subrout enerca **********************
c
c   calculation of statics in one point 
c
      subroutine enerca(q,pr,enej,time)
c
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      dimension q(ndv,m1,m2,m3),pr(m1,m2,m3)
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      common/meshu/d1x1,d1x2,d1x3
      common/qmean/qm(ndv),qrms(4)
      common/d1/re
      common/d123/alx1,alx2,alx3
      common/veltot/vit(ndv)
      common/eneav/enav,enavo
      common/vmean/vm(ndv),vrms(4)
      common/skfl/ske(4),fla(4)
      common/prrm/prm,prms
      common/skflo/skeo(4),flao(4)
      common/prrmo/prmo,prmso
      common/skder/skd1,skd2,skd3,fld1,fld2,fld3
      common/ledat/ics0,cvisc,ifiltr,ibox
c
      pi=2.*asin(1.)
      vol=(2.*pi)**3.
      vl13=1./float(n1m*n3m)
      vl123=1./float(n1m*n2m*n3m)
c     vit(1)=0.
c     vit(2)=0.
c     vit(3)=0.
      vm1m=0.
      vm2m=0.
      vm3m=0.
      prm=0.
      do 411 k=1,n3m
      kp=kpv(k)
      do 411 j=1,n2m
      jp=jpv(j)
      do 411 i=1,n1m
      ip=ipv(i)
       vm1m=q(1,i,j,k)+vm1m
       vm2m=q(2,i,j,k)+vm2m
       vm3m=q(3,i,j,k)+vm3m
      prm=prm+pr(i,j,k)
c     vit(1)=q(1,i,j,k)+vit(1)
c     vit(2)=q(2,i,j,k)+vit(2)
c     vit(3)=q(3,i,j,k)+vit(3)
  411 continue
      vm(1)=vm1m*vl123
      vm(2)=vm2m*vl123
      vm(3)=vm3m*vl123
      prm=prm*vl123
  410 continue
c     vit(1)=vit(1)*vl123/vol
c     vit(2)=vit(2)*vl123/vol
c     vit(3)=vit(3)*vl123/vol
      enej=0.
      u1rmm=0.
      u2rmm=0.
      u3rmm=0.
      u23mm=0.
      ppmm=0.
      do 414 k=1,n3m
      kp=kpv(k)
      do 414 j=1,n2m
      jp=jpv(j)
      do 414 i=1,n1m
      ip=ipv(i)
       v1m=q(1,i,j,k)-vm(1)
       v2m=q(2,i,j,k)-vm(2)
       v3m=q(3,i,j,k)-vm(3)
      ppm=pr(i,j,k)-prm
      u1rmm=u1rmm+v1m**2
      u2rmm=u2rmm+v2m**2
      u3rmm=u3rmm+v3m**2
      u23mm=u23mm+v2m*v3m
      ppmm=ppmm+ppm**2
      enejp=(v1m**2+v2m**2+v3m**2)*vl123
      enej=enej+enejp
  414 continue
      vrms(1)=u1rmm*vl123
      vrms(2)=u2rmm*vl123
      vrms(3)=u3rmm*vl123
      vrms(4)=u23mm*vl123
      prms=ppmm*vl123
c
c    skewness and flatness
c
      sk1m=0.
      sk2m=0.
      sk3m=0.
      sk4m=0.
      fl1m=0.
      fl2m=0.
      fl3m=0.
      fl4m=0.
      do 514 k=1,n3m
      kp=kpv(k)
      do 514 j=1,n2m
      jp=jpv(j)
      do 514 i=1,n1m
      ip=ipv(i)
      v1m=(q(1,ip,j,k)+q(1,i,j,k))*0.5-vm(1)
      v2m=(q(2,i,jp,k)+q(2,i,j,k))*0.5-vm(2)
      v3m=(q(3,i,j,kp)+q(3,i,j,k))*0.5-vm(3)
      ppm=pr(i,j,k)-prm
      sk1m=sk1m+(v1m/sqrt(vrms(1)))**3
      sk2m=sk2m+(v2m/sqrt(vrms(2)))**3
      sk3m=sk3m+(v3m/sqrt(vrms(3)))**3
      sk4m=sk4m+(ppm/sqrt(prms))**3
      fl1m=fl1m+(v1m/sqrt(vrms(1)))**4
      fl2m=fl2m+(v2m/sqrt(vrms(2)))**4
      fl3m=fl3m+(v3m/sqrt(vrms(3)))**4
      fl4m=fl4m+(ppm/sqrt(prms))**4
  514 continue
      ske(1)=sk1m*vl123
      ske(2)=sk2m*vl123
      ske(3)=sk3m*vl123
      ske(4)=sk4m*vl123
      fla(1)=fl1m*vl123
      fla(2)=fl2m*vl123
      fla(3)=fl3m*vl123
      fla(4)=fl4m*vl123
  513 continue
      enav=enej
      do 428 l=1,3
      qrms(l)=vrms(l)
      qm(l)=vm(l)
  428 continue
      qrms(4)=vrms(4)
      write(6,*) 'vm=',vm(1),vm(2),vm(3),prm
      write(6,*) 'vrms=',vrms(1),vrms(2),vrms(3),vrms(4)
      write(6,*) 'ske=',ske(1),ske(2),ske(3),ske(4)
      write(6,*) 'fla=',fla(1),fla(2),fla(3),fla(4)
c
c     skewness of velocity derivative
c
      rmsdux=0.
      rmsduy=0.
      rmsduz=0.
      do 714 k=1,n3m
      kp=kpv(k)
      do 714 j=1,n2m
      jp=jpv(j)
      do 714 i=1,n1m
      ip=ipv(i)
      dux=(q(1,ip,j,k)-q(1,i,j,k))*dx1
      duy=(q(2,i,jp,k)-q(2,i,j,k))*dx2
      duz=(q(3,i,j,kp)-q(3,i,j,k))*dx3
      rmsdux=rmsdux+dux*dux
      rmsduy=rmsduy+duy*duy
      rmsduz=rmsduz+duz*duz
 714  continue
      rmsdux=rmsdux*vl123
      rmsduy=rmsduy*vl123
      rmsduz=rmsduz*vl123
c     write(*,*) 'rmsder=',rmsdux,rmsduy,rmsduz
c
      skd1=0.
      skd2=0.
      skd3=0.
      fld1=0.
      fld2=0.
      fld3=0.
      do 724 k=1,n3m
      kp=kpv(k)
      do 724 j=1,n2m
      jp=jpv(j)
      do 724 i=1,n1m
      ip=ipv(i)
      dux=(q(1,ip,j,k)-q(1,i,j,k))*dx1
      duy=(q(2,i,jp,k)-q(2,i,j,k))*dx2
      duz=(q(3,i,j,kp)-q(3,i,j,k))*dx3
      skd1=skd1+(dux/sqrt(rmsdux))**3
      skd2=skd2+(duy/sqrt(rmsduy))**3
      skd3=skd3+(duz/sqrt(rmsduz))**3
      fld1=fld1+(dux/sqrt(rmsdux))**4
      fld2=fld2+(duy/sqrt(rmsduy))**4
      fld3=fld3+(duz/sqrt(rmsduz))**4
 724  continue
      skd1=skd1*vl123
      skd2=skd2*vl123
      skd3=skd3*vl123
      fld1=fld1*vl123
      fld2=fld2*vl123
      fld3=fld3*vl123
      write(*,*) 'skeder=',skd1,skd2,skd3
      write(*,*) 'flader=',fld1,fld2,fld3
      return
      end

c************************************************************
c****************  pdiss   **********************************
c************************************************************
c calculate the dissipation. average over strain at 4 points 
c
      subroutine pdiss(q,enej,diss,flamb,eta)
      include 'param.f'
      dimension q(ndv,m1,m2,m3)
      dimension ell2(m2),ell1(m2)
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      common/d123/alx1,alx2,alx3
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/ledat/ics0,cvisc,ifiltr,ibox
      common/d2/nstop,nprint,ntst,npin,npstf
      common/cosma/csma
      dimension visout(m2)
      character*20 titfil
      character*3 icount
c      
c
c      computation of Sij
c
c     call strper
c                                                                       
c     this calculates the  strain rate tensor at the center of
c     the cell. Perodic box
c                                                                       
c     
c     print*,'in pdiss'
      vl123=1./float(n1m*n2m*n3m)
      diss=0.
      diss1=0.
      do 4 kc=1,n3m 
      kp=kpv(kc)                                                    
      km=kmv(kc)                                                    
      do 4 jc=1,n2m                                                     
      jp=jpv(jc)
      jm=jmv(jc)
      do 4 ic=1,n1m                                                     
      ip=ipv(ic)                                                        
      im=imv(ic)                                                        
c
c   stl  =Sll  at i+1/2,j+1/2,k+1/2
c   st2  =S22  at i+1/2,j+1/2,k+1/2
c   st3  =S33  at i+1/2,j+1/2,k+1/2
c
      st1=(q(1,ip,jc,kc)-q(1,ic,jc,kc))*dx1
      st2=(q(2,ic,jp,kc)-q(2,ic,jc,kc))*dx2
      st3=(q(3,ic,jc,kp)-q(3,ic,jc,kc))*dx3
c
c   st4 = S12    at  i,j,k+1/2
c   st6 = S32    at  i+1/2,j,k
c   st5 = S13    at  i,j+1/2,k
c
c
      stccc4=(q(1,ic,jc,kc)-q(1,ic,jm,kc))*dx2
     1               +(q(2,ic,jc,kc)-q(2,im,jc,kc))*dx1
      stccc5=(q(1,ic,jc,kc)-q(1,ic,jc,km))*dx3
     1              +(q(3,ic,jc,kc)-q(3,im,jc,kc))*dx1
      stccc6=(q(3,ic,jc,kc)-q(3,ic,jm,kc))*dx2
     1               +(q(2,ic,jc,kc)-q(2,ic,jc,km))*dx3
c
      stcpc4=(q(1,ic,jp,kc)-q(1,ic,jc,kc))*dx2
     1               +(q(2,ic,jp,kc)-q(2,im,jp,kc))*dx1
      stccp5=(q(1,ic,jc,kp)-q(1,ic,jc,kc))*dx3
     1              +(q(3,ic,jc,kp)-q(3,im,jc,kp))*dx1
      stcpc6=(q(3,ic,jp,kc)-q(3,ic,jc,kc))*dx2
     1               +(q(2,ic,jp,kc)-q(2,ic,jp,km))*dx3
c
      stpcc4=(q(1,ip,jc,kc)-q(1,ip,jm,kc))*dx2
     1               +(q(2,ip,jc,kc)-q(2,ic,jc,kc))*dx1
      stpcc5=(q(1,ip,jc,kc)-q(1,ip,jc,km))*dx3
     1              +(q(3,ip,jc,kc)-q(3,ic,jc,kc))*dx1
      stccp6=(q(3,ic,jc,kp)-q(3,ic,jm,kp))*dx2
     1               +(q(2,ic,jc,kp)-q(2,ic,jc,kc))*dx3
c
      stppc4=(q(1,ip,jp,kc)-q(1,ip,jc,kc))*dx2
     1               +(q(2,ip,jp,kc)-q(2,ic,jp,kc))*dx1
      stpcp5=(q(1,ip,jc,kp)-q(1,ip,jc,kc))*dx3
     1              +(q(3,ip,jc,kp)-q(3,ic,jc,kp))*dx1
      stcpp6=(q(3,ic,jp,kp)-q(3,ic,jc,kp))*dx2
     1               +(q(2,ic,jp,kp)-q(2,ic,jp,kc))*dx3
c
c   st4 = S12    at  center
c   st5 = S13    at  center
c   st6 = S32    at  center
c
      st4=.5*0.25*(stccc4+stcpc4+stpcc4+stppc4)
      st5=.5*0.25*(stccc5+stccp5+stpcc5+stpcp5)
      st6=.5*0.25*(stccc6+stcpc6+stccp6+stcpp6)
        app1=.5*(st1*st1+st2*st2+ st3*st3+
     1   2.*(st4*st4+st5*st5+st6*st6))*vl123*8.
        diss=diss+app1*cvisc
    4 continue                                                          
c      
c     print *,'nu*diss=',cvisc*diss
      eta=(cvisc**3./diss)**(1./4.)
      flamb=sqrt(20*cvisc*enej/diss)
 
c     print*,' before retun in pdiss'
      return
      end
c
c    calculation of the energy spectrum
c
      subroutine spectre(q,qtil)
c     energy spectrum 
      include 'param.f'
      dimension q(ndv,m1,m2,m3)
      complex qtil(m1,(m2-1)/2+1,m3)
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/fft001/ifxx1(13),trigxx1(2*(m1-1))
      common/fft002/ifxx2(13),trigxx2(3*(m2-1)/2+1)
      common/fft003/ifxx3(13),trigxx3(2*(m3-1))
      common/numond/kx(m1),ky(m2),kz(m3),kkmax
      common/speene/e(ndv,0:m3)
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
       do kk=1,kkmax
       e(n,kk)=0.0
       enddo
       do k=1,n3m
        do j=2,n2mh
         do i=1,n1m
         akk=sqrt(float(kx(i)*kx(i)+ky(j)*ky(j)+kz(k)*kz(k)))
         kk=nint(akk)
         urea=real(qtil(i,j,k))
         uimm=aimag(qtil(i,j,k))
         e(n,kk)=e(n,kk)+2.*(urea*urea+uimm*uimm)
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
         e(n,kk)=e(n,kk)+(urea*urea+uimm*uimm)
         end do
       end do
  1    continue
       return
       end 
