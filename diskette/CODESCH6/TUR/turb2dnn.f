c
c   ****** subroutine vorqua  ****************
c
c     computes max,min values and corresponding coordinates
c  this subroutine is not so important for 2D turbulence
c
      subroutine vorqua(ntime,vor,psi,time,vormax,vormin)
      include 'param.f'
      dimension vor(m1,m2),psi(m1,m2)
      common/coor/yp1(m1),yp2(m2)
      common/dim/n1,n1m,n2,n2m
      common/mesh/dx1,dx1q,dx2,dx2q
      vorp=0.
      vorm=0.
      vormax=-10000.
      vormin=10000.
      do 328 j=1,n2m
      do 328 i=1,n1m
      vormax=max(vor(i,j),vormax)
      vormin=min(vor(i,j),vormin)
      if(vor(i,j).gt.0.) vorp=vorp+vor(i,j)
      if(vor(i,j).lt.0.) vorm=vorm+vor(i,j)
  328 continue
      return
      end
c
c    subroutine outth
c    output in file time history
c
      subroutine outth(ntime,time,vor,psi,vormax,vormin)
      include 'param.f'
      dimension vor(m1,m2),psi(m1,m2)
      common/camo/enerv,vorip,vorim,ensti,cflm
      common/hioo/oor3,oor4,oor6,oor8,oor10   
      common/hioo0/oo0r3,oo0r4,oo0r6,oo0r8,oo0r10   
      common/camo0/vmi0,ens0,ene0,circm0,pal0
      common/pale/palen
      rvma=vormax/vmi0
      rvmi=vormin/vmi0
      rvip=vorip/abs(circm0)
      rvim=vorim/abs(circm0)
      roo1=rvip+rvim
c     roo1=vorip+vorim
      rpal=palen/pal0
      rene=enerv/ene0
      rens=ensti/ens0
      roor4=oor4/oo0r4
      roor6=oor6/oo0r6
      roor8=oor8/oo0r8
      roor10=oor10/oo0r10
  158 format(1x,e10.4,1x,8(1x,e11.5),1x,e10.4)
      write(6,158)time,vormax,vormin,rvma,rvmi,rvip,rvim
     1            ,ensti,enerv,cflm
      write(22,158)time,roo1,oor3,rene,rens,rpal,roor4,roor6 
     1            ,roor8,cflm
      write(20,158)time,vormax,vormin,rvma,rvmi,rvip,rvim
     1            ,ensti,enerv,cflm
      return
      end
c
c  ****************************** subrout cfield **********************
c
c  this subroutine calculate the integral  quantities
c  vorticity , enstrophy and energy
c  Higher moments up to 10th of the vorticity are
c  calculated to check conservation properties
c  In addition the palenstrophy is calculated
c   that is an important quantity in 2D turbulence
c
      subroutine cfield(vor,psi)
      include 'param.f'
      common/mesh/dx1,dx1q,dx2,dx2q
      dimension vor(m1,m2),psi(m1,m2)
      common/dim/n1,n1m,n2,n2m
      common/camo/enerv,vorip,vorim,ensti,cflm
      common/pale/palen
      common/hioo/oor3,oor4,oor6,oor8,oor10   
      common/indbo/imv(m1),ipv(m1)
      common/indx2/jmv(m2),jpv(m2)
      common/coor/yp1(m1),yp2(m2)
      common/d1/tfin,alx1f,alx1i,alx2i,alx2f
      pi=2.*asin(1.)
      area=(alx1f-alx1i)*(alx2f-alx2i)
      uare=1./area
      anor=1./(dx1*dx2)*uare
      palen=0.
      vorip=0.
      vorim=0.
      ensti=0.
      enerv=0.
      oor4=0.
      oor3=0.
      oor6=0.
      oor8=0.
      oor10=0.
      do 410 i=1,n1m
      ip=ipv(i)
      im=imv(i)
      do 410 j=1,n2m
      jp=jpv(j)
      jm=jmv(j)
      dox=(vor(ip,j)-vor(im,j))*dx1*0.5
      doy=(vor(i,jp)-vor(i,jm))*dx2*0.5
      palen=palen+dox**2+doy**2
      enevc=vor(i,j)*psi(i,j)
      vorc=vor(i,j)
      enstc=vor(i,j)**2
      oor3c=vor(i,j)**3
      oor4c=vor(i,j)**4
      oor6c=vor(i,j)**6
      oor8c=vor(i,j)**8
      oor10c=vor(i,j)**10
      enerv=enerv+enevc*0.5
      if(vorc.gt.0.) vorip=vorip+vorc
      if(vorc.lt.0.) vorim=vorim+vorc
      ensti=ensti+enstc
      oor3=oor3+oor3C
      oor4=oor4+oor4C
      oor6=oor6+oor6C
      oor8=oor8+oor8C
      oor10=oor10+oor10C
  410 continue
      palen=palen*anor
      enerv=enerv*anor
      ensti=ensti*anor
      oor3=oor3*anor
      oor4=oor4*anor
      oor6=oor6*anor
      oor8=oor8*anor
      oor10=oor10*anor
      return
      end
c
c  ****************************** subrout cfl  **********************
c
c  in this subroutine is calculated the maximum courant number
c
c
      subroutine cfl(psi,cflm)
      include 'param.f'
      common/dim/n1,n1m,n2,n2m
      common/mesh/dx1,dx1q,dx2,dx2q
      dimension psi(m1,m2)
      common/indx1/imv(m1),ipv(m1)
      common/indx2/jmv(m2),jpv(m2)
c
      cflm=0.
      do 10 ic=1,n1m
      ip=ipv(ic)
      im=imv(ic)
      do 10 jc=1,n2m
      jm=jmv(jc)
      jp=jpv(jc)
      h22a=((psi(ic,jp)-psi(ic,jm))-
     1      (psi(ip,jc)-psi(im,jc)))
     1     *dx2*dx1*0.5
      cflm=amax1(abs(h22a),cflm)
   10 continue
      return
      end
c
c
c
c  ****************************** subrout meshes **********************
c
      subroutine meshes
      include 'param.f'
      common/dim/n1,n1m,n2,n2m
      common/mesh/dx1,dx1q,dx2,dx2q
      common/d1/tfin,alx1f,alx1i,alx2i,alx2f
c
c  evaluates grid size and then  the inverse
c
      dx1=(alx1f-alx1i)/float(n1m)
      dx2=(alx2f-alx2i)/float(n2m)
      dx1=1./dx1
      dx2=1./dx2
      dx1q=dx1*dx1
      dx2q=dx2*dx2
      return
      end
c
c  ****************************** subrout indic **********************
c
c  in this subroutine the indices ip,im and jp ,jm  are calculated
c  these are necessary for the non linear and viscous terms near the
c  points  at the boundaries where periodicity is assumed.
c
      subroutine indic
      include 'param.f'
      common/dim/n1,n1m,n2,n2m
      common/indx1/imv(m1),ipv(m1)
      common/indx2/jmv(m2),jpv(m2)
      common/in4x1/inv(m1),itv(m1)
      common/in4x2/jnv(m2),jtv(m2)
c
c   periodic direction x1
c
      do 11 ic=1,n1m
      ipv(ic)=ic+1
      if(ic.eq.n1m) ipv(ic)=1
      imv(ic)=ic-1
      if(ic.eq.1) imv(ic)=n1m
      itv(ic)=ic+2
      if(ic.eq.n1m-1) itv(ic)=1
      if(ic.eq.n1m) itv(ic)=2
      inv(ic)=ic-2
      if(ic.eq.1) inv(ic)=n1m-1
      if(ic.eq.2) inv(ic)=n1m
   11 continue
c
c   periodic direction  x2
c
      do 12 jc=1,n2m
      jpv(jc)=jc+1
      if(jc.eq.n2m) jpv(jc)=1
      jmv(jc)=jc-1
      if(jc.eq.1) jmv(jc)=n2m
      jtv(jc)=jc+2
      if(jc.eq.n2m-1) jtv(jc)=1
      if(jc.eq.n2m) jtv(jc)=2
      jnv(jc)=jc-2
      if(jc.eq.1) jnv(jc)=n2m-1
      if(jc.eq.2) jnv(jc)=n2m
   12 continue
      return
      end
c
c
c  ****************************** subrout coordi  **********************
c    this subroutine assign the physical coordinates as function
c    of the computational ones x1,x2
c  The two grids are uniform
c
      subroutine coordi
      include 'param.f'
      common/dim/n1,n1m,n2,n2m
      common/mesh/dx1,dx1q,dx2,dx2q
      common/coor/yp1(m1),yp2(m2)
      common/d1/tfin,alx1f,alx1i,alx2i,alx2f
      common/prinl/n1i,n1f,n2i,n2f
      pi=2.*asin(1.)
      do 65 i=1,n1
      x1=(i-1)/float(n1m)
      yp1(i)=alx1i+x1*(alx1f-alx1i)
   65 continue
      do 66 j=1,n2
      x2=(j-1)/float(n2m)
      yp2(j)=alx2i+x2*(alx2f-alx2i)
   66 continue
      call pricor(yp1,yp2)
      return
      end
c
c  ****************************** subrout pricor **********************3
c  the coordinates are printed for flow visualizations by TURB3D
c
      subroutine pricor(yp1,yp2)
      include 'param.f'
      dimension yp1(m1),yp2(m2)
      common/dim/n1,n1m,n2,n2m
  200 format(3x,'in pricor n1=',i3,3x,'n2=',i3)
      write(6,200) n1,n2
      write(18) n1,n2
      write(18) ((yp1(i),i=1,n1),j=1,n2),
     1          ((yp2(j),i=1,n1),j=1,n2)
      close(18)
      return
      end
c
c     ***  subr wrispe  ****************************************
c   in this routine the energy and enstrophy spectra are 
c   calculated and written in files with names initiating
c   with en for energy , es for enstrophy and pes for
c   palenstrophy
c
      subroutine wrispe(ntime,time,vor,psi)
      include 'param.f'
      dimension vor(m1,m2),psi(m1,m2)
      common/dim/n1,n1m,n2,n2m
      common/camo/enerv,vorip,vorim,ensti,cflm
      common/pspect/pen(1025),pes(1025),pak(1025),ppes(1025)
      common/tspect/enett,enstt,penstt
      common/hioo/oor3,oor4,oor6,oor8,oor10   
      common/visct/re
      common/tstep/dt
      common/wavek/akm(m1m,m2mh),kmax
      character*20 filena
      character*4 ptim
      itime=time+0.5
      write(ptim,88) itime
   88 format(i4.4)
      call calspe(vor,time)
      write(16,797)time,enett,enstt,penstt,enerv,ensti,oor4
  787 format(1x,' t =',e8.2,2x,'sum en(k) =',e12.5,2x,'tot ens=',2e12.5)
  797 format(1x,7e12.5)
      filena='en'//ptim//'.out'
      open(11,file=filena)
      kma=kmax/sqrt(2.)
      nfisp=11
      do 17 k=1,kma
          write(nfisp,*) pak(k),pen(k)
   17 continue
      close(11)
      filena='es'//ptim//'.out'
      open(11,file=filena)
      do 27 k=1,kma
          write(nfisp,*) pak(k),pes(k)
   27 continue
      close(11)
      filena='pes'//ptim//'.out'
      open(11,file=filena)
      nfisp=11
      do 37 k=1,kma
          write(nfisp,*) pak(k),ppes(k)
   37 continue
      close(11)
  219 format(3x,e12.4,6(3x,e12.4))
  220 format(1x,i5,1x,3(e12.4,1x))
c     close(nfisp)
      return
      end
c
c     ***  subr writfi  ****************************************
c   in this routine the vorticity field in the physical space are
c    written for flow visualizations by TURB3D package
c
      subroutine writfi(ntime,time,vor,psi,ru)
      include 'param.f'
      dimension vor(m1,m2),psi(m1,m2)
      dimension ru(m1,m2)
      common/dim/n1,n1m,n2,n2m
      common/camo/enerv,vorip,vorim,ensti,cflm
      common/visct/re
      common/tstep/dt
      character*20 filena
      character*4 ptim
      do 310 j=1,n2m
      psi(n1,j)=psi(1,j)
      vor(n1,j)=vor(1,j)
      ru(n1,j)=ru(1,j)
  310 continue
      do 311 i=1,n1
      psi(i,n2)=psi(i,1)
      vor(i,n2)=vor(i,1)
      ru(i,n2)=ru(i,1)
  311 continue
      itime=time+0.5
      write(ptim,88) itime
   88 format(i4.4)
      filena='vo'//ptim//'.dat'
      open(12,file=filena,form='unformatted')
      nfil=12
      rewind(nfil)
      write(nfil) n1,n2
      write(nfil) re,dt,enett,enstt,time
      write(nfil) ((vor(i,j),i=1,n1),j=1,n2)
c
c   here there is the possibility to visualize the distribution
c   of the non-linera terms
c
c    1           ,((ru(i,j),i=1,n1),j=1,n2)
      close(nfil)
      return
      end
c
c
c  ****************************** subrout calspe **********************
c  In this roitine the two-dimensional spectra of enstrophy and
c  palenstrophy are calculated
c
      subroutine calspe(vor,time)
      include 'param.f'
      common/pspect/pen(1025),pes(1025),pak(1025),ppes(1025)
      common/tspect/enett,enstt,penstt
      common/spetrc/ene(1025),ens(1025),akp(1025),pens(1025)
      dimension paen(1025),paes(1025)
      common/waves/an(m2mh),ap(m1m),ak2(m2mh),ak1(m1m)
      dimension vor(m1,m2)
      character*10 teste
      common/dim/n1,n1m,n2,n2m
      common/wavek/akm(m1m,m2mh),kmax
      complex psk(m1m,m2mh)
      common/coor/yp1(m1),yp2(m2)
      common/n2mhi/n2mh
c
c   fft trasform of vorticity
c
      call phwam(vor,psk)
c     anor=1./(4*asin(1.))
      anor=1.
      penstt=0.
      enstt=0.
      do 15 k=1,kmax
      ens(k)=0.1e-012
      pens(k)=0.1e-012
      rkm=k
   15 continue
      akp(1)=1.
      do 16 i=1,n1m
      do 16 j=1,n2mh
      rk2=ak1(i)**2+ak2(j)**2
      sk=sqrt(rk2)
      k=sk+0.5
      if(k.le.kmax.and.k.gt.0) then
      pens(k)=pens(k)+compmag(psk(i,j),psk(i,j))*rk2
      ens(k)=ens(k)+compmag(psk(i,j),psk(i,j))
      pes(k)=alog10(ens(k))
      ppes(k)=alog10(pens(k))
                    endif
   16 continue
      do k=2,kmax
      enstt=enstt+(ens(k)+ens(k-1))*0.5
      penstt=penstt+(pens(k)+pens(k-1))*0.5
      enddo
      psk(1,1)=0.
      do 694 i=1,n1m
      do 694 j=1,n2mh
      if(akm(i,j).gt.0.) then
      psk(i,j)=psk(i,j)/akm(i,j)**2
                         endif
 694  continue
      call spectr(psk,time)
      return
      end
c
c     ********************* spectr *******************************
c   Here the energy spectrum is evaluated by the vor(k1.k2)
c   evaluated in the previous subroutine
c
      subroutine spectr(psk,t)
      include 'param.f'
      common/pspect/pen(1025),pes(1025),pak(1025),ppes(1025)
      common/tspect/enett,enstt,penstt
      common/spetrc/ene(1025),ens(1025),akp(1025),pens(1025)
      dimension paen(1025),paes(1025)
      common/waves/an(m2mh),ap(m1m),ak2(m2mh),ak1(m1m)
      common/dim/n1,n1m,n2,n2m
      common/wavek/akm(m1m,m2mh),kmax
      complex psk(m1m,m2mh)
      common/n2mhi/n2mh
      common/spet0/ak0,en0,gam0
      common/inisp/ityp
c
      enett=0.
      do 15 k=1,kmax
      ene(k)=0.1e-012
   15 continue
      akp(1)=1.
      do 16 i=1,n1m
      do 16 j=1,n2mh
      rk2=ak1(i)**2+ak2(j)**2
      sk=sqrt(rk2)
      k=sk+0.5
      if(k.le.kmax.and.k.gt.0) then
      enel=compmag(psk(i,j),psk(i,j))*rk2*0.5
      ene(k)=ene(k)+enel
      akp(k)=k
      pen(k)=alog10(ene(k))
      pak(k)=alog10(akp(k))
                    endif
   16 continue
      do k=2,kmax
      enett=enett+(ene(k)+ene(k-1))*0.5
      enddo
      return
      end
