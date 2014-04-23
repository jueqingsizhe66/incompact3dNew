c
c
c  ****************************** subrout initia **********************
c
c   initial conditions  the energy spectra is given and
c   from this the vorticity distribution in physical space
c   is obtained
c
      subroutine initia(vor,psi,ru)
      include 'param.f'
      dimension vor(m1,m2),psi(m1,m2)
      common/dim/n1,n1m,n2,n2m
      dimension ru(m1,m2)
      common/mesh/dx1,dx1q,dx2,dx2q
      common/d2/nstop,nprint,ntst,npin,nread,nwrit
      common/idue/timrea
      common/tspect/enett,enstt,penstt
      common/spetrc/ene(1025),ens(1025),akp(1025),pens(1025)
      common/wavek/akm(m1m,m2mh),kmax
      common/spet0/ak0,en0,gam0
      common/coeen/scps
      common/icci/icut
      common/dimrat/factpi
      integer rand
      character*20 filena
      character*3 ptim
      pi=acos(-1.)
      do 121 k=2,kmax
      rkm=k
      enq=specib(rkm)
      enm=specib(rkm-1)
      ene0=ene0+(enq+enm)*0.5
  121 continue
      scps=sqrt(en0/ene0)
          write(6,222) en0,ene0,scps
  222 format(3x,' en0, ene0 ',2e12.3,'scps =',e12.3)
      if(nread.eq.0) then
c
c   here the vorticity is calculated with the field in the
c   whole computational domain
c
      write(6,*)' initial field inside the domain icut=', icut
      if(icut.eq.0)       then
      call vorand(vor,vorpm)
      filena='voincu.dat'
      open(12,file=filena,form='unformatted')
      nfil=12
      rewind(nfil)
      write(nfil) n1m,n2m
      write(nfil) ((vor(i,j),i=1,n1m),j=1,n2m)
                          else
c
c    here the vorticity field generated in a certain
c    box is read and inserted in a larger  computational 
c    box 
c
      filena='voincu.dat'
      open(12,file=filena,form='unformatted')
      nfil=12
      rewind(nfil)
      read(nfil) n1mi,n2mi
      write(6,*)' read initial vort for n1mi,m2mi',n1mi,n2mi
      write(6,*)' from ',filena 
      read(nfil) ((vor(i,j),i=1,n1mi),j=1,n2mi)
      n1mq=n1m/2-n1mi/2 
      n2mq=n2m/2-n2mi/2 
      do i=1,n1mi
      do j=1,n2mi
      ru(i,j)=vor(i,j)
      enddo
      enddo
      do i=1,n1m
      do j=1,n2m
      vor(i,j)=0.
      enddo
      enddo
      do ii=1,n1mi
      do jj=1,n2mi
      i=ii+n1mq
      j=jj+n2mq
      vor(i,j)=ru(ii,jj)
      enddo
      enddo
      open(12,file='vortcut.dat',form='unformatted')
      nfil=12
      rewind(nfil)
      write(nfil) n1,n2
      write(nfil) re,re,re,re
      write(nfil) ((vor(i,j),i=1,n1),j=1,n2)
     1           ,((ru(i,j),i=1,n1),j=1,n2)
      close(nfil)

                           endif
c
c  the energy spectra is evaluated from the vorticity
c  field in the physical space
c
      t=0.
      write(6,*)'  solo turbolenza'
      call calspe(vor,t)
      kma=kmax/sqrt(2.)
c
c   the spectra is compared with the assigned one
c
      do l=1,kma
      al=l
      eni=specin(al)
      reni=ene(l)/eni
      write(23,123)al,ene(l),eni
  123 format(3x,5e12.4)
      enddo
      close(23)
                     endif
      if(nread.eq.1)  then
      itime=timrea
      write(ptim,88) itime
   88 format(i3.3)
      filena='vo'//ptim//'.dat'
      open(12,file=filena,form='unformatted')
      nfil=12
      rewind(nfil)
      read(nfil) n1,n2
      read(nfil) re,dt,enett,enstt,time
      read(nfil) ((vor(i,j),i=1,n1),j=1,n2)
      call calspe(vor,t)
      write(6,*)'  turbolenza da file  en turb=   ',enett
      enetu=enett
      write(16,798)svor,vopma
 798  format(1x, ' vort svor=',e11.4,2x,'max vor dipolo ',e12.4)
                   endif
c
c   evaluates streamfunction
c
      call phcal(vor,psi)
      if(icut.eq.0)       then
      open(12,file='vortcut.dat',form='unformatted')
      nfil=12
      rewind(nfil)
      write(nfil) n1,n2
      write(nfil) re,re,re,re
      write(nfil) ((vor(i,j),i=1,n1),j=1,n2)
     1           ,((ru(i,j),i=1,n1),j=1,n2)
     1           ,((psi(i,j),i=1,n1),j=1,n2)
      close(nfil)

                           endif
      vomax=-10000.
      vomin=10000.
      psmin=10000.
      psmax=-10000.
      do 4 j=1,n2m
      do 4 i=1,n1m
      if(vor(i,j).gt.vomax) vomax=vor(i,j)
      if(psi(i,j).gt.psmax) psmax=psi(i,j)
      if(vor(i,j).lt.vomin) vomin=vor(i,j)
      if(psi(i,j).lt.psmin) psmin=psi(i,j)
      ru(i,j)=0.
    4 continue
      write(6,791)vomax,vomin,psmax,psmin
      write(16,791)vomax,vomin,psmax,psmin
  791 format(3x,'in init vor max min=',2e11.4,3x,'ps max min=',2e11.4)
      return
      end
c
c   ********************* function specib
c  the desired spectrum is given
c
      function specib(rkm)
      common/spet0/ak0,en0,gam0
      common/inisp/ityp
      if(ityp.eq.2) then
      specib=rkm**4*exp(-(rkm/ak0)**2)
      if(specib.lt..1e-05) specib=.1e-05
                    endif
      if(ityp.eq.1) then
      specib=rkm*(1.+(rkm/ak0)**(gam0+1))**(-1)*ak0**3*exp(-1.)
                    endif
      if(ityp.eq.0) then
      specib=rkm*exp(-(rkm/ak0)**2)*ak0**3
      if(specib.lt..1e-05) specib=.1e-05
                    endif
c     specib=rkm**(-3)
      return
      end
c
c   ********************* function specin
c  the desired spectrum is given
c
      function specin(rkm)
      common/spet0/ak0,en0,gam0
      common/coeen/scps
      common/inisp/ityp
      if(ityp.eq.2) then
      specin=scps*rkm**4*exp(-(rkm/ak0)**2)
      if(specin.lt..1e-05) specin=.1e-05
                    endif
      if(ityp.eq.1) then
      specin=scps*rkm*(1.+(rkm/ak0)**(gam0+1))**(-1)*ak0**3*exp(-1.)
                    endif
      if(ityp.eq.0) then
      specin=scps*rkm*exp(-(rkm/ak0)**2)*ak0**3
      if(specin.lt..1e-05) specin=.1e-05
                    endif
c     specin=rkm**(-3)
      return
      end
c   ********************* function gasdev
c random phases
c
      subroutine gasdes(ised,gasdev)
      iset=0
      if(iset.eq.0) then
    1   continue
        v1=2.*rand()-1.
        v2=2.*rand()-1.
        r=v1**2+v2**2
        if(r.ge.1.) go to 1
        fac=sqrt(-2.*log(r)/r)
        gset=v1*fac
        gasdev=v2*fac
        iset=1
                     else
        gasdev=gset
        iset=0
                     endif
      return
      end
c
c     ********************* vorand *******************************
c   from energy spectra the streamfunction is calculated
c  in the wave number space and then by FFTs the distribution
c  in the physical space is obtained
c
      subroutine vorand(vor,vorpm)
      include 'param.f'
      common/dim/n1,n1m,n2,n2m
      common/waves/an(m2mh),ap(m1m),ak2(m2mh),ak1(m1m)
      common/wavek/akm(m1m,m2mh),kmax
      dimension vor(m1,m2),psi(m1,m2),x(97),enec(1025)
      dimension enei(1025),eneco(1025)
      complex psk(m1m,m2mh)
      complex ai,eth1,eth2
      common/spet0/ak0,en0,gam0
      common/spetrc/ene(1025),ens(1025),akp(1025),pens(1025)
      common/spetrp/pen(1025),pes(1025),pak(1025),ppes(1025)
      common/coor/yp1(m1),yp2(m2)
      common/n2mhi/n2mh
      common/inisp/ityp
      common/coeen/scps
      common/ened/enede
      character*10 teste
      character*20 filena
      character*3  ptim
      open(69,file='inispe.out',form='formatted')
      itit=0
      icont=0
      write(ptim,88) itit
  88  format(i3.3)
      ene0=0.
      pi=2.*asin(1.)
      nran=97
      ai=cmplx(0.,1.)
      do i=1,n1m
      iseda=i
      do j=1,n2mh
      k=sqrt(ak1(i)**2+ak2(j)**2)+0.5
      psk(i,j)=cmplx(0.,0.)
      if(akm(i,j).gt.0.) then
      rkm=k
      enq=specin(rkm)
      enei(k)=enq
c
c  following the Rogallo definition
c
      etaq=(enq/(pi*rkm**3))**(0.5)
c
      call gasdes(iseda,xra)
      phi=2.*pi*rand()
      aar=etaq*xra*cos(phi)
      bbr=etaq*xra*sin(phi)
      psk(i,j)=cmplx(aar,bbr)
                         endif
      enddo
      enddo
      enet=0.
      do k=1,kmax
      eneco(k)=0.1e-012
      enec(k)=0.1e-012
      enddo
      do 21 i=1,n1m
      do 21 j=1,n2mh
      rk2=ak1(i)**2+ak2(j)**2
      k=sqrt(rk2)+0.5
      enel=compmag(psk(i,j),psk(i,j))*rk2*0.5
      eneco(k)=eneco(k)+enel
   21 continue
      do i=1,n1m
      do j=1,n2mh
      rk2=ak1(i)**2+ak2(j)**2
      k=sqrt(rk2)+0.5
      rkm=k
      if(specin(rkm).gt..1e-05) then
      rene=sqrt(specin(rkm)/eneco(k))
      psk(i,j)=psk(i,j)*rene 
                                endif
      enel=compmag(psk(i,j),psk(i,j))*rk2*0.5
      enec(k)=enec(k)+enel
      enddo
      enddo
      enema=0.
      do k=2,kmax
      enet=enet+(enec(k)+enec(k-1))*0.5
      enema=max(enec(k),enema)
      enddo
          write(6,*) 'enet=',enet,'    enema=',enema
  169 format(3x,5e11.4)
      write(69,169)alog10(1.),alog10(enei(1))
     1             ,alog10(enec(1))
      do k=2,kmax
      ak=k
      write(69,169)alog10(ak),alog10(enei(k))
     1             ,alog10(enec(k))
      enddo
      t=0.
c
c  psk from wave number to psi in physical space
c
      call phwap(psi,psk)
      re=0.
      pra=0.
      time=0.
c
c    vorticity  from psi in physical space  
c
      call psvo(vor,psi)
      pspma=0.
      vopma=0.
      vorpm=0.
      do 31 i=1,n1m
      do 31 j=1,n2m
      vorpm=vorpm+vor(i,j)
      pspma=amax1(abs(psi(i,j)),pspma)
      vopma=amax1(abs(vor(i,j)),vopma)
   31 continue
      vorpm=vorpm/(n1m*n2m)
      write(6,781)vopma,pspma,vorpm
      write(16,781)vopma,pspma,vorpm
  781 format(2x,'in vorand max vor',e12.4,2x,'ps max'
     1          ,e12.4,' vorpm=' ,e12.4)
      return
      end
c  ****************************** subrout psvo **********************
      subroutine psvo(vor,psi)
c  ****************************** subrout psvo **********************
      include 'param.f'
      dimension vor(m1,m2),psi(m1,m2)
      common/dim/n1,n1m,n2,n2m
      common/indx1/imv(m1),ipv(m1)
      common/mesh/dx1,dx1q,dx2,dx2q
      common/indx2/jmv(m2),jpv(m2)
c
c   vorticity from stream function in phys. space at t=0.
c
      do 1 jc=1,n2m
      do 1 ic=1,n1m
      ip=ipv(ic)
      im=imv(ic)
      jm=jmv(jc)
      jp=jpv(jc)
      visder=dx1q*(psi(ip,jc)+psi(im,jc))+dx2q*(psi(ic,jm)+psi(ic,jp))
     1      -2.*(dx1q+dx2q)*psi(ic,jc)
      vor(ic,jc)=-visder
    1 continue
      return
      end
c   ********************* function compmag
c  this subroutine perform the calculation of product of
c  a1 by the complex conjugate of a2
      function compmag(a1,a2)
      complex a1,a2
      compmag=real(a1)*real(a2)+aimag(a1)*aimag(a2)
      return
      end
c  ****************************** subrout psini  **********************
c
c   in this subr the coefficients of the poisson eq. for psi
c   are calculated this subr. is called only at the beginning
c
      subroutine psini
      include 'param.f'
      common/dim/n1,n1m,n2,n2m
      common/mesh/dx1,dx1q,dx2,dx2q
      common/waves/an(m2mh),ap(m1m),ak2(m2mh),ak1(m1m)
      common/wavek/akm(m1m,m2mh),kmax
      common/fftrc/ifr(13),trigr(3*m2m/2+1)
      common/fftcc/ifc(13),trigc(2*m1m)
      common/d1/tfin,alx1i,alx1f,alx2i,alx2f
      common/n2mhi/n2mh
      call cftfax(n1m,ifc,trigc)
      call fftfax(n2m,ifr,trigr)
      pi=2.*asin(1.)
      n2mh=n2m/2+1
      n1mh=n1m/2+1
      n1mp=n1mh+1
c
c     wave number definition
c
      do 16 k=1,n2mh
   16 an(k)=(k-1)
      do 18 i=1,n1mh
   18 ap(i)=(i-1)
      do 19 i=n1mp,n1m
   19 ap(i)=-(n1m-i+1)
c
c   modified wave number
c
      do 26 k=1,n2mh
      ak2(k)=an(k)
   26 continue
      do 28 i=1,n1m
      ak1(i)=ap(i)
   28 continue
c
      kmax=0
      do 594 i=1,n1m
      do 594 j=1,n2mh
      akm(i,j)=sqrt(ak1(i)**2+ak2(j)**2)
      km=sqrt(ak1(i)**2+ak2(j)**2)+0.5
      kmax=max(kmax,km)
  594 continue
      write(6,*)'    kmax=',kmax
  788 format(1x,20f6.1)
      return
      end
