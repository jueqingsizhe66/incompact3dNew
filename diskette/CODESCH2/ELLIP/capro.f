      program pssfftns
      include 'param.f'
      common/mesh/dx1,dx1q,dx2,dx2q
      dimension vor(m1,m2),psi(m1,m2),psa(m1,m2)
      dimension psi1(m1,m2),psi2(m1,m2)
      common/dim/n1,n1m,n2,n2m
      common/d1/alx1f,alx1i,alx2i,alx2f
      common/vort/y1c,y2c
      common/pert/amp,alpha,crad
      common/phco/iph
      common/boupa/n1b1,n1bn1,n2b1,n2bn2
      common/itype/ifish
      common/ncfft/nufft
      open(15,file='propro.d')
      read(15,301) dummy
      read(15,*) n1,n2
      read(15,301) dummy
      read(15,*) alx1i,alx1f,alx2i,alx2f
      read(15,301) dummy
      read(15,*) y1c,y2c
      read(15,301) dummy
      read(15,*) amp,alpha,crad
      read(15,301) dummy
      read(15,*) nfca,ifish,nufft
  301 format(a4)
      n1m=n1-1
      n2m=n2-1
      iph=0
      pi=2.*asin(1.)
c
c     mesh sizes calculations
c
      call meshes
      write(6,754)n1,n2
  754 format(10x,'n1=',i3,2x,'n2=',i3)
      write(6,755) dx1,dx2
  755 format(3x,'dx1=',e10.3,3x,'dx2=',e10.3,3x)
      call coordi
c
c     initial conditions
c     the vorticity and the analytical streamfunction are given
c
      do ncas=1,nfca
      read(15,301) dummy
      read(15,*) n1b1,n1bn1,n2b1,n2bn2
      call initia(vor,psa)
c
c    in the following routines with names starting with 
c    phin    the quantities necessary to the  different
c    FFT possibilities are calculated
c    
c
c   the following coefficients allow different boundary conditions
c   if =0 periodicity
c   if =1 Dirichlet
c   if =2 Neuwman
c
      if(ifish.eq.1) then
      call phini
                     else
      if(n2b1.eq.0.and.n1b1.eq.0.and.ifish.le.0) then
!      if(ifish.eq.-3) then
!      call phines(vor,psi)
!                      endif
      if(ifish.eq.-2) then
      call phin2p
                      endif
      if(ifish.eq.-1) then
      call phinip
                      endif
                                  endif
      if(n2b1.ge.1.and.n1b1.eq.1.and.ifish.eq.-1) then
      write(6,*)'  phinii for pnssc '
      call phinii
                                  endif
      if(n2b1.eq.1.and.n1b1.ge.1.and.ifish.eq.0) then
      write(6,*)'  phinij for pscns '
      call phinij
                                  endif
                     endif
c
c    the psi  (streamfunction) is calculated by
c    inverting the Laplacian
c
      call calc(vor,psi,psi1,psi2,psa)
c
c    the results are written
c
      call outpf(vor,psi,psi1,psi2,psa)
      enddo
      stop
      end
c
c  **************  subrout calc
c
      subroutine calc(vor,psi,psi1,psi2)
      include 'param.f'
      dimension vor(m1,m2),psi(m1,m2)
      dimension psi1(m1,m2),psi2(m1,m2)
      common/dim/n1,n1m,n2,n2m
      common/boupa/n1b1,n1bn1,n2b1,n2bn2
      common/itype/ifish
      common/ncfft/nufft
c
c  ********* calculation of the stream functio periodic in x1
c            and tridiag in vertical x2
      if(ifish.eq.1) then
      call phcal(vor,psi)
      write(6,*)'  fishpack used',nufft
                     else
      if(n2b1.eq.0.and.n1b1.eq.0.and.ifish.le.0) then
!      if(ifish.eq.-3) then
!      do n=1,nufft
!      call pcaes(vor,psi)
!      enddo
!      write(6,*)'  pcaes used',nufft
!                      endif
      if(ifish.eq.-2) then
      do n=1,nufft
      call p2fft(vor,psi)
      enddo
      write(6,*)'  p2fft used',nufft
                      endif
      if(ifish.eq.-1) then
      do n=1,nufft
      call psftpe(vor,psi)
      enddo
      write(6,*)'  psftpe used',nufft
                      endif
                                  endif
      if(n2b1.eq.1.and.n1b1.ge.1.and.ifish.eq.0) then
      do n=1,nufft
      call pscns(vor,psi)
      enddo
      write(6,*)'  pscns used',nufft
                                  endif
      if(n2b1.ge.1.and.n1b1.eq.1.and.ifish.lt.0) then
      do n=1,nufft
      call pnssc(vor,psi)
      enddo
      write(6,*)'  pnssc used',nufft
                                  endif
                     endif
      return
      end
c
c  ****************************** subrout meshes **********************
c
      subroutine meshes
      include 'param.f'
      common/dim/n1,n1m,n2,n2m
      common/mesh/dx1,dx1q,dx2,dx2q
      common/d1/alx1f,alx1i,alx2i,alx2f
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
c
c  ****************************** subrout coordi  **********************
c    this subroutine assign the physical coordinates as function
c    of the computational ones x1,x2
c    in this case are considered uniform grids
c
      subroutine coordi
      include 'param.f'
      common/dim/n1,n1m,n2,n2m
      common/mesh/dx1,dx1q,dx2,dx2q
      common/coor/yp1(m1),yp2(m2)
      common/d1/alx1f,alx1i,alx2i,alx2f
      pi=2.*asin(1.)
      do 65 i=1,n1
      x1=(i-1)/float(n1m)
      yp1(i)=alx1i+x1*(alx1f-alx1i)
   65 continue
      do 66 j=1,n2
      x2=(j-1)/float(n2m)
      yp2(j)=alx2i+x2*(alx2f-alx2i)
   66 continue
c
c   this file is necessary for contour plots following
c   the format for TURB3D package
c
!      open(18,file='coordi.dat',form='unformatted')
      open(18,file='coordi.dat')
      write(18,*) n1,n2

      write(18,19) ((yp1(i),i=1,n1),j=1,n2),
     1          ((yp2(j),i=1,n1),j=1,n2)
   19 format(f8.5,1x,f8.5)
      close(18)
      return
      end
c     ********************* initia *******************************
c  
      subroutine initia(vor,psa)
      include 'param.f'
      common/dim/n1,n1m,n2,n2m
      common/mesh/dx1,dx1q,dx2,dx2q
      common/coor/yp1(m1),yp2(m2)
      dimension vor(m1,m2),vorp(m1,m2),psa(m1,m2)
      common/vort/y1c,y2c
      common/pert/amp,alpha,crad
      common/boupa/n1b1,n1bn1,n2b1,n2bn2
      common/d1/alx1f,alx1i,alx2i,alx2f
      write(6,*)'in initia b.c.',n1b1,n1bn1,n2b1,n2bn2
      pi=2.*asin(1.)
      vorto=0.
      vortf=0.
      vorpm=0.
      ntvt=0
      nn1=2
      nn2=3
      ann=-(float(nn1**2)+float(nn2**2))
c
c   vorticity function of sin and cos
c
      do i=1,n1
      do j=1,n2
      y1r=(yp1(i)-yp1(1))/(alx1f-alx1i)*pi*2.
      y2r=(yp2(j)-yp2(1))/(alx2f-alx2i)*pi*2.
      if(n1b1.eq.0.and.n2b1.eq.1) then
      vor(i,j)=ann*sin(nn1*y1r)*cos(nn2*y2r)
      psa(i,j)=sin(nn1*y1r)*cos(nn2*y2r)
                                  endif
      if(n1b1.eq.1.and.n2b1.eq.0) then
      vor(i,j)=ann*cos(nn1*y1r)*sin(nn2*y2r)
      psa(i,j)=cos(nn1*y1r)*sin(nn2*y2r)
                                  endif
      if(n1b1.eq.0.and.n2b1.eq.0) then
      vor(i,j)=ann*sin(nn1*y1r)*sin(nn2*y2r)
      psa(i,j)=sin(nn1*y1r)*sin(nn2*y2r)
                                  endif
      if(n1b1.eq.1.and.n2b1.eq.1) then
      vor(i,j)=ann*sin(nn1*y1r)*sin(nn2*y2r)
      psa(i,j)=sin(nn1*y1r)*sin(nn2*y2r)
                                  endif
      if(n1b1.eq.1.and.n2b1.eq.2) then
      vor(i,j)=ann*sin(nn1*y1r)*cos(nn2*y2r)
      psa(i,j)=sin(nn1*y1r)*cos(nn2*y2r)
                                  endif
      if(n1b1.eq.2.and.n2b1.eq.2) then
      vor(i,j)=ann*cos(nn1*y1r)*cos(nn2*y2r)
      psa(i,j)=cos(nn1*y1r)*cos(nn2*y2r)
                                  endif
      if(n1b1.eq.2.and.n2b1.eq.1) then
      vor(i,j)=ann*cos(nn1*y1r)*sin(nn2*y2r)
      psa(i,j)=cos(nn1*y1r)*sin(nn2*y2r)
                                  endif
      vorto=vorto+vor(i,j)
c
c    eventual random disturbance 
c
c     xra=rand()*ampl*exp(-(rar**alpha*alpha-2.)**2)
      xra=0.
      vorp(i,j)=xra
      vorpm=vorpm+vorp(i,j)
      ntvt=ntvt+1
      enddo
      enddo
      vorpm=vorpm/float(ntvt)
      vorto=vorto/float(ntvt)
      do j=1,n2
      do i=1,n1
      if(abs(vor(i,j)).gt.0.) then
      vor(i,j)=vor(i,j)+vorp(i,j)-vorpm-vorto
                             endif
      vortf=vortf+vor(i,j)
      enddo
      enddo
      vortf=vortf/float(n1*n2)
      write(6,*)' ntvt, tot vor  , vorp,vortf ',ntvt,vorto,vorpm,vortf
      return
      end
      subroutine outpf(vor,psi,psi1,psi2,psa)
      include 'param.f'
      dimension vor(m1,m2),psi(m1,m2),psa(m1,m2)
      dimension psi1(m1,m2),psi2(m1,m2)
      common/dim/n1,n1m,n2,n2m
      common/boupa/n1b1,n1bn1,n2b1,n2bn2
      character*1 pn1b1,pn1bn1,pn2b1,pn2bn2
  101 format(i1.1)
      write(pn1b1,101) n1b1
      write(pn1bn1,101) n1bn1
      write(pn2b1,101) n2b1
      write(pn2bn2,101) n2bn2
      a1=1.
      ps1ma=-10.
      ps1mi=10.
      ps2ma=-10.
      ps2mi=10.
      psima=-10.
      psimi=10.
      vorma=-10.
      vormi=10.
c
c   in psi1 there is the error with respect the analytical 
c   solution
c
      do i=1,n1
      do j=1,n2
      psi1(i,j)=psi(i,j)-psa(i,j)
      ps1ma=amax1(ps1ma,psi1(i,j))
      ps1mi=amin1(ps1mi,psi1(i,j))
      psima=amax1(psima,psi(i,j))
      psimi=amin1(psimi,psi(i,j))
      vorma=amax1(vorma,vor(i,j))
      vormi=amin1(vormi,vor(i,j))
      enddo
      enddo
c
c   the results are written to perform contour plots
c   by the graphic routine TURB3d
c
      write(6,171)vormi,vorma
  171 format(3x,'max vor',2e12.4)
  172 format(3x,'max psi',6e12.4)
      write(6,172)psimi,psima,ps1ma,ps1mi
!      open(18,file='fu.'//pn1b1//pn1bn1//pn2b1//pn2bn2
!     1     ,form='unformatted')
      open(18,file='fu.'//pn1b1//pn1bn1//pn2b1//pn2bn2)
      write(18,*) n1,n2
      write(18,*) a1,a1,a1,a1
      write(18,20) ((vor(i,j),i=1,n1),j=1,n2),
     1          ((psi(i,j),i=1,n1),j=1,n2),
     1          ((psi1(i,j),i=1,n1),j=1,n2),
     1          ((psa(i,j),i=1,n1),j=1,n2)

   20 format (f8.5,1x,f8.5,1x,f8.5,1x,f8.5)
      close(18)
      return
      end

