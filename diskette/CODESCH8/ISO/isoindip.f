c
c  ****************************** subrout inidip **********************
c
c   initial conditions for the simulation of three-dimensional
c   vortices in a  triperiodic box.
c    The vorticity distribution is assigned in the plane x1-x2 the
c    coorrespective velocity field obtained by the streamfunction
c    is perturbed in the x3 direction
c
      subroutine inidip(q,pr,qcap,dph)
c
      include 'param.f'
      parameter (m1m=m1-1,m2m=m2-1)
      dimension q(ndv,m1,m2,m3),pr(m1,m2,m3)
      dimension qcap(m1,m2,m3),dph(m1,m2,m3)
      dimension vor(m1,m2),psi(m1,m2)
      dimension vca1(m1,m2),vca2(m1,m2)
      dimension vme1(m1,m2),vme2(m1,m2),vme3(m1,m2)
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/corpri/yp1(m1),yp2(m2),yp3(m3)
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      common/tstep/dt,beta,ren
      common/rhs3p/dp3ns
      common/vpeini/vper,omtres
      common/speene/e(ndv,0:m3)
      common/numond/kx(m1),ky(m2),kz(m3),kkmax
      common/ispec/imic
      pi=4.*atan(1.)
      if(imic.eq.-1) then
c
c evaluation of vorticity streamfunction of Lamb dipole in
c the 2D plane x1,x2
c
      call dipol(vor,psi,vca1,vca2)
                     endif
      if(imic.eq.-3) then
c
c   a Gaussian vortex is assigned
c
      call circvo(vor,psi,vca1,vca2)
                     endif
      if(imic.eq.-2) then
      call vorgam(vor,psi,vca1,vca2)
                     endif
      if(n3m.eq.1) then
c
c   in 2D only the dipole withouth any perturbation is assigned
c
      do k=1,n3m
            do j=1,n2m
                  do i=1,n1m
       q(1,i,j,k)=vca1(i,j)
       q(2,i,j,k)=vca2(i,j)
       q(3,i,j,k)=0.
                  enddo
            enddo
      enddo
                     else
c
c   in 3D the dipole is perturbed by adding to the velocity distribution
c   a random perturbation on the 3 velocity componenets
c   Depending on the computer a different random number generator
c   should be used
c   here is for the digital
c
      isid=17
      call srand(isid)
c     call ranset(isid)
      v1max=0.
      v2max=0.
      v3max=0.
            do j=1,n2m
                  do i=1,n1m
       vme1(i,j)=0.
       vme2(i,j)=0.
       vme3(i,j)=0.
       do k=1,n3m
       q(1,i,j,k)=vca1(i,j)
       q(2,i,j,k)=vca2(i,j)
       if(abs(vor(i,j)).gt.omtres) then
       q1p=vper*rand()
       vme1(i,j)=vme1(i,j)+q1p
       q(1,i,j,k)=q(1,i,j,k)+q1p
       q2p=vper*rand()
       q(2,i,j,k)=q(2,i,j,k)+q2p
       vme2(i,j)=vme2(i,j)+q2p
       q3p=vper*rand()
       q(3,i,j,k)=0.+q3p
       vme3(i,j)=vme3(i,j)+q3p
c
c    here is the random number generator for g77 compiler
c    on personal computers with UNIX add fcv.o to make file
c
c      q(1,i,j,k)=q(1,i,j,k)+vper*rnd()
c      q(2,i,j,k)=q(2,i,j,k)+vper*rnd()
c      q(3,i,j,k)=0.+vper*rnd()
                              endif
       v1max=max(q(1,i,j,k),v1max)
       v2max=max(q(2,i,j,k),v2max)
       v3max=max(q(3,i,j,k),v3max)
       enddo
       vme3(i,j)=vme3(i,j)/n3m
       vme1(i,j)=vme1(i,j)/n3m
       vme2(i,j)=vme2(i,j)/n3m
                  enddo
            enddo
            do j=1,n2m
                  do i=1,n1m
       do k=1,n3m
       q(2,i,j,k)=q(2,i,j,k)-vme2(i,j)
       q(1,i,j,k)=q(1,i,j,k)-vme1(i,j)
       q(3,i,j,k)=q(3,i,j,k)-vme3(i,j)
       enddo
                  enddo
            enddo
      write(6,*)'vimax',v1max,v2max,v3max
                     endif
        al=1.
c
c   this part up to the end of the rotine can be also 
c   cancelled infact the previous real step in time
c   evalauates a free-divergent velocity field
c
      call divgck(q,qmax)
      print *,'assigned velocity divergenza massima=',qmax
c
c  this field is not free-divergent
c
      call divg(qcap,q,al)
c
c qcap now contains the divergence of the velocity field q
c
c  ********* calculation of  dph by fft in two
c            directions and tridiag in vertical
c
      call phcalc(qcap,dph)
c
c  ********* calculation of solenoidal vel field
c
      call updvp(dph,q,1.)
      call divgck(q,qmax)
      print *,'divergenza massima=',qmax
      v1max=0.
      v2max=0.
      v3max=0.
      do k=1,n3m
            do j=1,n2m
                  do i=1,n1m
       v1max=max(q(1,i,j,k),v1max)
       v2max=max(q(2,i,j,k),v2max)
       v3max=max(q(3,i,j,k),v3max)
                  enddo
            enddo
      enddo
      write(6,*)'dopo div=0 vimax',v1max,v2max,v3max
      return
      end
c
c
c  ****************************** subrout vorgam  **********************
c
c   initial vort. and psi obtanied by a distribution of
c   circulation
c   vor=smgam(x2)*exp(-(x1-yc1mo)/velmo)^2
c   where smgam(x2)=d cagam(x2)/dx2
c   the distribution of cagam is assigned
c   This vorticity distribution can be usefull to study
c    the stability of trailing vortices generated
c   by airfoils with a certain circulation profile
c
      subroutine vorgam(vor,psi,vca1,vca2)
      include 'param.f'
      parameter (m1m=m1-1,m2m=m2-1)
      dimension vor(m1,m2),psi(m1,m2),vock(m1,m2)
      dimension vca1(m1,m2),vca2(m1,m2),voin(m1,m2)
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/corpri/yp1(m1),yp2(m2),yp3(m3)
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/pardip/thet0,vsi,yc1mo,yc2mo,akmo,velmo 
      common/pargam/ell1,ell2,bw1,bw2
      common/anmod/ar1,ar2
      dimension cagam1(m2),cagam2(m2),cagam(m2),smgam(m2),vegam(m3)
      dimension y2s(m2)
      pi=2.*asin(1.)
c
      bw1p=bw1/pi
      bwha1=bw1p*0.5
      bw2p=bw2/pi
      bwha2=bw2p*0.5
      write(6,*)'circulation distrbution  parabolic'
     1     , yc1mo,yc2mo,velmo,vsi
      do j=1,n2m
      y2s(j)=(yp2(j)+yp2(j+1))*0.5
      enddo
c
c    cagam   distribution
c
      do j=1,n2m
      yhal=y2s(j)
      if(abs(yhal).lt.bwha1) then
      cagam1(j)=sqrt(1.-4.*(yhal/bw1p)**2)
                            else
      cagam1(j)=0.
                            endif
      if(abs(yhal).lt.bwha2) then
      cagam2(j)=sqrt(1.-4.*(yhal/bw2p)**2)
                            else
      cagam2(j)=0.
                            endif
      enddo
      smgam(n2)=0.
      smgam(1)=0.
      open(65,file='cagam.out',status='unknown')
      do j=1,n2m
      cagam(j)=ell1*cagam1(j)+ell2*cagam2(j)
      write(65,*)y2s(j),cagam(j)
      enddo
      close(65)
c
c   smgam calculation
c
      do j=2,n2m
      smgam(j)=(cagam(j)-cagam(j-1))*dx2
      enddo
      do i=1,n3m
      exx1=(yp1(i)-yc1mo)/velmo
      if(abs(exx1).lt.10) then
      vegam(i)=exp(-exx1**2)
                          else
      vegam(i)=0.
                          endif
      enddo
      open(65,file='vegam.out',status='unknown')
      do i=1,n3m
      write(65,*)yp3(i),vegam(i)
      enddo
      close(65)
      open(65,file='smgam.out',status='unknown')
      do i=1,n2m
      write(65,*)yp2(i),smgam(i)
      enddo
      close(65)
      vorma=-100.
      vormi=+100.
      totgao=0.
      totgan=0.
      aint=1./(dx1*dx2)
c
c    vorticity calculation in a x1-x2 plane   
c
      do i=1,n1m
      do j=1,n2m
      vor(i,j)=vegam(i)*smgam(j)
      if(vor(i,j).gt.vorma) then
      vorma=vor(i,j)
      imap=i
      jmap=j
                            endif
      if(vor(i,j).lt.vormi) then
      vormi=vor(i,j)
      imin=i
      jmin=j
                            endif
      if(vor(i,j).ge.0.) then
      totgao=totgao+vor(i,j)*aint
                        endif
      voin(i,j)=vor(i,j)
      enddo
      enddo
      vomax=max(abs(vorma),abs(vormi))
      write(6,791) vorma,vormi,totgao
  791 format(3x,'***  in  vormax =',e12.4,'  vormi = ',e12.4
     1      ,3x,'totgao =',e12.4)
      do i=1,n1m
      do j=1,n2m
c     vor(i,j)=voin(i,j)/totgao
      vor(i,j)=voin(i,j)
      totgan=totgan+vor(i,j)*aint
      enddo
      enddo
      write(6,795) totgan
  795 format('totgan =',e12.4)
c
c   streamfunction calculation as in code in Chapt. 4
c
      call pscalc(vor,psi)
c
c   2D velocity field from streamfunction
c
      v1max=0.
      v2max=0.
      vockma=0.
      do i=1,n1m
            do j=1,n2m
      vca1(i,j)=(psi(i,jpv(j))-psi(i,j))*dx2
      vca2(i,j)=-(psi(ipv(i),j)-psi(i,j))*dx1
      v1max=max(vca1(i,j),v1max)
      v2max=max(vca2(i,j),v2max)
c 
c   here the vorticity is evaluated by the velocity for check
c
      vorck=(psi(ipv(i),j)-2*psi(i,j)+psi(imv(i),j))*dx1q
     1     +(psi(i,jpv(j))-2*psi(i,j)+psi(i,jmv(j)))*dx2q
      vock(i,j)=vorck
      vockma=max(vorck,vockma)
            enddo
      enddo
      write(6,*)'in dipol vimax',v1max,v2max,' ma vor',vockma

      do i=1,n1m
      vor(i,n2)=vor(i,1)
      vock(i,n2)=vock(i,1)
      psi(i,n2)=psi(i,1)
      enddo
      do j=1,n2
      vor(n1,j)=vor(1,j)
      vock(n1,j)=vock(1,j)
      psi(n1,j)=psi(1,j)
      enddo
      tt=0.
c
c   initial 2D  vorticity and streamfunction written in 
c   file for visualizations by TURB3D
c  
      open(57,file='voini.dat',form='unformatted')
      write(57) n1,n2
      write(57) tt,tt,tt,tt
      write(57) ((psi(i,j),i=1,n1),j=1,n2)
     1         ,((vor(i,j),i=1,n1),j=1,n2)
     1         ,((voin(i,j),i=1,n1),j=1,n2)
     1         ,((vock(i,j),i=1,n1),j=1,n2)
      close(57)
      return
      end
c
c
c  ****************************** subrout dipol  **********************
c
c   initial vort. and psi for Lamb dipole  batchelor pg 535
c
      subroutine dipol(vor,psi,vca1,vca2)
      include 'param.f'
      parameter (m1m=m1-1,m2m=m2-1)
      dimension vor(m1,m2),psi(m1,m2),vock(m1,m2)
      dimension vca1(m1,m2),vca2(m1,m2),voin(m1,m2)
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/corpri/yp1(m1),yp2(m2),yp3(m3)
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/pardip/thet0,vsi,yc1mo,yc2mo,akmo,velmo 
      common/anmod/ar1,ar2
      pi=2.*asin(1.)
c
      write(6,*)'Lamb dipole', yc1mo,yc2mo,akmo,velmo,vsi
      vorma=-100.
      vormi=+100.
      do i=1,n1m
      do j=1,n2m
      y1d=yp1(i)-yc1mo
      y2d=yp2(j)-yc2mo
      ramod=sqrt(y1d**2+y2d**2)
      if(ramod.lt..1e-06)then
      vor(i,j)=0.
      else
      stho=y2d/ramod
      if(abs(stho).gt.1.) write(6,719) i,j,y2d,ramod
  719 format(3x,2i4,2x,'err in asin',2e12.4)
      astho=(abs(stho)-1.)
c
c     definizione theta
c
      if(astho.gt..0) stho=stho-.1e-06*stho/abs(stho)
      thet=asin(stho)
      if(y2d.lt.0..and.y1d.gt.0.) thet=pi*2.-abs(thet)
      if(y2d.lt.0..and.y1d.lt.0.) thet=pi+abs(thet)
      if(y2d.ge.0..and.y1d.lt.0.) thet=-abs(thet)+pi
c
      sth=sin(thet)*cos(thet0)-sin(thet0)*cos(thet)
      ramo=1.
c  
c   here there is the possibility to change the shape of the
c   patch where the vorticity is concentrated
c   e.g. an ellipse with axis aa1 and aa2  remember to
c   enter the values in the routine
c
c     ramo=sqrt(aa1**2*aa2**2/(aa2**2+(aa1**2-aa2**2)*sth**2))
      radmod=ramod/ramo
      if(radmod.ge.1.) then
      vor(i,j)=0.
      else
c
c   
c   Bessel functions **********************
c   
      arbesj=radmod*akmo
      vor(i,j)=2.*vsi*velmo*sth*bessj1(arbesj)/bessj0(akmo)*akmo
      if(vor(i,j).gt.vorma) then
      vorma=vor(i,j)
      imap=i
      jmap=j
                            endif
      if(vor(i,j).lt.vormi) then
      vormi=vor(i,j)
      imin=i
      jmin=j
                            endif
      endif
      endif
      voin(i,j)=vor(i,j)
      enddo
      enddo
      vomax=max(abs(vorma),abs(vormi))
      write(6,791) vorma,vormi
  791 format(3x,'***  in dipol vormax =',e12.4,'  vormi = ',e12.4)
      icdip=(imap+imin)/2
      jcdip=(jmap+jmin)/2
      write(6,*)'positions icdip,jdip,imap.jmap,imin,jmin',
     1          icdip,jcdip,imap,jmap,imin,jmin
c
c  here the two patches are separated each other by the distance
c  ar1 along x1
c  ar2  along x2
c
      idii=ar1*dx1 
      jdij=ar2*dx2 
      write(6,*)idii,jdij,' dipol center displaced ar1='
     1          ,ar1,'   ar2=',ar2
      do i=1,n1m
      do j=1,n2m
      vca1(i,j)=0.
      vca2(i,j)=0.
      enddo
      enddo
      imapn=imap+idii
      iminn=imin-idii
      jmapn=jmap+jdij
      jminn=jmin-jdij
      icdipn=(imapn+iminn)/2
      jcdipn=(jmapn+jminn)/2
      if(icdipn.ne.icdip) then
      idii=-idii
                          endif
      if(jcdipn.ne.jcdip) then
      jdij=-jdij
                          endif
      write(6,*)' points dipol  displaced idii,jdij ',idii,jdij
      do i=1,n1m
      do j=1,n2m
      if(vor(i,j).gt.0) then
      vca1(i+idii,j+jdij)=vor(i,j)
                        endif
      enddo
      enddo
      do i=1,n1m
      do j=1,n2m
      if(vor(i,j).lt.0) then
      vca2(i-idii,j-jdij)=vor(i,j)
                        endif
      enddo
      enddo
      vorma=-100.
      vormi=+100.
      do i=1,n1m
      do j=1,n2m
      vor(i,j)=vca1(i,j)+vca2(i,j)
      if(vor(i,j).gt.vorma) then
      vorma=vor(i,j)
      imap=i
      jmap=j
                            endif
      if(vor(i,j).lt.vormi) then
      vormi=vor(i,j)
      imin=i
      jmin=j
                            endif
      enddo
      enddo
      icdip=(imap+imin)/2
      jcdip=(jmap+jmin)/2
c
c   the new position of the minimum and maximum vorticity is printed
c
      write(6,*)'NEW positions icdip,jdip,imap.jmap,imin,jmin',
     1          icdip,jcdip,imap,jmap,imin,jmin
c
c    the streamfunction is calculated
c
      call pscalc(vor,psi)                                               
c
c   here the velocity field
c   and the vorticity from velocity for check
c
      v1max=0.
      v2max=0.
      vockma=0.
      do i=1,n1m
            do j=1,n2m
      vca1(i,j)=(psi(i,jpv(j))-psi(i,j))*dx2
      vca2(i,j)=-(psi(ipv(i),j)-psi(i,j))*dx1
      v1max=max(vca1(i,j),v1max)
      v2max=max(vca2(i,j),v2max)
      vorck=(psi(ipv(i),j)-2*psi(i,j)+psi(imv(i),j))*dx1q
     1     +(psi(i,jpv(j))-2*psi(i,j)+psi(i,jmv(j)))*dx2q
      vock(i,j)=vorck
      vockma=max(vorck,vockma)
            enddo
      enddo
      write(6,*)'in dipol vimax',v1max,v2max,' ma vor',vockma
      do i=1,n1m
      vor(i,n2)=vor(i,1)
      vock(i,n2)=vock(i,1)
      psi(i,n2)=psi(i,1)
      enddo
      do j=1,n2
      vor(n1,j)=vor(1,j)
      vock(n1,j)=vock(1,j)
      psi(n1,j)=psi(1,j)
      enddo
      tt=0.
c
c  vorticity and straemfunction written in a file for visualizations
c  by TURB3D graphycs package
c
      open(57,file='voini.dat',form='unformatted')
      write(57) n1,n2
      write(57) tt,tt,tt,tt
      write(57) ((psi(i,j),i=1,n1),j=1,n2)
     1         ,((vor(i,j),i=1,n1),j=1,n2)
     1         ,((voin(i,j),i=1,n1),j=1,n2)
     1         ,((vock(i,j),i=1,n1),j=1,n2)
      close(57)
      return
      end
c
c
c  ****************************** subrout circvo  **********************
c
c   initial vort. for a Lamb dipole with an assigned circuletion.
c   The Lamb dipole is evaluated at the begin
c
      subroutine circvo(vor,psi,vca1,vca2)
      include 'param.f'
      parameter (m1m=m1-1,m2m=m2-1)
      dimension vor(m1,m2),psi(m1,m2),vock(m1,m2)
      dimension vca1(m1,m2),vca2(m1,m2),voin(m1,m2)
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/corpri/yp1(m1),yp2(m2),yp3(m3)
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/pardip/thet0,vsi,yc1mo,yc2mo,akmo,velmo 
      common/anmod/ar1,ar2
      pi=2.*asin(1.)
c
      write(6,*)'Lamb dipole', yc1mo,yc2mo,akmo,velmo,vsi
      vorma=-100.
      vormi=+100.
      circu=0.
      do i=1,n1m
      do j=1,n2m
      y1d=yp1(i)-yc1mo
      y2d=yp2(j)-yc2mo
      ramod=sqrt(y1d**2+y2d**2)
      if(ramod.lt..1e-06)then
      vor(i,j)=0.
      else
      stho=y2d/ramod
      if(abs(stho).gt.1.) write(6,719) i,j,y2d,ramod
  719 format(3x,2i4,2x,'err in asin',2e12.4)
      astho=(abs(stho)-1.)
c
c     definizione theta
c
      if(astho.gt..0) stho=stho-.1e-06*stho/abs(stho)
      thet=asin(stho)
      if(y2d.lt.0..and.y1d.gt.0.) thet=pi*2.-abs(thet)
      if(y2d.lt.0..and.y1d.lt.0.) thet=pi+abs(thet)
      if(y2d.ge.0..and.y1d.lt.0.) thet=-abs(thet)+pi
c
      sth=sin(thet)*cos(thet0)-sin(thet0)*cos(thet)
      ramo=1.
      radmod=ramod/ramo
      if(radmod.ge.1.) then
      vor(i,j)=0.
      else
c
c   
c   Bessel functions **********************
c   
      arbesj=radmod*akmo
      vor(i,j)=2.*vsi*velmo*sth*bessj1(arbesj)/bessj0(akmo)*akmo
      if(vor(i,j).gt.0.) then
      circu=circu+vor(i,j)/(dx1*dx2)
                            endif
      if(vor(i,j).gt.vorma) then
      vorma=vor(i,j)
      imap=i
      jmap=j
                            endif
      if(vor(i,j).lt.vormi) then
      vormi=vor(i,j)
      imin=i
      jmin=j
                            endif
      endif
      endif
      voin(i,j)=vor(i,j)
      enddo
      enddo
      vomax=max(abs(vorma),abs(vormi))
      write(6,791) vorma,circu
  791 format(3x,'Lamb in dipol vormax =',e12.4,'  circu = ',e12.4)
      n2mh=n2m/2+1
      circg=0.
      yc2mo=yc2mo-1.
      akmoc=2.40482
      avg=1./(dx1*dx2)
      write(6,*)'circu data first half',n2mh,yc2mo,avg
      do i=1,n1m
      do j=1,n2mh
      y1d=yp1(i)-yc1mo
      y2d=yp2(j)-yc2mo
      ramod=sqrt(y1d**2+y2d**2)
      radmod=ramod/ramo
      if(radmod.ge.1.) then
      vor(i,j)=0.
      else
c
c   
c   Bessel functions **********************
c   
      arbesj=radmod*akmoc
      vor(i,j)=bessj0(arbesj)
                        endif
      circg=circg+abs(vor(i,j))*avg
      enddo
      enddo
      write(6,*)'   circg',circg
      circf1=0.
      do i=1,n1m
      do j=1,n2mh
      vor(i,j)=vor(i,j)*circu/circg
      circf1=circf1+abs(vor(i,j))*avg
      if(vor(i,j).gt.0.) then
                         endif
      if(vor(i,j).gt.vorma) then
      vorma=vor(i,j)
      imap=i
      jmap=j
                            endif
      if(vor(i,j).lt.vormi) then
      vormi=vor(i,j)
      imin=i
      jmin=j
                            endif
      enddo
      enddo
      vomax=max(abs(vorma),abs(vormi))
      write(6,*)'  first half vorma,vormi',vorma,vormi
      yc2mo1=yc2mo
      yc2mo=yc2mo+ar2
      write(6,*)'circu data second half',n2mh,yc2mo,yc2mo1
      circf2=0.
      do i=1,n1m
      do j=n2mh,n2m
      y1d=yp1(i)-yc1mo
      y2d=yp2(j)-yc2mo
      ramod=sqrt(y1d**2+y2d**2)
      radmod=ramod/ramo
      if(radmod.ge.1.) then
      vor(i,j)=0.
      else
c
c   
c   Bessel functions **********************
c   
      arbesj=radmod*akmoc
      vor(i,j)=-circu/circg*bessj0(arbesj)
                         endif
      circf2=circf2+abs(vor(i,j))*avg
      if(vor(i,j).gt.vorma) then
      vorma=vor(i,j)
      imap=i
      jmap=j
                            endif
      if(vor(i,j).lt.vormi) then
      vormi=vor(i,j)
      imin=i
      jmin=j
                            endif
      enddo
      enddo
      write(6,*)'  second half vorma,vormi',vorma,vormi
      write(6,793) vorma,circf1,circf2
  793 format(3x,'CIRCULAR  vormax =',e12.4,'  circf1,circf2 = ',2e12.4)
      call pscalc(vor,psi)                                               
      v1max=0.
      v2max=0.
      vockma=0.
      do i=1,n1m
            do j=1,n2m
      vca1(i,j)=(psi(i,jpv(j))-psi(i,j))*dx2
      vca2(i,j)=-(psi(ipv(i),j)-psi(i,j))*dx1
      v1max=max(vca1(i,j),v1max)
      v2max=max(vca2(i,j),v2max)
      vorck=(psi(ipv(i),j)-2*psi(i,j)+psi(imv(i),j))*dx1q
     1     +(psi(i,jpv(j))-2*psi(i,j)+psi(i,jmv(j)))*dx2q
      vock(i,j)=vorck
      vockma=max(vorck,vockma)
            enddo
      enddo
      write(6,*)'in dipol vimax',v1max,v2max,' ma vor',vockma

      do i=1,n1m
      vor(i,n2)=vor(i,1)
      vock(i,n2)=vock(i,1)
      psi(i,n2)=psi(i,1)
      enddo
      do j=1,n2
      vor(n1,j)=vor(1,j)
      vock(n1,j)=vock(1,j)
      psi(n1,j)=psi(1,j)
      enddo
      tt=0.
      open(57,file='voini.dat',form='unformatted')
      write(57) n1,n2
      write(57) tt,tt,tt,tt
      write(57) ((psi(i,j),i=1,n1),j=1,n2)
     1         ,((vor(i,j),i=1,n1),j=1,n2)
     1         ,((voin(i,j),i=1,n1),j=1,n2)
     1         ,((vock(i,j),i=1,n1),j=1,n2)
      close(57)
      return
      end
