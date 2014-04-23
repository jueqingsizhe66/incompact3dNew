c   ********************* subr rearuu
c  this subroutine reads  the velocity correlations   
c  evaluated in a previous run to restart the time average
c
      subroutine rearuu(xrko,nmedr)
      include 'param.f'
      parameter (m3m=m3-1,m1m=m1-1,m12=2*m1m,m32=2*m3m)
      parameter (m2k=13)
      dimension xrko(3,m1m,m2k)
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/tauwal/cfnp
      common/jrrr/jri,jrf,djr,irejr,iruuca
      common/vscpd/ipol,vscmi,timp,tvsf,stresh
      character*10 ruuwfi
      character*1 ptyp         
      ptyp='s'
      ruuwfi='ruukw'//ptyp//'.out'
      n1mh=n1m/2+1
      open(70,file=ruuwfi,form='unformatted')
      rewind(70)
      read(70) nmedr
      jj=0    
      do 20 j=jri,jrf,djr
      jj=jj+1
      do 20 i=1,n1mh
      read(70)(xrko(l,i,jj),l=1,3)            
   20 continue
      close(70)
      return
      end
c
c   ********************* subr ruucal
c  this subroutine perform the calculation of velocity correlations   
c
      subroutine ruucal(qcap,dph,xrko,nmedr,q,nav)
      include 'param.f'
      parameter (m2k=13)
      parameter (m3m=m3-1,m1m=m1-1,m12=2*m1m,m32=2*m3m)
      dimension qcap(m1,m2,m3),q(3,m1,m2,m3)
      dimension dph(m1,m2,m3)
      dimension ruu(3,m1,m2k),xrko(3,m1m,m2k),ruun(3)
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/tauwal/cfo 
      common/y2sta/y2s(m2)
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      common/tstep/dt,beta,ren
      common/jrrr/jri,jrf,djr,irejr,iruuca
      common/vmean/vm(3,m2),vrms(6,m2)
      common/metria/caj(m2),cac(m2)
      character*10 ruuwfi,ruupl
      character*2 pypl         
      character*1 ptyp         
      n1mh=n1m/2+1
      ptyp='s'
c
c   the correlations are evaluated at the cell center
c
      ruuwfi='ruukw'//ptyp//'.out'
      do 10 l=1,3
      if(l.eq.1)  then
      do 13 j=1,n2m  
      do 13 k=1,n3m
      do 13 i=1,n1m
      q1m=(q(1,ipv(i),j,k)+q(1,i,j,k))*0.5
      vm1m=q1m-vm(1,j)
      qcap(i,j,k)=vm1m      
   13 continue
                  endif
      if(l.eq.2)  then
      do 11 j=1,n2m  
      jp=j+1
      do 11 k=1,n3m
      do 11 i=1,n1m
      q2m=(q(2,i,jp,k)+q(2,i,j,k))*0.5
      vm2m=q2m-vm(2,j)
      qcap(i,j,k)=vm2m      
   11 continue
                  endif
      if(l.eq.3)  then
      do 15 j=1,n2m  
      do 15 k=1,n3m
      kp=kpv(k)
      do 15 i=1,n1m
      vm3m=(q(3,i,j,kp)+q(3,i,j,k))*0.5
      qcap(i,j,k)=vm3m-vm(3,j)
   15 continue
                  endif
c
c   the spanwise correlations are evaluated in this routine
c
      call ruphkp(qcap,dph,xrko,ruu,nmedr,l)
   10 continue
      cfnp=cfo/nav
c
c    here the correlations are written in a file
c
      open(79,file=ruuwfi,form='unformatted')
      rewind(79)
      write(79) nmedr
      write(36,*)cfnp,ren
      cplus=sqrt(cfnp)*ren
      jj=0    
      do 20 j=jri,jrf,djr
      jj=jj+1
      if(y2s(j).le.0) yplus=(1.+y2s(j))*cplus
      if(y2s(j).ge.0) yplus=(1.-y2s(j))*cplus
      iplus=yplus
      write(pypl,98)j
   98 format(i2.2)
      do 21 i=1,n1mh
      write(79)(xrko(l,i,jj),l=1,3)            
   21 continue
      ruupl='rup'//ptyp//pypl//'.dat'
      write(36,*)ruupl
      open(81,file=ruupl)                    
      do 34 l=1,3
      ruun(l)=ruu(l,1,jj)
   34 continue
      do 22 i=1,n1m
      x1plus=(i-0.5)/dx1*cplus 
      do 32 l=1,3
      ruu(l,i,jj)=ruu(l,i,jj)/ruun(l)
   32 continue
      write(81,171)x1plus,(ruu(l,i,jj),l=1,3)            
  171 format(1x,4e12.4)
   22 continue
   20 continue
      close(79)
      close(81)
      return
      end
c
c   ********************* subr ruphkp
c  this subroutine perform the calculation of the correlations by
c  using the FFT. The direct calculations in the physical space
c  is much more time consuming.
c  A real FFT in the streamwise x3 direction is followed
c  by a coplex FFT in x1
c
      subroutine ruphkp(qcap,dph,xrko,ruu,nmedr,l)
      include 'param.f'
      parameter (m2k=13)
      parameter (m3m=m3-1,m1m=m1-1,m12=2*m1m,m32=2*m3m)
      dimension qcap(m1,m2,m3)
      dimension dph(m1,m2,m3)
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/fftcmk/ifxk(13),trigxk(3*m1m/2+1)
      common/fftcm1/ifx1(13),trigx1(3*m1m/2+1)
      common/fftcm3/ifx3(13),trigx3(m32)
      common/waves/an(m3),ap(m1),ak3(m3),ak1(m1)
      complex xa(m3m,m1m),wor(m3m,m1m)
      real xr(m1m+2,m3m),work(m1m+1,m3m)
      real xrk(m1m+2,m2k),workk(m1m+1,m2k)
      dimension ruu(3,m1,m2k),xrko(3,m1m,m2k)
      common/rhsc/rhs(m1,m2,m3)
      common/jrrr/jri,jrf,djr,irejr,iruuca
      n1mh=n1m/2+1
      n1md=n1m+2
      do 1 j=jri,jrf,djr
      do 11 k=1,n3m
      xr(1,k)=qcap(n1m,j,k)
      xr(n1md,k)=qcap(1,j,k)
   11 continue
      do 12 i=1,n1m
      is=i+1
      do 12 k=1,n3m
      xr(is,k)=qcap(i,j,k)
   12 continue
c
c   2-d   fft applied to the each velocity component  
c   from physical to wave number space
c
      call fft99(xr,work,trigx1,ifx1,1,m1+1,n1m,n3m,-1)
      do 33 i=1,n1mh
      ip=2*i
      id=2*i-1
      do 33 k=1,n3m
      xa(k,i)=cmplx(xr(id,k),xr(ip,k))
  33  continue
c
c  2-d complex fft
c
      call cfft99(xa,wor,trigx3,ifx3,1,m3m,n3m,n1mh,-1)
      do 13 k=1,n3m
      do 13 i=1,n1mh
      qcap(i,j,k)=real(xa(k,i)/(n3m))
      rhs(i,j,k)=aimag(xa(k,i)/(n3m))
   13 continue
    1 continue
c
c   summation on wave numbers k3 to have correlations 
c   in spanwise direction x1
c
      jj=0    
      do 2 j=jri,jrf,djr
      jj=jj+1
      do 21 i=1,n1mh
      xrkijj=0.
      do 22 k=1,n3m
      xrkijj=qcap(i,j,k)**2+rhs(i,j,k)**2+xrkijj
   22 continue
      ip=2*i
      id=2*i-1
      xrk(id,jj)=xrkijj/n3m+xrko(l,i,jj)
      xrk(ip,jj)=0.           
      xrko(l,i,jj)=xrk(id,jj)  
   21 continue
    2 continue
      jjm=jj
      do 44 j=1,jjm
      do 44 i=1,n1m
      xrk(i,j)=xrk(i,j)/nmedr
  44  continue

c
c   correlation     in  physical space
c
c
      call fft99(xrk,workk,trigx1,ifx1,1,m1+1,n1m,jjm,+1)
      do 42 i=1,n1m
      is=i+1
      do 42 j=1,jjm
      ruu(l,i,j)=xrk(is,j)
   42 continue
      return
      end
