c
c   ********************* subr fftqua
c  this subroutine perform the calculation of reduced wave numbers
c   initiate the fft routines in the ESSL IBM libraries
c
      subroutine fftqua
      include 'param.f' 
      common/ifftin/iftin
      pi=2.*asin(1.)
      n3mh=n3m/2
      n1mh=n1m/2
      n3mp=n3mh+1
      n1mp=n1mh+1
c
c     wave number definition
c
      do 16 k=1,n3mh
   16 an(k)=(k-1)*2.*pi
      do 17 k=n3mp,n3m
   17 an(k)=-(n3m-k+1)*2.*pi
      do 18 i=1,n1mh
   18 ap(i)=(i-1)*2.*pi
      do 19 i=n1mp,n1m
   19 ap(i)=-(n1m-i+1)*2.*pi
c
c   modified wave number
c
      do 26 k=1,n3m
      ak3(k)=2.*(1.-cos(an(k)/n3m))*dx3q
   26 continue
      do 28 i=1,n1m
      ak1(i)=2.*(1.-cos(ap(i)/n1m))*dx1q
   28 continue
c
c   initiation of fft routine
c
      iftin=1
      call phcalc
      iftin=0
      return
      end
c
c   ********************* subr phcalc
c  this subroutine perform the calculation of dph , periodic direction
c  along x3 and x1 
c
      subroutine phcalc
      include 'param.f' 
      parameter (m3mh=m3m/2+1,m3p=m3+1)
      real*8 aux1(22000),aux2(10000),aux3(m3m*m1m/2) 
      real*8 aaux1(22000),aaux2(10000),aaux3(m3m*m1m/2)
      real*4 xa(m3m,m1m),xb(m3m+2,m1m) 
      complex*8 y(m3m/2+1,m1m)
      common/ifftin/iftin
      n3mh=n3m/2+1
      n3fft=n3m
      n1fft=n1m
      do j=1,n2m
      do i=1,n1m
      do k=1,n3m
      xa(k,i)=qcap(i,j,k)
      enddo
      enddo
c
c   2-d   fft applied to the divg(qhat) by essl IBM
c   from physical to wave number space
c    
c     isign=+1 from physical to wave numbers
c     isign=-1 from wave to physical
c     scale=1. from phys to wav
c     scale=1./n1m from  wav to phys
c
c
      if(iftin.eq.1) then
      naux1=22000       
      naux2=10000      
      naux3=n3m*n1m/2 
      call srcft2 (1,xa,m3m,y,m3mh,n3m,n1m,1,1.0,aux1,naux1,   
     #aux2,naux2,aux3,naux3)  
                      else
C                           
      call srcft2 (0,xa,m3m,y,m3mh,n3m,n1m,1,1.0,aux1,naux1,  
     #aux2,naux2,aux3,naux3)  
                      endif                 
      do  k=1,n3mh
      do  i=1,n1m
      qcap(i,j,k)=real(y(k,i))
      dq(i,j,k)=aimag(y(k,i)) 
      enddo
      enddo
      enddo
c
c   solution of poisson equation real part
c
      call dsolv(qcap)
c
c   solution of poisson equation immag. part
c
      call dsolv(dq)
c
c   phi in wavenumber space
c
      do  j=1,n2m
      do  i=1,n1m
      do  k=1,n3mh
      y(k,i)=cmplx(qcap(i,j,k),dq(i,j,k))
      enddo
      enddo
c
c   2-d fft applied to the phi by ESSL IBM
c   from wave number space to physical space
c
       sca=1.0/float(n1m*n3m)             
      if(iftin.eq.1) then
C                                        
       call scrft2 (1,y,m3mh,xb,m3p,n3m,n1m,-1,sca,  
     # aaux1,naux1,aaux2,naux2,aaux3,naux3)
                     else
C                                         
       call scrft2 (0,y,m3mh,xb,m3p,n3m,n1m,-1,sca,
     # aaux1,naux1,aaux2,naux2,aaux3,naux3)  
C                     
                     endif     
      do  k=1,n3m
      do  i=1,n1m
      dph(i,j,k)=xb(k,i)
      enddo
      enddo
      enddo
      return
      end
c
c   ********************* subr dsolv
c
c   this subroutine performs the solution of the poisson equation
c   by solving a tridigonal matrix at each wave number k1 and k3
      subroutine dsolv(qkk)
      include 'param.f' 
      dimension qkk(m1,m2,m3)
      do 11 k=1,n3m/2+1
      do 15 j=1,n2m
      do 15 i=1,n1m
      fphj(i,j)=qkk(i,j,k)
      acphj(i,j)=acph(j)-ak1(i)/rm(j)**2-ak3(k)
      apphj(j)=apph(j)
      amphj(j)=amph(j)
   15 continue
c
c   zero wave number
c
      if(k.eq.1) then
      fphj(1,1)=0.
      acphj(1,1)=1.
      apphj(1)=0.
      amphj(1)=0.
      endif
c
c   tridiagonal inversion
c
      call tribj(n2m,n1m)
c
c  solution
c
      do 14 j=1,n2m
      do 14 i=1,n1m
      qkk(i,j,k)=qsbph(i,j)
   14 continue
   11 continue
      return
      end
c
c
c  ****************************** subrout tribj  **********************
c
      subroutine tribj(n,m)
      include 'param.f' 
      do 10 i=1,m
      bet(i)=acphj(i,1)
      qsbph(i,1)=fphj(i,1)/bet(i)
   10 continue
      do 11 j=2,n
      do 21 i=1,m
      gm(i,j)=apphj(j-1)/bet(i)
      bet(i)=acphj(i,j)-amphj(j)*gm(i,j)
      qsbph(i,j)=(fphj(i,j)-amphj(j)*qsbph(i,j-1))/bet(i)
   21 continue
   11 continue
      do 12 j=n-1,1,-1
      do 22 i=1,m
      qsbph(i,j)=qsbph(i,j)-gm(i,j+1)*qsbph(i,j+1)
   22 continue
   12 continue
      return
      end
c
c  ****************************** subrout phini  **********************
c
c   in this subr the coefficients of the poisson eq. for dph
c   are calculated this subr. is called only at the beginning
c
      subroutine phini
      include 'param.f' 
      call fftqua
c
c   tridiagonal matrix coefficients due to the radial finite-difference
c   discretization
c
      do 1 jc=1,n2m
c
c
      jm=jmv(jc)
      jp=jpv(jc)
      a22icc=rc(jc)*jmc(jc)*dx2q/g2rc(jc)
      a22icp=rc(jp)*jpc(jc)*dx2q/g2rc(jp)
      ac2=-(a22icc+a22icp)
      ugmmm=1./rm(jc)/g2rm(jc)
      amph(jc)=a22icc*ugmmm
      apph(jc)=a22icp*ugmmm
      acph(jc)=ac2*ugmmm
      write(66,166)jc,amph(jc),acph(jc),apph(jc)
  166 format(3x,i4,2x,3e12.5)
    1 continue
c
      return
      end
