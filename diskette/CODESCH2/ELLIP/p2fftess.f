c
c   ********************* subr fftes
c  this subroutine perform the calculation of trigz for temperton fft
c
      subroutine phines(vor,psi)
      include 'param.f' 
      common/dim/n1,n1m,n2,n2m
      common/mesh/dx1,dx1q,dx2,dx2q
      common/ifftin/iftin
      common/wav2f/an(m2),ap(m1),ak2(m2),ak1(m1)
      dimension vor(m1,m2),psi(m1,m2)
      pi=2.*asin(1.)
      n2mh=n2m/2
      n1mh=n1m/2
      n2mp=n2mh+1
      n1mp=n1mh+1
c
c     wave number definition
c
      do 16 k=1,n2mh
   16 an(k)=(k-1)*2.*pi
      do 17 k=n2mp,n2m
   17 an(k)=-(n2m-k+1)*2.*pi
      do 18 i=1,n1mh
   18 ap(i)=(i-1)*2.*pi
      do 19 i=n1mp,n1m
   19 ap(i)=-(n1m-i+1)*2.*pi
c
c   modified wave number
c
      do 26 k=1,n2m
      ak2(k)=2.*(1.-cos(an(k)/n2m))*dx2q
   26 continue
      do 28 i=1,n1m
      ak1(i)=2.*(1.-cos(ap(i)/n1m))*dx1q
   28 continue
      iftin=1
      call pcaes(vor,psi)
      iftin=0
      return
      end
c
c   ********************* subr pcaes
c  this subroutine perform the calculation of dph , periodic direction
c  along x3 and x1to use the real fourier transform
c
      subroutine pcaes(qcap,psi)
      include 'param.f' 
      real*8 aux1(42000),aux2(40000),aux3(m2m*m1m/2)
      real*8 aaux1(42000),aaux2(40000),aaux3(m2m*m1m/2)
      real*4 xa(m2m,m1m),xb(m2m+2,m1m)
      complex*8 y(m2m/2+1,m1m)
      common/ifftin/iftin
      dimension qcap(m1,m2),psi(m1,m2),rhs(m1,m2)
      common/dim/n1,n1m,n2,n2m
      common/wav2f/an(m2),ap(m1),ak2(m2),ak1(m1)
      n2mh=n2m/2+1
      n2md=n2m+2
      q1ma=0.
      q2ma=0.
      q3ma=0.
      q4ma=0.
      q5ma=0.
      do j=1,n2m
      do i=1,n1m
      xa(j,i)=qcap(i,j)
      q1ma=max(abs(qcap(i,j)),q1ma)
      enddo
      enddo
c
c   2-d   fft applied to the RHS of Laplacian by essl IBM
c
c     isign=+1 from physical to wave numbers
c     isign=-1 from wave to physical
c     scale=1. from phys to wav
c     scale=1./n1m from  wav to phys
c
c
c   from physical to wave number space
c
      if(iftin.eq.1) then
      naux1=42000
      naux2=40000
      naux3=n2m*n1m/2
      call srcft2 (1,xa,m2m,y,m2mh,n2m,n1m,1,1.0,aux1,naux1,
     #aux2,naux2,aux3,naux3)
                      else
C
      call srcft2 (0,xa,m2m,y,m2mh,n2m,n1m,1,1.0,aux1,naux1,
     #aux2,naux2,aux3,naux3)
                      endif
      do  i=2,n1m
      do  k=1,n2mh
      den=-1./(ak1(i)+ak2(k))
      rhs(i,k)=real(y(k,i)*den)
      psi(i,k)=aimag(y(k,i)*den)
      enddo
      enddo
      rhs(1,1)=0.
      psi(1,1)=0.
c
c   phi in wavenumber space
c
      do  i=1,n1m
      do  k=1,n2mh
      y(k,i)=cmplx(rhs(i,k),psi(i,k))
      q2ma=max(abs(rhs(i,k)),q2ma)
      enddo
      enddo
c
c   2-d fft applied to the psi by ESSL IBM
c
       sca=1.0/float(n1m*n2m)
      if(iftin.eq.1) then
C
c   from wave number space to physical space
c
       call scrft2 (1,y,m2mh,xb,m2p,n2m,n1m,-1,sca,
     # aaux1,naux1,aaux2,naux2,aaux3,naux3)
                     else
C
       call scrft2 (0,y,m2mh,xb,m2p,n2m,n1m,-1,sca,
     # aaux1,naux1,aaux2,naux2,aaux3,naux3)
C
                     endif
      do  k=1,n2m
      do  i=1,n1m
      psi(i,k)=xb(k,i)
      q4ma=max(abs(psi(i,k)),q4ma)
      enddo
      enddo
      write(6,*)'p2fftes ',q1ma,q2ma,q4ma
      return
      end
