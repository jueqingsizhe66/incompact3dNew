c
c   ********************* subr fftqua
c  this subroutine perform the calculation of trigz for temperton fft
c
      subroutine fftqua
      include 'param.f'
      parameter (m3m=m3-1,m1m=m1-1,m12=2*m1m,m32=2*m3m)
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/fftcm1/ifx1(13),trigx1(3*m1m/2+1)
      common/fftcm3/ifx3(13),trigx3(m32)
      common/waves/an(m3),ap(m1),ak3(m3),ak1(m1)
      common/numond/kx(m1),ky(m2),kz(m3),kkmax
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      pi=2.*asin(1.)
      n1mh=n1m/2
      n2mh=n2m/2
      n3mh=n3m/2
      n1mp=n1mh+1
      n2mp=n2mh+1
      n3mp=n3mh+1
c
c     wave number definition
c
      if(n3m.gt.1) then
      do 16 k=1,n3mh
       kz(k)=(k-1)
   16 an(k)=(k-1)*2.*pi
      do 17 k=n3mp,n3m
       kz(k)=-(n3m-k+1)
   17 an(k)=-(n3m-k+1)*2.*pi
      nx3fft=n3m
      call cftfax(nx3fft,ifx3,trigx3)
                   else 
      an(1)=0.
                   endif
      if(n1m.gt.1) then
      do 18 i=1,n1mh
      kx(i)=(i-1)
   18 ap(i)=(i-1)*2.*pi
      do 19 i=n1mp,n1m
      kx(i)=-(n1m-i+1)
   19 ap(i)=-(n1m-i+1)*2.*pi
                   else 
      ap(1)=0.
                   endif
      nx1fft=n1m
      if(n1m.gt.1) call fftfax(nx1fft,ifx1,trigx1)
c     numeri d'onda in y
      do j=1,n2mh
        ky(j)=(j-1)
      end do
      do j=n2mp,n2m
       ky(j)=-(n2m-j+1)
      end do
      akk=sqrt(float(n1mh*n1mh+n2mh*n2mh+n3mh*n3mh))
      kkmax=nint(akk)
      return
      end
c
c   ********************* subr phcalc
c  this subroutine perform the calculation of dph , periodic direction
c  along x3 and x1to use the real fourier transform
c
      subroutine phcalc(qcap,dph)
      include 'param.f'
      parameter (m3m=m3-1,m1m=m1-1,m12=2*m1m,m32=2*m3m)
      dimension qcap(m1,m2,m3)
      dimension dph(m1,m2,m3)
      real xr(m1m+2,m3m),work(m1m+1,m3m)
      complex xa(m3m,m1m),wor(m3m,m1m)
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/fftcm1/ifx1(13),trigx1(3*m1m/2+1)
      common/fftcm3/ifx3(13),trigx3(m32)
      common/waves/an(m3),ap(m1),ak3(m3),ak1(m1)
      common/rhsc/rhs(m1,m2,m3)
      n1mh=n1m/2+1
      n1md=n1m+2
      do 1 j=1,n2m
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
c   2-d   fft applied to the divg(qhat) by cfft99
c   from physical to wave number space
c
      if(n1m.gt.1) then
      call fft99(xr,work,trigx1,ifx1,1,m1+1,n1m,n3m,-1)
                   endif
      if(n3m.gt.1) then
      do 33 i=1,n1mh
      ip=2*i
      id=2*i-1
      do 33 k=1,n3m
      xa(k,i)=cmplx(xr(id,k),xr(ip,k))
  33  continue
c
c
      call cfft99(xa,wor,trigx3,ifx3,1,m3m,n3m,n1mh,-1)
      do 13 k=1,n3m
      do 13 i=1,n1mh
      qcap(i,j,k)=real(xa(k,i)/(n3m))
      rhs(i,j,k)=aimag(xa(k,i)/(n3m))
   13 continue
                   else 
      do 43 i=1,n1mh
      ip=2*i
      id=2*i-1
      do 43 k=1,n3m
      qcap(i,j,k)=xr(id,k)
      rhs(i,j,k)=xr(ip,k)
   43 continue
                   endif
    1 continue
      qcap(1,1,1)=0.
      rhs(1,1,1)=0.
c
c   solution of poisson equation real part
c
      call dsolv(qcap)
c
c   solution of poisson equation immag. part
c
      call dsolv(rhs)
c
c   phi in wavenumber space
c
      do 2 j=1,n2m
      if(n3m.gt.1) then
      do 21 k=1,n3m
      do 21 i=1,n1mh
      xa(k,i)=cmplx(qcap(i,j,k),rhs(i,j,k))
   21 continue
c
c   2-d fft applied to the phi by cfft99
c   from wave number space to physical space
c
      call cfft99(xa,wor,trigx3,ifx3,1,m3m,n3m,n1mh,+1)
      do 34 i=1,n1mh
      ip=2*i
      id=2*i-1
      do 34 k=1,n3m
      xr(id,k)=real(xa(k,i))
      xr(ip,k)=aimag(xa(k,i))
  34  continue
                   else
      do 44 i=1,n1mh
      ip=2*i
      id=2*i-1
      do 44 k=1,n3m
      xr(id,k)=qcap(i,j,k)
      xr(ip,k)=rhs(i,j,k) 
  44  continue
                   endif
c
      if(n1m.gt.1) then
      call fft99(xr,work,trigx1,ifx1,1,m1+1,n1m,n3m,+1)
                   endif
      do 22 i=1,n1m
      is=i+1
      do 22 k=1,n3m
      dph(i,j,k)=xr(is,k)
   22 continue
    2 continue
      return
      end
c
c   ********************* subr dsolv
c
c   this subroutine performs the solution of the poisson equation
c   by solving a tridigonal matrix at each wave number n and p
      subroutine dsolv(qk)
      include 'param.f'
      parameter (m1m=m1-1,m2m=m2-1,m3m=m3-1)
      dimension qk(m1,m2,m3)
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/waves/an(m3),ap(m1),ak3(m3),ak1(m1)
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      dimension amj(m1,m2),acj(m1,m2),apj(m1,m2),fj(m1,m2)
      n1mh=n1m/2+1
      do k=1,n3m
           do i=1,n1mh 
                do j=1,n2m
      acc = -2.*dx2q+(-ak1(i)-ak3(k))
      fj(i,j)=qk(i,j,k)/acc
      acj(i,j)= 1.
      apj(i,j)=dx2q/acc
      amj(i,j)=dx2q/acc
                enddo
      if(i.eq.1.and.k.eq.1) then
      acj(1,1)=1.
      apj(1,1)=0.
      amj(1,1)=0.
      fj(1,1)=0.
                            endif
           enddo
c
c   tridiagonal inversion
c
      call tripvv(amj,acj,apj,fj,1,n2m,1,n1mh)
c
c  solution
c
           do i=1,n1mh
                do j=1,n2m
      qk(i,j,k)=fj(i,j)
                enddo
           enddo
      enddo
      return
      end
c
c
c  ****************************** subrout phini  **********************
c
c   in this subr the coefficients of the poisson eq. for dph
c   are calculated this subr. is called only at the beginning
c
      subroutine phini
      include 'param.f'
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      common/tstep/dt,beta,ren
      common/indbo/imv(m1),ipv(m1),jmmv(m2),jppv(m2),kmv(m3),kpv(m3)
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/waves/an(m3),ap(m1),ak3(m3),ak1(m1)
      call fftqua
c
c   modified wave number x3 direct.
c
      if(n3m.gt.1) then
      do 26 k=1,n3m
      ak3(k)=2.*(1.-cos(an(k)/n3m))*dx3q
   26 continue
                   else 
      ak3(1)=0.
                   endif
c
c   modified wave number x1 direct.
c
      if(n1m.gt.1) then
      do 28 i=1,n1m
      ak1(i)=2.*(1.-cos(ap(i)/n1m))*dx1q
   28 continue
                   else 
      ak1(i)=0.
                   endif
      return
      end
c************************************************************************
c                                                                       *
c  ****************************** subrout tripvv  ***********************
c                                                                       *
c************************************************************************
      subroutine tripvv(ami,aci,api,fi,ji,jf,ni,nf)
      include 'param.f'
      dimension ami(m1,m2),aci(m1,m2),api(m1,m2),fi(m1,m2)
      dimension fei(m1,m2),qq(m1,m2),ss(m1,m2)
      dimension fnn(m1),p(m1)
c                                                                       
c     vectorized for right hand side and coefficients                   
c                                                                       
      ja = ji + 1                                                       
      jj = ji + jf                                                      
      do 20 k=ni,nf                                                     
      qq(k,ji) = -api(k,ji)/aci(k,ji)                                        
      ss(k,ji) = -ami(k,ji)/aci(k,ji)                                        
      fnn(k) = fi(k,jf)                                                   
      fi(k,ji) = fi(k,ji)/aci(k,ji)                                         
   20 continue                                                          
c                                                                       
c     forward elimination sweep                                         
c                                                                       
      do 10 j=ja,jf                                                     
      do 21 k=ni,nf                                                     
      p(k) =1./( aci(k,j) + ami(k,j)*qq(k,j-1))                              
      qq(k,j) = - api(k,j)*p(k)                                            
      ss(k,j) = - ami(k,j)*ss(k,j-1)*p(k)                                   
      fi(k,j) = ( fi(k,j) - ami(k,j)*fi(k,j-1))*p(k)                         
   21 continue                                                          
   10 continue                                                          
c                                                                       
c     backward pass                                                     
c                                                                       
      do 22 k=ni,nf                                                     
      ss(k,jf) = 1.                                                      
      fei(k,jf) = 0.                                                     
   22 continue                                                          
      do 11 i=ja,jf                                                     
      j = jj - i                                                        
      do 23 k=ni,nf                                                     
      ss(k,j) = ss(k,j) + qq(k,j)*ss(k,j+1)                                 
      fei(k,j) = fi(k,j) + qq(k,j)*fei(k,j+1)                               
   23 continue                                                          
   11 continue                                                          
      do 24 k=ni,nf                                                     
      fi(k,jf)=(fnn(k) - api(k,jf)*fei(k,ji) - ami(k,jf)*fei(k,jf-1))           
     &       /(api(k,jf)*ss(k,ji) + ami(k,jf)*ss(k,jf-1)  +aci(k,jf))           
   24 continue                                                          
c                                                                       
c     backward elimination pass                                         
c                                                                       
      do 12 i=ja,jf                                                     
      j = jj -i                                                         
      do 25 k=ni,nf                                                     
      fi(k,j) = fi(k,jf)*ss(k,j) + fei(k,j)                                 
   25 continue                                                          
   12 continue                                                          
      return                                                            
      end                                                               
