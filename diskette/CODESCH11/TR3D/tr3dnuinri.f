c ************************************************
c  
c    here the initial conditions for a vortex ring
c    impinging a solid wall are given
c  
c ***********************************************                       
      subroutine inring
      include 'param.f'
      dimension psi(m2,m3)                                              
      pi=2.*asin(1.)                                                    
      sig2=sig**2                                                       
      vmx=vmx/(pi*sig2)                                                 
      vorth = 0.
      write(6,*)'ring  ',vmx,sig,yc2mo,yc3mo,vper,epsil,lamb
c
c   a Gaussian azimuthal vorticity is distributed along
c   a radius perturbed in the radial direction
c   the azimuthal perturbation has a period given
c   by lamb and the simulation can be performed 
c         for   0 < theta < 2/pi  when irid.eq.0
c    or   for   0 < theta < 2/pi/lamb  when irid.eq1
c    to get a better resolution
c    remember that when irid = 1
c    dx1=2.*pi/(float(n1m)*float(lamb))
c
c    in the perturbation the azimuthal vorticity has
c    been projected on the axis
c
      do i=1,n1    
        if(irid.eq.0) then
          anper=(i-1)*2.*pi/float(n1m)*float(lamb)              
          app=-epsil*sin(2.*pi*(i-1)/float(n1m)*float(lamb))            
        else
          anper=(i-1)*2.*pi/float(n1m)
          app=-epsil*sin(2.*pi*(i-1)/float(n1m))            
        end if
        yteta=yc2mo+epsil*cos(anper)
        alpha=atan(app)                                                    
c                                                                       
        do j=1,n2                                                      
          do k=1,n3                                                     
            y1dm=rc(j)-yteta                                                 
            y1dp=rc(j)+yteta                                                 
            y2d=zz(k)-yc3mo                                                 
            raqm= y1dm**2 + y2d**2                                           
            raqp= y1dp**2 + y2d**2                                           
            ru1(i,j,k)=vmx*(
     %                   +exp(-raqm/sig2)                                     
     %                   -exp(-raqp/sig2)
     %                                    )*cos(alpha)
            if(abs(ru1(i,j,k)).gt.vorth) then                                
              vorth=abs(ru1(i,j,k))                                            
              imth=i                                                            
              jmth=j                                                            
              kmth=k                                                            
            endif                               
          end do
        end do
      end do
C
C     RADIAL VORTICITY IS COMPUTED FROM   DIV ( OME ) = 0
C     (OM_Z IS SET EQUAL TO 0)
C
      vorra = 0.
              imra=0
              jmra=0
              kmra=0
      do k=1,n3m
        do i=1,n1m
          im=imv(i)
          ru2(i,1,k)=0.
          do j=2,n2m
            jm=j-1
            ru2(i,j,k)=(  rm(jm)*ru2(i,jm,k)
     %        -dx1/dx2*g2rc(j)*(ru1(i,j,k)-ru1(im,j,k)))/rm(j)
            if(abs(ru2(i,j,k)).gt.vorra) then
              vorra=abs(ru2(i,j,k))
              imra=i
              jmra=j
              kmra=k
            endif
          end do
        end do
      end do
      do j=1,n2m
        do i=1,n1m
          ru2(i,j,n3)=0.
        end do
      end do
      do k=1,n3
        do i=1,n1m
          ru2(i,n2,k)=0.
        end do
      end do
      do k=1,n3
        do j=1,n2
          ru2(n1,j,k)=ru2(1,j,k)
        end do
      end do
c
c     THE AXIAL VORTICITY WAS SET = 0
C
C     THE PASSIVE SCALAR IS ASSIGNED AS THE VORTICITY MODULUS
C
      pscma = -10.
      do k=1,n3
        do j=1,n2
          do i=1,n1
            ru3(i,j,k)=0.
            psc(i,j,k)=sqrt(ru1(i,j,k)**2+ru2(i,j,k)**2)
            if(abs(psc(i,j,k)).gt.pscma) then                                
              pscma=abs(psc(i,j,k))                                            
              imps=i                                                            
              jmps=j                                                            
              kmps=k                                                            
            endif                               
          end do
        end do
      end do
      write(6,*)' OMt = ',vorth,' at i= ',imth,' j= ',jmth,' k= ',kmth
      write(6,*)' OMr = ',vorra,' at i= ',imra,' j= ',jmra,' k= ',kmra
      write(6,*)' Sca = ',pscma,' at i= ',imps,' j= ',jmps,' k= ',kmps
C
C
C     THE VELOCITY FIELD IS COMPUTED FROM  LAPL ( VEL ) = - CURL ( OME )
C
C     COMPUTE THE AXIAL COMPONENT OF CURL ( OME )
C
      do k=2,n3m
        do j=1,n2m
        jp=j+1
        do i=1,n1m
          ip=ipv(i)
          qcap(i,j,k)=-(
     %           (rc(jp)*ru1(i,jp,k)-rc(j)*ru1(i,j,k) )*dx2/g2rm(j)
     %               - ( ru2(ip,j,k)-ru2(i,j,k) )*dx1  
     %                                                )/rm(j)
          end do
        end do
      end do
      do j=1,n2m
        do i=1,n1m
          qcap(i,j,1)  = 0.
          qcap(i,j,n3) = 0.
        end do
      end do
C
C     AXIAL VELOCITY COMPONENT
C
C
c   Calculation of quantities necessary to the fft for
c   the velocity calculation
c
      call phinv
      call phcalv
C
      q3max=0.
      do k=1,n3
        do j=1,n2m
          do i=1,n1m
            q3(i,j,k)=dph(i,j,k)
            if(abs(q3(i,j,k)).gt.q3max) then
              iq3=i
              jq3=j
              kq3=k
              q3max=abs(q3(i,j,k))
            endif
          end do
        end do
      end do
C
C     COMPUTE THE AZIMUTHAL VELOCITY COMPONENT FROM OME_r definition
C
      q1max=0.
      do j=1,n2m
        do i=1,n1m
        im=imv(i)
        q1(i,j,1)=0.
        do k=2,n3m
          km=kmv(k)
          q1(i,j,k)=( 
     %    1./rm(j)*(q3(i,j,k)-q3(im,j,k))*dx1-ru2(i,j,k)
     %                     )/dx3*g3rc(k) + q1(i,j,km)
             if(abs(q1(i,j,k)).gt.q1max) then
               iq1=i
               jq1=j
               kq1=k
               q1max=abs(q1(i,j,k))
             endif
          end do
        end do
      end do
C
C     COMPUTE THE RADIAL VELOCITY COMPONENT FROM DIV ( VEL ) = 0
C
      do k=1,n3m
        do i=1,n1m
          q2(i,1,k)=0.
          q2(i,n2,k)=0.
        end do
      end do
      q2max=0.
      do k=1,n3m
        kp=kpv(k)
        do i=1,n1m
          q2(i,1,k)=0.
          ip=ipv(i)
          do j=1,n2m
            jp=j+1
            q2(i,jp,k)=q2(i,j,k)-(
     %      (q1(ip,j,k)-q1(i,j,k))*dx1+
     %      rm(j)*(q3(i,j,kp)-q3(i,j,k))*dx3/g3rm(k) )/dx2*g2rm(j)
            if(abs(q2(i,j,k)).gt.q2max) then
              iq2=i
              jq2=j
              kq2=k
              q2max=abs(q2(i,j,k))
            endif
          end do
        end do
      end do
      write(6,*)  '  '
      write(6,*)  ' V E L O C I T I E S  '
      write(6,*)  '  '
      write(6,321) q1max,iq1,jq1,kq1
  321 format(' max q1  ',e11.4,' at  ',3(i3,2x))
      write(6,320) q2max,iq2,jq2,kq2
  320 format(' max q2  ',e11.4,' at  ',3(i3,2x))
      write(6,322) q3max,iq3,jq3,kq3
  322 format(' max q3  ',e11.4,' at  ',3(i3,2x))
      call vmaxv
      call divgck(qmax,qtot)                                   
      write(6,706) epsil,lamb                                           
      write(32,706) epsil,lamb                                          
  706 format(3x,'from inqpr  epsil lamb =',e11.4,2x,i3)                 
      write(6,705)imxq,jmxq,kmxq,qmax,qtot,dq3m                                            
      write(32,705) qmax,qtot,d3m                                           
  705 format(3x,'from inqpr i,j,k qmax and qtot dq3m =',3(i3,1x),3e11.4)                   
      write(6,700) (vmax(l),l=1,3)                                      
      write(32,700) (vmax(l),l=1,3)                                     
  700 format(1x,'vmx1=',e11.4,2x,'vmx2=',e11.4,2x,'vmx3=',e11.4)        
      if(n1m.eq.1) then
      aa=1.
      open(59,file='inifiel.dat',form='unformatted')
      rewind 59
      write(59) n3m,n2m,1
      write(59) aa,aa,aa,aa
      write(59)
     1   ((q1(1,j,k),k=1,n3m),j=1,n2m),
     1   ((q2(1,j,k),k=1,n3m),j=1,n2m),
     1   ((q3(1,j,k),k=1,n3m),j=1,n2m),
     1   ((psc(1,j,k),k=1,n3m),j=1,n2m),
     1   ((ru1(1,j,k),k=1,n3m),j=1,n2m),
     1   ((ru2(1,j,k),k=1,n3m),j=1,n2m),
     1   ((ru3(1,j,k),k=1,n3m),j=1,n2m)
      close(59)
                   endif
      return                                                            
      end                                                               
c
c
c   ********************* subr phcalv
c  this subroutine perform the calculation of the
c   velocity field by solving the Laplacian V = - curl omega
c
      subroutine phcalv
      include 'param.f' 
      common /fftcm1/ ifx1(13),trigx1(m12)
      real xr(m1+1,m2m),work(m1,m2m)
      dimension acphvj(m2)

      my=m2
      pi=2.*asin(1.)
      n1mh=n1m/2+1
      if(n1m.gt.1) then
c
c   fft applied to the x1 direction to the RHS  of Laplacian V
c
      do k=1,n3
      do j=1,n2m
      xr(1,j)=qcap(n1m,j,k) 
      do i=1,n1m
      is=i+1
      xr(is,j)=qcap(i,j,k)
      enddo
      enddo
c   
c
c   real  fft 
c   from physical to wave number space
c
 
      call fft99(xr,work,trigx1,ifx1,1,m1+1,n1m,n2m,-1)
c
      do j=1,n2m
      do i=1,n1mh
      ip=2*i
      id=2*i-1
      qcap(i,j,k)=xr(id,j)
      dq(i,j,k)=xr(ip,j)
      enddo
      enddo
      enddo
                 endif
c
c   Now the fishpack is used to solve in the r-z plane at     
c      each wave number k1
c      Real part
c
            do i=1,n1mh
      do j=1,n2m
      acphvj(j)=acphjv(i,j)
      enddo
      do j=1,n2m
      do k=1,n3
      yfis(j,k)=qcap(i,j,k)
      enddo
      enddo
      call blktri(1,np,n3,amphkv,acphkv,apphkv
     1             ,mp,n2m,amphjv,acphvj,apphjv,my,yfis,ierror,w)
      do j=1,n2m
      do k=1,n3
      qcap(i,j,k)=yfis(j,k)
      enddo
      enddo
      if(n1m.gt.1) then
c
c      Immaginary part
c
      do j=1,n2m
      do k=1,n3
      yfis(j,k)=dq(i,j,k)
      enddo
      enddo
      call blktri(1,np,n3,amphkv,acphkv,apphkv
     1             ,mp,n2m,amphjv,acphvj,apphjv,my,yfis,ierror,w)
      do j=1,n2m
      do k=1,n3
      dq(i,j,k)=yfis(j,k)
      enddo
      enddo
                else
      do j=1,n2m
      do k=1,n3m
      dph(i,j,k)=qcap(i,j,k)
      enddo
      enddo
                endif
                enddo
c
c   phi in k1 wavenumber spacef
c
      if(n1m.gt.1) then
      n1mu=n1m-1
      do k=1,n3
      do j=1,n2m
      do i=1,n1mh
      ip=2*i
      id=2*i-1
      xr(id,j)=qcap(i,j,k)
      xr(ip,j)=dq(i,j,k)
      enddo
      enddo
c
c   from wave number to physical space
c
      call fft99(xr,work,trigx1,ifx1,1,m1+1,n1m,n2m,+1)
      do j=1,n2m
      dph(n1m,j,k)=xr(1,j)
      do i=1,n1mu
      is=i+1
      dph(i,j,k)=xr(is,j)
      enddo
      enddo
      enddo
                endif
      return
      end
c
c  ****************************** subrout phinv  **********************
c
c   in this subr the coefficients of the poisson eq. for Laplacian V
c   are calculated this routine is called only at the beginning
c
      subroutine phinv
      include 'param.f' 
      dimension acphjj(m2)
      n1mh=n1m/2+1
c
c   tridiagonal matrix coefficients to be used by FISHPACK
c   radial direction               
c
      do jc=1,n2m
      jm=jmv(jc)
      jp=jpv(jc)
      a22icc=rc(jc)*jmc(jc)*dx2q/g2rc(jc)
      a22icp=rc(jp)*jpc(jc)*dx2q/g2rc(jp)
      ac2=-(a22icc+a22icp)
      ugmmm=1./rm(jc)/g2rm(jc)
      amphjv(jc)=a22icc*ugmmm
      apphjv(jc)=a22icp*ugmmm
      acphjj(jc)=-(amphjv(jc)+apphjv(jc))
      enddo
      do ic=1,n1mh
      do jc=1,n2m
      acphjv(ic,jc)=acphjj(jc)-ak1(ic)/rm(jc)**2
      enddo
      enddo
      do ic=1,n1mh
c     write(6,*)' end phini coeff j=1,n2m'
      jc=1
c     write(6,101)jc,amphjv(jc),acphjj(ic,jc),apphjv(jc)
      jc=n2m
c     write(6,101)jc,amphjv(jc),acphjj(ic,jc),apphjv(jc)
  101 format(3x,i4,3x,3e12.5)
      enddo
c
c   tridiagonal matrix coefficients to be used by FISHPACK
c   axial direction               
c
      do kc=2,n3m
      km=kmv(kc)
      kp=kpv(kc)
      a33icc=dx3q/g3rm(km)
      a33icp=dx3q/g3rm(kc)
      ugmmm=1./g3rc(kc)
      amphkv(kc)=a33icc*ugmmm
      apphkv(kc)=a33icp*ugmmm
      acphkv(kc)=-(amphkv(kc)+apphkv(kc))
      enddo
      amphkv(1)= 0.
      apphkv(1)= 0.
      acphkv(1)= 1.
      amphkv(n3)= 0.
      apphkv(n3)= 0.
      acphkv(n3)= 1.
c 
      an3=n3
      ax=alog(an3)/alog(2.)
      k=ax+1
      L=2**(K+1)
      mw=(K-2)*L+K+5+MAX(2*n3,6*n2m)
      np=1
      mp=1
      do j=1,n2m
      do k=1,n3
      yfis(j,k)=1.
      enddo
      enddo
      my=m2
c
c    FISHPACK initiation
c
      call blktri(0,np,n3,amphkv,acphkv,apphkv
     1             ,mp,n2m,amphjv,acphjj,apphjv,my,yfis,ierror,w)
      return
      end
