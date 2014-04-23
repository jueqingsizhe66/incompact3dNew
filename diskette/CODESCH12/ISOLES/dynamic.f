c***********************************************
c     dynamic subgrid-scale stress model 
c     this model can be used only without any passive scalar
c     this is the model dveloped at the CTR with the contraction
c     suggested by Lilley
c***********************************************
      subroutine dynamic(q,pr,ntime,time)
      include 'param.f'
      dimension q(ndv,m1,m2,m3),pr(m1,m2,m3)
      common/strain/st(m1,m2,m3,6)
      common/mijdy/mij(m1,m2,m3,6)
      real*4 mij,lijmij,mijmij
      common/qdyn/g(m1,m2,m3)
      common/smapra/gold(m1,m2,m3)
      dimension accdyn(m2),acclij(m2),accmij(m2)
      common/qles1/deltax1,deltay1,deltaz1,ell1c
      common/qles2/deltax2,deltay2,deltaz2,ell2c
      common/rhsc/rhs(m1,m2,m3)
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      common/d123/alx1,alx2,alx3
      common/sma/vis(m1,m2,m3)
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/ledat/ics0,cvisc,ifiltr,ibox
      common/outd1/cdyn(m2),alij(m2),amij(m2)
      common/outd2/cs(m2),lijmij(m2),mijmij(m2)
      common/d2/nstop,nprint,ntst,npin,npstf
      common/cflco/icfl,cflc,tpin,tprin,tfin
      common/tstep/dt,beta,ren
      common/resca/iresca,tresca,rlamas
      common/dycopr/cdynav,visav,avg,avmij,nif,cist,flback,flint
      dimension visout(m2)
      common/strat2/igrad,bvais
      common/strat1/istrat,rho0,gg,schm
      character*20 titfil
      character*3 icount
c      
      if(ntime.eq.1.and.ibox.eq.1) print *,'Costante c(t) '
      if(ntime.eq.1.and.ibox.eq.0) then
       print *,'Costante c(x,y,z,t)'
       print *,'*******************'
      end if
      ratio=2.
c
c      computation of Sij=0.5*(dui/dxj + duj/dxi)
c
      call strper(q)
c      
      do k=1,n3m
          do j=1,n2m
              do i=1,n1m
      g(i,j,k)=2.*(st(i,j,k,1)**2+st(i,j,k,2)**2+st(i,j,k,3)**2)+ 
     1         4.*(st(i,j,k,4)**2+st(i,j,k,5)**2+st(i,j,k,6)**2)
      g(i,j,k)=sqrt(g(i,j,k)) 
              enddo
          enddo
      enddo
c
c     now g contains |S|=sqrt(2*Sij*Sij) 
c     gold stores |S| for the computation of viscosity at the end
c     of the subroutine
c
      do k=1,n3m
          do j=1,n2m
              do i=1,n1m
      gold(i,j,k)= g(i,j,k)
              enddo
          enddo
      enddo
c
c     pr=|S|Sij
c
            do n=1,6
              do k=1,n3m
          do j=1,n2m
      do i=1,n1m
      pr(i,j,k)=st(i,j,k,n)*g(i,j,k)
              enddo
          enddo
      enddo
c     ifiltr=0 spatial filter 
c     ifiltr=1 spectral  filter 
      if (ifiltr.eq.1) then
      call cut3d(pr,mij(1,1,1,n))
      else
      call filper(pr,mij(1,1,1,n))
      end if
            enddo
c
c     now mij contains (|S|Sij)^
c
            do n=1,6
              do k=1,n3m
          do j=1,n2m
      do i=1,n1m
      pr(i,j,k)=st(i,j,k,n)
              enddo
          enddo
      enddo
      if (ifiltr.eq.1) then
      call cut3d(pr,st(1,1,1,n))
      else
      call filper(pr,st(1,1,1,n))
      end if
            enddo
c
c     now st contains Sij^
c
              do k=1,n3m
          do j=1,n2m
      do i=1,n1m
      rhs(i,j,k)=2.*(st(i,j,k,1)**2+st(i,j,k,2)**2+st(i,j,k,3)**2)+
     1           4.*(st(i,j,k,4)**2+st(i,j,k,5)**2+st(i,j,k,6)**2)
      rhs(i,j,k)=sqrt(rhs(i,j,k)) 
              enddo
          enddo
      enddo
c
c     now rhs contains |S|^=sqrt(2*Sij^*Sij^)
c
      do 33 n=1,6
      do 33 k=1,n3m
      do 33 j=1,n2m
      do 33 i=1,n1m
      mij(i,j,k,n)=ell2c*rhs(i,j,k)*st(i,j,k,n)-ell1c*mij(i,j,k,n)
 33   continue
c
c   now mij=l2*|S|^Sij^-l1*(|S|Sij)^
c   l2 and l1  filter  grid and simulation grid
c   contraction of  mij with itself
c
      do 41 k=1,n3m
      do 41 j=1,n2m
      do 41 i=1,n1m
      pr(i,j,k)=(mij(i,j,k,1)**2+mij(i,j,k,2)**2+mij(i,j,k,3)**2)+
     1       2.*(mij(i,j,k,4)**2+mij(i,j,k,5)**2+mij(i,j,k,6)**2)
      rhs(i,j,k)=pr(i,j,k)
 41   continue
c     volume average (for isotropic turbulence)
      avmij=0.
      do k=1,n3m
       do j=1,n2m
        do i=1,n1m
          avmij=avmij+pr(i,j,k)
        end do
       end do
      end do
      avmij=avmij/float(n1m*n2m*n3m)
c
c     we have to contract lij with mij
c    apply the test cutoff to q
c    we have to define the q at the center
c     st(.,4)<------q1
c     st(.,5)<------q2
c     st(.,6)<------q3
      do 48 k=1,n3m
      kp=kpv(k)
      do 48 j=1,n2m
      jp=jpv(j)
      do 48 i=1,n1m
      ip=ipv(i)
      st(i,j,k,4)=.5*(q(1,ip,j,k)+q(1,i,j,k))
      st(i,j,k,5)=.5*(q(2,i,jp,k)+q(2,i,j,k))
      st(i,j,k,6)=.5*(q(3,i,j,kp)+q(3,i,j,k))
 48   continue
      if (ifiltr.eq.1) then
      call cut3d(st(1,1,1,4),st(1,1,1,1))        
      call cut3d(st(1,1,1,5),st(1,1,1,2))        
      call cut3d(st(1,1,1,6),st(1,1,1,3))        
      else
      call filper(st(1,1,1,4),st(1,1,1,1))        
      call filper(st(1,1,1,5),st(1,1,1,2))        
      call filper(st(1,1,1,6),st(1,1,1,3))        
      end if
c**************************** component 11 ********************** 
      do 50 k=1,n3m
      do 50 j=1,n2m
      do 50 i=1,n1m
      pr(i,j,k)=st(i,j,k,4)*st(i,j,k,4)
 50   continue
      if (ifiltr.eq.1) then
      call cut3d(pr,vis)        
      else
      call filper(pr,vis)        
      end if
c
c  g=(u_1u_1)^-u_1^u_1^)*m11
c
      do 52 k=1,n3m
      do 52 j=1,n2m
      do 52 i=1,n1m
      g(i,j,k)=(vis(i,j,k)-st(i,j,k,1)*st(i,j,k,1))*mij(i,j,k,1)
 52   continue
c**************************** component 22 ********************** 
      do 60 k=1,n3m
      do 60 j=1,n2m
      do 60 i=1,n1m
      pr(i,j,k)=st(i,j,k,5)*st(i,j,k,5)
 60   continue
      if (ifiltr.eq.1) then
      call cut3d(pr,vis)        
      else
      call filper(pr,vis)        
      end if
c
c  g=(u_2u_2)^-u_2^u_2^)*m22
c
      do 62 k=1,n3m
      do 62 j=1,n2m
      do 62 i=1,n1m
      g(i,j,k)=(vis(i,j,k)-st(i,j,k,2)*st(i,j,k,2))*mij(i,j,k,2)
     1          +g(i,j,k)        
 62   continue
c**************************** component 33 ********************** 
      do 70 k=1,n3m
      do 70 j=1,n2m
      do 70 i=1,n1m
      pr(i,j,k)=st(i,j,k,6)*st(i,j,k,6)
 70   continue
      if (ifiltr.eq.1) then
      call cut3d(pr,vis)        
      else
      call filper(pr,vis)        
      end if
c
c  g=(u_3u_3)^-u_3^u_3^)*m33
c
      do 72 k=1,n3m
      do 72 j=1,n2m
      do 72 i=1,n1m
      g(i,j,k)=(vis(i,j,k)-st(i,j,k,3)*st(i,j,k,3))*mij(i,j,k,3)
     1          +g(i,j,k)        
 72   continue
c**************************** component 12 ********************** 
      do 80 k=1,n3m
      do 80 j=1,n2m
      do 80 i=1,n1m
      pr(i,j,k)=st(i,j,k,4)*st(i,j,k,5)
 80   continue
      if (ifiltr.eq.1) then
      call cut3d(pr,vis)        
      else
      call filper(pr,vis)        
      end if
c
c  g=(u_1u_2)^-u_1^u_2^)*m12
c
      do 82 k=1,n3m
      do 82 j=1,n2m
      do 82 i=1,n1m
      g(i,j,k)=2.*(vis(i,j,k)-st(i,j,k,1)*st(i,j,k,2))*mij(i,j,k,4)
     1        +g(i,j,k)       
 82   continue
c**************************** component 13 ********************** 
      do 90 k=1,n3m
      do 90 j=1,n2m
      do 90 i=1,n1m
      pr(i,j,k)=st(i,j,k,4)*st(i,j,k,6)
 90   continue
      if (ifiltr.eq.1) then
      call cut3d(pr,vis)        
      else
      call filper(pr,vis)        
      end if
c
c  g=(u_1u_3)^-u_1^u_3^)*m13
c
      do 92 k=1,n3m
      do 92 j=1,n2m
      do 92 i=1,n1m
      g(i,j,k)=2.*(vis(i,j,k)-st(i,j,k,1)*st(i,j,k,3))*mij(i,j,k,5)
     1         +g(i,j,k)       
 92   continue
c**************************** component 23 ********************** 
      do 100 k=1,n3m
      do 100 j=1,n2m
      do 100 i=1,n1m
      pr(i,j,k)=st(i,j,k,5)*st(i,j,k,6)
100   continue
      if (ifiltr.eq.1) then
      call cut3d(pr,vis)        
      else
      call filper(pr,vis)        
      end if
c
c  g=(u_2u_3)^-u_2^u_3^)*m23
c
      do 102 k=1,n3m
      do 102 j=1,n2m
      do 102 i=1,n1m
      g(i,j,k)=2.*(vis(i,j,k)-st(i,j,k,2)*st(i,j,k,3))*mij(i,j,k,6)
     1         +g(i,j,k)       
102   continue
c     now g contains lij*mij
c      volume average (for isotropic turbulence)
      avg=0.
      do k=1,n3m
       do j=1,n2m
        do i=1,n1m
          avg=avg+g(i,j,k)
        end do
       end do
      end do
      avg=avg/float(n1m*n2m*n3m)
c       calculate smagorinsky constant
      cdynav=-.5*avg/avmij
c     print *,'time=',time,'C(volume average)=',cdynav
c    1       ,'C_s(volume average)=',sqrt(cdynav)
c
       nptot=n1m*n2m*n3m
       nback=0
       nif=0
c       ibox=1 media su tutto il bob C(t)
c       ibox=0 C(x,y,z,t) con viscosita' turbolenta>0
       if(ibox.eq.1) then
       visav=0.
        do k=1,n3m
          do j=1,n2m
            do i=1,n1m
            vistu=cdynav*ell1c*gold(i,j,k)
           visav=visav+vistu
            vis(i,j,k)=vistu+cvisc
          if(vis(i,j,k).lt.0.) then 
            vis(i,j,k)=0.
            nif=nif+1
          end if
            end do
          end do
        end do
       visav=visav/float(n1m*n2m*n3m)
       else
       visav=0.
        do 300 k=1,n3m
         do 300 j=1,n2m
          do 300 i=1,n1m
          cnew=-.5*g(i,j,k)/rhs(i,j,k)
          cist=cist+cnew
          vis(i,j,k)=cnew*ell1c*gold(i,j,k)+cvisc
          if(cnew.lt.0.) nback=nback+1
          if(cnew.gt.0.) nif=nif+1
          if(vis(i,j,k).lt.0.) then 
            vis(i,j,k)=0.
           visav=visav+vis(i,j,k)-cvisc
          end if
 300   continue
       visav=visav/float(n1m*n2m*n3m)
       cist=cist/float(n1m*n2m*n3m)
       flback=float(nback)/float(nptot)*100.
       flinf=float(nif)/float(nptot)*100.
       end if
      if(istrat.eq.1.and.igrad.eq.0) then
      do kc=1,n3m
      do jc=1,n2m
      do ic=1,n1m
      gold(ic,jc,kc)=vis(ic,jc,kc)
      enddo
      enddo
      enddo
                                     endif
       return
       end 
c
c
      subroutine filper(u,uf)
      include 'param.f'
c
c     calculates filtered function using a box filter in physical
c     space. Periodic box
c
      dimension u(m1,m2,m3),uf(m1,m2,m3)
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/dim/n1,n1m,n2,n2m,n3,n3m
      suv=1./64.
      do 10 kc=1,n3m
       kp=kpv(kc)
       km=kmv(kc)
      do 10 jc=1,n2m
       jp=jpv(jc)
       jm=jmv(jc)
      do 10 ic=1,n1m
       ip=ipv(ic)
       im=imv(ic)
      uf(ic,jc,kc)=suv*(
     1                 2.*u(im,jc,kp)+4.*u(ic,jc,kp)+2.*u(ip,jc,kp)+
     1                 4.*u(im,jc,kc)+8.*u(ic,jc,kc)+4.*u(ip,jc,kc)+
     1                 2.*u(im,jc,km)+4.*u(ic,jc,km)+2.*u(ip,jc,km)
     1                   +u(im,jm,kp)+2.*u(ic,jm,kp)+u(ip,jm,kp)+
     1                 2.*u(im,jm,kc)+4.*u(ic,jm,kc)+2.*u(ip,jm,kc)+
     1                    u(im,jm,km)+2.*u(ic,jm,km)+u(ip,jm,km) 
     1                   +u(im,jp,kp)+2.*u(ic,jp,kp)+u(ip,jp,kp)+
     1                 2.*u(im,jp,kc)+4.*u(ic,jp,kc)+2.*u(ip,jp,kc)+
     1                    u(im,jp,km)+2.*u(ic,jp,km)+u(ip,jp,km)
     1                 )
  10  continue
      return
      end
      subroutine cut3d(u,uf)
c
c     3d sharp Fourier cutoff filter
c     u is the resolved field
c     uf is the filtered field
c
      include 'param.f'
      dimension u(m1,m2,m3),uf(m1,m2,m3)
      complex qtil(m1,(m2-1)/2+1,m3)
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/fft001/ifxx1(13),trigxx1(2*(m1-1))
      common/fft002/ifxx2(13),trigxx2(3*(m2-1)/2+1)
      common/fft003/ifxx3(13),trigxx3(2*(m3-1))
      complex xa(m3-1,m1-1),wor(m3-1,m1-1)
      complex xa2(m1-1,m3-1),wor2(m1-1,m3-1)
      real xr(m2+1,m1-1),work(m2,m1-1)
c
      n2mh=n2m/2+1
      nxcut=n1m/2
      nycut=n2m/2
      nzcut=n3m/2

c
      do  10 k=1,n3m
c
        do i=1,n1m
         xr(1,i)=u(i,n2m,k)
         xr(n2m+2,i)=u(i,1,k)
         do j=1,n2m
          js=j+1
          xr(js,i)=u(i,j,k)
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
c
c     qtil is the velocity component in Fourier space.
c     Now apply a  sharp-Fourier cutoff
c     if kx>kxcut q=0.
c     if ky>kycut q=0.
c     if kz>kzcut q=0.
c
c     filter in z direction
c
      do k=nzcut/2+2,n3m-nzcut/2
       do j=1,n2mh
        do i=1,n1m
        qtil(i,j,k)=(0.,0.)
        end do
       end do
      end do
c
c     filter in x direction
c
      do k=1,n3m
       do j=1,n2mh
        do i=nxcut/2+2,n1m-nxcut/2
        qtil(i,j,k)=(0.,0.)
        end do
       end do
      end do
c
c     filter in y direction
c
      do k=1,n3m
       do j=nycut/2+2,n2mh
        do i=1,n1m
        qtil(i,j,k)=(0.,0.)
        end do
       end do
      end do

c
c   back to the physical space
c
c
c   2-d  cfft applied (twice) from
c   wave number to physical space
c
       do 30 j=1,n2mh
c
       do i=1,n1m
        do k=1,n3m
         xa(k,i)=qtil(i,j,k)
        enddo
       enddo
c
       call cfft99(xa,wor,trigxx3,ifxx3,1,m3-1,n3m,n1m,+1)
c
       do i=1,n1m
        do k=1,n3m
         qtil(i,j,k)=xa(k,i)
        enddo
       enddo
  30   continue
c
       do 40 j=1,n2mh
       do k=1,n3m
        do i=1,n1m
         xa2(i,k)=qtil(i,j,k)
        enddo
       enddo
       call cfft99(xa2,wor2,trigxx1,ifxx1,1,m1-1,n1m,n3m,+1)
c
       do i=1,n1m
        do k=1,n3m
         qtil(i,j,k)=xa2(i,k)
        enddo
       enddo
c
   40   continue
c
c   2-d  fft applied from
c   wave number to physical space
c
      do 50 k=1,n3m
c
       do j=1,n2mh
       jp=2*j
       jd=2*j-1
        do i=1,n1m
         xr(jd,i)=real(qtil(i,j,k))
         xr(jp,i)=aimag(qtil(i,j,k))
        enddo
       enddo
c
       call fft99(xr,work,trigxx2,ifxx2,1,m2+1,n2m,n1m,+1)
c
       do i=1,n1m
        do j=1,n2m
         js=j+1
         uf(i,j,k)=xr(js,i)
        enddo
       enddo
c
  50  continue
c
       return
       end 
