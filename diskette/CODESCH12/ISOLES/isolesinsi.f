c
c  ****************************** subrout inlero **********************
c
c   initial passive scalar conditions with spectrum as that by
c   Lesieur Rogallo
c
      subroutine inlero(rho,qtil4,qcap)
c
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/numond/kx(m1),ky(m2),kz(m3),kkmax
      dimension qcap(m1,m2,m3),rho(m1,m2,m3)
      complex qtil4(m1,(m2-1)/2+1,m3)
      complex cexp4
      common/strat1/istrat,rho0,g,schm
      common/strat2/igrad,bvais
      common/eispe/ekk(m1)
      common/sprh/srho(0:m3)
      common/kpmasp/kmax
      n2mh=n2m/2+1
      pi=4.*atan(1.)
       isd=19
       call srand(isd)
c      call ranset(isd)
       do k=1,n3m
        do j=1,n2mh
         do i=1,n1m
           if(i.eq.1.and.j.eq.1)then
            qtil4(i,j,k)=0.0
           else
            fk2=float(kx(i)*kx(i)+ky(j)*ky(j)+kz(k)*kz(k))
            fk=sqrt(fk2)
            kk=fk
            ampk=sqrt(ekk(kk)/(4.*pi*fk2))
c           x4=rnd()
            x4=rand()
c           x4=ranf()
            the4=2.*pi*x4
            cexp4=cmplx(cos(the4),sin(the4))
            qtil4(i,j,k)=ampk*cexp4
          endif
         end do
        end do
       end do
      call fft3d(rho,qtil4,+1)
       rhoma=-100.
       rhomi=100.
       do k=1,n3m
        do j=1,n2m
         do i=1,n1m
         rhoma=max(rhoma,rho(i,j,k))
         rhomi=min(rhoma,rho(i,j,k))
         end do
        end do
       end do
      call sperho(rho,qcap)
       eto=0.
       etro=0.
       do kk=1,kkmax
       eto=eto+ekk(kk)
       etro=etro+srho(kk)
       enddo
       write(6,*)'initial Lesieur Rogallo vel spec',eto,etro
       write(6,*)'initial Les Rog minmax pass ',rhoma,rhomi
      open(58,file='inico.plo')
      rewind(58)
       do k=1,kmax
       write (58,*)k,srho(k),ekk(k)
       end do
      close(58)
      return
      end
c
c  ****************************** subrout inchas **********************
c
c   initial passive scalar conditions with spectrum equal to that
c   of Chasnov  
c
      subroutine inchas(rho,qtil4,qcap,q)
c
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/numond/kx(m1),ky(m2),kz(m3),kkmax
      dimension qcap(m1,m2,m3),rho(m1,m2,m3)
      dimension q(ndv,m1,m2,m3)
      complex qtil4(m1,(m2-1)/2+1,m3)
      complex cexp4,cexp2
      common/strat1/istrat,rho0,g,schm
      common/spemw/akkpp,qq,sig
      common/strat2/igrad,bvais
      dimension ekk(m1)
      common/sprh/srho(0:m3)
      common/kpmasp/kmax
      n2mh=n2m/2+1
      pi=4.*atan(1.)
c
c  here we try  the initial spectrum of  Chasnov (JFM 235)
c
       isig=sig
       facn=1.
       lf=isig/2+1
       do ll=2,lf
       facn=facn*((ll-1)*2-1)
       enddo
       anu=isig**((isig+1)*0.5)
       akrap=akkpp
       aa1=qq*(2./pi)**0.5*anu/(facn*(akrap))
       write(6,*)'isig, facn, akrap, aa1 ',isig,facn,akrap,aa1
       eto=0.
       do kk=1,kkmax
       ak=kk
       rk=ak/akrap
       ek1kl=aa1*rk**sig*exp(-sig/2.*rk**2)
       ekk(kk)=ek1kl
       eto=eto+ekk(kk)
       enddo
       write(6,*)'   akpp=',akkpp,'   aa1=',aa1
       write(6,*)'   tt=',qq,'   rhoto=',eto

       isd=19
       call srand(isd)
c      call ranset(isd)
       do k=1,n3m
        do j=2,n2mh
         do i=1,n1m
            fk2=float(kx(i)*kx(i)+ky(j)*ky(j)+kz(k)*kz(k))
            fk=sqrt(fk2)
            kk=fk
            ampk=sqrt(ekk(kk)/(4.*pi*fk2))
c           x4=rnd()
            x4=rand()
c           x4=ranf()
            fke=float(kx(i)+ky(j)+kz(k))
            cexp2=cmplx(cos(fke),sin(fke))
            the4=2.*pi*x4
            cexp4=cmplx(cos(the4),sin(the4))
            qtil4(i,j,k)=ampk*cexp4*cexp2
         end do
        end do
       end do
       j=1
       do k=1,n3m
         do i=2,n1m
            fk2=float(kx(i)*kx(i)+ky(j)*ky(j)+kz(k)*kz(k))
            fk=sqrt(fk2)
            kk=fk
            ampk=sqrt(ekk(kk)/(4.*pi*fk2))
c           x4=rnd()
            x4=rand()
c           x4=ranf()
            fke=float(kx(i)+ky(j)+kz(k))
            cexp2=cmplx(cos(fke),sin(fke))
            the4=2.*pi*x4
            cexp4=cmplx(cos(the4),sin(the4))
            qtil4(i,j,k)=ampk*cexp4*cexp2
         end do
        end do
      call fft3d(rho,qtil4,+1)
       rhoma=-100.
       rhomi=100.
       do k=1,n3m
        do j=1,n2m
         do i=1,n1m
         rhoma=max(rhoma,rho(i,j,k))
         rhomi=min(rhoma,rho(i,j,k))
         end do
        end do
       end do
      call sperho(rho,qcap)
       eto=0.
       etro=0.
       do kk=1,kkmax
       eto=eto+ekk(kk)
       etro=etro+srho(kk)
       enddo
       write(6,*)'initial Chasnov dens spec',eto,etro, 
     1    'the velocity is now set = 0'
       write(6,*)'initial Les Rog minmax pass ',rhoma,rhomi
      open(58,file='inico.plo')
      rewind(58)
       do k=1,kmax
       ak=float(k)/akkpp
       write (58,*)ak,srho(k),ekk(k)
       end do
      close(58)
      do k=1,n3m
      do j=1,n2m
      do i=1,n1m
       q(1,i,j,k) =0.
       q(2,i,j,k) =0.
       q(3,i,j,k) =0.
      enddo
      enddo
      enddo
      return
      end
