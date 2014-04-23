c***********************************************
c     dynamic subgrid-scale stress model for
c     flows with density variations
c     This is the dynamic model with constant eddy viscosity and diffusivity
c***********************************************
      subroutine dywolil(q,pr,rho,ntime,time)
      include 'param.f'
      dimension q(ndv,m1,m2,m3),pr(m1,m2,m3)
      dimension rho(m1,m2,m3)
      common/strain/st(m1,m2,m3,6)
      common/mijdy/mij(m1,m2,m3,6)
      real*4 mij
      common/qdyn/g(m1,m2,m3)
      common/rhsc/rhs(m1,m2,m3)
      common/qles1/deltax1,deltay1,deltaz1,ell1c
      common/qles2/deltax2,deltay2,deltaz2,ell2c
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      common/sma/vis(m1,m2,m3)
      common/smapra/gold(m1,m2,m3)
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/ledat/ics0,cvisc,ifiltr,ibox
      common/strat1/istrat,rho0,gg,schm
      common/tstep/dt,beta,ren
      common/dycopr/cdynav,visav,avg,avmij,nif,cist,flback,flint
c      
c
c      computation of Sij=0.5*(dui/dxj + duj/dxi)
c
      call strper(q)
c
c     pr=Sij
c
            do n=1,6
              do k=1,n3m
          do j=1,n2m
      do i=1,n1m
      pr(i,j,k)=st(i,j,k,n)
              enddo
          enddo
      enddo
c
c    filtering in physical space
c
      call filper(pr,st(1,1,1,n))
            enddo
c
c     now st contains Sij^
c
              do k=1,n3m
          do j=1,n2m
      do i=1,n1m
      rhs(i,j,k)=2.*(st(i,j,k,1)**2+st(i,j,k,2)**2+st(i,j,k,3)**2)+
     1           4.*(st(i,j,k,4)**2+st(i,j,k,5)**2+st(i,j,k,6)**2)
      rhs(i,j,k)=rhs(i,j,k)*(1.-(ell2c/ell1c)**(2./3.))
              enddo
          enddo
      enddo
c
c     now Sij^ in mij to contract Lij
c
      do n=1,6
      do k=1,n3m
      do j=1,n2m
      do i=1,n1m
      mij(i,j,k,n)=st(i,j,k,n)
              enddo
          enddo
      enddo
      enddo

c
c  remember ell*=(dx1*dx2*dx3)^2/3
c
c  ell2  with 2*dx_i    ell1 with dx_i
c
c     now rhs contains  the denominator of eddy viscosity
c            |S|^=(2*Sij^*Sij^)
c
c
c     volume average (for isotropic turbulence)
c
      avmij=0.
      do k=1,n3m
       do j=1,n2m
        do i=1,n1m
          avmij=avmij+rhs(i,j,k)
        end do
       end do
      end do
      avmij=avmij/float(n1m*n2m*n3m)
      denvis=avmij
c     print *,'avmij=',avmij
c     we have to contract lij with sij
c    apply the test cutoff to q
c    we have to define the q at the center
c     st(.,4)<------q1
c     st(.,5)<------q2
c     st(.,6)<------q3
      do k=1,n3m
      kp=kpv(k)
      do j=1,n2m
      jp=jpv(j)
      do i=1,n1m
      ip=ipv(i)
      st(i,j,k,4)=.5*(q(1,ip,j,k)+q(1,i,j,k))
      st(i,j,k,5)=.5*(q(2,i,jp,k)+q(2,i,j,k))
      st(i,j,k,6)=.5*(q(3,i,j,kp)+q(3,i,j,k))
        end do
       end do
      end do
      call filper(st(1,1,1,4),st(1,1,1,1))        
      call filper(st(1,1,1,5),st(1,1,1,2))        
      call filper(st(1,1,1,6),st(1,1,1,3))        
c**************************** component 11 ********************** 
      do k=1,n3m
      do j=1,n2m
      do i=1,n1m
      pr(i,j,k)=st(i,j,k,4)*st(i,j,k,4)
        end do
       end do
      end do
      call filper(pr,vis)        
c
c  g=(u_1u_1)^-u_1^u_1^)*s11
c
      do k=1,n3m
      do j=1,n2m
      do i=1,n1m
      g(i,j,k)=(vis(i,j,k)-st(i,j,k,1)*st(i,j,k,1))*mij(i,j,k,1)
        end do
       end do
      end do
c**************************** component 22 ********************** 
      do k=1,n3m
      do j=1,n2m
      do i=1,n1m
      pr(i,j,k)=st(i,j,k,5)*st(i,j,k,5)
        end do
       end do
      end do
      call filper(pr,vis)        
c
c  g=(u_2u_2)^-u_2^u_2^)*m22
c
      do k=1,n3m
      do j=1,n2m
      do i=1,n1m
      g(i,j,k)=(vis(i,j,k)-st(i,j,k,2)*st(i,j,k,2))*mij(i,j,k,2)
     1          +g(i,j,k)        
        end do
       end do
      end do
c**************************** component 33 ********************** 
      do k=1,n3m
      do j=1,n2m
      do i=1,n1m
      pr(i,j,k)=st(i,j,k,6)*st(i,j,k,6)
        end do
       end do
      end do
      call filper(pr,vis)        
c
c  g=(u_3u_3)^-u_3^u_3^)*m33
c
      do k=1,n3m
      do j=1,n2m
      do i=1,n1m
      g(i,j,k)=(vis(i,j,k)-st(i,j,k,3)*st(i,j,k,3))*mij(i,j,k,3)
     1          +g(i,j,k)        
        end do
       end do
      end do
c**************************** component 12 ********************** 
      do k=1,n3m
      do j=1,n2m
      do i=1,n1m
      pr(i,j,k)=st(i,j,k,4)*st(i,j,k,5)
        end do
       end do
      end do
      call filper(pr,vis)        
c
c  g=(u_1u_2)^-u_1^u_2^)*m12
c
      do k=1,n3m
      do j=1,n2m
      do i=1,n1m
      g(i,j,k)=2.*(vis(i,j,k)-st(i,j,k,1)*st(i,j,k,2))*mij(i,j,k,4)
     1        +g(i,j,k)       
        end do
       end do
      end do
c**************************** component 13 ********************** 
      do k=1,n3m
      do j=1,n2m
      do i=1,n1m
      pr(i,j,k)=st(i,j,k,4)*st(i,j,k,6)
        end do
       end do
      end do
      call filper(pr,vis)        
c
c  g=(u_1u_3)^-u_1^u_3^)*m13
c
      do k=1,n3m
      do j=1,n2m
      do i=1,n1m
      g(i,j,k)=2.*(vis(i,j,k)-st(i,j,k,1)*st(i,j,k,3))*mij(i,j,k,5)
     1         +g(i,j,k)       
        end do
       end do
      end do
c**************************** component 23 ********************** 
      do k=1,n3m
      do j=1,n2m
      do i=1,n1m
      pr(i,j,k)=st(i,j,k,5)*st(i,j,k,6)
        end do
       end do
      end do
      call filper(pr,vis)        
c
c  g=(u_2u_3)^-u_2^u_3^)*m23
c
      do k=1,n3m
      do j=1,n2m
      do i=1,n1m
      g(i,j,k)=2.*(vis(i,j,k)-st(i,j,k,2)*st(i,j,k,3))*mij(i,j,k,6)
     1         +g(i,j,k)       
        end do
       end do
      end do
c     now g contains lij*Sij^
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
      anuvis=avg
c
c
c   the turbulent Prandtl number is now calculated
c
c
c      computation of heat fluxes Fi= drho/dxi)
c      F1=m(1) , F2=m(2)  , F3=m(3)
c
      do k=1,n3m
          do j=1,n2m
            do i=1,n1m
      mij(i,j,k,1)=((rho(ipv(i),j,kpv(k))-rho(imv(i),j,kpv(k)))*dx1*0.5
     1             +(rho(ipv(i),j,k)-rho(imv(i),j,k))*dx1*0.5 )*0.5
      mij(i,j,k,2)=((rho(i,jpv(j),kpv(k))-rho(i,jmv(j),kpv(k)))*dx2*0.5
     1             +(rho(i,jpv(j),k)-rho(i,jmv(j),k))*dx2*0.5 )*0.5
      mij(i,j,k,3)=(rho(i,j,kpv(k))-rho(i,j,k))*dx3
            end do
          end do
        end do
            do n=1,3
              do k=1,n3m
          do j=1,n2m
      do i=1,n1m
      pr(i,j,k)=mij(i,j,k,n)
              enddo
          enddo
      enddo
      call filper(pr,mij(1,1,1,n))
            enddo
              do k=1,n3m
          do j=1,n2m
      do i=1,n1m
      rhs(i,j,k)=(mij(i,j,k,1)**2+mij(i,j,k,2)**2+mij(i,j,k,3)**2)
      rhs(i,j,k)=rhs(i,j,k)*(1.-(ell2c/ell1c)**(2./3.))
              enddo
          enddo
      enddo
c
c     volume average denominator (for isotropic turbulence)
c
c       calculate rate of dissipation   constant
      if(avmij.eq.0.) then
      cdydis=0.
                      else
      cdydis=avg/avmij
                      endif
      avmij=0.
      do k=1,n3m
       do j=1,n2m
        do i=1,n1m
          avmij=avmij+rhs(i,j,k)
        end do
       end do
      end do
      avmij=avmij/float(n1m*n2m*n3m)
      dendif=avmij
c
c     now mij(n contains F_n^
c
c    Rememmber that
c   in st(.,4) q1  ,st(.,5) q2  ,st(.,6) q3
c   in st(.,1) q1^ ,st(.,2) q2^ ,st(.,3) q3^
c
c   the filter is applied to rho at the cell center
c
      do k=1,n3m
          do j=1,n2m
            do i=1,n1m
      rhs(i,j,k)=(rho(i,j,kpv(k))+rho(i,j,k))*0.5
            end do
          end do
        end do
      call filper(rhs,g)        
c
c   now rhs contains rho   and g rho^
c**************************** component u1rho ********************** 
      do k=1,n3m
      do j=1,n2m
      do i=1,n1m
      pr(i,j,k)=st(i,j,k,4)*rhs(i,j,k)
        end do
       end do
      end do
      call filper(pr,vis)        
c
c  gold=(u_1rho)^-u_1^rho^)*F_1^  
c
      do k=1,n3m
      do j=1,n2m
      do i=1,n1m
      gold(i,j,k)=(vis(i,j,k)-st(i,j,k,1)*g(i,j,k))*mij(i,j,k,1)
        end do
       end do
      end do
c**************************** component u2rho ********************** 
      do k=1,n3m
      do j=1,n2m
      do i=1,n1m
      pr(i,j,k)=st(i,j,k,5)*rhs(i,j,k)
        end do
       end do
      end do
      call filper(pr,vis)        
c
c  gold=(u_2rho)^-u_2^rho^)*F_2^
c
      do k=1,n3m
      do j=1,n2m
      do i=1,n1m
      gold(i,j,k)=(vis(i,j,k)-st(i,j,k,2)*g(i,j,k))*mij(i,j,k,2)
     1         +gold(i,j,k)       
        end do
       end do
      end do
c**************************** component u3rho ********************** 
      do k=1,n3m
      do j=1,n2m
      do i=1,n1m
      pr(i,j,k)=st(i,j,k,6)*rhs(i,j,k)
        end do
       end do
      end do
      call filper(pr,vis)        
c
c  gold=(rhou_3)^-rho^u_3^)*m13
c
      do k=1,n3m
      do j=1,n2m
      do i=1,n1m
      gold(i,j,k)=(vis(i,j,k)-st(i,j,k,3)*g(i,j,k))*mij(i,j,k,3)
     1         +gold(i,j,k)       
        end do
       end do
      end do
c     now gold contains R_ri*F_i^
c      volume average (for isotropic turbulence)
      avg=0.
      do k=1,n3m
       do j=1,n2m
        do i=1,n1m
          avg=avg+gold(i,j,k)
        end do
       end do
      end do
      avg=avg/float(n1m*n2m*n3m)
      anudif=avg
      cdydif=anudif/dendif
       visav=0.
       difav=0.
        do k=1,n3m
          do j=1,n2m
            do i=1,n1m
            vistu=cdydis
           visav=visav+vistu
            vis(i,j,k)=vistu+cvisc
          if(vis(i,j,k).lt.0.) then 
            vis(i,j,k)=cvisc
            nif=nif+1
          end if
            diftu=cdydif
           difav=difav+diftu
            gold(i,j,k)=diftu+cvisc/schm
          if(gold(i,j,k).lt.0.) then 
            gold(i,j,k)=cvisc
          end if
            end do
          end do
        end do
       visav=visav/float(n1m*n2m*n3m)
       difav=difav/float(n1m*n2m*n3m)
       write(6,161)time,denvis,anuvis,dendif,anudif,visav,difav,cvisc
  161 format(2x,'visc dyn',8e12.4)
            do n=1,6
              do k=1,n3m
          do j=1,n2m
      do i=1,n1m
      st(i,j,k,n)=0.
              enddo
          enddo
      enddo
            enddo
       return
       end 
