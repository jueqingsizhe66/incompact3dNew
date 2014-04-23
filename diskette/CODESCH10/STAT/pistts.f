c************************************************************************
c                                                                       *
c                                                                       *
c     this code perform a postprocessing of the fields generated        *
c     by the code in the directory CODESCH10/PIP  
c     the simulation was done  with
c     All variables in a staggered grid:                 *
c                                                                       *
c        flux variables were introduced 
c        q1=vtheta*r, q2= r*vr, q3= vz                                    *      
c                                                                       *
c                                                                       *
c************************************************************************
      program main                                                      
      include 'param.f'
      common/timavg/timav
      common/timini/timei
      common/npjet/n2t
      common/nbal/nba,ibudg
      open(15,file='pipstat.d')
  101 format(a38)
      read(15,301) dummy                                                
      read(15,*) n1,n2,n3,nsst,nwrit,nread
      read(15,301) dummy                                                
      read(15,*) n1p,n2p,n3p,nprde                                      
      read(15,301) dummy                                                
      read(15,*) ntst,nprint,npin,ireset ,timei,npouth
      read(15,301) dummy                                                
      read(15,*) alx3p                                                  
      read(15,301) dummy                                                
      read(15,*) re,vper,dt
      read(15,301) dummy                                                
      read(15,*) r0
      read(15,301) dummy                                                
      read(15,*) cflc,scfl,t0med1,t0pres                       
      read(15,301) dummy                                                
      read(15,*) ros                                                    
      read(15,301) dummy                                                
      read(15,*)nini,nfin,nstri,irid 
c$$$$$ parameters for non uniform grid (r distribution) $$$$$$$$$$$$$$$$$
      read(15,301) dummy                                                
      read(15,*)strr,rext,rint,rmed1
      read(15,301) dummy                                                
      read(15,*)istr,rmed,etdp,strb,n2t
c$$$$$$$parameters for boundary conditions $$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      read(15,301) dummy
      read(15,*) iav     
      read(15,301) dummy
      read(15,*) ialvo,icorspe,timav,ibudg
      read(15,301) dummy
      read(15,*) ichrc
      read(15,301) dummy
      read(15,*)strro,rexto,rinto,rmed1o
      read(15,301) dummy
      read(15,*)istro,rmedo,etdpo,strbo,n2to
301   format(a4)                                                        
      icfl=scfl
      ime1t0=t0med1
      ipr0=t0pres
c
c  when newvar=1 changes initial passive scalar mean profile
c                 
c
      pi=2.*asin(1.)                                                    
      alx3d=2.*pi*alx3p
      tfini=dt*ntst                                                     
      tfin=tfini                                                        
      n1m=n1-1                                                          
      n2m=n2-1                                                          
      n3m=n3-1                                                          
      n3mh=n3m/2+1                                                      
      if(n1m.eq.1) then
      iaxsy = 0
       else
      iaxsy=1
       endif
c                                                                       
c                                                                       
      alx3=alx3d                                                        
c                                                                       
      call openfi
c                                                                       
c     assign coefficients for time marching schemes                     
c
      write(6,*)'*****************************************************'
      write(32,*)'*****************************************************'
      write(6,*)'*                                                    *'
      write(32,*) '*                                                  *'
      write(6,*)'*               TURBULENT PIPE WITH ROTATION         *'
      write(32,*)'*              TURBULENT PIPE WITH ROTATION         *'
      write(6,*)'*      staistics from fields                         *'
      write(32,*)'*     staistics from fields                         *'
      write(6,*)'*****************************************************'
      write(32,*)'*****************************************************'
      write(6,*)'  '
      write(32,*) '  '
      write(6,*)'  '
      write(32,*) '  '
      write(6,112)alx3d,r0                                              
      write(32,112)alx3d,r0                                             
  112 format(10x,'domain dimensions L_z =',f8.5,' L_r = ',f8.5)                
      write(6,*)'  '
      write(32,*) '  '
      write(6,*)'  '
      write(32,*) '  '
c                                                                       
      call gcurv                                                        
                                                                        
      stop                                                              
      end                                                               
c************************************************************************
c************************************************************************
c                                                                       *
c                                                                       * 
c     q1=v(theta)*r  q2=v(r)*r     q3=v(zeta)                           *
c                                                                       *
c                                                                       *
c************************************************************************
c************************************************************************
      subroutine gcurv                                                  
      include 'param.f'
      common/tima3/tiax3d
      common/timavg/timav
      common/timini/timei
      common/nbal/nba,ibudg
      common/iolfil/iovco,iohel
c
c                                                                       
c     grid definition, indices and mesh size calculation                          
c                                                                       
      call meshes
      call indic                                                        
      call cordin                                                       
c
c
c
      read(15,301) dummy
301   format(a4)                                                        
       read(15,*) i3dou,njumk,njumj
      read(15,301) dummy
       read(15,*) npq,(jprq(n),n=1,npq)
       write(6,*) 'njumk,  njumj ',njumk,njumj 
       do n=1,npq
       write(6,*)'n,jprq(n) ',n,jprq(n)
       enddo
       ioldf=0
       iovco=0
       iohel=0
c
c*******************************************************
c
c     print some informations on the run
c
      write(6,*) '  '
      write(32,*) '  '
      write(6,754)n1,n2,n3                                              
      write(32,754)n1,n2,n3                                             
  754 format(10x,'number of grid points :'/                             
     1      5x,'n1=',i4,2x,'n2=',i4,2x,'n3=',i4/)                       
      write(6,*) '  '
      write(32,*) '  '
      write(6,756)n1p,n2p,n3p                                              
      write(32,756)n1p,n2p,n3p                                             
  756 format(10x,'plotting stride  :'/                             
     1      5x,'n1p=',i4,2x,'n2p=',i4,2x,'n3p=',i4/)                       
      write(6,*) '  '
      write(32,*) '  '
      write(6,755) re,ros
      write(32,755) re,ros
  755 format(3x,' Parameters of the flow: ',/,
     1 ' Reynolds number = ',e10.3,3x,' Rossby number = ',e10.3) 
      write(6,*) '  '
      write(32,*) '  '
      ren=re                                                            
      time=0.                                                           
      ntii=0                                                            
      beta=dt/re*0.5                                                    
            nap=0                                                           
            nav=0                                                           
            nvv=0                                                           
            nba=0
c
c  evaluation of metric quantities for the inversion
c
      call coetar
      call fftqua
c
c  here the initial intervals for the pdf calcualtion are assigned
c
      call pdfini
  159 format(1x,i4,2x,e10.4,3e10.3,3(1x,e10.4,1x,i3,1x,i3),e10.3)
c
      ncount=0
      ntii=1
       time=timei
      write(6,*)'ntii,ntst',ntii,ntst
      do 350 ntime=ntii,ntst                                           
c                                                                       
c     in tscheme the fields are read 
c     
      call tschem(ntime,time,ncount)
       ncount=ncount+1
c                                                                       
        ntt=ntt+1
c                                                                       
c*******************************************************
c                                                                       
c     write the flow field                                              
c                                                                       
c*******************************************************
           call outh(time,nav,ntime,cflm,nvv,navbu)            
         if(ntime.eq.ntst) go to 351
        time=time+dt
  350 continue                                                          
  351 continue                                                          
c
c  here the pdf are written
c  The pdf for the pipe are not reported in the book
c
      call pdffin(time)
      return                                                            
      end                                                               
c************************************************************************
c                                                                       *
c ****************************** subrout coetar  **********************  *
c                                                                       *
c    this subroutine calculates the coefficients for the              *
c    scretization of the viscous terms of the equation for q3 are evaluated
c    these coefficients are required to evaluated the pressure gradient 
c    and hence the friction velocity
c                                                                       *
c************************************************************************
      subroutine coetar
      include 'param.f'
      common/cor3j/ap3j(m2),ac3j(m2),am3j(m2)
c
c  *******   coefficients in several eq. funct of j
c
      do jc=1,n2m
      udx1q(jc)=dx1q/rm(jc)**2
      volz(jc)=rm(jc)*g2rm(jc)/(dx1*dx2*dx3)
      enddo
c
c  ***********  coefficients for q3   inner points
c
      do jc=2,n2m-1
      jp=jc+1
      a22=dx2q/g2rm(jc)/rm(jc)
      a22p= +a22*rc(jp)/g2rc(jp)
      a22m= +a22*rc(jc)/g2rc(jc)
      ap3j(jc)=a22p
      am3j(jc)=a22m
      ac3j(jc)=-(a22p+a22m)
      enddo
c
c    r=0 gives the following b.c. at axis  equiv. to dq3/dr=0
c
      jc=1
      jp=jc+1
      ugmm2=dx2q/g2rm(jc)/rm(jc)
      am3j(jc)=0.
      ac3j(jc)=ugmm2*rc(jp)/g2rc(jp)
      ap3j(jc)=(rc(jp)/g2rc(jp))*ugmm2
c
c    q3=0 has been assumed at the wall boundary
c
      jc=n2m
      jp=jc+1
      ugmm2=dx2q/g2rm(jc)/rm(jc)
      am3j(jc)=rc(jc)/g2rc(jc)*ugmm2
      ac3j(jc)=ugmm2*rc(jc)/g2rc(jc)
      ap3j(jc)=ugmm2*rc(jp)/g2rc(jp)*2.
      return
      end
c************************************************************************
c                                                                       *
c ****************************** subrout prgqso  ********************** *
c  this subroutine performs the calculation of pressure gradient and    *
c                                                                       *
c************************************************************************
      subroutine prgqso              
      include 'param.f'
      common/cor3j/ap3j(m2),ac3j(m2),am3j(m2)
      pi=2.*asin(1.)
      alre=1./ren
c
c
c    Pressure gradient
c
      s3tot=0.
      do kc=1,n3m
      km=kmv(kc)
      kp=kpv(kc)
            do jc=1,n2m
                  do ic=1,n1m
                  im=imv(ic)
                  ip=ipv(ic)
c
c   11 second derivatives of q3
c
      dq31=(q3(ip,jc,kc)
     1     -q3(ic,jc,kc)*2.
     1     +q3(im,jc,kc))*udx1q(jc)
c
c   33 second derivatives of q3
c
      dq33=(q3(ic,jc,kp)
     1     -q3(ic,jc,kc)*2.
     1     +q3(ic,jc,km))*dx3q
      dcq3=dq31+dq33
      s3tot=s3tot+dcq3*volz(jc)/ren
                  enddo

            enddo
      enddo
c     c
c   add second derivatie in r
c
      do kc=1,n3m
            do jc=2,n2m-1
            jm=jmv(jc)
            jp=jpv(jc)
                  do ic=1,n1m
c
c   22 second derivatives of q3
c
      dq32= q3(ic,jp,kc)*ap3j(jc)
     1     -q3(ic,jc,kc)*(ap3j(jc)+am3j(jc))
     1     +q3(ic,jm,kc)*am3j(jc)
      s3tot=s3tot+dq32*volz(jc)/ren
                  enddo
            enddo
      enddo
c
c   22 second derivatives of q3  at r=n2m
c
      jc=n2m
      jm=jc-1
      do kc=1,n3m
                  do ic=1,n1m
c
c   22 second derivatives of q3
c
      dq32=-q3(ic,jc,kc)*(ap3j(jc)+ac3j(jc))
     1     +q3(ic,jm,kc)*am3j(jc)
      s3tot=s3tot+dq32*volz(jc)/ren
                  enddo

      enddo
c     c
c   add second derivatie in r at r=0
c
      jc=1
      jp=jc+1
      do kc=1,n3m
                  do ic=1,n1m
c
c   22 second derivatives of q3
c
      dq32= q3(ic,jp,kc)*ap3j(jc)
     1     -q3(ic,jc,kc)*ap3j(jc)
      s3tot=s3tot+dq32*volz(jc)/ren
                  enddo
      enddo
       dp3ns=s3tot/(pi*alx3d)
      return
      end
c                                                                       
c  **************  subrout tschem                                       
c                                                                       
      subroutine tschem(ntime,time,ncount)
      include 'param.f'
      parameter (m3mh=m3m/2+1,m3p=m3+1)
      common/spck/rsp(m1m,m3mh),asp(m1m,m3mh),esp(m1m,m3mh)
      dimension rspv(6,m1m,m3mh),ak(m3mh),ai(m1m)
      common/timavg/timav
      common/timini/timei
      common/nbal/nba,ibudg
      dimension e1t(m2),e2t(m2),e3t(m2)
      dimension e1stt(m2),e2stt(m2),e3stt(m2)
      dimension e1stz(m2),e2stz(m2),e3stz(m2)
      dimension ene1(m2,m1m),ene2(m2,m1m),ene3(m2,m1m)
      dimension enep(m2,m1m)
      common/wavin/dlx1,dlx3,dkk1,dkk3
      common/ruudi/ruuthd(m2,m1m),ruuaxd(m2,m3m)

c
c   read the fields    
c
          call inirea(ntime,time,ntt,ncount,nap)           
c
c  evaluates the divergence
c
          call divgck(dmax,dtot)                                 
c         write(6,900)dmax,dtot                 ,ntime
  900   format(3x,'maxima local and global divergence of the read field',
     1     /,' dmax = ',e11.4,' dtot = ',e11.4,'ntime=',i5)
c
c   evaluates the pressure gradient
c 
          call prgqso
c
c   evaluates the mean velocity profile averaging in x1 and x3
c
      call velc
c
c   evaluates the mean voticity profile averaging in x1 and x3
c
      call vorc
      if(icorspe.eq.1) then
c
c   evaluates the azimuthal and axial
c   one-dimensional energy spectra in the
c   even the correlations in the azimuthal and in the
c   axial directions are calculated
c   remembers that here the correlations are evaluated
c   in the physiacal space. Very costly
c
      nba=nba+1
      dvotz=1./float(n1m*n3m)
      n2rm=n2m
      n2r=n2
      do j=1,n2m
      e1t(j)=0.
      do i=1,n1m
      do k=1,n3m
      rhs(i,j,k)=q1(i,j,k)/rm(j)-vmepo(1,j)
      e1t(j)=e1t(j)+rhs(i,j,k)**2
      enddo
      enddo
      e1t(j)=e1t(j)*dvotz
      enddo
c
c   axial spectrum  q1 velocity
c
           call ruuaxi(n2rm)
           call correl(n2rm)
c
c   q1 velocity correlation
c
             do jc=1,n2m
      do k=1,n3m
      r11axi(jc,k)=r11axi(jc,k)+ruuaxd(jc,k)
      enddo
              enddo
      do jc=1,n2m
      e1stz(jc)=0.
      do k=1,k3max
      e1stz(jc)=e1stz(jc)+enez(jc,k)*dkk3
      enddo
      enddo
      do jc=1,n2m
      e1stzo(jc)=e1stzo(jc)+e1stz(jc)
      do k=1,k3max
      ene1zo(jc,k)=ene1zo(jc,k)+enez(jc,k)
      enddo
      enddo
c
c   azimuthal spectrum   q1 velocity
c
           call ruuthe(n2rm)
             do jc=1,n2m
      do i=1,n1m
      r11the(jc,i)=r11the(jc,i)+ruuthd(jc,i)
      enddo
              enddo
      do jc=1,n2m
      e1stt(jc)=0.
      do k=1,k1max
      e1stt(jc)=e1stt(jc)+enet(jc,k)*dkk1
      enddo
      enddo
      do jc=1,n2m
      e1to(jc)=e1to(jc)+e1t(jc)
      e1stto(jc)=e1stto(jc)+e1stt(jc)
      do k=1,k1max
      ene1to(jc,k)=ene1to(jc,k)+enet(jc,k)
      enddo
      enddo
c
c   radial velocity
c
      do j=1,n2m
      e2t(j)=0.
      jp=j+1
      do k=1,n3m
      do i=1,n1m
      if(j.eq.1) then
      q2pos= (q2(i,jp,k) - q2(isym(i),jp,k))*0.5/rc(jp)
       else
      q2pos= q2(i,j,k)/rc(j)
       endif
      rhs(i,j,k)=q2pos-vmepo(2,j)
      q2med=(q2pos-vmepo(2,j)+q2(i,jp,k)/rc(jp)-vmepo(2,jp))*0.5
      e2t(j)=e2t(j)+q2med**2
      enddo
      enddo
      e2t(j)=e2t(j)*dvotz
      enddo
c
c   axial spectrum  q2 velocity
c
           call ruuaxi(n2rm)
           call correl(n2rm)
      n3mh=n3m/2+1
             jc=n2
      do k=1,n3m
      ruuaxd(jc,k)=0.
      enddo
             do jc=1,n2m
      do k=1,n3m
      ruuaxm=(ruuaxd(jc,k)+ruuaxd(jc+1,k))*0.5
      r22axi(jc,k)=r22axi(jc,k)+ruuaxm
      enddo
              enddo
      do jc=1,n2m
      e2stz(jc)=0.
      do k=1,k3max
      ene2m=(enez(jc,k)+enez(jc+1,k))*0.5
      e2stz(jc)=e2stz(jc)+ene2m*dkk3
      enddo
      enddo
      do jc=1,n2m
      e2stzo(jc)=e2stzo(jc)+e2stz(jc)
      do k=1,k3max
      ene2m=(enez(jc,k)+enez(jc+1,k))*0.5*2.
      ene2zo(jc,k)=ene2zo(jc,k)+ene2m
      enddo
      enddo
c
c   azimuthal spectrum   q2 velocity
c
c
c   q2 velocity correlation
c
c
           call ruuthe(n2rm)
             jc=n2
      do i=1,n1m
      ruuthd(jc,i)=0.
      enddo
             do jc=1,n2m
      do i=1,n1m
      ruuthm=(ruuthd(jc,i)+ruuthd(jc+1,i))*0.5
      r22the(jc,i)=r22the(jc,i)+ruuthm
      enddo
              enddo
      do k=1,k1max
      enet(n2,k)=0.
      enddo
      do jc=1,n2m
      e2stt(jc)=0.
      do k=1,k1max
      ene2m=(enet(jc,k)+enet(jc+1,k))*0.5
      e2stt(jc)=e2stt(jc)+ene2m*dkk1
      enddo
      enddo
      do jc=1,n2m
      e2to(jc)=e2to(jc)+e2t(jc)
      e2stto(jc)=e2stto(jc)+e2stt(jc)
      do k=1,k1max
      ene2m=(enet(jc,k)+enet(jc+1,k))*0.5
      ene2to(jc,k)=ene2to(jc,k)+ene2m
      enddo
      enddo
c
c    axial velocity q3
c
      do j=1,n2m
      e3t(j)=0.
      do i=1,n1m
      do k=1,n3m
      rhs(i,j,k)=q3(i,j,k)-vmepo(3,j)
      e3t(j)=e3t(j)+rhs(i,j,k)**2
      enddo
      enddo
      e3t(j)=e3t(j)*dvotz 
      enddo
c
c   azimuthal spectrum   q3 velocity
c
           call ruuthe(n2rm)
           call correl(n2rm)
              do jc=1,n2m
      do i=1,n1m
      r33the(jc,i)=r33the(jc,i)+ruuthd(jc,i)
      enddo
              enddo
      do jc=1,n2m
      e3stt(jc)=0.
      do k=1,k1max
      e3stt(jc)=e3stt(jc)+enet(jc,k)*dkk1
      enddo
      enddo
      do jc=1,n2m
      e3to(jc)=e3to(jc)+e3t(jc)
      e3stto(jc)=e3stto(jc)+e3stt(jc)
      do k=1,k1max
      ene3to(jc,k)=ene3to(jc,k)+enet(jc,k)
      enddo
      enddo
c
c   axial spectrum  q3 velocity
c
c
c   q3 velocity correlation
c
           call ruuaxi(n2rm)
              do jc=1,n2m
      do k=1,n3m
      r33axi(jc,k)=r33axi(jc,k)+ruuaxd(jc,k)
      enddo
              enddo
      do jc=1,n2m
      e3stz(jc)=0.
      do k=1,k3max
      e3stz(jc)=e3stz(jc)+enez(jc,k)*dkk3
      enddo
      enddo
      do jc=1,n2m
      e3stzo(jc)=e3stzo(jc)+e3stz(jc)
      do k=1,k3max
      ene3zo(jc,k)=ene3zo(jc,k)+enez(jc,k)
      enddo
      enddo
      do j=1,n2m
      do i=1,n1m
      do k=1,n3m
      rhs(i,j,k)=pr(i,j,k)-pmed(j)
      enddo
      enddo
      enddo
c
c   azimuthal spectrum  pressure
c
c
c   pressure correlation
c
           call ruuthe(n2rm)
           call correl(n2rm)
              do jc=1,n2m
      do i=1,n1m
      pcothe(jc,i)=pcothe(jc,i)+ruuthd(jc,i)
      enddo
              enddo
      do jc=1,n2m
      do k=1,k1max
      enepto(jc,k)=enepto(jc,k)+enet(jc,k)*dkk1
      enddo
      enddo
c
c   axial spectrum  pressure  
c
           call ruuaxi(n2rm)
              do jc=1,n2m
      do k=1,n3m
      pcoaxi(jc,k)=pcoaxi(jc,k)+ruuaxd(jc,k)
      enddo
              enddo
      do jc=1,n2m
      do k=1,k3max
      enepzo(jc,k)=enepzo(jc,k)+enez(jc,k)
      enddo
      enddo
c
c  ***********  compute the axial vorticity component
c               at         i,j,k+1/2
c
      do jc=2,n2m
      jm=jc-1
      do kc=1,n3m
      do ic=1,n1m
      im=imv(ic)
      dq1x2=(q1(ic,jc,kc)-q1(ic,jm,kc))*dx2/g2rc(jc)
      dq2x1=(q2(ic,jc,kc)-q2(im,jc,kc))*dx1/rc(jc)
      vorz=(dq1x2-dq2x1)/rc(jc)
      rhs(ic,jc,kc)=vorz-vompo(3,jc)
      enddo
      enddo
      enddo
      jc=n2
      jm=n2m
      do kc=1,n3m
      do ic=1,n1m
      dq1x2=-1./rc(jc)*q1(ic,jm,kc)*dx2*2./g2rc(jc)
      rhs(ic,jc,kc)=dq1x2-vompo(3,jc)
      enddo
      enddo
      jc=1
      jp=jc+1
      do kc=1,n3m
      vozcm=0.
      do ic=1,n1m
      im=imv(ic)
      dq1x2=(q1(ic,jp,kc)-q1(ic,jc,kc))*dx2/g2rc(jp)
      dq2x1=(q2(ic,jp,kc)-q2(im,jp,kc))*dx1/rc(jp)
      vorz=(dq1x2-dq2x1)/rc(jp)
      vozcm=vozcm+vorz
      enddo
      do ic=1,n1m
      rhs(ic,jc,kc)=vozcm/n1m-vompo(3,jc)
      enddo
      enddo
c
c   axial vorticity correlation
c
           call ruuthe(n2r)
           call correl(n2rm)
             do jc=1,n2m
      do i=1,n1m
      v33the(jc,i)=v33the(jc,i)+(ruuthd(jc,i)+ruuthd(jc+1,i))*0.5
      env33t(jc,i)=env33t(jc,i)+(enet(jc,i)+enet(jc+1,i))*0.5
      enddo
              enddo
           call ruuaxi(n2r)
             do jc=1,n2m
      do k=1,n3m
      v33axi(jc,k)=v33axi(jc,k)+(ruuaxd(jc,k)+ruuaxd(jc+1,k))*0.5
      env33z(jc,k)=env33z(jc,k)+(enez(jc,k)+enez(jc+1,k))*0.5
      enddo
              enddo
c
c  ***********  compute the radial  vorticity component
c                           at  i,j+1/2,k
c
      do jc=1,n2m
      do ic=1,n1m
      im=imv(ic)
      do kc=1,n3m
      km=kmv(kc)
      dq3x1=(q3(ic,jc,kc)-q3(im,jc,kc))*dx1/rm(jc)
      dq1x3=(q1(ic,jc,kc)-q1(ic,jc,km))*dx3/rm(jc)
      rhs(ic,jc,kc)=(dq3x1-dq1x3)-vompo(2,jc)
      enddo
      enddo
      enddo
c
c   radial vorticity correlation
c
           call ruuthe(n2rm)
           call correl(n2rm)
             do jc=1,n2m
      do i=1,n1m
      v22the(jc,i)=v22the(jc,i)+ruuthd(jc,i)
      env22t(jc,i)=env22t(jc,i)+enet(jc,i)
      enddo
              enddo
           call ruuaxi(n2rm)
             do jc=1,n2m
      do k=1,n3m
      v22axi(jc,k)=v22axi(jc,k)+ruuaxd(jc,k)
      env22z(jc,k)=env22z(jc,k)+enez(jc,k)
      enddo
              enddo
c
c  ***********  compute the azimuthal vorticity component
c               at         i+1/2,j,k
c
c  inside the field
      do jc=2,n2m
      jm=jmv(jc)
      do ic=1,n1m
      im=imv(ic)
      do kc=1,n3m
      km=kmv(kc)
      dq2x3=(q2(ic,jc,kc)-q2(ic,jc,km))*dx3/rc(jc)
      dq3x2=(q3(ic,jc,kc)-q3(ic,jm,kc))*dx2/g2rc(jc)
      rhs(ic,jc,kc)=(dq2x3-dq3x2)-vompo(1,jc)
      enddo
      enddo
      enddo
      jc=1
      jp=jc+1
      do ic=1,n1m
      do kc=1,n3m
      km=kmv(kc)
c  At  the axis
      q2akc= (q2(ic,jp,kc) - q2(isym(ic),jp,kc))*0.5/rc(jp)
      q2akm= (q2(ic,jp,km) - q2(isym(ic),jp,km))*0.5/rc(jp)
      q3axi= (q3(ic,jc,kc) + q3(isym(ic),jc,kc))*0.5
      dq2x3=(q2akc-q2akm)*dx3
      dq3x2=(q3(ic,jc,kc)-q3axi)*dx2*2./g2rc(jc)
      rhs(ic,jc,kc)=dq2x3-dq3x2-vompo(1,jc)
      enddo
      enddo
      jc=n2
      jm=jc-1
      do ic=1,n1m
      do kc=1,n3m
c  At the wall  (no-slip)
      dq3x2=(-q3(ic,jm,kc))*dx2*2./g2rc(jc)
      rhs(ic,jc,kc)=-dq3x2-vompo(1,jc)
      enddo
      enddo
c
c   azimuthal vorticity correlation
c
           call ruuthe(n2rm)
           call correl(n2rm)
             do jc=1,n2m
      do i=1,n1m
      v11the(jc,i)=v11the(jc,i)+(ruuthd(jc,i)+ruuthd(jc+1,i))*0.5
      env11t(jc,i)=env11t(jc,i)+(enet(jc,i)+enet(jc+1,i))*0.5
      enddo
              enddo
           call ruuaxi(n2r)
             do jc=1,n2m
      do k=1,n3m
      v11axi(jc,k)=v11axi(jc,k)+(ruuaxd(jc,k)+ruuaxd(jc+1,k))*0.5
      env11z(jc,k)=env11z(jc,k)+(enez(jc,k)+enez(jc+1,k))*0.5
      enddo
              enddo
          endif
      return                                                            
      end                                                               
c
c   ********************* subr inspru
c
      subroutine inspru(nba)
      include 'param.f' 
            nba=0
      do jc=1,n2m
      e3to(jc)=0.
      e2to(jc)=0.
      e1to(jc)=0.
      e3stto(jc)=0.
      e2stto(jc)=0.
      e1stto(jc)=0.
      e3stzo(jc)=0.
      e2stzo(jc)=0.
      e1stzo(jc)=0.
      do i=1,n1m
      v11the(jc,i)=0.
      v22the(jc,i)=0.
      v33the(jc,i)=0.
      enddo
      do i=1,n1m
      r11the(jc,i)=0.
      r22the(jc,i)=0.
      r33the(jc,i)=0.
      pcothe(jc,i)=0.
      scothe(jc,i)=0.
      enddo
      do k=1,n3m
      v11axi(jc,i)=0.
      v22axi(jc,i)=0.
      v33axi(jc,i)=0.
      enddo
      do k=1,n3m
      r11axi(jc,i)=0.
      r22axi(jc,i)=0.
      r33axi(jc,i)=0.
      pcoaxi(jc,i)=0.
      scoaxi(jc,i)=0.
      enddo
      do k=1,k1max
      ene1to(jc,k)=0.
      ene2to(jc,k)=0.
      ene3to(jc,k)=0.
      enepto(jc,k)=0.
      enddo
      enddo
      do jc=1,n2m
      do k=1,k3max
      ene1zo(jc,k)=0.
      ene2zo(jc,k)=0.
      ene3zo(jc,k)=0.
      enepzo(jc,k)=0.
      enddo
      enddo
      return
      end
c
