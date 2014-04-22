c************************************************************************
c                                                                       *
c                                                                       *
c     This code permits to perform velocity vorticity and pressure      *
c     visualizations in two dimensional planes .                        *
c     Even the full 3D flows can be visualized by in this case a       *
c     very large field with the 3 velocity and vorticity components     *
c     together with the pressure is generated.                          *
c                                                                       *
c     The quantities at the cell center are visualized                  *
c                                                                       *
c                                                                       *
c************************************************************************
      program main                                                      
      include 'param.f'
      common/timavg/timav
      common/timini/timei
      common/tstep/dt,beta,ren
      common/d2/nstop,nprint,ntst,npin,npstf
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/wrre/nwrit,nread
      common/strpar/str2
      common/tscoe/ga(3),ro(3),nsst
      common/d13/alx1,alx3
      common/vperin/vper
      common/averou/iav
      common/islwal/islv1s,islv1n,islv3s,islv3n
      common/jrrr/jri,jrf,djr,irejr,iruuca
      common/newdat/icost,timeav
      common/movwal/tosc,uosc
      common/slotfl/flowq2,tau2
      common/slotin/tim0sl
      common/slotpa/y1gsl,y1ssl,y3gsl,y3ssl
      common/slotdi/y1disl,y3disl
      common/oldso/ifield
      common/cflco/icfl,cflc,tpin,tprin,tfin
      open(15,file='chapn.d')
      read(15,*) n1,n2,n3,nsst
      read(15,*) nwrit,nread,iav,iprfi
      read(15,*) alx3d,alx1d,str2
      read(15,*) ren,vper
      read(15,*) dt,ntst,nprint,npin,npstf
      read(15,*) nstop
      read(15,*) icfl,cflc,tpin,tprin,tfin
      read(15,*) jri,jrf,djr,irejr,iruuca
      read(15,*) timeav
      read(15,*) islv1s,islv1n,islv3s,islv3n
      read(15,*) tosc,uosc
      read(15,*) flowq2,tau2,tim0sl
      read(15,*) y1gsd,y1ssd,y3gsd,y3ssd
      read(15,*) y1disd,y3disd
      read(15,*) ifield
      pi=2.*asin(1.)
      n1m=n1-1
      n2m=n2-1
      n3m=n3-1
      alx1=alx1d*pi
      alx3=alx3d*pi
      y1gsl=y1gsd*pi
      y1ssl=y1ssd*pi
      y3gsl=y3gsd*pi
      y3ssl=y3ssd*pi
      y1disl=y1disd*pi
      y3disl=y3disd*pi
      write(6,*)'*****************************************************'
      write(32,*)'*****************************************************'
      write(6,*)'*                                                    *'
      write(32,*) '*                                                  *'
      write(6,*)'*               CHANNEL  WITH TRANSPIRATION          *'
      write(32,*)'*              CHANNEL  WITH TRANSPIRATION          *'
      write(6,*)'*****************************************************'
      write(32,*)'*****************************************************'
      write(6,*)'  '
      write(32,*) '  '
      write(6,*)'  '
      write(32,*) '  '
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
c   This code uses routines just used for the code in CODESCH9/CHA
c   and for the code for the postprocessing in CODESCH9/POST2
c                                                                       *
c************************************************************************
c************************************************************************
      subroutine gcurv                                                  
      include 'param.f'
      common/timini/timei
      common/tstep/dt,beta,ren
      dimension y(m2)
      dimension ru(ndv,m1,m2,m3),sij(2*ndv,m1,m2,m3)
      dimension q(ndv,m1,m2,m3),dq(m1,m2,m3)
     1         ,pr(m1,m2,m3)
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/d13/alx1,alx3
      common/i3dpr/i3dou,iprq(15),jprq(15),kprq(15)
      write(6,*) 'enter  timei'
      read(5,*) timei
c                                                                       
c     grid definition, indices and mesh size calculation                          
c                                                                       
      call meshes
      call indic
      call coordi(y)
      call metric(y)
      read(15,301) dummy
301   format(a4)                                                        
c
c   select the location of the planes where to perform the
c   visualizations
c
c
c   i3do =1 the full field is visualized
c   the output file should be very large
c
       read(15,*) i3do
      read(15,301) dummy
       read(15,*) npj,(jprq(n),n=1,npj)
      npi=3
c
c   here only 3  x2-x3 planes are selected 
c
      iprq(3)=n1m/4+1
      iprq(2)=n1m/2+1
      iprq(1)=3*n1m/4+1
      npk=3
c
c   here only 3  x2-x1 planes are selected 
c
      kprq(3)=n3m/4+1
      kprq(2)=n3m/2+1
      kprq(1)=3*n3m/4+1
       do n=1,npi
       write(6,*)'n,iprq(n) ',n,iprq(n)
       enddo
       do n=1,npj
       write(6,*)'n,jprq(n) ',n,jprq(n)
       enddo
       do n=1,npk
       write(6,*)'n,kprq(n) ',n,kprq(n)
       enddo
c
c     print the grid
c
      call pricor
c
c
c     print some informations on the run
c
      write(6,754)n1,n2,n3
  754 format(10x,'centered velocities',2x
     1      ,5x,'n1=',i3,2x,'n2=',i3,2x,'n3=',i3)
      write(6,755) ren,dx1,dx2,dx3,ntst
  755 format(3x,'ren=',e10.3,3x,'dx1=',e10.3,3x,'dx2=',e10.3,3x
     1      ,3x,'dx3=',e10.3,3x,
     1      'ntst=',i5,3x)
  159 format(1x,i4,2x,e10.4,3e10.3,3(1x,e10.4,1x,i3,1x,i3),e10.3)
c
c                                                                       
      ncount=0
      ntii=1
       time=timei
      ntime=ntii
c
c   the field is read
c                                                 
      call tschem(q,pr,ru,dq,time,sij)
c
c   the output fields are written
c
      call outpf(time,q,ru,pr,sij)
c                                                                       
c                                                                       
      return                                                            
      end                                                               
c                                                                       
c
c  **************  subrout tschem ***************************************
c
c
      subroutine tschem(q,pr,ru,dq,time,sij)
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      dimension q(ndv,m1,m2,m3),dq(m1,m2,m3),pr(m1,m2,m3)
      dimension ru(ndv,m1,m2,m3),sij(2*ndv,m1,m2,m3)
          call inirea(ntii,time,q,ru,pr)
      call velc(q,pr)
      call vorc(q,ru,dq,sij,pr)
      return                                                            
      end                                                               
c                                                                       
c  ************************* subrout,velc  **********************  
c                                                                       
c     this subroutine calculates the velocity at the centre of the cell
c     and then the fluctuations
c                                                                       
      subroutine velc(q,pr)
      include 'param.f'
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      dimension pr(m1,m2,m3)
      dimension q(ndv,m1,m2,m3)
      common/corpri/yp1(m1),yp2(m2),yp3(m3)
      dimension volz(m2)
      common/d13/alx1,alx3
      common/vmean/vmed(3,m2),pmed(m2)
      common/vmeao/vmeo(3,m2),pmeo(m2)
      common/vmpos/vmepo(3,m2),vompo(3,m2)
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/veltot/vit(4)
      avgn=1./(float(n1m*n3m))
      vl123=1./(n1m*n2m*n3m)
      volto=1./(2*alx3*alx1)
c
c  velocities at the center of the cell
c
      vit(1)=0.
      vit(2)=0.
      vit(3)=0.
      vit(4)=0.
      do jc=1,n2m                                                     
      jm=jmv(jc)                                                        
      jp=jc+1
      pmed(jc)=0.
      do l=1,3
      vmed(l,jc)=0.
      vmepo(l,jc)=0.
      enddo
      do kc=1,n3m 
      do ic=1,n1m                                                     
      dphc=(q(2,ic,jc,kc)+q(2,ic,jp,kc))*0.5
      dqc=( q(1,ic,jc,kc)+q(1,ipv(ic),jc,kc) )*0.5
      qcapc=(q(3,ic,jc,kc)+q(3,ic,jc,kpv(kc)) )*0.5 
      vmed(2,jc)=vmed(2,jc)+dphc
      vmed(1,jc)=vmed(1,jc)+dqc
      vmed(3,jc)=vmed(3,jc)+qcapc
      vmepo(1,jc)=vmepo(1,jc)+q(1,ic,jc,kc)
      vmepo(3,jc)=vmepo(3,jc)+q(3,ic,jc,kc)
      volz(jc)=(yp1(ic+1)-yp1(ic))*(yp2(jc+1)-yp2(jc))*yp3(kc+1)-yp3(kc))
      vit(1)=dqc*volz(jc)+vit(1)
      vit(2)=dphc*volz(jc)+vit(2)
      vit(3)=qcapc*volz(jc)+vit(3)
      vit(4)=pr(ic,jc,kc)*volz(jc)+vit(4)
      pmed(jc)=pmed(jc)+pr(ic,jc,kc)
      enddo
      enddo
      do l=1,3
      vmed(l,jc)=vmed(l,jc)*avgn
      vmepo(l,jc)=vmepo(l,jc)*avgn
      enddo
      pmed(jc)=pmed(jc)*avgn
      enddo
      do jc=1,n2                                                     
      vmepo(2,jc)=0.
      do kc=1,n3m 
      do ic=1,n1m                                                     
      vmepo(2,jc)=vmepo(2,jc)+q(2,ic,jc,kc)
      enddo
      enddo
      vmepo(2,jc)=vmepo(2,jc)*avgn
      enddo
      do l=1,4
      vit(l)=vit(l)*volto
      enddo
      return
      end
c
c                                                                       
c  ************************ subrout vorc  **********************  
c                                                                       
c     this subroutine calculates the vorticity components
c     and the six components  of strain rate  tensor
c                                                                       
      subroutine vorc(q,voc,rhs,sij,pr)
      include 'param.f'
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/corpri/yp1(m1),yp2(m2),yp3(m3)
      common/y2sta/y2s(m2)
      common/y13sta/y1s(m1),y3s(m3)
      dimension volz(m2)
      common/d13/alx1,alx3
      common/qurms/verms(3,m2),vorms(3,m2),prrms(m2)
      common/vmean/vmed(3,m2),pmed(m2)
      common/vmeao/vmeo(3,m2),pmeo(m2)
      common/vmpos/vmepo(3,m2),vompo(3,m2)
      common/vomean/vorv(3,m2)
      common/dim/n1,n1m,n2,n2m,n3,n3m
      dimension q(ndv,m1,m2,m3),rhs(m1,m2,m3),qapp(m1,m2,m3)
      dimension voc(ndv,m1,m2,m3)
      dimension pr(m1,m2,m3)
      dimension sij(2*ndv,m1,m2,m3)
      dimension nj1(2),nj2(2)
      common/stmean/strac(6,m2),strai(6,m2)
      avgn=1./(float(n1m*n3m))
c                                                                       
c  ***********  compute the spanwise vorticity component               
c               at         i+1/2,j,k 

c                                                                       
c  inside the field
      do jc=2,n2m
      jm=jmv(jc)
      vompo(1,jc)=0.
      strai(4,jc)=0.
            do ic=1,n1m
      im=imv(ic)
                  do kc=1,n3m
      km=kmv(kc)
      udx2l=1./(y2s(jc)-y2s(jm))
      udx1l=1./(y1s(ic)-y1s(im))
      udx3l=1./(y3s(kc)-y3s(km))
      dq3x2=(q(3,ic,jc,kc)-q(3,ic,jm,kc))*udx2l
      dq2x3=(q(2,ic,jc,kc)-q(2,ic,jc,km))*udx3l
      rhs(ic,jc,kc)=dq3x2-dq2x3
      vompo(1,jc)=vompo(1,jc)+rhs(ic,jc,kc)
      qapp(ic,jc,kc)=(dq3x2+dq2x3)*0.5
      strai(4,jc)=strai(4,jc)+qapp(ic,jc,kc)
                  enddo
            enddo
      vompo(1,jc)=vompo(1,jc)*avgn
      strai(4,jc)=strai(4,jc)*avgn
      enddo
c
c  At the walls  (no-slip)
c
      jc=1
      vompo(1,jc)=0.
      strai(4,jc)=0.
      jp=jc+1
            do ic=1,n1m
                  do kc=1,n3m
      km=kmv(kc)
      udx2l=1./(y2s(jc)-yp2(jc))
      dq3x2=+q(3,ic,jc,kc)*udx2l
      rhs(ic,jc,kc)=dq3x2
      vompo(1,jc)=vompo(1,jc)+rhs(ic,jc,kc)
      qapp(ic,jc,kc)=dq3x2*0.5
      strai(4,jc)=strai(4,jc)+qapp(ic,jc,kc)
                  enddo
            enddo
      vompo(1,jc)=vompo(1,jc)*avgn
      strai(4,jc)=strai(4,jc)*avgn
      jc=n2
      vompo(1,jc)=0.
      strai(4,jc)=0.
      jm=jc-1
            do ic=1,n1m
                  do kc=1,n3m
      km=kmv(kc)
      udx2l=1./(yp2(jc)-y2s(jm))
      dq3x2=-q(3,ic,jm,kc)*udx2l
      rhs(ic,jc,kc)=dq3x2
      vompo(1,jc)=vompo(1,jc)+rhs(ic,jc,kc)
      qapp(ic,jc,kc)=dq3x2*0.5
      strai(4,jc)=strai(4,jc)+qapp(ic,jc,kc)
                  enddo
            enddo
      vompo(1,jc)=vompo(1,jc)*avgn
      strai(4,jc)=strai(4,jc)*avgn
c
c   spanwise vorticity at the cell centre 
c
      do jc=1,n2m
      vorv(1,jc)=0.
      strac(4,jc)=0.
            do kc=1,n3m     
      kp=kpv(kc)                                              
                  do ic=1,n1m
      dqc=(rhs(ic,jc,kc)+rhs(ic,jc+1,kc)+
     1     rhs(ic,jc,kp)+rhs(ic,jc+1,kp))*0.25
      vorv(1,jc)=vorv(1,jc)+dqc
      voc(1,ic,jc,kc)=dqc
      dsc=(qapp(ic,jc,kc)+qapp(ic,jc+1,kc)+
     1     qapp(ic,jc,kp)+qapp(ic,jc+1,kp))*0.25
      sij(4,ic,jc,kc)=dsc
      strac(4,jc)=strac(4,jc)+dsc
                  enddo
            enddo
      vorv(1,jc)=vorv(1,jc)*avgn
      strac(4,jc)=strac(4,jc)*avgn
      enddo
c                                                                  
c  ***********  compute the normal  vorticity component            
c                           at  i,j+1/2,k
c                                                                  
      do jc=1,n2m
      vompo(2,jc)=0.
      strai(5,jc)=0.
            do ic=1,n1m
      im=imv(ic)
                  do kc=1,n3m
      km=kmv(kc)
      udx1l=1./(y1s(ic)-y1s(im))
      udx3l=1./(y3s(kc)-y3s(km))
      dq1x3=(q(1,ic,jc,kc)-q(1,ic,jc,km))*udx3l
      dq3x1=(q(3,ic,jc,kc)-q(3,im,jc,kc))*udx1l
      rhs(ic,jc,kc)=dq1x3-dq3x1
      vompo(2,jc)=vompo(2,jc)+rhs(ic,jc,kc)
      qapp(ic,jc,kc)=(dq1x3+dq3x1)*0.5
      strai(5,jc)=strai(5,jc)+qapp(ic,jc,kc)
                  enddo
            enddo
      vompo(2,jc)=vompo(2,jc)*avgn
      strai(5,jc)=strai(5,jc)*avgn
      enddo

c
c   vorticity at the cell centre inner field
c
      do jc=1,n2m
      vorv(2,jc)=0.
      strac(5,jc)=0.
                  do kc=1,n3m                                                   
      kp=kpv(kc)
                  do ic=1,n1m
      ip=ipv(ic)                                                  
      dphc=(rhs(ic,jc,kc)+rhs(ip,jc,kc)+
     1      rhs(ic,jc,kp)+rhs(ip,jc,kp))*0.25
      vorv(2,jc)=vorv(2,jc)+dphc
      voc(2,ic,jc,kc)=dphc
      dsc=(qapp(ic,jc,kc)+qapp(ip,jc,kc)+
     1      qapp(ic,jc,kp)+qapp(ip,jc,kp))*0.25
      sij(5,ic,jc,kc)=dsc
      strac(5,jc)=strac(5,jc)+dsc
                  enddo
                  enddo
      vorv(2,jc)=vorv(2,jc)*avgn
      strac(5,jc)=strac(5,jc)*avgn
      enddo
c                                                                  
c  ***********  compute the streamwise vorticity component        
c               at         i,j,k+1/2        
c
      do jc=2,n2m                                               
      vompo(3,jc)=0.
      strai(6,jc)=0.
      jm=jc-1                                                  
            do kc=1,n3m                                                
                  do ic=1,n1m                                               
      im=imv(ic)                                                  
      udx2l=1./(y2s(jc)-y2s(jm))
      udx1l=1./(y1s(ic)-y1s(im))
      dq1x2=(q(1,ic,jc,kc)-q(1,ic,jm,kc))*udx2l
      dq2x1=(q(2,ic,jc,kc)-q(2,im,jc,kc))*udx1l
      rhs(ic,jc,kc)=dq2x1-dq1x2
      vompo(3,jc)=vompo(3,jc)+rhs(ic,jc,kc)
      qapp(ic,jc,kc)=(dq1x2+dq2x1)*0.5
      strai(6,jc)=strai(6,jc)+qapp(ic,jc,kc)
                  enddo
            enddo
      vompo(3,jc)=vompo(3,jc)*avgn
      strai(6,jc)=strai(6,jc)*avgn
      enddo
      jc=n2
      vompo(3,jc)=0.
      strai(6,jc)=0.
      jm=n2m
            do kc=1,n3m                                                
                  do ic=1,n1m                                               
      im=imv(ic)                                                  
      udx2l=1./(yp2(jc)-y2s(jm))
      udx1l=1./(y1s(ic)-y1s(im))
      dq1x2=-q(1,ic,jm,kc)*udx2l
      rhs(ic,jc,kc)=-dq1x2
      vompo(3,jc)=vompo(3,jc)+rhs(ic,jc,kc)
      qapp(ic,jc,kc)=(dq1x2)*0.5
      strai(6,jc)=strai(6,jc)+qapp(ic,jc,kc)
                  enddo
            enddo
      vompo(3,jc)=vompo(3,jc)*avgn
      strai(6,jc)=strai(6,jc)*avgn
      jc=1
      vompo(3,jc)=0.
      strai(6,jc)=0.
            do kc=1,n3m
                  do ic=1,n1m
      im=imv(ic)                                                  
      udx2l=1./(y2s(jc)-yp2(jc))
      udx1l=1./(y1s(ic)-y1s(im))
      dq1x2=q(1,ic,jc,kc)*udx2l
      rhs(ic,jc,kc)=-dq1x2
      vompo(3,jc)=vompo(3,jc)+rhs(ic,jc,kc)
      qapp(ic,jc,kc)=(dq1x2)*0.5
      strai(6,jc)=strai(6,jc)+qapp(ic,jc,kc)
                  enddo
            enddo
      vompo(3,jc)=vompo(3,jc)*avgn
      strai(6,jc)=strai(6,jc)*avgn

c
c   vorticity at the cell centre 
c
      do jc=1,n2m
      vorv(3,jc)=0.
      strac(6,jc)=0.
            do kc=1,n3m                                                   
                  do ic=1,n1m
      ip=ipv(ic)                                                  
      qcapc=(rhs(ic,jc,kc)+rhs(ic,jc+1,kc)+
     1       rhs(ip,jc,kc)+rhs(ip,jc+1,kc))*0.25
      voc(3,ic,jc,kc)=qcapc
      vorv(3,jc)=vorv(3,jc)+qcapc
      dsc=(qapp(ic,jc,kc)+qapp(ic,jc+1,kc)+
     1       qapp(ip,jc,kc)+qapp(ip,jc+1,kc))*0.25
      sij(6,ic,jc,kc)=dsc
      strac(6,jc)=strac(6,jc)+dsc
                  enddo
            enddo
      vorv(3,jc)=vorv(3,jc)*avgn
      strac(6,jc)=strac(6,jc)*avgn
      enddo
c
c   normal strain at the cell centre 
c
      do jc=1,n2m
      jp=jc+1
      strac(1,jc)=0.
      strac(2,jc)=0.
      strac(3,jc)=0.
            do kc=1,n3m                                                   
      kp=kpv(kc)                                                  
                  do ic=1,n1m
      ip=ipv(ic)                                                  
      udx1l=1./(yp1(ic+1)-yp1(ic))
      udx2l=1./(yp2(jc+1)-yp2(jc))
      udx3l=1./(yp3(kc+1)-yp3(kc))
      dq1x1=(q(1,ip,jc,kc)-q(1,ic,jc,kc))*udx1l
      dq2x2=(q(2,ic,jp,kc)-q(2,ic,jc,kc))*udx2l
      dq3x3=(q(3,ic,jc,kp)-q(3,ic,jc,kc))*udx3l
      sij(1,ic,jc,kc)=dq1x1
      sij(2,ic,jc,kc)=dq2x2
      sij(3,ic,jc,kc)=dq3x3
      strac(1,jc)=strac(1,jc)+dq1x1
      strac(2,jc)=strac(2,jc)+dq2x2
      strac(3,jc)=strac(3,jc)+dq3x3
                  enddo
            enddo
      strac(1,jc)=strac(1,jc)*avgn
      strac(2,jc)=strac(2,jc)*avgn
      strac(3,jc)=strac(3,jc)*avgn
      enddo
c
c   evaluation of rms quantities
c   vel vor and press
c
      do j=1,n2m
      jp=j+1
      prrms(j)=0.
      do n=1,3
      verms(n,j)=0.
      vorms(n,j)=0.
      enddo
      avn1n3=1./float(n1m*n3m)
      do k=1,n3m                                                   
      do i=1,n1m
      q2c=
     1  (    q(2,i,j,k)+q(2,i,jp,k))*0.5
      q1c=
     1  (    q(1,i,j,k)+q(1,ipv(i),j,k) )*0.5
      q3c=(q(3,i,j,k)+q(3,i,j,kpv(k)) )*0.5 
      prp=pr(i,j,k)-pmed(j)
      q1p=q1c-vmed(1,j)
      q2p=q2c-vmed(2,j)
      q3p=q3c-vmed(3,j)
      vo1p=voc(1,i,j,k)-vorv(1,j)
      vo2p=voc(2,i,j,k)-vorv(2,j)
      vo3p=voc(3,i,j,k)-vorv(3,j)
      verms(1,j)=verms(1,j)+q1p**2
      verms(2,j)=verms(2,j)+q2p**2
      verms(3,j)=verms(3,j)+q3p**2
      vorms(1,j)=vorms(1,j)+vo1p**2
      vorms(2,j)=vorms(2,j)+vo2p**2
      vorms(3,j)=vorms(3,j)+vo3p**2
      prrms(j)=prrms(j)+prp**2
      enddo
      enddo
      verms(1,j)=sqrt(verms(1,j)*avn1n3)
      verms(2,j)=sqrt(verms(2,j)*avn1n3)
      verms(3,j)=sqrt(verms(3,j)*avn1n3)
      vorms(1,j)=sqrt(vorms(1,j)*avn1n3)
      vorms(2,j)=sqrt(vorms(2,j)*avn1n3)
      vorms(3,j)=sqrt(vorms(3,j)*avn1n3)
      prrms(j)=sqrt(prrms(j)*avn1n3)
      enddo
      do nsub=1,3
      if(nsub.eq.1) then
      nj1(1)=1
      nj2(1)=20
      nj1(2)=n2m-20
      nj2(2)=n2m
                    endif
      if(nsub.eq.2) then
      nj1(1)=21 
      nj2(1)=53
      nj1(2)=n2m-53
      nj2(2)=n2m-21
                    endif
      if(nsub.eq.3) then
      nj1(1)=54 
      nj2(1)=75
      nj1(2)=n2m-76
      nj2(2)=n2m-54
                    endif
      enema=-1000.
      ensma=-1000.
      disma=-1000.
      espma=-1000.
      plama=-1000.
      helma=-1000.
      detma=-1000.
      enemi=+1000.
      ensmi=+1000.
      dismi=+1000.
      espmi=+1000.
      plami=+1000.
      helmi=+1000.
      detmi=+1000.
      do k=1,n3m                                                   
      do i=1,n1m
      do npar=1,2
      do j=nj1(npar),nj2(npar)
      jp=j+1
      q2c=
     1  (    q(2,i,j,k)+q(2,i,jp,k))*0.5
      q1c=
     1  (    q(1,i,j,k)+q(1,ipv(i),j,k) )*0.5
      q3c=(q(3,i,j,k)+q(3,i,j,kpv(k)) )*0.5 
      q1p=q1c-vmed(1,j)
      q2p=q2c-vmed(2,j)
      q3p=q3c-vmed(3,j)
      vo1p=voc(1,i,j,k)-vorv(1,j)
      vo2p=voc(2,i,j,k)-vorv(2,j)
      vo3p=voc(3,i,j,k)-vorv(3,j)
      s11p=sij(1,i,j,k)-strac(1,j)
      s22p=sij(2,i,j,k)-strac(2,j)
      s33p=sij(3,i,j,k)-strac(3,j)
      s23p=sij(4,i,j,k)-strac(4,j)
      s13p=sij(5,i,j,k)-strac(5,j)
      s12p=sij(6,i,j,k)-strac(6,j)
      enep=q1p**2+q2p**2+q3p**2
      ensp=vo1p**2+vo2p**2+vo3p**2
      helep=(q1p*vo1p+q2p*vo2p+q3p*vo3p)/sqrt(enep*ensp)
      dislo=(s11p**2+s22p**2+s33p**2+
     1      2.*(s12p**2+s13p**2+s23p**2))
      enspro=vo1p**2*s11p+vo2p**2*s22p+vo3p**2*s33p+
     1      2.*(vo1p*vo2p*s12p+vo2p*vo3p*s23p+vo1p*vo3p*s13p)
      detsij=s11p*(s22p*s33p-s23p*s23p)
     1      -s12p*(s12p*s33p-s23p*s13p)
     1      +s13p*(s12p*s23p-s22p*s13p)
      prlap=ensp/4.-dislo
      enema=max(enema,enep)
      enemi=min(enemi,enep)
      ensma=max(ensma,ensp)
      ensmi=min(ensmi,ensp)
      disma=max(disma,dislo)
      dismi=min(dismi,dislo)
      helma=max(helma,helep)
      helmi=min(helmi,helep)
      detma=max(detma,detsij)
      detmi=min(detmi,detsij)
      espma=max(espma,enspro)
      espmi=min(espmi,enspro)
      plama=max(plama,prlap)
      plami=min(plami,prlap)
      enddo
      enddo
      enddo
      enddo
      write(6,*)nsub,nj1(1),nj2(1),nj1(2),nj2(2)
      write(18,*)nsub,nj1(1),nj2(1),nj1(2),nj2(2)
      write(6,*)'max qua in vorc'
     1    ,enema,ensma,disma,espma,plama,helma,detma
      write(6,*)'min qua in vorc'
     1    ,enemi,ensmi,dismi,espmi,plami,helmi,detmi
      write(18,*)'max qua in vorc'
     1    ,enema,ensma,disma,espma,plama,helma,detma
      write(18,*)'min qua in vorc'
     1    ,enemi,ensmi,dismi,espmi,plami,helmi,detmi
      enddo
      return                                                            
      end                                                               
************************************************************************
c                                                                       *
c     ********* subrout pricor ******************                       *
c                                                                       *
c************************************************************************
      subroutine pricor
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      common/y2sta/y2s(m2)
      common/y13sta/y1s(m1),y3s(m3)
      common/i3dpr/i3dou,iprq(15),jprq(15),kprq(15)
      character*60 namfile
      namfile='cord13.dat'
      open(18,file=namfile,form='unformatted')
      aaa=1.
      write(6,*) 'in prico 13' ,n3m,n1m
      write(18) n3m,n1m,1
      write(18)
     1            ((y3s(k),k=1,n3m),i=1,n1m),
     1            ((y1s(i),k=1,n3m),i=1,n1m),
     1            ((aaa,k=1,n3m),j=1,n1m)
      close(18)
      namfile='cord12.dat'
      open(18,file=namfile,form='unformatted')
      write(6,*) 'in prico 12' ,n1,n2m
      write(18) n1m,n2m,1
      write(18)
     1   ((y1s(i),i=1,n1m),j=1,n2m),
     1   ((y2s(j),i=1,n1m),j=1,n2m),
     1   ((aaa,i=1,n1m),j=1,n2m)
      close(18)
      namfile='cord32.dat'
      open(18,file=namfile,form='unformatted')
      write(18) n3m,n2m,1
      write(18)
     1            ((y3s(k),k=1,n3m),j=1,n2m),
     1            ((y2s(j),k=1,n3m),j=1,n2m),
     1            ((aaa,k=1,n3m),j=1,n2m)
      close(18)
      if(i3dou.eq.1) then
      namfile='cordto.dat'
      open(18,file=namfile,form='unformatted')
      write(18) n1m,n2m,n3m
      write(18)
     1   (((y1s(i),i=1,n1m),j=1,n2m),k=1,n3m),
     1   (((y2s(j),i=1,n1m),j=1,n2m),k=1,n3m),
     1   (((y3s(k),i=1,n1m),j=1,n2m),k=1,n3m)
      close(18)
                        endif
      return
      end
c************************************************************************
c                                                                       *
c  ****************************** subrout vmaxv **********************  *
c                                                                       *
c************************************************************************
      subroutine vmaxv(q,n1mv,n2mv,n3mv)
      include 'param.f'
      common/velmax/vmax(ndv),vmaxo(ndv)
      dimension q(ndv,m1,m2,m3)
c
c     find the maximum velocities in order to check convergence or
c     to derive stability conditions.
c
      vca=0.
c
      do l=1,ndv
      vmax(l)=0.
      do k=1,n3mv
      do j=1,n2mv
      do i=1,n1mv
      if(l.eq.1) vca=q(1,i,j,k)
      if(l.eq.2) vca=q(2,i,j,k)
      if(l.eq.3) vca=q(3,i,j,k)
      if(abs(vca).ge.vmax(l)) then
      vmax(l)=abs(vca)
                              endif
      enddo
      enddo
      enddo
      enddo
      return
      end
c************************************************************************
c                                                                       *
c  ****************************** subrout outpf  ********************** *
c                                                                       *
c************************************************************************
      subroutine outpf(time,q,voc,pr,sij)
      include 'param.f'
      dimension q1rz(m2,m3),q2rz(m2,m3),q3rz(m2,m3)
      dimension vo1rz(m2,m3),vo2rz(m2,m3),vo3rz(m2,m3)
      dimension prrz(m2,m3)
      dimension q1tr(m1,m2),q2tr(m1,m2),q3tr(m1,m2)
      dimension vo1tr(m1,m2),vo2tr(m1,m2),vo3tr(m1,m2)
      dimension prtr(m1,m2)
      dimension q1tz(m1,m3),q2tz(m1,m3),q3tz(m1,m3)
      dimension vo1tz(m1,m3),vo2tz(m1,m3),vo3tz(m1,m3)
      dimension prtz(m1,m3),dpr1tz(m1,m3),dpr3tz(m1,m3)
      dimension q(ndv,m1,m2,m3),pr(m1,m2,m3)
      dimension qua(7,m3,m3),quarm(7,m3,m3)
      dimension voc(ndv,m1,m2,m3)
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      common/y2sta/y2s(m2)
      common/y13sta/y1s(m1),y3s(m3)
      common/qurms/verms(3,m2),vorms(3,m2),prrms(m2)
      common/vmean/vmed(3,m2),pmed(m2)
      common/vmeao/vmeo(3,m2),pmeo(m2)
      common/vmpos/vmepo(3,m2),vompo(3,m2)
      common/vomean/vorv(3,m2)
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/i3dpr/i3dou,iprq(15),jprq(15),kprq(15)
      dimension sij(2*ndv,m1,m2,m3)
      common/stmean/strac(6,m2),strai(6,m2)

      character*60 namfi3
      character*3 nkpse,njpse
      character*4 navps
      npq=1
c
c     form the name of the file
c
      nptz=0
      irep=re
      itime=nint(time)
   98 format(i4.4)
   82 format(i4.4)
   83 format(i3.3)
c
c   write the profiles of the mean and the rms
c   these are profiles obtained by only one field
c   these should be compared with those averaged also
c   in time
c
      open(49,file='profirms.out')
      open(29,file='profiles.out')
      open(39,file='profisij.out')
      do j=1,n2m
      write(49,137)y2s(j),(verms(l,j),l=1,3),(vorms(l,j),l=1,3),
     1                     prrms(j)
      write(29,137)y2s(j),(vmed(l,j),l=1,3),(vorv(l,j),l=1,3),
     1                     pmed(j)
      write(39,137)y2s(j),(strac(l,j),l=1,6)
      enddo
  137 format(1x,e13.5,2x,7e14.7)
      close(29)
      close(39)
      close(49)
c
c  cross-stream  section x1-x2
c
      npk=1
      do k=n3m,1,-1
      if(k.eq.kprq(npk)) then
      write(nkpse,83)k
      do i=1,n1m
      do j=1,n2m
      jp=j+1
      q2c=
     1  (    q(2,i,j,k)+q(2,i,jp,k))*0.5
      q1c=
     1  (    q(1,i,j,k)+q(1,ipv(i),j,k) )*0.5
      q3c=(q(3,i,j,k)+q(3,i,j,kpv(k)) )*0.5 
      q1p=q1c-vmed(1,j)
      q2p=q2c-vmed(2,j)
      q3p=q3c-vmed(3,j)
      vo1p=voc(1,i,j,k)-vorv(1,j)
      vo2p=voc(2,i,j,k)-vorv(2,j)
      vo3p=voc(3,i,j,k)-vorv(3,j)
      q1tr(i,j)=q1p
      q2tr(i,j)=q2p
      q3tr(i,j)=q3p
      prtr(i,j)=pr(i,j,k)-pmed(j)
      vo1tr(i,j)=vo1p
      vo2tr(i,j)=vo2p
      vo3tr(i,j)=vo3p
      s11p=sij(1,i,j,k)-strac(1,j)
      s22p=sij(2,i,j,k)-strac(2,j)
      s33p=sij(3,i,j,k)-strac(3,j)
      s23p=sij(4,i,j,k)-strac(4,j)
      s13p=sij(5,i,j,k)-strac(5,j)
      s12p=sij(6,i,j,k)-strac(6,j)
      enep=q1p**2+q2p**2+q3p**2
      ensp=vo1p**2+vo2p**2+vo3p**2
      helep=(q1p*vo1p+q2p*vo2p+q3p*vo3p)/sqrt(enep*ensp)
      dislo=(s11p**2+s22p**2+s33p**2+
     1      2.*(s12p**2+s13p**2+s23p**2))
      enspro=vo1p**2*s11p+vo2p**2*s22p+vo3p**2*s33p+
     1      2.*(vo1p*vo2p*s12p+vo2p*vo3p*s23p+vo1p*vo3p*s13p)
      prlap=ensp/4.-dislo
      qua(1,i,j)=enep
      qua(2,i,j)=ensp/4.
      qua(3,i,j)=helep
      qua(4,i,j)=dislo
      qua(5,i,j)=enspro
      qua(6,i,j)=prlap
      detsij=s11p*(s22p*s33p-s23p*s23p)
     1      -s12p*(s12p*s33p-s23p*s13p)
     1      +s13p*(s12p*s23p-s22p*s13p)
      qua(7,i,j)=detsij
      enddo
      enddo
c
c   here velocity vorticity components and pressure are written
c   in a file compatible with TURB3D package
c
      namfi3='verseupo'//nkpse//'.dat'
      open(59,file=namfi3,form='unformatted')
      rewind 59
      write(59) n1m,n2m,1
      write(59) re,re,re,time
      write(59)
     1   ((q1tr(i,j),i=1,n1m),j=1,n2m),
     1   ((q2tr(i,j),i=1,n1m),j=1,n2m),
     1   ((q3tr(i,j),i=1,n1m),j=1,n2m),
     1   ((prtr(i,j),i=1,n1m),j=1,n2m),
     1   ((vo1tr(i,j),i=1,n1m),j=1,n2m),
     1   ((vo2tr(i,j),i=1,n1m),j=1,n2m),
     1   ((vo3tr(i,j),i=1,n1m),j=1,n2m)
      close(59)
c
c  a bunch of other statistics are visualized
c
      namfi3='verseeeh'//nkpse//'.dat'
      open(59,file=namfi3,form='unformatted')
      rewind 59
      write(59) n1m,n2m,1
      write(59) re,re,re,time
      write(59)
     1   ((qua(1,i,j),i=1,n1m),j=1,n2m),
     1   ((qua(2,i,j),i=1,n1m),j=1,n2m),
     1   ((qua(3,i,j),i=1,n1m),j=1,n2m),
     1   ((qua(4,i,j),i=1,n1m),j=1,n2m),
     1   ((qua(5,i,j),i=1,n1m),j=1,n2m),
     1   ((qua(6,i,j),i=1,n1m),j=1,n2m)
     1  ,((qua(7,i,j),i=1,n1m),j=1,n2m)
      close(59)
      npk=npk+1
                            endif
      enddo
c
c  horizonthal section x3-x1 planes
c
      npj=1
      do j=n2m,1,-1
      if(j.eq.jprq(npj)) then
      jp=j+1
      write(njpse,83)j
      do k=1,n3m
      do i=1,n1m
      q2c=
     1  (    q(2,i,j,k)+q(2,i,jp,k))*0.5
      q1c=
     1  (    q(1,i,j,k)+q(1,ipv(i),j,k) )*0.5
      q3c=(q(3,i,j,k)+q(3,i,j,kpv(k)) )*0.5 
      q1p=q1c-vmed(1,j)
      q2p=q2c-vmed(2,j)
      q3p=q3c-vmed(3,j)
      vo1p=voc(1,i,j,k)-vorv(1,j)
      vo2p=voc(2,i,j,k)-vorv(2,j)
      vo3p=voc(3,i,j,k)-vorv(3,j)
      q1tz(i,k)=q1p
      q2tz(i,k)=q2p
      q3tz(i,k)=q3p
      prp=pr(i,j,k)-pmed(j)
      prtz(i,k)=prp
      dpr1tz(i,k)=(pr(i,j,k)-pr(imv(i),j,k))*dx1
      dpr3tz(i,k)=(pr(i,j,k)-pr(i,j,kmv(k)))*dx3
      vo1tz(i,k)=vo1p
      vo2tz(i,k)=vo2p
      vo3tz(i,k)=vo3p
      quarm(1,i,k)=q1p/verms(1,j)
      quarm(2,i,k)=q2p/verms(2,j)
      quarm(3,i,k)=q3p/verms(3,j)
      quarm(4,i,k)=prp/prrms(j)
      quarm(5,i,k)=vo1p/vorms(1,j)
      quarm(6,i,k)=vo2p/vorms(2,j)
      quarm(7,i,k)=vo3p/vorms(3,j)
      s11p=sij(1,i,j,k)-strac(1,j)
      s22p=sij(2,i,j,k)-strac(2,j)
      s33p=sij(3,i,j,k)-strac(3,j)
      s23p=sij(4,i,j,k)-strac(4,j)
      s13p=sij(5,i,j,k)-strac(5,j)
      s12p=sij(6,i,j,k)-strac(6,j)
      enep=q1p**2+q2p**2+q3p**2
      enep=q1p**2+q2p**2+q3p**2
      ensp=vo1p**2+vo2p**2+vo3p**2
      helep=(q1p*vo1p+q2p*vo2p+q3p*vo3p)/sqrt(enep*ensp)
      dislo=(s11p**2+s22p**2+s33p**2+
     1      2.*(s12p**2+s13p**2+s23p**2))
      enspro=vo1p**2*s11p+vo2p**2*s22p+vo3p**2*s33p+
     1      2.*(vo1p*vo2p*s12p+vo2p*vo3p*s23p+vo1p*vo3p*s13p)
      prlap=ensp/4.-dislo
      qua(1,i,k)=enep
      qua(2,i,k)=ensp/4.
      qua(3,i,k)=helep
      qua(4,i,k)=dislo
      qua(5,i,k)=enspro
      qua(6,i,k)=prlap
      detsij=s11p*(s22p*s33p-s23p*s23p)
     1      -s12p*(s12p*s33p-s23p*s13p)
     1      +s13p*(s12p*s23p-s22p*s13p)
      qua(7,i,k)=detsij
      enddo
      enddo
c
c   here velocity vorticity components and pressure are written
c   in a file compatible with TURB3D package
c
      namfi3='horseupo'//njpse//'.dat'
      open(59,file=namfi3,form='unformatted')
      write(59) n3m,n1m,1
      write(59) re,re,re,time
      write(59)
     1   ((q1tz(i,k),k=1,n3m),i=1,n1m),
     1   ((q2tz(i,k),k=1,n3m),i=1,n1m),
     1   ((q3tz(i,k),k=1,n3m),i=1,n1m),
     1   ((prtz(i,k),k=1,n3m),i=1,n1m),
     1   ((vo1tz(i,k),k=1,n3m),i=1,n1m),
     1   ((vo2tz(i,k),k=1,n3m),i=1,n1m),
     1   ((vo3tz(i,k),k=1,n3m),i=1,n1m)
      close(59)
c
c  the previous fluctuating quantities normalized
c  with respect to their rms value at this j location
c
      namfi3='horuoprms'//njpse//'.dat'
      open(59,file=namfi3,form='unformatted')
      rewind 59
      write(59) n3m,n1m,1
      write(59) re,re,re,time
      write(59)
     1   ((quarm(1,i,k),k=1,n3m),i=1,n1m),
     1   ((quarm(2,i,k),k=1,n3m),i=1,n1m),
     1   ((quarm(3,i,k),k=1,n3m),i=1,n1m),
     1   ((quarm(4,i,k),k=1,n3m),i=1,n1m),
     1   ((quarm(5,i,k),k=1,n3m),i=1,n1m),
     1   ((quarm(6,i,k),k=1,n3m),i=1,n1m)
     1  ,((quarm(7,i,k),k=1,n3m),i=1,n1m)
      close(59)
c
c  a bunch of other statistics are visualized
c
      namfi3='horseeeh'//njpse//'.dat'
      open(59,file=namfi3,form='unformatted')
      rewind 59
      write(59) n3m,n1m,1
      write(59) re,re,re,time
      write(59)
     1   ((qua(1,i,k),k=1,n3m),i=1,n1m),
     1   ((qua(2,i,k),k=1,n3m),i=1,n1m),
     1   ((qua(3,i,k),k=1,n3m),i=1,n1m),
     1   ((qua(4,i,k),k=1,n3m),i=1,n1m),
     1   ((qua(5,i,k),k=1,n3m),i=1,n1m),
     1   ((qua(6,i,k),k=1,n3m),i=1,n1m)
     1  ,((qua(7,i,k),k=1,n3m),i=1,n1m)
      close(59)
      if(j.eq.1.or.j.eq.n2m) then
      namfi3='horsepre'//njpse//'.dat'
      open(59,file=namfi3,form='unformatted')
      write(59) n3m,n1m,1
      write(59) re,re,re,time
      write(59)
     1   ((prtz(i,k),k=1,n3m),i=1,n1m),
     1   ((dpr1tz(i,k),k=1,n3m),i=1,n1m),
     1   ((dpr3tz(i,k),k=1,n3m),i=1,n1m)
      close(59)
                              endif
      npj=npj+1
                            endif
               enddo
c
c  side-view section  x2-x3 planes
c
      npi=1
      do i=n3m,1,-1
      if(i.eq.iprq(npi)) then
      write(nkpse,83)i
      do k=1,n3m
      do j=1,n2m
      jp=j+1
      q2c=
     1  (    q(2,i,j,k)+q(2,i,jp,k))*0.5
      q1c=
     1  (    q(1,i,j,k)+q(1,ipv(i),j,k) )*0.5
      q3c=(q(3,i,j,k)+q(3,i,j,kpv(k)) )*0.5
      q1p=q1c-vmed(1,j)
      q2p=q2c-vmed(2,j)
      q3p=q3c-vmed(3,j)
      vo1p=voc(1,i,j,k)-vorv(1,j)
      vo2p=voc(2,i,j,k)-vorv(2,j)
      vo3p=voc(3,i,j,k)-vorv(3,j)
      q1rz(j,k)=q1p
      q2rz(j,k)=q2p
      q3rz(j,k)=q3p
      prrz(j,k)=pr(i,j,k)-pmed(j)
      vo1rz(j,k)=vo1p
      vo2rz(j,k)=vo2p
      vo3rz(j,k)=vo3p
      s11p=sij(1,i,j,k)-strac(1,j)
      s22p=sij(2,i,j,k)-strac(2,j)
      s33p=sij(3,i,j,k)-strac(3,j)
      s23p=sij(4,i,j,k)-strac(4,j)
      s13p=sij(5,i,j,k)-strac(5,j)
      s12p=sij(6,i,j,k)-strac(6,j)
      enep=q1p**2+q2p**2+q3p**2
      ensp=vo1p**2+vo2p**2+vo3p**2
      helep=(q1p*vo1p+q2p*vo2p+q3p*vo3p)/sqrt(enep*ensp)
      dislo=(s11p**2+s22p**2+s33p**2+
     1      2.*(s12p**2+s13p**2+s23p**2))
      enspro=vo1p**2*s11p+vo2p**2*s22p+vo3p**2*s33p+
     1      2.*(vo1p*vo2p*s12p+vo2p*vo3p*s23p+vo1p*vo3p*s13p)
      prlap=ensp/4.-dislo
      qua(1,j,k)=enep
      qua(2,j,k)=ensp/4.
      qua(3,j,k)=helep
      qua(4,j,k)=dislo
      qua(5,j,k)=enspro
      qua(6,j,k)=prlap
      detsij=s11p*(s22p*s33p-s23p*s23p)
     1      -s12p*(s12p*s33p-s23p*s13p)
     1      +s13p*(s12p*s23p-s22p*s13p)
      qua(7,j,k)=detsij
      enddo
      enddo
      namfi3='sidseupo'//nkpse//'.dat'
      write(6,*)' written sections aty1=',y1s(i)
      open(59,file=namfi3,form='unformatted')
      write(59) n3m,n2m,1
      write(59) re,re,re,time
      write(59)
     1   ((q1rz(j,k),k=1,n3m),j=1,n2m),
     1   ((q2rz(j,k),k=1,n3m),j=1,n2m),
     1   ((q3rz(j,k),k=1,n3m),j=1,n2m),
     1   ((prrz(j,k),k=1,n3m),j=1,n2m),
     1   ((vo1rz(j,k),k=1,n3m),j=1,n2m),
     1   ((vo2rz(j,k),k=1,n3m),j=1,n2m),
     1   ((vo3rz(j,k),k=1,n3m),j=1,n2m)
          close(59)
      rewind 59
      namfi3='sidseeeh'//nkpse//'.dat'
      write(6,*)' written sections aty1=',y1s(i)
      open(59,file=namfi3,form='unformatted')
      write(59) n3m,n2m,1
      write(59) re,re,re,time
      write(59)
     1   ((qua(1,j,k),k=1,n3m),j=1,n2m),
     1   ((qua(2,j,k),k=1,n3m),j=1,n2m),
     1   ((qua(3,j,k),k=1,n3m),j=1,n2m),
     1   ((qua(4,j,k),k=1,n3m),j=1,n2m),
     1   ((qua(5,j,k),k=1,n3m),j=1,n2m),
     1   ((qua(6,j,k),k=1,n3m),j=1,n2m)
     1  ,((qua(7,j,k),k=1,n3m),j=1,n2m)
      close(59)
      npi=npi+1
                            endif
                 enddo

      return
      end
c
c  ****************************** subrout coordi  **********************
c    this subroutine assign the physical coordinates as function
c    of the computational ones x2
c
      subroutine coordi(y)
c
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      dimension x(m1),eta(m2)
      common/strpar/str2
      dimension y(m2)
      common/d13/alx1,alx3
      common/y2sta/y2s(m2)
      common/y13sta/y1s(m1),y3s(m3)
      common/corpri/yp1(m1),yp2(m2),yp3(m3)
      common/jprts/jpri,ipri,kpri
c
      tstr2=tanh(str2*0.5)
      do 63 j=1,n2
      x2=(j-1)/float(n2m)
      eta(j)=0.5*(1.+tanh(str2*(x2-0.5))/tstr2)
   63 continue
      do 65 j=1,n2
      yp2(j)=(-0.5+eta(j))*2.
      y(j)=yp2(j)
   65 continue
c     print *,'y al centro',y(n2/2+1)
      do 67 j=1,n2m
      y2s(j)=(yp2(j)+yp2(j+1))*0.5
   67 continue
      if(mod(n2,2).eq.1) then
      print *,'n2 dispari forza la simmetria nella coordinata y' 
      y(n2/2+1)=0.
      do j=1,n2/2
      y(n2-j+1)=-y(j)
      end do
      do j=1,n2/2
      y2s(n2-j)=-y2s(j)
      end do
      end if
      open(47,file='ymesh')
      do j=1,n2m
      dy=y(j+1)-y(j) 
      write(47,*) y2s(j),dy
      end do
      close(47)
      do j=1,n2m/2+1
      ydp=(y(j)-y(1) )*180
      if(ydp.gt.20) then
      jpri=j
      go to 23
                    endif
      end do
   23 continue
      kpri=n3m/2+1
      ipri=n1m/2+1
      jin=jpri-3
      jfi=jpri+3
      iin=ipri-3
      ifi=ipri+3
      kin=kpri-3
      kfi=kpri+3
      npoin=0
      do i =iin,ifi
      do j =jin,jfi
      do k =kin,kfi
      npoin=npoin+1
      enddo
      enddo
      enddo
      write(6,*)' print npoin=',npoin,'   around j=',jpri
c
      do  i=1,n1
          yp1(i) = float(i-1)/float(n1m) * alx1
      enddo
      do  k=1,n3
          yp3(k) = float(k-1)/float(n3m) * alx3
      enddo
      do  i=1,n1m
          y1s(i) = 0.5*(yp1(i) + yp1(i+1) )
      enddo
      do  k=1,n3m
          y3s(k) = 0.5*(yp3(k) + yp3(k+1) )
      enddo
      return
      end

c  ****************************** subrout divgck  **********************
c
c  this subroutine perform the calculation of divg(q)
c
      subroutine divgck(vq,qmax)
c
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      dimension vq(ndv,m1,m2,m3)
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/metria/caj(m2),cac(m2)
c
c  ***** compute the divg(q)
      qmax=0.
      do 11 kc=1,n3m
      kp=kpv(kc)
      do 11 jc=1,n2m
      jp=jc+1
      sucaj=1./caj(jc)
      do 11 ic=1,n1m
      ip=ipv(ic)
      dqcap=(vq(1,ip,jc,kc)-vq(1,ic,jc,kc))*dx1
     1     +(vq(2,ic,jp,kc)-vq(2,ic,jc,kc))*dx2*sucaj
     1     +(vq(3,ic,jc,kp)-vq(3,ic,jc,kc))*dx3
      qmax=max(abs(dqcap),qmax)
   11 continue
      return
      end
c
c  ****************************** subrout indic **********************
c
c  in this subroutine the indices ip,im,jp,jm,kp,km are calculated
c
      subroutine indic
c
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/injmp/jpppv(m2),jmmmv(m2)
      common/inimp/ipppv(m1),immmv(m1)
c
c
c   periodic directions
c
      n1mm=n1m-1
      do 11 ic=1,n1m
      ipv(ic)=ic+1
      if(ic.eq.n1m) ipv(ic)=1
      imv(ic)=ic-1
      if(ic.eq.1) imv(ic)=n1m
      ip=ipv(ic)
      ipppv(ic)=ip+1
      if(ic.eq.n1mm) ipppv(ic)=1
      im=imv(ic)
      immmv(ic)=im-1
      if(ic.eq.2) immmv(ic)=n1m
   11 continue
      do 2 kc=1,n3m
      kmv(kc)=kc-1
      kpv(kc)=kc+1
      if(kc.eq.1) kmv(kc)=n3m
      if(kc.eq.n3m) kpv(kc)=1
    2 continue
c
c   direction normal to non-slip walls
c
      do 3 jc=1,n2m
      jp=jc+1
      jpppv(jc)=jp+1
      if(jc.eq.n2m) jpppv(jc)=n2
      jmv(jc)=jc-1
      jpv(jc)=jc+1
      if(jc.eq.1) jmv(jc)=jc
      if(jc.eq.n2m) jpv(jc)=jc
    3 continue
      do 5 jc=2,n2m
      jmmmv(jc)=jc-2
      if(jc.eq.2) jmmmv(jc)=1
    5 continue
      return
      end
c
c
c  ****************************** subrout meshes **********************
c
c  generates the mesh inverse of the spatial steps dx1, dx2, dx3,
c  and the squares of the dx1 and dx2.
c
      subroutine meshes
c
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      common/d13/alx1,alx3
      common/meshu/d1x1,d1x2,d1x3
c
      d1x2=1./float(n2m)
      d1x1=alx1/float(n1m)
      d1x3=alx3/float(n3m)
      dx1=1./d1x1
      dx2=1./d1x2
      dx3=1./d1x3
      dx1q=dx1*dx1
      dx2q=dx2*dx2
      dx3q=dx3*dx3
      return
      end
c
c
c  ****************************** subrout metric  **********************
c
c  this subroutine performs the calculation of the metric quantities
c  at j ; j+1/2 
c  ** the quantities with  j(  evaluated at .,j+1/2,.     positiion
c  ** the quantities with  c(  evaluated at .,j,.         positiion
c  ** the metric are evaluated by centered differences.
c
      subroutine metric(y)
c
      include 'param.f'
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      common/metria/caj(m2),cac(m2)
      dimension y(m2)
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/d13/alx1,alx3
c
c  *********                 i,j+1/2
c
      do 1 j=1,n2m
c
c  ********* derivatives of cartesian corrdinates interior points
c  cn1 deriv respect to x1 , cn2 deriv respect to x2
c
      cn22=(y(j+1)-y(j))*dx2
   10 continue
      caj(j)=cn22
    1 continue
c
c  *********                 i,j
c
      do 6 j=1,n2
      if(j.gt.1.and.j.lt.n2) go to 31
c
c  ********* derivatives of cartesian corrdinates horizontal boundary
c            three points backw.
      if(j.eq.1) js=1
      if(j.eq.n2) js=-1
      cn22=-js*(y(j)-y(j+js))*dx2
      go to 34
   31 continue
c
c  ********* derivatives of cartesian corrdinates interior points
c
      cn22=(y(j+1)-y(j-1))*dx2*0.5
   60 continue
   34 continue
      cac(j)=cn22
    6 continue
      return
      end
c
c  ****************************** subrout inirea **********************
c
c   reading initial  conditions from file
c
      subroutine inirea(ntii,time,q,ru,pr)
c
      include 'param.f'
      dimension pr(m1,m2,m3)
      common/velmax/vmax(ndv),vmaxo(ndv)
      common/dim/n1,n1m,n2,n2m,n3,n3m
      dimension q(ndv,m1,m2,m3)
      dimension ru(ndv,m1,m2,m3)
      common/veltot/vit(4)
      common/omegas/omr,omi
      common/inener/ene0
      common/rhs3p/dp3ns
      common/eneav/enav,enavo,cfn,cfo
      common/metria/caj(m2),cac(m2)
      common/d13/alx1,alx3
      common/oldso/ifield
      common/qwallo/q1s(m1,m3),q2s(m1,m3),q3s(m1,m3)
      common/qwalup/q1n(m1,m3),q2n(m1,m3),q3n(m1,m3)
      character*60 namfile
c
      character*4 ipfi
      nfil=23
      itime=nint(time)
      write(ipfi,82)itime
   82 format(i4.4)
c
c   this file was written during the simulation
c   by the code in the directory CODESCH9/CHA
c   This file is a restarting file
c   of large dimensions since the pressure and
c   q2 are given. These are not necessary
c   since q2 can be obtained by the divg if
c   only q1 and q3 re read the pressure can be
c   obtained by performing one time step advancement
c   by the code in CODESCH9/CHA
c
      namfile='../field'//ipfi//'.dat'
      open(23,file=namfile,form='unformatted')
      rewind(nfil)
      read(nfil)n1lm,n2,n3lm
      read(nfil)ntime,time,ene0,dp3ns,enavo,cfo,dtold
      read(nfil) (((q(1,i,j,k),i=1,n1lm),j=1,n2),k=1,n3lm),
     1            (((q(2,i,j,k),i=1,n1lm),j=1,n2),k=1,n3lm),
     1            (((q(3,i,j,k),i=1,n1lm),j=1,n2),k=1,n3lm),
     1            (((pr(i,j,k),i=1,n1lm),j=1,n2),k=1,n3lm)
      read(nfil)ntime,time,vit(1),vit(2),vit(3)
     1         ,dp3ns,cfuw,enen,vmax(2),vmax(3),cflm
      if (ifield.eq.1) then
      read(nfil)  ((q1s(i,k),i=1,n1lm),k=1,n3lm),
     1            ((q2s(i,k),i=1,n1lm),k=1,n3lm),
     1            ((q3s(i,k),i=1,n1lm),k=1,n3lm),                        
     1            ((q1n(i,k),i=1,n1lm),k=1,n3lm),                        
     1            ((q2n(i,k),i=1,n1lm),k=1,n3lm),                        
     1            ((q3n(i,k),i=1,n1lm),k=1,n3lm)                         
      end if
      close(nfil)
      call vmaxv(q,n1lm,n2,n3lm)
      write(6,*)' in inirea vmax(l=',(vmax(l),l=1,3)
      call divgck(q,qmax)
      write(6,*)' in inirea after read divg max=',qmax
      return
      end
