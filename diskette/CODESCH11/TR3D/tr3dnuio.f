c************************************************************************
c                                                                       *
c     ********* subrout pricor ******************                       *
c                                                                       *
c************************************************************************
      subroutine pricor
      include 'param.f'
      dimension zzmm(m3),xxm(m1,m2),yym(m1,m2),rh(2*m2)
      n3pp=(n3-1)/n3p
      n2pp=(n2-1)/n2p
      n1pp=(n1-1)/n1p+1
      aaa=1.
c
c   a coordinate  file is written to investigate
c   how the simulation evolves by looking at
c   contours in r-z planes obtained by averaging
c   in theta
c
      do k=1,n3m
      zzmm(k)=(zz(k)+zz(k+1))*0.5
      enddo
      namfile='cordaxsym.dat'
      open(18,file=namfile,form='unformatted')
      write(18) n3m,n2m,1
      write(18)
     1            ((zzmm(k),k=1,n3m),j=1,n2m),
     1            ((rm(j),k=1,n3m),j=1,n2m),
     1            ((aaa,k=1,n3m),j=1,n2m)
      close(18)
      return
      end
c************************************************************************
c                                                                       *
c  ****************************** subrout outpf  ********************** *
c                                                                       *
c************************************************************************
      subroutine outpf(time,nav)
      include 'param.f'
      character*3 nkpse,njpse
      character*4 navps
      dimension quap(7,m3,m2),qc(7)
c
c   averages velocity  in theta
c
      call veltur
c
c   averages vorticity  in theta
c
      call vortur
      itim=nint(time)
      write(nkpse,83)itim
   83 format(i3.3)
c
c   the averages are written in a file to be visualized
c          by TURB3D  package
c
      namfi3='qmed'//nkpse//'.dat'
      write(6,*)' after ntime=',nav,'  print mean quant. on',namfi3
      do k=1,n3m
      do j=1,n2m
      quap(1,k,j)=vmed(1,k,j)
      quap(2,k,j)=vmed(2,k,j)
      quap(3,k,j)=vmed(3,k,j)
      quap(4,k,j)=pscmed(k,j)
      quap(5,k,j)=vorv(1,k,j)
      quap(6,k,j)=vorv(2,k,j)
      quap(7,k,j)=vorv(3,k,j)
      enddo
      enddo
      open(59,file=namfi3,form='unformatted')
      rewind 59
      write(59) n3m,n2m,1
      write(59) re,re,re,time
      write(59)
     1   ((quap(1,k,j),k=1,n3m),j=1,n2m),
     1   ((quap(2,k,j),k=1,n3m),j=1,n2m),
     1   ((quap(3,k,j),k=1,n3m),j=1,n2m),
     1   ((quap(4,k,j),k=1,n3m),j=1,n2m),
     1   ((quap(5,k,j),k=1,n3m),j=1,n2m),
     1   ((quap(6,k,j),k=1,n3m),j=1,n2m),
     1   ((quap(7,k,j),k=1,n3m),j=1,n2m)
      close(59)
c
c    these commented instructions permit to write
c    files with the rms distributions
c     of fluctuating velocity and vorticity
c
c     namfi3='veturb'//nkpse//'.dat'
c     do k=1,n3m
c     do j=1,n2m
c     quap(1,k,j)=turstr(1,k,j)
c     quap(2,k,j)=turstr(2,k,j)
c     quap(3,k,j)=turstr(3,k,j)
c     quap(4,k,j)=turstr(4,k,j)
c     quap(5,k,j)=turstr(5,k,j)
c     quap(6,k,j)=turstr(6,k,j)
c     quap(7,k,j)=
c    1      (turstr(1,k,j)**2+turstr(2,k,j)**2+turstr(3,k,j)**2)
c     enddo
c     enddo
c     open(59,file=namfi3,form='unformatted')
c     rewind 59
c     write(59) n3m,n2m,1
c     write(59) re,re,re,time
c     write(59)
c    1   ((quap(1,k,j),k=1,n3m),j=1,n2m),
c    1   ((quap(2,k,j),k=1,n3m),j=1,n2m),
c    1   ((quap(3,k,j),k=1,n3m),j=1,n2m),
c    1   ((quap(4,k,j),k=1,n3m),j=1,n2m),
c    1   ((quap(5,k,j),k=1,n3m),j=1,n2m),
c    1   ((quap(6,k,j),k=1,n3m),j=1,n2m),
c    1   ((quap(7,k,j),k=1,n3m),j=1,n2m)
c     close(59)
c     namfi3='voturb'//nkpse//'.dat'
c     do k=1,n3m
c     do j=1,n2m
c     quap(1,k,j)=vorstr(1,k,j)
c     quap(2,k,j)=vorstr(2,k,j)
c     quap(3,k,j)=vorstr(3,k,j)
c     quap(4,k,j)=vorstr(4,k,j)
c     quap(5,k,j)=vorstr(5,k,j)
c     quap(6,k,j)=vorstr(6,k,j)
c     quap(7,k,j)=enpsc(k,j)
c     enddo
c     enddo
c     open(59,file=namfi3,form='unformatted')
c     rewind 59
c     write(59) nn3m,n2m,1
c     write(59) re,re,re,time
c     write(59)
c    1   ((quap(1,k,j),k=1,n3m),j=1,n2m),
c    1   ((quap(2,k,j),k=1,n3m),j=1,n2m),
c    1   ((quap(3,k,j),k=1,n3m),j=1,n2m),
c    1   ((quap(4,k,j),k=1,n3m),j=1,n2m),
c    1   ((quap(5,k,j),k=1,n3m),j=1,n2m),
c    1   ((quap(6,k,j),k=1,n3m),j=1,n2m),
c    1   ((quap(7,k,j),k=1,n3m),j=1,n2m)
c     close(59)
      return
      end
c************************************************************************
c                                                                       *
c    *****   subro outh   ********************************              *
c                                                                       *
c************************************************************************
      subroutine outh(time,nav,ntime,cflm,nvv,navbu)        
      include 'param.f'
      common/pscmm/pscmax,pscmin
      call vmaxv
      call divgck(qma,tma)
c
c   time history of global quantities written
c
      write(33,783)time,vomax(1),vomax(2)
      write(34,783)time,vmax(2)
      write(49,783)time,vomin(3)
      write(40,783)time,pscmin,pscmax
      write(32,783)time,vmax(1),vmax(2),vmax(3),vomax(1)
     1        ,vomax(2),vomax(3),vomin(3),cflm
      write(6,783)time,vmax(1),vmax(2),vmax(3),vomax(1)
     1        ,vomax(2),vomax(3),vomin(3),cflm,qma,pscmax
  783 format(3x,11(1x,e12.5))
      return                                                            
      end                                                               
c************************************************************************
c                                                                       *
c  ****************************** subrout vmaxv **********************  *
c                                                                       *
c************************************************************************
      subroutine vmaxv
      include 'param.f'
      common/pscmm/pscmax,pscmin
c
c     find the maximum velocities in order to check convergence or
c     to derive stability conditions.
c
      vca=0.
      pscmin=100.
      pscmax=-1000.
c
      do l=1,ndv
      vmax(l)=-100.
      vmin(l)=+100.
      enddo
      do k=1,n3m
      do j=2,n2m
      do i=1,n1m
      vca2=q2(i,j,k)/rc(j)
      vmax(2)=max(vca2,vmax(2))
      vmin(2)=min(vca2,vmin(2))
      enddo
      enddo
      enddo
      do k=1,n3m
      do j=1,n2m
      do i=1,n1m
      vca1=q1(i,j,k)
      vca3=q3(i,j,k)
      vmax(1)=max(vca1,vmax(1))
      vmin(1)=min(vca1,vmin(1))
      vmax(3)=max(vca3,vmax(3))
      vmin(3)=min(vca3,vmin(3))
      pscmax=max(psc(i,j,k),pscmax)
      pscmin=min(psc(i,j,k),pscmin)
      enddo
      enddo
      enddo
      return
      end
c************************************************************************
c                                                                       *
c  ****************************** subrout contwr ********************** *
c     print flow fields   to be used as restarting file
c     or to perform flow visualizations by the code in
c     the directory CODESCH11/VIS
c                                                                       *
c************************************************************************
      subroutine contwr(ntime,time,ntt,nap)             
      include 'param.f'
      character*20 filcnwp
      do i=1,n1m
      do j=1,n2
      q1(i,j,n3)=q1(i,j,1)
      q2(i,j,n3)=q2(i,j,1)
      q3(i,j,n3)=q3(i,j,1)
      pr(i,j,n3)=pr(i,j,1)
      psc(i,j,n3)=psc(i,j,1)
      enddo
      enddo
      do j=1,n2
      do k=1,n3
      q1(n1,j,k)=q1(1,j,k)
      q2(n1,j,k)=q2(1,j,k)
      q3(n1,j,k)=q3(1,j,k)
      pr(n1,j,k)=pr(1,j,k)
      psc(n1,j,k)=psc(1,j,k)
      enddo
      enddo
      itime=nint(time)
      write(ipfi,82)itime
   82 format(i4.4)
      namfi3='field'//ipfi//'.dat'
      open(13,file=namfi3,form='unformatted')
      nfil=13
      rewind(nfil)
      write(nfil) n1,n2,n3
      write(nfil) ros,alx3d,re,time
      write(nfil) (((q1(i,j,k),i=1,n1),j=1,n2),k=1,n3),           
     1            (((q2(i,j,k),i=1,n1),j=1,n2),k=1,n3),        
     1            (((q3(i,j,k),i=1,n1),j=1,n2),k=1,n3),
     1            (((psc(i,j,k),i=1,n1),j=1,n2),k=1,n3),
     1            (((pr(i,j,k),i=1,n1),j=1,n2),k=1,n3)     
      write(nfil) ntime,ntt,nap
      close(nfil)                                                       
      return                                                            
      end                                                               
c************************************************************************
c                                                                       *
c  ****************************** subrout inirea ********************** *
c     read restarting  file to continue the simulation                  *
c                                                                       *
c************************************************************************
      subroutine inirea(ntii,time,ntt,ncount,nap)
      include 'param.f'
      common/avgin/avpscn(m2),vnew(m2)
      open(13,file=filcnr,form='unformatted')
      nfil=13                                                           
      read(nfil) n1l,n2l,n3l                                        
      read(nfil) epsil,lamb,rele,time 
      write(6,*)' legge da ',filcnr,'  n1l,n2l,n3l ',n1l,n2l,n3l
      read(nfil)  (((q1(i,j,k),i=1,n1l),j=1,n2l),k=1,n3l),           
     1            (((q2(i,j,k),i=1,n1l),j=1,n2l),k=1,n3l),        
     1            (((q3(i,j,k),i=1,n1l),j=1,n2l),k=1,n3l),
     1            (((psc(i,j,k),i=1,n1l),j=1,n2l),k=1,n3l),
     1            (((pr(i,j,k),i=1,n1l),j=1,n2l),k=1,n3l)      
      call veltur
      call vortur
      return                                                            
      end                                                               
c
c  ****************************** subrout initia  **********************
c
      subroutine initia
      include 'param.f'
      do j=1,n2
      do i=1,n1
      do k=1,n3
      pr(i,j,k)=0.
      q1(i,j,k)=0.
      q2(i,j,k)=0.
      q3(i,j,k)=0.
      ru1(i,j,k)=0.
      ru2(i,j,k)=0.
      ru3(i,j,k)=0.
      enddo
      enddo
      enddo
      return
      end

