c
c  ****************************** subrout vorcal **********************
c
c   from the streamfunction the vorticity is calculated
c
      subroutine vorcal(vor,psi)
      include 'param.f'
      dimension psi(m1,m2),vor(m1,m2)
      common/dim/n1,n1m,n2,n2m
      common/tstep/dt
      common/cnlev/n1g(0:nlev),n2g(0:nlev),n12g(0:nlev)
      common/comg/co(10,nij)
      igjg(i,j,l)=i+(j-1)*n1g(l)+n12g(l-1)
      do  ic=2,n1m
      ip=ic+1
      im=ic-1
      do  jc=2,n2m
      jm=jc-1
      jp=jc+1
      ij=igjg(ic,jc,1)
c
c    cross terms
c
      psidcr= +co(1,ij)*psi(ic+1,jc+1)+
     1         co(3,ij)*psi(ic+1,jc-1)+
     1         co(5,ij)*psi(ic-1,jc-1)+
     1         co(7,ij)*psi(ic-1,jc+1)
c
c    diagonal terms
c
      psider=+ co(2,ij)*psi(ic+1,jc)+
     1         co(4,ij)*psi(ic,jc-1)+
     1         co(6,ij)*psi(ic-1,jc)+
     1         co(8,ij)*psi(ic,jc+1)+
     1         co(9,ij)*psi(ic,jc)
      vor(ic,jc)=-(psidcr+psider)
      enddo
      enddo
      return
      end
c
c  ****************************** subrout carmod **********************
c
c  this subroutine calculate the chjaracteristic quantities of the dipole.
c  total vorticity and enstrophy and the position of maxima vorticity
c
      subroutine carmod(vor,psi,y)
      include 'param.f'
      common/mesh/dx1,dx1q,dx2,dx2q
      common/metrst/gccc(m1,m2)
      dimension y(ndd,m1,m2)
      dimension vor(m1,m2),psi(m1,m2)
      common/dim/n1,n1m,n2,n2m
      common/camo/y1ma,y2ma,y1mi,y2mi,vori,ensti
      common/circul/vorip,vorim
      vorip=0.
      vorim=0.
      ensti=0.
      i=1
      do 410 i=1,n1m
      do 410 j=2,n2m
      if(vor(i,j).lt.0.) then
      vorim=vorim+vor(i,j)/(dx1*dx2)*gccc(i,j)
                         endif
      if(vor(i,j).gt.0.) then
      vorip=vorip+vor(i,j)/(dx1*dx2)*gccc(i,j)
                         endif
      vori=vorip+vorim
      ensti=ensti+vor(i,j)**2/(dx1*dx2)*gccc(i,j)
  410 continue
      vormi=1.e+05
      vorma=-1.e+05
      do 411 i=1,n1
      do 411 j=1,n2
      if(vor(i,j).gt.vorma) then
      ima=i
      jma=j
      vorma=vor(i,j)
      endif
      if(vor(i,j).lt.vormi) then
      imi=i
      jmi=j
      vormi=vor(i,j)
      endif
  411 continue
      y1ma=y(1,ima,jma)
      y2ma=y(2,ima,jma)
      y1mi=y(1,imi,jmi)
      y2mi=y(2,imi,jmi)
      return
      end
c
c
c  ****************************** subrout outpf  **********************
c   the solution is written in a file for visualizations
c   format as requested by TURB3D
c
      subroutine outpf(vor,psi,time,psc)
      include 'param.f'
      common/d2/nstop,nprint,ntst,npin,nmolp,ntra,mnpi
      common/mesh/dx1,dx1q,dx2,dx2q
      dimension vor(m1,m2),psi(m1,m2)
      dimension psc(m1,m2)
      common/dim/n1,n1m,n2,n2m
      character*4 pntim
      character*60 namfile
            itim=nint(time*nmolp)
      write(pntim,77) itim
   77 format(i4.4)
      namfile='field'//pntim//'.dat'
      open(17,file=namfile,form='formatted')
      write(17,*)n1,n2
      write(17,*) time,time,time,time
      write(17,*)   ((psi(i,j),i=1,n1),j=1,n2)
     1           ,((psc(i,j),i=1,n1),j=1,n2)
     1           ,((psc(i,j),i=1,n1),j=1,n2)
     1           ,((vor(i,j),i=1,n1),j=1,n2)
      close(17)
      return
      end
c
c
c  ****************************** subrout enerca **********************
c
c   calculation of total energy =sum psi(i,j)*vor(i,j)
c
      subroutine enerca(vor,psi,enej,y,psc)
      include 'param.f'
      common/mesh/dx1,dx1q,dx2,dx2q
      dimension vor(m1,m2),psi(m1,m2),y(ndd,m1,m2),psc(m1,m2)
      common/dim/n1,n1m,n2,n2m
      common/metrst/gccc(m1,m2)
      common/vpmima/vorma,vormi,pscma,pscmi
      enej=0.
      pscma=-100.
      vorma=-100.
      pscmi=+100.
      vormi=+100.
      do j=2,n2m
      do i=2,n1m
      enej=enej+psi(i,j)*vor(i,j)/(dx1*dx2)*gccc(i,j)
      enddo
      enddo
      do j=1,n2
      do i=1,n1
      vorma=max(vorma,vor(i,j))
      pscma=max(pscma,psc(i,j))
      vormi=min(vormi,vor(i,j))
      pscmi=min(pscmi,psc(i,j))
      enddo
      enddo
      return
      end
c
c  ****************************** subrout inirea **********************
c
c   initial conditions in the whole field for psi and vor
c  
      subroutine inirea(vor,psi,time,ntime,psc)
      include 'param.f'
      dimension vor(m1,m2),psi(m1,m2),psc(m1,m2)
      common/dim/n1,n1m,n2,n2m
      open(17,file='../inirea.dat',form='formatted')                           
      read(17,*)n1,n2                                        
      read(17,*) time,time,time,time                             
      read(17,*)   ((psi(i,j),i=1,n1),j=1,n2)       
     1           ,((psc(i,j),i=1,n1),j=1,n2)
     1           ,((psc(i,j),i=1,n1),j=1,n2)                      
     1           ,((vor(i,j),i=1,n1),j=1,n2)                         
      close(17)                                                               
      return
      end
c
c  ****************************** subrout meshes **********************
c
      subroutine meshes
      include 'param.f'
      common/dim/n1,n1m,n2,n2m
      common/mesh/dx1,dx1q,dx2,dx2q
      dx1=1./float(n1m)
      dx2=1./float(n2m)
      dx1=1./dx1
      dx2=1./dx2
      dx1q=dx1*dx1
      dx2q=dx2*dx2
      return
      end
c  ****************************** subrout cfl  **********************
c
c  in this subroutine is calculated the maximum courant number
c
      subroutine cfl(psi,cflm)
      include 'param.f'
      common/dim/n1,n1m,n2,n2m
      common/mesh/dx1,dx1q,dx2,dx2q
      dimension psi(m1,m2)
      common/chnlc/chal,chbe,chga
      common/metrst/gccc(m1,m2)
      common/metrih/cac(ndd,ndd,m1,m2)
      cflm=0.
      do 10 ic=2,n1m
      ip=ic+1
      im=ic-1
      do 10 jc=2,n2m
      jm=jc-1
      jp=jc+1
c
c     terms as Arakawa 
c
      h22a= ((psi(ic,jp)-psi(ic,jm))-
     1       (psi(ip,jc)-psi(im,jc)))
     1     *dx2/gccc(ic,jc)*dx1*0.5
      cflm=max(abs(h22a),cflm)
   10 continue
      return
      end
