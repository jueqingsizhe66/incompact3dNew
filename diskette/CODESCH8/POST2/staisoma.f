c
c   in thic code the field generated in the simulation are read and
c   from this many statistics are evaluated
c
      program main
      common/d1/re
      common/d2/nstop,nprint,ntst,npin,npstf
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/wrre/nwrit,nread
      common/tscoe/ga(3),ro(3),nsst
      common/d123/alx1,alx2,alx3
      common/averou/iav
      common/inior/indrea,icosma,icont
      common/ledat/ics0,cvisc,ifiltr,ibox
      common/cosma/csma
      common/cflco/icfl,cflc,tpin,tprin,tfin
      common/rot/f0
      common/ispec/imic
      common/vpeini/vper
      common/spemw/akkpp,qq,sig
      common/tstep/dt,beta,ren
      common/pardip/thet0,vsi,yc1mo,yc2mo,akmo,velmo
      common/anmod/ar1,ar2
      common/papdf/nlr,nlrju,lipdf,lipdh,lipdn
      common/paqdf/sfuma,qquma,sfoma,qqoma
      common/alpdf/alasfu(201),alaqu(201),alahe(201),aladi(201)
      common/itystf/inostf
      pi=2.*asin(1.)
c
      open(15,file='isostat.d')
      read(15,*) n1,n2,n3,nsst,alx1d,alx2d,alx3d
      read(15,*) nwrit,nread,iav
      read(15,*) re,vper,dt,ntst,nprint,npin ,npstf
      read(15,*) nstop
      read(15,*) icfl,cflc,tpin,tprin,tfin
      read(15,*) ics0,ifiltr,csma,ibox
      read(15,*) f0
c
c rotation to add to the momentum equation
c coriolis parameter f0, rotation about z-axis i.e.  axis 3
c
      read(15,*) imic
      if(imic.ge.0) then
      read(15,*) akkpp,qq,sig
                    else
      read(15,*) yc1mo,yc2mo,ramo,velmo,vsi
      read(15,*) thet0
      read(15,*) ar1,ar2
      akmo=3.83711/ramo
      pi=2.*asin(1.)
      thet0=thet0*pi/180.
      yc1mo=yc1mo*pi
      yc2mo=yc2mo*pi
      write(6,201)thet0
  201 format(3x,'modone vort stream funct. tht0=',e10.4)
                    endif
      read(15,*)nlr,nlrju,lipdf,sfuma,qquma,sfoma,qqoma
      read(15,*)inostf

      mpq=npq
      cvisc=1./re
      n1m=n1-1
      n2m=n2-1
      n3m=n3-1
      alx1=alx1d*pi
      alx2=alx2d*pi
      alx3=alx3d*pi
      call openfi
      call pospro
      stop
      end
c
      subroutine openfi
      common/names/filcnw,filcnr,filth,filou,filuu,filsf
      character*27 filcnw,filcnr,filth,filou,filuu,filsf
      open(46,file='nfisodyn')
c     continuation file
      read(46,'(a)')filcnw
c     restart file
      read(46,'(a)')filcnr
c     time history
      read(46,'(a)')filth
c     rms file
      read(46,'(a)')filou
c     skewness
      read(46,'(a)')filuu
c     flatness
      read(46,'(a)')filsf
      open(23,file=filcnr,form='unformatted')
      open(32,file=filth)
      rewind 23
      rewind 32
      return
      end
c
c  ************************************************************
c  ************subrout pospro **********************************
c  ************************************************************
c
      subroutine pospro
c
c     code to post-process the velocity field to get statistics
c     several routines used into the code to generate the
c     data are used.
c
      include 'param.f'
      parameter (m1m=m1-1)
c
      common/cflco/icfl,cflc,tpin,tprin,tfin
      common/wrre/nwrit,nread
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      dimension q(ndv,m1,m2,m3),pr(m1,m2,m3)
      dimension qcap(m1,m2,m3)
      common/velmax/vmax(ndv),vmaxo(ndv)
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/d1/re
      common/d123/alx1,alx2,alx3
      common/d2/nstop,nprint,ntst,npin,npstf
      common/tstep/dt,beta,ren
      common/tscoe/ga(3),ro(3),nsst
      common/inener/ene0
      common/averou/iav
      common/sc/sc,sc1
      common/names/filcnw,filcnr,filth,filou,filuu,filsf
      common/jrrr/jri,jrf,djr,irejr,iruuca
      common/ledat/ics0,cvisc,ifiltr,ibox
      common/newdat/icost,timeav
      common/cosma/csma
      common/rot/f0
      common/speene/e(ndv,0:m3)
      common/numond/kx(m1),ky(m2),kz(m3),kkmax
      common/ispec/imic
      dimension vor(m1,m2)
      character*27 filcnw,filcnr,filth,filou,filuu,filsf
      character*4 pntim
c
c
      pi=2.*asin(1.)
c
c     step  and mesh sizes calculations
c
c
      call meshes
      call indic
      call coordi
c
      write(6,754)n1,n2,n3
  754 format(10x,'centered velocities',2x
     1      ,5x,'n1=',i4,2x,'n2=',i4,2x,'n3=',i4)
      write(6,755) re,dx1,dx2,dx3,dt,ntst
  755 format(3x,'re=',e10.3,3x,'dx1=',e10.3,3x,'dx2=',e10.3,3x
     1      ,3x,'dx3=',e10.3,3x,
     1      'dt=',e10.3,3x,'ntst=',i5)
      write(6,756)  f0
  756 format(3x,'f0=',e10.3)
      ren=re
c
c   here the velocity field is read
c
      call inirea(ntii,time,q,pr)
      print *,'the field was correctly read'
c
c    check of divergence 
c
      call divgck(q,qmax)
      print *,'divergence of the field is =',qmax
      call vmaxv(q)
      write(6,*)'  vmax',(vmax(l),l=1,3)    
      write(57,*)'  vmax',(vmax(l),l=1,3)    
c
c  calculation of quantities for FFT
c
      call speini
      call pdfini   
c
c   OUTPUT
c
      call outpf(time,enen,nav,q,pr,qcap)
      return
      end
