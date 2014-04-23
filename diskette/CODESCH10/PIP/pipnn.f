c************************************************************************
c                                                                       *
      subroutine openfi                                                 *
c                                                                       *
c************************************************************************
      include 'param.f'
c
c   a bunch of files are opened to write
c   the outputs
c
      open(46,file='nfpipe')
      read(46,'(a)')filcnw
      read(46,'(a)')filcnr
      read(46,'(a)')filth
      read(46,'(a)')filvm
      read(46,'(a)')filpo
      read(46,'(a)')filen
      read(46,'(a)')filet
      read(46,'(a)')filer
      read(46,'(a)')filez
      read(46,'(a)')filed
      read(46,'(a)')filev
      open(32,file=filth)
      open(33,file=filvm)
      open(34,file=filpo)
      open(39,file=filen)
      open(40,file=filet)
      open(41,file=filer)
      open(42,file=filez)
      open(49,file=filed)
      open(50,file='piqm.out')
      open(59,file='piav.out')
      rewind 12
      rewind 33
      rewind 34
      rewind 32
      rewind 39
      rewind 40
      rewind 41
      rewind 42
      rewind 49
      rewind 48
      if(iwlop.gt.0) then
      numf=90
      kpr=n3m/12
      kpi=n3m/12+1
      ipr=n1m/4
      do jp=1,njprs
      j=npjp(jp)
      write(6,*) jp,j
c     do k=kpi,n3m,kpr
c     write(ipfk,82)k
      write(ipfj,82)j
   82 format(i3.3)
c     namfi3='fk'//ipfk//'j'//ipfj//'.dat'
      namfi3='fj'//ipfj//'.dat'
      open(numf,file=namfi3,form='unformatted')
      numf=numf+1
c     enddo
      enddo
      write(6,*)numf
                      endif
      return
      end   
