c
c  ****************************** subrout itsolv **********************
c this is an iterative solver for the multigrid
c here there is the POINT SOR
c  this is not appropriate for CRAY computers
c
c
      subroutine itsolv(nl,ph,res)
      include 'param.f'
      common/cnlev/n1g(0:nlev),n2g(0:nlev),n12g(0:nlev)
      dimension ph(nij),res(nij)
      common/comg/co(10,nij)
      common/maxit/nmaxi
      common/omsor/omeg
      igjg(i,j,l)=i+(j-1)*n1g(l)+n12g(l-1)
      nit=0
   50 continue
      do jc=2,n2g(nl)-1
      do ic=2,n1g(nl)-1
      ij=igjg(ic,jc,nl)
      icjp=igjg(ic,jc+1,nl)
      icjm=igjg(ic,jc-1,nl)
      ipjp=icjp+1
      ipjm=icjm+1
      imjm=icjm-1
      imjp=icjp-1
      ipjc=ij+1
      imjc=ij-1
      acj=co(9,ij)
      fj=(-(co(1,ij)*ph(ipjp)+co(3,ij)*ph(ipjm)+co(2,ij)*ph(ipjc)
     1    +co(5,ij)*ph(imjm)+co(7,ij)*ph(imjp)+co(6,ij)*ph(imjc)
     1    +co(4,ij)*ph(icjm)+co(8,ij)*ph(icjp))
     1    +res(ij))
      fej=fj/acj
      pho=ph(ij)
      ph(ij)=omeg*fej+(1.-omeg)*ph(ij)
      enddo
      enddo
      nit=nit+1
      if(nit.ge.nmaxi) go to 51
      go to 50
   51 continue
      return
      end
