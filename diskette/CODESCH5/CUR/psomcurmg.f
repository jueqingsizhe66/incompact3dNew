c
c   ********************* subr phini
c  this subroutine calculate the  number of computational
c  points in each level of the multigrid solver
c in matsb the metric quantities at each level are calculated
c
      subroutine phini
      include 'param.f'
      common/dim/n1,n1m,n2,n2m
      common/cnlev/n1g(0:nlev),n2g(0:nlev),n12g(0:nlev)
      common/mesh/dx1,dx1q,dx2,dx2q
      common/mlv/mlev,iwmg,mlw,maxcmg,epsm
      common/dak/ndm
      common/phkin/phk(mdm)
      igjg(i,j,l)=i+(j-1)*n1g(l)+n12g(l-1)
      ndm=n1*n2
      n1g(0)=0
      n2g(0)=0
      n12g(0)=0
      do nl=1,mlev
      n1g(nl)=n1m*2**float(-nl+1)+1
      n2g(nl)=n2m*2**float(-nl+1)+1
      n12g(nl)=n1g(nl)*n2g(nl)+n12g(nl-1)
      write(6,333)nl,n1g(nl),n2g(nl),n12g(nl)
  333 format(3x,'nl=',i3,3(2x,i5))
      enddo
      call matsb
      do j=1,n2
      do i=1,n1
      ij=igjg(i,j,1)
      phk(ij)=0.
      enddo
      enddo
      return
      end
c
c   ********************* subr dsolv
c    this is the principal routine for the
c    mutigrid solver
c
      subroutine dsolv(qk,phk)
      include 'param.f'
      dimension qk(mdm),phk(mdm),phm(nij),res(nij),rest(nij)
      common/cnlev/n1g(0:nlev),n2g(0:nlev),n12g(0:nlev)
      common/dak/ndm
      common/mlv/mlev,iwmg,mlw,maxcmg,epsm
      common/mgoui/ncm
      common/mgout/reml1,reml2,rrm
      common/psbou/psbi(2,m1),psbj(2,m2)
      igjg(i,j,l)=i+(j-1)*n1g(l)+n12g(l-1)
c
c the finel level
c
      ic=1
      ip=n1g(1)
      do jc=1,n2g(1)
      phm(igjg(ic,jc,1))=psbj(1,jc)
      phm(igjg(ip,jc,1))=psbj(2,jc)
      enddo
      jc=1
      jp=n2g(1)
      do ic=1,n1g(1)
      phm(igjg(ic,jc,1))=psbi(1,ic)
      phm(igjg(ic,jp,1))=psbi(2,ic)
      enddo
      do ic=2,n1g(1)-1
      do jc=2,n2g(1)-1
      ij=igjg(ic,jc,1)
      phm(ij)=phk(ij)
      rest(ij)=qk(ij)
      enddo
      enddo
      ncm=0
      call resma(1,phm,rest)
  338 format(3x,'nl=1',2x,2e12.5,2x,2i4)
c     write(6,338) reml1
      reool1=reml1
    9 continue
c
c   cfila routine  where the levels from fine to coarse
c   are performed
c
      call cfila(1,mlev,phm,rest,res)
      if(iwmg.gt.0) then
c
c   W cicle with iwmg the number of
c   W cicles
c
      do 30 mw=1,iwmg
      mlwl=mlev-mw
      if(iwmg.eq.1) mlwl=mlw
c
c   clafi routine  where the levels from coarse
c   to fine are performed
c
      call clafi(mlwl,mlev,phm,rest)
      call cfila(mlwl,mlev,phm,rest,res)
   30 continue
                     endif
      call clafi(1,mlev,phm,rest)
      call resma(1,phm,rest)
      ncm=ncm+1
c     write(6,339) ncm,reml1
  339 format(3x,'ncm=',i4,2x,2e12.5)
      alrem1=log(reml1)
      if(ncm.gt.maxcmg) go to 12
      if(reml1.gt.epsm) then
      go to 9
      endif
   12 continue
      aml1=log(reml1/reool1)/ncm
      rrm=exp(-aml1)
      do ic=1,n1g(1)
      do jc=1,n2g(1)
      ij=igjg(ic,jc,1)
      phk(ij)=phm(ij)
      enddo
      enddo
      return
      end
c
c   ********************* subr cfila
c
c  From fine to coarse cicle
c  it consist in a neumber of iterative solvers by (itsolv) 
c   followed by the evaluation of the residual calculation
c   in resca
c   when the coarsest level is reached the residual is 
c   transfered to a fine grid.
c
      subroutine cfila(mli,mlev,phm,rest,res)
      include 'param.f'
      dimension phm(nij),res(nij),rest(nij)
      common/mgout/reml1,reml2,rrm
      common/cnlev/n1g(0:nlev),n2g(0:nlev),n12g(0:nlev)
      igjg(i,j,l)=i+(j-1)*n1g(l)+n12g(l-1)
      do 20 nl=mli,mlev
      call itsolv(nl,phm,rest)
      call resca(nl,phm,res,rest)
      do ic=1,n1g(nl+1)
      do jc=1,n2g(nl+1)
      ij=igjg(ic,jc,nl+1)
      phm(ij)=0.
      enddo
      enddo
      if(nl.lt.mlev) then
      call tfila(res,nl,rest)
      endif
   20 continue
      return
      end
c
c
c   ********************* subr clafi
c
c  from coarse to fine
c   iterative solver at the corsest followed
c   by the tlafi where the  solution psi is 
c   transfered by a coarse to a fine grid.
c   this solution is the initial solution of 
c   the iterative solver
c
      subroutine clafi(mli,mlev,phm,rest)
      include 'param.f'
      common/mgout/reml1,reml2,rrm
      dimension phm(nij),rest(nij)
      call itsolv(mlev,phm,rest)
      do 30 nl=mlev,mli+1,-1
      call tlafi(phm,nl)
      call itsolv(nl-1,phm,rest)
   30 continue
      return
      end
c
c   ********************* subr phcalc
c  this subroutine perform the calculation of psi by a multigrid
c  solver the vorticity field is tyransfered into a one dimensional
c  array qk( and the psi returns in the one-dimensional array phk(
c
      subroutine phcalc(qcap,dph)
      include 'param.f'
      common/dim/n1,n1m,n2,n2m
      dimension qk(mdm),qcap(m1,m2),dph(m1,m2)
      common/phkin/phk(mdm)
      common/mesh/dx1,dx1q,dx2,dx2q
      common/dak/ndm
      common/cnlev/n1g(0:nlev),n2g(0:nlev),n12g(0:nlev)
      igjg(i,j,l)=i+(j-1)*n1g(l)+n12g(l-1)
      do j=2,n2m
      do i=2,n1m
      qk(igjg(i,j,1))=-qcap(i,j)
      phk(igjg(i,j,1))=dph(i,j)
      enddo
      enddo
c
c    the multigrid is in this routine
c
      call dsolv(qk,phk)
      do j=1,n2
      do i=1,n1
      dph(i,j)=phk(igjg(i,j,1))
      enddo
      enddo
      return
      end
c
c  ****************************** subrout resca **********************
c   residual calculation of the tentative solution
c
      subroutine resca(nl,ph,rs,rst)
      include 'param.f'
      common/cnlev/n1g(0:nlev),n2g(0:nlev),n12g(0:nlev)
      common/comg/co(10,nij)
      dimension ph(nij),rs(nij),rst(nij)
      igjg(i,j,l)=i+(j-1)*n1g(l)+n12g(l-1)
      do jc=2,n2g(nl)-1
      do ic=2,n1g(nl)-1
      ij=igjg(ic,jc,nl)
      d2phij=+(co(1,ij)*ph(igjg(ic+1,jc+1,nl))+
     1         co(2,ij)*ph(igjg(ic+1,jc,nl))+
     1         co(3,ij)*ph(igjg(ic+1,jc-1,nl))+
     1         co(4,ij)*ph(igjg(ic,jc-1,nl))+
     1         co(5,ij)*ph(igjg(ic-1,jc-1,nl))+
     1         co(6,ij)*ph(igjg(ic-1,jc,nl))+
     1         co(7,ij)*ph(igjg(ic-1,jc+1,nl))+
     1         co(8,ij)*ph(igjg(ic,jc+1,nl))+
     1         co(9,ij)*ph(ij))
      rs(ij)=(-d2phij+rst(ij))
      enddo
      enddo
      return
      end
c
c  ****************************** subrout resma **********************
c
c  calculation of the maximum residual in the whole computational box
c
      subroutine resma(nl,ph,res)
      include 'param.f'
      common/cnlev/n1g(0:nlev),n2g(0:nlev),n12g(0:nlev)
      common/comg/co(10,nij)
      dimension ph(nij),res(nij)
      common/mgout/reml1,reml2,rrm
      igjg(i,j,l)=i+(j-1)*n1g(l)+n12g(l-1)
      reml1=0.
      do jc=2,n2g(nl)-1
      do ic=2,n1g(nl)-1
      ij=igjg(ic,jc,nl)
      d2ph=+(co(1,ij)*ph(igjg(ic+1,jc+1,nl))+
     1       co(2,ij)*ph(igjg(ic+1,jc,nl))+
     1       co(3,ij)*ph(igjg(ic+1,jc-1,nl))+
     1       co(4,ij)*ph(igjg(ic,jc-1,nl))+
     1       co(5,ij)*ph(igjg(ic-1,jc-1,nl))+
     1       co(6,ij)*ph(igjg(ic-1,jc,nl))+
     1       co(7,ij)*ph(igjg(ic-1,jc+1,nl))+
     1       co(8,ij)*ph(igjg(ic,jc+1,nl))+
     1       co(9,ij)*ph(ij))
      rsij=(-d2ph+res(ij))
      reml1=max(abs(rsij),reml1)
      enddo
      enddo
      return
      end
c
c  ****************************** subrout tfila **********************
c
c  the residual in a fine grid is transfered in a coarser grid
c
      subroutine tfila(res,nl,rest)
      include 'param.f'
      common/cnlev/n1g(0:nlev),n2g(0:nlev),n12g(0:nlev)
      common/comg/co(10,nij)
      dimension res(nij),rest(nij)
      igjg(i,j,l)=i+(j-1)*n1g(l)+n12g(l-1)
      do jcl=2,n2g(nl+1)-1
      do icl=2,n1g(nl+1)-1
      ijl=igjg(icl,jcl,nl+1)
      icjcf=igjg(2*icl-1,2*jcl-1,nl)
      ipjcf=icjcf+1
      imjcf=icjcf-1
      icjpf=igjg(2*icl-1,2*jcl,nl)
      icjmf=igjg(2*icl-1,2*jcl-2,nl)
      ipjpf=icjpf+1
      imjmf=icjmf-1
      imjpf=icjpf-1
      ipjmf=icjmf+1
      rest(ijl)=((res(ipjcf)*co(10,ipjcf)+res(imjcf)*co(10,imjcf)
     1           +res(icjpf)*co(10,icjpf)+res(icjmf)*co(10,icjmf))*2.
     1           +res(ipjpf)*co(10,ipjpf)+res(imjpf)*co(10,imjpf)
     1           +res(imjmf)*co(10,imjmf)+res(ipjmf)*co(10,ipjmf)
     1        +4.*res(icjcf)*co(10,icjcf))/(16.*co(10,ijl))
      enddo
      enddo
      return
      end
c
c  ****************************** subrout tlafi **********************
c
c  the solution in a coarse grid is transfered in a finer grid
c
      subroutine tlafi(ph,nl)
      include 'param.f'
      common/cnlev/n1g(0:nlev),n2g(0:nlev),n12g(0:nlev)
      dimension ph(nij)
      igjg(i,j,l)=i+(j-1)*n1g(l)+n12g(l-1)
      do jc=2,n2g(nl)-1,2
      do ic=2,n1g(nl)-1,2
      icjcl=igjg(ic,jc,nl)
      ipjcl=icjcl+1
      imjcl=icjcl-1
      icjpl=igjg(ic,jc+1,nl)
      ipjpl=icjpl+1
      imjpl=icjpl-1
      icjml=igjg(ic,jc-1,nl)
      ipjml=icjpl+1
      imjml=icjpl-1
      icjcf=igjg(2*ic-1,2*jc-1,nl-1)
      imjcf=icjcf-1
      ipjcf=icjcf+1
      icjmf=igjg(2*ic-1,2*jc-2,nl-1)
      imjmf=icjmf-1
      ipjmf=icjmf+1
      icjpf=igjg(2*ic-1,2*jc,nl-1)
      imjpf=icjpf-1
      ipjpf=icjpf+1
      ph(icjcf)=ph(icjcl)+ph(icjcf)
      ph(icjpf)=(ph(icjpl)+ph(icjcl))*0.5+ph(icjpf)
      ph(ipjpf)=(ph(ipjpl)+ph(ipjcl)+ph(icjcl)+ph(icjpl))*0.25+ph(ipjpf)
      ph(ipjcf)=(ph(ipjcl)+ph(icjcl))*0.5+ph(ipjcf)
      ph(ipjmf)=(ph(ipjcl)+ph(ipjml)+ph(icjml)+ph(icjcl))*0.25+ph(ipjmf)
      ph(icjmf)=(ph(icjcl)+ph(icjml))*0.5+ph(icjmf)
      ph(imjmf)=(ph(icjcl)+ph(icjml)+ph(imjml)+ph(imjcl))*0.25+ph(imjmf)
      ph(imjcf)=(ph(icjcl)+ph(imjcl))*0.5+ph(imjcf)
      ph(imjpf)=(ph(icjpl)+ph(icjcl)+ph(imjcl)+ph(imjpl))*0.25+ph(imjpf)
      enddo
      enddo
      return
      end
c
c  ****************************** subrout matsb  **********************
c
c   in this subr the coefficients of the poisson eq. at each
c   grid level  are calculated. this subr. is called only 
c   once at the beginning of the simulation
c
      subroutine matsb
      include 'param.f'
      common/cnlev/n1g(0:nlev),n2g(0:nlev),n12g(0:nlev)
      common/metrip/alfi(ndd,ndd,m1,m2),alfj(ndd,ndd,m1,m2)
      common/metriv/alfc(ndd,ndd,m1,m2)
      common/metrst/gccc(m1,m2)
      common/comg/co(10,nij)
      common/mlv/mlev,iwmg,mlw,maxcmg,epsm
      igjg(i,j,l)=i+(j-1)*n1g(l)+n12g(l-1)
      do 100 nl=1,mlev
      dx1q=float((n1g(nl)-1)**2)
      dx2q=float((n2g(nl)-1)**2)
      dxcd=float(n2g(nl)-1)*float(n1g(nl)-1)/4.
      write(6,335)nl,n1g(nl),n2g(nl),dx1q,dx2q,dxcd
  335 format(3x,'nl=',i3,2x,2(2x,i3),2x,'mesh',3e12.4)
      do jc=2,n2g(nl)-1
      do ic=2,n1g(nl)-1
c
c   coefficients for the nine point stencil of the poisson eq.
c   1=(ip,jp);2=(ip,jc);3=(ip,jm);4=(ic,jm);5=(im,jm);6=(im,jc)
c   7=(im,jp);8=(ic,jp);9=(ic,jc)
      ip=ic+1
      im=ic-1
      jm=jc-1
      jp=jc+1
      jci=(jc-1)*2**(nl-1)+1
      ici=ic*2**(nl-1)-(nl-2)
      imi=(ic-1)*2**(nl-1)-(nl-2)
      icj=(ic-1)*2**(nl-1)+1
      jcj=jc*2**(nl-1)-(nl-2)
      jmj=(jc-1)*2**(nl-1)-(nl-2)
      if(nl.eq.1) a22jcc=alfj(2,2,ic,jc)
      if(nl.gt.1) a22jcc=alfc(2,2,icj,jcj)
      if(nl.eq.1) a21jcc=alfj(2,1,ic,jc)
      if(nl.gt.1) a21jcc=alfc(2,1,icj,jcj)
      if(nl.eq.1) a22jcm=alfj(2,2,ic,jm)
      if(nl.gt.1) a22jcm=alfc(2,2,icj,jmj)
      if(nl.eq.1) a21jcm=alfj(2,1,ic,jm)
      if(nl.gt.1) a21jcm=alfc(2,1,icj,jmj)
      if(nl.eq.1) a11icc=alfi(1,1,ic,jc)
      if(nl.gt.1) a11icc=alfc(1,1,ici,jci)
      if(nl.eq.1) a12icc=alfi(1,2,ic,jc)
      if(nl.gt.1) a12icc=alfc(1,2,ici,jci)
      if(nl.eq.1) a11imc=alfi(1,1,im,jc)
      if(nl.gt.1) a11imc=alfc(1,1,imi,jci)
      if(nl.eq.1) a12imc=alfi(1,2,im,jc)
      if(nl.gt.1) a12imc=alfc(1,2,imi,jci)
      if(nl.eq.1) ugccc=1./gccc(ic,jc)
      if(nl.gt.1) ugccc=1./gccc(icj,jci)
      ij=igjg(ic,jc,nl)
      co(1,ij)=(a12icc+a21jcc)*dxcd*ugccc
      co(2,ij)=(a11icc*dx1q+(+a21jcc-a21jcm)*dxcd)*ugccc
      co(3,ij)=(-a12icc-a21jcm)*dxcd*ugccc
      co(4,ij)=(a22jcm*dx2q+(-a12icc+a12imc)*dxcd)*ugccc
      co(5,ij)=(+a12imc+a21jcm)*dxcd*ugccc
      co(6,ij)=(a11imc*dx1q+(-a21jcc+a21jcm)*dxcd)*ugccc
      co(7,ij)=(-a12imc-a21jcc)*dxcd*ugccc
      co(8,ij)=(a22jcc*dx2q+(+a12icc-a12imc)*dxcd)*ugccc
      co(9,ij)=-((a11icc+a11imc)*dx1q+(a22jcc+a22jcm)*dx2q)*ugccc
      co(10,ij)=1./ugccc
      enddo
      enddo
  100 continue
      close(37)
      return
      end
