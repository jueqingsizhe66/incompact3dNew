c
c   ********************* subr phcal
c  this subroutine perform the calculation of dph , periodic direction
c  along x3 allows to use the cfft99
c
      subroutine phcal(qcap,psi)
      include 'param.f'
      common/dim/n1,n1m,n2,n2m
      common/mesh/dx1,dx1q,dx2,dx2q
      common/ctrdii/ampi(m1),acpi(m1),appi(m1)
      common/ctrdjj/ampj(m2),acpj(m2),appj(m2)
      common/boupa/n1b1,n1bn1,n2b1,n2bn2
      common/boupn/np,mp
      common/ncfft/nufft
      dimension qcap(m1,m2),psi(m1,m2)
      if(np.eq.0) then
      an1m=n1m
      ax=alog(an1m)/alog(2.)
      k=ax+1
      L=2**(K+1)
      mw=(K-2)*L+K+5+2*n1+MAX(2*n1,6*n2)
                  else
      an1=n1
      ax=alog(an1)/alog(2.)
      k=ax+1
      L=2**(K+1) 
      mw=(K-2)*L+K+5+MAX(2*n1,6*n2)
                  endif
      npp=np
      mpp=mp
      if(mp.eq.0)  then
      nx=n1m
                   else
      nx=n1
                   endif
      if(np.eq.0)  then
      ny=n2m
                   else
      ny=n2
                   endif
      my=m1
c
c   here the entuality to multiply the 
c   matrices by Delta x1^2 Delta x2^2 is possible
c
c     den=1./(dx1q*dx2q)
      den=1.
      do j=1,n2
      ampj(j)=ampj(j)*den
      acpj(j)=acpj(j)*den
      appj(j)=appj(j)*den
      enddo
      do i=1,n1
      ampi(i)=ampi(i)*den
      acpi(i)=acpi(i)*den
      appi(i)=appi(i)*den
      do j=1,n2
      psi(i,j)=qcap(i,j)*den
      enddo
      enddo
      call sofish(npp,mpp,mw,my,psi,0,ny,nx)
      do n=1,nufft
      do i=1,n1
      do j=1,n2
      psi(i,j)=qcap(i,j)*den
      enddo
      enddo
      call sofish(npp,mpp,mw,my,psi,1,ny,nx)
      if(mp.eq.0.and.np.gt.0)  then
      do j=1,n2
      psi(n1,j)=psi(1,j)
      enddo
                   endif
      if(np.eq.0.and.mp.gt.0)  then
      do i=1,n1
      psi(i,n2)=psi(i,1)
      enddo
                   endif
      if(mp.eq.0.and.np.eq.0)  then
      do i=1,n1m
      psi(i,n2)=psi(i,1)
      enddo
      do j=1,n2
      psi(n1,j)=psi(1,j)
      enddo
                   endif
      enddo
      return
      end
c
c   ********************* subr sofish
c  this subroutine perform the calculation of dph , periodic direction
c  along x3 allows to use the cfft99
c
      subroutine sofish(np,mp,mw,idimy,y,iflg,n,m)
      include 'param.f'
      parameter (mwd=m1m*m2m)
      common/ctrdii/am(m1),bm(m1),cm(m1)
      common/ctrdjj/an(m2),bn(m2),cn(m2)
      dimension y(m1,m2),w(mwd)
c
c   this is the routine in the fishpack.f package
c
      call blktri(iflg,np,n,an,bn,cn,mp,m,am,bm,cm,idimy,y,ierror,w)
      if(iflg.eq.0) then
      do j=2,n
      prod=an(j)*cn(j-1)
      if(prod.le.0.) then
      write(6,133)j,an(j),bn(j),cn(j),prod
                     endif
      enddo
      write(6,*)'np,mp,ierror   ',np,mp,ierror
      i=1
      write(6,132)i,am(i),bm(i),cm(i)
      i=(m-1)/2+1
      write(6,132)i,am(i),bm(i),cm(i)
      i=m
      write(6,132)i,am(i),bm(i),cm(i)
  132 format(1x,'i=',i4,3x,4e12.5)
      j=1
      write(6,133)j,an(j),bn(j),cn(j)
      j=(n-1)/2+1
      write(6,133)j,an(j),bn(j),cn(j)
      j=n
      write(6,133)j,an(j),bn(j),cn(j)
  133 format(1x,'j=',i4,3x,4e12.5)
      write(6,*)'dimension w ',mw,w(1)
                    endif
      return
      end
c
c   ********************* subr phini
c  this subroutine calculate the l,u matrices to solve the poisson
c  eq. for dph this subroutine is called only once.
c  correspond to the operation ijob=1 in the leqt1b imsl routine.
c
      subroutine phini
      include 'param.f'
      common/dim/n1,n1m,n2,n2m
      common/mesh/dx1,dx1q,dx2,dx2q
      common/ctrdii/ampi(m1),acpi(m1),appi(m1)
      common/ctrdjj/ampj(m2),acpj(m2),appj(m2)
      common/boupa/n1b1,n1bn1,n2b1,n2bn2
      common/boupn/np,mp
c
      write(6,*)n1b1,n1bn1,n2b1,n2bn2
c
c  coefficients for the x1 direction
c
      do  i=1,n1
        ampi(i)=dx1q
        appi(i)=dx1q
        acpi(i)=-dx1q*2.
      enddo
      mp=0
c
c  coefficients for the x2 direction
c
      do  j=1,n2
        ampj(j)=dx2q
        appj(j)=dx2q
        acpj(j)=-dx2q*2.
      enddo
c
c   coefficients at the boundaries depending on
c   the chosen case
c
      !!!分为边界为1  和边界n的情况
      np=0
        i=1                          !!!!!!!!!i=1
        if(n1b1.eq.1) then
      mp=1
        ampi(i)=0.
        appi(i)=0.
        acpi(i)=1.
                      endif
        if(n1b1.eq.2) then
      mp=1
        ampi(i)=0.
        appi(i)=dx1
        acpi(i)=-dx1
                      endif
        i=n1                        !!!!!!!!!!!n
        if(n1bn1.eq.1) then
      mp=1
        ampi(i)=0.
        appi(i)=0.
        acpi(i)=1.
                      endif
        if(n1bn1.eq.2) then
      mp=1
        ampi(i)=-dx1
        appi(i)=0.
        acpi(i)=dx1
                      endif
        j=1                         !!!!!!!!!!!!j=1
        if(n2b1.eq.1) then
      np=1
        ampj(j)=0.
        appj(j)=0.
        acpj(j)=1.
                      endif
        if(n2b1.eq.2) then
      np=1
        ampj(j)=0.
        appj(j)=dx2
        acpj(j)=-dx2
                      endif
        j=n2                        !!!!!!!!!!!!! j=n2
        if(n2bn2.eq.1) then
      np=1
        ampj(j)=0.
        appj(j)=0.
        acpj(j)=1.
                      endif
        if(n2bn2.eq.2) then
      np=1
        ampj(j)=dx2
        appj(j)=0.
        acpj(j)=-dx2
                      endif
      return
      end
