c
c  ****************************** subrout meshes **********************
c
      subroutine meshes
      include 'param.f'
      common/dim/n1,n1m,n2,n2m
      common/mesh/dx1,dx1q,dx2,dx2q
      common/d1/alx1i,alx1f,alx2i,alx2f
      dx1=1./float(n1m)
      dx2=(alx2f-alx2i)/float(n2m)
      dx1=1./dx1
      dx2=1./dx2
      dx1q=dx1*dx1
      dx2q=dx2*dx2
      return
      end
c
c  ****************************** subrout indic **********************
c
c  in this subroutine the indices ip,im and jp jm are calculated
c  these are used in the evaluation of the derivatives to
c  avoid if statements for the boundary points
c 
      subroutine indic
      include 'param.f'
      common/dim/n1,n1m,n2,n2m
      common/indbo/imv(m1),ipv(m1)
      common/indwa/jmv(m2),jpv(m2)
      common/iperi/ib2per
c
c   x1 radiative  direction
c
      do 11 ic=1,n1
      ipv(ic)=ic+1
      if(ic.eq.n1) ipv(ic)=ic
      imv(ic)=ic-1
      if(ic.eq.1) imv(ic)=ic
   11 continue
      if(ib2per.eq.1) then
c
c   periodic directions
c
      do  jc=1,n2m
      jpv(jc)=jc+1
      if(jc.eq.n2m) jpv(jc)=1
      jmv(jc)=jc-1
      if(jc.eq.1) jmv(jc)=n2m
      enddo
                       else
c
c   free-slip and no-slip walls directions
c
      do jc=1,n2m
      jpv(jc)=jc+1
      jmv(jc)=jc-1
      if(jc.eq.1) jmv(jc)=jc
      enddo
                       endif
      return
      end
c
c
c  ****************************** subrout coordi  **********************
c in this routine the coordinate transformations
c are used these are described in the Sect.2.2.3
c
      subroutine coordi
      include 'param.f'
      common/dim/n1,n1m,n2,n2m
      common/mesh/dx1,dx1q,dx2,dx2q
      common/coor/yp1(m1),yp2(m2)
      common/coorh/ym1(m1)
      common/d1/alx1i,alx1f,alx2i,alx2f
      common/ydif/ydi(m2)
      dimension eta(m1),xtp(m1),xts(m1)
      common/metri/g2m(m1),g2c(m1)
      common/strpar/str1,xcra,etra,istr
      common/cwalv/inbcvs,inbcvn,inbcvw,inbcve
      if(istr.eq.1) then
      nc=n1m*etra
      xc=alx1f*xcra
      xf=alx1f
      alin=str1
      csic=(nc-1)/float(n1m)*0.5+0.5
      etac=abs(csic-0.5)
      etaf=0.5
      zitaf=etaf-etac
      al=alin
      beta=xc/etac
      iter=1
      val=beta/(xf-xc)
c     write(6,*)val,zitaf
   33 continue
      res=gal(al,zitaf)-val
      if(abs(res).lt..1e-04) then
      go to 35
                             else
      tang=dgal(al,zitaf)
      dal=res/tang
      write(6,*)iter,res,al,tang,dal
      al=al-dal
      iter=iter+1
      if(iter.gt.21) go to 35
      go to 33
                             endif
   35 continue
      write(6,*)'convergence',iter,res,al
      gam=tanh(al*zitaf)                                                      
      do i=1,n1                                                      
      csi=(i-1)/float(n1m)                                   
      etal=csi-0.5
      if(abs(etal).le.etac)  then
      eta(i)=etal*xc/etac
                            else 
      siet=etal/abs(etal) 
      if(siet.ge.0) then
      zita=(etal-etac)
      azit=al*(zita-zitaf)
      eta(i)=(xf+tanh(azit)/gam*(xf-xc))
                            else 
      zita=(etal+etac)
      azit=al*(zita+zitaf)
      eta(i)=(-xf+tanh(azit)/gam*(xf-xc))
                            endif                                   
                            endif                                   
c     write(18,110)csi,eta(i)
  110 format(3x,2e12.5)
      enddo
      pi=2.*asin(1.)
      do 65 i=1,n1
      yp1(i)=eta(i)
c     write(68,*)i,yp1(i)
   65 continue
                       endif
      if(istr.eq.-2) then
      al=str1
      gam=tanh(al*0.5)
      do i=1,n1
      csi=(i-1)/float(n1m)
      etal=al*(csi-0.5)
      xt1=0.5*(1.+tanh(etal)/gam)
      yp1(i)=(-0.5+xt1)*2.*alx1f
      enddo
                     endif
      if(istr.eq.-1) then
      al=str1
      gam=tanh(al)
      do i=1,n1
      csi=(i-1)/float(n1m)
      etal=al*(csi-1.)
      xt1=(1.+tanh(etal)/gam)
      yp1(i)=xt1*alx1f
      enddo
                     endif
      if(istr.eq.0) then
      do i=1,n1                                                      
      csi=(i-1)/float(n1m)                                   
      yp1(i)=alx1i+csi*(alx1f-alx1i)
      enddo
                       endif
      if(istr.eq.2) then
      al=str1
      gam=tanh(al)
      do i=1,n1
      csi=(i-1)/float(n1m)
      etal=al*csi
      xt1=tanh(etal)/gam
      yp1(i)=alx1i+xt1*(alx1f-alx1i)
      enddo
                       endif
      if(istr.eq.3) then
      nc=n1m*etra
      xc=alx1f*xcra
      xf=alx1f
      al1=str1
      al2=al1+1.
      csic=(nc-1)/float(n1m)
      gam=tanh(al1*csic)                                                      
      do i=1,n1
      csi=(i-1)/float(n1m)
      xtp(i)=xc/xf*tanh(al1*csi)/gam
      enddo
      do i=1,n1
      csi=(i-1)/float(n1m)
      xts(i)=1./xtp(n1)+(1.-1./xtp(n1))
     1     *tanh(al2*(csi-1.))/tanh(al2*(csic-1.))
      yp1(i)=xtp(i)*xts(i)*xf
      enddo
                           endif
      do i=1,n1m
      ym1(i)=(yp1(i+1)+yp1(i))*0.5
      enddo
      do 66 j=1,n2
      x2=(j-1)/float(n2m)
      yp2(j)=alx2i+x2*(alx2f-alx2i)
      ydi(j)=x2
   66 continue
      do i=1,n1m
      g2m(i)=(yp1(i+1)-yp1(i))*dx1
      enddo
      open(79,file='strcor.out')
      do i=2,n1m
      g2c(i)=(yp1(i+1)-yp1(i-1))*dx1*0.5
      enddo
      g2c(1)=(yp1(2)-yp1(1))*dx1
      g2c(n1)= (yp1(n1)-yp1(n1m))*dx1
c thus alx1i<yp1<alx1f
c and  alx2i<yp2<alx2f
c next just put yp1 and yp2 in an unformatted file
      do i=1,n1
      csi=(i-1)/float(n1m)                                   
      write(79,133)csi,yp1(i)
  133 format(3x,5e12.5)
      enddo
      do i=1,n1m
      csi=(i-1)/float(n1m)                                   
      write(89,133)csi,yp1(i),g2c(i),g2m(i)
      enddo
      close(79)
      call pricor(yp1,yp2)
      return
      end
c
c
c  ****************************** subrout tribk  **********************
c  this routine inverts a tridiagonal matrix
c
      subroutine tribk(a,b,c,f,n,mi,mf)
      include 'param.f'
      dimension gam(md2),u(md2,md2)
      dimension c(md2),b(md2),a(md2),f(md2,md2)
      bet=b(1)
      do k=mi,mf
      u(k,1)=f(k,1)/bet
      do j=2,n
      gam(j)=c(j-1)/bet
      bet=b(j)-a(j)*gam(j)
      u(k,j)=(f(k,j)-a(j)*u(k,j-1))/bet
      enddo
      enddo
      do j=n-1,1,-1
      do k=mi,mf
      u(k,j)=u(k,j)-gam(j+1)*u(k,j+1)
      enddo
      enddo
      do j=1,n
      do k=mi,mf
      f(k,j)=u(k,j)
      enddo
      enddo
      return
      end
c
c  ****************************** subrout tripv  **********************
c this routine inverts atridiagonal matrix
c with periodic conditions
c
      subroutine tripv( a,b,c,f,j2,mi,mf )
      include 'param.f'
c
      dimension a(md2),b(md2),c(md2),f(md2,md2)
     1         ,q(md2),s(md2),qe(md2,md2),fn(md2)
      j1=1
      ja = j1 + 1
      jj = j1 + j2
      q(j1) = -c(j1)/b(j1)
      s(j1) = - a(j1)/b(j1)
      do 20 k=mi,mf
      fn(k) = f(k,j2)
      f(k,j1) = f(k,j1)/b(j1)
   20 continue
c
c     forward elimination sweep
c
      do 10 j=ja,j2
      p =1./( b(j) + a(j)*q(j-1))
      q(j) = - c(j)*p
      s(j) = - a(j)*s(j-1)*p
      do 21 k=mi,mf
      f(k,j) = ( f(k,j) - a(j)*f(k,j-1))*p
   21 continue
   10 continue
c
c     backward pass
c
      s(j2) = 1.
      do 22 k=mi,mf
      qe(k,j2) = 0.
   22 continue
      do 11 i=ja,j2
      j = jj - i
      s(j) = s(j) + q(j)*s(j+1)
      do 23 k=mi,mf
      qe(k,j) = f(k,j) + q(j)*qe(k,j+1)
   23 continue
   11 continue
      do 24 k=mi,mf
      f(k,j2)=(fn(k) - c(j2)*qe(k,j1) - a(j2)*qe(k,j2-1))
     &       /(c(j2)*s(j1) + a(j2)*s(j2-1)  +b(j2))
   24 continue
c
c     backward elimination pass
c
      do 12 i=ja,j2
      j = jj -i
      do 25 k=mi,mf
      f(k,j) = f(k,j2)*s(j) + qe(k,j)
   25 continue
   12 continue
      return
      end
c  **************************************************
      subroutine minmax(w,nm1,nm2,imn,jmn,imx,jmx)
      include 'param.f'
      dimension w(m1,m2)
      wmn=w(1,1)
      imn=1
      jmn=1
      wmx=w(1,1)
      imx=1
      jmx=1
      do 200 j=1,nm2
      do 200 i=1,nm1
      wt=w(i,j)
      testn=wmn-wt
      if(testn.le.0.)go to 20
      imn=i
      jmn=j
      wmn=wt
20    continue
      testx=wmx-wt
      if(testx.ge.0.)go to 200
      imx=i
      jmx=j
      wmx=wt
200   continue
      write(6,*)'  in minmax min = ',wmn,' max = ',wmx
      return
      end
c
c  ****************************** subrout cfl  **********************
c
c  in this subroutine is calculated the maximum courant number
c
      subroutine cfl(psi,cflm)
      include 'param.f'
      common/dim/n1,n1m,n2,n2m
      common/mesh/dx1,dx1q,dx2,dx2q
      common/metri/g2m(m1),g2c(m1)
      dimension psi(m1,m2)
      common/indbo/imv(m1),ipv(m1)
      common/indpe2/n2i,n2f
      common/indwa/jmv(m2),jpv(m2)
c
      cflm=0.
      udx12=dx2*dx1*0.5
      do 10 ic=2,n1m   
      ip=ipv(ic)
      im=imv(ic)
      ud12=udx12/g2c(ic)
      do 10 jc=n2i,n2m
      jm=jmv(jc)
      jp=jpv(jc)
      h22a=((psi(ic,jp)-psi(ic,jm))-
     1      (psi(ip,jc)-psi(im,jc)))*ud12
      cflm=max(abs(h22a),cflm)
   10 continue
      return
      end
c
c  ****************************** subrout inirea  **********************
c   in this routine the initail fields are read,
c   these fields were written in a previous run
c   in the routine  OUTPF
c
      subroutine inirea(vor,psi,timei,vorq,psc,ru,rt)
      include 'param.f'
      character*22 namfile
      character*4 ipre
      character*3 ipfi
      dimension vor(m1,m2),psi(m1,m2),vorq(m1,m2)
      dimension psc(mpsc,m1,m2),ru(m1,m2),rt(mpsc,m1,m2)
      common/coor/yp1(m1),yp2(m2)
      common/mesh/dx1,dx1q,dx2,dx2q
      common/dim/n1,n1m,n2,n2m
      common/inmodo/inmod
      common/velcho/choe(2+mpsc,m2),chow(2+mpsc,m2)
      common/velcve/cven(2+mpsc,m1),cves(2+mpsc,m1)
      common/dvpsol/dvesto(m2),dpesto(m2),dvweso(m2),dpweso(m2)
     1         ,dsesto(mpsc,m2),dsweso(mpsc,m2)
      common/dvpsvl/dvnoro(m1),dpnoro(m1),dvsudo(m1),dpsudo(m1)
     1         ,dsnoro(mpsc,m1),dssudo(mpsc,m1)
      common/pscqu/sch(mpsc),scla,npscf
      common/iperi/ib2per
      nfilp=nint(timei)
      write(ipfi,96)nfilp
      namfile='tobe'//ipfi//'.dat'
   96 format(i3.3)
c
c   write potential vort,bigpsi(flow u0)  psi and relative vorticity
c
      nfilo=19
      write(6,202)timei,namfile
  202 format(1x,'a t=',e10.3,3x,3x,
     1       'read quantities on file=',a22)
      open(nfilo,file=namfile,form='unformatted')
      rmrc=0.
      read(nfilo)n1j,n2j,nju
      read(nfilo)aa,aa,aa,timei
      if(inmod.ne.10) then
      if(npscf.eq.2) then
      read(nfilo) ((psi(i,j),i=1,n1),j=1,n2)
     1            ,((psc(2,i,j),i=1,n1),j=1,n2)
     1            ,((psc(1,i,j),i=1,n1),j=1,n2)
     1            ,((vor(i,j),i=1,n1),j=1,n2)
                     endif
      if(npscf.eq.1) then
      read(nfilo) ((psi(i,j),i=1,n1),j=1,n2)
     1            ,((vorq(i,j),i=1,n1),j=1,n2)
     1            ,((psc(1,i,j),i=1,n1),j=1,n2)
     1            ,((vor(i,j),i=1,n1),j=1,n2)
                     endif
                     else
      read(nfilo) ((psi(i,j),i=1,n1),j=1,n2)
     1            ,((vorq(i,j),i=1,n1),j=1,n2)
     1            ,((psc(1,i,j),i=1,n1),j=1,n2)
     1            ,((vor(i,j),i=1,n1),j=1,n2)
     1            ,((fo,i=1,n1),j=1,n2)
     1            ,((vc,i=1,n1),j=1,n2)
     1            ,((vc,i=1,n1),j=1,n2)
                     endif
      close(nfilo)
      namfile='bcpsco'//ipfi//'.dat'
      open(nfilo,file=namfile,form='unformatted')
      read(nfilo)choe,chow,cves,cven,
     1          dvesto,dpesto,dvweso,dpweso
     1         ,dsesto,dsweso
     1         ,dvnoro,dpnoro,dvsudo,dpsudo
     1         ,dsnoro,dssudo
      close(nfilo)
      do i=1,n1
      do j=1,n2
      ru(i,j)=0.
      vorq(i,j)=0.
      do npsc=1,npscf
      rt(npsc,i,j)=0.
      enddo
      enddo
      enddo
      return
      end
c
c  ****************************** subrout outpf  **********************
c  in this routine the fields are wriiten for visualizations
c  by TUR3D . The same field can be used as restarting file
c  to be read in the routine inirea
c
      subroutine outpf(vor,psi,time,vorq,psc)
      include 'param.f'
      character*22 namfile
      character*4 ipre
      character*3 ipfi
      dimension vor(m1,m2),psi(m1,m2),vorq(m1,m2)
      dimension vca1(m1,m2),vca2(m1,m2)
      dimension psc(mpsc,m1,m2)
      common/forci/forc(m1,m2),fff
      common/inmodo/inmod
      common/coor/yp1(m1),yp2(m2)
      common/mesh/dx1,dx1q,dx2,dx2q
      common/metri/g2m(m1),g2c(m1)
      common/dim/n1,n1m,n2,n2m
      common/visct/re,ekmn
      common/nmolt/nmolp
      common/bcvedp/dpest(m2),dpwes(m2)
      common/bvedph/dphest(m2),dphwes(m2)
      common/bchodp/dpsud(m1),dpnor(m1)
      common/betef/beta
      common/ydif/ydi(m2)
      common/velcho/choe(2+mpsc,m2),chow(2+mpsc,m2)
      common/velcve/cven(2+mpsc,m1),cves(2+mpsc,m1)
      common/dvpsol/dvesto(m2),dpesto(m2),dvweso(m2),dpweso(m2)
     1         ,dsesto(mpsc,m2),dsweso(mpsc,m2)
      common/dvpsvl/dvnoro(m1),dpnoro(m1),dvsudo(m1),dpsudo(m1)
     1         ,dsnoro(mpsc,m1),dssudo(mpsc,m1)
      common/bcvev/vwes(m2),vest(m2)
      common/bcvep/pwes(m2),pest(m2)
      common/bcves/swes(mpsc,m2),sest(mpsc,m2)
      common/bchov/vsud(m1),vnor(m1)
      common/bchop/psud(m1),pnor(m1)
      common/bchos/ssud(mpsc,m1),snor(mpsc,m1)
      common/angmd/thet0
      common/pscqu/sch(mpsc),scla,npscf
      common/topog/htop(m1,m2)
      common/quainf/qinp(6,md2)
      common/cwalv/inbcvs,inbcvn,inbcvw,inbcve
      dimension pwal(m2)
      common/iperi/ib2per
      if(ib2per.eq.1) then
      do i=1,n1
      vor(i,n2)=vor(i,1)
      psi(i,n2)=psi(i,1)
      enddo
      do l=1,npscf
      do i=1,n1
      psc(l,i,n2)=psc(l,i,1)
      enddo
      enddo
                      endif
c
c   the Cartesian velocities are evaluated
c
      do i=1,n1
      do j=2,n2m
      vca1(i,j)=(psi(i,j+1)-psi(i,j-1))*dx2*0.5
      enddo
      enddo
      do i=1,n1
      vca1(i,n2)=(psi(i,n2)-psi(i,n2m))*dx2
      vca1(i,1)=(psi(i,2)-psi(i,1))*dx2
      enddo
      do j=1,n2
      vca2(1,j)=-(psi(2,j)-psi(1,j))/(yp1(2)-yp1(1))
      vca2(n1,j)=-(psi(n1,j)-psi(n1m,j))/(yp1(n1)-yp1(n1m))
      enddo
      do i=2,n1m
      do j=1,n2
      vca2(i,j)=-(psi(i+1,j)-psi(i-1,j))/(yp1(i+1)-yp1(i-1))
      enddo
      enddo
c
c   potential vorticity
c
      do i=1,n1
      do j=1,n2
      vorq(i,j)=vor(i,j)+beta*yp2(j)+htop(i,j)
      enddo
      enddo
      psmal=0.
      vomal=0.
      vqmal=0.
      do i=1,n1
      do j=1,n2
      psmal=max(abs(psi(i,j)),psmal)
      vomal=max(abs(vor(i,j)),vomal)
      vqmal=max(abs(vorq(i,j)),vqmal)
      enddo
      enddo
      write(6,137)psmal,vomal,vqmal 
  137 format(3x,'  max quant  ***** in outp   ',3e12.4,2x)
      nfilp=nint(time*nmolp)
      write(ipfi,96)nfilp
      namfile='tobe'//ipfi//'.dat'
   96 format(i3.3)
c
c   write potential vort,bigpsi(flow u0)  psi and relative vorticity
c
      nfilo=19
  202 format(1x,'a t=',e10.3,3x,3x,
     1       'write quantities on file=',a22)
      open(nfilo,file=namfile,form='unformatted')
      nju=1
      rmrc=0.
      write(nfilo)n1,n2,nju
      write(nfilo)thet0,beta,re,time
      if(inmod.ne.10) then
      if(npscf.eq.2) then
      write(6,202)time,namfile
      write(nfilo) ((psi(i,j),i=1,n1),j=1,n2)
     1            ,((psc(2,i,j),i=1,n1),j=1,n2)
     1            ,((psc(1,i,j),i=1,n1),j=1,n2)
     1            ,((vor(i,j),i=1,n1),j=1,n2)
     1            ,((vca1(i,j),i=1,n1),j=1,n2)
     1            ,((vca2(i,j),i=1,n1),j=1,n2)
                     endif
      if(npscf.eq.1) then
      write(6,202)time,namfile
      write(nfilo) ((psi(i,j),i=1,n1),j=1,n2)
     1            ,((vorq(i,j),i=1,n1),j=1,n2)
     1            ,((psc(1,i,j),i=1,n1),j=1,n2)
     1            ,((vor(i,j),i=1,n1),j=1,n2)
     1            ,((vca1(i,j),i=1,n1),j=1,n2)
     1            ,((vca2(i,j),i=1,n1),j=1,n2)
                     endif
                     else
      write(6,203)time,namfile
  203 format(1x,'a t=',e10.3,3x,' forcing',3x,
     1       'write quantities on file=',a22)
      write(nfilo) ((psi(i,j),i=1,n1),j=1,n2)
     1            ,((vorq(i,j),i=1,n1),j=1,n2)
     1            ,((psc(1,i,j),i=1,n1),j=1,n2)
     1            ,((vor(i,j),i=1,n1),j=1,n2)
     1            ,((forc(i,j),i=1,n1),j=1,n2)
     1            ,((vca1(i,j),i=1,n1),j=1,n2)
     1            ,((vca2(i,j),i=1,n1),j=1,n2)
                     endif
      close(nfilo)
c
c   here the profiles of vor, at
c   the north south boundaries are wriiten
c
      namfile='bcnosu'//ipfi//'.dat'
      open(nfilo,file=namfile,form='formatted')
      do i=1,n1
      write(nfilo,*)yp1(i),vor(i,1),vor(i,n2)
      enddo
      close(nfilo)
c
c   here the profiles of vor, psi and passive scalar at
c   the boundaries are wriiten
c
      namfile='bcps'//ipfi//'.dat'
      open(nfilo,file=namfile,form='formatted')
      write(nfilo,*)'  est west boundaries'
      do j=1,n2,2
      write(nfilo,136)yp2(j)
     1                      ,vor(1,j),psi(1,j),psc(1,1,j)
     1                      ,vor(n1,j),psi(n1,j),psc(1,n1,j)
      enddo
      write(nfilo,*)'  sud nord boundaries'
      do i=1,n1,2
      write(nfilo,136)yp1(i)
     1                      ,vor(i,1),psi(i,1),psc(1,i,1)
     1                      ,vor(i,n2),psi(i,n2),psc(1,i,n2)
      enddo
  136 format(13e10.3)
      close(nfilo)
      namfile='bcpsco'//ipfi//'.dat'
      open(nfilo,file=namfile,form='unformatted')
      rewind(nfilo)
      write(nfilo)choe,chow,cves,cven,
     1          dvesto,dpesto,dvweso,dpweso
     1         ,dsesto,dsweso
     1         ,dvnoro,dpnoro,dvsudo,dpsudo
     1         ,dsnoro,dssudo
      close(nfilo)
      if(inbcvs.eq.-2) then
      namfile='insud'//ipfi//'.dat'
      open(47,file=namfile)
      do i=1,n1
      write(47,133) yp1(i),(qinp(l,i),l=1,6)
      enddo
      close(47)
                       endif
      if(inbcvw.eq.-2) then
      namfile='inwest'//ipfi//'.dat'
      open(47,file=namfile)
      do j=1,n2
      write(47,133) yp2(j),(qinp(l,j),l=1,6)
      enddo
  133 format(3x,8e12.4)
      close(47)
                       endif
      if(inbcve.eq.-1) then
      namfile='pwal'//ipfi//'.dat'
      open(47,file=namfile)
      pwal(1)=0.
      do  j=2,n2
      pwal(j)=pwal(j-1)-( (vor(n1,j)+vor(n1,j-1))
     1                   -(vor(n1m,j)+vor(n1m,j-1)) )
     1        *0.5/dx1/(g2m(n1m)*re)
      enddo
      do  j=1,n2
      write(47,101) yp2(j),pwal(j)
      enddo
  101 format(3e13.6)
      close(47)
      namfile='vowal'//ipfi//'.dat'
      open(47,file=namfile)
      do  j=1,n2
      write(47,101) yp2(j),vor(n1,j)
      enddo
      close(47)
                   endif
      return
      end
c
c  ****************************** subrout pricor **********************3
c here the coordiantes are wriiten the file has the
c  format compatible with TURB3D
c
      subroutine pricor(yp1,yp2)
      include 'param.f'
      character*12 namfile
      character*3 ngrid
      dimension yp1(m1),yp2(m2)
      common/dim/n1,n1m,n2,n2m
      write(ngrid,80) n1
   80 format(i3.3)
      nju=1
      aaa=1.
      namfile='coordi.dat'
      open(18,file=namfile,form='unformatted')
      write(18) n1,n2,nju
      write(18) ((yp1(i),i=1,n1),j=1,n2),
     1          ((yp2(j),i=1,n1),j=1,n2),
     1          ((aaa,i=1,n1),j=1,n2)
      close(18)
      return
      end
c
c    subroutine outth                ***************************************
c    output in file time history of certain global quantities
c    calculated in cfield
c
      subroutine outth(ntime,time,vor,psi)
      include 'param.f'
      common/dim/n1,n1m,n2,n2m
      common/d2/nfield,ntst,nhist,nscrn
      common/circul/vort,vortp,vortm,vorw
      common/enes/enerv,ensti,pssqu,enqui
      common/camo/cflm
      dimension vor(m1,m2),psi(m1,m2)
      common/camoc/vocmax,y1cma,y2cma,vocmin,y1cmi,y2cmi
      common/psmoc/pscmax,y1pma,y2pma,pscmin,y1pmi,y2pmi
      common/sfmoc/sfcmax,y1sma,y2sma,sfcmin,y1smi,y2smi
      write(20,158)ntime,time,y1cma,y2cma,vocmax,y1cmi,y2cmi,vocmin,cflm
c     write(6,158)ntime,time,y1cma,y2cma,vocmax,y1cmi,y2cmi,vocmin,cflm
      write(6,158)ntime,time,sfcmax,pscmax,vocmax,
     1                       sfcmin,pscmin,vocmin,cflm
c     write(6,158)time,y1mi,y2mi,vmi,y1ma,vma,vorim,vorip
c    1            ,ensti,enerv,cflm
      write(33,159)y1cma,y2cma
      write(34,159)y1cmi,y2cmi
      write(35,159)time,vort,vortp,vortm,vorw
      write(38,159)time,enerv,ensti,enqui,pssqu
      write(36,159)time,y1cma,y2cma,vocmax,y1cmi,y2cmi,vocmin
      write(37,159)time,y1pma,y2pma,pscmax,y1pmi,y2pmi,pscmin
  158 format(1x,i5,1x,e10.3,1x,10(1x,e10.4))
  159 format(10(1x,e12.5))
      return
      end
c
c  ****************************** subrout cfield **********************
c
c  this subroutine calculate the integral  quantities
c  vorticity and enstrophy and energy
c
c
c     controllare il modo di calcolare l'energia
c
c
      subroutine cfield(vor,psi,psc)
      include 'param.f'
      common/d1/alx1i,alx1f,alx2i,alx2f
      dimension vor(m1,m2),psi(m1,m2),psc(mpsc,m1,m2)
      common/vely/v2(m1,m2)
      common/dim/n1,n1m,n2,n2m
      common/mesh/dx1,dx1q,dx2,dx2q
      common/metri/g2m(m1),g2c(m1)
      common/camoc/vocmax,y1cma,y2cma,vocmin,y1cmi,y2cmi
      common/psmoc/pscmax,y1pma,y2pma,pscmin,y1pmi,y2pmi
      common/sfmoc/sfcmax,y1sma,y2sma,sfcmin,y1smi,y2smi
      common/indbo/imv(m1),ipv(m1)
      common/coor/yp1(m1),yp2(m2)
      common/enes/enerv,ensti,pssqu,enqui
      common/circul/vort,vortp,vortm,vorw
      common/cwalv/inbcvs,inbcvn,inbcvw,inbcve
      common/coorh/ym1(m1)
      common/indpe2/n2i,n2f
      n1m=n1-1
      if(inbcvw.eq.0) then
      n1mt=n1m
                      else
      n1mt=n1m/2+1
                      endif
      vort=0.
      vortp=0.
      vortm=0.
      udx=1./(dx2)
      vorw=0.
      do i=2,n1mt
      vorw=vorw+vor(i,1)
      enddo
      vorw=vorw/float(n1mt)
      do j=2,n2m
      do i=2,n1mt
      undx=udx*(ym1(i)-ym1(i-1))
      if(vor(i,j).gt.0.) vortp=vortp+vor(i,j)*undx
      if(vor(i,j).le.0.) vortm=vortm+vor(i,j)*undx
      vort=vort+vor(i,j)*undx
      enddo
      enddo
      j=n2
      do i=2,n1mt
      undx=udx*(ym1(i)-ym1(i-1))*0.5
      if(vor(i,j).gt.0.) vortp=vortp+vor(i,j)*undx
      if(vor(i,j).le.0.) vortm=vortm+vor(i,j)*undx
      vort=vort+vor(i,j)*undx
      enddo
      j=1
      do i=2,n1mt
      undx=udx*(ym1(i)-ym1(i-1))*0.5
      if(vor(i,j).gt.0.) vortp=vortp+vor(i,j)*undx
      if(vor(i,j).le.0.) vortm=vortm+vor(i,j)*undx
      vort=vort+vor(i,j)*undx
      enddo
      i=1
      do j=2,n2m
      undx=udx*(ym1(i)-yp1(i))
      if(vor(i,j).gt.0.) vortp=vortp+vor(i,j)*undx
      if(vor(i,j).le.0.) vortm=vortm+vor(i,j)*undx
      vort=vort+vor(i,j)*undx
      enddo
      i=n1mt
      do j=2,n2m
      undx=udx*(yp1(i+1)-ym1(i))
      if(vor(i,j).gt.0.) vortp=vortp+vor(i,j)*undx
      if(vor(i,j).le.0.) vortm=vortm+vor(i,j)*undx
      vort=vort+vor(i,j)*undx
      enddo
      vma=-.1e+04
      vmi=+.1e+04
      n1mh=n1m/2+1
      do i=2,n1m
      do j=2,n2m
      if(vor(i,j).lt.vmi) then
      vmi=vor(i,j)
      y1mi=yp1(i)
      y2mi=yp2(j)
      ixmin=i
      jymin=j
      endif
      enddo
      enddo
      vocmin=vmi
      exmin=0.
      eymin=0.
      y1cmi=y1mi+exmin
      y2cmi=y2mi+eymin
      do i=2,n1m
      do j=2,n2m
      if(vor(i,j).gt.vma) then
      vma=vor(i,j)
      y1ma=yp1(i)
      y2ma=yp2(j)
      ixmax=i
      jymax=j
      endif
      enddo
      enddo
      vocmax=vma
      exmax=0.
      eymax=0.
      y1cma=y1ma+exmax
      y2cma=y2ma+eymax
      vma=-1.e+04
      vmi=+.1e+04
      do i=2,n1m
      do j=2,n2m
      if(psc(1,i,j).lt.vmi) then
      vmi=psc(1,i,j)
      y1mi=yp1(i)
      y2mi=yp2(j)
      ixmin=i
      jymin=j
      endif
      enddo
      enddo
      y1pmi=y1mi
      y2pmi=y2mi
      pscmin=vmi
      do i=2,n1m
      do j=2,n2m
      if(psc(1,i,j).gt.vma) then
      vma=psc(1,i,j)
      y1ma=yp1(i)
      y2ma=yp2(j)
      ixmax=i
      jymax=j
      endif
      enddo
      enddo
      y1pma=y1ma
      y2pma=y2ma
      pscmax=vma
      vma=-.1e+04
      vmi=+.1e+04
      n1mh=n1m/2+1
      do i=2,n1m
      do j=2,n2m
      if(psi(i,j).lt.vmi) then
      vmi=psi(i,j)
      y1mi=yp1(i)
      y2mi=yp2(j)
      ixmin=i
      jymin=j
      endif
      enddo
      enddo
      sfcmin=vmi
      y1smi=y1mi+exmin
      y2smi=y2mi+eymin
      do i=2,n1m
      do j=2,n2m
      if(psi(i,j).gt.vma) then
      vma=psi(i,j)
      y1ma=yp1(i)
      y2ma=yp2(j)
      ixmax=i
      jymax=j
      endif
      enddo
      enddo
      sfcmax=vma
      y1sma=y1ma+exmax
      y2sma=y2ma+eymax
      area=(alx1f-alx1i)*(alx2f-alx2i)
      pssqu=0.
      enqui=0.
      ensti=0.
      enerv=0.
      udx=1./(dx1*dx2)
      do i=1,n1 
      undx=g2c(i)*udx
      do j=1,n2f
      enevc=(vor(i,j)*psi(i,j))*undx
      vorqu=(vor(i,j)**4)*undx
      enstc=(vor(i,j)**2)*undx
      psscc=(psc(1,i,j)**2)*undx
      enerv=enerv+enevc
      ensti=ensti+enstc
      enqui=enqui+vorqu
      pssqu=pssqu+psscc
      enddo
      enddo
      enerv=enerv/area
      ensti=ensti/area
      enqui=enqui/area
      pssqu=pssqu/area
      return
      end
c  **************************************************
       subroutine fext(x,y,ze,z,i,j)
       include 'param.f'
      common/mesh/dx1,dx1q,dx2,dx2q
      common/metri/g2m(m1),g2c(m1)

c      see page 296 numerical recipes
c      f(x,y)=f(0,0)+gradf(0).r+.5r.A.r
c      where A is matrix of second derivatives.
c      then note grad f at peak is zero
c      so solve A.r=-gradf(0)
c
c      this subroutine finds the extrema of the modon
c      by iterpolation from the 4 nearest neighbor points centerc      
c      peak on the grid. all distances are measured
c      in terma of d=1. the distance between grid points
c      and must be rescaled for application (i.e.
c      multiply x and y after output by bl/n).

c      dimension z(3,3)
      dimension z(m1,m2)

c grid parameters and central point

c indices for surounding points


c rename points
       zmm=z(i-1,j-1)
       z0m=z(i  ,j-1)
       zpm=z(i+1,j-1)
       zm0=z(i-1,j  )
       z00=z(i  ,j  )
       zp0=z(i+1,j  )
       zmp=z(i-1,j+1)
       z0p=z(i  ,j+1)
       zpp=z(i+1,j+1)

c first derivatives
       fx=(zp0-zm0)*0.5/g2c(i)*dx1
       fy=(z0p-z0m)/2.*dx2

c second derivatives
       fxx=((zp0-z00)/g2m(i)-(z00-zm0)/g2m(i-1))*dx1q/g2c(i)
       fxy=((zpp-zpm)-(zmp-zmm))/4.*dx1*dx2/g2c(i)
       fyy=(z0p-2.*z00+z0m)*dx2q

c hessian (determinant of matrix second derivatives)
       det=fxx*fyy-fxy**2

c position of extremum relative to center at i,j
      if(det.eq.0.) then
         x=0.
         y=0.
      else
       x=(fxy*fy-fyy*fx)/det
       y=(fxy*fx-fxx*fy)/det
      endif
      if( (abs(x).gt.1.).or.(abs(y).gt.1.) ) then
      x=0.
      y=0.
      endif

c value of extremum
       ze=z00+fx*x+fy*y+.5*fxx*x**2+fxy*x*y+.5*fyy*y**2
       return
       end
      function dgal(al,zitaf)
      alzitf=al*zitaf
      dg1=1./(sinh(alzitf)*cosh(alzitf))
      dg2=-alzitf/sinh(alzitf)**2
      dg3=-alzitf/cosh(alzitf)**2
      dgal=dg1+dg2+dg3
      return
      end
c
c
      function gal(al,zitaf)
      alzitf=al*zitaf
      gal=al/(cosh(alzitf)*sinh(alzitf))
      return
      end
