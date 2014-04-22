c
c  ****************************** subrout coordi  **********************
c    this subroutine assign the physical coordinates as function
c    of the computational ones x1,x2
 
      subroutine coordi(y)
      include 'param.f'
      common/dim/n1,n1m,n2,n2m
      common/corbum/bum1,sigb1,bum2,sigb2,xp1bu,xp2bu 
      common/cordat/alx1i,alx1f,alx2i,alx2f
      common/strx1/str1,xcr1,etr1,istr1
      common/strx2/str2,xcr2,etr2,istr2
      common/lewai/xl1,xl2
      dimension y(ndd,m1,m2),ycb1(m1),ycb2(m1),x1(m1)
      dimension xt11(m1),xt12(m1),et1(m1)
      dimension xt21(m2),xt22(m2),et2(m2)
      common/slewa/xls(m1)
      common/itypco/itypc
      pi=2.*asin(1.)
c
c  non-uniform grid in x1
c
c  here a first coordinate transformation is used to stretch
c  or cluster the grid points in the region of interest
c
      if(istr1.eq.1) then
c
c  uniform in the center stretched on the sides
c
      nc=n1m*etr1
      xc=xcr1
      xf=1.
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
   43 continue
      res=gal1(al,zitaf)-val
      if(abs(res).lt..1e-04) then
      go to 45
                             else
      tang=dgal1(al,zitaf)
      dal=res/tang
      write(6,*)iter,res,al,tang,dal
      al=al-dal
      iter=iter+1
      if(iter.gt.21) go to 45
      go to 43
                             endif
   45 continue
      write(6,*)'convergence',iter,res,al
      gam=tanh(al*zitaf)                                                      
      do i=1,n1                                                      
      csi=(i-1)/float(n1m)                                   
      etal=csi-0.5
      if(abs(etal).le.etac)  then
      et1(i)=etal*xc/etac
                            else 
      siet=etal/abs(etal) 
      if(siet.ge.0) then
      zita=(etal-etac)
      azit=al*(zita-zitaf)
      et1(i)=(xf+tanh(azit)/gam*(xf-xc))
                            else 
      zita=(etal+etac)
      azit=al*(zita+zitaf)
      et1(i)=(-xf+tanh(azit)/gam*(xf-xc))
                            endif                                   
                            endif                                   
c     write(18,110)csi,et1(i)
  110 format(3x,2e12.5)
      enddo
      do  i=1,n1
      do  j=1,n2
      y(1,i,j)=alx1f*et1(i)
      enddo
      enddo
                     endif
      if(istr1.eq.0) then
c
c  uniform 
c
      do  i=1,n1
      et1(i)=float(i-1)/n1m
      enddo
      do  i=1,n1
      do  j=1,n2
      y(1,i,j)=alx1i+(alx1f-alx1i)*et1(i)
      enddo
      enddo
                    endif
      if(istr1.eq.2) then
c
c  stretched after the bump
c
      tstr=tanh(str1*etr1)
      xmed1=xcr1
      strb=str1+2.
      do  i=1,n1
      x1(i)=float(i-1)/n1m
      xt11(i)=xmed1/alx1f*tanh(str1*x1(i))/tstr
      enddo
      do i=1,n1
      xt12(i)=1./xt11(n1)+(1.-1./xt11(n1))
     1     *tanh(strb*(x1(i)-1.))/tanh(strb*(etr1-1.))
       et1(i)= xt11(i)*xt12(i)
      do  j=1,n2
      y(1,i,j)=alx1i+(alx1f-alx1i)*et1(i)
      enddo
      enddo
                    endif
c
c  irregular boundaries at j=1 and j=n2
c
      if(itypc.eq.0) then
c   here the coplex boundary is given by an gaussian distribution
c   to simulate a bump or a depression
c
      bum1c=bum1*bum2
      bum1t=bum1-bum1c
      if(bum2.gt.0.) then 
      cost0=bum1t/bum1
      the0=acos(cost0)
      d0=sigb1*bum1*sin(the0)
      al0=(alx2i-d0)*tan(the0)/d0
      al=1.
      raf=abs(y(1,1,1)-xp1bu)
      arf=(raf-bum1t)
      iter=1
      val=al0
   33 continue
      res=gal(al,arf)-val
      if(abs(res).lt..1e-04) then
      go to 35
                             else
      tang=dgal(al,arf)
      dal=res/tang
      write(6,*)iter,res,al,tang,dal
      al=al+dal
      iter=iter+1
      if(iter.gt.21) go to 35
      go to 33
                             endif
   35 continue
      write(6,*)'convergence',iter,res,al
      red=d0+(alx2i-d0)
      raf=abs(y(1,1,1)-xp1bu)
      arf=(raf-bum1t)*al
      write(6,*)' coord. smooth ,d0,al,red ',d0,al,red
      icor=1
                             else
      icor=0
                             endif
      do  i=1,n1
      rad=abs(y(1,i,1)-xp1bu)
      if(rad.ge.bum1t) then
      if(icor.eq.1) then
      arg=(rad-bum1t)*al
      ycb1(i)=d0+(alx2i-d0)*tanh(arg)/tanh(arf)
                    else
      ycb1(i)=alx2i
                    endif
                      else
      cost=rad/bum1
      the=acos(cost)
      ycb1(i)=alx2i+sigb1*bum1*sin(the)
                      endif
      ycb2(i)=alx2f
      enddo
                       else
      do  i=1,n1
      ycb1(i)=alx2i+bum1*exp(-((y(1,i,1)-xp1bu)/sigb1)**2)
      ycb2(i)=alx2f+bum2*exp(-((y(1,i,n2)-xp2bu)/sigb2)**2)
      enddo
                      endif

c
c   non-uniform in x2 or j
c
      if(istr2.eq.-2) then
c
c  clustered near north and south bound.
c
      al=str2
      gam=tanh(al*0.5)
      do j=1,n2
      csi=(j-1)/float(n2m)
      etal=al*(csi-0.5)
      xt1=0.5*(1.+tanh(etal)/gam)
      et2(j)=(-0.5+xt1)*2.
      enddo
                     endif
      if(istr2.eq.0) then
c
c   uniform 
c
      do j=1,n2                                                      
      csi=(j-1)/float(n2m)                                   
      et2(j)=csi
      enddo
                       endif
      if(istr2.eq.2) then
c
c  clustered near the south boundary
c
      al=str2
      gam=tanh(al)
      do j=1,n2
      csi=(j-1)/float(n2m)
      etal=al*(csi-1.)
      xt1=1.+tanh(etal)/gam
      et2(j)=xt1
      write(18,110)csi,et2(j)
   
      enddo
                       endif
      if(istr2.eq.3) then
c
c   stretched in one region and clustered in a different one
c    coord. 2.4.f
c
      nc=n2m*etr2
      xc=alx2f*xcr2
      xf=alx2f
      al1=str2
      al2=al1+1.
      csic=(nc-1)/float(n2m)
      gam=tanh(al1*csic)                                                      
      do j=1,n2
      csi=(j-1)/float(n2m)
      xt21(j)=xc/xf*tanh(al1*csi)/gam
      enddo
      do j=1,n2
      csi=(j-1)/float(n2m)
      xt22(j)=1./xt21(n2)+(1.-1./xt21(n2))
     1     *tanh(al2*(csi-1.))/tanh(al2*(csic-1.))
      et2(j)=xt22(j)*xt21(j)*xf
      enddo
                           endif
      do  i=1,n1
      do  j=1,n2
      y(2,i,j)=ycb1(i)+(ycb2(i)-ycb1(i))*et2(j)
      enddo
      enddo
c
c    here the coordiante is written in a file to be used in
c    graphics packages here the format if for TURB3D 
c
      write(6,*)' extrema points of ccord'
      write(6,102) y(1,1,1),y(2,1,1),y(1,1,n2),y(2,1,n2)
      write(6,102) y(1,n1,1),y(2,n1,1),y(1,n1,n2),y(2,n1,n2)
  102 format(3x,4e11.4)
      open(18,file='coordi.dat',form='formatted')
      write(18,*)n1,n2
      write(18,*) ((y(1,i,j),i=1,n1),j=1,n2)
     1         ,((y(2,i,j),i=1,n1),j=1,n2)
      close(18)
      xls(1)=0.
      iin=0
      il1=0
      il2=0
      do i =2,n1
c     xse=sqrt((y(1,i,1)-y(1,i-1,1))**2
c    1        +(y(2,i,1)-y(2,i-1,1))**2)
c     xls(i)=xls(i-1)+xse
      xls(i)=y(1,i,1)
      if(xls(i).gt.xl1.and.iin.eq.0) then
      iin=1
      il1=i
                        endif
      if(xls(i).ge.xl2.and.iin.eq.1) then
      il2=i
      go to 133
                        endif
      enddo
  133 continue
      write(6,*)'xl1,xl2   ',xl1,xl2,'  il1,il2 ',il1,il2
      write(20,*)'xl1,xl2   ',xl1,xl2,'  il1,il2 ',il1,il2
      return
      end
      function dgal(al,arf)
      alzitf=al*arf
      dg1=-tanh(alzitf)
      dg2=alzitf/cosh(alzitf)**2
      dgal=(dg1+dg2)/al**2
      return
      end
c
c
      function gal(al,arf)
      alzitf=al*arf
      gal=tanh(alzitf)/al
      return
      end
      function dgal1(al,zitaf)
      alzitf=al*zitaf
      dg1=1./(sinh(alzitf)*cosh(alzitf))
      dg2=-alzitf/sinh(alzitf)**2
      dg3=-alzitf/cosh(alzitf)**2
      dgal1=dg1+dg2+dg3
      return
      end
c
c
      function gal1(al,zitaf)
      alzitf=al*zitaf
      gal1=al/(cosh(alzitf)*sinh(alzitf))
      return
      end
c
c  ****************************** subrout pmetri  **********************
c   here the metric quantities are written in a file
c   to evaluate their distribution in the computational 
c   space
c  This routine is only for check and is not important
c
      subroutine pmetri(y)
      include 'param.f'
      common/dim/n1,n1m,n2,n2m
      common/mesh/dx1,dx1q,dx2,dx2q
      common/metrip/alfi(ndd,ndd,m1,m2),alfj(ndd,ndd,m1,m2)
      common/metrst/gccc(m1,m2)
      dimension y(ndd,m1,m2)
      common/metriv/alfc(ndd,ndd,m1,m2)
      dimension yp1(m1,m2),yp2(m1,m2)
      do i=1,n1
      do j=1,n2m
      yp1(i,j)=(y(1,i,j)+y(1,i,j+1))*0.5
      yp2(i,j)=(y(2,i,j)+y(2,i,j+1))*0.5
      enddo
      enddo
      open(18,file='cordjh.dat',form='unformatted')
      write(18)n1,n2m
      write(18) ((yp1(i,j),i=1,n1),j=1,n2m)
     1         ,((yp2(i,j),i=1,n1),j=1,n2m)
      close(18)
      time=0.
      open(17,file='metrj.dat',form='unformatted')
      write(17)n1,n2m
      write(17) time,time,time,time
      write(17)   ((alfj(1,1,i,j),i=1,n1),j=1,n2m)
     1           ,((alfj(1,2,i,j),i=1,n1),j=1,n2m)
     1           ,((alfj(2,1,i,j),i=1,n1),j=1,n2m)
     1           ,((alfj(2,2,i,j),i=1,n1),j=1,n2m)
      close(17)
      do i=1,n1m
      do j=1,n2
      yp1(i,j)=(y(1,i+1,j)+y(1,i,j))*0.5
      yp2(i,j)=(y(2,i+1,j)+y(2,i,j))*0.5
      enddo
      enddo
      open(18,file='cordih.dat',form='unformatted')
      write(18)n1m,n2
      write(18) ((yp1(i,j),i=1,n1m),j=1,n2)
     1         ,((yp2(i,j),i=1,n1m),j=1,n2)
      close(18)
      time=0.
      open(17,file='metri.dat',form='unformatted')
      write(17)n1m,n2
      write(17) time,time,time,time
      write(17)   ((alfi(1,1,i,j),i=1,n1m),j=1,n2)
     1           ,((alfi(1,2,i,j),i=1,n1m),j=1,n2)
     1           ,((alfi(2,1,i,j),i=1,n1m),j=1,n2)
     1           ,((alfi(2,2,i,j),i=1,n1m),j=1,n2)
      close(17)
      time=0.
      open(17,file='metrc.dat',form='unformatted')
      write(17)n1,n2
      write(17) time,time,time,time
      write(17)   ((alfc(1,1,i,j),i=1,n1),j=1,n2)
     1           ,((alfc(1,2,i,j),i=1,n1),j=1,n2)
     1           ,((alfc(2,1,i,j),i=1,n1),j=1,n2)
     1           ,((alfc(2,2,i,j),i=1,n1),j=1,n2)
      close(17)
      return
      end
c
c  ****************************** subrout metric  **********************
c  this subroutine performs the calculation of the metric quantities
c  at i,j ; i+1/2,j ; i,j+1/2 ; 1+1/2,j+1/2
c  ***** the quantities with  j(  evaluated at i,j+1/2      positiion
c  ***** the quantities with  i(  evaluated at i+1/2,j      positiion
c  ***** the quantities with  c(  evaluated at i,j          positiion
c  ***** the quantities with  m(  evaluated at i+1/2,j+1/2  positiion
c  ***** the metric are evaluated by centered differences.
c  g are the jacobians
c
      subroutine metric(y)
      include 'param.f'
      common/dim/n1,n1m,n2,n2m
      common/mesh/dx1,dx1q,dx2,dx2q
      common/metrip/alfi(ndd,ndd,m1,m2),alfj(ndd,ndd,m1,m2)
      common/metrih/cac(ndd,ndd,m1,m2)
      common/metrst/gccc(m1,m2)
      dimension y(ndd,m1,m2),cn2(ndd),cn1(ndd),gn12x1(ndd)
      common/metrbw/gwa112(2,m1)
      common/metriv/alfc(ndd,ndd,m1,m2)
c
c  *********                 i,j+1/2
c
      do i=2,n1m
      ip=i+1
      im=i-1
      do j=1,n2m
      do nd=1,ndd
c
c  ********* derivatives of cartesian corrdinates interior points
c  cn1 deriv respect to x1 , cn2 deriv respect to x2
c
      cn1(nd)=((y(nd,ip,j+1)+y(nd,ip,j))
     1       -(y(nd,im,j+1)+y(nd,im,j)))*dx1/4.
      cn2(nd)=(y(nd,i,j+1)-y(nd,i,j))*dx2
      enddo
      gjj=(cn1(1)*cn2(2)-cn1(2)*cn2(1))
      if(gjj.lt..1e-10)gjj=.1e-10
      cbj11=cn2(2)/gjj
      cbj12=-cn2(1)/gjj
      cbj21=-cn1(2)/gjj
      cbj22=cn1(1)/gjj
      alfj(2,1,i,j)=gjj*(cbj21*cbj11+cbj22*cbj12)
      alfj(2,2,i,j)=gjj*(cbj21**2+cbj22**2)
      alfj(1,2,i,j)=gjj*(cbj11*cbj21+cbj12*cbj22)
      alfj(1,1,i,j)=gjj*(cbj11**2+cbj12**2)
      enddo
      enddo
c
c  *********    i,j+1/2  at i=1 and i=n1
c
      do i=1,n1,n1m
      if(i.eq.1) is=1
      if(i.eq.n1) is=-1
      ip=i+1
      im=i-1
      do j=1,n2m
      do nd=1,ndd
c
c  ********* derivatives of cartesian corrdinates interior points
c  cn1 deriv respect to x1 , cn2 deriv respect to x2
c
      cn1(nd)=-is*((y(nd,i,j)-y(nd,i+is,j))
     1            +(y(nd,i,j+1)-y(nd,i+is,j+1)))*dx1*0.5
      cn2(nd)=(y(nd,i,j+1)-y(nd,i,j))*dx2
      enddo
      gjj=(cn1(1)*cn2(2)-cn1(2)*cn2(1))
      if(gjj.lt..1e-10)gjj=.1e-10
      cbj11=cn2(2)/gjj
      cbj12=-cn2(1)/gjj
      cbj21=-cn1(2)/gjj
      cbj22=cn1(1)/gjj
      alfj(2,1,i,j)=gjj*(cbj21*cbj11+cbj22*cbj12)
      alfj(2,2,i,j)=gjj*(cbj21**2+cbj22**2)
      alfj(1,2,i,j)=gjj*(cbj11*cbj21+cbj12*cbj22)
      alfj(1,1,i,j)=gjj*(cbj11**2+cbj12**2)
      enddo
      enddo
c
c  *********                 i+1/2,j
c
      do i=1,n1m
      ip=i+1
      do j=2,n2m
      do nd=1,ndd
c
c  ********* derivatives of cartesian corrdinates interior points
c
      cn2(nd)=((y(nd,ip,j+1)+y(nd,i,j+1))
     1       -(y(nd,ip,j-1)+y(nd,i,j-1)))*dx2/4.
      cn1(nd)=(y(nd,ip,j)-y(nd,i,j))*dx1
      enddo
      gji=(cn1(1)*cn2(2)-cn1(2)*cn2(1))
      if(gji.lt..1e-10)gji=.1e-10
      cbi11=cn2(2)/gji
      cbi12=-cn2(1)/gji
      cbi21=-cn1(2)/gji
      cbi22=cn1(1)/gji
      alfi(1,1,i,j)=gji*(cbi11**2+cbi12**2)
      alfi(1,2,i,j)=gji*(cbi11*cbi21+cbi12*cbi22)
      alfi(2,1,i,j)=gji*(cbi21*cbi11+cbi22*cbi12)
      alfi(2,2,i,j)=gji*(cbi21**2+cbi22**2)
      enddo
      enddo
c
c  *********   i+1/2,j  at j=1 and n2
c
      do j=1,n2,n2m
      if(j.eq.1) js=1
      if(j.eq.n2) js=-1
      do i=1,n1m
      ip=i+1
      do nd=1,ndd
c
c  ********* derivatives of cartesian corrdinates interior points
c
      cn2(nd)=-js*((y(nd,i,j)-y(nd,i,j+js))
     1            +(y(nd,ip,j)-y(nd,ip,j+js)))*dx2*0.5
      cn1(nd)=(y(nd,ip,j)-y(nd,i,j))*dx1
      enddo
      gji=(cn1(1)*cn2(2)-cn1(2)*cn2(1))
      if(gji.lt..1e-10)gji=.1e-10
      cbi11=cn2(2)/gji
      cbi12=-cn2(1)/gji
      cbi21=-cn1(2)/gji
      cbi22=cn1(1)/gji
      alfi(1,1,i,j)=gji*(cbi11**2+cbi12**2)
      alfi(1,2,i,j)=gji*(cbi11*cbi21+cbi12*cbi22)
      alfi(2,1,i,j)=gji*(cbi21*cbi11+cbi22*cbi12)
      alfi(2,2,i,j)=gji*(cbi21**2+cbi22**2)
      enddo
      enddo
c
c  *********                 i,j
c
      do i=1,n1,n1m
      if(i.eq.1) is=1
      if(i.eq.n1) is=-1
      do j=1,n2
c
c  ********* derivatives of cartesian corrdinates vertical boundary
c            three points backw.
      do nd=1,ndd
      cn1(nd)=-is*(y(nd,i,j)-y(nd,i+is,j))*dx1
      if(j.gt.1.and.j.lt.n2) then
      cn2(nd)=(y(nd,i,j+1)-y(nd,i,j-1))*dx2*0.5
      else
      if(j.eq.1) js=1
      if(j.eq.n2) js=-1
      cn2(nd)=-js*(y(nd,i,j)-y(nd,i,j+js))*dx2
      endif
      enddo
      cac(2,2,i,j)=cn2(2)
      cac(1,2,i,j)=cn2(1)
      cac(2,1,i,j)=cn1(2)
      cac(1,1,i,j)=cn1(1)
      gjc=(cn1(1)*cn2(2)-cn1(2)*cn2(1))
      if(gjc.lt..1e-10)write(6,714)i,j,gjc
      if(gjc.lt..1e-10)gjc=.1e-10
      gccc(i,j)=gjc
      cbc11=cn2(2)/gjc
      cbc12=-cn2(1)/gjc
      cbc21=-cn1(2)/gjc
      cbc22=cn1(1)/gjc
      alfc(1,1,i,j)=gjc*(cbc11**2+cbc12**2)
      alfc(1,2,i,j)=gjc*(cbc11*cbc21+cbc12*cbc22)
      alfc(2,1,i,j)=gjc*(cbc21*cbc11+cbc22*cbc12)
      alfc(2,2,i,j)=gjc*(cbc21**2+cbc22**2)
      enddo
      enddo
c
      do j=1,n2,n2m
      if(j.eq.1) js=1
      if(j.eq.n2) js=-1
      do i=2,n1m
c
c  ********* derivatives of cartesian corrdinates horizontal boundary
c            three points backw.
      do nd=1,ndd
      cn1(nd)=(y(nd,i+1,j)-y(nd,i-1,j))*dx1*0.5
      cn2(nd)=-js*(y(nd,i,j)-y(nd,i,j+js))*dx2
      enddo
      cac(2,2,i,j)=cn2(2)
      cac(1,2,i,j)=cn2(1)
      cac(2,1,i,j)=cn1(2)
      cac(1,1,i,j)=cn1(1)
      gjc=(cn1(1)*cn2(2)-cn1(2)*cn2(1))
      if(gjc.lt..1e-10)write(6,714)i,j,gjc
      if(gjc.lt..1e-10)gjc=.1e-10
      gccc(i,j)=gjc
      cbc11=cn2(2)/gjc
      cbc12=-cn2(1)/gjc
      cbc21=-cn1(2)/gjc
      cbc22=cn1(1)/gjc
      alfc(1,1,i,j)=gjc*(cbc11**2+cbc12**2)
      alfc(1,2,i,j)=gjc*(cbc11*cbc21+cbc12*cbc22)
      alfc(2,1,i,j)=gjc*(cbc21*cbc11+cbc22*cbc12)
      alfc(2,2,i,j)=gjc*(cbc21**2+cbc22**2)
      enddo
      enddo
c
      do i=2,n1m
      do j=2,n2m
c
c  ********* derivatives of cartesian corrdinates interior points
c
      do nd=1,ndd
      cn1(nd)=(y(nd,i+1,j)-y(nd,i-1,j))*dx1*0.5
      cn2(nd)=(y(nd,i,j+1)-y(nd,i,j-1))*dx2*0.5
      enddo
      cac(2,2,i,j)=cn2(2)
      cac(1,2,i,j)=cn2(1)
      cac(2,1,i,j)=cn1(2)
      cac(1,1,i,j)=cn1(1)
      gjc=(cn1(1)*cn2(2)-cn1(2)*cn2(1))
      if(gjc.lt..1e-10)write(6,714)i,j,gjc
  714 format(3x,'i,j',2i4,3x,e10.4)
      if(gjc.lt..1e-10)gjc=.1e-10
      gccc(i,j)=gjc
      cbc11=cn2(2)/gjc
      cbc12=-cn2(1)/gjc
      cbc21=-cn1(2)/gjc
      cbc22=cn1(1)/gjc
      alfc(1,1,i,j)=gjc*(cbc11**2+cbc12**2)
      alfc(1,2,i,j)=gjc*(cbc11*cbc21+cbc12*cbc22)
      alfc(2,1,i,j)=gjc*(cbc21*cbc11+cbc22*cbc12)
      alfc(2,2,i,j)=gjc*(cbc21**2+cbc22**2)
      enddo
      enddo
c
c     metric quantities in the wall vorticity
c
c
      do j=1,n2,n2m
      if(j.eq.1) ja=1
      if(j.eq.n2) ja=2
      do i=2,n1m
c
c  ********* derivatives of c12 and c21 along horizontal boundary
c
      do nd=1,ndd
      gn12x1(nd)=(cac(nd,2,i+1,j)-cac(nd,2,i-1,j))*dx1*0.5
     1             *(cac(nd,1,i,j)*alfc(1,1,i,j)
     1              +cac(nd,2,i,j)*alfc(1,2,i,j))
     1             /(gccc(i,j)**2*alfc(1,1,i,j))
      enddo
      gwa112(ja,i)=gn12x1(1)+gn12x1(2)
      enddo
      enddo
c
      return
      end
