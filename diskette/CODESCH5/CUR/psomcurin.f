c
c  ****************************** subrout dipole  **********************
c
c   initial vort. and psi for Lamb dipole
c
      subroutine dipole(vor,y,ramo,vsi,yc1mo,yc2mo,akmo,velmo,psc)
      include 'param.f'
      dimension vor(m1,m2),psc(m1,m2)
      common/dim/n1,n1m,n2,n2m
      dimension y(ndd,m1,m2)
      common/angmd/thet0
      pi=2.*asin(1.)
      vorpar=2.*vsi*velmo/bessj0(akmo)*akmo
      vorma=-100.
      vormi=+100.
      do i=1,n1
      do j=1,n2
      y1d=y(1,i,j)-yc1mo
      y2d=y(2,i,j)-yc2mo
      ramod=sqrt(y1d**2+y2d**2)
      if(ramod.lt..1e-06)then
      vor(i,j)=0.
      else
      radmod=ramod/ramo
      if(radmod.ge.1.) then
      vor(i,j)=0.
      else
      stho=y2d/ramod
      if(abs(stho).gt.1.) write(6,719) i,j,y2d,ramod
  719 format(3x,2i4,2x,'err in asin',2e12.4)
      astho=(abs(stho)-1.)
c
c     definizione theta
c
      if(astho.gt..0) stho=stho-.1e-06*stho/abs(stho)
      thet=asin(stho)
      if(y2d.lt.0..and.y1d.gt.0.) thet=pi*2.-abs(thet)
      if(y2d.lt.0..and.y1d.lt.0.) thet=pi+abs(thet)
      if(y2d.ge.0..and.y1d.lt.0.) thet=-abs(thet)+pi
c
      sth=sin(thet)*cos(thet0)-sin(thet0)*cos(thet)
c
c   
c   Bessel functions **********************
c   
      arbesj=radmod*akmo
c     write(6,729) i,j,y1d,y2d,ramod,arbesj,thet,sth
  729 format(2i4,3x,2e10.3,3x,4e12.4)
      vor(i,j)=vorpar*sth*bessj1(arbesj)
      vorma=max(vor(i,j),vorma)
      vormi=min(vor(i,j),vormi)
      endif
      endif
      enddo
      enddo
      write(6,791) vorma,vormi
  791 format(3x,'***  in dipol vormax =',e12.4,'  vormi = ',e12.4)
      do i=1,n1
      do j=1,n2
      psc(i,j)=abs(vor(i,j))/vorma
      enddo
      enddo
      return
      end
c
c  ****************************** subrout initia **********************
c
c   initial conditions in the whole field for psi and vor
c
      subroutine initia(vor,psi,y,dvor,psc)
      include 'param.f'
      dimension vor(m1,m2),psi(m1,m2),dvor(m1,m2),psc(m1,m2)
      common/dvowal/dvorbs(m1),dvorbn(m1)
      common/dim/n1,n1m,n2,n2m
      dimension y(ndd,m1,m2),cn2(ndd),cn1(ndd)
      common/metrst/gccc(m1,m2)
      common/metrip/alfi(ndd,ndd,m1,m2),alfj(ndd,ndd,m1,m2)
      common/nonlt/rt(m1,m2),ru(m1,m2)
      common/parmod/yc1mo,yc2mo,akmo,ramo,velmo,vsi
      common/mesh/dx1,dx1q,dx2,dx2q
      common/psbou/psbi(2,m1),psbj(2,m2)
      common/vorwal/vorlo(m1),vorup(m1)
      common/psiwal/psins(m1)
      common/vopsin/psinf(m2),vobj(m2)
      common/uiwai/uinfi(m1),vinfi(m1)
      common/uinfl/uinfl(m2),vinfl(m2)
      common/bdqout/dvest(m2),dpest(m2),dpsest(m2)
      common/slewa/xls(m1)
      common/lewai/xl1,xl2
      common/sliwal/insls,insln
      common/vistr/vampu,vampv
      common/icase/imod,istrea
      ga0=1.
      ro0=0.
      fti=0.
      do j=1,n2
      do i=1,n1
      ru(i,j)=0.
      rt(i,j)=0.
      psc(i,j)=0.
      enddo
      enddo
      if(imod.eq.1) then
c
c    vorticity distribution of a Lamb dipole
c
      arbi=akmo*ramo
      write(6,106)yc1mo,yc2mo,velmo,ramo,akmo,vsi
  106 format(2x,'first  dipole,y1c,y2c,v,r,k,vs',6e10.3)
      call dipole(vor,y,ramo,vsi,yc1mo,yc2mo,akmo,velmo,psc)
      do i =1,n1
      psbi(1,i)=0.
      psbi(2,i)=0.
      uinfi(i)=0.
      vinfi(i)=0.
      enddo
      do j =1,n2
      psbj(1,j)=0.
      psbj(2,j)=0.
      enddo
c
c   from vorticity the streamfunction is calculated
c
      call phcalc(vor,psi)
      t0=0.
      call outpf(vor,psi,t0,psc)
      do i=2,n1m
      vorlo(i)=0.
      vorup(i)=0.
      enddo
      call bcnsvi(ga0,ro0,vor,psi,y,dvor,psc,fti)
      call bcwest(ga0,ro0,vor,psi,psc,y,dvor,psc,fti)
      call bcest(ga0,ro0,vor,psi,psc)
      do i=2,n1m
      vor(i,1)=dvorbs(i)+vor(i,2)
      vor(i,n2)=dvorbn(i)+vor(i,n2m)
      enddo
c    endif  imod
                    endif
      if(istrea.eq.1) then
c
c   here the an inlet profile on the west
c   boundary is assigned
c
      do j=1,n2m
      ym=(y(2,1,j+1)+y(2,1,j))*0.5
      yd=(ym-y(2,1,1))/(y(2,1,n2)-y(2,1,1))
      if(insln.eq.1) then
      uinfl(j)=yd*(1.-yd)/0.25
                     else
      uinfl(j)=yd**(1./7.) 
                     endif
c  endif insln
      vinfl(j)=0.
      enddo
      utrma=0.
      vtrma=0.
c
c    here a inlet velocity profile can be assigned
c    on th sout boundary to get a rough idea of the
c    effects of a discharge on a main stream
c
      do i =1,n1
      uinfi(i)=0.
      vinfi(i)=0.
      if(xls(i).ge.xl1.and.xls(i).le.xl2) then
      xd=(xls(i)-xl1)/(xl2-xl1)
      vinfi(i)=xd*(1.-xd)/0.25*vampv
      uinfi(i)=xd*(1.-xd)/0.25*vampu
c    endif  xls(i
                                endif
      utrma=max(utrma,abs(uinfi(i)))
      vtrma=max(vtrma,abs(vinfi(i)))
      enddo
      write(6,*)' traspirat vel max',utrma,vtrma
      is=1
      i=1
      psbj(1,1)=0.
      psbj(2,1)=0.
      do j =1,n2m
      do nd=1,ndd
      cn1(nd)=-is*((y(nd,i,j)-y(nd,i+is,j))
     1            +(y(nd,i,j+1)-y(nd,i+is,j+1)))*dx1*0.5
      cn2(nd)=(y(nd,i,j+1)-y(nd,i,j))*dx2
      enddo
      q1=cn1(1)*uinfl(j)+cn1(2)*vinfl(j)
      psbj(1,j+1)=psbj(1,j)+q1/dx1
      psbj(2,j+1)=psbj(1,j+1)    
      enddo
      do j =1,n2
      psinf(j)=0.
      enddo
      do i =1,n1
      psins(i)=0.
      psbi(1,i)=0.
      psbi(2,i)=psbj(1,n2)
      enddo
      do j =1,n2
      do i =1,n1
      psi(i,j)=psbj(1,j)
      enddo
      enddo
      do j =2,n2m
      vobj(j)=-(alfj(2,2,1,j)*(psbj(1,j+1)-psbj(1,j))
     1         -alfj(2,2,1,j-1)*(psbj(1,j)-psbj(1,j-1)))
     1        *dx2q/gccc(1,j)
      enddo
      i=1
      ual2=-1./gccc(i,2)
      ual1=-1./gccc(i,1)
      dq12=-(ual2*(alfj(2,2,i,2)*(psbj(i,3)-psbj(i,2))
     1            -alfj(2,2,i,1)*(psbj(i,2)-psbj(i,1)))
     1    -2.*ual1*alfj(2,2,i,1)*(psbj(i,2)-psbj(i,1)))*dx2q
      vorjb=dq12
      vobj(1)=vobj(2)+vorjb 
      ual2=-1./gccc(i,n2m)
      ual1=-1./gccc(i,n2)
      dq12=-(ual2*(alfj(2,2,i,n2m-1)*(psbj(1,n2m)-psbj(1,n2m-1))
     1            -alfj(2,2,i,n2m)*(psbj(1,n2)-psbj(1,n2m)))
     1    -2.*ual1*alfj(2,2,i,n2m)*(psbj(1,n2)-psbj(1,n2m)))*dx2q
      vorjb=dq12
      vobj(n2)=vobj(n2m)+vorjb 
      call vorcal(vor,psi)
      call bcnsvi(ga0,ro0,vor,psi,y,dvor,psc,fti)
      call bcwest(ga0,ro0,vor,psi,psc,y,dvor,psc,fti)
      call bcest(ga0,ro0,vor,psi,psc)
c    endif istrea
                    endif
      do i=2,n1m
      vor(i,1)=vor(i,2)+dvorbs(i)
      vor(i,n2)=vor(i,n2m)+dvorbn(i)
      enddo
      do j=1,n2
      vor(n1,j)=vor(n1m,j)
      enddo
      open(27,file='inico.out')
      write(27,*)'  inlet'
      do j=1,n2
      write(27,132)j,vobj(j),vor(2,j),psbj(1,j),uinfl(j)
  132 format(3x,i4,2x,5e12.5)
      enddo
      write(27,*)'  outlet'
      do j=1,n2
      write(27,132)j,vor(n1,j),vor(n1m,j),psi(n1,j)
      enddo
      close(27)
      vorma=0.
      psima=0.
      pscma=0.
      do i=1,n1
      do j=1,n2
      pscma=max(pscma,abs(psc(i,j)))
      psima=max(psima,abs(psi(i,j)))
      vorma=max(vorma,abs(vor(i,j)))
      enddo
      enddo
      write(6,*)'***** initia pscma vorma,psima  '
     1          ,pscma,vorma,psima
      return
      end
c
c  BESSEL FUNCTIONS
c
      function bessj1(x)
      include 'param.f'
      real*8 y,p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,
     *    s1,s2,s3,s4,s5,s6
      data r1,r2,r3,r4,r5,r6/72362614232.d0,-7895059235.d0,242396853.1d0
     *,
     *    -2972611.439d0,15704.48260d0,-30.16036606d0/,
     *    s1,s2,s3,s4,s5,s6/144725228442.d0,2300535178.d0,
     *    18583304.74d0,99447.43394d0,376.9991397d0,1.d0/
      data p1,p2,p3,p4,p5/1.d0,.183105d-2,-.3516396496d-4,.2457520174d-5
     *,
     *    -.240337019d-6/, q1,q2,q3,q4,q5/.04687499995d0,-.2002690873d-3
     *,
     *    .8449199096d-5,-.88228987d-6,.105787412d-6/
      if(abs(x).lt.8.)then
        y=x**2
        bessj1=x*(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))
     *      /(s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6)))))
      else
        ax=abs(x)
        z=8./ax
        y=z**2
        xx=ax-2.356194491
        bessj1=sqrt(.636619772/ax)*(cos(xx)*(p1+y*(p2+y*(p3+y*(p4+y
     *      *p5))))-z*sin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))
     *      *sign(1.,x)
      endif
      return
      end
      function bessj0(x)
      include 'param.f'
      real*8 y,p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,
     *    s1,s2,s3,s4,s5,s6
      data p1,p2,p3,p4,p5/1.d0,-.1098628627d-2,.2734510407d-4,
     *    -.2073370639d-5,.2093887211d-6/, q1,q2,q3,q4,q5/-.1562499995d-
     *1,
     *    .1430488765d-3,-.6911147651d-5,.7621095161d-6,-.934945152d-7/
      data r1,r2,r3,r4,r5,r6/57568490574.d0,-13362590354.d0,651619640.7d
     *0,
     *    -11214424.18d0,77392.33017d0,-184.9052456d0/,
     *    s1,s2,s3,s4,s5,s6/57568490411.d0,1029532985.d0,
     *    9494680.718d0,59272.64853d0,267.8532712d0,1.d0/
      if(abs(x).lt.8.)then
        y=x**2
        bessj0=(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))
     *      /(s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6)))))
      else
        ax=abs(x)
        z=8./ax
        y=z**2
        xx=ax-.785398164
        bessj0=sqrt(.636619772/ax)*(cos(xx)*(p1+y*(p2+y*(p3+y*(p4+y
     *      *p5))))-z*sin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))
      endif
      return
      end
