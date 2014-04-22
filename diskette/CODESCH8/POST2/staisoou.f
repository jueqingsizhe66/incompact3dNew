c
c  ****************************** subrout outpf  **********************
c
      subroutine outpf(time,enen,nav,q,pr,qcap)
c
      include 'param.f'
      parameter (m1m=m1-1,m2m=m2-1,m3m=m3-1)
      parameter (mlr=m3m/2+1)
      common/dim/n1,n1m,n2,n2m,n3,n3m
      dimension pr(m1,m2,m3),sca(m1,m2,m3)
      dimension qcap(m1,m2,m3)
      dimension q(ndv,m1,m2,m3)
      dimension ru(ndv,m1,m2,m3)
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      common/tstep/dt,beta,ren
      dimension vprms(4),vmp(3),skp(4),flp(4)
      common/filep/ifilp
      common/vmean/vm(ndv),vrms(4)
      common/omean/vo(ndv),vorms(3)
      common/eneav/enav,enavo
      common/skfl/ske(4),fla(4)
      common/prrm/prm,prms
      common/vmeao/vmo(ndv),vrmso(4)
      common/names/filcnw,filcnr,filth,filou,filuu,filsf
      common/corpri/yp1(m1),yp2(m2),yp3(m3)
      common/speene/e(ndv,0:m3)
      common/corre/corf(ndv,m3)
      dimension ruu(ndv,m3),roo(ndv,m3)
      common/numond/kx(m1),ky(m2),kz(m3),kkmax
      common/ledat/ics0,cvisc,ifiltr,ibox
      common/kolqua/diss,eta
      common/ispec/imic
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/papdf/nlr,nlrju,lipdf,lipdh,lipdn
      common/paqdf/sfuma,qquma,sfoma,qqoma
      common/itystf/inostf
      common/alpdf/alasf(201),alaqu(201),alahe(201),aladi(201)
      common/nllpdf/nllsf(9,mlr,201),nllqu(7,201),nllsc(7,201)
      common/noupdf/nousf(9,mlr),nouqu(7),nousc(7)
      common/nouto/ntqll(7),ntqou(7),ntsll(7),ntsou(7)
      common/nsfto/ntsfl(9,mlr),ntsfu(9,mlr)
      dimension nllq(201)
      common/pdfqdf/pdfsf(9,mlr,201),pdfqu(7,201),pdfsc(7,201)
      common/ncopdf/nllco(6,201),ntcoll(6)
      common/pdfcog/pdfco(6,201),pdfcot(6)
c      
      character*27 filcnw,filcnr,filth,filou,filuu,filsf
      character*17 filcos
      character*4 pntim
      character*6 slp
      character*1 is1s,is3s
      character*1 is1n,is3n
      character*60 namfi3
      character*3 nkpse,njpse
      itim=time+0.3
      write(pntim,77) itim
   77 format(i4.4)
      nav=1
      if(imic.gt.0) then
c
c  statistics calculation
c
      call enerca(q,pr,enen,time)
      call pdiss(q,diss,sca)
      call vort(q,ru)
c     print*,'in outpf'
c
      dissve=2*cvisc*diss
      enavp=enen
      filcos='spenc.'//pntim
      open(57,file=filcos)
      rewind(57)
      filcos='spenko.'//pntim
      open(59,file=filcos)
      rewind(59)
      do 615 l=1,3
      skp(l)=ske(l)
      flp(l)=fla(l)
      vmp(l)=vm(l)
      vrms(l)=sqrt(vrms(l))
      vorms(l)=sqrt(vorms(l))
  615 continue
      skp(4)=ske(4)
      flp(4)=fla(4)
      vprms(4)=vrms(4)
      ppmp=prms
  611 continue
  612 format(1x,i3,2x,e14.7,2x,9(1x,e11.4)) 
  613 format(e14.7,1x,9(1x,e11.4)) 
      close(42)
      close(62)
c
c  evaluation of the velocity spectra
c
      call spectre(q,qcap)
c
c   Kolmogorov scaling
c
      speska=(cvisc**5./dissve)**(1./4.)
      eta=(cvisc**3./dissve)**(1./4.)
       ee1t0=0.
       ee2t0=0.
       ee3t0=0.
       en1t0=0.
       en2t0=0.
       en3t0=0.
       do k=1,kkmax
       write (57,*)k,e(1,k),e(2,k),e(3,k)
       ee1to=ee1to+e(1,k)
       ee2to=ee2to+e(2,k)
       ee3to=ee3to+e(3,k)
       en1to=en1to+e(1,k)*k**2
       en2to=en2to+e(2,k)*k**2
       en3to=en3to+e(3,k)*k**2
       akol=k*eta
       se1=e(1,k)/speska
       se2=e(2,k)/speska
       se3=e(3,k)/speska
       write (59,*)akol,se1,se2,se3          
       end do
       eeto=ee1to+ee2to+ee3to
       ento=en1to+en2to+en3to
      close(57)
      close(59)
      write(6,133)time,ee1to,ee2to,ee3to,eeto,enen
  133 format('energy',e12.4,3x,3e12.4,3x,'etot sp,phy',2e12.5)
      write(6,113)time,en1to,en2to,en3to,ento
  113 format('en k^2',e12.4,3x,3e12.4,3x,'entot sp,phy',2e12.5)
       disssp=2.*cvisc*ento
      flamb=sqrt(20*cvisc*enen/disssp)
      rlamb=enen*sqrt(20./(3.*cvisc*disssp))
      write(6,*)'eta=',eta,'   flamb=',flamb,'  rlamb',rlamb
      write(6,*)'dissve=',dissve,'   disssp=',disssp
c
c  evaluation of the vorticity spectra
c
      call spectre(ru,qcap)
      filcos='spesc.'//pntim
      open(57,file=filcos)
      rewind(57)
      filcos='spesko.'//pntim
      open(59,file=filcos)
      rewind(59)
c     speska=(cvisc**5./diss)**(1./4.)
       es1t0=0.
       es2t0=0.
       es3t0=0.
       do k=1,kkmax
       write (57,*)k,e(1,k),e(2,k),e(3,k)
       es1to=es1to+e(1,k)
       es2to=es2to+e(2,k)
       es3to=es3to+e(3,k)
       akol=k*eta
       se1=e(1,k)/speska
       se2=e(2,k)/speska
       se3=e(3,k)/speska
       write (59,*)akol,se1,se2,se3          
       end do
      close(57)
      close(59)
      esto=es1to+es2to+es3to
      write(6,134)time,es1to,es2to,es3to,esto,diss
  134 format('enstrophy',e12.4,3x,3e12.4,3x,'diss sp,phy',2e12.5)
             endif
      filcos='corvel.'//pntim
      open(57,file=filcos)
      rewind(57)
      filcos='corvor.'//pntim
      open(59,file=filcos)
      rewind(59)
c
c   calculation of velocity correlations
c
      call correl(q,qcap)
      n3mh=n3m/2+1
      do n=1,3
      do k=1,n3mh
      ruu(n,k)=corf(n,k)/corf(n,1)
      enddo
      enddo
c
c   calculation of vorticity correlations
c
      call correl(ru,qcap)
      do n=1,3
      do k=1,n3mh
      roo(n,k)=corf(n,k)/corf(n,1)
      enddo
      enddo
       do k=1,n3mh
       kk=k+n3mh-1
       yd=yp3(kk)
       write (57,*)yd,ruu(1,k),ruu(2,k),ruu(3,k)
       write (59,*)yd,roo(1,k),roo(2,k),roo(3,k)
       enddo
      close(57)
      close(59)
      qmax=qquma
      write(6,*)' enter pdf qmax=',qmax
c
c    PDF of velocity u_1
c
      qnorm=vrms(1)
      call pdfpdf(q,qcap,1,qnorm,qmax,0)
c
c    PDF of velocity u_2
c
      qnorm=vrms(2)
      call pdfpdf(q,qcap,2,qnorm,qmax,0)
c
c    PDF of velocity u_3
c
      qnorm=vrms(3)
      call pdfpdf(q,qcap,3,qnorm,qmax,0)
      qmax=qqoma
c
c    PDF of vorticity om_1
c
      qnorm=vorms(1)
      call pdfpdf(ru,qcap,1,qnorm,qmax,3)
c
c    PDF of vorticity om_2
c
      qnorm=vorms(2)
      call pdfpdf(ru,qcap,2,qnorm,qmax,3)
c
c    PDF of vorticity om_3
c
      qnorm=vorms(3)
      call pdfpdf(ru,qcap,3,qnorm,qmax,3)
c
c    PDF of pressure pr   
c
      qmax=qqoma
      n=7
      qnorm=sqrt(prms)
      ntqll(n)=0
      ntqou(n)=0
      nouq=0
      do ll=1,lipdf
      nllq(ll)=0
      enddo
      call pdfqua(pr,qnorm,nllq,nouq,qmax)
      nouqu(n)=nouq
      ntqou(n)=ntqou(n)+nouq
      do ll=1,lipdf
      nllqu(n,ll)=nllq(ll)
      ntqll(n)=ntqll(n)+nllq(ll)
      enddo
      write(18,121)n,ntqll(n),ntqou(n)
  121 format(3x,'in outp press',3x,5i10)
      do ll=1,lipdf
      pdfqu(n,ll)=nllqu(n,ll)/float(ntqll(n))
c     write(18,123)ll,nllqu(n,ll),pdfqu(n,ll)
  123 format(3x,i4,3x,i10,3x,e12.5)
      enddo
c
c    PDF of helicity density  
c
      do n=1,6
      ntcoll(n)=0
      enddo
      n=1
      qnorm=1.
      qmax=1.
      ntsll(n)=0
      ntsou(n)=0
      nousc(1)=0
      do ll=1,lipdf
      nllsc(1,ll)=0
      do n=1,6
      nllco(n,ll)=0
      enddo
      enddo
      helmed=0.
      vla1me=0.
      vla2me=0.
      vla3me=0.
      helrms=0.
      vla1rm=0.
      vla2rm=0.
      vla3rm=0.
      vl123=1./float(n1m*n2m*n3m)
      do k=1,n3m
      kp=kpv(k)
      do j=1,n2m
      jp=jpv(j)
      do i=1,n1m
      ip=ipv(i)
      v1p=(q(1,ip,j,k)+q(1,i,j,k))*0.5-vm(1)
      v2p=(q(2,i,jp,k)+q(2,i,j,k))*0.5-vm(2)
      v3p=(q(3,i,j,kp)+q(3,i,j,k))*0.5-vm(3)
      o1p=(ru(1,i,jp,kp)+ru(1,i,j,kp)
     1    +ru(1,i,jp,k)+ru(1,i,j,k))*0.25-vo(1)
      o2p=(ru(2,ip,j,kp)+ru(2,i,j,kp)
     1    +ru(2,ip,j,k)+ru(2,i,j,k))*0.25-vo(2)
      o3p=(ru(3,ip,jp,k)+ru(3,i,jp,k)
     1    +ru(3,ip,j,k)+ru(3,i,j,k))*0.25-vo(3)
      voplo=sqrt(o1p**2+o2p**2+o3p**2)
      veplo=sqrt(v1p**2+v2p**2+v3p**2)
      denhe=veplo*voplo
      dhel1=v1p*o1p/denhe
      dhel2=v2p*o2p/denhe
      dhel3=v3p*o3p/denhe
      hede=dhel1+dhel2+dhel3
      helmed=helmed+hede
      helrms=helrms+hede**2
c
c   total helicity density
c
      all=lipdh*(qmax+hede)/qmax+1.5
      ll=all
      nllsc(1,ll)=nllsc(1,ll)+1
c
c   total helicity density componenents
c
      all=lipdh*(qmax+dhel1)/qmax+1.5
      ll=all
      nllco(1,ll)=nllco(1,ll)+1
      all=lipdh*(qmax+dhel2)/qmax+1.5
      ll=all
      nllco(2,ll)=nllco(2,ll)+1
      all=lipdh*(qmax+dhel3)/qmax+1.5
      ll=all
      nllco(3,ll)=nllco(3,ll)+1
c
c   componenents of Lamb vector V X omega    
c
      vela1=(o3p*v2p-o2p*v3p)/denhe
      vela2=(o1p*v3p-o3p*v1p)/denhe
      vela3=(o2p*v1p-o1p*v2p)/denhe
      vla1me=vla1me+vela1
      vla2me=vla2me+vela2
      vla3me=vla3me+vela3
      vla1rm=vla1rm+vela1**2
      vla2rm=vla2rm+vela2**2
      vla3rm=vla3rm+vela3**2
      all=lipdh*(qmax+vela1)/qmax+1.5
      ll=all
      nllco(4,ll)=nllco(4,ll)+1
      all=lipdh*(qmax+vela2)/qmax+1.5
      ll=all
      nllco(5,ll)=nllco(5,ll)+1
      all=lipdh*(qmax+vela3)/qmax+1.5
      ll=all
      nllco(6,ll)=nllco(6,ll)+1
      enddo
      enddo
      enddo
      helmed=helmed*vl123
      vla1me=vla1me*vl123
      vla2me=vla2me*vl123
      vla3me=vla3me*vl123
      helrms=helrms*vl123
      vla1rm=vla1rm*vl123
      vla2rm=vla2rm*vl123
      vla3rm=vla3rm*vl123
      write(6,*)' heli med=',helmed
      write(6,*)' heli rms=',helrms
      write(6,*)' lam med=',vla1me,vla2me,vla3me
      write(6,*)' lam rms=',vla1rm,vla2rm,vla3rm
      ntsou(1)=ntsou(1)+nousc(1)
      do ll=1,lipdf
      ntsll(1)=ntsll(1)+nllsc(1,ll)
      do n=1,6
      ntcoll(n)=ntcoll(n)+nllco(n,ll)
      enddo
      enddo
      write(18,141)ntsll(1),(ntcoll(n),n=1,6)
  141 format(3x,'heli',1x,i10,2x,' components ',6i10)
      do ll=1,lipdf
      pdfsc(1,ll)=nllsc(1,ll)/float(ntsll(1))
c     write(18,123)ll,nllsc(n,ll),pdfsc(n,ll)
      do n=1,6
      pdfco(n,ll)=nllco(n,ll)/float(ntcoll(n))
      enddo
      enddo
c
c    PDF of dissipation  local
c
      n=2
      qnorm=diss
      qmax=qqoma
      ntsll(n)=0
      ntsou(n)=0
      nouq=0
      do ll=1,lipdf
      nllq(ll)=0
      enddo
      call pdfqlo(sca,qnorm,nllq,nouq,qmax)
      nousc(n)=nouq
      ntsou(n)=ntsou(n)+nouq
      do ll=1,lipdf
      nllsc(n,ll)=nllq(ll)
      ntsll(n)=ntsll(n)+nllq(ll)
      enddo
      write(18,143)n,ntsll(n),ntsou(n)
  143 format(3x,'in outp diss',3x,5i10)
      do ll=1,lipdf
      pdfsc(n,ll)=nllsc(n,ll)/float(ntsll(n))
c     write(18,123)ll,nllsc(n,ll),pdfsc(n,ll)
      enddo
      call pdffiq(time)
c
c    pdf  of structure function
c
      do n=1,6
           do lr=1,nlr,nlrju
      ntsfl(n,lr)=0
      ntsfu(n,lr)=0
           enddo
      enddo
      if(inostf.eq.0) then
c
c  velocity struct. funct.
c
      call stfun(q,0,sfuma)
c
c  vorticity struct. funct.
c
      call stfun(ru,3,sfoma)
c
c  pressure struct. funct.
c
      call sfprn(pr,sfoma)
                      else
c
c  velocity struct. funct.
c
      call stfunn(q,0,sfuma)
c
c  vorticity struct. funct.
c
      call stfunn(ru,3,sfoma)
c
c  pressure struct. funct.
c
      call sfprnn(pr,sfoma)
                      endif
      do n=1,9
           do lr=1,nlr,nlrju
      do ll=1,lipdf
      pdfsf(n,lr,ll)=nllsf(n,lr,ll)/float(ntsfl(n,lr))
      enddo
           enddo
      enddo
c
c    the pdf are written
c
      call pdffis(time)
      return
      end
