c
c  
c   In this code the energy spectra from the simulation 
c   are read and are rewritten in kolmogorov scaling
c
      program main  
c
      include 'param.f'
      common/speene/e(ndv,0:m3)
      dimension tl(50),enl(50),esl(50),rla(50)
      dimension tt(50),enn(50),ens(50)
      character*17 filcos
      character*4 pntim
      pi=2.*asin(1.)
      filcos='post.d'
      open(5,file=filcos)
c     write(6,*)' enter ispew=1 to write the spectra'
      read(5,*)ispew
c     write(6,*) 'enter Re,nn'
      read(5,*) re,nn,f0
      if(ispew.eq.1) then
      write(6,*) 'enter itime'
      read(5,*) itim
      ntfi=1
                     else
c     write(6,*) 'enter ntfi,dti,timei'
      read(5,*) ntfi,dti,timei
      filcos='b33vv.out'
      open(45,file=filcos)
      filcos='b33oo.out'
      open(44,file=filcos)
      filcos='b33ooro.out'
      open(43,file=filcos)
      filcos='rvv3.out'
      open(54,file=filcos)
      filcos='roo3.out'
      open(55,file=filcos)
      filcos='lescal.out'
      open(75,file=filcos)
      filcos='anirms.out'
      open(74,file=filcos)
      filcos='anivorms.out'
      open(72,file=filcos)
      filcos='vormsk.out'
      open(71,file=filcos)
      filcos='rmsk.out'
      open(73,file=filcos)
      filcos='slorla.out'
      open(63,file=filcos)
      filcos='ckslen.out'
      open(64,file=filcos)
      filcos='sloes.out'
      open(62,file=filcos)
      filcos='quato.out'
      open(61,file=filcos)
      filcos='slofo.out'
      open(79,file=filcos)
      filcos='sloen.out'
      open(92,file=filcos)
      filcos='c2esc.out'
      open(93,file=filcos)
      filcos='c2rla.out'
      open(95,file=filcos)
      filcos='ene.out'
      open(25,file=filcos)
      filcos='altu.out'
      open(23,file=filcos)
      filcos='tsca.out'
      open(24,file=filcos)
      filcos='ens.out'
      open(26,file=filcos)
      filcos='enetau.out'
      open(27,file=filcos)
      filcos='enstau.out'
      open(28,file=filcos)
       if(f0.gt.0.) then
      filcos='rolato.out'
      open(65,file=filcos)
      filcos='rosoto.out'
      open(66,file=filcos)
                    endif
                     endif
      cvisc=1./re
      kkmax=nn/2
      deltx=2.*pi/float(nn)
  202 continue
      do nnt=1,ntfi
      if(ispew.eq.0) itim=(nnt-1)*dti+timei
      time=itim
      write(pntim,77) itim
   77 format(i4.4)
      filcos='../spec.'//pntim
      open(57,file=filcos)
      rewind(57)
c     write(6,*)'read from',filcos
       ee1to=0.
       ee2to=0.
       ee3to=0.
       en1to=0.
       en2to=0.
       en3to=0.
       al1=0.
       al2=0.
       al3=0.
       alto=0.
       do k=1,512
       read (57,*,end=121)lk,e(1,k),e(2,k),e(3,k)
       ee1to=ee1to+e(1,k)
       ee2to=ee2to+e(2,k)
       ee3to=ee3to+e(3,k)
       en1to=en1to+e(1,k)*k**2
       en2to=en2to+e(2,k)*k**2
       en3to=en3to+e(3,k)*k**2
       alto=alto+(e(1,k)+e(2,k)+e(3,k))/k
       al1=al1+e(1,k)/k
       al2=al2+e(2,k)/k
       al3=al3+e(3,k)/k
       end do
  121  continue
      close(57)
       eeto=ee1to+ee2to+ee3to
       ento=en1to+en2to+en3to
       diss=2.*cvisc*ento
       tken=eeto
       rlam=tken*sqrt(20./(3.*cvisc*diss))
       eta=(cvisc**3/diss)**(1./4.)
       speska=(cvisc**5.*diss)**(1./4.)
       upri=4.*eeto/3.
       all=2.*asin(1.)*alto/tken
       tle=all/sqrt(tken)
       ala=sqrt(20.*tken*cvisc/diss)
       write(6,133)time,tken,diss,eta,speska,rlam,all,ala,tle
       eps=.1e-05
       tl(nnt)=alog(time+eps)
       tt(nnt)=time
       enn(nnt)=tken
       ens(nnt)=diss
       enl(nnt)=alog(tken)
       esl(nnt)=alog(diss)
       rla(nnt)=rlam
       ralle=all/eta
       ralae=ala/eta
       rdelx=deltx/eta
       write(61,133)time,tken,diss,eta,speska,rlam,all,ralle,ralae
     1             ,rdelx
       if(ispew.eq.0) then
       ru3=0.5*(ee1to+ee2to)/ee3to
       ro3=0.5*(en1to+en2to)/en3to
       all1=2.*asin(1.)*al1/ee1to
       all2=2.*asin(1.)*al2/ee2to
       all3=2.*asin(1.)*al3/ee3to
       ru11=ee1to/eeto-1./3.
       ru22=ee2to/eeto-1./3.
       ru33=ee3to/eeto-1./3.
       ro11=en1to/ento-1./3.
       ro22=en2to/ento-1./3.
       ro33=en3to/ento-1./3.
      tau=ens(1)/enn(1)*(tt(nnt)-tt(1))
      write(54,133)tau,ru3
      write(55,133)tau,ro3
      write(75,133)time,all1,all2,all3
      write(45,133)tau,ru33
      write(44,133)tau,ro33
      write(74,133)time,ru11,ru22,ru33
      write(72,133)time,ro11,ro22,ro33
      write(71,133)time,en1to,en2to,en3to
      write(73,133)time,ee1to,ee2to,ee3to
       if(f0.gt.0.) then
       rosla=2./f0*re/rlam
       rosvo=sqrt(ento)*2./f0
       rosll=sqrt(ee3to)/all3/f0
      write(66,133)time,rosvo
      write(65,133)time,rosla
      write(43,133)rosll,ro33
                     endif
                     endif
  133  format(3x,10e15.6)
      if(ispew.eq.1) then
      filcos='spek1'//pntim//'.plo'
      open(11,file=filcos)
      filcos='spek2'//pntim//'.plo'
      open(12,file=filcos)
      filcos='spek3'//pntim//'.plo'
      open(13,file=filcos)
      filcos='speko'//pntim//'.plo'
      open(59,file=filcos)
      filcos='speto'//pntim//'.plo'
      open(57,file=filcos)
       do k=1,kkmax
       akol=k*eta
       ak=k
       se1=e(1,k)/speska
       se2=e(2,k)/speska
       se3=e(3,k)/speska
       write (59,*)akol,se1+se2+se3          
       write (11,*)akol,se1
       write (12,*)akol,se2
       write (13,*)akol,se3
       write (57,*)ak,e(1,k)+e(2,k)+e(3,k)
       end do
      close(59)
      close(57)
      write(6,*)' enter imor=1 for a new time'
      read(5,*)imor
      if(imor.eq.1) go to 202
                          endif
       enddo
      if(ispew.eq.0) then
      do nnt=1,ntfi
      tdif=(tt(nnt)-tt(1))
      if(f0.eq.0)  then
      tau=ens(1)/enn(1)*(tt(nnt)-tt(1))
                   else
      tau=ens(1)/enn(1)*(tt(nnt)-tt(1))
                   endif
      tsca=enn(nnt)/ens(nnt)*ens(1)/enn(1)
      altu=enn(nnt)**(3./2.)/ens(nnt)
      ren=enn(1)/enn(nnt)
      rdi=ens(1)/ens(nnt)
      write(27,133)tau,ren
      write(28,133)tau,rdi
      write(25,133)tdif,ren
      write(26,133)tdif,rdi
      write(24,133)tau,tsca
      write(23,133)tau,altu
      enddo
      nth=(ntfi-1)/2
      do nnt=2,ntfi-1
      dnen=-(enl(nnt-1)-enl(nnt+1))/(tl(nnt-1)-tl(nnt+1))
      tau=ens(nnt)/enn(nnt)*tt(nnt)
      dnes=-(esl(nnt-1)-esl(nnt+1))/(tl(nnt-1)-tl(nnt+1))
      des=-(ens(nnt-1)-ens(nnt+1))/(tt(nnt-1)-tt(nnt+1))
      desc=dnes/tt(nnt)*ens(nnt)
      if(f0.eq.0)  then
      c2esc=dnes/tau
      c2en=(dnen+1)/dnen
      c2es=dnes/tt(nnt)*enn(nnt)/ens(nnt)
      write(62,*)tt(nnt),dnes,des,c2en,c2es,c2esc
      write(93,*)tt(nnt),c2esc
      write(95,*)rla(nnt),c2esc
                   else
      if(nnt.eq.nth.or.nnt.eq.ntfi-1) then
      c2es=1.70
      cres=-(-dnes+c2es*tau)*2./f0/tt(nnt)
      write(62,*)f0,cres
                                    endif
                   endif
      write(63,*)rla(nnt),dnen
      write(64,*)tt(nnt),dnen,tau
      write(92,*)tt(nnt),dnen
      if(nnt.eq.nth.or.nnt.eq.ntfi-1) then
      write(79,*)f0,dnen
                                    endif
      write(6,*)rla(nnt),dnen,dnes,c2en,c2es,c2esc,des,desc
      enddo
                     endif
      stop
      end
