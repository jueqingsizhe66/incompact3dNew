      include'paramdi.f'
      parameter (m1=mm1,m2=mm2,m3=mm3)
      parameter (ndd=2,ndv=3)
      parameter (m1m=m1-1,m2m=m2-1,m3m=m3-1)
      parameter (m3c=2*m3m,m3p=m3c+2)
      parameter (m12=3*m1m/2+1,m32=2*m3c)
c************************************************************************
      common/qqq/q1(m1,m2,m3),q2(m1,m2,m3),q3(m1,m2,m3)
      common/rrr/ru1(m1,m2,m3),ru2(m1,m2,m3),ru3(m1,m2,m3)
      common/ddq/dq(m1,m2,m3)
      common/coo/yp1(m1),yp2(m2),yp3(m3)
      common/ppr/pr(m1,m2,m3),qcap(m1,m2,m3),dph(m1,m2,m3)
c************************************************************************
      common/passc/psc(m1,m2,m3),rupsc(m1,m2,m3)
      common /bcpsc/ pscbs(m1,m2),pscbn(m1,m2)
      character*4 ipre2,ipfi
      character*4 ipre
      character*4 dummy, ri3d
      character*60 namfi3
      character*60 namfile
      character*60 filcnw,filcnr,filth,filvm,filpo,filen
      character*60 filet,filer,filez,filed,filev
      common/area/a12
      common/cordv2/rc(m2),rm(m2),g2rc(m2),g2rm(m2)
      common/cordv3/zz(m3),zzm(m3),g3rc(m3),g3rm(m3)
      common/corrt/thetac(m1),ragc(m2),thetam(m1),ragm(m2)
      common/ctrdph/amph(m2),acph(m2),apph(m2)
      common/util/x(m1)
      common/d1/re,tfin,eps
      common/d13/alx3
      common/pigr/pi
      common/timp/timew
      common/d2/nstop,nprint,ntst,npin,npstf,ireset
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/dimens/alx3d
      common/nonunz/str3,rmed31,istr3,etdp3,strb3
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common /inwa2/  jpc(m2),jmc(m2)
      common /inwa3/  kpc(m3),kmc(m3)
      common/inqca/imxq,jmxq,kmxq
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      common/namefi/filcnw,filcnr,filth,filvm,filpo
      common/njump/n1p,n2p,n3p
      common/posma/radm(m1),zetm(m1)
      common/pscpar/pran,pscwal
      common/rhsc/rhs(m1*m2*m3)
      common/rota/ros
      common/nonunif/strr,rext,rint,rmed,rmed1,istr,etdp,strb
      common/npjeto/n2to
      common/olddim/n1l,n2l,n3l
      common/tstep/dt,beta,ren
      common/tscoe/gam(3),rom(3),alm(3),nsst
      common/velmax/vmax(ndv),vomax(ndv)
      common/velmin/vmin(ndv),vomin(ndv)
      common/wrre/nwrit,nread
      common/waves/an(m3),ap(m1),ak3(m3),ak1(m1)
c
c   time dependent inflow cond.
c

      common/indax/isym(m1)
      common /inslip/ inslws,inslwn,inslwr

c
c     quantities for pipe
c
      common/vmean/vmed(3,m3,m2),turstr(6,m3,m2),pmed(m3,m2)
      common/vomean/vorv(3,m3,m2),vorstr(6,m3,m2)
      common/qflux/qsourc,qsouro
      common/sctuo/pscmed(m3,m2),enpsc(m3,m2)
      common/enpsav/enpsv,enpsvo,disstps,disstpo
c
c     dt variabile
c
      common/cfcost/icfl,cflc

      character*60 filqm,filve,filvo,filvn,filvt,filvmp
      character*60 filns,filts,filsk,filfl,filuvp,filnsp
      character*60 filvsp,filnu,filba1,filba2,filba3
      character*60 filvtb,filval,filpsc,filpco,filpsp
      character*4 pntim
      character*2 jsp
      character*1 jzo
      character*3 jrp
c
c
c   spectra and correlations
c

      common/wavn/apik2(m1m),ankk2(m3m)
      common/kma/k1max,k3max

      common /cor3ck/  ap3ck(m3),ac3ck(m3),am3ck(m3)
      common /cor3sk/  ap3sk(m3),ac3sk(m3),am3sk(m3)
      common /cor3ss/  ap3ssk(m3),ac3ssk(m3),am3ssk(m3)
      common /cor1j/  ap1j(m2),ac1j(m2),am1j(m2)
      common /cor2j/  ap2j(m2),ac2j(m2),am2j(m2)
      common /cor2je/ ap2je(m2),ac2je(m2),am2je(m2)
      common /cor3j/  ap3j(m2),ac3j(m2),am3j(m2)
      common /corscj/ apscj(m2),acscj(m2),amscj(m2)
      common /bcdens/ denbs(m1,m2),denbn(m1,m2)
      common /stfft/  coecos(m3),coesin(m3)
      common /trici/  ami(m2),aci(m2),api(m2),
     %                fei(m2,m1),q(m2,m1),s(m2,m1)
      common /tricj/  amj(m2),acj(m2),apj(m2)
      common /trick/  amk(m3),ack(m3),apk(m3)
      common /cft/    nx3fft
      common /fftvel/ ifxv(13),trigxv(3*m1m/2+1)
      common/phcoj/amphj(m2),
     1   acphj(m1,m2),apphj(m2),qsbph(m1,m2),fphj(m1,m2)
      common/phcok/amphk(m3),acphk(m3),apphk(m3)
      common/phcojv/amphjv(m2),acphjv(m1,m2),apphjv(m2)
      common/phcokv/amphkv(m3),acphkv(m3),apphkv(m3)
      common/quafi/yfis(m2,m3),w(m2*m3)
            common/cofisn/mw,np,mp
      common/timavg/timav
      common/npjet/n2t,n2v
      common/pert/amp,alpha,crad,amp3
      common/nwav/nwa,nwa3,iran1,iran3,sigr,rap,rapsc
      common/inipr/nprde
      common/riduz/irid
      common/typvor/iring,itrip
      common/pricfl/tpin,tprin,tchpr,tprich
      common/vperin/vper,epsil,lamb                                
      common/parrin/vmx,sig,yc2mo,yc3mo


