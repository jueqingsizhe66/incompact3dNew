      parameter (m1=97,m2=193,m3=193,ndd=2,ndv=3)
      parameter (m1m=m1-1,m2m=m2-1,m3m=m3-1)
      parameter (m3md=m3m+2,m3mp=m3m+4)
      parameter (m12=2*m1m,m32=2*m3m)
      parameter (m1p=m1m+2,msc=1)
c************************************************************************
      common/qqq/q1(m1,m2,m3),q2(m1,m2,m3),q3(m1,m2,m3)
      common/voc/dq(m1,m2,m3),dph(m1,m2,m3),qcap(m1,m2,m3)
      common/coo/yp1(m1),yp2(m2),yp3(m3)
      common/ppr/pr(m1,m2,m3)
c************************************************************************
      common/passc/psc(m1,m2,m3)
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
      common/cordoo/rco(m2),rmo(m2),g2rco(m2),g2rmo(m2)
      common/corrt/thetac(m1),ragc(m2),thetam(m1),ragm(m2)
      common/ctrdph/amph(m2),acph(m2),apph(m2)
      common/ctrdph2/amphp(m2),acphp(m2),apphp(m2)
      common/util/bet(m1),x(m1),gm(m1,m2)
      common/d1/re,tfin,eps
      common/d13/alx3
      common/pigr/pi
      common/timp/timew
      common/d2/nstop,nprint,ntst,npin,npstf,ireset
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/dimens/alx3d
      common/npjet/n2t,n2v
      common/nonunz/str3,rmed31,istr3,etdp3,strb3
      common/stfft/coecos(m3),coesin(m3)
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/inqca/imxq,jmxq,kmxq
      common/inwal/jpc(m2),jup(m2),jmc(m2),jum(m2)
      common/invoma/iazm,jazm,kazm,irrm,jrrm,krrm
     1             ,irzm,jrzm,krzm
      common/jump/jum1,jum2,jum3,jum4
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      common/namefi/filcnw,filcnr,filth,filvm,filpo
      common/njump/n1p,n2p,n3p
      common/numfil/nfil
      common/parcoo/r0
      common/posma/radm(m1),zetm(m1)
      common/posmax/jmaxv(ndv),kmaxv(ndv)
      common/pscpar/pran,pscwal
      common/rhsc/rhs(m1,m2,m3)
      common/rota/ros
      common/section/nini,nfin,nstri,irid
      common/nonunif/strr,rext,rint,rmed,rmed1,istr,etdp,strb
      common/nounifo/strro,rexto,rinto,rmedo,rmed1o,istro,etdpo,strbo
      common/npjeto/n2to
      common/olddim/n1l,n2l,n3l
      common/tstep/dt,beta,ren
      common/tscoe/gam(3),rom(3),alm(3),nsst
      common/velmax/vmax(ndv),vomax(ndv),vomin(ndv)
      common/vperin/vper,epsil,lamb
c     common/frequi/alamk,alamk2,alamk3,alamk4
      common/vore/omamax,omamin,omrmax,omrmin,omzmax,omzmin
      common/wrre/nwrit,nread
      common/waves/an(m3),ap(m1),ak3(m3),ak1(m1)
      common/wavesp/anp(m3),app(m1),ak3p(m3),ak1p(m1)
      common/nprdel/nprde
      common/velaxi/vax(m2),vaxk(m1,m2)
      common/zewa/waze,iaxsy
c
c   time dependent inflow cond.
c
      common/tinfl/uinf(m2),deinf(m2)
      common/tprfim/tprfi
      common/vinfl/strh,vpin
      common/inflth/uinfth(m1,m2),deinth(m1,m2)
      common/lesin/ialvo,icorspe,npouth
c
c   quantities in coetar for hdnl and invtr
c

      common/coefj/upd1(m2),upd2(m2),udx1q(m2),udxi(m2),udx1(m2)
     1            ,udx2(m2),vdx1(m2),vd1d2(m2),vd1d3(m2),vdxi(m2)
     1            ,udx1c(m2),udx2c(m2),vdyh2(m2),vdyh3(m2),ugmv(m2)
     1            ,udh2(m2),vh23(m2),vh13(m2),a11(m2),uvdx1(m2)
     1            ,vdxh2(m2),volz(m2),udx1n(m2)
      common/indax/isym(m1)
      common /inslip/ inslws,inslwn,inslwr
      common /cor1j/  ap1j(m2),ac1j(m2),am1j(m2)
      common /cor2j/  ap2j(m2),ac2j(m2),am2j(m2)
      common /cor2je/ ap2je(m2),ac2je(m2),am2je(m2)
      common /cor3j/  ap3j(m2),ac3j(m2),am3j(m2)
      common /corscj/ apscj(m2),acscj(m2),amscj(m2)
      common /bcdens/ denbs(m1,m2),denbn(m1,m2)
c
c     quantities for pipe
c
      common/it0q1/ime1t0,ipr0,ipsc0
      common/vmean/vmed(3,m3,m2),turstr(6,m3,m2),pmed(m3,m2)
      common/vomean/vorv(3,m3,m2),vorstr(6,m3,m2)
      common/sctuo/pscmed(m3,m2),enpsc(m3,m2)
      common/quamp/psimed(m3,m2),vmep(3,m3,m2),vormed(m3,m2)

      common/axv/q3ax(m3)
      common/veltot/vit(5)
      common/dismax/dissma,dissmn
      common/eneav/enav,enavo,cfn,cfo,dp3mo,disst,dissto
      common/quath/dp3th,qsouth
      common/iavg/iav
      common/tauwal/cfnp
      common/presgr/dp3ns
      common/qflux/qsourc,qsouro
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
      common/i3dpr/i3dou,njumk,njumj,ioldf,jprq(15),n3mol,kprq(15)

