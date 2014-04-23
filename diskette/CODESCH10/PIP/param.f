      parameter (m1=129,m2=129,m3=257,ndd=2,ndv=3)
      parameter (m1m=m1-1,m2m=m2-1,m3m=m3-1)
      parameter (m3md=m3m+2,m3mp=m3m+4)
      parameter (m12=2*m1m,m32=2*m3m)
      parameter (m1p=m1m+2,msc=1)
c************************************************************************
      common/qqq/q1(m1,m2,m3),q2(m1,m2,m3),q3(m1,m2,m3)
      common/rrr/ru1(m1,m2,m3),ru2(m1,m2,m3),ru3(m1,m2,m3)
      common/ddq/dq(m1,m2,m3),qcap(m1,m2,m3),dph(m1,m2,m3)
      common/coo/yp1(m1),yp2(m2),yp3(m3)
      common/ppr/pr(m1,m2,m3)
      common/dumma/dum(m1,m2,m3)
      common/rhsc/rhs(m1*m2*m3)
c************************************************************************
      common/dipovp/vor(m1,m2),psi(m1,m2)
      common/ipwlo/iwlop,ipr,kpi,kpr,njprs,npjp(10)
      character*3 ipfk,ipfj
      character*4 ipre2,ipfi
      character*4 ipre
      character*4 dummy, ri3d
      character*60 namfi3
      character*60 namfile
      character*60 filcnw,filcnr,filth,filvm,filpo,filen
      character*60 filet,filer,filez,filed,filev
      character*60 namdir
      common/namdi/namdir
      common/iflow/ipipe
      common/iwall/islip
      common/area/a12
      common/cordvo/rc(m2),rm(m2),g2rc(m2),g2rm(m2),zz(m3)
      common/cordoo/rco(m2),rmo(m2),g2rco(m2),g2rmo(m2)
      common/corrt/thetac(m1),ragc(m2),thetam(m1),ragm(m2)
      common/util/bet(m1),x(m1),gm(m1,m2)
      common/d1/re
      common/d13/alx3
      common/pigr/pi
      common/timp/timav
      common/d2/ntst,ireset
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/dimens/alx3d
      common/stfft/coecos(m3),coesin(m3)
      common/indbo/imv(m1),ipv(m1),jmv(m2),jpv(m2),kmv(m3),kpv(m3)
      common/inqca/imxq,jmxq,kmxq
      common/inwal/jpc(m2),jup(m2),jmc(m2),jum(m2)
      common/bbou/b1bou,b3bou 
      common/mesh/dx1,dx1q,dx2,dx2q,dx3,dx3q
      common/namefi/filcnw,filcnr,filth,filvm,filpo
      common/njump/n1p,n2p,n3p
      common/parcoo/r0
      common/posma/radm(m1),zetm(m1)
      common/posmax/jmaxv(ndv),kmaxv(ndv)
      common/rota/ros,pran
      common/section/nini,nfin,nstri,irid
      common/nonunif/strr,rext,rint,rmed,rmed1,istr,etdp,strb
      common/nounifo/strro,rexto,rinto,rmedo,rmed1o,istro,etdpo,strbo
      common/npjeto/n2to
      common/olddim/n1l,n2l,n3l
      common/tstep/dt,beta,ren
      common/tscoe/gam(3),rom(3),alm(3),nsst
      common/velmax/bmax(ndv),vmax(ndv)
      common/vperin/vper
      common/wrre/nwrit,nread,ipsc0
      common/waves/an(m3),ap(m1),ak3(m3),ak1(m1)
      common/wavesp/anp(m3),app(m1),ak3p(m3),ak1p(m1)
      common/phcoe/amphj(m2), 
     1   acphj(m1,m2),apphj(m2),qsbph(m1,m2),fphj(m1,m2)
      common/ctpkji/amk(m3),ack(m3),apk(m3)
      common/ctpijk/ami(m2),aci(m2),api(m2)
      common/ctrdph/amph(m2),acph(m2),apph(m2)
      common/intpr/npouth,nprint,npin
c
c   time dependent inflow cond.
c
      common/tinfl/uinf(m2)
      common/tprfim/tpin,tprin,tpouth,tfin,dtfin
      common/vinfl/strh,vpin
      common/inflth/uinfth(m1,m2)
c
c   coeff. per tridiagonals in common
c
c
c   quantities in coetar for hdnl and invtr
c

      common/coefj/upd1(m2),upd2(m2),udx1q(m2),udxi(m2),udx1(m2)
     1            ,udx2(m2),vdx1(m2),vd1d2(m2),vd1d3(m2),vdxi(m2)
     1            ,udx1c(m2),udx2c(m2),vdyh2(m2),vdyh3(m2),ugmv(m2)
     1            ,udh2(m2),vh23(m2),vh13(m2),a11(m2),uvdx1(m2)
     1            ,vdxh2(m2),volz(m2),udx1n(m2),volz2(m2)
      common/indax/isym(m1)
c
c     quantities for pipe
c
      common/vmean/vmed(3,m2),turstr(6,m2),pmed(m2),pstr(m2)
      common/vmeao/vmeo(3,m2),tursto(6,m2),pmeo(m2),psto(m2)
      common/veltot/vit(5)
      common/eneav/enav,enavo,cfn,cfo,dp3mo
      common/iavg/iav
      common/skfl/ske(4,m2),fla(4,m2)
      common/skflo/skeo(4,m2),flao(4,m2)
      common/tauwal/cfnp
      common/presgr/dp3ns
c
c     dt variabile
c
      common/cfcost/icfl,cflc
c
c   2D dipole calculation
c
      common/dippar/yc1mo,yc2mo,velmo
C 
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
c
c   tridiagonal coefficients for i direction
cc
c   coeff. per tridiagonals in common
c
      common/tricis/fi(m2,m1),fei(m2,m1),fnn(m2)
      common/tricif/amif(m2,m1,3),acif(m2,m1,3),apif(m2,m1,3),
     1   qqf(m2,m1,3),ssf(m2,m1,3),pf(m2)
      common/triciv/amiv(m2,m1,3),aciv(m2,m1,3),apiv(m2,m1,3),
     1   qqv(m2,m1,3),ssv(m2,m1,3),pv(m2)






