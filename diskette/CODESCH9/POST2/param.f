      parameter (m1=65,m2=151,m3=257,ndd=2,ndv=3)
      parameter (m1m=m1-1,m2m=m2-1,m3m=m3-1)
      parameter (m3md=m3m+2,m3mp=m3m+4)
      parameter (m12=2*m1m,m32=2*m3m)
      parameter (m1p=m1m+2,msc=1)
      parameter (mpq=13)
c************************************************************************
      common/qqq/q1(m1,m2,m3),q2(m1,m2,m3),q3(m1,m2,m3)
      common/coo/yp1(m1),yp2(m2),yp3(m3)
      common/ppr/pr(m1,m2,m3)
c************************************************************************
      character*4 ipre2,ipfi
      character*4 ipre
      character*4 dummy, ri3d
      character*60 namfi3
      character*60 namfile
      character*60 filcnw,filcnr,filth,filvm,filpo,filen
      character*60 filet,filer,filez,filed,filev
      common/area/a12
      common/cordvo/y(m2),y2s(m2),cac(m2),caj(m2),y1s(m1),y3s(m3)
      common/jprts/ipri(5),kpri(5)
      common/ctrdph/amph(m2),acph(m2),apph(m2)
      common/ctrdph2/amphp(m2),acphp(m2),apphp(m2)
      common/util/bet(m1),xxx(m1),gm(m1,m2)
      common/d1/re,tfin,eps
      common/d13/alx3,alx1
      common/pigr/pi
      common/timp/timew
      common/d2/nstop,nprint,ntst,npin,npstf,ireset
      common/dim/n1,n1m,n2,n2m,n3,n3m
      common/dimens/alx3d,alx1d
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
      common/posmax/jmaxv(ndv),kmaxv(ndv)
      common/rhsc/rhs(m1,m2,m3)
      common/section/nini,nfin,nstri,irid
      common/nonunif/strr,rext,rint,rmed,rmed1,istr,etdp,strb
      common/olddim/n1l,n2l,n3l
      common/tstep/dt,beta,ren
      common/velmax/vmax(ndv),vmaxo(ndv)
      common/vperin/vper,epsil,lamb
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
      common/tprfim/tprfi
      common/vinfl/strh,vpin
      common/lesin/itot,icorspe,icorr,npouth
c
c   quantities in coetar for hdnl and invtr
c

      common/coefj/upd1(m2),upd2(m2),udx1q(m2),udxi(m2),udx1(m2)
     1            ,udx2(m2),vdx1(m2),vd1d2(m2),vd1d3(m2),vdxi(m2)
     1            ,udx1c(m2),udx2c(m2),vdyh2(m2),vdyh3(m2),ugmv(m2)
     1            ,udh2(m2),vh23(m2),vh13(m2),a11(m2),uvdx1(m2)
     1            ,vdxh2(m2),volz(m2),udx1n(m2)
c
c     quantities for pipe
c
      common/it0q1/ime1t0,ipr0,imhd0
      common/vmean/vmed(3,m2),turstr(6,m2),pmed(m2)
      common/vmeao/vmeo(3,m2),tursto(6,m2),pmeo(m2)
      common/vmpos/vmepo(3,m2),vompo(3,m2)
      common/vapos/vpmeo(3,m2),vpome(3,m2)
      common/vcclwa/vocl(3,m1,m3),vowa(3,m1,m3)
      common/axv/q3ax(m3)
      common/veltot/vit(5)
      common/eneav/enav,enavo,cfn,cfo,dp3mo,disst,dissto
      common/quath/dp3th,qsouth
      common/iavg/iav
      common/skfl/ske(4,m2),fla(4,m2),pvc(4,m2)
      common/skflo/skeo(4,m2),flao(4,m2),pvmo(4,m2)
      common/skflqo/skequo(12,m2),flaquo(12,m2),sququo(12,m2)
      common/skflq/skequ(12,m2),flaqu(12,m2),sququ(12,m2)
      common/skflqq/quaav(12,m2),quasq(12,m2),quacu(12,m2),quafo(12,m2)
      common/tauwal/cfnp
      common/vomean/vorv(3,m2),vorstr(6,m2),vdtomt(3,m2),vcromt(6,m2)
      common/vomeao/voro(3,m2),vorsto(6,m2),vdtomo(3,m2),vcromo(6,m2)
      common/prmeao/prgro(3,m2)
      common/prmean/prgrd(3,m2)
      common/vomeao/vdtoma(3,m2),vcroma(6,m2),pgrada(3,m2)
      common/diss/dissj(m2),stmed(6,m2),enspr(m2),ensy(m2)
      common/disso/dissjo(m2),stmedo(6,m2),enspo(m2),ensyo(m2)
      common/vtrv/vtv(5,m2),vtvo(5,m2),dvto(5,m2)
      common/presgr/dp3ns
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

      common/cospez/en12zo(m2,m3m),en13zo(m2,m3m),en23zo(m2,m3m)
     1   ,env12z(m2,m3m),env13z(m2,m3m),env23z(m2,m3m)
     1   ,evo11z(m2,m3m),evo12z(m2,m3m),evo13z(m2,m3m)
     1   ,evo21z(m2,m3m),evo22z(m2,m3m),evo23z(m2,m3m)
     1   ,evo31z(m2,m3m),evo32z(m2,m3m),evo33z(m2,m3m)
      common/specz/ene1zo(m2,m3m),ene2zo(m2,m3m),ene3zo(m2,m3m)
     1   ,enepzo(m2,m3m)
     1   ,env33z(m2,m3m),env22z(m2,m3m),env11z(m2,m3m)
      common/spect/ene1to(m2,m1m),ene2to(m2,m1m),ene3to(m2,m1m)
     1   ,enepto(m2,m1m)
     1   ,env33t(m2,m1m),env22t(m2,m1m),env11t(m2,m1m)
      common/cospet/en12to(m2,m1m),en13to(m2,m1m),en23to(m2,m1m)
     1   ,env12t(m2,m1m),env13t(m2,m1m),env23t(m2,m1m)
     1   ,evo11t(m2,m1m),evo12t(m2,m1m),evo13t(m2,m1m)
     1   ,evo21t(m2,m1m),evo22t(m2,m1m),evo23t(m2,m1m)
     1   ,evo31t(m2,m1m),evo32t(m2,m1m),evo33t(m2,m1m)
      common/wavn/apik2(m1m),ankk2(m3m)
      common/kma/k1max,k3max
      common/enec/ene1(m2,m1m),ene3(m2,m3m)
      common/ent/e1to(m2),e2to(m2),e3to(m2)
     1   ,e1stto(m2),e2stto(m2),e3stto(m2)
     1   ,e1stzo(m2),e2stzo(m2),e3stzo(m2)
      common/cor3/corx3(m2,m3)
      common/cor1/corx1(m2,m1)
      common/corx3/r11x3(m2,m3m),r22x3(m2,m3m),r33x3(m2,m3m)
     1            ,r12x3(m2,m3m),r23x3(m2,m3m),r31x3(m2,m3m)
      common/covx3/v11x3(m2,m3m),v22x3(m2,m3m),v33x3(m2,m3m)
     1            ,v12x3(m2,m3m),v23x3(m2,m3m),v31x3(m2,m3m)
      common/covox3/vo11x3(m2,m3m),vo12x3(m2,m3m),vo13x3(m2,m3m)
     1             ,vo21x3(m2,m3m),vo22x3(m2,m3m),vo23x3(m2,m3m)
     1             ,vo31x3(m2,m3m),vo32x3(m2,m3m),vo33x3(m2,m3m)
      common/corx1/r11x1(m2,m1m),r22x1(m2,m1m),r33x1(m2,m1m)
     1            ,r12x1(m2,m1m),r23x1(m2,m1m),r31x1(m2,m1m)
      common/covx1/v11x1(m2,m1m),v22x1(m2,m1m),v33x1(m2,m1m)
     1            ,v12x1(m2,m1m),v23x1(m2,m1m),v31x1(m2,m1m)
      common/covox1/vo11x1(m2,m1m),vo12x1(m2,m1m),vo13x1(m2,m1m)
     1             ,vo21x1(m2,m1m),vo22x1(m2,m1m),vo23x1(m2,m1m)
     1             ,vo31x1(m2,m1m),vo32x1(m2,m1m),vo33x1(m2,m1m)
      common/copx3/pcox3(m2,m3m)
      common/copx1/pcox1(m2,m1m)
      common/ruut/ruuth(m2,m1m)
      common/ruux/ruuax(m2,m3m)
      common/i3dpr/i3dou,njumk,njumj,ioldf,jprq(15)
c
      common/dismax/dissma(mpq)
      common/ensmax/enssma(mpq),esprma(mpq)
      common/espmm/espmi(mpq),espma(mpq)
      common/nfpq/npqf
      common/x3metria/caj3(m3),cac3(m3)





