      parameter (m1=513)
      parameter (m1m=m1-1)
      parameter (m12=2*m1m)
      parameter (m1p=m1m+2)
c************************************************************************
      common/qqq/q1(m1)
      common/rrr/ru1(m1)
      common/ddq/dq(m1)
      common/coo/yp1(m1),ym1(m1)
c************************************************************************
      common/rhsc/rhs(m1)
      character*4 ipre2,ipfi
      character*4 ipre
      character*4 dummy, ri3d
      character*60 namfi3
      character*60 namfile
      character*60 filcnw,filcnr,filth,filvm,filpo,filen
      character*60 filet,filer,filez,filed,filev
      common/cordvo/g1m(m1),g1c(m1)
      common/d1/re,pnu
      common/pigr/pi
      common/d2/ntst
      common/dim/n1,n1m
      common/indbo/imv(m1),ipv(m1)
      common/mesh/dx1,dx1q
      common/namefi/filcnw,filcnr,filth,filvm,filpo
      common/nonunif/strr,rext,istr
      common/tstep/dt,beta,ren,gam,rho,sig,thet,zita
      common/parou/devma,tidem,duxc
      common/vperin/vper
c
c   coeff. per tridiagonals in common   
c
      common/intpr/npouth,nprint,npin
c
c   quantities in coetar for hdnl and invtr
c

      common/veltot/vit(1)
c
c     dt variabile
c







