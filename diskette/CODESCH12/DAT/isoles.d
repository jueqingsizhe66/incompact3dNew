 65 65 65 3 2. 2. 2.  n1  n2   n3  nsst  alx1  alx2   alx3
 1      0      0  1 nwrit   nread   iav mutim
1.000+06 1.0E-01  .05 501   50  5    5 Re  vper  dt  ntst   nprint  npin  npstf
 501             nstop 
 1    1.5   1.  20.  20.          icfl cflc ,tpin,tprin,tfin
 3   1    .18  0.55  0   iles(isc0)  ifilter csma pratu ibox 
 0.                       f0
 2          imic=1  for C.B. ,  imic=2 Mansour Wray imic=3 Les. Rog.  4 Chasnov
   25. 1.     10.        akkpp,qq,sig
 1 1.0 1.00 1.            istrat,rho0(mks),g(mks),schm
 0 0.000                   igrad,bvais(cph)
0      1.25    960.        iresca,tresca,rlamas
1                          irunpc
1    0                      iturb,idipol
0.5   0.   1.   1.0   1.0     yc1mo,yc2mo,ramo,velmo,vsi
0.    4.                     thet0,omtres
0.    0.                      ar1,ar2
2.25   -1.25   4.  3.0       gell1,gell2,bw1,bw2



ILES=0  No les
ILES=1  LES: dynamic model
ILES=2  LES: equilibrium model
ILES=3  LES: Smagorinsky model with csma constant

IFILTER=0  Filtro fisico
IFILTER=1  Filtro spettrale

IBOX=0  C(x,y,z,t)
IBOX=1  C(t)  for isotropic turbulence

ICFL=0   constant dt
ICFL=1   constant cfl=cflcost


ics0 (0 for normal viscosity, 1 dynamic, 2 equil., 3 smagorinsky)
cmsa smagorinsky constant
f0 coriolis parameter (for the earth=2*Omega*sin (latitude)) where Omega=2*pi/day
