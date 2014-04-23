 161 161 161 3 2. 2. 2.  n1  n2   n3  nsst  alx1  alx2   alx3
 1      1      0   nwrit   nread   iav
 3000. 1.0E-01  .05 1001     20  2    5 Re  vper  tfin  ntst   nprint  npin  npstf
 1001             nstop 
 1    1.5 .2  2.  30.          icfl cflc ,tpin,tprin,tfin
 0   0    .15  1   iles(isc0)  ifilter csma ibox
 0.                 f0
 2                  imic=1  for C.B. spectrum  imic=2 Mansour Wray
   13. 3.    2.     akkpp,qq,sig
80  1  81  6.0  5.0 12.0 7.5    nlr,lipdf,sfuma,qquma,sfoma,qqoma
1

ILES=0  No les
ILES=1  LES: dynamic model
ILES=2  LES: equilibrium model
ILES=3  LES: Smagorinsky model with csma constant

IFILTER=0  Filtro fisico
IFILTER=1  Filtro spettrale

IBOX=0  C(x,y,z,t)
IBOX=1  C(t)

ICFL=0   constant dt
ICFL=1   constant cfl=cflcost

ISTRAT=0 do not evolve density rho
ISTRAT=1 evolve density rho

IGRAD=-1 read background density gradient from file grbar.in
IGRAD=+1 calculate background density gradient brunt-vaisala frequency (bvais) 

ics0 (0 for normal viscosity, 1 dynamic, 2 equil., 3 smagorinsky)
cmsa smagorinsky constant
f0 coriolis parameter (for the earth=2*Omega*sin (latitude)) where Omega=2*pi/day
rho0 vertically averaged density
g    acceleration due to gravity
schm scmidt number (nu/kappa)
bvais brunt-vaisala frequency used when IGRAD=+1 
