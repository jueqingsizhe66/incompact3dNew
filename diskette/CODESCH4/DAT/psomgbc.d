129 129  3    n1  n2  nsst(=3 for Runge-Kutta 3rd order,=1 A-B)
0     nread  (0 generate ic, otherwise read ic) 
1    nmolp(multiplier for making names of output files)
.05 401  40 4  0  dt  ntst  nfield nhist nscrn(read ic from file nscrn)
1     1.7       .05            icflm,cflma,dtl
20.00    0.2      400.           tprin,tpin,tfin
-40.00  40.00  4.25e-00 0.2  0.7  1 alx1i  alx1f (xmin,xmax) stret xc et istr
0.0   120.0       alx2i  alx2f (ymin,ymax)
0.1e+03  .0  .0     re   beta  ekmn
.3333333 .3333333 .3333333   chal  chbe  chga(Arakawa each=1/3)
-2 1 1 1 inbcvs,vn,vw,ve (vort bc)  -1 no slip only south -1 only west inlet
-2 1 1 1 inbcps,pn,pw,pe (psi)    0 free slip  inbcvs=2 omega=0
-2 1 1 1 inbcss,sn,sw,se (scalar)  1 radiation  inbc*e=-1 for om , psi no slip
0      ib2per
1    1.   npscf   scla 
1.    10.   sch(1)    sch(2)
-5 -20. -20.0 3.5 5.0   1. inm(1dip,0mod,-1alpha,-2gauss) yc1 yc2 ra vel vsi
0.      thet0 (angle for dipole and modon)
1.   1.   ar1   ar2
0    ipert (=1 normal, =2 sin perturb, =3 1+eps*sin(2*pi*om*t))
0  3. 1.0  .0  0.3  -0.0 +0.25  itop,dim1i,dim1f,dim2i,dim2f,hmax,dhmax
1.00   0.2  2.5      u0, blve, cpp (initial velocity, inlet width)
.5  .25   0.5 .5  0.  tau1,ampltp,timepe,dtau,time0
1.    .1           om,eps
1      0.862      1.            nvor akx,alx1d