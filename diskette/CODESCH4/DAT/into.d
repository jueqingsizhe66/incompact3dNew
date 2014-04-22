193 129  3  n1  n2  nsst(=3 for Runge-Kutta 3rd order,=1 A-B)
0  nread  (0 generate ic, otherwise read ic) 
1                         nmolp(multiplier for making names of output files)
.04 3  5 1  0  dt  ntst  nfield nhist nscrn(read ic from file nscrn)
1     1.7    .05            icflm,cflma,dtl
10.00    0.25   80. 0          tprin,tpin,tfin,tpstar
-0.0  20.0  3.00 0.55  0.7 3  alx1i  alx1f (xmin,xmax)  stret  xc et istr
0.0   10.0    alx2i  alx2f (ymin,ymax)
3.0e+03  .0  .0   re   beta  ekmn
.3333333 .3333333 .3333333   chal  chbe  chga(for Arakawa  each=1/3)
2  1  -2 1 inbcvs,vn,vw,ve (vort bc)  -1 no slip only south -1 only west inlet
0  1  -2 1 inbcps,pn,pw,pe (psi)   0 free slip  inbcvs=2 omega=0
0  1  -2 1 inbcss,sn,sw,se (scalar)  1 rad  inbc*e=-1 for om , psi no-slip
0                            ib2per
2    1.                 npscf   scla 
1.   1.         sch(1)    sch(2)
-3 -0. -2.5 1.0 1.0   1.  inmod(1dip,0mod,-1alpha,-2gauss) yc1 yc2 ra vel vsi
38.5     thet0 (angle for dipole and modon)
1.   1.           ar1   ar2
1    ipert (=1 normal, =2 sin perturb, =3 1+eps*sin(2*pi*om*t))
1  9. 1.0  0.  0.3  2. -1.0   itop,dim1i,dim1f,dim2i,dim2f,hmax,dhmax
1.0   1.0  2.5          u0    , blve,cpp
.5    .25   0.5e+09 .5     1.   tau1,ampltp,timepe,dtau,time0
1.    .1           om,eps
1      0.862      1.            nvor akx,alx1d
