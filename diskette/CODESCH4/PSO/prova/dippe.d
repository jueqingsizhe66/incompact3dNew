96 129    3    n1  n2  nsst(=3 for Runge-Kutta 3rd order,=1 A-B)
0              nread  (0 generate ic, otherwise read ic) 
1          nmolp(multiplier for making names of output files)
.1 401  40 4  0    dt  ntst  nfield nhist nscrn(read ic from file nscrn)
1     1.7       .05            icflm,cflma,dtl
4.00    0.5  40. 0.        tprin,tpin,tfin,tpstar
-5.00  5.00  1.50e-00 0.2  0.7  1 alx1i  alx1f (xmin,xmax)  stretching  xcra  etra   istr
-5.   5.0                       alx2i  alx2f (ymin,ymax)
8.055e+04  .0   0.0e-03                 re   beta ekmn
.3333333 .3333333 .3333333      chal  chbe  chga(for arakawa scheme each=1/3)
1    1      1    1  inbcvs,vn,vw,ve (vort bc)  -1 no slip only south -1 only west inlet
1    1      1    1  inbcps,pn,pw,pe (psi)       0 free slip  inbcvs=2 omega=0
1    1      1    1  inbcss,sn,sw,se (scalar)  1 radiation  inbc*e=-1 for om and psi no slip
1                    ib2per
1    1.       npscf   scla 
1.    10.               sch(1)    sch(2)
1 0. 0. 1.0 1.0  -1.0   inmod yc1mo1 yc2mo ramo velmo vsi(1dip,0mod,-1alpha,-2gauss) 
90.      thet0 (angle for dipole and modon)
1.   1.      ar1   ar2
0  3. 1.0 .0  0.3  -0.0 +0.25  itop,dim1i,dim1f,dim2i,dim2f,hmax,dhmax
1.00   1.0  2.5      u0, blve, cpp (initial velocity, inlet width)
.5    .25   0.5e+00 .5     0.   tau1,ampltp,timepe,dtau,time0
1.    .1              om,eps
1      0.862  1   1.            nvor akx,alx1d
