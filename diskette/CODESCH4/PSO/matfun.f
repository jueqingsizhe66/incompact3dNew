c***************************************************************
c     
c    THIS SET OF ROUTINES HAVE BEEN COPIED BY THE NUMERICAL RECEIPT
C    BOOK
C**************************************************************
      subroutine spline(x,y,n,yp1,ypn,y2)
      parameter (nmax=500)
      dimension x(nmax),y(nmax),y2(nmax),u(nmax)
      if (yp1.gt..99e30) then
        y2(1)=0.
        u(1)=0.
      else
        y2(1)=-0.5
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))
     *      /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
11    continue
      if (ypn.gt..99e30) then
        qn=0.
        un=0.
      else
        qn=0.5
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue
      return
      end
      subroutine splint(xa,ya,y2a,n,x,y)
      dimension xa(n),ya(n),y2a(n)
      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
      goto 1
      endif
      h=xa(khi)-xa(klo)
      if (h.eq.0.) pause 'bad xa input.'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+
     *      ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
      return
      end

c***************************************************
      function func(flambda)
      common/c1/c,rho,a
      fj1=bessj1(flambda*a)
      fk1=bessk1(rho*a)
      fj2=bessj(2,flambda*a)
      fk2=bessk(2,rho*a)
      func=(fj2/fj1)/flambda+(fk2/fk1)/rho
      return
      end
c***************************************************
      FUNCTION BESSJ(N,X) 
      PARAMETER (IACC=40,BIGNO=1.E10,BIGNI=1.E-10) 
      IF(N.LT.2)PAUSE 'bad argument N in BESSJ' 
      TOX=2./X 
      IF(X.GT.FLOAT(N))THEN 
        BJM=BESSJ0(X) 
        BJ=BESSJ1(X) 
        DO 11 J=1,N-1 
          BJP=J*TOX*BJ-BJM 
          BJM=BJ 
          BJ=BJP 
11      CONTINUE 
        BESSJ=BJ 
      ELSE 
        M=2*((N+INT(SQRT(FLOAT(IACC*N))))/2) 
        BESSJ=0. 
        JSUM=0 
        SUM=0. 
        BJP=0. 
        BJ=1. 
        DO 12 J=M,1,-1 
          BJM=J*TOX*BJ-BJP 
          BJP=BJ 
          BJ=BJM 
          IF(ABS(BJ).GT.BIGNO)THEN 
            BJ=BJ*BIGNI 
            BJP=BJP*BIGNI 
            BESSJ=BESSJ*BIGNI 
            SUM=SUM*BIGNI 
          ENDIF 
          IF(JSUM.NE.0)SUM=SUM+BJ 
          JSUM=1-JSUM 
          IF(J.EQ.N)BESSJ=BJP 
12      CONTINUE 
        SUM=2.*SUM-BJ 
        BESSJ=BESSJ/SUM 
      ENDIF 
      RETURN 
      END 
      FUNCTION BESSK(N,X) 
      IF (N.LT.2) PAUSE 'bad argument N in BESSK' 
      TOX=2.0/X 
      BKM=BESSK0(X) 
      BK=BESSK1(X) 
      DO 11 J=1,N-1 
        BKP=BKM+J*TOX*BK 
        BKM=BK 
        BK=BKP 
11    CONTINUE 
      BESSK=BK 
      RETURN 
      END 
      FUNCTION BESSK0(X) 
      REAL*8 Y,P1,P2,P3,P4,P5,P6,P7, 
     *    Q1,Q2,Q3,Q4,Q5,Q6,Q7 
      DATA P1,P2,P3,P4,P5,P6,P7/-0.57721566D0,0.42278420D0,0.23069756D0, 
     *    0.3488590D-1,0.262698D-2,0.10750D-3,0.74D-5/ 
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7/1.25331414D0,-0.7832358D-1,0.2189568D-1, 
     *    -0.1062446D-1,0.587872D-2,-0.251540D-2,0.53208D-3/ 
      IF (X.LE.2.0) THEN 
        Y=X*X/4.0 
        BESSK0=(-LOG(X/2.0)*BESSI0(X))+(P1+Y*(P2+Y*(P3+ 
     *        Y*(P4+Y*(P5+Y*(P6+Y*P7)))))) 
      ELSE 
        Y=(2.0/X) 
        BESSK0=(EXP(-X)/SQRT(X))*(Q1+Y*(Q2+Y*(Q3+ 
     *        Y*(Q4+Y*(Q5+Y*(Q6+Y*Q7)))))) 
      ENDIF 
      RETURN 
      END 
      FUNCTION BESSK1(X)
      REAL*8 Y,P1,P2,P3,P4,P5,P6,P7,
     *    Q1,Q2,Q3,Q4,Q5,Q6,Q7
      DATA P1,P2,P3,P4,P5,P6,P7/1.0D0,0.15443144D0,-0.67278579D0,
     *    -0.18156897D0,-0.1919402D-1,-0.110404D-2,-0.4686D-4/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7/1.25331414D0,0.23498619D0,-0.3655620D-1,
     *    0.1504268D-1,-0.780353D-2,0.325614D-2,-0.68245D-3/
      IF (X.LE.2.0) THEN
        Y=X*X/4.0
        BESSK1=(LOG(X/2.0)*BESSI1(X))+(1.0/X)*(P1+Y*(P2+
     *      Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
      ELSE
        Y=2.0/X
        BESSK1=(EXP(-X)/SQRT(X))*(Q1+Y*(Q2+Y*(Q3+
     *      Y*(Q4+Y*(Q5+Y*(Q6+Y*Q7))))))
      ENDIF
      RETURN
      END
      FUNCTION BESSI(N,X) 
      PARAMETER(IACC=40,BIGNO=1.0E10,BIGNI=1.0E-10) 
      IF (N.LT.2) PAUSE 'bad argument N in BESSI' 
      TOX=2.0/X 
      BIP=0.0 
      BI=1.0 
      BESSI=0. 
      M=2*((N+INT(SQRT(FLOAT(IACC*N))))) 
      DO 11 J=M,1,-1 
        BIM=BIP+FLOAT(J)*TOX*BI 
        BIP=BI 
        BI=BIM 
        IF (ABS(BI).GT.BIGNO) THEN 
          BESSI=BESSI*BIGNI 
          BI=BI*BIGNI 
          BIP=BIP*BIGNI 
        ENDIF 
        IF (J.EQ.N) BESSI=BIP 
11    CONTINUE 
      BESSI=BESSI*BESSI0(X)/BI 
      RETURN 
      END 
      FUNCTION BESSI0(X) 
      REAL*8 Y,P1,P2,P3,P4,P5,P6,P7, 
     *    Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9 
      DATA P1,P2,P3,P4,P5,P6,P7/1.0D0,3.5156229D0,3.0899424D0,1.2067492D 
     *0, 
     *    0.2659732D0,0.360768D-1,0.45813D-2/ 
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,0.1328592D-1, 
     *    0.225319D-2,-0.157565D-2,0.916281D-2,-0.2057706D-1, 
     *    0.2635537D-1,-0.1647633D-1,0.392377D-2/ 
      IF (ABS(X).LT.3.75) THEN 
        Y=(X/3.75)**2 
        BESSI0=P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))) 
      ELSE 
        AX=ABS(X) 
        Y=3.75/AX 
        BESSI0=(EXP(AX)/SQRT(AX))*(Q1+Y*(Q2+Y*(Q3+Y*(Q4 
     *      +Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9)))))))) 
      ENDIF 
      RETURN 
      END 
      FUNCTION BESSI1(X)
      REAL*8 Y,P1,P2,P3,P4,P5,P6,P7,
     *    Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9
      DATA P1,P2,P3,P4,P5,P6,P7/0.5D0,0.87890594D0,0.51498869D0,
     *    0.15084934D0,0.2658733D-1,0.301532D-2,0.32411D-3/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,-0.3988024D-1,
     *    -0.362018D-2,0.163801D-2,-0.1031555D-1,0.2282967D-1,
     *    -0.2895312D-1,0.1787654D-1,-0.420059D-2/
      IF (ABS(X).LT.3.75) THEN
        Y=(X/3.75)**2
        BESSI1=X*(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
      ELSE
        AX=ABS(X)
        Y=3.75/AX
        BESSI1=(EXP(AX)/SQRT(AX))*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+
     *      Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9))))))))
      ENDIF
      RETURN
      END
      FUNCTION BESSJ0(X)
      REAL*8 Y,P1,P2,P3,P4,P5,Q1,Q2,Q3,Q4,Q5,R1,R2,R3,R4,R5,R6,
     *    S1,S2,S3,S4,S5,S6
      DATA P1,P2,P3,P4,P5/1.D0,-.1098628627D-2,.2734510407D-4,
     *    -.2073370639D-5,.2093887211D-6/, Q1,Q2,Q3,Q4,Q5/-.1562499995D-
     *1,
     *    .1430488765D-3,-.6911147651D-5,.7621095161D-6,-.934945152D-7/
      DATA R1,R2,R3,R4,R5,R6/57568490574.D0,-13362590354.D0,651619640.7D
     *0,
     *    -11214424.18D0,77392.33017D0,-184.9052456D0/,
     *    S1,S2,S3,S4,S5,S6/57568490411.D0,1029532985.D0,
     *    9494680.718D0,59272.64853D0,267.8532712D0,1.D0/
      IF(ABS(X).LT.8.)THEN
        Y=X**2
        BESSJ0=(R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6)))))
     *      /(S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6)))))
      ELSE
        AX=ABS(X)
        Z=8./AX
        Y=Z**2
        XX=AX-.785398164
        BESSJ0=SQRT(.636619772/AX)*(COS(XX)*(P1+Y*(P2+Y*(P3+Y*(P4+Y
     *      *P5))))-Z*SIN(XX)*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)))))
      ENDIF
      RETURN
      END
      FUNCTION BESSJ1(X)
      REAL*8 Y,P1,P2,P3,P4,P5,Q1,Q2,Q3,Q4,Q5,R1,R2,R3,R4,R5,R6,
     *    S1,S2,S3,S4,S5,S6
      DATA R1,R2,R3,R4,R5,R6/72362614232.D0,-7895059235.D0,242396853.1D0
     *,
     *    -2972611.439D0,15704.48260D0,-30.16036606D0/,
     *    S1,S2,S3,S4,S5,S6/144725228442.D0,2300535178.D0,
     *    18583304.74D0,99447.43394D0,376.9991397D0,1.D0/
      DATA P1,P2,P3,P4,P5/1.D0,.183105D-2,-.3516396496D-4,.2457520174D-5
     *,
     *    -.240337019D-6/, Q1,Q2,Q3,Q4,Q5/.04687499995D0,-.2002690873D-3
     *,
     *    .8449199096D-5,-.88228987D-6,.105787412D-6/
      IF(ABS(X).LT.8.)THEN
        Y=X**2
        BESSJ1=X*(R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6)))))
     *      /(S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6)))))
      ELSE
        AX=ABS(X)
        Z=8./AX
        Y=Z**2
        XX=AX-2.356194491
        BESSJ1=SQRT(.636619772/AX)*(COS(XX)*(P1+Y*(P2+Y*(P3+Y*(P4+Y
     *      *P5))))-Z*SIN(XX)*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)))))
     *      *SIGN(1.,X)
      ENDIF
      RETURN
      END
c******************************************************
      FUNCTION ZBRENT(FUNC,X1,X2,TOL) 
      PARAMETER (ITMAX=100,EPS=3.E-8) 
      A=X1 
      B=X2 
      FA=FUNC(A) 
      FB=FUNC(B) 
      IF(FB*FA.GT.0.) PAUSE 'Root must be bracketed for ZBRENT.' 
      FC=FB 
      DO 11 ITER=1,ITMAX 
        IF(FB*FC.GT.0.) THEN 
          C=A 
          FC=FA 
          D=B-A 
          E=D 
        ENDIF 
        IF(ABS(FC).LT.ABS(FB)) THEN 
          A=B 
          B=C 
          C=A 
          FA=FB 
          FB=FC 
          FC=FA 
        ENDIF 
        TOL1=2.*EPS*ABS(B)+0.5*TOL 
        XM=.5*(C-B) 
        IF(ABS(XM).LE.TOL1 .OR. FB.EQ.0.)THEN 
          ZBRENT=B 
          RETURN 
        ENDIF 
        IF(ABS(E).GE.TOL1 .AND. ABS(FA).GT.ABS(FB)) THEN 
          S=FB/FA 
          IF(A.EQ.C) THEN 
            P=2.*XM*S 
            Q=1.-S 
          ELSE 
            Q=FA/FC 
            R=FB/FC 
            P=S*(2.*XM*Q*(Q-R)-(B-A)*(R-1.)) 
            Q=(Q-1.)*(R-1.)*(S-1.) 
          ENDIF 
          IF(P.GT.0.) Q=-Q 
          P=ABS(P) 
          IF(2.*P .LT. MIN(3.*XM*Q-ABS(TOL1*Q),ABS(E*Q))) THEN 
            E=D 
            D=P/Q 
          ELSE 
            D=XM 
            E=D 
          ENDIF 
        ELSE 
          D=XM 
          E=D 
        ENDIF 
        A=B 
        FA=FB 
        IF(ABS(D) .GT. TOL1) THEN 
          B=B+D 
        ELSE 
          B=B+SIGN(TOL1,XM) 
        ENDIF 
        FB=FUNC(B) 
11    CONTINUE 
      PAUSE 'ZBRENT exceeding maximum iterations.' 
      ZBRENT=B 
      RETURN 
      END 
c                                                                       
c    subroutine generating random numbers
c
      subroutine gerand(idum,npa,x)
      dimension r(10000),x(1)
      data iff /0/
      m1=259200
      ia1=7141
      ic1=54773
      rm1=3.8580247e-6
      m2=134456
      ia2=8121
      ic2=28411
      rm2=7.4373773e-6
      m3=243000
      ia3=4561
      ic3=51349
      if (idum.lt.0.or.iff.eq.0) then
        iff=1
        ix1=mod(ic1-idum,m1)
        ix1=mod(ia1*ix1+ic1,m1)
        ix2=mod(ix1,m2)
        ix1=mod(ia1*ix1+ic1,m1)
        ix3=mod(ix1,m3)
        do 11 j=1,npa
          ix1=mod(ia1*ix1+ic1,m1)
          ix2=mod(ia2*ix2+ic2,m2)
          r(j)=(float(ix1)+float(ix2)*rm2)*rm1
11      continue
        idum=1
      endif
      ix1=mod(ia1*ix1+ic1,m1)
      ix2=mod(ia2*ix2+ic2,m2)
      ix3=mod(ia3*ix3+ic3,m3)
      j=1+(npa*ix3)/m3
      if(j.gt.npa) j=npa
      if(j.lt.1) j=1 
c     ran1=r(j)
      ll=0
      do 12 l=j,npa
      ll=ll+1
      x(ll)=r(l)
   12 continue
      do 13 l=1,j-1
      ll=ll+1
      x(ll)=r(l)
   13 continue
      r(j)=(float(ix1)+float(ix2)*rm2)*rm1
      return
      end
