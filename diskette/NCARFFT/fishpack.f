
      SUBROUTINE BLKTRI (IFLG,NP,N,AN,BN,CN,MP,M,AM,BM,CM,IDIMY,Y,
     1                   IERROR,W)
C
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C     *                                                               *
C     *                        F I S H P A K                          *
C     *                                                               *
C     *                                                               *
C     *     A PACKAGE OF FORTRAN SUBPROGRAMS FOR THE SOLUTION OF      *
C     *                                                               *
C     *      SEPARABLE ELLIPTIC PARTIAL DIFFERENTIAL EQUATIONS        *
C     *                                                               *
C     *                  (VERSION 3.1 , OCTOBER 1980)                  *
C     *                                                               *
C     *                             BY                                *
C     *                                                               *
C     *        JOHN ADAMS, PAUL SWARZTRAUBER AND ROLAND SWEET         *
C     *                                                               *
C     *                             OF                                *
C     *                                                               *
C     *         THE NATIONAL CENTER FOR ATMOSPHERIC RESEARCH          *
C     *                                                               *
C     *                BOULDER, COLORADO  (80307)  U.S.A.             *
C     *                                                               *
C     *                   WHICH IS SPONSORED BY                       *
C     *                                                               *
C     *              THE NATIONAL SCIENCE FOUNDATION                  *
C     *                                                               *
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C     SUBROUTINE BLKTRI SOLVES A SYSTEM OF LINEAR EQUATIONS OF THE FORM
C
C          AN(J)*X(I,J-1) + AM(I)*X(I-1,J) + (BN(J)+BM(I))*X(I,J)
C
C          + CN(J)*X(I,J+1) + CM(I)*X(I+1,J) = Y(I,J)
C
C               FOR I = 1,2,...,M  AND  J = 1,2,...,N.
C
C     I+1 AND I-1 ARE EVALUATED MODULO M AND J+1 AND J-1 MODULO N, I.E.,
C
C          X(I,0) = X(I,N),  X(I,N+1) = X(I,1),
C          X(0,J) = X(M,J),  X(M+1,J) = X(1,J).
C
C     THESE EQUATIONS USUALLY RESULT FROM THE DISCRETIZATION OF
C     SEPARABLE ELLIPTIC EQUATIONS.  BOUNDARY CONDITIONS MAY BE
C     DIRICHLET, NEUMANN, OR PERIODIC.
C
C
C     * * * * * * * * * *     ON INPUT     * * * * * * * * * *
C
C     IFLG
C       = 0  INITIALIZATION ONLY.  CERTAIN QUANTITIES THAT DEPEND ON NP,
C            N, AN, BN, AND CN ARE COMPUTED AND STORED IN THE WORK
C            ARRAY  W.
C       = 1  THE QUANTITIES THAT WERE COMPUTED IN THE INITIALIZATION ARE
C            USED TO OBTAIN THE SOLUTION X(I,J).
C
C       NOTE   A CALL WITH IFLG=0 TAKES APPROXIMATELY ONE HALF THE TIME
C              TIME AS A CALL WITH IFLG = 1  .  HOWEVER, THE
C              INITIALIZATION DOES NOT HAVE TO BE REPEATED UNLESS NP, N,
C              AN, BN, OR CN CHANGE.
C
C     NP
C       = 0  IF AN(1) AND CN(N) ARE NOT ZERO, WHICH CORRESPONDS TO
C            PERIODIC BOUNARY CONDITIONS.
C       = 1  IF AN(1) AND CN(N) ARE ZERO.
C
C     N
C       THE NUMBER OF UNKNOWNS IN THE J-DIRECTION. N MUST BE GREATER
C       THAN 4. THE OPERATION COUNT IS PROPORTIONAL TO MNLOG2(N), HENCE
C       N SHOULD BE SELECTED LESS THAN OR EQUAL TO M.
C
C     AN,BN,CN
C       ONE-DIMENSIONAL ARRAYS OF LENGTH N THAT SPECIFY THE COEFFICIENTS
C       IN THE LINEAR EQUATIONS GIVEN ABOVE.
C
C     MP
C       = 0  IF AM(1) AND CM(M) ARE NOT ZERO, WHICH CORRESPONDS TO
C            PERIODIC BOUNDARY CONDITIONS.
C       = 1  IF AM(1) = CM(M) = 0  .
C
C     M
C       THE NUMBER OF UNKNOWNS IN THE I-DIRECTION. M MUST BE GREATER
C       THAN 4.
C
C     AM,BM,CM
C       ONE-DIMENSIONAL ARRAYS OF LENGTH M THAT SPECIFY THE COEFFICIENTS
C       IN THE LINEAR EQUATIONS GIVEN ABOVE.
C
C     IDIMY
C       THE ROW (OR FIRST) DIMENSION OF THE TWO-DIMENSIONAL ARRAY Y AS
C       IT APPEARS IN THE PROGRAM CALLING BLKTRI.  THIS PARAMETER IS
C       USED TO SPECIFY THE VARIABLE DIMENSION OF Y.  IDIMY MUST BE AT
C       LEAST M.
C
C     Y
C       A TWO-DIMENSIONAL ARRAY THAT SPECIFIES THE VALUES OF THE RIGHT
C       SIDE OF THE LINEAR SYSTEM OF EQUATIONS GIVEN ABOVE.  Y MUST BE
C       DIMENSIONED AT LEAST M*N.
C
C     W
C       A ONE-DIMENSIONAL ARRAY THAT MUST BE PROVIDED BY THE USER FOR
C       WORK SPACE.
C             IF NP=1 DEFINE K=INT(LOG2(N))+1 AND SET L=2**(K+1) THEN
C                     W MUST HAVE DIMENSION (K-2)*L+K+5+MAX(2N,6M)
C
C             IF NP=0 DEFINE K=INT(LOG2(N-1))+1 AND SET L=2**(K+1) THEN
C                     W MUST HAVE DIMENSION (K-2)*L+K+5+2N+MAX(2N,6M)
C
C       **IMPORTANT** FOR PURPOSES OF CHECKING, THE REQUIRED DIMENSION
C                     OF W IS COMPUTED BY BLKTRI AND STORED IN W(1)
C                     IN FLOATING POINT FORMAT.
C
C     * * * * * * * * * *     ON OUTPUT     * * * * * * * * * *
C
C     Y
C       CONTAINS THE SOLUTION X.
C
C     IERROR
C       AN ERROR FLAG THAT INDICATES INVALID INPUT PARAMETERS.  EXCEPT
C       FOR NUMBER ZERO, A SOLUTION IS NOT ATTEMPTED.
C
C       = 0  NO ERROR.
C       = 1  M IS LESS THAN 5
C       = 2  N IS LESS THAN 5
C       = 3  IDIMY IS LESS THAN M.
C       = 4  BLKTRI FAILED WHILE COMPUTING RESULTS THAT DEPEND ON THE
C            COEFFICIENT ARRAYS AN, BN, CN.  CHECK THESE ARRAYS.
C       = 5  AN(J)*CN(J-1) IS LESS THAN 0 FOR SOME J. POSSIBLE REASONS
C            FOR THIS CONDITION ARE
C            1. THE ARRAYS AN AND CN ARE NOT CORRECT
C            2. TOO LARGE A GRID SPACING WAS USED IN THE DISCRETIZATION
C               OF THE ELLIPTIC EQUATION
C            3. THE LINEAR EQUATIONS RESULTED FROM A PARTIAL
C               DIFFERENTIAL EQUATION WHICH WAS NOT ELLIPTIC
C
C     W
C       CONTAINS INTERMEDIATE VALUES THAT MUST NOT BE DESTROYED IF
C       BLKTRI WILL BE CALLED AGAIN WITH IFLG=1. W(1) CONTAINS THE
C       NUMBER OF LOCATIONS REQUIRED BY W IN FLOATING POINT FORMAT.
C
C     * * * * * * *   PROGRAM SPECIFICATIONS    * * * * * * * * * * * *
C
C     DIMENSION OF   AN(N),BN(N),CN(N),AM(M),BM(M),CM(M),Y(IDIMY,N)
C     ARGUMENTS      W(SEE ARGUMENT LIST)
C
C     LATEST         JUNE 1979
C     REVISION
C
C     REQUIRED       BLKTRI,BLKTRI,PROD,PRODP,CPROD,CPRODP,COMPB,INDXA,
C     SUBPROGRAMS    INDXB,INDXC,PPADD,PSGF,PPSGF,PPSPF,BSRH,TEVLS,
C                    EPMACH,STORE
C
C     SPECIAL        THE ALGORITHM MAY FAIL IF ABS(BM(I)+BN(J)) IS LESS
C     CONDITIONS     THAN ABS(AM(I))+ABS(AN(J))+ABS(CM(I))+ABS(CN(J))
C                    FOR SOME I AND J. THE ALGORITHM WILL ALSO FAIL IF
C                    AN(J)*CN(J-1) IS LESS THAN ZERO FOR SOME J
C                    SEE THE DISCRIPTION OF THE OUTPUT PARAMETER IERROR.
C
C     COMMON         CBLKT,VALUE
C     BLOCKS
C
C     I/O            NONE
C
C     PRECISION      SINGLE
C
C     SPECIALIST     PAUL SWARZTRAUBER
C
C     LANGUAGE       FORTRAN
C
C     HISTORY        VERSION 1 SEPTEMBER 1973
C                    VERSION 2 APRIL     1976
C                    VERSION 3 JUNE      1979
C
C     ALGORITHM      GENERALIZED CYCLIC REDUCTION (SEE REFERENCE BELOW)
C
C     SPACE
C     REQUIRED       CONTROL DATA 7600
C
C     PORTABILITY    AMERICAN NATIONAL STANDARDS INSTITUTE FORTRAN.
C                    THE APPROXIMATE MACHINE ACCURACY IS COMPUTED IN
C                    FUNCTION EPMACH
C
C     REQUIRED       NONE
C     RESIDENT
C     ROUTINES
C
C     REFERENCES     SWARZTRAUBER,P. AND R. SWEET, 'EFFICIENT FORTRAN
C                    SUBPROGRAMS FOR THE SOLUTION OF ELLIPTIC EQUATIONS'
C                    NCAR TN/IA-109, JULY, 1975, 138 PP.
C
C                    SWARZTRAUBER P. N.,A DIRECT METHOD FOR THE DISCRETE
C                    SOLUTION OF SEPARABLE ELLIPTIC EQUATIONS, S.I.A.M.
C                    J. NUMER. ANAL.,11(1974) PP. 1136-1150.
C
C     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C
C
      DIMENSION       AN(1)      ,BN(1)      ,CN(1)      ,AM(1)      ,
     1                BM(1)      ,CM(1)      ,Y(IDIMY,1) ,W(1)
      EXTERNAL        PROD       ,PRODP      ,CPROD      ,CPRODP
      COMMON /CBLKT/  NPP        ,K          ,EPS        ,CNV        ,
     1                NM         ,NCMPLX     ,IK
C
C TEST M AND N FOR THE PROPER FORM
C
      NM = N
      IERROR = 0
      IF (M-5) 101,102,102
  101 IERROR = 1
      GO TO 119
  102 IF (NM-3) 103,104,104
  103 IERROR = 2
      GO TO 119
  104 IF (IDIMY-M) 105,106,106
  105 IERROR = 3
      GO TO 119
  106 NH = N
      NPP = NP
      IF (NPP) 107,108,107
  107 NH = NH+1
  108 IK = 2
      K = 1
  109 IK = IK+IK
      K = K+1
      IF (NH-IK) 110,110,109
  110 NL = IK
      IK = IK+IK
      NL = NL-1
      IWAH = (K-2)*IK+K+6
      IF (NPP) 111,112,111
C
C     DIVIDE W INTO WORKING SUB ARRAYS
C
  111 IW1 = IWAH
      IWBH = IW1+NM
      W(1) = FLOAT(IW1-1+MAX0(2*NM,6*M))
      GO TO 113
  112 IWBH = IWAH+NM+NM
      IW1 = IWBH
      W(1) = FLOAT(IW1-1+MAX0(2*NM,6*M))
      NM = NM-1
C
C SUBROUTINE COMP B COMPUTES THE ROOTS OF THE B POLYNOMIALS
C
  113 IF (IERROR) 119,114,119
  114 IW2 = IW1+M
      IW3 = IW2+M
      IWD = IW3+M
      IWW = IWD+M
      IWU = IWW+M
      IF (IFLG) 116,115,116
  115 CALL COMPB (NL,IERROR,AN,BN,CN,W(2),W(IWAH),W(IWBH))
      GO TO 119
  116 IF (MP) 117,118,117
C
C SUBROUTINE BLKTR1 SOLVES THE LINEAR SYSTEM
C
  117 CALL BLKTR1 (NL,AN,BN,CN,M,AM,BM,CM,IDIMY,Y,W(2),W(IW1),W(IW2),
     1             W(IW3),W(IWD),W(IWW),W(IWU),PROD,CPROD)
      GO TO 119
  118 CALL BLKTR1 (NL,AN,BN,CN,M,AM,BM,CM,IDIMY,Y,W(2),W(IW1),W(IW2),
     1             W(IW3),W(IWD),W(IWW),W(IWU),PRODP,CPRODP)
  119 CONTINUE
      RETURN
      END
      subroutine cprodp (nd,bd,nm1,bm1,nm2,bm2,na,aa,x,yy,m,a,b,c,d,u,y)
c
c prodp applies a sequence of matrix operations to the vector x and
c stores the result in yy       periodic boundary conditions
c and  complex  case
c
c bd,bm1,bm2 are arrays containing roots of certian b polynomials
c nd,nm1,nm2 are the lengths of the arrays bd,bm1,bm2 respectively
c aa   array containing scalar multipliers of the vector x
c na is the length of the array aa
c x,yy the matrix operations are applied to x and the result is yy
c a,b,c  are arrays which contain the tridiagonal matrix
c m  is the order of the matrix
c d,u,y are working arrays
c isgn  determines whether or not a change in sign is made
c
      complex         y          ,d          ,u          ,v          ,
     1                den        ,bh         ,ym         ,am         ,
     2                y1         ,y2         ,yh         ,bd         ,
     3                crt
      dimension       a(1)       ,b(1)       ,c(1)       ,x(1)       ,
     1                y(1)       ,d(1)       ,u(1)       ,bd(1)      ,
     2                bm1(1)     ,bm2(1)     ,aa(1)      ,yy(1)
      do 101 j=1,m
         y(j) = cmplx(x(j),0.)
  101 continue
      mm = m-1
      mm2 = m-2
      id = nd
      m1 = nm1
      m2 = nm2
      ia = na
  102 iflg = 0
      if (id) 111,111,103
  103 crt = bd(id)
      id = id-1
      iflg = 1
c
c begin solution to system
c
      bh = b(m)-crt
      ym = y(m)
      den = b(1)-crt
      d(1) = c(1)/den
      u(1) = a(1)/den
      y(1) = y(1)/den
      v = cmplx(c(m),0.)
      if (mm2-2) 106,104,104
  104 do 105 j=2,mm2
         den = b(j)-crt-a(j)*d(j-1)
         d(j) = c(j)/den
         u(j) = -a(j)*u(j-1)/den
         y(j) = (y(j)-a(j)*y(j-1))/den
         bh = bh-v*u(j-1)
         ym = ym-v*y(j-1)
         v = -v*d(j-1)
  105 continue
  106 den = b(m-1)-crt-a(m-1)*d(m-2)
      d(m-1) = (c(m-1)-a(m-1)*u(m-2))/den
      y(m-1) = (y(m-1)-a(m-1)*y(m-2))/den
      am = a(m)-v*d(m-2)
      bh = bh-v*u(m-2)
      ym = ym-v*y(m-2)
      den = bh-am*d(m-1)
      if (cabs(den)) 107,108,107
  107 y(m) = (ym-am*y(m-1))/den
      go to 109
  108 y(m) = (1.,0.)
  109 y(m-1) = y(m-1)-d(m-1)*y(m)
      do 110 j=2,mm
         k = m-j
         y(k) = y(k)-d(k)*y(k+1)-u(k)*y(m)
  110 continue
  111 if (m1) 112,112,114
  112 if (m2) 123,123,113
  113 rt = bm2(m2)
      m2 = m2-1
      go to 119
  114 if (m2) 115,115,116
  115 rt = bm1(m1)
      m1 = m1-1
      go to 119
  116 if (abs(bm1(m1))-abs(bm2(m2))) 118,118,117
  117 rt = bm1(m1)
      m1 = m1-1
      go to 119
  118 rt = bm2(m2)
      m2 = m2-1
c
c matrix multiplication
c
  119 yh = y(1)
      y1 = (b(1)-rt)*y(1)+c(1)*y(2)+a(1)*y(m)
      if (mm-2) 122,120,120
  120 do 121 j=2,mm
         y2 = a(j)*y(j-1)+(b(j)-rt)*y(j)+c(j)*y(j+1)
         y(j-1) = y1
         y1 = y2
  121 continue
  122 y(m) = a(m)*y(m-1)+(b(m)-rt)*y(m)+c(m)*yh
      y(m-1) = y1
      iflg = 1
      go to 102
  123 if (ia) 126,126,124
  124 rt = aa(ia)
      ia = ia-1
      iflg = 1
c
c scalar multiplication
c
      do 125 j=1,m
         y(j) = rt*y(j)
  125 continue
  126 if (iflg) 127,127,102
  127 do 128 j=1,m
         yy(j) = real(y(j))
  128 continue
      return
      end
      subroutine prodp (nd,bd,nm1,bm1,nm2,bm2,na,aa,x,y,m,a,b,c,d,u,w)
c
c prodp applies a sequence of matrix operations to the vector x and
c stores the result in y        periodic boundary conditions
c
c bd,bm1,bm2 are arrays containing roots of certian b polynomials
c nd,nm1,nm2 are the lengths of the arrays bd,bm1,bm2 respectively
c aa   array containing scalar multipliers of the vector x
c na is the length of the array aa
c x,y  the matrix operations are applied to x and the result is y
c a,b,c  are arrays which contain the tridiagonal matrix
c m  is the order of the matrix
c d,u,w are working arrays
c is  determines whether or not a change in sign is made
c
      dimension       a(1)       ,b(1)       ,c(1)       ,x(1)       ,
     1                y(1)       ,d(1)       ,u(1)       ,bd(1)      ,
     2                bm1(1)     ,bm2(1)     ,aa(1)      ,w(1)
      do 101 j=1,m
         y(j) = x(j)
         w(j) = y(j)
  101 continue
      mm = m-1
      mm2 = m-2
      id = nd
      ibr = 0
      m1 = nm1
      m2 = nm2
      ia = na
  102 if (ia) 105,105,103
  103 rt = aa(ia)
      if (nd .eq. 0) rt = -rt
      ia = ia-1
      do 104 j=1,m
         y(j) = rt*w(j)
  104 continue
  105 if (id) 128,128,106
  106 rt = bd(id)
      id = id-1
      if (id .eq. 0) ibr = 1
c
c begin solution to system
c
      bh = b(m)-rt
      ym = y(m)
      den = b(1)-rt
      d(1) = c(1)/den
      u(1) = a(1)/den
      w(1) = y(1)/den
      v = c(m)
      if (mm2-2) 109,107,107
  107 do 108 j=2,mm2
         den = b(j)-rt-a(j)*d(j-1)
         d(j) = c(j)/den
         u(j) = -a(j)*u(j-1)/den
         w(j) = (y(j)-a(j)*w(j-1))/den
         bh = bh-v*u(j-1)
         ym = ym-v*w(j-1)
         v = -v*d(j-1)
  108 continue
  109 den = b(m-1)-rt-a(m-1)*d(m-2)
      d(m-1) = (c(m-1)-a(m-1)*u(m-2))/den
      w(m-1) = (y(m-1)-a(m-1)*w(m-2))/den
      am = a(m)-v*d(m-2)
      bh = bh-v*u(m-2)
      ym = ym-v*w(m-2)
      den = bh-am*d(m-1)
      if (den) 110,111,110
  110 w(m) = (ym-am*w(m-1))/den
      go to 112
  111 w(m) = 1.
  112 w(m-1) = w(m-1)-d(m-1)*w(m)
      do 113 j=2,mm
         k = m-j
         w(k) = w(k)-d(k)*w(k+1)-u(k)*w(m)
  113 continue
      if (na) 116,116,102
  114 do 115 j=1,m
         y(j) = w(j)
  115 continue
      ibr = 1
      go to 102
  116 if (m1) 117,117,118
  117 if (m2) 114,114,123
  118 if (m2) 120,120,119
  119 if (abs(bm1(m1))-abs(bm2(m2))) 123,123,120
  120 if (ibr) 121,121,122
  121 if (abs(bm1(m1)-bd(id))-abs(bm1(m1)-rt)) 114,122,122
  122 rt = rt-bm1(m1)
      m1 = m1-1
      go to 126
  123 if (ibr) 124,124,125
  124 if (abs(bm2(m2)-bd(id))-abs(bm2(m2)-rt)) 114,125,125
  125 rt = rt-bm2(m2)
      m2 = m2-1
  126 do 127 j=1,m
         y(j) = y(j)+rt*w(j)
  127 continue
      go to 102
  128 return
      end
      subroutine cprod (nd,bd,nm1,bm1,nm2,bm2,na,aa,x,yy,m,a,b,c,d,w,y)
c
c prod applies a sequence of matrix operations to the vector x and
c stores the result in yy           (complex case)
c aa   array containing scalar multipliers of the vector x
c nd,nm1,nm2 are the lengths of the arrays bd,bm1,bm2 respectively
c bd,bm1,bm2 are arrays containing roots of certian b polynomials
c na is the length of the array aa
c x,yy the matrix operations are applied to x and the result is yy
c a,b,c  are arrays which contain the tridiagonal matrix
c m  is the order of the matrix
c d,w,y are working arrays
c isgn  determines whether or not a change in sign is made
c
      complex         y          ,d          ,w          ,bd         ,
     1                crt        ,den        ,y1         ,y2
      dimension       a(1)       ,b(1)       ,c(1)       ,x(1)       ,
     1                y(1)       ,d(1)       ,w(1)       ,bd(1)      ,
     2                bm1(1)     ,bm2(1)     ,aa(1)      ,yy(1)
      do 101 j=1,m
         y(j) = cmplx(x(j),0.)
  101 continue
      mm = m-1
      id = nd
      m1 = nm1
      m2 = nm2
      ia = na
  102 iflg = 0
      if (id) 109,109,103
  103 crt = bd(id)
      id = id-1
c
c begin solution to system
c
      d(m) = a(m)/(b(m)-crt)
      w(m) = y(m)/(b(m)-crt)
      do 104 j=2,mm
         k = m-j
         den = b(k+1)-crt-c(k+1)*d(k+2)
         d(k+1) = a(k+1)/den
         w(k+1) = (y(k+1)-c(k+1)*w(k+2))/den
  104 continue
      den = b(1)-crt-c(1)*d(2)
      if (cabs(den)) 105,106,105
  105 y(1) = (y(1)-c(1)*w(2))/den
      go to 107
  106 y(1) = (1.,0.)
  107 do 108 j=2,m
         y(j) = w(j)-d(j)*y(j-1)
  108 continue
  109 if (m1) 110,110,112
  110 if (m2) 121,121,111
  111 rt = bm2(m2)
      m2 = m2-1
      go to 117
  112 if (m2) 113,113,114
  113 rt = bm1(m1)
      m1 = m1-1
      go to 117
  114 if (abs(bm1(m1))-abs(bm2(m2))) 116,116,115
  115 rt = bm1(m1)
      m1 = m1-1
      go to 117
  116 rt = bm2(m2)
      m2 = m2-1
  117 y1 = (b(1)-rt)*y(1)+c(1)*y(2)
      if (mm-2) 120,118,118
c
c matrix multiplication
c
  118 do 119 j=2,mm
         y2 = a(j)*y(j-1)+(b(j)-rt)*y(j)+c(j)*y(j+1)
         y(j-1) = y1
         y1 = y2
  119 continue
  120 y(m) = a(m)*y(m-1)+(b(m)-rt)*y(m)
      y(m-1) = y1
      iflg = 1
      go to 102
  121 if (ia) 124,124,122
  122 rt = aa(ia)
      ia = ia-1
      iflg = 1
c
c scalar multiplication
c
      do 123 j=1,m
         y(j) = rt*y(j)
  123 continue
  124 if (iflg) 125,125,102
  125 do 126 j=1,m
         yy(j) = real(y(j))
  126 continue
      return
      end
      subroutine prod (nd,bd,nm1,bm1,nm2,bm2,na,aa,x,y,m,a,b,c,d,w,u)
c
c prod applies a sequence of matrix operations to the vector x and
c stores the result in y
c bd,bm1,bm2 are arrays containing roots of certian b polynomials
c nd,nm1,nm2 are the lengths of the arrays bd,bm1,bm2 respectively
c aa   array containing scalar multipliers of the vector x
c na is the length of the array aa
c x,y  the matrix operations are applied to x and the result is y
c a,b,c  are arrays which contain the tridiagonal matrix
c m  is the order of the matrix
c d,w,u are working arrays
c is  determines whether or not a change in sign is made
c
      dimension       a(1)       ,b(1)       ,c(1)       ,x(1)       ,
     1                y(1)       ,d(1)       ,w(1)       ,bd(1)      ,
     2                bm1(1)     ,bm2(1)     ,aa(1)      ,u(1)
      do 101 j=1,m
         w(j) = x(j)
         y(j) = w(j)
  101 continue
      mm = m-1
      id = nd
      ibr = 0
      m1 = nm1
      m2 = nm2
      ia = na
  102 if (ia) 105,105,103
  103 rt = aa(ia)
      if (nd .eq. 0) rt = -rt
      ia = ia-1
c
c scalar multiplication
c
      do 104 j=1,m
         y(j) = rt*w(j)
  104 continue
  105 if (id) 125,125,106
  106 rt = bd(id)
      id = id-1
      if (id .eq. 0) ibr = 1
c
c begin solution to system
c
      d(m) = a(m)/(b(m)-rt)
      w(m) = y(m)/(b(m)-rt)
      do 107 j=2,mm
         k = m-j
         den = b(k+1)-rt-c(k+1)*d(k+2)
         d(k+1) = a(k+1)/den
         w(k+1) = (y(k+1)-c(k+1)*w(k+2))/den
  107 continue
      den = b(1)-rt-c(1)*d(2)
      w(1) = 1.
      if (den) 108,109,108
  108 w(1) = (y(1)-c(1)*w(2))/den
  109 do 110 j=2,m
         w(j) = w(j)-d(j)*w(j-1)
  110 continue
      if (na) 113,113,102
  111 do 112 j=1,m
         y(j) = w(j)
  112 continue
      ibr = 1
      go to 102
  113 if (m1) 114,114,115
  114 if (m2) 111,111,120
  115 if (m2) 117,117,116
  116 if (abs(bm1(m1))-abs(bm2(m2))) 120,120,117
  117 if (ibr) 118,118,119
  118 if (abs(bm1(m1)-bd(id))-abs(bm1(m1)-rt)) 111,119,119
  119 rt = rt-bm1(m1)
      m1 = m1-1
      go to 123
  120 if (ibr) 121,121,122
  121 if (abs(bm2(m2)-bd(id))-abs(bm2(m2)-rt)) 111,122,122
  122 rt = rt-bm2(m2)
      m2 = m2-1
  123 do 124 j=1,m
         y(j) = y(j)+rt*w(j)
  124 continue
      go to 102
  125 return
      end
      subroutine blktr1 (n,an,bn,cn,m,am,bm,cm,idimy,y,b,w1,w2,w3,wd,
     1                   ww,wu,prdct,cprdct)
c
c blktr1 solves the linear system
c
c b  contains the roots of all the b polynomials
c w1,w2,w3,wd,ww,wu  are all working arrays
c prdct  is either prodp or prod depending on whether the boundary
c conditions in the m direction are periodic or not
c cprdct is either cprodp or cprod which are the complex versions
c of prodp and prod. these are called in the event that some
c of the roots of the b sub p polynomial are complex
c
c
      dimension       an(1)      ,bn(1)      ,cn(1)      ,am(1)      ,
     1                bm(1)      ,cm(1)      ,b(1)       ,w1(1)      ,
     2                w2(1)      ,w3(1)      ,wd(1)      ,ww(1)      ,
     3                wu(1)      ,y(idimy,1)
      common /cblkt/  npp        ,k          ,eps        ,cnv        ,
     1                nm         ,ncmplx     ,ik
c
c begin reduction phase
c
      kdo = k-1
      do 109 l=1,kdo
         ir = l-1
         i2 = 2**ir
         i1 = i2/2
         i3 = i2+i1
         i4 = i2+i2
         irm1 = ir-1
         call indxb (i2,ir,im2,nm2)
         call indxb (i1,irm1,im3,nm3)
         call indxb (i3,irm1,im1,nm1)
         call prdct (nm2,b(im2),nm3,b(im3),nm1,b(im1),0,dum,y(1,i2),w3,
     1               m,am,bm,cm,wd,ww,wu)
         if = 2**k
         do 108 i=i4,if,i4
            if (i-nm) 101,101,108
  101       ipi1 = i+i1
            ipi2 = i+i2
            ipi3 = i+i3
            call indxc (i,ir,idxc,nc)
            if (i-if) 102,108,108
  102       call indxa (i,ir,idxa,na)
            call indxb (i-i1,irm1,im1,nm1)
            call indxb (ipi2,ir,ip2,np2)
            call indxb (ipi1,irm1,ip1,np1)
            call indxb (ipi3,irm1,ip3,np3)
            call prdct (nm1,b(im1),0,dum,0,dum,na,an(idxa),w3,w1,m,am,
     1                  bm,cm,wd,ww,wu)
            if (ipi2-nm) 105,105,103
  103       do 104 j=1,m
               w3(j) = 0.
               w2(j) = 0.
  104       continue
            go to 106
  105       call prdct (np2,b(ip2),np1,b(ip1),np3,b(ip3),0,dum,
     1                  y(1,ipi2),w3,m,am,bm,cm,wd,ww,wu)
            call prdct (np1,b(ip1),0,dum,0,dum,nc,cn(idxc),w3,w2,m,am,
     1                  bm,cm,wd,ww,wu)
  106       do 107 j=1,m
               y(j,i) = w1(j)+w2(j)+y(j,i)
  107       continue
  108    continue
  109 continue
      if (npp) 132,110,132
c
c     the periodic case is treated using the capacitance matrix method
c
  110 if = 2**k
      i = if/2
      i1 = i/2
      call indxb (i-i1,k-2,im1,nm1)
      call indxb (i+i1,k-2,ip1,np1)
      call indxb (i,k-1,iz,nz)
      call prdct (nz,b(iz),nm1,b(im1),np1,b(ip1),0,dum,y(1,i),w1,m,am,
     1            bm,cm,wd,ww,wu)
      izr = i
      do 111 j=1,m
         w2(j) = w1(j)
  111 continue
      do 113 ll=2,k
         l = k-ll+1
         ir = l-1
         i2 = 2**ir
         i1 = i2/2
         i = i2
         call indxc (i,ir,idxc,nc)
         call indxb (i,ir,iz,nz)
         call indxb (i-i1,ir-1,im1,nm1)
         call indxb (i+i1,ir-1,ip1,np1)
         call prdct (np1,b(ip1),0,dum,0,dum,nc,cn(idxc),w1,w1,m,am,bm,
     1               cm,wd,ww,wu)
         do 112 j=1,m
            w1(j) = y(j,i)+w1(j)
  112    continue
         call prdct (nz,b(iz),nm1,b(im1),np1,b(ip1),0,dum,w1,w1,m,am,
     1               bm,cm,wd,ww,wu)
  113 continue
      do 118 ll=2,k
         l = k-ll+1
         ir = l-1
         i2 = 2**ir
         i1 = i2/2
         i4 = i2+i2
         ifd = if-i2
         do 117 i=i2,ifd,i4
            if (i-i2-izr) 117,114,117
  114       if (i-nm) 115,115,118
  115       call indxa (i,ir,idxa,na)
            call indxb (i,ir,iz,nz)
            call indxb (i-i1,ir-1,im1,nm1)
            call indxb (i+i1,ir-1,ip1,np1)
            call prdct (nm1,b(im1),0,dum,0,dum,na,an(idxa),w2,w2,m,am,
     1                  bm,cm,wd,ww,wu)
            do 116 j=1,m
               w2(j) = y(j,i)+w2(j)
  116       continue
            call prdct (nz,b(iz),nm1,b(im1),np1,b(ip1),0,dum,w2,w2,m,
     1                  am,bm,cm,wd,ww,wu)
            izr = i
            if (i-nm) 117,119,117
  117    continue
  118 continue
  119 do 120 j=1,m
         y(j,nm+1) = y(j,nm+1)-cn(nm+1)*w1(j)-an(nm+1)*w2(j)
  120 continue
      call indxb (if/2,k-1,im1,nm1)
      call indxb (if,k-1,ip,np)
      if (ncmplx) 121,122,121
  121 call cprdct (nm+1,b(ip),nm1,b(im1),0,dum,0,dum,y(1,nm+1),
     1             y(1,nm+1),m,am,bm,cm,w1,w3,ww)
      go to 123
  122 call prdct (nm+1,b(ip),nm1,b(im1),0,dum,0,dum,y(1,nm+1),
     1            y(1,nm+1),m,am,bm,cm,wd,ww,wu)
  123 do 124 j=1,m
         w1(j) = an(1)*y(j,nm+1)
         w2(j) = cn(nm)*y(j,nm+1)
         y(j,1) = y(j,1)-w1(j)
         y(j,nm) = y(j,nm)-w2(j)
  124 continue
      do 126 l=1,kdo
         ir = l-1
         i2 = 2**ir
         i4 = i2+i2
         i1 = i2/2
         i = i4
         call indxa (i,ir,idxa,na)
         call indxb (i-i2,ir,im2,nm2)
         call indxb (i-i2-i1,ir-1,im3,nm3)
         call indxb (i-i1,ir-1,im1,nm1)
         call prdct (nm2,b(im2),nm3,b(im3),nm1,b(im1),0,dum,w1,w1,m,am,
     1               bm,cm,wd,ww,wu)
         call prdct (nm1,b(im1),0,dum,0,dum,na,an(idxa),w1,w1,m,am,bm,
     1               cm,wd,ww,wu)
         do 125 j=1,m
            y(j,i) = y(j,i)-w1(j)
  125    continue
  126 continue
c
      izr = nm
      do 131 l=1,kdo
         ir = l-1
         i2 = 2**ir
         i1 = i2/2
         i3 = i2+i1
         i4 = i2+i2
         irm1 = ir-1
         do 130 i=i4,if,i4
            ipi1 = i+i1
            ipi2 = i+i2
            ipi3 = i+i3
            if (ipi2-izr) 127,128,127
  127       if (i-izr) 130,131,130
  128       call indxc (i,ir,idxc,nc)
            call indxb (ipi2,ir,ip2,np2)
            call indxb (ipi1,irm1,ip1,np1)
            call indxb (ipi3,irm1,ip3,np3)
            call prdct (np2,b(ip2),np1,b(ip1),np3,b(ip3),0,dum,w2,w2,m,
     1                  am,bm,cm,wd,ww,wu)
            call prdct (np1,b(ip1),0,dum,0,dum,nc,cn(idxc),w2,w2,m,am,
     1                  bm,cm,wd,ww,wu)
            do 129 j=1,m
               y(j,i) = y(j,i)-w2(j)
  129       continue
            izr = i
            go to 131
  130    continue
  131 continue
c
c begin back substitution phase
c
  132 do 144 ll=1,k
         l = k-ll+1
         ir = l-1
         irm1 = ir-1
         i2 = 2**ir
         i1 = i2/2
         i4 = i2+i2
         ifd = if-i2
         do 143 i=i2,ifd,i4
            if (i-nm) 133,133,143
  133       imi1 = i-i1
            imi2 = i-i2
            ipi1 = i+i1
            ipi2 = i+i2
            call indxa (i,ir,idxa,na)
            call indxc (i,ir,idxc,nc)
            call indxb (i,ir,iz,nz)
            call indxb (imi1,irm1,im1,nm1)
            call indxb (ipi1,irm1,ip1,np1)
            if (i-i2) 134,134,136
  134       do 135 j=1,m
               w1(j) = 0.
  135       continue
            go to 137
  136       call prdct (nm1,b(im1),0,dum,0,dum,na,an(idxa),y(1,imi2),
     1                  w1,m,am,bm,cm,wd,ww,wu)
  137       if (ipi2-nm) 140,140,138
  138       do 139 j=1,m
               w2(j) = 0.
  139       continue
            go to 141
  140       call prdct (np1,b(ip1),0,dum,0,dum,nc,cn(idxc),y(1,ipi2),
     1                  w2,m,am,bm,cm,wd,ww,wu)
  141       do 142 j=1,m
               w1(j) = y(j,i)+w1(j)+w2(j)
  142       continue
            call prdct (nz,b(iz),nm1,b(im1),np1,b(ip1),0,dum,w1,y(1,i),
     1                  m,am,bm,cm,wd,ww,wu)
  143    continue
  144 continue
      return
      end
      subroutine indxa (i,ir,idxa,na)
      common /cblkt/  npp        ,k          ,eps        ,cnv        ,
     1                nm         ,ncmplx     ,ik
      na = 2**ir
      idxa = i-na+1
      if (i-nm) 102,102,101
  101 na = 0
  102 return
      end
      subroutine indxc (i,ir,idxc,nc)
      common /cblkt/  npp        ,k          ,eps        ,cnv        ,
     1                nm         ,ncmplx     ,ik
      nc = 2**ir
      idxc = i
      if (idxc+nc-1-nm) 102,102,101
  101 nc = 0
  102 return
      end
      subroutine compb (n,ierror,an,bn,cn,b,ah,bh)
c
c     compb computes the roots of the b polynomials using subroutine
c     tevls which is a modification the eispack program tqlrat.
c     ierror is set to 4 if either tevls fails or if a(j+1)*c(j) is
c     less than zero for some j.  ah,bh are temporary work arrays.
c
      dimension       an(1)      ,bn(1)      ,cn(1)      ,b(1)       ,
     1                ah(1)      ,bh(1)
      common /cblkt/  npp        ,k          ,eps        ,cnv        ,
     1                nm         ,ncmplx     ,ik
      eps = epmach(dum)
      bnorm = abs(bn(1))
      do 102 j=2,nm
         bnorm = amax1(bnorm,abs(bn(j)))
         arg = an(j)*cn(j-1)
         if (arg) 119,101,101
  101    b(j) = sign(sqrt(arg),an(j))
  102 continue
      cnv = eps*bnorm
      if = 2**k
      kdo = k-1
      do 108 l=1,kdo
         ir = l-1
         i2 = 2**ir
         i4 = i2+i2
         ipl = i4-1
         ifd = if-i4
         do 107 i=i4,ifd,i4
            call indxb (i,l,ib,nb)
            if (nb) 108,108,103
  103       js = i-ipl
            jf = js+nb-1
            ls = 0
            do 104 j=js,jf
               ls = ls+1
               bh(ls) = bn(j)
               ah(ls) = b(j)
  104       continue
            call tevls (nb,bh,ah,ierror)
            if (ierror) 118,105,118
  105       lh = ib-1
            do 106 j=1,nb
               lh = lh+1
               b(lh) = -bh(j)
  106       continue
  107    continue
  108 continue
      do 109 j=1,nm
         b(j) = -bn(j)
  109 continue
      if (npp) 117,110,117
  110 nmp = nm+1
      nb = nm+nmp
      do 112 j=1,nb
         l1 = mod(j-1,nmp)+1
         l2 = mod(j+nm-1,nmp)+1
         arg = an(l1)*cn(l2)
         if (arg) 119,111,111
  111    bh(j) = sign(sqrt(arg),-an(l1))
         ah(j) = -bn(l1)
  112 continue
      call tevls (nb,ah,bh,ierror)
      if (ierror) 118,113,118
  113 call indxb (if,k-1,j2,lh)
      call indxb (if/2,k-1,j1,lh)
      j2 = j2+1
      lh = j2
      n2m2 = j2+nm+nm-2
  114 d1 = abs(b(j1)-b(j2-1))
      d2 = abs(b(j1)-b(j2))
      d3 = abs(b(j1)-b(j2+1))
      if ((d2 .lt. d1) .and. (d2 .lt. d3)) go to 115
      b(lh) = b(j2)
      j2 = j2+1
      lh = lh+1
      if (j2-n2m2) 114,114,116
  115 j2 = j2+1
      j1 = j1+1
      if (j2-n2m2) 114,114,116
  116 b(lh) = b(n2m2+1)
      call indxb (if,k-1,j1,j2)
      j2 = j1+nmp+nmp
      call ppadd (nm+1,ierror,an,cn,b(j1),b(j1),b(j2))
  117 return
  118 ierror = 4
      return
  119 ierror = 5
      return
      end
      subroutine ppadd (n,ierror,a,c,cbp,bp,bh)
c
c     ppadd computes the eigenvalues of the periodic tridiagonal matrix
c     with coefficients an,bn,cn
c
c n is the order of the bh and bp polynomials
c     on output bp contians the eigenvalues
c cbp is the same as bp except type complex
c bh is used to temporarily store the roots of the b hat polynomial
c which enters through bp
c
      complex         cf         ,cx         ,fsg        ,hsg        ,
     1                dd         ,f          ,fp         ,fpp        ,
     2                cdis       ,r1         ,r2         ,r3         ,
     3                cbp
      dimension       a(1)       ,c(1)       ,bp(1)      ,bh(1)      ,
     1                cbp(1)
      common /cblkt/  npp        ,k          ,eps        ,cnv        ,
     1                nm         ,ncmplx     ,ik
      external        psgf       ,ppspf      ,ppsgf
      scnv = sqrt(cnv)
      iz = n
      izm = iz-1
      izm2 = iz-2
      if (bp(n)-bp(1)) 101,142,103
  101 do 102 j=1,n
         nt = n-j
         bh(j) = bp(nt+1)
  102 continue
      go to 105
  103 do 104 j=1,n
         bh(j) = bp(j)
  104 continue
  105 ncmplx = 0
      modiz = mod(iz,2)
      is = 1
      if (modiz) 106,107,106
  106 if (a(1)) 110,142,107
  107 xl = bh(1)
      db = bh(3)-bh(1)
  108 xl = xl-db
      if (psgf(xl,iz,c,a,bh)) 108,108,109
  109 sgn = -1.
      cbp(1) = cmplx(bsrh(xl,bh(1),iz,c,a,bh,psgf,sgn),0.)
      is = 2
  110 if = iz-1
      if (modiz) 111,112,111
  111 if (a(1)) 112,142,115
  112 xr = bh(iz)
      db = bh(iz)-bh(iz-2)
  113 xr = xr+db
      if (psgf(xr,iz,c,a,bh)) 113,114,114
  114 sgn = 1.
      cbp(iz) = cmplx(bsrh(bh(iz),xr,iz,c,a,bh,psgf,sgn),0.)
      if = iz-2
  115 do 136 ig=is,if,2
         xl = bh(ig)
         xr = bh(ig+1)
         sgn = -1.
         xm = bsrh(xl,xr,iz,c,a,bh,ppspf,sgn)
         psg = psgf(xm,iz,c,a,bh)
         if (abs(psg)-eps) 118,118,116
  116    if (psg*ppsgf(xm,iz,c,a,bh)) 117,118,119
c
c     case of a real zero
c
  117    sgn = 1.
         cbp(ig) = cmplx(bsrh(bh(ig),xm,iz,c,a,bh,psgf,sgn),0.)
         sgn = -1.
         cbp(ig+1) = cmplx(bsrh(xm,bh(ig+1),iz,c,a,bh,psgf,sgn),0.)
         go to 136
c
c     case of a multiple zero
c
  118    cbp(ig) = cmplx(xm,0.)
         cbp(ig+1) = cmplx(xm,0.)
         go to 136
c
c     case of a complex zero
c
  119    it = 0
         icv = 0
         cx = cmplx(xm,0.)
  120    fsg = (1.,0.)
         hsg = (1.,0.)
         fp = (0.,0.)
         fpp = (0.,0.)
         do 121 j=1,iz
            dd = 1./(cx-bh(j))
            fsg = fsg*a(j)*dd
            hsg = hsg*c(j)*dd
            fp = fp+dd
            fpp = fpp-dd*dd
  121    continue
         if (modiz) 123,122,123
  122    f = (1.,0.)-fsg-hsg
         go to 124
  123    f = (1.,0.)+fsg+hsg
  124    i3 = 0
         if (cabs(fp)) 126,126,125
  125    i3 = 1
         r3 = -f/fp
  126    i2 = 0
         if (cabs(fpp)) 132,132,127
  127    i2 = 1
         cdis = csqrt(fp**2-2.*f*fpp)
         r1 = cdis-fp
         r2 = -fp-cdis
         if (cabs(r1)-cabs(r2)) 129,129,128
  128    r1 = r1/fpp
         go to 130
  129    r1 = r2/fpp
  130    r2 = 2.*f/fpp/r1
         if (cabs(r2) .lt. cabs(r1)) r1 = r2
         if (i3) 133,133,131
  131    if (cabs(r3) .lt. cabs(r1)) r1 = r3
         go to 133
  132    r1 = r3
  133    cx = cx+r1
         it = it+1
         if (it .gt. 50) go to 142
         if (cabs(r1) .gt. scnv) go to 120
         if (icv) 134,134,135
  134    icv = 1
         go to 120
  135    cbp(ig) = cx
         cbp(ig+1) = conjg(cx)
  136 continue
      if (cabs(cbp(n))-cabs(cbp(1))) 137,142,139
  137 nhalf = n/2
      do 138 j=1,nhalf
         nt = n-j
         cx = cbp(j)
         cbp(j) = cbp(nt+1)
         cbp(nt+1) = cx
  138 continue
  139 ncmplx = 1
      do 140 j=2,iz
         if (aimag(cbp(j))) 143,140,143
  140 continue
      ncmplx = 0
      do 141 j=2,iz
         bp(j) = real(cbp(j))
  141 continue
      go to 143
  142 ierror = 4
  143 continue
      return
      end
      function ppsgf (x,iz,c,a,bh)
      dimension       a(1)       ,c(1)       ,bh(1)
      sum = 0.
      do 101 j=1,iz
         sum = sum-1./(x-bh(j))**2
  101 continue
      ppsgf = sum
      return
      end
      function ppspf (x,iz,c,a,bh)
      dimension       a(1)       ,c(1)       ,bh(1)
      sum = 0.
      do 101 j=1,iz
         sum = sum+1./(x-bh(j))
  101 continue
      ppspf = sum
      return
      end
      function bsrh (xll,xrr,iz,c,a,bh,f,sgn)
      dimension       a(1)       ,c(1)       ,bh(1)
      common /cblkt/  npp        ,k          ,eps        ,cnv        ,
     1                nm         ,ncmplx     ,ik
      xl = xll
      xr = xrr
      dx = .5*abs(xr-xl)
  101 x = .5*(xl+xr)
      if (sgn*f(x,iz,c,a,bh)) 103,105,102
  102 xr = x
      go to 104
  103 xl = x
  104 dx = .5*dx
      if (dx-cnv) 105,105,101
  105 bsrh = .5*(xl+xr)
      return
      end
      function psgf (x,iz,c,a,bh)
      dimension       a(1)       ,c(1)       ,bh(1)
      fsg = 1.
      hsg = 1.
      do 101 j=1,iz
         dd = 1./(x-bh(j))
         fsg = fsg*a(j)*dd
         hsg = hsg*c(j)*dd
  101 continue
      if (mod(iz,2)) 103,102,103
  102 psgf = 1.-fsg-hsg
      return
  103 psgf = 1.+fsg+hsg
      return
      end
      subroutine tevls (n,d,e2,ierr)
c
      integer         i          ,j          ,l          ,m          ,
     1                n          ,ii         ,l1         ,mml        ,
     2                ierr
      real            d(n)       ,e2(n)
      real            b          ,c          ,f          ,g          ,
     1                h          ,p          ,r          ,s          ,
     2                machep
c
c     real sqrt,abs,sign
c
      common /cblkt/  npp        ,k          ,machep     ,cnv        ,
     1                nm         ,ncmplx     ,ik
c
c     this subroutine is a modification of the eispack subroutine tqlrat
c     algorithm 464, comm. acm 16, 689(1973) by reinsch.
c
c     this subroutine finds the eigenvalues of a symmetric
c     tridiagonal matrix by the rational ql method.
c
c     on input-
c
c        n is the order of the matrix,
c
c        d contains the diagonal elements of the input matrix,
c
c        e2 contains the                subdiagonal elements of the
c          input matrix in its last n-1 positions.  e2(1) is arbitrary.
c
c      on output-
c
c        d contains the eigenvalues in ascending order.  if an
c          error exit is made, the eigenvalues are correct and
c          ordered for indices 1,2,...ierr-1, but may not be
c          the smallest eigenvalues,
c
c        e2 has been destroyed,
c
c        ierr is set to
c          zero       for normal return,
c          j          if the j-th eigenvalue has not been
c                     determined after 30 iterations.
c
c     questions and comments should be directed to b. s. garbow,
c     applied mathematics division, argonne national laboratory
c
c
c     ********** machep is a machine dependent parameter specifying
c                the relative precision of floating point arithmetic.
c
c                **********
c
      ierr = 0
      if (n .eq. 1) go to 115
c
      do 101 i=2,n
         e2(i-1) = e2(i)*e2(i)
  101 continue
c
      f = 0.0
      b = 0.0
      e2(n) = 0.0
c
      do 112 l=1,n
         j = 0
         h = machep*(abs(d(l))+sqrt(e2(l)))
         if (b .gt. h) go to 102
         b = h
         c = b*b
c
c     ********** look for small squared sub-diagonal element **********
c
  102    do 103 m=l,n
            if (e2(m) .le. c) go to 104
c
c     ********** e2(n) is always zero, so there is no exit
c                through the bottom of the loop **********
c
  103    continue
c
  104    if (m .eq. l) go to 108
  105    if (j .eq. 30) go to 114
         j = j+1
c
c     ********** form shift **********
c
         l1 = l+1
         s = sqrt(e2(l))
         g = d(l)
         p = (d(l1)-g)/(2.0*s)
         r = sqrt(p*p+1.0)
         d(l) = s/(p+sign(r,p))
         h = g-d(l)
c
         do 106 i=l1,n
            d(i) = d(i)-h
  106    continue
c
         f = f+h
c
c     ********** rational ql transformation **********
c
         g = d(m)
         if (g .eq. 0.0) g = b
         h = g
         s = 0.0
         mml = m-l
c
c     ********** for i=m-1 step -1 until l do -- **********
c
         do 107 ii=1,mml
            i = m-ii
            p = g*h
            r = p+e2(i)
            e2(i+1) = s*r
            s = e2(i)/r
            d(i+1) = h+s*(h+d(i))
            g = d(i)-e2(i)/g
            if (g .eq. 0.0) g = b
            h = g*p/r
  107    continue
c
         e2(l) = s*g
         d(l) = h
c
c     ********** guard against underflowed h **********
c
         if (h .eq. 0.0) go to 108
         if (abs(e2(l)) .le. abs(c/h)) go to 108
         e2(l) = h*e2(l)
         if (e2(l) .ne. 0.0) go to 105
  108    p = d(l)+f
c
c     ********** order eigenvalues **********
c
         if (l .eq. 1) go to 110
c
c     ********** for i=l step -1 until 2 do -- **********
c
         do 109 ii=2,l
            i = l+2-ii
            if (p .ge. d(i-1)) go to 111
            d(i) = d(i-1)
  109    continue
c
  110    i = 1
  111    d(i) = p
  112 continue
c
      if (abs(d(n)) .ge. abs(d(1))) go to 115
      nhalf = n/2
      do 113 i=1,nhalf
         ntop = n-i
         dhold = d(i)
         d(i) = d(ntop+1)
         d(ntop+1) = dhold
  113 continue
      go to 115
c
c     ********** set error -- no convergence to an
c                eigenvalue after 30 iterations **********
c
  114 ierr = l
  115 return
c
c     ********** last card of tqlrat **********
c
      end
      subroutine indxb (i,ir,idx,idp)
c
c b(idx) is the location of the first root of the b(i,ir) polynomial
c
      common /cblkt/  npp        ,k          ,eps        ,cnv        ,
     1                nm         ,ncmplx     ,ik
      idp = 0
      if (ir) 107,101,103
  101 if (i-nm) 102,102,107
  102 idx = i
      idp = 1
      return
  103 izh = 2**ir
      id = i-izh-izh
      idx = id+id+(ir-1)*ik+ir+(ik-i)/izh+4
      ipl = izh-1
      idp = izh+izh-1
      if (i-ipl-nm) 105,105,104
  104 idp = 0
      return
  105 if (i+ipl-nm) 107,107,106
  106 idp = nm+ipl-i+1
  107 return
      end
      function epmach (dum)
c
c     this program computes an approximate machiine epsilon (accuracy)
c
      common /value/  v
      eps = 1.
  101 eps = eps/10.
      call store (eps+1.)
      if (v-1.) 102,102,101
  102 epmach = 100.*eps
      return
      end
      subroutine store (x)
      common /value/  v
      v = x
      return
      end
