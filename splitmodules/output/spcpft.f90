  SUBROUTINE SPCPFT( A, B, C, UIN, UOUT, ISIGN )

    IMPLICIT NONE

    INTEGER :: ISIGN,A,B,C,IA,IB,IC,JCR,JC

    DOUBLE PRECISION :: ANGLE

    COMPLEX(mytype) :: UIN(B,C,A),UOUT(B,A,C),DELTA,OMEGA,SUM

    ANGLE = 6.28318530717958_mytype / REAL( A * C, kind=mytype )
    OMEGA = CMPLX( 1.0, 0.0, kind=mytype )

    IF( ISIGN .EQ. 1 ) THEN
       DELTA = CMPLX( DCOS(ANGLE), DSIN(ANGLE), kind=mytype )
    ELSE
       DELTA = CMPLX( DCOS(ANGLE), -DSIN(ANGLE), kind=mytype )
    END IF

    DO IC = 1, C
       DO IA = 1, A
          DO IB = 1, B
             SUM = UIN( IB, C, IA )
             DO JCR = 2, C
                JC = C + 1 - JCR
                SUM = UIN( IB, JC, IA ) + OMEGA * SUM
             END DO
             UOUT( IB, IA, IC ) = SUM
          END DO
          OMEGA = DELTA * OMEGA
       END DO
    END DO

    RETURN
  END SUBROUTINE SPCPFT
