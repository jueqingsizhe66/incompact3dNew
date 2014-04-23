module glassman

  use decomp_2d, only : mytype

  implicit none

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Following is a FFT implementation based on algorithm proposed by
  ! Glassman, a general FFT algorithm supporting arbitrary input length.
  !
  ! W. E. Ferguson, Jr., "A simple derivation of Glassman general-n fast
  ! Fourier transform," Comput. and Math. with Appls., vol. 8, no. 6, pp.
  ! 401-411, 1982.
  !
  ! Original implemtation online at http://www.jjj.de/fft/fftpage.html
  !
  ! Updated
  !  -  to handle double-precision as well
  !  -  unnecessary scaling code removed
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE SPCFFT(U,N,ISIGN,WORK)

    IMPLICIT NONE

    LOGICAL :: INU
    INTEGER :: A,B,C,N,I,ISIGN
    COMPLEX(mytype) :: U(*),WORK(*)

    A = 1
    B = N
    C = 1
    INU = .TRUE.

    DO WHILE ( B .GT. 1 )
       A = C * A
       C = 2
       DO WHILE ( MOD(B,C) .NE. 0 )
          C = C + 1
       END DO
       B = B / C
       IF ( INU ) THEN
          CALL SPCPFT (A,B,C,U,WORK,ISIGN)
       ELSE
          CALL SPCPFT (A,B,C,WORK,U,ISIGN)
       END IF
       INU = ( .NOT. INU )
    END DO

    IF ( .NOT. INU ) THEN
       DO I = 1, N
          U(I) = WORK(I)
       END DO
    END IF

    RETURN
  END SUBROUTINE SPCFFT
