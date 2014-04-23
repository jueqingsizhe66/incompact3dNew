program verif

  implicit none

  integer, parameter :: nx1=129 , ny1=33, nz1=16
  real(4),dimension(nx1,ny1,nz1) :: umean,uumean
  integer :: i,j,k,count
  real(4) :: x,y,xitime

 xitime=1000.

 OPEN(11,FILE='umean.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=4, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) umean(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
  ENDDO
  CLOSE(11)

OPEN(11,FILE='uumean.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=4, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) uumean(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
  ENDDO
  CLOSE(11)

umean=umean/xitime
uumean=uumean/xitime

uumean=uumean-umean*umean




 ! DO K=1,nx1
 ! DO J=1,ny1
 ! DO I=1,nx1
 !    print *,i,umean(i,j,k),ux(4*i-3,4*j-3,4*k-3),umean(i,j,k)-ux(4*i-3,4*j-3,4*k-3)
 ! ENDDO
 ! ENDDO
 ! ENDDO

  open(10,file='test1.dat',status='unknown',form='formatted')
  do j=1,ny1
     write(10,*) (j-1.)*12./32,uumean(38,j,1),uumean(45,j,1),uumean(50,j,1)
  enddo
  close(10)



    end program verif
