program visu

  implicit none

  integer, parameter :: nx=17 , ny=17, nz=17
  real(4),dimension(nx,ny,nz) :: ux
  integer :: i,j,k,count
  real(4) :: x,y

  do k=1,nz
  do j=1,ny
  do i=1,nx
     ux(i,j,k)=i+nx*(j-1)+nx*ny*(k-1)
  enddo
  enddo
  enddo

  OPEN(10,FILE='epsilon.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=4)
  COUNT = 1
  DO K=1,nz
     DO J=1,ny
        DO I=1,nx
           WRITE(10,REC=COUNT) ux(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
  ENDDO
  CLOSE(10)


    end program visu
