!BE CAREFUL WITH FORTRAN COMPILER

!USE OPTION ifort -assume byterecl

program visu

  implicit none

  integer, parameter :: nx=17 , ny=16, nz=16, nlength=100
  real(8),dimension(nx,ny,nz,nlength) :: ux
  integer :: i,j,k,count,itime,ii

  OPEN(10,FILE='probe_phi',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
 COUNT = 1
 DO II=1,nlength
    DO K=1,nz
       DO J=1,ny
          DO I=1,nx
             READ(10,REC=COUNT) ux(I,J,K,II)
             COUNT = COUNT + 1
          ENDDO
       ENDDO
        print *,ii,k,ux(nx/2,ny/2,k,ii)
    ENDDO
!    print *,ii,ux(nx/2,ny/2,nz/2,ii)
 ENDDO
  CLOSE(10)


  open(20,file='tampon.avs',form='formatted')
  do ii=1,nlength
     write(20,*) ii,ux(nx/2,ny/2,nz/2,ii)
  enddo
  close(20)

    end program visu
   
