subroutine channel (ux)
!
!********************************************************************

USE decomp_2d
USE decomp_2d_poisson
USE variables
USE param
USE var
USE MPI

implicit none

real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ux

integer :: j,i,k,code
real(mytype) :: can,ut3,ut,ut4

ut3=0.
do k=1,ysize(3)
do i=1,ysize(1)
   ut=0.
   do j=1,ny-1
      if (istret.ne.0) ut=ut+(yp(j+1)-yp(j))*(ux(i,j+1,k)-0.5*(ux(i,j+1,k)-ux(i,j,k)))
      if (istret.eq.0) ut=ut+(yly/(ny-1))*(ux(i,j+1,k)-0.5*(ux(i,j+1,k)-ux(i,j,k)))
   enddo
   ut=ut/yly
   ut3=ut3+ut
enddo
enddo
ut3=ut3/ysize(1)/ysize(3)

call MPI_ALLREDUCE(ut3,ut4,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
ut4=ut4/nproc

!can=-2.*xnu*gdt(itr) ! Poisseuille
can=-(2./3.-ut4) ! constant flow rate

if (nrank==0) print *,nrank,'UT',ut4,can

do k=1,ysize(3)
do i=1,ysize(1)
do j=2,ny-1
   ux(i,j,k)=-can+ux(i,j,k)
enddo
enddo
enddo

return
end subroutine channel
