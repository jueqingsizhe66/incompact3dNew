subroutine test_speed_min_max(ux,uy,uz)
!
!*******************************************************************

USE decomp_2d
USE decomp_2d_poisson
USE variables
USE param
USE var
USE MPI


implicit none

integer :: i,j,k
real(mytype) :: uxmax,uymax,uzmax,uxmin,uymin,uzmin,cfl
real(mytype) :: uxmax1,uymax1,uzmax1,uxmin1,uymin1,uzmin1

real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
integer :: code

uxmax=-1609.
uymax=-1609.
uzmax=-1609.
uxmin=1609.
uymin=1609.
uzmin=1609.
do k=1,xsize(3)
do j=1,xsize(2)
do i=1,xsize(1)
   if (ux(i,j,k).gt.uxmax) uxmax=ux(i,j,k)
   if (uy(i,j,k).gt.uymax) uymax=uy(i,j,k)
   if (uz(i,j,k).gt.uzmax) uzmax=uz(i,j,k)
   if (ux(i,j,k).lt.uxmin) uxmin=ux(i,j,k)
   if (uy(i,j,k).lt.uymin) uymin=uy(i,j,k)
   if (uz(i,j,k).lt.uzmin) uzmin=uz(i,j,k)
enddo
enddo
enddo


call MPI_REDUCE(uxmax,uxmax1,1,real_type,MPI_MAX,0,MPI_COMM_WORLD,code)
call MPI_REDUCE(uymax,uymax1,1,real_type,MPI_MAX,0,MPI_COMM_WORLD,code)
call MPI_REDUCE(uzmax,uzmax1,1,real_type,MPI_MAX,0,MPI_COMM_WORLD,code)
call MPI_REDUCE(uxmin,uxmin1,1,real_type,MPI_MIN,0,MPI_COMM_WORLD,code)
call MPI_REDUCE(uymin,uymin1,1,real_type,MPI_MIN,0,MPI_COMM_WORLD,code)
call MPI_REDUCE(uzmin,uzmin1,1,real_type,MPI_MIN,0,MPI_COMM_WORLD,code)
if (nrank==0) then
   print *,'U,V,W max=',uxmax1,uymax1,uzmax1
   print *,'U,V,W min=',uxmin1,uymin1,uzmin1
endif


return
end subroutine test_speed_min_max
