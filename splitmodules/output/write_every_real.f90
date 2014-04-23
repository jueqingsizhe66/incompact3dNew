  subroutine write_every_real(ipencil,var,iskip,jskip,kskip, &
       filename, from1)

    implicit none

    integer, intent(IN) :: ipencil !(x-pencil=1; y-pencil=2; z-pencil=3)
    real(mytype), dimension(:,:,:), intent(IN) :: var
    integer, intent(IN) :: iskip,jskip,kskip
    character(len=*), intent(IN) :: filename
    logical, intent(IN) :: from1  ! .true.  - save 1,n+1,2n+1...
                                  ! .false. - save n,2n,3n...

    real(mytype), allocatable, dimension(:,:,:) :: wk, wk2
    integer(kind=MPI_OFFSET_KIND) :: filesize, disp
    integer, dimension(3) :: sizes, subsizes, starts
    integer :: i,j,k, ierror, newtype, fh, key,color,newcomm, data_type
    integer, dimension(3) :: xsz,ysz,zsz,xst,yst,zst,xen,yen,zen,skip

    data_type = real_type

#include "io_write_every.f90"

    return
  end subroutine write_every_real
