  subroutine read_one_real(ipencil,var,filename,opt_decomp)

    implicit none

    integer, intent(IN) :: ipencil
    real(mytype), dimension(:,:,:), intent(INOUT) :: var
    character(len=*), intent(IN) :: filename
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    TYPE(DECOMP_INFO) :: decomp
    integer(kind=MPI_OFFSET_KIND) :: disp
    integer, dimension(3) :: sizes, subsizes, starts
    integer :: ierror, newtype, fh, data_type

    data_type = real_type

#include "io_read_one.f90"

    return
  end subroutine read_one_real
