  subroutine write_one_complex(ipencil,var,filename,opt_decomp)

    implicit none

    integer, intent(IN) :: ipencil
    complex(mytype), dimension(:,:,:), intent(IN) :: var
    character(len=*), intent(IN) :: filename
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    TYPE(DECOMP_INFO) :: decomp
    integer(kind=MPI_OFFSET_KIND) :: filesize, disp
    integer, dimension(3) :: sizes, subsizes, starts
    integer :: ierror, newtype, fh, data_type, info, gs

    data_type = complex_type

#include "io_write_one.f90"

    return
  end subroutine write_one_complex
