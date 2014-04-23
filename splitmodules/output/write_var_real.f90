  subroutine write_var_real(fh,disp,ipencil,var,opt_decomp)

    implicit none

    integer, intent(IN) :: fh
    integer(KIND=MPI_OFFSET_KIND), intent(INOUT) :: disp
    integer, intent(IN) :: ipencil
    real(mytype), dimension(:,:,:), intent(IN) :: var
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    TYPE(DECOMP_INFO) :: decomp
    integer, dimension(3) :: sizes, subsizes, starts
    integer :: ierror, newtype, data_type

    data_type = real_type

#include "io_write_var.f90"

    return
  end subroutine write_var_real
