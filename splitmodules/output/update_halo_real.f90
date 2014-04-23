  subroutine update_halo_real(in, out, level, opt_decomp, opt_global)

    implicit none

    integer, intent(IN) :: level      ! levels of halo cells required
    real(mytype), dimension(:,:,:), intent(IN) :: in
    real(mytype), allocatable, dimension(:,:,:), intent(OUT) :: out
    TYPE(DECOMP_INFO), optional :: opt_decomp
    logical, optional :: opt_global

    TYPE(DECOMP_INFO) :: decomp
    logical :: global

    ! starting/ending index of array with halo cells
    integer :: xs, ys, zs, xe, ye, ze

    integer :: i, j, k, s1, s2, s3, ierror
    integer :: data_type

    integer :: icount, ilength, ijump
    integer :: halo12, halo21, halo31, halo32
    integer, dimension(4) :: requests
    integer, dimension(MPI_STATUS_SIZE,4) :: status
    integer :: tag_e, tag_w, tag_n, tag_s, tag_t, tag_b

    data_type = real_type

#include "halo_common.f90"

    return
  end subroutine update_halo_real
