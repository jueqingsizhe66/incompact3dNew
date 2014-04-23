  subroutine transpose_y_to_z_real_wait(handle, src, dst, sbuf, rbuf, &
       opt_decomp)

    implicit none

    integer :: handle
    real(mytype), dimension(:,:,:) :: src, dst, sbuf, rbuf
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    TYPE(DECOMP_INFO) :: decomp

    integer :: ierror

    if (present(opt_decomp)) then
       decomp = opt_decomp
    else
       decomp = decomp_main
    end if

    call NBC_WAIT(handle, ierror)

    dst = rbuf

    return
  end subroutine transpose_y_to_z_real_wait
