  subroutine transpose_y_to_x_real_wait(handle, src, dst, sbuf, rbuf, &
       opt_decomp)

    implicit none

    integer :: handle
    real(mytype), dimension(:,:,:) :: src, dst, sbuf, rbuf
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    TYPE(DECOMP_INFO) :: decomp

    integer :: d1,d2,d3
    integer :: ierror

    if (present(opt_decomp)) then
       decomp = opt_decomp
    else
       decomp = decomp_main
    end if

    d1 = SIZE(dst,1)
    d2 = SIZE(dst,2)
    d3 = SIZE(dst,3)

    call NBC_WAIT(handle, ierror)

    ! rearrange receive buffer
    call mem_merge_yx_real(rbuf, d1, d2, d3, dst, dims(1), &
         decomp%x1dist, decomp)

    return
  end subroutine transpose_y_to_x_real_wait
