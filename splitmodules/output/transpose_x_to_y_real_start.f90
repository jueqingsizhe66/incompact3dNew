  subroutine transpose_x_to_y_real_start(handle, src, dst, sbuf, rbuf, &
       opt_decomp)

    implicit none

    integer :: handle
    real(mytype), dimension(:,:,:) :: src, dst, sbuf, rbuf
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    TYPE(DECOMP_INFO) :: decomp

    integer :: s1,s2,s3
    integer :: ierror

    if (present(opt_decomp)) then
       decomp = opt_decomp
    else
       decomp = decomp_main
    end if

    s1 = SIZE(src,1)
    s2 = SIZE(src,2)
    s3 = SIZE(src,3)

    ! rearrange source array as send buffer
    call mem_split_xy_real(src, s1, s2, s3, sbuf, dims(1), &
         decomp%x1dist, decomp)

#ifdef EVEN
    call NBC_IALLTOALL(sbuf, decomp%x1count, real_type, &
         rbuf, decomp%y1count, real_type, &
         DECOMP_2D_COMM_COL, handle, ierror)
#else
    call NBC_IALLTOALLV(sbuf, decomp%x1cnts, decomp%x1disp, real_type, &
         rbuf, decomp%y1cnts, decomp%y1disp, real_type, &
         DECOMP_2D_COMM_COL, handle, ierror)
#endif

    return
  end subroutine transpose_x_to_y_real_start
