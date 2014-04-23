  subroutine transpose_y_to_x_real_start(handle, src, dst, sbuf, rbuf, &
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
    call mem_split_yx_real(src, s1, s2, s3, sbuf, dims(1), &
         decomp%y1dist, decomp)

#ifdef EVEN
    call NBC_IALLTOALL(sbuf, decomp%y1count, real_type, &
         rbuf, decomp%x1count, real_type, &
         DECOMP_2D_COMM_COL, handle, ierror)
#else
    call NBC_IALLTOALLV(sbuf, decomp%y1cnts, decomp%y1disp, real_type, &
         rbuf, decomp%x1cnts, decomp%x1disp, real_type, &
         DECOMP_2D_COMM_COL, handle, ierror)
#endif

    return
  end subroutine transpose_y_to_x_real_start
