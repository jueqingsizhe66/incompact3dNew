  subroutine transpose_z_to_y_complex_start(handle, src, dst, sbuf, &
       rbuf, opt_decomp)

    implicit none

    integer :: handle
    complex(mytype), dimension(:,:,:) :: src, dst, sbuf, rbuf
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    TYPE(DECOMP_INFO) :: decomp

    integer :: ierror

    if (present(opt_decomp)) then
       decomp = opt_decomp
    else
       decomp = decomp_main
    end if

    sbuf = src

#ifdef EVEN
    call NBC_IALLTOALL(sbuf, decomp%z2count, &
         complex_type, rbuf, decomp%y2count, &
         complex_type, DECOMP_2D_COMM_ROW, handle, ierror)
#else
    call NBC_IALLTOALLV(sbuf, decomp%z2cnts, decomp%z2disp, &
         complex_type, rbuf, decomp%y2cnts, decomp%y2disp, &
         complex_type, DECOMP_2D_COMM_ROW, handle, ierror)
#endif

    return
  end subroutine transpose_z_to_y_complex_start
