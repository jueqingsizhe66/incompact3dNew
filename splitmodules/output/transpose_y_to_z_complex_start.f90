  subroutine transpose_y_to_z_complex_start(handle, src, dst, sbuf, &
       rbuf, opt_decomp)

    implicit none

    integer :: handle
    complex(mytype), dimension(:,:,:) :: src, dst, sbuf, rbuf
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
    call mem_split_yz_complex(src, s1, s2, s3, sbuf, dims(2), &
         decomp%y2dist, decomp)

#ifdef EVEN
    call NBC_IALLTOALL(sbuf, decomp%y2count, &
         complex_type, rbuf, decomp%z2count, &
         complex_type, DECOMP_2D_COMM_ROW, handle, ierror)
#else
    call NBC_IALLTOALLV(sbuf, decomp%y2cnts, decomp%y2disp, &
         complex_type, rbuf,decomp%z2cnts, decomp%z2disp, &
         complex_type, DECOMP_2D_COMM_ROW, handle, ierror)
#endif

    return
  end subroutine transpose_y_to_z_complex_start
