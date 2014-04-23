  subroutine transpose_y_to_z_complex(src, dst, opt_decomp)

    implicit none

    complex(mytype), dimension(:,:,:), intent(IN) :: src
    complex(mytype), dimension(:,:,:), intent(OUT) :: dst
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    TYPE(DECOMP_INFO) :: decomp

#ifdef SHM
    complex(mytype) :: work1(*), work2(*)
    POINTER  (work1_p, work1), (work2_p, work2)  ! Cray pointers
#endif

    integer :: s1,s2,s3,d1,d2,d3
    integer :: ierror

    if (present(opt_decomp)) then
       decomp = opt_decomp
    else
       decomp = decomp_main
    end if

    s1 = SIZE(src,1)
    s2 = SIZE(src,2)
    s3 = SIZE(src,3)
    d1 = SIZE(dst,1)
    d2 = SIZE(dst,2)
    d3 = SIZE(dst,3)

    ! rearrange source array as send buffer
#ifdef SHM
    work1_p = decomp%ROW_INFO%SND_P_c
    call mem_split_yz_complex(src, s1, s2, s3, work1, dims(2), &
         decomp%y2dist, decomp)
#else
    call mem_split_yz_complex(src, s1, s2, s3, work1_c, dims(2), &
         decomp%y2dist, decomp)
#endif

    ! define receive buffer
#ifdef SHM
    work2_p = decomp%ROW_INFO%RCV_P_c
    call MPI_BARRIER(decomp%ROW_INFO%CORE_COMM, ierror)
#endif

#ifdef SHM
    if (decomp%ROW_INFO%CORE_ME==1) THEN
       call MPI_ALLTOALLV(work1, decomp%y2cnts_s, decomp%y2disp_s, &
            complex_type, work2, decomp%z2cnts_s, decomp%z2disp_s, &
            complex_type, decomp%ROW_INFO%SMP_COMM, ierror)
    end if
#else
#ifdef EVEN
    if (decomp%even) then
       call MPI_ALLTOALL(work1_c, decomp%y2count, &
            complex_type, dst, decomp%z2count, &
            complex_type, DECOMP_2D_COMM_ROW, ierror)
    else
       call MPI_ALLTOALL(work1_c, decomp%y2count, &
            complex_type, work2_c, decomp%z2count, &
            complex_type, DECOMP_2D_COMM_ROW, ierror)
    end if
#else
    call MPI_ALLTOALLV(work1_c, decomp%y2cnts, decomp%y2disp, &
         complex_type, dst, decomp%z2cnts, decomp%z2disp, &
         complex_type, DECOMP_2D_COMM_ROW, ierror)
#endif
#endif

    ! rearrange receive buffer
#ifdef SHM
    call MPI_BARRIER(decomp%ROW_INFO%CORE_COMM, ierror)
    call mem_merge_yz_complex(work2, d1, d2, d3, dst, dims(2), &
         decomp%z2dist, decomp)
#else
#ifdef EVEN
    if (.not. decomp%even) then
       call mem_merge_yz_complex(work2_c, d1, d2, d3, dst, dims(2), &
         decomp%z2dist, decomp)
    end if
#else
    ! note the receive buffer is already in natural (i,j,k) order
    ! so no merge operation needed
#endif
#endif

    return
  end subroutine transpose_y_to_z_complex
