  subroutine transpose_y_to_x_real(src, dst, opt_decomp)

    implicit none

    real(mytype), dimension(:,:,:), intent(IN) :: src
    real(mytype), dimension(:,:,:), intent(OUT) :: dst
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    TYPE(DECOMP_INFO) :: decomp

#ifdef SHM
    real(mytype) :: work1(*), work2(*)
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
    work1_p = decomp%COL_INFO%SND_P
    call mem_split_yx_real(src, s1, s2, s3, work1, dims(1), &
         decomp%y1dist, decomp)
#else
    call mem_split_yx_real(src, s1, s2, s3, work1_r, dims(1), &
         decomp%y1dist, decomp)
#endif

    ! define receive buffer
#ifdef SHM
    work2_p = decomp%COL_INFO%RCV_P
    call MPI_BARRIER(decomp%COL_INFO%CORE_COMM, ierror)
#endif

    ! transpose using MPI_ALLTOALL(V)
#ifdef SHM
    if (decomp%COL_INFO%CORE_ME==1) THEN
       call MPI_ALLTOALLV(work1, decomp%y1cnts_s, decomp%y1disp_s, &
            real_type, work2, decomp%x1cnts_s, decomp%x1disp_s, &
            real_type, decomp%COL_INFO%SMP_COMM, ierror)
    end if
#else
#ifdef EVEN
    call MPI_ALLTOALL(work1_r, decomp%y1count, &
         real_type, work2_r, decomp%x1count, &
         real_type, DECOMP_2D_COMM_COL, ierror)
#else
    call MPI_ALLTOALLV(work1_r, decomp%y1cnts, decomp%y1disp, &
         real_type, work2_r, decomp%x1cnts, decomp%x1disp, &
         real_type, DECOMP_2D_COMM_COL, ierror)
#endif
#endif

    ! rearrange receive buffer
#ifdef SHM
    call MPI_BARRIER(decomp%COL_INFO%CORE_COMM, ierror)
    call mem_merge_yx_real(work2, d1, d2, d3, dst, dims(1), &
         decomp%x1dist, decomp)
#else
    call mem_merge_yx_real(work2_r, d1, d2, d3, dst, dims(1), &
         decomp%x1dist, decomp)
#endif

    return
  end subroutine transpose_y_to_x_real
