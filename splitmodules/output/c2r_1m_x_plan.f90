  subroutine c2r_1m_x_plan(plan1, decomp_sp, decomp_ph)

    implicit none

    integer*8, intent(OUT) :: plan1
    TYPE(DECOMP_INFO), intent(IN) :: decomp_sp
    TYPE(DECOMP_INFO), intent(IN) :: decomp_ph

    complex(mytype), allocatable, dimension(:,:,:) :: a1
    real(mytype), allocatable, dimension(:,:,:) :: a2

    allocate(a1(decomp_sp%xsz(1),decomp_sp%xsz(2),decomp_sp%xsz(3)))
    allocate(a2(decomp_ph%xsz(1),decomp_ph%xsz(2),decomp_ph%xsz(3)))
#ifdef DOUBLE_PREC
    call dfftw_plan_many_dft_c2r(plan1, 1, decomp_ph%xsz(1), &
         decomp_ph%xsz(2)*decomp_ph%xsz(3), a1, decomp_sp%xsz(1), 1, &
         decomp_sp%xsz(1), a2, decomp_ph%xsz(1), 1, decomp_ph%xsz(1), &
         plan_type)
#else
    call sfftw_plan_many_dft_c2r(plan1, 1, decomp_ph%xsz(1), &
         decomp_ph%xsz(2)*decomp_ph%xsz(3), a1, decomp_sp%xsz(1), 1, &
         decomp_sp%xsz(1), a2, decomp_ph%xsz(1), 1, decomp_ph%xsz(1), &
         plan_type)
#endif
    deallocate(a1,a2)

    return
  end subroutine c2r_1m_x_plan
