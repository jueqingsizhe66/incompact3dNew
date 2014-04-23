  subroutine c2c_1m_y_plan(plan1, decomp, isign)

    implicit none

    integer*8, intent(OUT) :: plan1
    TYPE(DECOMP_INFO), intent(IN) :: decomp
    integer, intent(IN) :: isign

    complex(mytype), allocatable, dimension(:,:) :: a1

    ! Due to memory pattern of 3D arrays, 1D FFTs along Y have to be
    ! done one Z-plane at a time. So plan for 2D data sets here.

    allocate(a1(decomp%ysz(1),decomp%ysz(2)))

#ifdef DOUBLE_PREC
    call dfftw_plan_many_dft(plan1, 1, decomp%ysz(2), decomp%ysz(1), &
         a1, decomp%ysz(2), decomp%ysz(1), 1, a1, decomp%ysz(2), &
         decomp%ysz(1), 1, isign, plan_type)
#else
    call sfftw_plan_many_dft(plan1, 1, decomp%ysz(2), decomp%ysz(1), &
         a1, decomp%ysz(2), decomp%ysz(1), 1, a1, decomp%ysz(2), &
         decomp%ysz(1), 1, isign, plan_type)
#endif

    deallocate(a1)

    return
  end subroutine c2c_1m_y_plan
