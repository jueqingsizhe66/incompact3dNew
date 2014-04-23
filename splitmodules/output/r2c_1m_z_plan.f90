  subroutine r2c_1m_z_plan(plan1, decomp_ph, decomp_sp)

    implicit none

    integer*8, intent(OUT) :: plan1
    TYPE(DECOMP_INFO), intent(IN) :: decomp_ph
    TYPE(DECOMP_INFO), intent(IN) :: decomp_sp

    real(mytype), allocatable, dimension(:,:,:) :: a1
    complex(mytype), allocatable, dimension(:,:,:) :: a2

    allocate(a1(decomp_ph%zsz(1),decomp_ph%zsz(2),decomp_ph%zsz(3)))
    allocate(a2(decomp_sp%zsz(1),decomp_sp%zsz(2),decomp_sp%zsz(3)))
#ifdef DOUBLE_PREC
    call dfftw_plan_many_dft_r2c(plan1, 1, decomp_ph%zsz(3), &
         decomp_ph%zsz(1)*decomp_ph%zsz(2), a1, decomp_ph%zsz(3), &
         decomp_ph%zsz(1)*decomp_ph%zsz(2), 1, a2, decomp_sp%zsz(3), &
         decomp_sp%zsz(1)*decomp_sp%zsz(2), 1, plan_type)
#else
    call sfftw_plan_many_dft_r2c(plan1, 1, decomp_ph%zsz(3), &
         decomp_ph%zsz(1)*decomp_ph%zsz(2), a1, decomp_ph%zsz(3), &
         decomp_ph%zsz(1)*decomp_ph%zsz(2), 1, a2, decomp_sp%zsz(3), &
         decomp_sp%zsz(1)*decomp_sp%zsz(2), 1, plan_type)
#endif
    deallocate(a1,a2)

    return
  end subroutine r2c_1m_z_plan
