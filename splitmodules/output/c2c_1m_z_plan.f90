  subroutine c2c_1m_z_plan(plan1, decomp, isign)

    implicit none

    integer*8, intent(OUT) :: plan1
    TYPE(DECOMP_INFO), intent(IN) :: decomp
    integer, intent(IN) :: isign

    complex(mytype), allocatable, dimension(:,:,:) :: a1

    allocate(a1(decomp%zsz(1),decomp%zsz(2),decomp%zsz(3)))

#ifdef DOUBLE_PREC
    call dfftw_plan_many_dft(plan1, 1, decomp%zsz(3), &
         decomp%zsz(1)*decomp%zsz(2), a1, decomp%zsz(3), &
         decomp%zsz(1)*decomp%zsz(2), 1, a1, decomp%zsz(3), &
         decomp%zsz(1)*decomp%zsz(2), 1, isign, plan_type)
#else
    call sfftw_plan_many_dft(plan1, 1, decomp%zsz(3), &
         decomp%zsz(1)*decomp%zsz(2), a1, decomp%zsz(3), &
         decomp%zsz(1)*decomp%zsz(2), 1, a1, decomp%zsz(3), &
         decomp%zsz(1)*decomp%zsz(2), 1, isign, plan_type)
#endif

    deallocate(a1)

    return
  end subroutine c2c_1m_z_plan
