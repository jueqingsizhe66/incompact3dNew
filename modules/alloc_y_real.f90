  subroutine alloc_y_real(var, opt_decomp, opt_global)

    implicit none

    real(mytype), allocatable, dimension(:,:,:) :: var
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp
    logical, intent(IN), optional :: opt_global

    TYPE(DECOMP_INFO) :: decomp
    logical :: global
    integer :: alloc_stat, errorcode

    if (present(opt_decomp)) then
       decomp = opt_decomp
    else
       decomp = decomp_main
    end if

    if (present(opt_global)) then
       global = opt_global
    else
       global = .false.
    end if

    if (global) then
       allocate(var(decomp%yst(1):decomp%yen(1), &
            decomp%yst(2):decomp%yen(2), decomp%yst(3):decomp%yen(3)), &
            stat=alloc_stat)
    else
       allocate(var(decomp%ysz(1),decomp%ysz(2),decomp%ysz(3)), &
            stat=alloc_stat)
    end if

    if (alloc_stat /= 0) then
       errorcode = 8
       call decomp_2d_abort(errorcode, &
            'Memory allocation failed when creating new arrays')
    end if

    return
  end subroutine alloc_y_real
