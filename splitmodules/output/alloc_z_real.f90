  subroutine alloc_z_real(var, opt_decomp, opt_global)

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
       allocate(var(decomp%zst(1):decomp%zen(1), &
            decomp%zst(2):decomp%zen(2), decomp%zst(3):decomp%zen(3)), &
            stat=alloc_stat)
    else
       allocate(var(decomp%zsz(1),decomp%zsz(2),decomp%zsz(3)), &
            stat=alloc_stat)
    end if

    if (alloc_stat /= 0) then
       errorcode = 8
       call decomp_2d_abort(errorcode, &
            'Memory allocation failed when creating new arrays')
    end if

    return
  end subroutine alloc_z_real
