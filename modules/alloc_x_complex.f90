  subroutine alloc_x_complex(var, opt_decomp, opt_global)

    implicit none

    complex(mytype), allocatable, dimension(:,:,:) :: var
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
       allocate(var(decomp%xst(1):decomp%xen(1), &
            decomp%xst(2):decomp%xen(2), decomp%xst(3):decomp%xen(3)), &
            stat=alloc_stat)
    else
       allocate(var(decomp%xsz(1),decomp%xsz(2),decomp%xsz(3)), &
            stat=alloc_stat)
    end if

    if (alloc_stat /= 0) then
       errorcode = 8
       call decomp_2d_abort(errorcode, &
            'Memory allocation failed when creating new arrays')
    end if

    return
  end subroutine alloc_x_complex
