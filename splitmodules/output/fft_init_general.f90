  subroutine fft_init_general(pencil, nx, ny, nz)

    implicit none

    integer, intent(IN) :: pencil
    integer, intent(IN) :: nx, ny, nz

    logical, dimension(2) :: dummy_periods
    integer, dimension(2) :: dummy_coords
    integer :: status, errorcode, ierror

    if (initialised) then
       errorcode = 4
       call decomp_2d_abort(errorcode, &
            'FFT library should only be initialised once')
    end if

    format = pencil
    nx_fft = nx
    ny_fft = ny
    nz_fft = nz

    ! determine the processor grid in use
    call MPI_CART_GET(DECOMP_2D_COMM_CART_X, 2, &
         dims, dummy_periods, dummy_coords, ierror)

    ! for c2r/r2c interface:
    ! if in physical space, a real array is of size: nx*ny*nz
    ! in spectral space, the complex array is of size:
    !         (nx/2+1)*ny*nz, if PHYSICAL_IN_X
    !      or nx*ny*(nz/2+1), if PHYSICAL_IN_Z

    call decomp_info_init(nx, ny, nz, ph)
    if (format==PHYSICAL_IN_X) then
       call decomp_info_init(nx/2+1, ny, nz, sp)
    else if (format==PHYSICAL_IN_Z) then
       call decomp_info_init(nx, ny, nz/2+1, sp)
    end if

    allocate(wk2_c2c(ph%ysz(1),ph%ysz(2),ph%ysz(3)), STAT=status)
    allocate(wk2_r2c(sp%ysz(1),sp%ysz(2),sp%ysz(3)), STAT=status)
    if (format==PHYSICAL_IN_X) then
       allocate(wk13(sp%xsz(1),sp%xsz(2),sp%xsz(3)), STAT=status)
    else if (format==PHYSICAL_IN_Z) then
       allocate(wk13(sp%zsz(1),sp%zsz(2),sp%zsz(3)), STAT=status)
    end if
    if (status /= 0) then
       errorcode = 3
       call decomp_2d_abort(errorcode, &
            'Out of memory when initialising FFT')
    end if

    call init_fft_engine

    initialised = .true.

    return
  end subroutine fft_init_general
