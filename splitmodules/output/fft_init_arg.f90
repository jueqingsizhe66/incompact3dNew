  subroutine fft_init_arg(pencil)     ! allow to handle Z-pencil input

    implicit none

    integer, intent(IN) :: pencil

    call fft_init_general(pencil, nx_global, ny_global, nz_global)

    return
  end subroutine fft_init_arg
