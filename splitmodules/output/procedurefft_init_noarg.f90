     module procedure fft_init_noarg
     module procedure fft_init_arg
     module procedure fft_init_general
  end interface

  interface decomp_2d_fft_3d
     module procedure fft_3d_c2c
     module procedure fft_3d_r2c
     module procedure fft_3d_c2r
  end interface


contains


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Initialise the FFT module
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine fft_init_noarg

    implicit none

    call fft_init_arg(PHYSICAL_IN_X)  ! default input is X-pencil data

    return
  end subroutine fft_init_noarg
