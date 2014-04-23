module decomp_2d_fft

  use decomp_2d  ! 2D decomposition module
  use glassman

  implicit none

  private        ! Make everything private unless declared public

  ! engine-specific global variables
  complex(mytype), allocatable, dimension(:) :: buf, scratch

  ! common code used for all engines, including global variables,
  ! generic interface definitions and several subroutines
#include "fft_common.f90"

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  This routine performs one-time initialisations for the FFT engine
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_fft_engine

    implicit none

    integer :: cbuf_size

    if (nrank==0) then
       write(*,*) ' '
       write(*,*) '***** Using the generic FFT engine *****'
       write(*,*) ' '
    end if

    cbuf_size = max(ph%xsz(1), ph%ysz(2))
    cbuf_size = max(cbuf_size, ph%zsz(3))
    allocate(buf(cbuf_size))
    allocate(scratch(cbuf_size))

    return
  end subroutine init_fft_engine
