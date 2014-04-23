module decomp_2d_io

  use decomp_2d
  use MPI
#ifdef T3PIO
  use t3pio
#endif

  implicit none

  private        ! Make everything private unless declared public

  public :: decomp_2d_write_one, decomp_2d_read_one, &
       decomp_2d_write_var, decomp_2d_read_var, &
       decomp_2d_write_scalar, decomp_2d_read_scalar, &
       decomp_2d_write_plane, decomp_2d_write_every

  ! Generic interface to handle multiple data types and decompositions

  interface decomp_2d_write_one
     module procedure write_one_real
     module procedure write_one_complex
     module procedure mpiio_write_real_coarse
  end interface decomp_2d_write_one

  interface decomp_2d_read_one
     module procedure read_one_real
     module procedure read_one_complex
  end interface decomp_2d_read_one

  interface decomp_2d_write_var
     module procedure write_var_real
     module procedure write_var_complex
  end interface decomp_2d_write_var

  interface decomp_2d_read_var
     module procedure read_var_real
     module procedure read_var_complex
  end interface decomp_2d_read_var

  interface decomp_2d_write_scalar
     module procedure write_scalar_real
     module procedure write_scalar_complex
     module procedure write_scalar_integer
  end interface decomp_2d_write_scalar

  interface decomp_2d_read_scalar
     module procedure read_scalar_real
     module procedure read_scalar_complex
     module procedure read_scalar_integer
  end interface decomp_2d_read_scalar

  interface decomp_2d_write_plane
     module procedure write_plane_3d_real
     module procedure write_plane_3d_complex
!     module procedure write_plane_2d
  end interface decomp_2d_write_plane

  interface decomp_2d_write_every
     module procedure write_every_real
     module procedure write_every_complex
  end interface decomp_2d_write_every

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Using MPI-IO library to write a single 3D array to a file
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine write_one_real(ipencil,var,filename,opt_decomp)

    implicit none

    integer, intent(IN) :: ipencil
    real(mytype), dimension(:,:,:), intent(IN) :: var
    character(len=*), intent(IN) :: filename
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    TYPE(DECOMP_INFO) :: decomp
    integer(kind=MPI_OFFSET_KIND) :: filesize, disp
    integer, dimension(3) :: sizes, subsizes, starts
    integer :: ierror, newtype, fh, data_type, info, gs

    data_type = real_type

#include "io_write_one.f90"

    return
  end subroutine write_one_real
