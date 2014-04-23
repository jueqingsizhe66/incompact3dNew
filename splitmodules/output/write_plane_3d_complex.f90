  subroutine write_plane_3d_complex(ipencil,var,iplane,n, &
       filename,opt_decomp)

    implicit none

    integer, intent(IN) :: ipencil !(x-pencil=1; y-pencil=2; z-pencil=3)
    complex(mytype), dimension(:,:,:), intent(IN) :: var
    integer, intent(IN) :: iplane !(x-plane=1; y-plane=2; z-plane=3)
    integer, intent(IN) :: n ! which plane to write (global coordinate)
    character(len=*), intent(IN) :: filename
    TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

    complex(mytype), allocatable, dimension(:,:,:) :: wk, wk2
    complex(mytype), allocatable, dimension(:,:,:) :: wk2d
    TYPE(DECOMP_INFO) :: decomp
    integer(kind=MPI_OFFSET_KIND) :: filesize, disp
    integer, dimension(3) :: sizes, subsizes, starts
    integer :: i,j,k, ierror, newtype, fh, data_type

    data_type = complex_type

#include "io_write_plane.f90"

    return
  end subroutine write_plane_3d_complex
