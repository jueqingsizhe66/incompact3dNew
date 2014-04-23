  subroutine read_scalar_complex(fh,disp,n,var)

    implicit none

    integer, intent(IN) :: fh
    integer(KIND=MPI_OFFSET_KIND), intent(INOUT) :: disp
    integer, intent(IN) :: n
    complex(mytype), dimension(n), intent(INOUT) :: var

    integer :: ierror

    call MPI_FILE_SET_VIEW(fh,disp,complex_type, &
         complex_type,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_READ_ALL(fh, var, n, complex_type, &
         MPI_STATUS_IGNORE, ierror)
    disp = disp + n*mytype_bytes*2

    return
  end subroutine read_scalar_complex
