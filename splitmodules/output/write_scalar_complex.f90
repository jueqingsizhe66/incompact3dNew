  subroutine write_scalar_complex(fh,disp,n,var)

    implicit none

    integer, intent(IN) :: fh
    integer(KIND=MPI_OFFSET_KIND), intent(INOUT) :: disp
    integer, intent(IN) :: n
    complex(mytype), dimension(n), intent(IN) :: var

    integer :: m, ierror

    call MPI_FILE_SET_VIEW(fh,disp,complex_type, &
         complex_type,'native',MPI_INFO_NULL,ierror)
    if (nrank==0) then
       m = n
    else
       m = 0
    end if
    call MPI_FILE_WRITE_ALL(fh, var, m, complex_type, &
         MPI_STATUS_IGNORE, ierror)
    disp = disp + n*mytype_bytes*2

    return
  end subroutine write_scalar_complex
