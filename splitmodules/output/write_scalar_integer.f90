  subroutine write_scalar_integer(fh,disp,n,var)

    implicit none

    integer, intent(IN) :: fh
    integer(KIND=MPI_OFFSET_KIND), intent(INOUT) :: disp
    integer, intent(IN) :: n
    integer, dimension(n), intent(IN) :: var

    integer :: m, ierror

    call MPI_FILE_SET_VIEW(fh,disp,MPI_INTEGER, &
         MPI_INTEGER,'native',MPI_INFO_NULL,ierror)
    if (nrank==0) then
       m = n
    else
       m = 0
    end if
    call MPI_FILE_WRITE_ALL(fh, var, m, MPI_INTEGER, &
         MPI_STATUS_IGNORE, ierror)
    call MPI_TYPE_SIZE(MPI_INTEGER,m,ierror)
    disp = disp + n*m

    return
  end subroutine write_scalar_integer
