  subroutine read_scalar_logical(fh,disp,n,var)

    implicit none

    integer, intent(IN) :: fh
    integer(KIND=MPI_OFFSET_KIND), intent(INOUT) :: disp
    integer, intent(IN) :: n
    logical, dimension(n), intent(INOUT) :: var

    integer :: m, ierror

    call MPI_FILE_SET_VIEW(fh,disp,MPI_LOGICAL, &
         MPI_LOGICAL,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_READ_ALL(fh, var, n, MPI_LOGICAL, &
         MPI_STATUS_IGNORE, ierror)
    call MPI_TYPE_SIZE(MPI_LOGICAL,m,ierror)
    disp = disp + n*m

    return
  end subroutine read_scalar_logical
