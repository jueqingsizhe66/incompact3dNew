  subroutine read_scalar_real(fh,disp,n,var)

    implicit none

    integer, intent(IN) :: fh             ! file handle
    integer(KIND=MPI_OFFSET_KIND), &
         intent(INOUT) :: disp            ! displacement
    integer, intent(IN) :: n              ! number of scalars
    real(mytype), dimension(n), &
         intent(INOUT) :: var             ! array of scalars

    integer :: ierror

    call MPI_FILE_SET_VIEW(fh,disp,real_type, &
         real_type,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_READ_ALL(fh, var, n, real_type, &
         MPI_STATUS_IGNORE, ierror)
    disp = disp + n*mytype_bytes

    return
  end subroutine read_scalar_real
