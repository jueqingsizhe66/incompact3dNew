  subroutine write_scalar_real(fh,disp,n,var)

    implicit none

    integer, intent(IN) :: fh             ! file handle
    integer(KIND=MPI_OFFSET_KIND), &
         intent(INOUT) :: disp            ! displacement
    integer, intent(IN) :: n              ! number of scalars
    real(mytype), dimension(n), &
         intent(IN) :: var                ! array of scalars

    integer :: m, ierror

    call MPI_FILE_SET_VIEW(fh,disp,real_type, &
         real_type,'native',MPI_INFO_NULL,ierror)
    if (nrank==0) then
       m = n ! only one rank needs to write
    else
       m = 0
    end if
    call MPI_FILE_WRITE_ALL(fh, var, m, real_type, &
         MPI_STATUS_IGNORE, ierror)
    disp = disp + n*mytype_bytes

    return
  end subroutine write_scalar_real
