  subroutine mpiio_write_real_decomp(decomp,ipencil,var,filename)

    implicit none

    TYPE(DECOMP_INFO), intent(IN) :: decomp
    integer, intent(IN) :: ipencil !(x-pencil=1; y-pencil=2; z-pencil=3)
    real(mytype), dimension(:,:,:), intent(IN) :: var
    character(len=*) :: filename

    integer (kind=MPI_OFFSET_KIND) :: filesize, disp
    integer, dimension(3) :: sizes, subsizes, starts
    integer :: i,j,k, ierror, newtype, fh

    sizes(1) = decomp%xsz(1)
    sizes(2) = decomp%ysz(2)
    sizes(3) = decomp%zsz(3)

    if (ipencil == 1) then
       subsizes(1) = decomp%xsz(1)
       subsizes(2) = decomp%xsz(2)
       subsizes(3) = decomp%xsz(3)
       starts(1) = decomp%xst(1)-1  ! 0-based index
       starts(2) = decomp%xst(2)-1
       starts(3) = decomp%xst(3)-1
    else if (ipencil == 2) then
       subsizes(1) = decomp%ysz(1)
       subsizes(2) = decomp%ysz(2)
       subsizes(3) = decomp%ysz(3)
       starts(1) = decomp%yst(1)-1
       starts(2) = decomp%yst(2)-1
       starts(3) = decomp%yst(3)-1
    else if (ipencil == 3) then
       subsizes(1) = decomp%zsz(1)
       subsizes(2) = decomp%zsz(2)
       subsizes(3) = decomp%zsz(3)
       starts(1) = decomp%zst(1)-1
       starts(2) = decomp%zst(2)-1
       starts(3) = decomp%zst(3)-1
    endif

    call MPI_TYPE_CREATE_SUBARRAY(3, sizes, subsizes, starts,  &
         MPI_ORDER_FORTRAN, real_type, newtype, ierror)
    call MPI_TYPE_COMMIT(newtype,ierror)
    call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, &
         MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
         fh, ierror)
    filesize = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
    disp = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_VIEW(fh,disp,MPI_INTEGER, &
         newtype,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_WRITE_ALL(fh, var, &
         subsizes(1)*subsizes(2)*subsizes(3), &
         real_type, MPI_STATUS_IGNORE, ierror)
    call MPI_FILE_CLOSE(fh,ierror)
    call MPI_TYPE_FREE(newtype,ierror)

    return
  end subroutine mpiio_write_real_decomp