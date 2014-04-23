  subroutine read_var(fh,disp,nx,ny,nz,ipencil,var)

    implicit none

    integer, intent(IN) :: fh             ! file handle
    integer, intent(IN) :: nx,ny,nz       ! global size
    integer, intent(IN) :: ipencil        ! pencil orientation
    real(mytype), dimension(:,:,:), &
         intent(INOUT) :: var                ! distributed 3D data array
    integer(KIND=MPI_OFFSET_KIND), &
         intent(INOUT) :: disp            ! displacement

    integer, dimension(3) :: sizes, subsizes, starts
    integer :: i,j,k, ierror, newtype

    ! Create file type and set file view
    sizes(1) = nx
    sizes(2) = ny
    sizes(3) = nz
    if (ipencil == 1) then
       subsizes(1) = xsize(1)
       subsizes(2) = xsize(2)
       subsizes(3) = xsize(3)
       starts(1) = xstart(1)-1  ! 0-based index
       starts(2) = xstart(2)-1
       starts(3) = xstart(3)-1
    else if (ipencil == 2) then
       subsizes(1) = ysize(1)
       subsizes(2) = ysize(2)
       subsizes(3) = ysize(3)
       starts(1) = ystart(1)-1
       starts(2) = ystart(2)-1
       starts(3) = ystart(3)-1
    else if (ipencil == 3) then
       subsizes(1) = zsize(1)
       subsizes(2) = zsize(2)
       subsizes(3) = zsize(3)
       starts(1) = zstart(1)-1
       starts(2) = zstart(2)-1
       starts(3) = zstart(3)-1
    endif

    call MPI_TYPE_CREATE_SUBARRAY(3, sizes, subsizes, starts,  &
         MPI_ORDER_FORTRAN, real_type, newtype, ierror)
    call MPI_TYPE_COMMIT(newtype,ierror)
    call MPI_FILE_SET_VIEW(fh,disp,MPI_INTEGER, &
         newtype,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_READ_ALL(fh, var, &
         subsizes(1)*subsizes(2)*subsizes(3), &
         real_type, MPI_STATUS_IGNORE, ierror)

    call MPI_TYPE_FREE(newtype,ierror)

    ! update displacement for next write
    disp = disp + nx*ny*nz*mytype

    return
  end subroutine read_var
