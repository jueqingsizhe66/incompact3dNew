  subroutine mpiio_write_real_coarse(ipencil,var,filename,icoarse)

    USE param
    USE variables

    implicit none

    integer, intent(IN) :: ipencil !(x-pencil=1; y-pencil=2; z-pencil=3)
    integer, intent(IN) :: icoarse !(nstat=1; nvisu=2)
    real(mytype), dimension(:,:,:), intent(IN) :: var

    character(len=*) :: filename

    integer (kind=MPI_OFFSET_KIND) :: filesize, disp
    integer, dimension(3) :: sizes, subsizes, starts
    integer :: i,j,k, ierror, newtype, fh

    if (icoarse==1) then
       sizes(1) = xszS(1)
       sizes(2) = yszS(2)
       sizes(3) = zszS(3)

       if (ipencil == 1) then
          subsizes(1) = xszS(1)
          subsizes(2) = xszS(2)
          subsizes(3) = xszS(3)
          starts(1) = xstS(1)-1  ! 0-based index
          starts(2) = xstS(2)-1
          starts(3) = xstS(3)-1
       else if (ipencil == 2) then
          subsizes(1) = yszS(1)
          subsizes(2) = yszS(2)
          subsizes(3) = yszS(3)
          starts(1) = ystS(1)-1
          starts(2) = ystS(2)-1
          starts(3) = ystS(3)-1
       else if (ipencil == 3) then
          subsizes(1) = zszS(1)
          subsizes(2) = zszS(2)
          subsizes(3) = zszS(3)
          starts(1) = zstS(1)-1
          starts(2) = zstS(2)-1
          starts(3) = zstS(3)-1
       endif
    endif

    if (icoarse==2) then
       sizes(1) = xszV(1)
       sizes(2) = yszV(2)
       sizes(3) = zszV(3)

       if (ipencil == 1) then
          subsizes(1) = xszV(1)
          subsizes(2) = xszV(2)
          subsizes(3) = xszV(3)
          starts(1) = xstV(1)-1  ! 0-based index
          starts(2) = xstV(2)-1
          starts(3) = xstV(3)-1
       else if (ipencil == 2) then
          subsizes(1) = yszV(1)
          subsizes(2) = yszV(2)
          subsizes(3) = yszV(3)
          starts(1) = ystV(1)-1
          starts(2) = ystV(2)-1
          starts(3) = ystV(3)-1
       else if (ipencil == 3) then
          subsizes(1) = zszV(1)
          subsizes(2) = zszV(2)
          subsizes(3) = zszV(3)
          starts(1) = zstV(1)-1
          starts(2) = zstV(2)-1
          starts(3) = zstV(3)-1
       endif
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
    call MPI_FILE_SET_VIEW(fh,disp,real_type, &
         newtype,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_WRITE_ALL(fh, var, &
         subsizes(1)*subsizes(2)*subsizes(3), &
         real_type, MPI_STATUS_IGNORE, ierror)
    call MPI_FILE_CLOSE(fh,ierror)
    call MPI_TYPE_FREE(newtype,ierror)


    return
  end subroutine mpiio_write_real_coarse
