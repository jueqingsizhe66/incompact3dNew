  subroutine write_plane_3d(ipencil,var,iplane,n,filename)

    integer, intent(IN) :: ipencil !(x-pencil=1; y-pencil=2; z-pencil=3)
    real(mytype), dimension(:,:,:), intent(IN) :: var
    integer, intent(IN) :: iplane !(x-plane=1; y-plane=2; z-plane=3)
    integer, intent(IN) :: n ! which plane to write (global coordinate)
    character(len=*) :: filename

    real(mytype), allocatable, dimension(:,:,:) :: wk, wk2
    real(mytype), allocatable, dimension(:,:,:) :: wk2d

    integer (kind=MPI_OFFSET_KIND) :: filesize, disp
    integer, dimension(3) :: sizes, subsizes, starts
    integer :: i,j,k, ierror, newtype, fh

    ! if needed transpose 3D data so that all mpi rank participate I/O
    if (iplane==1) then
       allocate(wk(xsize(1),xsize(2),xsize(3)))
       if (ipencil==1) then
          wk = var
       else if (ipencil==2) then
          call transpose_y_to_x(var,wk)
       else if (ipencil==3) then
          allocate(wk2(ysize(1),ysize(2),ysize(3)))
          call transpose_z_to_y(var,wk2)
          call transpose_y_to_x(wk2,wk)
          deallocate(wk2)
       end if
       allocate(wk2d(1,xsize(2),xsize(3)))
       do k=1,xsize(3)
          do j=1,xsize(2)
             wk2d(1,j,k)=wk(n,j,k)
          end do
       end do
       sizes(1) = 1
       sizes(2) = ysize(2)
       sizes(3) = zsize(3)
       subsizes(1) = 1
       subsizes(2) = xsize(2)
       subsizes(3) = xsize(3)
       starts(1) = 0
       starts(2) = xstart(2)-1
       starts(3) = xstart(3)-1
    else if (iplane==2) then
       allocate(wk(ysize(1),ysize(2),ysize(3)))
       if (ipencil==1) then
          call transpose_x_to_y(var,wk)
       else if (ipencil==2) then
          wk = var
       else if (ipencil==3) then
          call transpose_z_to_y(var,wk)
       end if
       allocate(wk2d(ysize(1),1,ysize(3)))
       do k=1,ysize(3)
          do i=1,ysize(1)
             wk2d(i,1,k)=wk(i,n,k)
          end do
       end do
       sizes(1) = xsize(1)
       sizes(2) = 1
       sizes(3) = zsize(3)
       subsizes(1) = ysize(1)
       subsizes(2) = 1
       subsizes(3) = ysize(3)
       starts(1) = ystart(1)-1
       starts(2) = 0
       starts(3) = ystart(3)-1
    else if (iplane==3) then
       allocate(wk(zsize(1),zsize(2),zsize(3)))
       if (ipencil==1) then
          allocate(wk2(ysize(1),ysize(2),ysize(3)))
          call transpose_x_to_y(var,wk2)
          call transpose_y_to_z(wk2,wk)
          deallocate(wk2)
       else if (ipencil==2) then
          call transpose_y_to_z(var,wk)
       else if (ipencil==3) then
          wk = var
       end if
       allocate(wk2d(zsize(1),zsize(2),1))
       do j=1,zsize(2)
          do i=1,zsize(1)
             wk2d(i,j,1)=wk(i,j,n)
          end do
       end do
       sizes(1) = xsize(1)
       sizes(2) = ysize(2)
       sizes(3) = 1
       subsizes(1) = zsize(1)
       subsizes(2) = zsize(2)
       subsizes(3) = 1
       starts(1) = zstart(1)-1
       starts(2) = zstart(2)-1
       starts(3) = 0
    end if

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
    call MPI_FILE_WRITE_ALL(fh, wk2d, &
         subsizes(1)*subsizes(2)*subsizes(3), &
         real_type, MPI_STATUS_IGNORE, ierror)
    call MPI_FILE_CLOSE(fh,ierror)
    call MPI_TYPE_FREE(newtype,ierror)

    deallocate(wk,wk2d);

    return
  end subroutine write_plane_3d
