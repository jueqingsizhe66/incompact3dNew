  subroutine write_subdomain(ipencil,var,is,ie,js,je,ks,ke,filename)

    implicit none

    integer, intent(IN) :: ipencil !(x-pencil=1; y-pencil=2; z-pencil=3)
    real(mytype), dimension(:,:,:), intent(IN) :: var
    integer, intent(IN) :: is, ie, js, je, ks, ke
    character(len=*), intent(IN) :: filename

    real(mytype), allocatable, dimension(:,:,:) :: wk, wk2
    integer(kind=MPI_OFFSET_KIND) :: filesize, disp
    integer, dimension(3) :: sizes, subsizes, starts
    integer :: color, key, errorcode, newcomm, ierror
    integer :: newtype, fh, data_type, i, j, k
    integer :: i1, i2, j1, j2, k1, k2

    data_type = real_type

    ! validate the input paramters
    if (is<1 .OR. ie>nx_global .OR. js<1 .OR. je>ny_global .OR. &
         ks<1 .OR. ke>nz_global) then
       errorcode = 10
       call decomp_2d_abort(errorcode, &
            'Invalid subdomain specified in I/O')
    end if

    ! create a communicator for all those MPI ranks containing the subdomain
    color = 1
    key = 1
    if (ipencil==1) then
       if (xstart(1)>ie .OR. xend(1)<is .OR. xstart(2)>je .OR. xend(2)<js &
            .OR. xstart(3)>ke .OR. xend(3)<ks) then
          color = 2
       end if
    else if (ipencil==2) then
       if (ystart(1)>ie .OR. yend(1)<is .OR. ystart(2)>je .OR. yend(2)<js &
            .OR. ystart(3)>ke .OR. yend(3)<ks) then
          color = 2
       end if
    else if (ipencil==3) then
       if (zstart(1)>ie .OR. zend(1)<is .OR. zstart(2)>je .OR. zend(2)<js &
            .OR. zstart(3)>ke .OR. zend(3)<ks) then
          color = 2
       end if
    end if
    call MPI_COMM_SPLIT(MPI_COMM_WORLD,color,key,newcomm,ierror)

    if (color==1) then ! only ranks in this group do IO collectively

       ! generate MPI-IO subarray information

       ! global size of the sub-domain to write
       sizes(1) = ie - is + 1
       sizes(2) = je - js + 1
       sizes(3) = ke - ks + 1

       ! 'subsizes' & 'starts' as required by MPI_TYPE_CREATE_SUBARRAY
       ! note the special code whe subdomain only occupy part of the pencil
       if (ipencil==1) then

          subsizes(1) = xsize(1)
          starts(1) = xstart(1) - is
          if (xend(1)>ie .AND. xstart(1)<is) then
             subsizes(1) = ie - is + 1
             starts(1) = 0
          else if (xstart(1)<is) then
             subsizes(1) = xend(1) - is + 1
             starts(1) = 0
          else if (xend(1)>ie) then
             subsizes(1) = ie - xstart(1) + 1
          end if
          subsizes(2) = xsize(2)
          starts(2) = xstart(2) - js
          if (xend(2)>je .AND. xstart(2)<js) then
             subsizes(2) = je - js + 1
             starts(2) = 0
          else if (xstart(2)<js) then
             subsizes(2) = xend(2) - js + 1
             starts(2) = 0
          else if (xend(2)>je) then
             subsizes(2) = je - xstart(2) + 1
          end if
          subsizes(3) = xsize(3)
          starts(3) = xstart(3) - ks
          if (xend(3)>ke .AND. xstart(3)<ks) then
             subsizes(3) = ke - ks + 1
             starts(3) = 0
          else if (xstart(3)<ks) then
             subsizes(3) = xend(3) - ks + 1
             starts(3) = 0
          else if (xend(3)>ke) then
             subsizes(3) = ke - xstart(3) + 1
          end if

       else if (ipencil==2) then

          ! TODO

       else if (ipencil==3) then

          ! TODO

       end if


       ! copy data from orginal to a temp array
       ! pay attention to blocks only partially cover the sub-domain
       if (ipencil==1) then

          if (xend(1)>ie .AND. xstart(1)<is) then
             i1 = is
             i2 = ie
          else if (xend(1)>ie) then
             i1 = xstart(1)
             i2 = ie
          else if (xstart(1)<is) then
             i1 = is
             i2 = xend(1)
          else
             i1 = xstart(1)
             i2 = xend(1)
          end if

          if (xend(2)>je .AND. xstart(2)<js) then
             j1 = js
             j2 = je
          else if (xend(2)>je) then
             j1 = xstart(2)
             j2 = je
          else if (xstart(2)<js) then
             j1 = js
             j2 = xend(2)
          else
             j1 = xstart(2)
             j2 = xend(2)
          end if

          if (xend(3)>ke .AND. xstart(3)<ks) then
             k1 = ks
             k2 = ke
          else if (xend(3)>ke) then
             k1 = xstart(3)
             k2 = ke
          else if (xstart(3)<ks) then
             k1 = ks
             k2 = xend(3)
          else
             k1 = xstart(3)
             k2 = xend(3)
          end if

          allocate(wk(i1:i2, j1:j2, k1:k2))
          allocate(wk2(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)))
          wk2 = var
          do k=k1,k2
             do j=j1,j2
                do i=i1,i2
                   wk(i,j,k) = wk2(i,j,k)
                end do
             end do
          end do

       else if (ipencil==2) then

          ! TODO

       else if (ipencil==3) then

          ! TODO

       end if

       deallocate(wk2)

       ! MPI-IO
       call MPI_TYPE_CREATE_SUBARRAY(3, sizes, subsizes, starts,  &
            MPI_ORDER_FORTRAN, data_type, newtype, ierror)
       call MPI_TYPE_COMMIT(newtype,ierror)
       call MPI_FILE_OPEN(newcomm, filename, &
            MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
            fh, ierror)
       filesize = 0_MPI_OFFSET_KIND
       call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
       disp = 0_MPI_OFFSET_KIND
       call MPI_FILE_SET_VIEW(fh,disp,data_type, &
            newtype,'native',MPI_INFO_NULL,ierror)
       call MPI_FILE_WRITE_ALL(fh, wk, &
            subsizes(1)*subsizes(2)*subsizes(3), &
            data_type, MPI_STATUS_IGNORE, ierror)
       call MPI_FILE_CLOSE(fh,ierror)
       call MPI_TYPE_FREE(newtype,ierror)

       deallocate(wk)

    end if

    return
  end subroutine write_subdomain
