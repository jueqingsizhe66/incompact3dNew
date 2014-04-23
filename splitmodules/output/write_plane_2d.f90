  subroutine write_plane_2d(ipencil,var,filename)
    integer, intent(IN) :: ipencil !(x-pencil=1; y-pencil=2; z-pencil=3)
    real(mytype), dimension(:,:), intent(IN) :: var ! 2D array
    character(len=*) :: filename

    if (ipencil==1) then
       ! var should be defined as var(xsize(2)

    else if (ipencil==2) then

    else if (ipencil==3) then

    end if

    return
  end subroutine write_plane_2d
