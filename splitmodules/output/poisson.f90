  subroutine poisson(rhs, bcx, bcy, bcz)

    implicit none

    real(mytype), dimension(:,:,:), intent(INOUT) :: rhs
    integer, intent(IN) :: bcx, bcy, bcz  ! boundary conditions
    integer :: i

    if (bcx==0 .and. bcy==0 .and. bcz==0) then
       call poisson_000(rhs)
    else if (bcx==1 .and. bcy==0 .and. bcz==0) then
       call poisson_100(rhs)
    else if (bcx==0 .and. bcy==1 .and. bcz==0) then
       call poisson_010(rhs)
    else if (bcx==1 .and. bcy==1) then   ! 110 & 111
       call poisson_11x(rhs, bcz)
    else
       stop 'boundary condition not supported'
    end if

    return
  end subroutine poisson
