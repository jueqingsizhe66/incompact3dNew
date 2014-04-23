  subroutine abxyz(ax,ay,az,bx,by,bz,nx,ny,nz,bcx,bcy,bcz)
	
    use param

    implicit none

    integer, intent(IN) :: nx,ny,nz
    integer, intent(IN) :: bcx,bcy,bcz
    real(mytype), dimension(:), intent(OUT) :: ax,bx
    real(mytype), dimension(:), intent(OUT) :: ay,by
    real(mytype), dimension(:), intent(OUT) :: az,bz

    integer :: i,j,k

    if (bcx==0) then
       do i=1,nx
          ax(i) = sin(real(i-1, kind=mytype)*PI/real(nx, kind=mytype))
          bx(i) = cos(real(i-1, kind=mytype)*PI/real(nx, kind=mytype))
       end do
    else if (bcx==1) then
       do i=1,nx
          ax(i) = sin(real(i-1, kind=mytype)*PI/2.0_mytype/ &
               real(nx, kind=mytype))
          bx(i) = cos(real(i-1, kind=mytype)*PI/2.0_mytype/ &
               real(nx, kind=mytype))
       end do
    end if

    if (bcy==0) then
       do j=1,ny
          ay(j) = sin(real(j-1, kind=mytype)*PI/real(ny, kind=mytype))
          by(j) = cos(real(j-1, kind=mytype)*PI/real(ny, kind=mytype))
       end do
    else if (bcy==1) then
       do j=1,ny
          ay(j) = sin(real(j-1, kind=mytype)*PI/2.0_mytype/ &
               real(ny, kind=mytype))
          by(j) = cos(real(j-1, kind=mytype)*PI/2.0_mytype/ &
               real(ny, kind=mytype))
       end do
    end if

    if (bcz==0) then
       do k=1,nz
          az(k) = sin(real(k-1, kind=mytype)*PI/real(nz, kind=mytype))
          bz(k) = cos(real(k-1, kind=mytype)*PI/real(nz, kind=mytype))
       end do
    else if (bcz==1) then
       do k=1,nz
          az(k) = sin(real(k-1, kind=mytype)*PI/2.0_mytype/ &
               real(nz, kind=mytype))
          bz(k) = cos(real(k-1, kind=mytype)*PI/2.0_mytype/ &
               real(nz, kind=mytype))
       end do
    end if

    return
  end subroutine abxyz
