  subroutine glassman_3d_r2c(in_r,nx,ny,nz,out_c)

    implicit none

    integer, intent(IN) :: nx,ny,nz
    real(mytype), dimension(nx,ny,nz) :: in_r
    complex(mytype), dimension(nx/2+1,ny,nz) :: out_c

    complex(mytype), allocatable, dimension(:) :: buf, scratch
    integer :: maxsize, i,j,k

    maxsize = max(nx, max(ny,nz))
    allocate(buf(maxsize))
    allocate(scratch(maxsize))

    ! ===== 1D FFTs in X =====
    do k=1,nz
       do j=1,ny
          ! Glassman's 1D FFT is c2c only,
          ! needing some pre- and post-processing for r2c
          ! pack real input in complex storage
          do i=1,nx
             buf(i) = cmplx(in_r(i,j,k),0._mytype, kind=mytype)
          end do
          call spcfft(buf,nx,-1,scratch)
          ! simply drop the redundant part of the complex output
          do i=1,nx/2+1
             out_c(i,j,k) = buf(i)
          end do
       end do
    end do

    ! ===== 1D FFTs in Y =====
    do k=1,nz
       do i=1,nx/2+1
          do j=1,ny
             buf(j) = out_c(i,j,k)
          end do
          call spcfft(buf,ny,-1,scratch)
          do j=1,ny
             out_c(i,j,k) = buf(j)
          end do
       end do
    end do

    ! ===== 1D FFTs in Z =====
    do j=1,ny
       do i=1,nx/2+1
          do k=1,nz
             buf(k) = out_c(i,j,k)
          end do
          call spcfft(buf,nz,-1,scratch)
          do k=1,nz
             out_c(i,j,k) = buf(k)
          end do
       end do
    end do

    deallocate(buf,scratch)

    return
  end subroutine glassman_3d_r2c
