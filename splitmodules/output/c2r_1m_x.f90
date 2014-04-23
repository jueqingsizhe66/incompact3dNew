  subroutine c2r_1m_x(input, output)

    implicit none

    complex(mytype), dimension(:,:,:), intent(IN)  ::  input
    real(mytype), dimension(:,:,:), intent(OUT) :: output

    integer :: i,j,k, d1,d2,d3

    d1 = size(output,1)
    d2 = size(output,2)
    d3 = size(output,3)

    do k=1,d3
       do j=1,d2
          ! Glassman's FFT is c2c only,
          ! needing some pre- and post-processing for c2r
          do i=1,d1/2+1
             buf(i) = input(i,j,k)
          end do
          ! expanding to a full-size complex array
          ! For odd N, the storage is:
          !  1, 2, ...... N/2+1   integer division rounded down
          !     N, ...... N/2+2   => a(i) is conjugate of a(N+2-i)
          ! For even N, the storage is:
          !  1, 2, ...... N/2  , N/2+1
          !     N, ...... N/2+2  again a(i) conjugate of a(N+2-i)
          do i=d1/2+2,d1
             buf(i) =  conjg(buf(d1+2-i))
          end do
          call spcfft(buf,d1,1,scratch)
          do i=1,d1
             ! simply drop imaginary part
             output(i,j,k) = real(buf(i), kind=mytype)
          end do
       end do
    end do

    return

  end subroutine c2r_1m_x
