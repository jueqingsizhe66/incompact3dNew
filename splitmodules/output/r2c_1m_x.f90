  subroutine r2c_1m_x(input, output)

    implicit none

    real(mytype), dimension(:,:,:), intent(IN)  ::  input
    complex(mytype), dimension(:,:,:), intent(OUT) :: output

    integer :: i,j,k, s1,s2,s3, d1

    s1 = size(input,1)
    s2 = size(input,2)
    s3 = size(input,3)
    d1 = size(output,1)

    do k=1,s3
       do j=1,s2
          ! Glassman's FFT is c2c only,
          ! needing some pre- and post-processing for r2c
          ! pack real input in complex storage
          do i=1,s1
             buf(i) = cmplx(input(i,j,k),0._mytype, kind=mytype)
          end do
          call spcfft(buf,s1,-1,scratch)
          ! note d1 ~ s1/2+1
          ! simply drop the redundant part of the complex output
          do i=1,d1
             output(i,j,k) = buf(i)
          end do
       end do
    end do

    return

  end subroutine r2c_1m_x
