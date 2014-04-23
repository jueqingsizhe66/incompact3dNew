  subroutine c2c_1m_y(inout, isign, decomp)

    implicit none

    complex(mytype), dimension(:,:,:), intent(INOUT) :: inout
    integer, intent(IN) :: isign
    TYPE(DECOMP_INFO), intent(IN) :: decomp

    integer :: i,j,k

    do k=1,decomp%ysz(3)
       do i=1,decomp%ysz(1)
          do j=1,decomp%ysz(2)
             buf(j) = inout(i,j,k)
          end do
          call spcfft(buf,decomp%ysz(2),isign,scratch)
          do j=1,decomp%ysz(2)
             inout(i,j,k) = buf(j)
          end do
       end do
    end do

    return

  end subroutine c2c_1m_y
