  subroutine c2c_1m_z(inout, isign, decomp)

    implicit none

    complex(mytype), dimension(:,:,:), intent(INOUT) :: inout
    integer, intent(IN) :: isign
    TYPE(DECOMP_INFO), intent(IN) :: decomp

    integer :: i,j,k

    do j=1,decomp%zsz(2)
       do i=1,decomp%zsz(1)
          do k=1,decomp%zsz(3)
             buf(k) = inout(i,j,k)
          end do
          call spcfft(buf,decomp%zsz(3),isign,scratch)
          do k=1,decomp%zsz(3)
             inout(i,j,k) = buf(k)
          end do
       end do
    end do

    return

  end subroutine c2c_1m_z
