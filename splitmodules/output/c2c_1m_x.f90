  subroutine c2c_1m_x(inout, isign, decomp)

    implicit none

    complex(mytype), dimension(:,:,:), intent(INOUT) :: inout
    integer, intent(IN) :: isign
    TYPE(DECOMP_INFO), intent(IN) :: decomp

    integer :: i,j,k

    do k=1,decomp%xsz(3)
       do j=1,decomp%xsz(2)
          do i=1,decomp%xsz(1)
             buf(i) = inout(i,j,k)
          end do
          call spcfft(buf,decomp%xsz(1),isign,scratch)
          do i=1,decomp%xsz(1)
             inout(i,j,k) = buf(i)
          end do
       end do
    end do

    return

  end subroutine c2c_1m_x
