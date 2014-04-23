  subroutine fft_3d_r2c(in_r, out_c)

    implicit none

    real(mytype), dimension(:,:,:), intent(IN) :: in_r
    complex(mytype), dimension(:,:,:), intent(OUT) :: out_c

    if (format==PHYSICAL_IN_X) then

       ! ===== 1D FFTs in X =====
       call r2c_1m_x(in_r,wk13)

       ! ===== Swap X --> Y; 1D FFTs in Y =====
       if (dims(1)>1) then
          call transpose_x_to_y(wk13,wk2_r2c,sp)
          call c2c_1m_y(wk2_r2c,-1,plan(0,2))
       else
          call c2c_1m_y(wk13,-1,plan(0,2))
       end if

       ! ===== Swap Y --> Z; 1D FFTs in Z =====
       if (dims(1)>1) then
          call transpose_y_to_z(wk2_r2c,out_c,sp)
       else
          call transpose_y_to_z(wk13,out_c,sp)
       end if
       call c2c_1m_z(out_c,-1,plan(0,3))

    else if (format==PHYSICAL_IN_Z) then

       ! ===== 1D FFTs in Z =====
       call r2c_1m_z(in_r,wk13)

       ! ===== Swap Z --> Y; 1D FFTs in Y =====
       if (dims(1)>1) then
          call transpose_z_to_y(wk13,wk2_r2c,sp)
          call c2c_1m_y(wk2_r2c,-1,plan(0,2))
       else  ! out_c==wk2_r2c if 1D decomposition
          call transpose_z_to_y(wk13,out_c,sp)
          call c2c_1m_y(out_c,-1,plan(0,2))
       end if

       ! ===== Swap Y --> X; 1D FFTs in X =====
       if (dims(1)>1) then
          call transpose_y_to_x(wk2_r2c,out_c,sp)
       end if
       call c2c_1m_x(out_c,-1,plan(0,1))

    end if

    return
  end subroutine fft_3d_r2c
