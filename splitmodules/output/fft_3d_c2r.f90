  subroutine fft_3d_c2r(in_c, out_r)

    implicit none

    complex(mytype), dimension(:,:,:), intent(INOUT) :: in_c
    real(mytype), dimension(:,:,:), intent(OUT) :: out_r

#ifndef OVERWRITE
    complex(mytype), allocatable, dimension(:,:,:) :: wk1
#endif

    if (format==PHYSICAL_IN_X) then

       ! ===== 1D FFTs in Z =====
#ifdef OVERWRITE
       call c2c_1m_z(in_c,1,plan(2,3))
#else
       allocate (wk1(sp%zsz(1),sp%zsz(2),sp%zsz(3)))
       wk1 = in_c
       call c2c_1m_z(wk1,1,plan(2,3))
#endif

       ! ===== Swap Z --> Y; 1D FFTs in Y =====
#ifdef OVERWRITE
       call transpose_z_to_y(in_c,wk2_r2c,sp)
#else
       call transpose_z_to_y(wk1,wk2_r2c,sp)
#endif
       call c2c_1m_y(wk2_r2c,1,plan(2,2))

       ! ===== Swap Y --> X; 1D FFTs in X =====
       if (dims(1)>1) then
          call transpose_y_to_x(wk2_r2c,wk13,sp)
          call c2r_1m_x(wk13,out_r)
       else
          call c2r_1m_x(wk2_r2c,out_r)
       end if

    else if (format==PHYSICAL_IN_Z) then

       ! ===== 1D FFTs in X =====
#ifdef OVERWRITE
       call c2c_1m_x(in_c,1,plan(2,1))
#else
       allocate(wk1(sp%xsz(1),sp%xsz(2),sp%xsz(3)))
       wk1 = in_c
       call c2c_1m_x(wk1,1,plan(2,1))
#endif

       ! ===== Swap X --> Y; 1D FFTs in Y =====
       if (dims(1)>1) then
#ifdef OVERWRITE
          call transpose_x_to_y(in_c,wk2_r2c,sp)
#else
          call transpose_x_to_y(wk1,wk2_r2c,sp)
#endif
          call c2c_1m_y(wk2_r2c,1,plan(2,2))
       else  ! in_c==wk2_r2c if 1D decomposition
#ifdef OVERWRITE
          call c2c_1m_y(in_c,1,plan(2,2))
#else
          call c2c_1m_y(wk1,1,plan(2,2))
#endif
       end if

       ! ===== Swap Y --> Z; 1D FFTs in Z =====
       if (dims(1)>1) then
          call transpose_y_to_z(wk2_r2c,wk13,sp)
       else
#ifdef OVERWRITE
          call transpose_y_to_z(in_c,wk13,sp)
#else
          call transpose_y_to_z(wk1,wk13,sp)
#endif
       end if
       call c2r_1m_z(wk13,out_r)

    end if

#ifndef OVERWRITE
    deallocate (wk1)
#endif

    return
  end subroutine fft_3d_c2r
