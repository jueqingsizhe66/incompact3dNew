  subroutine fft_3d_c2c(in, out, isign)

    implicit none

    complex(mytype), dimension(:,:,:), intent(INOUT) :: in
    complex(mytype), dimension(:,:,:), intent(OUT) :: out
    integer, intent(IN) :: isign

#ifndef OVERWRITE
    complex(mytype), allocatable, dimension(:,:,:) :: wk1
#endif

    if (format==PHYSICAL_IN_X .AND. isign==DECOMP_2D_FFT_FORWARD .OR.  &
         format==PHYSICAL_IN_Z .AND. isign==DECOMP_2D_FFT_BACKWARD) then

       ! ===== 1D FFTs in X =====
#ifdef OVERWRITE
       call c2c_1m_x(in,isign,plan(isign,1))
#else
       allocate (wk1(ph%xsz(1),ph%xsz(2),ph%xsz(3)))
       wk1 = in
       call c2c_1m_x(wk1,isign,plan(isign,1))
#endif

       ! ===== Swap X --> Y; 1D FFTs in Y =====

       if (dims(1)>1) then
#ifdef OVERWRITE
          call transpose_x_to_y(in,wk2_c2c,ph)
#else
          call transpose_x_to_y(wk1,wk2_c2c,ph)
#endif
          call c2c_1m_y(wk2_c2c,isign,plan(isign,2))
       else
#ifdef OVERWRITE
          call c2c_1m_y(in,isign,plan(isign,2))
#else
          call c2c_1m_y(wk1,isign,plan(isign,2))
#endif
       end if

       ! ===== Swap Y --> Z; 1D FFTs in Z =====
       if (dims(1)>1) then
          call transpose_y_to_z(wk2_c2c,out,ph)
       else
#ifdef OVERWRITE
          call transpose_y_to_z(in,out,ph)
#else
          call transpose_y_to_z(wk1,out,ph)
#endif
       end if
       call c2c_1m_z(out,isign,plan(isign,3))

    else if (format==PHYSICAL_IN_X .AND. isign==DECOMP_2D_FFT_BACKWARD &
         .OR. &
         format==PHYSICAL_IN_Z .AND. isign==DECOMP_2D_FFT_FORWARD) then

       ! ===== 1D FFTs in Z =====
#ifdef OVERWRITE
       call c2c_1m_z(in,isign,plan(isign,3))
#else
       allocate (wk1(ph%zsz(1),ph%zsz(2),ph%zsz(3)))
       wk1 = in
       call c2c_1m_z(wk1,isign,plan(isign,3))
#endif

       ! ===== Swap Z --> Y; 1D FFTs in Y =====
       if (dims(1)>1) then
#ifdef OVERWRITE
          call transpose_z_to_y(in,wk2_c2c,ph)
#else
          call transpose_z_to_y(wk1,wk2_c2c,ph)
#endif
          call c2c_1m_y(wk2_c2c,isign,plan(isign,2))
       else  ! out==wk2_c2c if 1D decomposition
#ifdef OVERWRITE
          call transpose_z_to_y(in,out,ph)
#else
          call transpose_z_to_y(wk1,out,ph)
#endif
          call c2c_1m_y(out,isign,plan(isign,2))
       end if

       ! ===== Swap Y --> X; 1D FFTs in X =====
       if (dims(1)>1) then
          call transpose_y_to_x(wk2_c2c,out,ph)
       end if
       call c2c_1m_x(out,isign,plan(isign,1))

    end if

#ifndef OVERWRITE
    deallocate (wk1)
#endif

    return
  end subroutine fft_3d_c2c
