  subroutine init_fft_engine

    implicit none

    if (nrank==0) then
       write(*,*) ' '
       write(*,*) '***** Using the FFTW (version 3.x) engine *****'
       write(*,*) ' '
    end if

    if (format == PHYSICAL_IN_X) then

       ! For C2C transforms
       call c2c_1m_x_plan(plan(-1,1), ph, FFTW_FORWARD )
       call c2c_1m_y_plan(plan(-1,2), ph, FFTW_FORWARD )
       call c2c_1m_z_plan(plan(-1,3), ph, FFTW_FORWARD )
       call c2c_1m_z_plan(plan( 1,3), ph, FFTW_BACKWARD)
       call c2c_1m_y_plan(plan( 1,2), ph, FFTW_BACKWARD)
       call c2c_1m_x_plan(plan( 1,1), ph, FFTW_BACKWARD)

       ! For R2C/C2R tranforms
       call r2c_1m_x_plan(plan(0,1), ph, sp)
       call c2c_1m_y_plan(plan(0,2), sp, FFTW_FORWARD )
       call c2c_1m_z_plan(plan(0,3), sp, FFTW_FORWARD )
       call c2c_1m_z_plan(plan(2,3), sp, FFTW_BACKWARD)
       call c2c_1m_y_plan(plan(2,2), sp, FFTW_BACKWARD)
       call c2r_1m_x_plan(plan(2,1), sp, ph)

    else if (format == PHYSICAL_IN_Z) then

       ! For C2C transforms
       call c2c_1m_z_plan(plan(-1,3), ph, FFTW_FORWARD )
       call c2c_1m_y_plan(plan(-1,2), ph, FFTW_FORWARD )
       call c2c_1m_x_plan(plan(-1,1), ph, FFTW_FORWARD )
       call c2c_1m_x_plan(plan( 1,1), ph, FFTW_BACKWARD)
       call c2c_1m_y_plan(plan( 1,2), ph, FFTW_BACKWARD)
       call c2c_1m_z_plan(plan( 1,3), ph, FFTW_BACKWARD)

       ! For R2C/C2R tranforms
       call r2c_1m_z_plan(plan(0,3), ph, sp)
       call c2c_1m_y_plan(plan(0,2), sp, FFTW_FORWARD )
       call c2c_1m_x_plan(plan(0,1), sp, FFTW_FORWARD )
       call c2c_1m_x_plan(plan(2,1), sp, FFTW_BACKWARD)
       call c2c_1m_y_plan(plan(2,2), sp, FFTW_BACKWARD)
       call c2r_1m_z_plan(plan(2,3), sp, ph)

    end if

    return
  end subroutine init_fft_engine
