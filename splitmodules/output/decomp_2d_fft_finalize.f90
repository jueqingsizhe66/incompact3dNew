  subroutine decomp_2d_fft_finalize

    implicit none

    call decomp_info_finalize(ph)
    call decomp_info_finalize(sp)

    deallocate(wk2_c2c, wk2_r2c, wk13)

    call finalize_fft_engine

    initialised = .false.

    return
  end subroutine decomp_2d_fft_finalize
