  subroutine finalize_fft_engine

    implicit none

    deallocate(buf,scratch)

    return
  end subroutine finalize_fft_engine
