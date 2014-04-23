  subroutine decomp_2d_fft_get_size(istart, iend, isize)

    implicit none
    integer, dimension(3), intent(OUT) :: istart, iend, isize

    if (format==PHYSICAL_IN_X) then
       istart = sp%zst
       iend   = sp%zen
       isize  = sp%zsz
    else if (format==PHYSICAL_IN_Z) then
       istart = sp%xst
       iend   = sp%xen
       isize  = sp%xsz
    end if

    return
  end subroutine decomp_2d_fft_get_size
