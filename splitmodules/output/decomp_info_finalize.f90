  subroutine decomp_info_finalize(decomp)

    implicit none

    TYPE(DECOMP_INFO), intent(INOUT) :: decomp

    deallocate(decomp%x1dist,decomp%y1dist,decomp%y2dist,decomp%z2dist)
    deallocate(decomp%x1cnts,decomp%y1cnts,decomp%y2cnts,decomp%z2cnts)
    deallocate(decomp%x1disp,decomp%y1disp,decomp%y2disp,decomp%z2disp)

#ifdef SHM
    deallocate(decomp%x1disp_o,decomp%y1disp_o,decomp%y2disp_o, &
         decomp%z2disp_o)
    deallocate(decomp%x1cnts_s,decomp%y1cnts_s,decomp%y2cnts_s, &
         decomp%z2cnts_s)
    deallocate(decomp%x1disp_s,decomp%y1disp_s,decomp%y2disp_s, &
         decomp%z2disp_s)
#endif

    return
  end subroutine decomp_info_finalize
