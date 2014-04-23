  subroutine decomp_info_init_shm(decomp)

    implicit none

    TYPE(DECOMP_INFO), intent(INOUT) :: decomp

    ! a copy of old displacement array (will be overwritten by shm code)
    allocate(decomp%x1disp_o(0:dims(1)-1),decomp%y1disp_o(0:dims(1)-1), &
         decomp%y2disp_o(0:dims(2)-1),decomp%z2disp_o(0:dims(2)-1))
    decomp%x1disp_o = decomp%x1disp
    decomp%y1disp_o = decomp%y1disp
    decomp%y2disp_o = decomp%y2disp
    decomp%z2disp_o = decomp%z2disp

    call prepare_shared_buffer(decomp%ROW_INFO,MPI_COMM_ROW,decomp)
    call prepare_shared_buffer(decomp%COL_INFO,MPI_COMM_COL,decomp)

    return
  end subroutine decomp_info_init_shm
