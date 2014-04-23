  subroutine get_dist(nx,ny,nz,decomp)

    integer, intent(IN) :: nx, ny, nz
    TYPE(DECOMP_INFO), intent(INOUT) :: decomp
    integer, allocatable, dimension(:) :: st,en,sz

    allocate(st(0:dims(1)-1))
    allocate(en(0:dims(1)-1))
    call distribute(nx,dims(1),st,en,decomp%x1dist)
    call distribute(ny,dims(1),st,en,decomp%y1dist)
    deallocate(st,en)

    allocate(st(0:dims(2)-1))
    allocate(en(0:dims(2)-1))
    call distribute(ny,dims(2),st,en,decomp%y2dist)
    call distribute(nz,dims(2),st,en,decomp%z2dist)
    deallocate(st,en)

    return
  end subroutine get_dist
