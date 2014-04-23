  subroutine prepare_buffer(decomp)

  implicit none

  TYPE(DECOMP_INFO), intent(INOUT) :: decomp
  integer :: i

  do i=0, dims(1)-1
    decomp%x1cnts(i) = decomp%x1dist(i)*decomp%xsz(2)*decomp%xsz(3)
    decomp%y1cnts(i) = decomp%ysz(1)*decomp%y1dist(i)*decomp%ysz(3)
    if (i==0) then
      decomp%x1disp(i) = 0  ! displacement is 0-based index
      decomp%y1disp(i) = 0
    else
      decomp%x1disp(i) = decomp%x1disp(i-1) + decomp%x1cnts(i-1)
      decomp%y1disp(i) = decomp%y1disp(i-1) + decomp%y1cnts(i-1)
    endif
  enddo

  do i=0, dims(2)-1
    decomp%y2cnts(i) = decomp%ysz(1)*decomp%y2dist(i)*decomp%ysz(3)
    decomp%z2cnts(i) = decomp%zsz(1)*decomp%zsz(2)*decomp%z2dist(i)
    if (i==0) then
      decomp%y2disp(i) = 0  ! displacement is 0-based index
      decomp%z2disp(i) = 0
    else
      decomp%y2disp(i) = decomp%y2disp(i-1) + decomp%y2cnts(i-1)
      decomp%z2disp(i) = decomp%z2disp(i-1) + decomp%z2cnts(i-1)
    endif
  enddo

  ! simpler information for ALLTOALL
  decomp%x1count = decomp%xsz(1)*decomp%xsz(2)*decomp%xsz(3)/dims(1)
  decomp%y1count = decomp%ysz(1)*decomp%ysz(2)*decomp%ysz(3)/dims(1)
  decomp%y2count = decomp%ysz(1)*decomp%ysz(2)*decomp%ysz(3)/dims(2)
  decomp%z2count = decomp%zsz(1)*decomp%zsz(2)*decomp%zsz(3)/dims(2)

  return
  end subroutine prepare_buffer
