  subroutine fine_to_coarseV(ipencil,var_fine,var_coarse)

    implicit none

    real(mytype), dimension(:,:,:) :: var_fine
    real(mytype), dimension(:,:,:) :: var_coarse
    integer, intent(IN) :: ipencil

    real(mytype), allocatable, dimension(:,:,:) :: wk, wk2
    integer :: i,j,k

    if (ipencil==1) then
       allocate(wk(xstV(1):xenV(1),xstV(2):xenV(2),xstV(3):xenV(3)))
       allocate(wk2(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)))
       wk2=var_fine
       if (coarse_mesh_starts_from_1) then
          do k=xstV(3),xenV(3)
             do j=xstV(2),xenV(2)
                do i=xstV(1),xenV(1)
                   wk(i,j,k) = wk2((i-1)*iskipV+1,(j-1)*jskipV+1,(k-1)*kskipV+1)
                end do
             end do
          end do
       else
          do k=xstV(3),xenV(3)
             do j=xstV(2),xenV(2)
                do i=xstV(1),xenV(1)
                   wk(i,j,k) = wk2(i*iskipV,j*jskipV,k*kskipV)
                end do
             end do
          end do
       end if
       var_coarse=wk
    else if (ipencil==2) then
       allocate(wk(ystV(1):yenV(1),ystV(2):yenV(2),ystV(3):yenV(3)))
       allocate(wk2(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)))
       wk2=var_fine
       if (coarse_mesh_starts_from_1) then
          do k=ystV(3),yenV(3)
             do j=ystV(2),yenV(2)
                do i=ystV(1),yenV(1)
                   wk(i,j,k) = wk2((i-1)*iskipV+1,(j-1)*jskipV+1,(k-1)*kskipV+1)
                end do
             end do
          end do
       else
          do k=ystV(3),yenV(3)
             do j=ystV(2),yenV(2)
                do i=ystV(1),yenV(1)
                   wk(i,j,k) = wk2(i*iskipV,j*jskipV,k*kskipV)
                end do
             end do
          end do
       end if
       var_coarse=wk
    else if (ipencil==3) then
       allocate(wk(zstV(1):zenV(1),zstV(2):zenV(2),zstV(3):zenV(3)))
       allocate(wk2(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3)))
       wk2=var_fine
       if (coarse_mesh_starts_from_1) then
          do k=zstV(3),zenV(3)
             do j=zstV(2),zenV(2)
                do i=zstV(1),zenV(1)
                   wk(i,j,k) = wk2((i-1)*iskipV+1,(j-1)*jskipV+1,(k-1)*kskipV+1)
                end do
             end do
          end do
       else
          do k=zstV(3),zenV(3)
             do j=zstV(2),zenV(2)
                do i=zstV(1),zenV(1)
                   wk(i,j,k) = wk2(i*iskipV,j*jskipV,k*kskipV)
                end do
             end do
          end do
       end if
       var_coarse=wk
    end if

    deallocate(wk,wk2)

    return
  end subroutine fine_to_coarseV
