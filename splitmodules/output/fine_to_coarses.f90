  subroutine fine_to_coarseS(ipencil,var_fine,var_coarse)

    implicit none

    real(mytype), dimension(:,:,:) :: var_fine
    real(mytype), dimension(:,:,:) :: var_coarse
    integer, intent(IN) :: ipencil

    real(mytype), allocatable, dimension(:,:,:) :: wk, wk2
    integer :: i,j,k

    if (ipencil==1) then
       allocate(wk(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3)))
       allocate(wk2(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)))
       wk2=var_fine
       if (coarse_mesh_starts_from_1) then
          do k=xstS(3),xenS(3)
             do j=xstS(2),xenS(2)
                do i=xstS(1),xenS(1)
                   wk(i,j,k) = wk2((i-1)*iskipS+1,(j-1)*jskipS+1,(k-1)*kskipS+1)
                end do
             end do
          end do
       else
          do k=xstS(3),xenS(3)
             do j=xstS(2),xenS(2)
                do i=xstS(1),xenS(1)
                   wk(i,j,k) = wk2(i*iskipS,j*jskipS,k*kskipS)
                end do
             end do
          end do
       end if
       var_coarse=wk
    else if (ipencil==2) then
       allocate(wk(ystS(1):yenS(1),ystS(2):yenS(2),ystS(3):yenS(3)))
       allocate(wk2(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)))
       wk2=var_fine
       if (coarse_mesh_starts_from_1) then
          do k=ystS(3),yenS(3)
             do j=ystS(2),yenS(2)
                do i=ystS(1),yenS(1)
                   wk(i,j,k) = wk2((i-1)*iskipS+1,(j-1)*jskipS+1,(k-1)*kskipS+1)
                end do
             end do
          end do
       else
          do k=ystS(3),yenS(3)
             do j=ystS(2),yenS(2)
                do i=ystS(1),yenS(1)
                   wk(i,j,k) = wk2(i*iskipS,j*jskipS,k*kskipS)
                end do
             end do
          end do
       end if
       var_coarse=wk
    else if (ipencil==3) then
       allocate(wk(zstS(1):zenS(1),zstS(2):zenS(2),zstS(3):zenS(3)))
       allocate(wk2(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3)))
       wk2=var_fine
       if (coarse_mesh_starts_from_1) then
          do k=zstS(3),zenS(3)
             do j=zstS(2),zenS(2)
                do i=zstS(1),zenS(1)
                   wk(i,j,k) = wk2((i-1)*iskipS+1,(j-1)*jskipS+1,(k-1)*kskipS+1)
                end do
             end do
          end do
       else
          do k=zstS(3),zenS(3)
             do j=zstS(2),zenS(2)
                do i=zstS(1),zenS(1)
                   wk(i,j,k) = wk2(i*iskipS,j*jskipS,k*kskipS)
                end do
             end do
          end do
       end if
       var_coarse=wk
    end if

    deallocate(wk,wk2)

    return
  end subroutine fine_to_coarseS
