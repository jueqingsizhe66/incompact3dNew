  subroutine init_coarser_mesh_statV(i_skip,j_skip,k_skip,from1)

    implicit none

    integer, intent(IN) :: i_skip,j_skip,k_skip
    logical, intent(IN) :: from1  ! .true.  - save 1,n+1,2n+1...
                                  ! .false. - save n,2n,3n...

    integer, dimension(3) :: skip
    integer :: i

    coarse_mesh_starts_from_1 = from1
    iskipV = i_skip
    jskipV = j_skip
    kskipV = k_skip

    skip(1)=iskipV
    skip(2)=jskipV
    skip(3)=kskipV

    do i=1,3
       if (from1) then
          xstV(i) = (xstart(i)+skip(i)-1)/skip(i)
          if (mod(xstart(i)+skip(i)-1,skip(i))/=0) xstV(i)=xstV(i)+1
          xenV(i) = (xend(i)+skip(i)-1)/skip(i)
       else
          xstV(i) = xstart(i)/skip(i)
          if (mod(xstart(i),skip(i))/=0) xstV(i)=xstV(i)+1
          xenV(i) = xend(i)/skip(i)
       end if
       xszV(i) = xenV(i)-xstV(i)+1
    end do

    do i=1,3
       if (from1) then
          ystV(i) = (ystart(i)+skip(i)-1)/skip(i)
          if (mod(ystart(i)+skip(i)-1,skip(i))/=0) ystV(i)=ystV(i)+1
          yenV(i) = (yend(i)+skip(i)-1)/skip(i)
       else
          ystV(i) = ystart(i)/skip(i)
          if (mod(ystart(i),skip(i))/=0) ystV(i)=ystV(i)+1
          yenV(i) = yend(i)/skip(i)
       end if
       yszV(i) = yenV(i)-ystV(i)+1
    end do

    do i=1,3
       if (from1) then
          zstV(i) = (zstart(i)+skip(i)-1)/skip(i)
          if (mod(zstart(i)+skip(i)-1,skip(i))/=0) zstV(i)=zstV(i)+1
          zenV(i) = (zend(i)+skip(i)-1)/skip(i)
       else
          zstV(i) = zstart(i)/skip(i)
          if (mod(zstart(i),skip(i))/=0) zstV(i)=zstV(i)+1
          zenV(i) = zend(i)/skip(i)
       end if
       zszV(i) = zenV(i)-zstV(i)+1
    end do

    return
  end subroutine init_coarser_mesh_statV
