  subroutine init_coarser_mesh_statP(i_skip,j_skip,k_skip,from1)

    implicit none

    integer, intent(IN) :: i_skip,j_skip,k_skip
    logical, intent(IN) :: from1  ! .true.  - save 1,n+1,2n+1...
                                  ! .false. - save n,2n,3n...

    integer, dimension(3) :: skip
    integer :: i

    coarse_mesh_starts_from_1 = from1
    iskipP = i_skip
    jskipP = j_skip
    kskipP = k_skip

    skip(1)=iskipP
    skip(2)=jskipP
    skip(3)=kskipP

    do i=1,3
       if (from1) then
          xstP(i) = (xstart(i)+skip(i)-1)/skip(i)
          if (mod(xstart(i)+skip(i)-1,skip(i))/=0) xstP(i)=xstP(i)+1
          xenP(i) = (xend(i)+skip(i)-1)/skip(i)
       else
          xstP(i) = xstart(i)/skip(i)
          if (mod(xstart(i),skip(i))/=0) xstP(i)=xstP(i)+1
          xenP(i) = xend(i)/skip(i)
       end if
       xszP(i) = xenP(i)-xstP(i)+1
    end do

    do i=1,3
       if (from1) then
          ystP(i) = (ystart(i)+skip(i)-1)/skip(i)
          if (mod(ystart(i)+skip(i)-1,skip(i))/=0) ystP(i)=ystP(i)+1
          yenP(i) = (yend(i)+skip(i)-1)/skip(i)
       else
          ystP(i) = ystart(i)/skip(i)
          if (mod(ystart(i),skip(i))/=0) ystP(i)=ystP(i)+1
          yenP(i) = yend(i)/skip(i)
       end if
       yszP(i) = yenP(i)-ystP(i)+1
    end do

    do i=1,3
       if (from1) then
          zstP(i) = (zstart(i)+skip(i)-1)/skip(i)
          if (mod(zstart(i)+skip(i)-1,skip(i))/=0) zstP(i)=zstP(i)+1
          zenP(i) = (zend(i)+skip(i)-1)/skip(i)
       else
          zstP(i) = zstart(i)/skip(i)
          if (mod(zstart(i),skip(i))/=0) zstP(i)=zstP(i)+1
          zenP(i) = zend(i)/skip(i)
       end if
       zszP(i) = zenP(i)-zstP(i)+1
    end do

    return
  end subroutine init_coarser_mesh_statP
