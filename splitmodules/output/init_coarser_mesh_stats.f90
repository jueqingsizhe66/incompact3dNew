  subroutine init_coarser_mesh_statS(i_skip,j_skip,k_skip,from1)

    implicit none

    integer, intent(IN) :: i_skip,j_skip,k_skip
    logical, intent(IN) :: from1  ! .true.  - save 1,n+1,2n+1...
                                  ! .false. - save n,2n,3n...

    integer, dimension(3) :: skip
    integer :: i

    coarse_mesh_starts_from_1 = from1
    iskipS = i_skip
    jskipS = j_skip
    kskipS = k_skip

    skip(1)=iskipS
    skip(2)=jskipS
    skip(3)=kskipS

    do i=1,3
       if (from1) then
          xstS(i) = (xstart(i)+skip(i)-1)/skip(i)
          if (mod(xstart(i)+skip(i)-1,skip(i))/=0) xstS(i)=xstS(i)+1
          xenS(i) = (xend(i)+skip(i)-1)/skip(i)
       else
          xstS(i) = xstart(i)/skip(i)
          if (mod(xstart(i),skip(i))/=0) xstS(i)=xstS(i)+1
          xenS(i) = xend(i)/skip(i)
       end if
       xszS(i) = xenS(i)-xstS(i)+1
    end do

    do i=1,3
       if (from1) then
          ystS(i) = (ystart(i)+skip(i)-1)/skip(i)
          if (mod(ystart(i)+skip(i)-1,skip(i))/=0) ystS(i)=ystS(i)+1
          yenS(i) = (yend(i)+skip(i)-1)/skip(i)
       else
          ystS(i) = ystart(i)/skip(i)
          if (mod(ystart(i),skip(i))/=0) ystS(i)=ystS(i)+1
          yenS(i) = yend(i)/skip(i)
       end if
       yszS(i) = yenS(i)-ystS(i)+1
    end do

    do i=1,3
       if (from1) then
          zstS(i) = (zstart(i)+skip(i)-1)/skip(i)
          if (mod(zstart(i)+skip(i)-1,skip(i))/=0) zstS(i)=zstS(i)+1
          zenS(i) = (zend(i)+skip(i)-1)/skip(i)
       else
          zstS(i) = zstart(i)/skip(i)
          if (mod(zstart(i),skip(i))/=0) zstS(i)=zstS(i)+1
          zenS(i) = zend(i)/skip(i)
       end if
       zszS(i) = zenS(i)-zstS(i)+1
    end do

    return
  end subroutine init_coarser_mesh_statS
