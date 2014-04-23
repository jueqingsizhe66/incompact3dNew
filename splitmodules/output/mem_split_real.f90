  subroutine mem_split_real(ndir,in,n1,n2,n3,out,iproc,dist,decomp)

    implicit none

    ! ndir -- x->y: 1; y->z: 2; z->y: 3; y->x: 4.
    integer, intent(IN) :: ndir
    integer, intent(IN) :: n1,n2,n3,iproc
    real(mytype), dimension(n1,n2,n3), intent(IN) :: in
    real(mytype), dimension(*), intent(OUT) :: out
    integer, dimension(0:iproc-1), intent(IN) :: dist
    TYPE(DECOMP_INFO), intent(IN) :: decomp

    integer :: i,j,k, m,i1,i2,pos

#ifndef SHM
    pos = 1
#endif
    do m=0,iproc-1
       if (m==0) then
          i1 = 1
          i2 = dist(0)
       else
          i1 = i2+1
          i2 = i1+dist(m)-1
       endif

       if (ndir==1) then
#ifdef SHM
          pos = decomp%x1disp_o(m) + 1
#endif
          do k=1,n3
             do j=1,n2
                do i=i1,i2
                   out(pos) = in(i,j,k)
                   pos = pos + 1
                enddo
             enddo
          enddo
       else if (ndir==2) then
#ifdef SHM
          pos = decomp%y2disp_o(m) + 1
#endif
          do k=1,n3
             do j=i1,i2
                do i=1,n1
                   out(pos) = in(i,j,k)
                   pos = pos + 1
                enddo
             enddo
          enddo
       else if (ndir==3) then
#ifdef SHM
          pos = decomp%z2disp_o(m) + 1
#endif
          do k=i1,i2
             do j=1,n2
                do i=1,n1
                   out(pos) = in(i,j,k)
                   pos = pos + 1
                enddo
             enddo
          enddo
       else if (ndir==4) then
#ifdef SHM
          pos = decomp%y1disp_o(m) + 1
#endif
          do k=1,n3
             do j=i1,i2
                do i=1,n1
                   out(pos) = in(i,j,k)
                   pos = pos + 1
                enddo
             enddo
          enddo

       endif
    enddo

    return
  end subroutine mem_split_real
