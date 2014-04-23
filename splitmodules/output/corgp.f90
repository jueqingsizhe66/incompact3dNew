subroutine corgp (ux,gx,uy,uz,px,py,pz)
!
!********************************************************************

USE decomp_2d
USE variables
USE param
USE var
USE MPI

implicit none

integer :: ijk,nxyz
real,dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz,px,py,pz
real,dimension(ysize(1),ysize(2),ysize(3)) :: gx

nxyz=xsize(1)*xsize(2)*xsize(3)

do ijk=1,nxyz
   ux(ijk,1,1)=-px(ijk,1,1)+ux(ijk,1,1)
   uy(ijk,1,1)=-py(ijk,1,1)+uy(ijk,1,1)
   uz(ijk,1,1)=-pz(ijk,1,1)+uz(ijk,1,1)
enddo

if (itype==2) then !channel flow
   call transpose_x_to_y(ux,gx)
   call channel(gx)
   call transpose_y_to_x(gx,ux)
endif

return
end subroutine corgp
