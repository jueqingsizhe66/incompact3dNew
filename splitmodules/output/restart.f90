subroutine restart(ux1,uy1,uz1,ep1,pp3,phi1,gx1,gy1,gz1,px1,py1,pz1,phis1,&
                   hx1,hy1,hz1,phiss1,phG,irestart)
!
!*******************************************************************

USE decomp_2d
USE decomp_2d_io
USE variables
USE param
USE MPI

implicit none

TYPE(DECOMP_INFO) :: phG
integer :: i,j,k,irestart,nzmsize,fh,ierror,code
real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,ep1
real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: gx1,gy1,gz1
real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: hx1,hy1,hz1
real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: px1,py1,pz1
real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: phi1, phis1, phiss1
real(mytype), dimension(phG%zst(1):phG%zen(1),phG%zst(2):phG%zen(2),phG%zst(3):phG%zen(3)) :: pp3
integer (kind=MPI_OFFSET_KIND) :: filesize, disp
real(mytype) :: xdt
integer, dimension(2) :: dims, dummy_coords
logical, dimension(2) :: dummy_periods

if (iscalar==0) then
   if (nscheme.ne.4) then
      if (irestart==1) then !write
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve.dat', &
              MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
              fh, ierror)
         filesize = 0_MPI_OFFSET_KIND
         call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_write_var(fh,disp,1,ux1)
         call decomp_2d_write_var(fh,disp,1,uy1)
         call decomp_2d_write_var(fh,disp,1,uz1)
         call decomp_2d_write_var(fh,disp,1,gx1)
         call decomp_2d_write_var(fh,disp,1,gy1)
         call decomp_2d_write_var(fh,disp,1,gz1)
         call decomp_2d_write_var(fh,disp,1,px1)
         call decomp_2d_write_var(fh,disp,1,py1)
         call decomp_2d_write_var(fh,disp,1,pz1)
         call decomp_2d_write_var(fh,disp,1,pp3,phG)
         call MPI_FILE_CLOSE(fh,ierror)
      else
         if (nrank==0) print *,'RESTART'
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve.dat', &
              MPI_MODE_RDONLY, MPI_INFO_NULL, &
              fh, ierror)
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_read_var(fh,disp,1,ux1)
         call decomp_2d_read_var(fh,disp,1,uy1)
         call decomp_2d_read_var(fh,disp,1,uz1)
         call decomp_2d_read_var(fh,disp,1,gx1)
         call decomp_2d_read_var(fh,disp,1,gy1)
         call decomp_2d_read_var(fh,disp,1,gz1)
         call decomp_2d_read_var(fh,disp,1,px1)
         call decomp_2d_read_var(fh,disp,1,py1)
         call decomp_2d_read_var(fh,disp,1,pz1)
         call decomp_2d_read_var(fh,disp,1,pp3,phG)
         call MPI_FILE_CLOSE(fh,ierror)
      endif
   else !AB3
     if (irestart==1) then !write
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve.dat', &
              MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
              fh, ierror)
         filesize = 0_MPI_OFFSET_KIND
         call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_write_var(fh,disp,1,ux1)
         call decomp_2d_write_var(fh,disp,1,uy1)
         call decomp_2d_write_var(fh,disp,1,uz1)
         call decomp_2d_write_var(fh,disp,1,gx1)
         call decomp_2d_write_var(fh,disp,1,gy1)
         call decomp_2d_write_var(fh,disp,1,gz1)
         call decomp_2d_write_var(fh,disp,1,hx1)
         call decomp_2d_write_var(fh,disp,1,hy1)
         call decomp_2d_write_var(fh,disp,1,hz1)
         call decomp_2d_write_var(fh,disp,1,px1)
         call decomp_2d_write_var(fh,disp,1,py1)
         call decomp_2d_write_var(fh,disp,1,pz1)
         call decomp_2d_write_var(fh,disp,1,pp3,phG)
         call MPI_FILE_CLOSE(fh,ierror)
      else
         if (nrank==0) print *,'RESTART'
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve.dat', &
              MPI_MODE_RDONLY, MPI_INFO_NULL, &
              fh, ierror)
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_read_var(fh,disp,1,ux1)
         call decomp_2d_read_var(fh,disp,1,uy1)
         call decomp_2d_read_var(fh,disp,1,uz1)
         call decomp_2d_read_var(fh,disp,1,gx1)
         call decomp_2d_read_var(fh,disp,1,gy1)
         call decomp_2d_read_var(fh,disp,1,gz1)
         call decomp_2d_read_var(fh,disp,1,hx1)
         call decomp_2d_read_var(fh,disp,1,hy1)
         call decomp_2d_read_var(fh,disp,1,hz1)
         call decomp_2d_read_var(fh,disp,1,px1)
         call decomp_2d_read_var(fh,disp,1,py1)
         call decomp_2d_read_var(fh,disp,1,pz1)
         call decomp_2d_read_var(fh,disp,1,pp3,phG)
         call MPI_FILE_CLOSE(fh,ierror)
      endif
   endif
else !SCALAR
if (nscheme.ne.4) then
      if (irestart==1) then !write
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve.dat', &
              MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
              fh, ierror)
         filesize = 0_MPI_OFFSET_KIND
         call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_write_var(fh,disp,1,ux1)
         call decomp_2d_write_var(fh,disp,1,uy1)
         call decomp_2d_write_var(fh,disp,1,uz1)
         call decomp_2d_write_var(fh,disp,1,gx1)
         call decomp_2d_write_var(fh,disp,1,gy1)
         call decomp_2d_write_var(fh,disp,1,gz1)
         call decomp_2d_write_var(fh,disp,1,px1)
         call decomp_2d_write_var(fh,disp,1,py1)
         call decomp_2d_write_var(fh,disp,1,pz1)
         call decomp_2d_write_var(fh,disp,1,pp3,phG)
         call decomp_2d_write_var(fh,disp,1,phi1)
         call decomp_2d_write_var(fh,disp,1,phis1)
         call MPI_FILE_CLOSE(fh,ierror)
      else
         if (nrank==0) print *,'RESTART'
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve.dat', &
              MPI_MODE_RDONLY, MPI_INFO_NULL, &
              fh, ierror)
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_read_var(fh,disp,1,ux1)
         call decomp_2d_read_var(fh,disp,1,uy1)
         call decomp_2d_read_var(fh,disp,1,uz1)
         call decomp_2d_read_var(fh,disp,1,gx1)
         call decomp_2d_read_var(fh,disp,1,gy1)
         call decomp_2d_read_var(fh,disp,1,gz1)
         call decomp_2d_read_var(fh,disp,1,px1)
         call decomp_2d_read_var(fh,disp,1,py1)
         call decomp_2d_read_var(fh,disp,1,pz1)
         call decomp_2d_read_var(fh,disp,1,pp3,phG)
         call decomp_2d_read_var(fh,disp,1,phi1)
         call decomp_2d_read_var(fh,disp,1,phis1)
         call MPI_FILE_CLOSE(fh,ierror)
      endif
   else !SCALAR + AB3
     if (irestart==1) then !write
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve.dat', &
              MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
              fh, ierror)
         filesize = 0_MPI_OFFSET_KIND
         call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_write_var(fh,disp,1,ux1)
         call decomp_2d_write_var(fh,disp,1,uy1)
         call decomp_2d_write_var(fh,disp,1,uz1)
         call decomp_2d_write_var(fh,disp,1,gx1)
         call decomp_2d_write_var(fh,disp,1,gy1)
         call decomp_2d_write_var(fh,disp,1,gz1)
         call decomp_2d_write_var(fh,disp,1,hx1)
         call decomp_2d_write_var(fh,disp,1,hy1)
         call decomp_2d_write_var(fh,disp,1,hz1)
         call decomp_2d_write_var(fh,disp,1,px1)
         call decomp_2d_write_var(fh,disp,1,py1)
         call decomp_2d_write_var(fh,disp,1,pz1)
         call decomp_2d_write_var(fh,disp,1,pp3,phG)
         call decomp_2d_write_var(fh,disp,1,phi1)
         call decomp_2d_write_var(fh,disp,1,phis1)
         call decomp_2d_write_var(fh,disp,1,phiss1)
         call MPI_FILE_CLOSE(fh,ierror)
      else
         if (nrank==0) print *,'RESTART'
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve.dat', &
              MPI_MODE_RDONLY, MPI_INFO_NULL, &
              fh, ierror)
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_read_var(fh,disp,1,ux1)
         call decomp_2d_read_var(fh,disp,1,uy1)
         call decomp_2d_read_var(fh,disp,1,uz1)
         call decomp_2d_read_var(fh,disp,1,gx1)
         call decomp_2d_read_var(fh,disp,1,gy1)
         call decomp_2d_read_var(fh,disp,1,gz1)
         call decomp_2d_read_var(fh,disp,1,hx1)
         call decomp_2d_read_var(fh,disp,1,hy1)
         call decomp_2d_read_var(fh,disp,1,hz1)
         call decomp_2d_read_var(fh,disp,1,px1)
         call decomp_2d_read_var(fh,disp,1,py1)
         call decomp_2d_read_var(fh,disp,1,pz1)
         call decomp_2d_read_var(fh,disp,1,pp3,phG)
         call decomp_2d_read_var(fh,disp,1,phi1)
         call decomp_2d_read_var(fh,disp,1,phis1)
         call decomp_2d_read_var(fh,disp,1,phiss1)
         call MPI_FILE_CLOSE(fh,ierror)
      endif
   endif
endif

if (irestart==0) then
! reconstruction of the dp/dx, dp/dy and dp/dz from px1,py1 and pz1
! Temporal scheme (1:AB2, 2: RK3, 3:RK4, 4:AB3)
   if (nscheme==1) xdt=adt(1)+bdt(1)
   if (nscheme==2) xdt=(3./4.)*dt + (-5./12.)*dt
   if (nscheme==3) xdt=0.041717869325*dt
   if (nscheme==4) xdt=adt(1)+bdt(1)+cdt(1)

   do k=1,xsize(3)
   do j=1,xsize(2)
      dpdyx1(j,k)=py1(1,j,k)/xdt
      dpdzx1(j,k)=pz1(1,j,k)/xdt
      dpdyxn(j,k)=py1(nx,j,k)/xdt
      dpdzxn(j,k)=pz1(nx,j,k)/xdt
   enddo
   enddo

   if (xsize(3)==1) then
      do j=1,xsize(2)
      do i=1,xsize(1)
         dpdxz1(i,j)=px1(i,j,1)/xdt
         dpdyz1(i,j)=py1(i,j,1)/xdt
      enddo
      enddo
   endif
   if (xsize(3)==nz) then
      do j=1,xsize(2)
      do i=1,xsize(1)
         dpdxzn(i,j)=px1(i,j,nz)/xdt
         dpdyzn(i,j)=py1(i,j,nz)/xdt
      enddo
      enddo
   endif

   ! determine the processor grid in use
   call MPI_CART_GET(DECOMP_2D_COMM_CART_X, 2, &
         dims, dummy_periods, dummy_coords, code)

   if (dims(1)==1) then
      do k=1,xsize(3)
      do i=1,xsize(1)
      dpdxy1(i,k)=px1(i,1,k)/xdt
      dpdzy1(i,k)=pz1(i,1,k)/xdt
      enddo
      enddo
      do k=1,xsize(3)
      do i=1,xsize(1)
      dpdxyn(i,k)=px1(i,xsize(2),k)/xdt
      dpdzyn(i,k)=pz1(i,xsize(2),k)/xdt
      enddo
      enddo
   else
!find j=1 and j=ny
      if (xstart(2)==1) then
         do k=1,xsize(3)
         do i=1,xsize(1)
         dpdxy1(i,k)=px1(i,1,k)/xdt
         dpdzy1(i,k)=pz1(i,1,k)/xdt
         enddo
         enddo
      endif
!      print *,nrank,xstart(2),ny-(nym/p_row)
       if (ny-(nym/dims(1))==xstart(2)) then
         do k=1,xsize(3)
         do i=1,xsize(1)
         dpdxyn(i,k)=px1(i,xsize(2),k)/xdt
         dpdzyn(i,k)=pz1(i,xsize(2),k)/xdt
         enddo
         enddo
      endif

   endif


   if (nrank==0) print *,'reconstruction pressure gradients done!'
endif

if (irestart==0) then
if (ivirt==2) then
   call MPI_FILE_OPEN(MPI_COMM_WORLD, 'epsilon.dat', &
        MPI_MODE_RDONLY, MPI_INFO_NULL, &
        fh, ierror)
   disp = 0_MPI_OFFSET_KIND
   call decomp_2d_read_var(fh,disp,1,ep1)
   call MPI_FILE_CLOSE(fh,ierror)
   if (nrank==0) print *,'read epsilon file done from restart'
endif
endif

end subroutine restart
