subroutine VISU_INSTA (ux1,uy1,uz1,phi1,ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1,&
     ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2,di2,&
     ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3,phG,uvisu)
!
!############################################################################

USE param
USE variables
USE decomp_2d
USE decomp_2d_io

implicit none

TYPE(DECOMP_INFO) :: phG
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,phi1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2,di2
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3
real(mytype),dimension(xszV(1),xszV(2),xszV(3)) :: uvisu

integer :: code,icomplet
integer :: ijk,nvect1,nvect2,nvect3,i,j,k
character(len=20) nfichier,nfichier1
character(len=20) :: filename

nvect1=xsize(1)*xsize(2)*xsize(3)
!x-derivatives
call derx (ta1,ux1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
call derx (tb1,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
call derx (tc1,uz1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
!y-derivatives
call transpose_x_to_y(ux1,td2)
call transpose_x_to_y(uy1,te2)
call transpose_x_to_y(uz1,tf2)
call dery (ta2,td2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
call dery (tb2,te2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
call dery (tc2,tf2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
!!z-derivatives
call transpose_y_to_z(td2,td3)
call transpose_y_to_z(te2,te3)
call transpose_y_to_z(tf2,tf3)
call derz (ta3,td3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
call derz (tb3,te3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
call derz (tc3,tf3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)
!!all back to x-pencils
call transpose_z_to_y(ta3,td2)
call transpose_z_to_y(tb3,te2)
call transpose_z_to_y(tc3,tf2)
call transpose_y_to_x(td2,tg1)
call transpose_y_to_x(te2,th1)
call transpose_y_to_x(tf2,ti1)
call transpose_y_to_x(ta2,td1)
call transpose_y_to_x(tb2,te1)
call transpose_y_to_x(tc2,tf1)
!du/dx=ta1 du/dy=td1 and du/dz=tg1
!dv/dx=tb1 dv/dy=te1 and dv/dz=th1
!dw/dx=tc1 dw/dy=tf1 and dw/dz=ti1


!############################################################################
!VORTICITY
di1=0.
do ijk=1,nvect1
   di1(ijk,1,1)=sqrt((tf1(ijk,1,1)-th1(ijk,1,1))**2+&
        (tg1(ijk,1,1)-tc1(ijk,1,1))**2+&
        (tb1(ijk,1,1)-td1(ijk,1,1))**2)
enddo
uvisu=0.
call fine_to_coarseV(1,di1,uvisu)
990 format('vort',I3.3)
write(filename, 990) itime/imodulo
call decomp_2d_write_one(1,uvisu,filename,2)
!call decomp_2d_write_one(nx_global,ny_global,nz_global,&
!     1,di1,filename)
!############################################################################

!############################################################################
!VELOCITY
uvisu=0.
call fine_to_coarseV(1,ux1,uvisu)
993 format('ux',I3.3)
      write(filename, 993) itime/imodulo
call decomp_2d_write_one(1,uvisu,filename,2)
!call decomp_2d_write_one(nx_global,ny_global,nz_global,&
!           1,ux1,filename)
uvisu=0.
call fine_to_coarseV(1,uy1,uvisu)
994 format('uy',I3.3)
      write(filename, 994) itime/imodulo
call decomp_2d_write_one(1,uvisu,filename,2)
!call decomp_2d_write_one(nx_global,ny_global,nz_global,&
!           1,uy1,filename)
uvisu=0.
call fine_to_coarseV(1,uz1,uvisu)
995 format('uz',I3.3)
      write(filename, 995) itime/imodulo
call decomp_2d_write_one(1,uvisu,filename,2)
!call decomp_2d_write_one(nx_global,ny_global,nz_global,&
!           1,uz1,filename)
!############################################################################

!############################################################################
!PASSIVE SCALAR
if (iscalar==1) then
uvisu=0.
call fine_to_coarseV(1,phi1,uvisu)
996 format('phi',I3.3)
   write(filename, 996) itime/imodulo
   call decomp_2d_write_one(1,uvisu,filename,2)
!   call decomp_2d_write_one(nx_global,ny_global,nz_global,&
!        1,phi1,filename)
endif
!############################################################################

!############################################################################
!PRESSURE
!IT IS IN A SEPARATE SUBROUTINE
!############################################################################
end subroutine VISU_INSTA
