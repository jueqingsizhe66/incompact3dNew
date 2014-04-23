subroutine VISU_PRE (pp3,ta1,tb1,di1,ta2,tb2,di2,&
     ta3,di3,nxmsize,nymsize,nzmsize,phG,ph2,ph3,uvisu)
!
!############################################################################

USE param
USE variables
USE decomp_2d
USE decomp_2d_io

implicit none

integer :: nxmsize,nymsize,nzmsize
TYPE(DECOMP_INFO) :: phG,ph2,ph3
real(mytype),dimension(xszV(1),xszV(2),xszV(3)) :: uvisu

real(mytype),dimension(ph3%zst(1):ph3%zen(1),ph3%zst(2):ph3%zen(2),nzmsize) :: pp3
!Z PENCILS NXM NYM NZM-->NXM NYM NZ
real(mytype),dimension(ph3%zst(1):ph3%zen(1),ph3%zst(2):ph3%zen(2),zsize(3)) :: ta3,di3
!Y PENCILS NXM NYM NZ -->NXM NY NZ
real(mytype),dimension(ph3%yst(1):ph3%yen(1),nymsize,ysize(3)) :: ta2
real(mytype),dimension(ph3%yst(1):ph3%yen(1),ysize(2),ysize(3)) :: tb2,di2
!X PENCILS NXM NY NZ  -->NX NY NZ
real(mytype),dimension(nxmsize,xsize(2),xsize(3)) :: ta1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: tb1,di1

integer :: code,icomplet
integer :: ijk,nvect1,nvect2,nvect3,i,j,k
character(len=20) nfichier,nfichier1
character(len=20) :: filename

!WORK Z-PENCILS
call interiz6(ta3,pp3,di3,sz,cifip6z,cisip6z,ciwip6z,cifz6,cisz6,ciwz6,&
     (ph3%zen(1)-ph3%zst(1)+1),(ph3%zen(2)-ph3%zst(2)+1),nzmsize,zsize(3),1)
!WORK Y-PENCILS
call transpose_z_to_y(ta3,ta2,ph3) !nxm nym nz
call interiy6(tb2,ta2,di2,sy,cifip6y,cisip6y,ciwip6y,cify6,cisy6,ciwy6,&
     (ph3%yen(1)-ph3%yst(1)+1),nymsize,ysize(2),ysize(3),1)
!WORK X-PENCILS
call transpose_y_to_x(tb2,ta1,ph2) !nxm ny nz
call interi6(tb1,ta1,di1,sx,cifip6,cisip6,ciwip6,cifx6,cisx6,ciwx6,&
     nxmsize,xsize(1),xsize(2),xsize(3),1)
!The pressure field on the main mesh is in tb1
!PRESSURE
uvisu=0.
call fine_to_coarseV(1,tb1,uvisu)
990 format('pp',I3.3)
      write(filename, 990) itime/imodulo
call decomp_2d_write_one(1,uvisu,filename,2)
!call decomp_2d_write_one(nx_global,ny_global,nz_global,&
!           1,tb1,filename)

end subroutine VISU_PRE
