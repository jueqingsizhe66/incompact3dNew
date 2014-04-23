subroutine STATISTIC(ux1,uy1,uz1,phi1,ta1,umean,vmean,wmean,phimean,uumean,vvmean,wwmean,&
     uvmean,uwmean,vwmean,phiphimean,tmean)
!
!############################################################################

USE param
USE variables
USE decomp_2d
USE decomp_2d_io

implicit none

real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,phi1
real(mytype),dimension(xszS(1),xszS(2),xszS(3)) :: umean,vmean,wmean,uumean,vvmean,wwmean,uvmean,uwm
real(mytype),dimension(xszS(1),xszS(2),xszS(3)) :: phimean, phiphimean
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ta1

!umean=ux1
call fine_to_coarseS(1,ux1,tmean)
umean(:,:,:)=umean(:,:,:)+tmean(:,:,:)

!vmean=uy1
call fine_to_coarseS(1,uy1,tmean)
vmean(:,:,:)=vmean(:,:,:)+tmean(:,:,:)

!wmean=uz1
call fine_to_coarseS(1,uz1,tmean)
wmean(:,:,:)=wmean(:,:,:)+tmean(:,:,:)

if (iscalar==1) then
   !phimean=phi1
   call fine_to_coarseS(1,phi1,tmean)
   phimean(:,:,:)=phimean(:,:,:)+tmean(:,:,:)
endif

!uumean=ux1*ux1
ta1(:,:,:)=ux1(:,:,:)*ux1(:,:,:)
call fine_to_coarseS(1,ta1,tmean)
uumean(:,:,:)=uumean(:,:,:)+tmean(:,:,:)

!vvmean=uy1*uy1
ta1(:,:,:)=uy1(:,:,:)*uy1(:,:,:)
call fine_to_coarseS(1,ta1,tmean)
vvmean(:,:,:)=vvmean(:,:,:)+tmean(:,:,:)

!wwmean=uz1*uz1
ta1(:,:,:)=uz1(:,:,:)*uz1(:,:,:)
call fine_to_coarseS(1,ta1,tmean)
wwmean(:,:,:)=wwmean(:,:,:)+tmean(:,:,:)

!uvmean=ux1*uy1
ta1(:,:,:)=ux1(:,:,:)*uy1(:,:,:)
call fine_to_coarseS(1,ta1,tmean)
uvmean(:,:,:)=uvmean(:,:,:)+tmean(:,:,:)

!uwmean=ux1*uz1
ta1(:,:,:)=ux1(:,:,:)*uz1(:,:,:)
call fine_to_coarseS(1,ta1,tmean)
uwmean(:,:,:)=uwmean(:,:,:)+tmean(:,:,:)

!vwmean=uy1*uz1
ta1(:,:,:)=uy1(:,:,:)*uz1(:,:,:)
call fine_to_coarseS(1,ta1,tmean)
vwmean(:,:,:)=vwmean(:,:,:)+tmean(:,:,:)

if (iscalar==1) then
   !phiphimean=phi1*phi1
   ta1(:,:,:)=phi1(:,:,:)*phi1(:,:,:)
   call fine_to_coarseS(1,ta1,tmean)
   phiphimean(:,:,:)=phiphimean(:,:,:)+tmean(:,:,:)
endif

!for a verification
!call decomp_2d_write_one(nx_global,ny_global,nz_global,&
!           1,ta1,'compa.dat')

if (mod(itime,isave)==0) then
   call decomp_2d_write_one(1,umean,'umean.dat',1)
   call decomp_2d_write_one(1,vmean,'vmean.dat',1)
   call decomp_2d_write_one(1,wmean,'wmean.dat',1)
   call decomp_2d_write_one(1,uumean,'uumean.dat',1)
   call decomp_2d_write_one(1,vvmean,'vvmean.dat',1)
   call decomp_2d_write_one(1,wwmean,'wwmean.dat',1)
   call decomp_2d_write_one(1,uvmean,'uvmean.dat',1)
   call decomp_2d_write_one(1,uwmean,'uwmean.dat',1)
   call decomp_2d_write_one(1,vwmean,'vwmean.dat',1)
   if (nrank==0) print *,'write stat arrays velocity done!'
   if (iscalar==1) then
      call decomp_2d_write_one(1,phimean,'phimean.dat',1)
      call decomp_2d_write_one(1,phiphimean,'phiphimean.dat',1)
      if (nrank==0) print *,'write stat arrays scalar done!'
   endif
!   call decomp_2d_write_one(nx_global,ny_global,nz_global,1,ux1,'compa.dat')

endif

end subroutine STATISTIC
