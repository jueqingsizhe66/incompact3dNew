subroutine divergence (ux1,uy1,uz1,ep1,ta1,tb1,tc1,di1,td1,te1,tf1,&
     td2,te2,tf2,di2,ta2,tb2,tc2,ta3,tb3,tc3,di3,td3,te3,tf3,pp3,&
     nxmsize,nymsize,nzmsize,ph1,ph3,ph4,nlock)
!
!********************************************************************

USE param
USE IBM
USE decomp_2d
USE variables
USE MPI

implicit none

TYPE(DECOMP_INFO) :: ph1,ph3,ph4

!X PENCILS NX NY NZ  -->NXM NY NZ
real,dimension(xsize(1),xsize(2),xsize(3)) :: ta1,tb1,tc1,di1,ux1,uy1,uz1,ep1
real,dimension(nxmsize,xsize(2),xsize(3)) :: td1,te1,tf1
!Y PENCILS NXM NY NZ  -->NXM NYM NZ
real,dimension(ph1%yst(1):ph1%yen(1),ysize(2),ysize(3)) :: td2,te2,tf2,di2
real,dimension(ph1%yst(1):ph1%yen(1),nymsize,ysize(3)) :: ta2,tb2,tc2
!Z PENCILS NXM NYM NZ  -->NXM NYM NZM
real,dimension(ph1%zst(1):ph1%zen(1),ph1%zst(2):ph1%zen(2),zsize(3)) :: ta3,tb3,tc3,di3
real,dimension(ph1%zst(1):ph1%zen(1),ph1%zst(2):ph1%zen(2),ph1%zst(3):ph1%zen(3)) :: td3,te3,tf3,pp3



integer :: ijk,nvect1,nvect2,nvect3,i,j,k,nlock
integer :: nxmsize,nymsize,nzmsize,code
real :: tmax,tmoy,tmax1,tmoy1


nvect1=xsize(1)*xsize(2)*xsize(3)
nvect2=ysize(1)*ysize(2)*ysize(3)


if (nlock==1) then
   if (ivirt.eq.0) ep1(:,:,:)=0.
   do ijk=1,nvect1
      ta1(ijk,1,1)=(1.-ep1(ijk,1,1))*ux1(ijk,1,1)
      tb1(ijk,1,1)=(1.-ep1(ijk,1,1))*uy1(ijk,1,1)
      tc1(ijk,1,1)=(1.-ep1(ijk,1,1))*uz1(ijk,1,1)
   enddo
else
   ta1(:,:,:)=ux1(:,:,:)
   tb1(:,:,:)=uy1(:,:,:)
   tc1(:,:,:)=uz1(:,:,:)
endif



!WORK X-PENCILS
call decx6(td1,ta1,di1,sx,cfx6,csx6,cwx6,xsize(1),nxmsize,xsize(2),xsize(3),0)
call inter6(te1,tb1,di1,sx,cifxp6,cisxp6,ciwxp6,xsize(1),nxmsize,xsize(2),xsize(3),1)
call inter6(tf1,tc1,di1,sx,cifxp6,cisxp6,ciwxp6,xsize(1),nxmsize,xsize(2),xsize(3),1)

call transpose_x_to_y(td1,td2,ph4)!->NXM NY NZ
call transpose_x_to_y(te1,te2,ph4)
call transpose_x_to_y(tf1,tf2,ph4)


!WORK Y-PENCILS

call intery6(ta2,td2,di2,sy,cifyp6,cisyp6,ciwyp6,(ph1%yen(1)-ph1%yst(1)+1),ysize(2),nymsize,ysize(3)
call decy6(tb2,te2,di2,sy,cfy6,csy6,cwy6,ppyi,(ph1%yen(1)-ph1%yst(1)+1),ysize(2),nymsize,ysize(3),0)
call intery6(tc2,tf2,di2,sy,cifyp6,cisyp6,ciwyp6,(ph1%yen(1)-ph1%yst(1)+1),ysize(2),nymsize,ysize(3)

call transpose_y_to_z(ta2,ta3,ph3)!->NXM NYM NZ
call transpose_y_to_z(tb2,tb3,ph3)
call transpose_y_to_z(tc2,tc3,ph3)


!WORK Z-PENCILS
call interz6(td3,ta3,di3,sz,cifzp6,ciszp6,ciwzp6,(ph1%zen(1)-ph1%zst(1)+1),&
     (ph1%zen(2)-ph1%zst(2)+1),zsize(3),(ph1%zen(3)-ph1%zen(3)+1),1)
call interz6(te3,tb3,di3,sz,cifzp6,ciszp6,ciwzp6,(ph1%zen(1)-ph1%zst(1)+1),&
     (ph1%zen(2)-ph1%zst(2)+1),zsize(3),(ph1%zen(3)-ph1%zen(3)+1),1)
call decz6(tf3,tc3,di3,sz,cfz6,csz6,cwz6,(ph1%zen(1)-ph1%zst(1)+1),&
     (ph1%zen(2)-ph1%zst(2)+1),zsize(3),(ph1%zen(3)-ph1%zen(3)+1),0)



do k=ph1%zst(3),ph1%zen(3)
do j=ph1%zst(2),ph1%zen(2)
do i=ph1%zst(1),ph1%zen(1)
pp3(i,j,k)=td3(i,j,k)+te3(i,j,k)+tf3(i,j,k)
enddo
enddo
enddo

nvect3=(ph1%zen(1)-ph1%zst(1)+1)*(ph1%zen(2)-ph1%zst(2)+1)*(ph1%zen(3)-ph1%zst(3)+1)


tmax=-1609.
tmoy=0.
do k=ph1%zst(3),ph1%zen(3)
do j=ph1%zst(2),ph1%zen(2)
do i=ph1%zst(1),ph1%zen(1)
   if (pp3(i,j,k).gt.tmax) tmax=pp3(i,j,k)
   tmoy=tmoy+abs(pp3(i,j,k))
enddo
enddo
enddo
tmoy=tmoy/nvect3
!print *,'TMOY PER RANK',nrank,tmoy,tmax

call MPI_REDUCE(tmax,tmax1,1,MPI_REAL,MPI_MAX,0,MPI_COMM_WORLD,code)
call MPI_REDUCE(tmoy,tmoy1,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,code)!

if (nrank==0) then
     if (nlock==2) then
        print *,'DIV U final Max=',tmax1
        print *,'DIV U final Moy=',tmoy1/real(nproc)
     else
        print *,'DIV U* Max=',tmax1
        print *,'DIV U* Moyy=',tmoy1/real(nproc)
     endif
endif

return
end subroutine divergence
