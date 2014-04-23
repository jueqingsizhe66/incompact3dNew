subroutine gradp(ta1,tb1,tc1,di1,td2,tf2,ta2,tb2,tc2,di2,&
     ta3,tc3,di3,pp3,nxmsize,nymsize,nzmsize,ph2,ph3)
!
!********************************************************************

USE param
USE decomp_2d
USE variables

implicit none

TYPE(DECOMP_INFO) :: ph2,ph3
integer :: i,j,k,ijk,nxmsize,nymsize,nzmsize

real,dimension(ph3%zst(1):ph3%zen(1),ph3%zst(2):ph3%zen(2),nzmsize) :: pp3
!Z PENCILS NXM NYM NZM-->NXM NYM NZ
real,dimension(ph3%zst(1):ph3%zen(1),ph3%zst(2):ph3%zen(2),zsize(3)) :: ta3,tc3,di3
!Y PENCILS NXM NYM NZ -->NXM NY NZ
real,dimension(ph3%yst(1):ph3%yen(1),nymsize,ysize(3)) :: ta2,tc2
real,dimension(ph3%yst(1):ph3%yen(1),ysize(2),ysize(3)) :: tb2,td2,tf2,di2
!X PENCILS NXM NY NZ  -->NX NY NZ
real,dimension(nxmsize,xsize(2),xsize(3)) :: td1,te1,tf1
real,dimension(xsize(1),xsize(2),xsize(3)) :: ta1,tb1,tc1,di1




!WORK Z-PENCILS

call interiz6(ta3,pp3,di3,sz,cifip6z,cisip6z,ciwip6z,cifz6,cisz6,ciwz6,&
     (ph3%zen(1)-ph3%zst(1)+1),(ph3%zen(2)-ph3%zst(2)+1),nzmsize,zsize(3),1)
call deciz6(tc3,pp3,di3,sz,cfip6z,csip6z,cwip6z,cfz6,csz6,cwz6,&
     (ph3%zen(1)-ph3%zst(1)+1),(ph3%zen(2)-ph3%zst(2)+1),nzmsize,zsize(3),1)

!     do k=1,zsize(3)
!        do j=ph3%zst(2),ph3%zen(2)
!            do i=ph3%zst(1),ph3%zen(1)
!               print *,i,j,k,ta3(i,j,k)-pp3(i,j,k)
!            enddo
!        enddo
!      enddo

!      if (nrank==1) then
!         open(10,file='testA.dat',status='unknown',form='formatted')
!         do k=1,nz
!         write(10,*) (real(k)-0.5)*dz,ta3(1,j,1),pp3(1,j,1),tc2(1,j,1)
!        enddo
!        close(10)
!     endif

!WORK Y-PENCILS
call transpose_z_to_y(ta3,ta2,ph3) !nxm nym nz
call transpose_z_to_y(tc3,tc2,ph3)



!print *,'nrank',nrank,ph3%yst(1),ph3%yen(1)
!
!do k=1,ysize(3)
!do j=1,nym
!do i=ph3%yst(1),ph3%yen(1)
!   ta2(i,j,k)=cos(2.*pi*ypi(j))
!enddo
!enddo
!enddo

call interiy6(tb2,ta2,di2,sy,cifip6y,cisip6y,ciwip6y,cify6,cisy6,ciwy6,&
     (ph3%yen(1)-ph3%yst(1)+1),nymsize,ysize(2),ysize(3),1)
call deciy6(td2,ta2,di2,sy,cfip6y,csip6y,cwip6y,cfy6,csy6,cwy6,ppy,&
     (ph3%yen(1)-ph3%yst(1)+1),nymsize,ysize(2),ysize(3),1)
call interiy6(tf2,tc2,di2,sy,cifip6y,cisip6y,ciwip6y,cify6,cisy6,ciwy6,&
     (ph3%yen(1)-ph3%yst(1)+1),nymsize,ysize(2),ysize(3),1)

!      if (nrank==0) then
!         open(10,file='test1.dat',status='unknown',form='formatted')
!         do j=1,ny
!         write(10,*) yp(j),tb2(1,j,1),cos(2.*pi*yp(j))
!         enddo
!        close(10)
!     endif!


!      if (nrank==1) then
!         open(10,file='testA.dat',status='unknown',form='formatted')
!         do j=1,ny
!         write(10,*) (real(j)-1.)*dy,tb2(1,j,1),td2(1,j,1),tf2(1,j,1)
!        enddo
!        close(10)
!     endif

!print *,tb2
!     do k=1,ysize(3)
!         do j=1,ysize(2)
!            do i=ph3%yst(1),ph3%yen(1)
!               print *,i,j,k,tb2(i,j,k)
!            enddo
!         enddo
!      enddo

!WORK X-PENCILS

call transpose_y_to_x(tb2,td1,ph2) !nxm ny nz
call transpose_y_to_x(td2,te1,ph2)
call transpose_y_to_x(tf2,tf1,ph2)

!print *,te1

!     do k=ph3%xst(3),ph3%xen(3)
!         do j=ph3%xst(2),ph3%xen(2)
!            do i=ph3%xst(1),ph3%xen(1)
!               print *,i,j,k,te1(i,j,k)
!            enddo
!         enddo
!      enddo

call deci6(ta1,td1,di1,sx,cfip6,csip6,cwip6,cfx6,csx6,cwx6,&
     nxmsize,xsize(1),xsize(2),xsize(3),1)
call interi6(tb1,te1,di1,sx,cifip6,cisip6,ciwip6,cifx6,cisx6,ciwx6,&
     nxmsize,xsize(1),xsize(2),xsize(3),1)
call interi6(tc1,tf1,di1,sx,cifip6,cisip6,ciwip6,cifx6,cisx6,ciwx6,&
     nxmsize,xsize(1),xsize(2),xsize(3),1)



!we are in X pencils:
do k=1,xsize(3)
do j=1,xsize(2)
   dpdyx1(j,k)=tb1(1,j,k)/gdt(itr)
   dpdzx1(j,k)=tc1(1,j,k)/gdt(itr)
   dpdyxn(j,k)=tb1(nx,j,k)/gdt(itr)
   dpdzxn(j,k)=tc1(nx,j,k)/gdt(itr)
enddo
enddo


if (xsize(3)==1) then
   do j=1,xsize(2)
   do i=1,xsize(1)
      dpdxz1(i,j)=ta1(i,j,1)/gdt(itr)
      dpdyz1(i,j)=tb1(i,j,1)/gdt(itr)
   enddo
   enddo
endif
if (xsize(3)==nz) then
   do j=1,xsize(2)
   do i=1,xsize(1)
      dpdxzn(i,j)=ta1(i,j,nz)/gdt(itr)
      dpdyzn(i,j)=tb1(i,j,nz)/gdt(itr)
   enddo
   enddo
endif

if (xsize(2)==1) then
   do k=1,xsize(3)
   do i=1,xsize(1)
      dpdxy1(i,k)=ta1(i,1,k)/gdt(itr)
      dpdzy1(i,k)=tc1(i,1,k)/gdt(itr)
   enddo
   enddo
endif
if (xsize(2)==ny) then
   do k=1,xsize(3)
   do i=1,xsize(1)
      dpdxyn(i,k)=ta1(i,ny,k)/gdt(itr)
      dpdzyn(i,k)=tc1(i,ny,k)/gdt(itr)
   enddo
   enddo
endif

return
end subroutine gradp
