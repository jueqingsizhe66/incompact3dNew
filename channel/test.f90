PROGRAM incompact3d

USE decomp_2d
USE decomp_2d_poisson
use decomp_2d_io
USE variables
USE param
USE var
USE MPI
USE IBM

implicit none

integer :: code,nlock,i,j,k,bcx,bcy,bcz
real(mytype) :: x,y,z,tmp1
double precision :: t1,t2
character(len=6) :: filename

TYPE(DECOMP_INFO) :: phG,ph1,ph2,ph3,ph4

CALL MPI_INIT(code)
call decomp_2d_init(nx,ny,nz,p_row,p_col)
call init_coarser_mesh(nstat,nstat,nstat,.true.) !start from 1 == true
call parameter()

call init_variables

call schemes()

if (ifilter.eq.1) call filter()

if (nclx==0) then
   bcx=0
else
   bcx=1
endif
if (ncly==0) then
   bcy=0
else
   bcy=1
endif
if (nclz==0) then
   bcz=0
else
   bcz=1
endif

call decomp_2d_poisson_init(bcx,bcy,bcz)

call decomp_info_init(nxm,nym,nzm,phG)


if (ilit==0) call init(ux1,uy1,uz1,phi1,gx1,gy1,gz1,phis1,hx1,hy1,hz1,phiss1)  
if (ilit==1) call restart(ux1,uy1,uz1,pp3,phi1,gx1,gy1,gz1,&
        px1,py1,pz1,phis1,hx1,hy1,hz1,phiss1,phG,0)
call test_speed_min_max(ux1,uy1,uz1)
if (iscalar==1) call test_scalar_min_max(phi1)

!array for stat to zero
umean=0.;vmean=0.;wmean=0.
uumean=0.;vvmean=0.;wwmean=0.
uvmean=0.;uwmean=0.;vwmean=0.


t1 = MPI_WTIME()

!div: nx ny nz --> nxm ny nz --> nxm nym nz --> nxm nym nzm
call decomp_info_init(nxm, nym, nzm, ph1)
call decomp_info_init(nxm, ny, nz, ph4)

!gradp: nxm nym nzm -> nxm nym nz --> nxm ny nz --> nx ny nz
call decomp_info_init(nxm, ny, nz, ph2)  !for grad p
call decomp_info_init(nxm, nym, nz, ph3) !for grad p


istret=0
do itime=ifirst,ilast
   t=(itime-1)*dt
   if (nrank==0) then
      write(*,1001) itime,t
1001  format('Time step =',i7,', Time unit =',F9.3)
   endif
   
   do itr=1,iadvance_time

      if (nclx.eq.2) then
         call inflow (ux1,uy1,uz1,phi1) !X PENCILS
         call outflow(ux1,uy1,uz1,phi1) !X PENCILS 
      endif

     !X-->Y-->Z-->Y-->X
      call convdiff(ux1,uy1,uz1,ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1,&
           ux2,uy2,uz2,ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2,di2,&
           ux3,uy3,uz3,ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3)
           
      

      if (iscalar==1) then
         call scalar(ux1,uy1,uz1,phi1,phis1,phiss1,di1,tg1,th1,ti1,td1,&
     uy2,uz2,phi2,di2,ta2,tb2,tc2,td2,uz3,phi3,di3,ta3,tb3,ep1) 
      endif


      !X PENCILS
      call intt (ux1,uy1,uz1,gx1,gy1,gz1,hx1,hy1,hz1,ta1,tb1,tc1) 


      call pre_correc(ux1,uy1,uz1)


      if (ivirt==1) then !solid body old school
         !we are in X-pencil
         call corgp_IBM(ux1,uy1,uz1,px1,py1,pz1,1)
         call body(ux1,uy1,uz1,ep1)
         call corgp_IBM(ux1,uy1,uz1,px1,py1,pz1,2)
      endif

      !X-->Y-->Z
      call divergence (ux1,uy1,uz1,ep1,ta1,tb1,tc1,di1,td1,te1,tf1,&
           td2,te2,tf2,di2,ta2,tb2,tc2,ta3,tb3,tc3,di3,td3,te3,tf3,pp3,&
           nxmsize,nymsize,nzmsize,ph1,ph3,ph4,1)       

!      do k=ph1%zst(3),ph1%zen(3)
!         z=(real(k)-0.5)/real(nzm)
!          do j=ph1%zst(2),ph1%zen(2)
!             y=ypi(j)!(real(j)-0.5)/real(nym)
!             do i=ph1%zst(1),ph1%zen(1)
!              x=(real(i)-0.5)/real(nzm)
!              pp3(i,j,k) =cos(2.*pi*x)*cos(2.*pi*y)*cos(2.*pi*z) 
!(cos(2.*pi*z)+sin(2.*pi*z)+cos(4.*pi*z)+sin(4.*pi*z)&
!                   +cos(6.*pi*z)+sin(6.*pi*z))*&
!                   (cos(2.*pi*x)+sin(2.*pi*x)+cos(4.*pi*x)+sin(4.*pi*x)&
!                   +cos(6.*pi*x)+sin(6.*pi*x)+sin(8.*pi*x))*&
!                   (cos(1.*pi*y)+cos(2.*pi*y)+cos(3.*pi*y)+cos(4.*pi*y)&
!!                   +cos(5.*pi*y)+cos(6.*pi*y)+cos(7.*pi*y))
! !(cos(2.*pi*z)+sin(2.*pi*z)+cos(4.*pi*z)+sin(4.*pi*z)&
! !                  +cos(6.*pi*z)+sin(6.*pi*z))*&
! !                  cos(2.*pi*y)
! !                  (cos(1.*pi*y)+cos(2.*pi*y)+cos(3.*pi*y)+cos(4.*pi*y)&
! !                  +cos(5.*pi*y)+cos(6.*pi*y)+cos(7.*pi*y))
!                po3(i,j,k)=pp3(i,j,k)
!            end do
!          end do
!       end do


      !POISSON Z-->Z 
      call decomp_2d_poisson_stg(pp3,bcx,bcy,bcz)

      !Z-->Y-->X
      call gradp(px1,py1,pz1,di1,td2,tf2,ta2,tb2,tc2,di2,&
           ta3,tc3,di3,pp3,nxmsize,nymsize,nzmsize,ph2,ph3)

!      call divergence (px1,py1,pz1,ep1,ta1,tb1,tc1,di1,td1,te1,tf1,&
!           td2,te2,tf2,di2,ta2,tb2,tc2,ta3,tb3,tc3,di3,td3,te3,tf3,pp3,&
!           nxmsize,nymsize,nzmsize,ph1,ph3,ph4,1)

 !      tmp2=0.
 !      tmp3=0.
 !      tmp2=pp3(1,1,1)
 !      tmp3=pp4(1,1,1)
 !      do k=ph1%zst(3),ph1%zen(3)
 !         do j=ph1%zst(2),ph1%zen(2)
 !            do i=ph1%zst(1),ph1%zen(1)
 !               pp3(i,j,k)=pp3(i,j,k)-tmp2
 !               pp4(i,j,k)=pp4(i,j,k)-tmp3
 !            enddo^M
 !         enddo^M
 !      enddo^M
 !      tmp1=0.^M
!      do k=ph1%zst(3),ph1%zen(3)
!          do j=ph1%zst(2),ph1%zen(2)
!             do i=ph1%zst(1),ph1%zen(1)
!                tmp1=tmp1+(pp3(i,j,k)-po3(i,j,k))**2
!  !               if (nrank==0) write(*,*) i,j,k, pp3(i,j,k),pp4(i,j,k), pp3(i,j,k)-pp4(i,j,k)
!             enddo
!          enddo
!       enddo
!       print *,'ERROR DIV FINAL',nrank,tmp1/(ph1%zst(3)-ph1%zen(3))/&
!            (ph1%zst(2)-ph1%zen(2))/(ph1%zst(1)-ph1%zen(1))


      !X PENCILS
      call corgp(ux1,ux2,uy1,uz1,px1,py1,pz1) 
 
      
     !does not matter -->output=DIV U=0 (in dv3)
      call divergence (ux1,uy1,uz1,ep1,ta1,tb1,tc1,di1,td1,te1,tf1,&
           td2,te2,tf2,di2,ta2,tb2,tc2,ta3,tb3,tc3,di3,td3,te3,tf3,dv3,&
           nxmsize,nymsize,nzmsize,ph1,ph3,ph4,2)

      call test_speed_min_max(ux1,uy1,uz1)
      if (iscalar==1) call test_scalar_min_max(phi1)

   enddo

   call STATISTIC(ux1,uy1,uz1,ta1,umean,vmean,wmean,uumean,vvmean,wwmean,&
        uvmean,uwmean,vwmean,tmean)

   if (mod(itime,5)==0) call restart(ux1,uy1,uz1,pp3,phi1,gx1,gy1,gz1,&
        px1,py1,pz1,phis1,hx1,hy1,hz1,phiss1,phG,1)
     
   if (mod(itime,100)==0) then
      call Q_R (ux1,uy1,uz1,ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1,&
           ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2,di2,&
           ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3,phG)
   endif
enddo

t2=MPI_WTIME()-t1
call MPI_ALLREDUCE(t2,t1,1,MPI_REAL8,MPI_SUM, &
                   MPI_COMM_WORLD,code)
if (nrank==0) print *,'time per time_step: ', &
     t1/float(nproc)/(ilast-ifirst+1),' seconds'
if (nrank==0) print *,'simulation with nx*ny*nz=',nx,ny,nz,'mesh nodes'
if (nrank==0) print *,'Mapping p_row*p_col=',p_row,p_col


!call decomp_2d_poisson_finalize
call decomp_2d_finalize
CALL MPI_FINALIZE(code)

end PROGRAM incompact3d
