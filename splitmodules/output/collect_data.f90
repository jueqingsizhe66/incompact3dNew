subroutine collect_data()
!
!*****************************************************************

USE decomp_2d
USE decomp_2d_poisson
USE variables
USE param
USE var
USE MPI

implicit none

integer :: i,j,imin,ii,code
real(mytype), dimension(200) :: ta
integer,dimension(200) :: tai,tbi

!TOTAL TIME FOR COLLECTION T=50000*0.01=500
!I want 100 files
!100 numbers from 0 to 500

if (nrank==0) then
call random_number(ta)
do i=1,200
   tai(i)=int(ta(i)*2000./dt)
enddo
do j=1,200
imin=999998
ii=0
do i=1,200
   if (tai(i).le.imin) then
      imin=tai(i)
      ii=i
   endif
enddo
tbi(j)=imin
tai(ii)=999999
enddo
idata=tbi
do i=1,200
print *,i,idata(i),'VALUE COLLECTION DATA'
enddo
endif

call MPI_BCAST(idata,200,MPI_INTEGER,0,MPI_COMM_WORLD,code)

end subroutine collect_data
