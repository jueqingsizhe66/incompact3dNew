subroutine pre_correc(ux,uy,uz)
!
!****************************************************************************

USE decomp_2d
USE variables
USE param
USE var
USE MPI


implicit none

real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
integer :: i,j,k,code
real(mytype) :: ut,ut1,utt,ut11

if (itime==1) then
   dpdyx1=0.
   dpdzx1=0.
   dpdyxn=0.
   dpdzxn=0.
endif


!we are in X pencils:
do k=1,xsize(3)
do j=1,xsize(2)
   dpdyx1(j,k)=dpdyx1(j,k)*gdt(itr)
   dpdzx1(j,k)=dpdzx1(j,k)*gdt(itr)
   dpdyxn(j,k)=dpdyxn(j,k)*gdt(itr)
   dpdzxn(j,k)=dpdzxn(j,k)*gdt(itr)
enddo
enddo

if (xsize(3)==1) then
   do j=1,xsize(2)
   do i=1,xsize(1)
      dpdxz1(i,j)=dpdxz1(i,j)*gdt(itr)
      dpdyz1(i,j)=dpdyz1(i,j)*gdt(itr)
   enddo
   enddo
endif
if (xsize(3)==nz) then
   do j=1,xsize(2)
   do i=1,xsize(1)
      dpdxzn(i,j)=dpdxzn(i,j)*gdt(itr)
      dpdyzn(i,j)=dpdyzn(i,j)*gdt(itr)
   enddo
   enddo
endif

if (xsize(2)==1) then
   do k=1,xsize(3)
   do i=1,xsize(1)
      dpdxy1(i,k)=dpdxy1(i,k)*gdt(itr)
      dpdzy1(i,k)=dpdzy1(i,k)*gdt(itr)
   enddo
   enddo
endif
if (xsize(2)==ny) then
   do k=1,xsize(3)
   do i=1,xsize(1)
      dpdxyn(i,k)=dpdxyn(i,k)*gdt(itr)
      dpdzyn(i,k)=dpdzyn(i,k)*gdt(itr)
   enddo
   enddo
endif

!Computatation of the flow rate Inflow/Outflow
!we are in X pencils:
if (nclx==2) then
   ut1=0.
   do k=1,xsize(3)
   do j=1,xsize(2)
      ut1=ut1+bxx1(j,k)
   enddo
   enddo
   ut1=ut1/xsize(2)/xsize(3)
   call MPI_ALLREDUCE(ut1,ut11,1,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,code)
   ut11=ut11/nproc
   ut=0.
   do k=1,xsize(3)
   do j=1,xsize(2)
      ut=ut+bxxn(j,k)
   enddo
   enddo
   ut=ut/xsize(2)/xsize(3)
   call MPI_ALLREDUCE(ut,utt,1,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,code)
   utt=utt/nproc

   if (nrank==0) print *,'FLOW RATE I/O',ut11,utt

   do k=1,xsize(3)
   do j=1,xsize(2)
      bxxn(j,k)=bxxn(j,k)-utt+ut11
   enddo
   enddo
endif

!********NCLX==2*************************************
!****************************************************
if (nclx.eq.2) then
   do k=1,xsize(3)
   do j=1,xsize(2)
      ux(1 ,j,k)=bxx1(j,k)
      uy(1 ,j,k)=bxy1(j,k)+dpdyx1(j,k)
      uz(1 ,j,k)=bxz1(j,k)+dpdzx1(j,k)
      ux(nx,j,k)=bxxn(j,k)
      uy(nx,j,k)=bxyn(j,k)+dpdyxn(j,k)
      uz(nx,j,k)=bxzn(j,k)+dpdzxn(j,k)
   enddo
   enddo
endif
!****************************************************
!********NCLY==2*************************************
!****************************************************
!WE ARE IN X PENCIL!!!!!!
if (ncly==2) then
   if (itype.eq.2) then
!find j=1 and j=ny
      if (xstart(2)==1) then
         do k=1,xsize(3)
         do i=1,xsize(1)
            ux(i,1,k)=0.+dpdxy1(i,k)
            uy(i,1,k)=0.
            uz(i,1,k)=0.+dpdzy1(i,k)
         enddo
         enddo
      endif
      if (xend(2)==ny) then
         do k=1,xsize(3)
         do i=1,xsize(1)
            ux(i,ny,k)=0.+dpdxyn(i,k)
            uy(i,ny,k)=0.
            uz(i,ny,k)=0.+dpdzyn(i,k)
         enddo
         enddo
      endif
   endif
endif
!****************************************************
!********NCLZ==2*************************************
!****************************************************
!****************************************************

!#####################################################

return
end subroutine pre_correc
