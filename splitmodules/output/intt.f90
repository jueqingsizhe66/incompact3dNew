subroutine  intt (ux,uy,uz,gx,gy,gz,hx,hy,hz,ta1,tb1,tc1)
!
!********************************************************************

USE param
USE variables
USE decomp_2d

implicit none

integer :: ijk,nxyz
real,dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
real,dimension(xsize(1),xsize(2),xsize(3)) :: gx,gy,gz
real,dimension(xsize(1),xsize(2),xsize(3)) :: hx,hy,hz
real,dimension(xsize(1),xsize(2),xsize(3)) :: ta1,tb1,tc1
!real, dimension(xsize() :: ux,uy,uz,hx,hy,hz,gx,gy,gz

nxyz=xsize(1)*xsize(2)*xsize(3)

if ((nscheme.eq.1).or.(nscheme.eq.2)) then
if ((nscheme.eq.1.and.itime.eq.1.and.ilit.eq.0).or.&
     (nscheme.eq.2.and.itr.eq.1)) then
   do ijk=1,nxyz
      ux(ijk,1,1)=gdt(itr)*ta1(ijk,1,1)+ux(ijk,1,1)
      uy(ijk,1,1)=gdt(itr)*tb1(ijk,1,1)+uy(ijk,1,1)
      uz(ijk,1,1)=gdt(itr)*tc1(ijk,1,1)+uz(ijk,1,1)
      gx(ijk,1,1)=ta1(ijk,1,1)
      gy(ijk,1,1)=tb1(ijk,1,1)
      gz(ijk,1,1)=tc1(ijk,1,1)
   enddo
else
   if (nz.gt.1) then
      do ijk=1,nxyz
         ux(ijk,1,1)=adt(itr)*ta1(ijk,1,1)+bdt(itr)*gx(ijk,1,1)+ux(ijk,1,1)
         uy(ijk,1,1)=adt(itr)*tb1(ijk,1,1)+bdt(itr)*gy(ijk,1,1)+uy(ijk,1,1)
         uz(ijk,1,1)=adt(itr)*tc1(ijk,1,1)+bdt(itr)*gz(ijk,1,1)+uz(ijk,1,1)
         gx(ijk,1,1)=ta1(ijk,1,1)
         gy(ijk,1,1)=tb1(ijk,1,1)
         gz(ijk,1,1)=tc1(ijk,1,1)
      enddo
   else
      do ijk=1,nxyz
         ux(ijk,1,1)=adt(itr)*ta1(ijk,1,1)+bdt(itr)*gx(ijk,1,1)+ux(ijk,1,1)
         uy(ijk,1,1)=adt(itr)*tb1(ijk,1,1)+bdt(itr)*gy(ijk,1,1)+uy(ijk,1,1)
         gx(ijk,1,1)=ta1(ijk,1,1)
         gy(ijk,1,1)=tb1(ijk,1,1)
      enddo
   endif
endif
endif

if (nscheme.eq.3) then
   if (nz.gt.1) then
      if (adt(itr)==0.) then
         do ijk=1,nxyz
            gx(ijk,1,1)=dt*ta1(ijk,1,1)
            gy(ijk,1,1)=dt*tb1(ijk,1,1)
            gz(ijk,1,1)=dt*tc1(ijk,1,1)
         enddo
      else
         do ijk=1,nxyz
            gx(ijk,1,1)=adt(itr)*gx(ijk,1,1)+dt*ta1(ijk,1,1)
            gy(ijk,1,1)=adt(itr)*gy(ijk,1,1)+dt*tb1(ijk,1,1)
            gz(ijk,1,1)=adt(itr)*gz(ijk,1,1)+dt*tc1(ijk,1,1)
         enddo
      endif
      do ijk=1,nxyz
         ux(ijk,1,1)=ux(ijk,1,1)+bdt(itr)*gx(ijk,1,1)
         uy(ijk,1,1)=uy(ijk,1,1)+bdt(itr)*gy(ijk,1,1)
         uz(ijk,1,1)=uz(ijk,1,1)+bdt(itr)*gz(ijk,1,1)
      enddo
   else
      if (adt(itr)==0.) then
         do ijk=1,nxyz
            gx(ijk,1,1)=dt*ta1(ijk,1,1)
            gy(ijk,1,1)=dt*tb1(ijk,1,1)
         enddo
      else
         do ijk=1,nxyz
            gx(ijk,1,1)=adt(itr)*gx(ijk,1,1)+dt*ta1(ijk,1,1)
            gy(ijk,1,1)=adt(itr)*gy(ijk,1,1)+dt*tb1(ijk,1,1)
         enddo
      endif
      do ijk=1,nxyz
         ux(ijk,1,1)=ux(ijk,1,1)+bdt(itr)*gx(ijk,1,1)
         uy(ijk,1,1)=uy(ijk,1,1)+bdt(itr)*gy(ijk,1,1)
      enddo
   endif
endif

if (nscheme==4) then
   if ((itime.eq.1).and.(ilit.eq.0)) then
      if (nrank==0) print *,'start with Euler',itime
      do ijk=1,nxyz !start with Euler
         ux(ijk,1,1)=dt*ta1(ijk,1,1)+ux(ijk,1,1)
         uy(ijk,1,1)=dt*tb1(ijk,1,1)+uy(ijk,1,1)
         uz(ijk,1,1)=dt*tc1(ijk,1,1)+uz(ijk,1,1)
         gx(ijk,1,1)=ta1(ijk,1,1)
         gy(ijk,1,1)=tb1(ijk,1,1)
         gz(ijk,1,1)=tc1(ijk,1,1)
      enddo
   else
      if  ((itime.eq.2).and.(ilit.eq.0)) then
          if (nrank==0) print *,'then with AB2',itime
         do ijk=1,nxyz
            ux(ijk,1,1)=1.5*dt*ta1(ijk,1,1)-0.5*dt*gx(ijk,1,1)+ux(ijk,1,1)
            uy(ijk,1,1)=1.5*dt*tb1(ijk,1,1)-0.5*dt*gy(ijk,1,1)+uy(ijk,1,1)
            uz(ijk,1,1)=1.5*dt*tc1(ijk,1,1)-0.5*dt*gz(ijk,1,1)+uz(ijk,1,1)
            hx(ijk,1,1)=gx(ijk,1,1)
            hy(ijk,1,1)=gy(ijk,1,1)
            hz(ijk,1,1)=gz(ijk,1,1)
            gx(ijk,1,1)=ta1(ijk,1,1)
            gy(ijk,1,1)=tb1(ijk,1,1)
            gz(ijk,1,1)=tc1(ijk,1,1)
         enddo
      else
         do ijk=1,nxyz
            ux(ijk,1,1)=adt(itr)*ta1(ijk,1,1)+bdt(itr)*gx(ijk,1,1)+&
                 cdt(itr)*hx(ijk,1,1)+ux(ijk,1,1)
            uy(ijk,1,1)=adt(itr)*tb1(ijk,1,1)+bdt(itr)*gy(ijk,1,1)+&
                 cdt(itr)*hy(ijk,1,1)+uy(ijk,1,1)
            uz(ijk,1,1)=adt(itr)*tc1(ijk,1,1)+bdt(itr)*gz(ijk,1,1)+&
                 cdt(itr)*hz(ijk,1,1)+uz(ijk,1,1)
            hx(ijk,1,1)=gx(ijk,1,1)
            hy(ijk,1,1)=gy(ijk,1,1)
            hz(ijk,1,1)=gz(ijk,1,1)
            gx(ijk,1,1)=ta1(ijk,1,1)
            gy(ijk,1,1)=tb1(ijk,1,1)
            gz(ijk,1,1)=tc1(ijk,1,1)
         enddo
      endif
   endif
endif


return
end subroutine intt
