subroutine cyclix(fiffx,fifx,ficx,fibx,fibbx,filax,fiz1x,fiz2x,nx)
!
!*********************************************************************

USE parfiX

implicit none

integer :: nx,i
real(mytype), dimension(nx) :: fiffx,fifx,ficx,fibx,fibbx,fiz1x,fiz2x
real(mytype), dimension(nx,2) :: filax
real(mytype) :: xcoef, xdet

do i=1,nx
   fiz1x(i)=0.
   fiz2x(i)=0.
enddo
xcoef=2.
fiz1x(1)=xcoef
fiz1x(nx-1)=-1
fiz2x(2)=-1.
fiz2x(nx)=xcoef
do i=1,nx-2
   fiz1x(i+1)=fiz1x(i+1)-filax(i,1)*fiz1x(i)
   fiz1x(i+2)=fiz1x(i+2)-filax(i,2)*fiz1x(i)
   fiz2x(i+1)=fiz2x(i+1)-filax(i,1)*fiz2x(i)
   fiz2x(i+2)=fiz2x(i+2)-filax(i,2)*fiz2x(i)
enddo

fiz1x(nx)=fiz1x(nx)-filax(nx-1,1)*fiz1x(nx-1)
fiz2x(nx)=fiz2x(nx)-filax(nx-1,1)*fiz2x(nx-1)

fiz1x(nx)=fiz1x(nx)/fiffx(nx)
fiz2x(nx)=fiz2x(nx)/fiffx(nx)
fiz1x(nx-1)=(fiz1x(nx-1)-fifx(nx-1)*fiz1x(nx))/fiffx(nx-1)
fiz2x(nx-1)=(fiz2x(nx-1)-fifx(nx-1)*fiz2x(nx))/fiffx(nx-1)
fiz1x(nx-2)=(fiz1x(nx-2)-fifx(nx-2)*fiz1x(nx-1)-&
     ficx(nx-2)*fiz1x(nx))/fiffx(nx-2)
fiz2x(nx-2)=(fiz2x(nx-2)-fifx(nx-2)*fiz2x(nx-1)-&
     ficx(nx-2)*fiz2x(nx))/fiffx(nx-2)
fiz1x(nx-3)=(fiz1x(nx-3)-fifx(nx-3)*fiz1x(nx-2)-&
     ficx(nx-3)*fiz1x(nx-1)-&
     fibx(nx-3)*fiz1x(nx))/fiffx(nx-3)
fiz2x(nx-3)=(fiz2x(nx-3)-fifx(nx-3)*fiz2x(nx-2)-&
     ficx(nx-3)*fiz2x(nx-1)-&
     fibx(nx-3)*fiz2x(nx))/fiffx(nx-3)
do i=nx-4,1,-1
   fiz1x(i)=(fiz1x(i)-fifx(i)*fiz1x(i+1)-&
        ficx(i)*fiz1x(i+2)-&
        fibx(i)*fiz1x(i+3)-&
        fibbx(i)*fiz1x(i+4))/fiffx(i)
   fiz2x(i)=(fiz2x(i)-fifx(i)*fiz2x(i+1)-&
        ficx(i)*fiz2x(i+2)-&
        fibx(i)*fiz2x(i+3)-&
        fibbx(i)*fiz2x(i+4))/fiffx(i)
enddo
xdet=(1.-fibex*fiz1x(1)+fiz1x(nx-1)*fibex/xcoef+&
     fialx*fiz1x(nx)/xcoef)*&
     (1.+fialx*fiz2x(1)/xcoef+fibex*fiz2x(2)/xcoef-&
     fibex*fiz2x(nx))-&
     (-fibex*fiz2x(1)+fiz2x(nx-1)*fibex/xcoef+&
     fialx*fiz2x(nx)/xcoef)*&
     (fialx*fiz1x(1)/xcoef+fibex*fiz1x(2)/xcoef-&
     fibex*fiz1x(nx))

fih1x=(1.+fialx*fiz2x(1)/xcoef+fibex*fiz2x(2)/xcoef-&
     fibex*fiz2x(nx))/xdet
fih2x=-(fialx*fiz1x(1)/xcoef+fibex*fiz1x(2)/xcoef-&
     fibex*fiz1x(nx))/xdet
fih3x=-(-fibex*fiz2x(1)+fiz2x(nx-1)*fibex/xcoef+&
     fialx*fiz2x(nx)/xcoef)/xdet
fih4x=(1.-fibex*fiz1x(1)+fiz1x(nx-1)*fibex/xcoef+&
     fialx*fiz1x(nx)/xcoef)/xdet


return
end subroutine cyclix
