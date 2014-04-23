subroutine cycliy(fiffy,fify,ficy,fiby,fibby,filay,fiz1y,fiz2y,ny)
!
!*********************************************************************

USE parfiY

implicit none
integer :: ny,j
real(mytype), dimension(ny) :: fiffy,fify,ficy,fiby,fibby,fiz1y,fiz2y
real(mytype), dimension(ny,2) :: filay
real(mytype) :: xcoef, ydet

do j=1,ny
   fiz1y(j)=0.
   fiz2y(j)=0.
enddo
xcoef=2.
fiz1y(1)=xcoef
fiz1y(ny-1)=-1.
fiz2y(2)=-1.
fiz2y(ny)=xcoef

do j=1,ny-2
   fiz1y(j+1)=fiz1y(j+1)-filay(j,1)*fiz1y(j)
   fiz1y(j+2)=fiz1y(j+2)-filay(j,2)*fiz1y(j)
   fiz2y(j+1)=fiz2y(j+1)-filay(j,1)*fiz2y(j)
   fiz2y(j+2)=fiz2y(j+2)-filay(j,2)*fiz2y(j)
enddo
fiz1y(ny)=fiz1y(ny)-filay(ny-1,1)*fiz1y(ny-1)
fiz2y(ny)=fiz2y(ny)-filay(ny-1,1)*fiz2y(ny-1)
fiz1y(ny)=fiz1y(ny)/fiffy(ny)
fiz2y(ny)=fiz2y(ny)/fiffy(ny)
fiz1y(ny-1)=(fiz1y(ny-1)-fify(ny-1)*fiz1y(ny))/fiffy(ny-1)
fiz2y(ny-1)=(fiz2y(ny-1)-fify(ny-1)*fiz2y(ny))/fiffy(ny-1)
fiz1y(ny-2)=(fiz1y(ny-2)-fify(ny-2)*fiz1y(ny-1)-&
     ficy(ny-2)*fiz1y(ny))/fiffy(ny-2)
fiz2y(ny-2)=(fiz2y(ny-2)-fify(ny-2)*fiz2y(ny-1)-&
     ficy(ny-2)*fiz2y(ny))/fiffy(ny-2)
fiz1y(ny-3)=(fiz1y(ny-3)-fify(ny-3)*fiz1y(ny-2)-&
     ficy(ny-3)*fiz1y(ny-1)-&
     fiby(ny-3)*fiz1y(ny))/fiffy(ny-3)
fiz2y(ny-3)=(fiz2y(ny-3)-fify(ny-3)*fiz2y(ny-2)-&
     ficy(ny-3)*fiz2y(ny-1)-&
     fiby(ny-3)*fiz2y(ny))/fiffy(ny-3)
do j=ny-4,1,-1
   fiz1y(j)=(fiz1y(j)-fify(j)*fiz1y(j+1)-&
        ficy(j)*fiz1y(j+2)-&
        fiby(j)*fiz1y(j+3)-&
        fibby(j)*fiz1y(j+4))/fiffy(j)
   fiz2y(j)=(fiz2y(j)-fify(j)*fiz2y(j+1)-&
        ficy(j)*fiz2y(j+2)-&
        fiby(j)*fiz2y(j+3)-&
                fibby(j)*fiz2y(j+4))/fiffy(j)
enddo
ydet=(1.-fibey*fiz1y(1)+fibey*fiz1y(ny-1)/xcoef+&
     fialy*fiz1y(ny)/xcoef)*&
     (1.+fialy*fiz2y(1)/xcoef+fibey*fiz2y(2)/xcoef-&
     fibey*fiz2y(ny))-&
     (-fibey*fiz2y(1)+fibey*fiz2y(ny-1)/xcoef+&
     fialy*fiz2y(ny)/xcoef)*&
     (fialy*fiz1y(1)/xcoef+fibey*fiz1y(2)/xcoef-&
     fibey*fiz1y(ny))

fih1y=(1.+fialy*fiz2y(1)/xcoef+fibey*fiz2y(2)/xcoef-&
     fibey*fiz2y(ny))/ydet
fih2y=-(fialy*fiz1y(1)/xcoef+fibey*fiz1y(2)/xcoef-&
     fibey*fiz1y(ny))/ydet
fih3y=-(-fibey*fiz2y(1)+fibey*fiz2y(ny-1)/xcoef+&
     fialy*fiz2y(ny)/xcoef)/ydet
fih4y=(1.-fibey*fiz1y(1)+fibey*fiz1y(ny-1)/xcoef+&
     fialy*fiz1y(ny)/xcoef)/ydet


return
end subroutine cycliy
