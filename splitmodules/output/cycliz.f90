subroutine cycliz(fiffz,fifz,ficz,fibz,fibbz,filaz,fiz1z,fiz2z,nz)
!
!*********************************************************************

USE parfiZ

implicit none

integer :: nz,k
real(mytype), dimension(nz) :: fiffz,fifz,ficz,fibz,fibbz,fiz1z,fiz2z
real(mytype), dimension(nz,2) :: filaz
real(mytype) :: xcoef, zdet

do k=1,nz
   fiz1z(k)=0.
   fiz2z(k)=0.
enddo
xcoef=2.
fiz1z(1)=xcoef
fiz1z(nz-1)=-1.
fiz2z(2)=-1.
fiz2z(nz)=xcoef
do k=1,nz-2
   fiz1z(k+1)=fiz1z(k+1)-filaz(k,1)*fiz1z(k)
   fiz1z(k+2)=fiz1z(k+2)-filaz(k,2)*fiz1z(k)
   fiz2z(k+1)=fiz2z(k+1)-filaz(k,1)*fiz2z(k)
   fiz2z(k+2)=fiz2z(k+2)-filaz(k,2)*fiz2z(k)
enddo
fiz1z(nz)=fiz1z(nz)-filaz(nz-1,1)*fiz1z(nz-1)
fiz2z(nz)=fiz2z(nz)-filaz(nz-1,1)*fiz2z(nz-1)
fiz1z(nz)=fiz1z(nz)/fiffz(nz)
fiz2z(nz)=fiz2z(nz)/fiffz(nz)
fiz1z(nz-1)=(fiz1z(nz-1)-fifz(nz-1)*fiz1z(nz))/fiffz(nz-1)
fiz2z(nz-1)=(fiz2z(nz-1)-fifz(nz-1)*fiz2z(nz))/fiffz(nz-1)
fiz1z(nz-2)=(fiz1z(nz-2)-fifz(nz-2)*fiz1z(nz-1)-&
     ficz(nz-2)*fiz1z(nz))/fiffz(nz-2)
fiz2z(nz-2)=(fiz2z(nz-2)-fifz(nz-2)*fiz2z(nz-1)-&
     ficz(nz-2)*fiz2z(nz))/fiffz(nz-2)
fiz1z(nz-3)=(fiz1z(nz-3)-fifz(nz-3)*fiz1z(nz-2)-&
     ficz(nz-3)*fiz1z(nz-1)-&
     fibz(nz-3)*fiz1z(nz))/fiffz(nz-3)
fiz2z(nz-3)=(fiz2z(nz-3)-fifz(nz-3)*fiz2z(nz-2)-&
     ficz(nz-3)*fiz2z(nz-1)-&
     fibz(nz-3)*fiz2z(nz))/fiffz(nz-3)
do k=nz-4,1,-1
   fiz1z(k)=(fiz1z(k)-fifz(k)*fiz1z(k+1)-&
        ficz(k)*fiz1z(k+2)-&
        fibz(k)*fiz1z(k+3)-&
        fibbz(k)*fiz1z(k+4))/fiffz(k)
   fiz2z(k)=(fiz2z(k)-fifz(k)*fiz2z(k+1)-&
        ficz(k)*fiz2z(k+2)-&
        fibz(k)*fiz2z(k+3)-&
        fibbz(k)*fiz2z(k+4))/fiffz(k)
enddo
zdet=(1.-fibez*fiz1z(1)+fiz1z(nz-1)*fibez/xcoef+&
     fialz*fiz1z(nz)/xcoef)*&
     (1.+fialz*fiz2z(1)/xcoef+fibez*fiz2z(2)/xcoef-&
     fibez*fiz2z(nz))-&
     (-fibez*fiz2z(1)+fiz2z(nz-1)*fibez/xcoef+&
     fialz*fiz2z(nz)/xcoef)*&
     (fialz*fiz1z(1)/xcoef+fibez*fiz1z(2)/xcoef-&
     fibez*fiz1z(nz))

fih1z=(1.+fialz*fiz2z(1)/xcoef+fibez*fiz2z(2)/xcoef-&
     fibez*fiz2z(nz))/zdet
fih2z=-(fialz*fiz1z(1)/xcoef+fibez*fiz1z(2)/xcoef-&
     fibez*fiz1z(nz))/zdet
fih3z=-(-fibez*fiz2z(1)+fiz2z(nz-1)*fibez/xcoef+&
     fialz*fiz2z(nz)/xcoef)/zdet
fih4z=(1.-fibez*fiz1z(1)+fiz1z(nz-1)*fibez/xcoef+&
     fialz*fiz1z(nz)/xcoef)/zdet


return
end subroutine cycliz
