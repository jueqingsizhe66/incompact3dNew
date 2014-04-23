subroutine inversion5_v1(aaa,eee,spI)
!
!*****************************************************************

USE decomp_2d
USE decomp_2d_poisson
USE variables
USE param
USE var
USE MPI

implicit none

! decomposition object for spectral space
TYPE(DECOMP_INFO) :: spI

#ifdef DOUBLE_PREC
  real(mytype), parameter :: epsilon = 1.e-16
#else
  real(mytype), parameter :: epsilon = 1.e-8
#endif

complex(mytype),dimension(spI%yst(1):spI%yen(1),ny/2,spI%yst(3):spI%yen(3),5) :: aaa
complex(mytype),dimension(spI%yst(1):spI%yen(1),spI%yst(2):spI%yen(2),spI%yst(3):spI%yen(3)) :: eee
integer :: i,j,k,m,mi,jc
integer,dimension(2) :: ja,jb
complex(mytype),dimension(spI%yst(1):spI%yen(1),spI%yst(3):spI%yen(3)) :: sr
complex(mytype),dimension(spI%yst(1):spI%yen(1),spI%yst(3):spI%yen(3)) :: a1,b1

real(mytype) :: tmp1,tmp2,tmp3,tmp4

do i=1,2
   ja(i)=4-i
   jb(i)=5-i
enddo
do m=1,ny/2-2
   do i=1,2
      mi=m+i
      do k=spI%yst(3),spI%yen(3)
      do j=spI%yst(1),spI%yen(1)
         if (real(aaa(j,m,k,3), kind=mytype).ne.0.) tmp1=real(aaa(j,mi,k,3-i), kind=mytype)/real(aaa
         if (aimag(aaa(j,m,k,3)).ne.0.)tmp2=aimag(aaa(j,mi,k,3-i))/aimag(aaa(j,m,k,3))
         sr(j,k)=cmplx(tmp1,tmp2, kind=mytype)
         eee(j,mi,k)=cmplx(real(eee(j,mi,k), kind=mytype)-tmp1*real(eee(j,m,k), kind=mytype),&
              aimag(eee(j,mi,k))-tmp2*aimag(eee(j,m,k)), kind=mytype)
      enddo
      enddo
      do jc=ja(i),jb(i)
      do k=spI%yst(3),spI%yen(3)
      do j=spI%yst(1),spI%yen(1)
         aaa(j,mi,k,jc)=cmplx(real(aaa(j,mi,k,jc), kind=mytype)-real(sr(j,k), kind=mytype)*real(aaa(
              aimag(aaa(j,mi,k,jc))-aimag(sr(j,k))*aimag(aaa(j,m,k,jc+i)), kind=mytype)
      enddo
      enddo
      enddo
   enddo
enddo


do k=spI%yst(3),spI%yen(3)
do j=spI%yst(1),spI%yen(1)
   if (abs(real(aaa(j,ny/2-1,k,3), kind=mytype)).gt.epsilon) then
      tmp1=real(aaa(j,ny/2,k,2), kind=mytype)/real(aaa(j,ny/2-1,k,3), kind=mytype)
   else
      tmp1=0.
   endif
   if (abs(aimag(aaa(j,ny/2-1,k,3))).gt.epsilon) then
      tmp2=aimag(aaa(j,ny/2,k,2))/aimag(aaa(j,ny/2-1,k,3))
   else
      tmp2=0.
   endif
   sr(j,k)=cmplx(tmp1,tmp2, kind=mytype)
   b1(j,k)=cmplx(real(aaa(j,ny/2,k,3), kind=mytype)-tmp1*real(aaa(j,ny/2-1,k,4), kind=mytype),&
        aimag(aaa(j,ny/2,k,3))-tmp2*aimag(aaa(j,ny/2-1,k,4)), kind=mytype)

   if (abs(real(b1(j,k), kind=mytype)).gt.epsilon) then
      tmp1=real(sr(j,k), kind=mytype)/real(b1(j,k), kind=mytype)
      tmp3=real(eee(j,ny/2,k), kind=mytype)/real(b1(j,k), kind=mytype)-tmp1*real(eee(j,ny/2-1,k), ki
   else
      tmp1=0.
      tmp3=0.
   endif
   if (abs(aimag(b1(j,k))).gt.epsilon) then
      tmp2=aimag(sr(j,k))/aimag(b1(j,k))
      tmp4=aimag(eee(j,ny/2,k))/aimag(b1(j,k))-tmp2*aimag(eee(j,ny/2-1,k))
   else
      tmp2=0.
      tmp4=0.
   endif
   a1(j,k)=cmplx(tmp1,tmp2, kind=mytype)
   eee(j,ny/2,k)=cmplx(tmp3,tmp4, kind=mytype)

   if (abs(real(aaa(j,ny/2-1,k,3), kind=mytype)).gt.epsilon) then
      tmp1=1./real(aaa(j,ny/2-1,k,3), kind=mytype)
   else
      tmp1=0.
   endif
   if (abs(aimag(aaa(j,ny/2-1,k,3))).gt.epsilon) then
      tmp2=1./aimag(aaa(j,ny/2-1,k,3))
   else
      tmp2=0.
   endif
   b1(j,k)=cmplx(tmp1,tmp2, kind=mytype)
   a1(j,k)=cmplx(real(aaa(j,ny/2-1,k,4), kind=mytype)*real(b1(j,k), kind=mytype),&
        aimag(aaa(j,ny/2-1,k,4))*aimag(b1(j,k)), kind=mytype)
   eee(j,ny/2-1,k)=cmplx(real(eee(j,ny/2-1,k))*real(b1(j,k))-real(a1(j,k))*real(eee(j,ny/2,k)),&
        aimag(eee(j,ny/2-1,k))*aimag(b1(j,k))-aimag(a1(j,k))*aimag(eee(j,ny/2,k)), kind=mytype)
enddo
enddo

do i=ny/2-2,1,-1
do k=spI%yst(3),spI%yen(3)
do j=spI%yst(1),spI%yen(1)
   if (abs(real(aaa(j,i,k,3), kind=mytype)).gt.epsilon) then
      tmp1=1./real(aaa(j,i,k,3), kind=mytype)
   else
      tmp1=0.
   endif
   if (abs(aimag(aaa(j,i,k,3))).gt.epsilon) then
      tmp2=1./aimag(aaa(j,i,k,3))
   else
      tmp2=0.
   endif
   sr(j,k)=cmplx(tmp1,tmp2, kind=mytype)
   a1(j,k)=cmplx(real(aaa(j,i,k,4), kind=mytype)*real(sr(j,k), kind=mytype),&
        aimag(aaa(j,i,k,4))*aimag(sr(j,k)), kind=mytype)
   b1(j,k)=cmplx(real(aaa(j,i,k,5), kind=mytype)*real(sr(j,k), kind=mytype),&
        aimag(aaa(j,i,k,5))*aimag(sr(j,k)), kind=mytype)
   eee(j,i,k)=cmplx(real(eee(j,i,k), kind=mytype)*real(sr(j,k), kind=mytype)-&
        real(a1(j,k), kind=mytype)*real(eee(j,i+1,k), kind=mytype)-&
        real(b1(j,k), kind=mytype)*real(eee(j,i+2,k), kind=mytype),&
        aimag(eee(j,i,k))*aimag(sr(j,k))-&
        aimag(a1(j,k))*aimag(eee(j,i+1,k))-aimag(b1(j,k))*aimag(eee(j,i+2,k)), kind=mytype)
enddo
enddo
enddo

return

end subroutine inversion5_v1
