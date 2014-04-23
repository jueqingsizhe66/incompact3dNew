  subroutine poisson_010(rhs)

    implicit none

    real(mytype), dimension(:,:,:), intent(INOUT) :: rhs

    complex(mytype) :: xyzk
    real(mytype) :: tmp1, tmp2, tmp3, tmp4
    real(mytype) :: xx1,xx2,xx3,xx4,xx5,xx6,xx7,xx8

    integer :: nx,ny,nz, i,j,k

100 format(1x,a8,3I4,2F12.6)

    nx = nx_global
    ny = ny_global - 1
    nz = nz_global

    ! rhs is in Z-pencil but requires global operations in Y
    call transpose_z_to_y(rhs,rw2,ph)
    do k=ph%yst(3),ph%yen(3)
       do i=ph%yst(1),ph%yen(1)
          do j=1,ny/2
             rw2b(i,j,k)=rw2(i,2*(j-1)+1,k)
          enddo
          do j=ny/2+1,ny
             rw2b(i,j,k)=rw2(i,2*ny-2*j+2,k)
          enddo
       enddo
    end do
    call transpose_y_to_z(rw2b,rhs,ph)

    if (.not. fft_initialised) then
       call decomp_2d_fft_init(PHYSICAL_IN_Z,nx,ny,nz)
       fft_initialised = .true.
    end if
    ! compute r2c transform
    call decomp_2d_fft_3d(rhs,cw1)

    ! normalisation
    cw1 = cw1 / real(nx, kind=mytype) /real(ny, kind=mytype) &
         / real(nz, kind=mytype)
#ifdef DEBUG
    do k = sp%xst(3),sp%xen(3)
       do j = sp%xst(2),sp%xen(2)
          do i = sp%xst(1),sp%xen(1)
             if (abs(cw1(i,j,k)) > 1.0e-4) then
                write(*,100) 'START',i,j,k,cw1(i,j,k)
             end if
          end do
       end do
    end do
#endif

    ! post-processing in spectral space

    ! POST PROCESSING IN Z
    do k = sp%xst(3),sp%xen(3)
       do j = sp%xst(2),sp%xen(2)
          do i = sp%xst(1),sp%xen(1)
             tmp1 = real(cw1(i,j,k), kind=mytype)
             tmp2 = aimag(cw1(i,j,k))
             cw1(i,j,k) = cmplx(tmp1*bz(k)+tmp2*az(k), &
                  tmp2*bz(k)-tmp1*az(k), kind=mytype)
#ifdef DEBUG
             if (abs(cw1(i,j,k)) > 1.0e-4) &
                  write(*,100) 'after z',i,j,k,cw1(i,j,k)
#endif
          end do
       end do
    end do

    ! POST PROCESSING IN X
    do k = sp%xst(3),sp%xen(3)
       do j = sp%xst(2),sp%xen(2)
          do i = sp%xst(1),sp%xen(1)
             tmp1 = real(cw1(i,j,k), kind=mytype)
             tmp2 = aimag(cw1(i,j,k))
             cw1(i,j,k) = cmplx(tmp1*bx(i)+tmp2*ax(i), &
                  tmp2*bx(i)-tmp1*ax(i), kind=mytype)
             if (i.gt.(nx/2+1)) cw1(i,j,k)=-cw1(i,j,k)
#ifdef DEBUG
             if (abs(cw1(i,j,k)) > 1.0e-4) &
                  write(*,100) 'after x',i,j,k,cw1(i,j,k)
#endif
          end do
       end do
    end do

    ! POST PROCESSING IN Y
    ! NEED TO BE IN Y PENCILS!!!!!!!!!!!!!!!
    call transpose_x_to_y(cw1,cw2,sp)

    do k = sp%yst(3), sp%yen(3)
       do i = sp%yst(1), sp%yen(1)
          cw2b(i,1,k)=cw2(i,1,k)
          do j = 2,ny
             tmp1 = real(cw2(i,j,k), kind=mytype)
             tmp2 = aimag(cw2(i,j,k))
             tmp3 = real(cw2(i,ny-j+2,k), kind=mytype)
             tmp4 = aimag(cw2(i,ny-j+2,k))
             xx1=tmp1*by(j)/2._mytype
             xx2=tmp1*ay(j)/2._mytype
             xx3=tmp2*by(j)/2._mytype
             xx4=tmp2*ay(j)/2._mytype
             xx5=tmp3*by(j)/2._mytype
             xx6=tmp3*ay(j)/2._mytype
             xx7=tmp4*by(j)/2._mytype
             xx8=tmp4*ay(j)/2._mytype
             cw2b(i,j,k) = cmplx(xx1+xx4+xx5-xx8,-xx2+xx3+xx6+xx7, &
                  kind=mytype)
          end do
       end do
    end do
#ifdef DEBUG
    do k = sp%yst(3), sp%yen(3)
       do j = sp%yst(2), sp%yen(2)
          do i = sp%yst(1), sp%yen(1)
             if (abs(cw2b(i,j,k)) > 1.0e-4) then
                write(*,100) 'after y',i,j,k,cw2b(i,j,k)
                print *,kxyz(i,j,k)
             end if
          end do
       end do
    end do
#endif

    if (istret==0) then

    ! Solve Poisson
    ! doing wave number division in Y-pencil
    do k = sp%yst(3), sp%yen(3)
       do j = sp%yst(2), sp%yen(2)
          do i = sp%yst(1), sp%yen(1)
             !tmp1=real(zk2(k)+yk2(j)+xk2(i), kind=mytype)
             !tmp2=aimag(zk2(k)+yk2(j)+xk2(i))
             tmp1=real(kxyz(i,j,k), kind=mytype)
             tmp2=aimag(kxyz(i,j,k))
             !xyzk=cmplx(tmp1,tmp2, kind=mytype)
             !CANNOT DO A DIVISION BY ZERO
             if ((abs(tmp1).lt.epsilon).and.(abs(tmp2).lt.epsilon)) then
                cw2b(i,j,k)=cmplx(0._mytype,0._mytype, kind=mytype)
             end if
             if ((abs(tmp1).lt.epsilon).and.(abs(tmp2).ge.epsilon)) then
                cw2b(i,j,k)=cmplx(0._mytype, &
                     aimag(cw2b(i,j,k))/(-tmp2), kind=mytype)
             end if
             if ((abs(tmp1).ge.epsilon).and.(abs(tmp2).lt.epsilon)) then
                cw2b(i,j,k)=cmplx( real(cw2b(i,j,k), kind=mytype) &
                     /(-tmp1), 0._mytype, kind=mytype)
             end if
             if ((abs(tmp1).ge.epsilon).and.(abs(tmp2).ge.epsilon)) then
                cw2b(i,j,k)=cmplx( real(cw2b(i,j,k), kind=mytype) &
                     /(-tmp1), &
                     aimag(cw2b(i,j,k))/(-tmp2), kind=mytype)
             end if
          end do
       end do
    end do

    else

       call matrice_refinement()
!       do k = sp%yst(3), sp%yen(3)
!          do j = 1,ny/2
!             do i = sp%yst(1), sp%yen(1)
!                print *,i,j,k,a(i,j,k,3)
!!                if (nrank.le.1) print *,i,j,k,a(i,j,k,3)
!!                if (nrank.gt.1) print *,i+4,j,k,a(i,j,k,3)
!             enddo
!          enddo
!       enddo


       if (istret.ne.3) then
          cw2(:,:,:)=0.;cw2c(:,:,:)=0.
          do k = sp%yst(3), sp%yen(3)
          do j = 1,ny/2
          do i = sp%yst(1), sp%yen(1)
             cw2(i,j,k)=cw2b(i,2*j-1,k)
             cw2c(i,j,k)=cw2b(i,2*j,k)
          enddo
          enddo
          enddo

 !   do k = sp%yst(3), sp%yen(3)
 !      do j = 1,ny/2
 !         do i = sp%yst(1), sp%yen(1)
 !            if (abs(cw2(i,j,k)) > 1.0e-4) then
 !               write(*,*) 'before IN',i,j,k,cw2(i,j,k)!*2.
 !!            end if
  !        end do
  !     end do
  !  end do

         call inversion5_v1(a,cw2,sp)
         call inversion5_v1(a2,cw2c,sp)

!         cw2(1,1,1)=cw2(1,1,1)*0.5


!   do k = sp%yst(3), sp%yen(3)
!       do j = 1,ny/2
!          do i = sp%yst(1), sp%yen(1)
!             if (abs(cw2c(i,j,k)) > 1.0e-4) then
!                write(*,*) 'after IN',i,j,k,cw2c(i,j,k)!*2.
!             end if
!          end do
!       end do
!    end do

          cw2b(:,:,:)=0.
          do k=sp%yst(3), sp%yen(3)
          do j=1,ny-1,2
          do i=sp%yst(1), sp%yen(1)
             cw2b(i,j,k)=cw2(i,(j+1)/2,k)
          enddo
          enddo
          do j=2,ny,2
          do i=sp%yst(1), sp%yen(1)
             cw2b(i,j,k)=cw2c(i,j/2,k)
          enddo
          enddo
          enddo
          !do k=sp%yst(3), sp%yen(3)
          !do i=sp%yst(1), sp%yen(1)
          !   if ((xkx(i)==0).and.(zkz(k)==0)) then
          !   !   cw2b(i,1,1)=0.
          !   !   cw2b(i,ny,1)=0.
          !   endif
          !enddo
          !enddo
       else
          do k = sp%yst(3), sp%yen(3)
          do j = 1,ny
          do i = sp%yst(1), sp%yen(1)
             cw2(i,j,k)=cw2b(i,j,k)
          enddo
          enddo
          enddo
          call inversion5_v2(a3,cw2,sp)
          do k = sp%yst(3), sp%yen(3)
          do j = 1,ny
          do i = sp%yst(1), sp%yen(1)
             cw2b(i,j,k)=cw2(i,j,k)
          enddo
          enddo
          enddo
       endif

    endif

!    print *,nrank, sp%yst(3),sp%yen(3),sp%yst(1),sp%yen(1)

!we are in Y pencil
    do k = sp%yst(3), sp%yen(3)
    do i = sp%yst(1), sp%yen(1)
       if ((i==nx/2+1).and.(k==nz/2+1)) then
          cw2b(i,:,k)=0.
       endif
    enddo
    enddo
#ifdef DEBUG
    do k = sp%yst(3), sp%yen(3)
       do j = sp%yst(2), sp%yen(2)
          do i = sp%yst(1), sp%yen(1)
             if (abs(cw2b(i,j,k)) > 1.0e-4) then
                write(*,100) 'AFTER',i,j,k,cw2b(i,j,k)
                print *,kxyz(i,j,k)
             end if
          end do
       end do
    end do
#endif

    ! post-processing backward

    ! POST PROCESSING IN Y
    do k = sp%yst(3), sp%yen(3)
       do i = sp%yst(1), sp%yen(1)
          cw2(i,1,k)=cw2b(i,1,k)
          do j = 2,ny
             tmp1 = real(cw2b(i,j,k), kind=mytype)
             tmp2 = aimag(cw2b(i,j,k))
             tmp3 = real(cw2b(i,ny-j+2,k), kind=mytype)
             tmp4 = aimag(cw2b(i,ny-j+2,k))
             xx1=tmp1*by(j)
             xx2=tmp1*ay(j)
             xx3=tmp2*by(j)
             xx4=tmp2*ay(j)
             xx5=tmp3*by(j)
             xx6=tmp3*ay(j)
             xx7=tmp4*by(j)
             xx8=tmp4*ay(j)
             cw2(i,j,k) = cmplx(xx1-xx4+xx6+xx7,-(-xx2-xx3+xx5-xx8), &
                  kind=mytype)
          end do
       end do
    end do

    ! Back to X-pencil
    call transpose_y_to_x(cw2,cw1,sp)
#ifdef DEBUG
    do k = sp%xst(3),sp%xen(3)
       do j = sp%xst(2),sp%xen(2)
          do i = sp%xst(1),sp%xen(1)
             if (abs(cw1(i,j,k)) > 1.0e-4) then
                write(*,100) 'AFTER Y',i,j,k,cw1(i,j,k)
             end if
          end do
       end do
    end do
#endif

    ! POST PROCESSING IN X
    do k = sp%xst(3),sp%xen(3)
       do j = sp%xst(2),sp%xen(2)
          do i = sp%xst(1),sp%xen(1)
             tmp1 = real(cw1(i,j,k), kind=mytype)
             tmp2 = aimag(cw1(i,j,k))
             cw1(i,j,k) = cmplx(tmp1*bx(i)-tmp2*ax(i), &
                  tmp2*bx(i)+tmp1*ax(i), kind=mytype)
             if (i.gt.(nx/2+1)) cw1(i,j,k)=-cw1(i,j,k)
#ifdef DEBUG
             if (abs(cw1(i,j,k)) > 1.0e-4) &
                  write(*,100) 'AFTER X',i,j,k,cw1(i,j,k)
#endif
          end do
       end do
    end do

    ! POST PROCESSING IN Z
    do k = sp%xst(3),sp%xen(3)
       do j = sp%xst(2),sp%xen(2)
          do i = sp%xst(1),sp%xen(1)
             tmp1 = real(cw1(i,j,k), kind=mytype)
             tmp2 = aimag(cw1(i,j,k))
             cw1(i,j,k) = cmplx(tmp1*bz(k)-tmp2*az(k), &
                  tmp2*bz(k)+tmp1*az(k), kind=mytype)
#ifdef DEBUG
             if (abs(cw1(i,j,k)) > 1.0e-4) &
                  write(*,100) 'END',i,j,k,cw1(i,j,k)
#endif
          end do
       end do
    end do

    ! compute c2r transform, back to physical space
    call decomp_2d_fft_3d(cw1,rhs)

    ! rhs is in Z-pencil but requires global operations in Y
    call transpose_z_to_y(rhs,rw2,ph)
    do k=ph%yst(3),ph%yen(3)
       do i=ph%yst(1),ph%yen(1)
          do j=1,ny/2
             rw2b(i,2*j-1,k)=rw2(i,j,k)
          enddo
          do j=1,ny/2
             rw2b(i,2*j,k)=rw2(i,ny-j+1,k)
          enddo
       enddo
    end do
    call transpose_y_to_z(rw2b,rhs,ph)

  !  call decomp_2d_fft_finalize

    return
  end subroutine poisson_010
