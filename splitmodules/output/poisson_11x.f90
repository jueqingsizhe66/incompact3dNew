  subroutine poisson_11x(rhs, nclz1)

    implicit none

    integer, intent(IN) :: nclz1	
    real(mytype), dimension(:,:,:), intent(INOUT) :: rhs

    complex(mytype) :: xyzk
    real(mytype) :: tmp1, tmp2, tmp3, tmp4
    real(mytype) :: xx1,xx2,xx3,xx4,xx5,xx6,xx7,xx8

    integer :: nx,ny,nz, i,j,k

100 format(1x,a8,3I4,2F12.6)

    nx = nx_global - 1
    ny = ny_global - 1

    if (nclz1==1) then	
       nz = nz_global - 1
    else if (nclz1==0) then
       nz = nz_global
    end if

    if (nclz1==1) then
       do j=1,ph%zsz(2)
          do i=1,ph%zsz(1)
             do k=1,nz/2
                rw3(i,j,k)=rhs(i,j,2*(k-1)+1)
             end do
             do k=nz/2+1,nz
                rw3(i,j,k)=rhs(i,j,2*nz-2*k+2)
             end do
          end do
       end do
       call transpose_z_to_y(rw3,rw2,ph)
    else if (nclz1==0) then
       call transpose_z_to_y(rhs,rw2,ph)
    end if


    do k=ph%yst(3),ph%yen(3)
       do i=ph%yst(1),ph%yen(1)
          do j=1,ny/2
             rw2b(i,j,k)=rw2(i,2*(j-1)+1,k)
          end do
          do j=ny/2+1,ny
             rw2b(i,j,k)=rw2(i,2*ny-2*j+2,k)
          end do
       end do
    end do

    ! the global operations in X
    call transpose_y_to_x(rw2b,rw1,ph)

    do k=ph%xst(3),ph%xen(3)
       do j=ph%xst(2),ph%xen(2)
          do i=1,nx/2
             rw1b(i,j,k)=rw1(2*(i-1)+1,j,k)
          end do
          do i=nx/2+1,nx
             rw1b(i,j,k)=rw1(2*nx-2*i+2,j,k)
          end do
       end do
    end do

    ! back to Z-pencil
    call transpose_x_to_y(rw1b,rw2,ph)
    call transpose_y_to_z(rw2,rhs,ph)

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

    ! POST PROCESSING IN Y
    ! WE HAVE TO BE IN Y PENCILS
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

    ! back to X-pencil
    call transpose_y_to_x(cw2b,cw1,sp)
#ifdef DEBUG
    do k = sp%xst(3),sp%xen(3)
       do j = sp%xst(2),sp%xen(2)
          do i = sp%xst(1),sp%xen(1)
             if (abs(cw1(i,j,k)) > 1.0e-4) then
                write(*,100) 'after y',i,j,k,cw1(i,j,k)
             end if
          end do
       end do
    end do
#endif

    ! POST PROCESSING IN X
    do k = sp%xst(3),sp%xen(3)
       do j = sp%xst(2),sp%xen(2)
          cw1b(1,j,k)=cw1(1,j,k)
          do i = 2,nx
             tmp1 = real(cw1(i,j,k), kind=mytype)
             tmp2 = aimag(cw1(i,j,k))
             tmp3 = real(cw1(nx-i+2,j,k), kind=mytype)
             tmp4 = aimag(cw1(nx-i+2,j,k))
             xx1=tmp1*bx(i)/2._mytype
             xx2=tmp1*ax(i)/2._mytype
             xx3=tmp2*bx(i)/2._mytype
             xx4=tmp2*ax(i)/2._mytype
             xx5=tmp3*bx(i)/2._mytype
             xx6=tmp3*ax(i)/2._mytype
             xx7=tmp4*bx(i)/2._mytype
             xx8=tmp4*ax(i)/2._mytype
             cw1b(i,j,k) = cmplx(xx1+xx4+xx5-xx8,-xx2+xx3+xx6+xx7, &
                  kind=mytype)
          end do
       end do
    end do

!#ifdef DEBUG
    do k = sp%xst(3),sp%xen(3)
       do j = sp%xst(2),sp%xen(2)
          do i = sp%xst(1),sp%xen(1)
             if (abs(cw1b(i,j,k)) > 1.0e-4) then
                write(*,*) 'BEFORE',i,j,k,cw1b(i,j,k)
             end if
          end do
       end do
    end do
!#endif

    if (istret==0) then

    ! Solve Poisson
    do k = sp%xst(3),sp%xen(3)
       do j = sp%xst(2),sp%xen(2)
          do i = sp%xst(1),sp%xen(1)
             !tmp1=real(zk2(k)+yk2(j)+xk2(i), kind=mytype)
             !tmp2=aimag(zk2(k)+yk2(j)+xk2(i))
             tmp1=real(kxyz(i,j,k), kind=mytype)
             tmp2=aimag(kxyz(i,j,k))
             !xyzk=cmplx(tmp1,tmp2, kind=mytype)
             !CANNOT DO A DIVISION BY ZERO
             if ((abs(tmp1).lt.epsilon).and.(abs(tmp2).lt.epsilon)) then
                cw1b(i,j,k)=cmplx(0._mytype,0._mytype, kind=mytype)
             end if
             if ((abs(tmp1).lt.epsilon).and.(abs(tmp2).ge.epsilon)) then
                cw1b(i,j,k)=cmplx(0._mytype, &
                     aimag(cw1b(i,j,k))/(-tmp2), kind=mytype)
             end if
             if ((abs(tmp1).ge.epsilon).and.(abs(tmp2).lt.epsilon)) then
                cw1b(i,j,k)=cmplx( real(cw1b(i,j,k), kind=mytype) &
                     /(-tmp1), 0._mytype, kind=mytype)
             end if
             if ((abs(tmp1).ge.epsilon).and.(abs(tmp2).ge.epsilon)) then
                cw1b(i,j,k)=cmplx( real(cw1b(i,j,k), kind=mytype) &
                     /(-tmp1), &
                     aimag(cw1b(i,j,k))/(-tmp2), kind=mytype)
             end if
          end do
       end do
    end do

    else
       call matrice_refinement()
! the stretching is only working in Y pencils
       call transpose_x_to_y(cw1b,cw2b,sp)
       !we are now in Y pencil

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



          call inversion5_v1(a,cw2,sp)
          call inversion5_v1(a2,cw2c,sp)

! do k = sp%yst(3), sp%yen(3)
!       do j = 1,ny/2
!          do i = sp%yst(1), sp%yen(1)
!             if (abs(cw2(i,j,k)) > 1.0e-4) then
!                write(*,*) 'before IN',i,j,k,cw2(i,j,k)!*2.
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
!we have to go back in X pencils
       call transpose_y_to_x(cw2b,cw1b,sp)
    endif

!#ifdef DEBUG
    do k = sp%xst(3),sp%xen(3)
       do j = sp%xst(2),sp%xen(2)
          do i = sp%xst(1),sp%xen(1)
             if (abs(cw1b(i,j,k)) > 1.0e-4) then
                write(*,*) 'AFTER',i,j,k,cw1b(i,j,k)
             end if
          end do
       end do
    end do
!#endif
stop
    ! post-processing backward

    do k = sp%xst(3),sp%xen(3)
       do j = sp%xst(2),sp%xen(2)
          cw1(1,j,k)=cw1b(1,j,k)
          do i = 2,nx
             tmp1 = real(cw1b(i,j,k), kind=mytype)
             tmp2 = aimag(cw1b(i,j,k))
             tmp3 = real(cw1b(nx-i+2,j,k), kind=mytype)
             tmp4 = aimag(cw1b(nx-i+2,j,k))
             xx1=tmp1*bx(i)
             xx2=tmp1*ax(i)
             xx3=tmp2*bx(i)
             xx4=tmp2*ax(i)
             xx5=tmp3*bx(i)
             xx6=tmp3*ax(i)
             xx7=tmp4*bx(i)
             xx8=tmp4*ax(i)
             cw1(i,j,k) = cmplx(xx1-xx4+xx6+xx7,-(-xx2-xx3+xx5-xx8), &
                  kind=mytype)
          end do
       end do
    end do
#ifdef DEBUG
    do k = sp%xst(3),sp%xen(3)
       do j = sp%xst(2),sp%xen(2)
          do i = sp%xst(1),sp%xen(1)
             if (abs(cw1(i,j,k)) > 1.0e-4) then
                write(*,100) 'AFTER X',i,j,k,cw1(i,j,k)
             end if
          end do
       end do
    end do
#endif

    ! POST PROCESSING IN Y
    ! NEED to be in Y-pencil
    call transpose_x_to_y(cw1,cw2,sp)
    do k = sp%yst(3), sp%yen(3)
       do i = sp%yst(1), sp%yen(1)
          cw2b(i,1,k)=cw2(i,1,k)
          do j = 2,ny
             tmp1 = real(cw2(i,j,k), kind=mytype)
             tmp2 = aimag(cw2(i,j,k))
             tmp3 = real(cw2(i,ny-j+2,k), kind=mytype)
             tmp4 = aimag(cw2(i,ny-j+2,k))
             xx1=tmp1*by(j)
             xx2=tmp1*ay(j)
             xx3=tmp2*by(j)
             xx4=tmp2*ay(j)
             xx5=tmp3*by(j)
             xx6=tmp3*ay(j)
             xx7=tmp4*by(j)
             xx8=tmp4*ay(j)
             cw2b(i,j,k) = cmplx(xx1-xx4+xx6+xx7,-(-xx2-xx3+xx5-xx8), &
                  kind=mytype)
          end do
       end do
    end do
#ifdef DEBUG
    do k = sp%yst(3), sp%yen(3)
       do j = sp%yst(2), sp%yen(2)
          do i = sp%yst(1), sp%yen(1)
             if (abs(cw2b(i,j,k)) > 1.0e-4) then
                write(*,100) 'AFTER Y',i,j,k,cw2b(i,j,k)
             end if
          end do
       end do
    end do
#endif
    ! back to X-pencil
    call transpose_y_to_x(cw2b,cw1,sp)

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

    if (nclz1==1) then
       do j=1,ph%zsz(2)
          do i=1,ph%zsz(1)
             do k=1,nz/2
                rw3(i,j,2*k-1)=rhs(i,j,k)
             end do
             do k=1,nz/2
                rw3(i,j,2*k)=rhs(i,j,nz-k+1)
             end do
          end do
       end do
       call transpose_z_to_y(rw3,rw2,ph)
    else if (nclz1==0) then
       call transpose_z_to_y(rhs,rw2,ph)
    end if

    do k=ph%yst(3),ph%yen(3)
       do i=ph%yst(1),ph%yen(1)
          do j=1,ny/2
             rw2b(i,2*j-1,k)=rw2(i,j,k)
          end do
          do j=1,ny/2
             rw2b(i,2*j,k)=rw2(i,ny-j+1,k)
          end do
       enddo
    end do
    call transpose_y_to_x(rw2b,rw1,ph)
    do k=ph%xst(3),ph%xen(3)
       do j=ph%xst(2),ph%xen(2)
          do i=1,nx/2
             rw1b(2*i-1,j,k)=rw1(i,j,k)
          enddo
          do i=1,nx/2
             rw1b(2*i,j,k)=rw1(nx-i+1,j,k)
          enddo
       enddo
    end do
    call transpose_x_to_y(rw1b,rw2,ph)
    call transpose_y_to_z(rw2,rhs,ph)

  !  call decomp_2d_fft_finalize



    return
  end subroutine poisson_11x
