  subroutine poisson_100(rhs)

    implicit none

    real(mytype), dimension(:,:,:), intent(INOUT) :: rhs

    complex(mytype) :: xyzk
    real(mytype) :: tmp1, tmp2, tmp3, tmp4
    real(mytype) :: xx1,xx2,xx3,xx4,xx5,xx6,xx7,xx8

    integer :: nx,ny,nz, i,j,k, itmp

100 format(1x,a8,3I4,2F12.6)

    nx = nx_global - 1
    ny = ny_global
    nz = nz_global

    ! rhs is in Z-pencil but requires global operations in X
    call transpose_z_to_y(rhs,rw2,ph)
    call transpose_y_to_x(rw2,rw1,ph)
    do k=ph%xst(3),ph%xen(3)
       do j=ph%xst(2),ph%xen(2)
          do i=1,nx/2
             rw1b(i,j,k)=rw1(2*(i-1)+1,j,k)
          enddo
          do i=nx/2+1,nx
             rw1b(i,j,k)=rw1(2*nx-2*i+2,j,k)
          enddo
       enddo
    end do

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
    do k = sp%xst(3),sp%xen(3)
       do j = sp%xst(2),sp%xen(2)
          do i = sp%xst(1),sp%xen(1)
             tmp1 = real(cw1(i,j,k), kind=mytype)
             tmp2 = aimag(cw1(i,j,k))
             cw1(i,j,k) = cmplx(tmp1*by(j)+tmp2*ay(j), &
                  tmp2*by(j)-tmp1*ay(j), kind=mytype)
             if (j.gt.(ny/2+1)) cw1(i,j,k)=-cw1(i,j,k)
#ifdef DEBUG
             if (abs(cw1(i,j,k)) > 1.0e-4) &
                  write(*,100) 'after y',i,j,k,cw1(i,j,k)
#endif
          end do
       end do
    end do

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
#ifdef DEBUG
    do k = sp%xst(3),sp%xen(3)
       do j = sp%xst(2),sp%xen(2)
          do i = sp%xst(1),sp%xen(1)
             if (abs(cw1b(i,j,k)) > 1.0e-4) then
                write(*,100) 'after x',i,j,k,cw1b(i,j,k)
             end if
          end do
       end do
    end do
#endif

    ! Solve Poisson
    do k = sp%xst(3),sp%xen(3)
       do j = sp%xst(2),sp%xen(2)
          do i = sp%xst(1),sp%xen(1)
             !tmp1=real(zk2(k)+yk2(j)+xk2(i), kind=mytype)
             !tmp2=aimag(zk2(k)+yk2(j)+xk2(i))
             tmp1=real(kxyz(i,j,k), kind=mytype)
             tmp2=aimag(kxyz(i,j,k))
             !xyzk=cmplx(tmp1,tmp2, kind=mytype)
             ! CANNOT DO A DIVISION BY ZERO
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
#ifdef DEBUG
             if (abs(cw1b(i,j,k)) > 1.0e-4) &
                  write(*,100) 'AFTER',i,j,k,cw1b(i,j,k)
#endif
          end do
       end do
    end do

    ! post-processing backward

    ! POST PROCESSING IN X
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
    do k = sp%xst(3),sp%xen(3)
       do j = sp%xst(2),sp%xen(2)
          do i = sp%xst(1),sp%xen(1)
             tmp1 = real(cw1(i,j,k), kind=mytype)
             tmp2 = aimag(cw1(i,j,k))
             cw1(i,j,k) = cmplx(tmp1*by(j)-tmp2*ay(j), &
                  tmp2*by(j)+tmp1*ay(j), kind=mytype)
             if (j.gt.(ny/2+1)) cw1(i,j,k)=-cw1(i,j,k)
#ifdef DEBUG
             if (abs(cw1(i,j,k)) > 1.0e-4) &
                  write(*,100) 'AFTER Y',i,j,k,cw1(i,j,k)
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

    ! compute c2r transform
    call decomp_2d_fft_3d(cw1,rhs)

    ! rhs is in Z-pencil but requires global operations in X
    call transpose_z_to_y(rhs,rw2,ph)
    call transpose_y_to_x(rw2,rw1,ph)
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
  end subroutine poisson_100
