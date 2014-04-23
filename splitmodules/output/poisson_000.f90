  subroutine poisson_000(rhs)

    use derivX
    use derivY
    use derivZ

    ! right-hand-side of Poisson as input
    ! solution of Poisson as output
    real(mytype), dimension(:,:,:), intent(INOUT) :: rhs

    integer, dimension(3) :: fft_start, fft_end, fft_size

    complex(mytype) :: xyzk

    complex(mytype) :: ytt,xtt,ztt,yt1,xt1,yt2,xt2
    complex(mytype) :: xtt1,ytt1,ztt1,zt1,zt2


    real(mytype) :: tmp1, tmp2,x ,y, z

    integer :: nx,ny,nz, i,j,k

    nx = nx_global
    ny = ny_global
    nz = nz_global

    if (.not. fft_initialised) then
       call decomp_2d_fft_init(PHYSICAL_IN_Z)
       fft_initialised = .true.
    end if

    ! compute r2c transform
    call decomp_2d_fft_3d(rhs,cw1)

    ! normalisation
    cw1 = cw1 / real(nx, kind=mytype) /real(ny, kind=mytype) &
         / real(nz, kind=mytype)

    do k = sp%xst(3),sp%xen(3)
       do j = sp%xst(2),sp%xen(2)
          do i = sp%xst(1),sp%xen(1)

             ! post-processing in spectral space

             ! POST PROCESSING IN Z
             tmp1 = real(cw1(i,j,k), kind=mytype)
             tmp2 = aimag(cw1(i,j,k))
             cw1(i,j,k) = cmplx(tmp1*bz(k)+tmp2*az(k), &
                  tmp2*bz(k)-tmp1*az(k), kind=mytype)

             ! POST PROCESSING IN Y
             tmp1 = real(cw1(i,j,k), kind=mytype)
             tmp2 = aimag(cw1(i,j,k))
             cw1(i,j,k) = cmplx(tmp1*by(j)+tmp2*ay(j), &
                  tmp2*by(j)-tmp1*ay(j), kind=mytype)
             if (j.gt.(ny/2+1)) cw1(i,j,k)=-cw1(i,j,k)

             ! POST PROCESSING IN X
             tmp1 = real(cw1(i,j,k), kind=mytype)
             tmp2 = aimag(cw1(i,j,k))
             cw1(i,j,k) = cmplx(tmp1*bx(i)+tmp2*ax(i), &
                  tmp2*bx(i)-tmp1*ax(i), kind=mytype)
             if (i.gt.(nx/2+1)) cw1(i,j,k)=-cw1(i,j,k)

             ! Solve Poisson
             tmp1=real(kxyz(i,j,k), kind=mytype)
             tmp2=aimag(kxyz(i,j,k))
             ! CANNOT DO A DIVISION BY ZERO
             if ((tmp1.lt.epsilon).or.(tmp2.lt.epsilon)) then
                cw1(i,j,k)=0._mytype
!                print *,'DIV 0',i,j,k,epsilon
             else
                cw1(i,j,k)=cmplx( real(cw1(i,j,k), kind=mytype) / (-tmp1), &
                     aimag(cw1(i,j,k))/(-tmp2), kind=mytype)
             end if

           !Print result in spectal space after Poisson
      !     if (abs(out(i,j,k)) > 1.0e-4) then
      !        write(*,*) 'AFTER',i,j,k,out(i,j,k),xyzk
      !     end if

             ! post-processing backward

             ! POST PROCESSING IN Z
             tmp1 = real(cw1(i,j,k), kind=mytype)
             tmp2 = aimag(cw1(i,j,k))
             cw1(i,j,k) = cmplx(tmp1*bz(k)-tmp2*az(k), &
                  -tmp2*bz(k)-tmp1*az(k), kind=mytype)

             ! POST PROCESSING IN Y
             tmp1 = real(cw1(i,j,k), kind=mytype)
             tmp2 = aimag(cw1(i,j,k))
             cw1(i,j,k) = cmplx(tmp1*by(j)+tmp2*ay(j), &
                  tmp2*by(j)-tmp1*ay(j), kind=mytype)
             if (j.gt.(ny/2+1)) cw1(i,j,k)=-cw1(i,j,k)

             ! POST PROCESSING IN X
             tmp1 = real(cw1(i,j,k), kind=mytype)
             tmp2 = aimag(cw1(i,j,k))
             cw1(i,j,k) = cmplx(tmp1*bx(i)+tmp2*ax(i), &
                  -tmp2*bx(i)+tmp1*ax(i), kind=mytype)
             if (i.gt.(nx/2+1)) cw1(i,j,k)=-cw1(i,j,k)

          end do
       end do
    end do

    ! compute c2r transform
    call decomp_2d_fft_3d(cw1,rhs)

 !   call decomp_2d_fft_finalize

    return
  end subroutine poisson_000
