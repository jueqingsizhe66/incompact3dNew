module decomp_2d_poisson

  use decomp_2d
  use decomp_2d_fft

  use param
  use variables

  implicit none

  private        ! Make everything private unless declared public

!  real(mytype), private, parameter :: PI = 3.14159265358979323846_mytype

#ifdef DOUBLE_PREC
  real(mytype), parameter :: epsilon = 1.e-16
#else
  real(mytype), parameter :: epsilon = 1.e-8
#endif

  ! boundary conditions
  integer, save :: bcx, bcy, bcz

  ! decomposition object for physical space
  TYPE(DECOMP_INFO), save :: ph

  ! decomposition object for spectral space
  TYPE(DECOMP_INFO), save :: sp

  ! store sine/cosine factors
  real(mytype), save, allocatable, dimension(:) :: az,bz
  real(mytype), save, allocatable, dimension(:) :: ay,by
  real(mytype), save, allocatable, dimension(:) :: ax,bx

  ! wave numbers
  complex(mytype), save, allocatable, dimension(:,:,:) :: kxyz
  !wave numbers for stretching in a pentadiagonal matrice
  complex(mytype), save, allocatable, dimension(:,:,:,:) :: a,a2,a3
  ! work arrays,
  ! naming convention: cw (complex); rw (real);
  !                    1 = X-pencil; 2 = Y-pencil; 3 = Z-pencil
  real(mytype), allocatable, dimension(:,:,:) :: rw1,rw1b,rw2,rw2b,rw3
  complex(mytype), allocatable, dimension(:,:,:) :: cw1,cw1b,cw2,cw22,cw2b,cw2c

  ! underlying FFT library only needs to be initialised once
  logical, save :: fft_initialised = .false.

  public :: decomp_2d_poisson_stg, decomp_2d_poisson_init, &
       decomp_2d_poisson_finalize

  ! For staggered mesh where main variables are defined in the centre of
  ! control volumes while boundary conditions are defined on interfaces
  interface decomp_2d_poisson_stg
     module procedure poisson
  end interface
contains



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Initialise Poisson solver for given boundary conditions
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine decomp_2d_poisson_init(bcx1, bcy1, bcz1)

    implicit none

    integer, intent(IN) :: bcx1, bcy1, bcz1
    integer :: nx, ny, nz, i

    bcx = bcx1
    bcy = bcy1
    bcz = bcz1

    nx = nx_global
    ny = ny_global
    nz = nz_global

    ! pressure-grid having 1 fewer point for non-periodic directions
    if (bcx==1) nx=nx-1
    if (bcy==1) ny=ny-1
    if (bcz==1) nz=nz-1

    allocate(ax(nx),bx(nx))
    allocate(ay(ny),by(ny))
    allocate(az(nz),bz(nz))
    call abxyz(ax,ay,az,bx,by,bz,nx,ny,nz,bcx,bcy,bcz)

    call decomp_info_init(nx, ny, nz, ph)
    call decomp_info_init(nx, ny, nz/2+1, sp)

    ! allocate work space
    if (bcx==0 .and. bcy==0 .and. bcz==0) then
       allocate(cw1(sp%xst(1):sp%xen(1),sp%xst(2):sp%xen(2), &
            sp%xst(3):sp%xen(3)))
       allocate(kxyz(sp%xst(1):sp%xen(1),sp%xst(2):sp%xen(2), &
            sp%xst(3):sp%xen(3)))
       allocate(a(sp%yst(1):sp%yen(1),ny/2,sp%yst(3):sp%yen(3),5))
       allocate(a2(sp%yst(1):sp%yen(1),ny/2,sp%yst(3):sp%yen(3),5))
       allocate(a3(sp%yst(1):sp%yen(1),ny,sp%yst(3):sp%yen(3),5))
    else if (bcx==1 .and. bcy==0 .and. bcz==0) then
       allocate(cw1(sp%xst(1):sp%xen(1),sp%xst(2):sp%xen(2), &
                 sp%xst(3):sp%xen(3)))
       allocate(cw1b(sp%xst(1):sp%xen(1),sp%xst(2):sp%xen(2), &
                 sp%xst(3):sp%xen(3)))
       allocate(rw1(ph%xst(1):ph%xen(1),ph%xst(2):ph%xen(2), &
            ph%xst(3):ph%xen(3)))
       allocate(rw1b(ph%xst(1):ph%xen(1),ph%xst(2):ph%xen(2), &
            ph%xst(3):ph%xen(3)))
       allocate(rw2(ph%yst(1):ph%yen(1),ph%yst(2):ph%yen(2), &
            ph%yst(3):ph%yen(3)))
       allocate(kxyz(sp%xst(1):sp%xen(1),sp%xst(2):sp%xen(2), &
            sp%xst(3):sp%xen(3)))
       allocate(a(sp%yst(1):sp%yen(1),ny/2,sp%yst(3):sp%yen(3),5))
       allocate(a2(sp%yst(1):sp%yen(1),ny/2,sp%yst(3):sp%yen(3),5))
       allocate(a3(sp%yst(1):sp%yen(1),ny,sp%yst(3):sp%yen(3),5))
    else if (bcx==0 .and. bcy==1 .and. bcz==0) then
       allocate(rw2(ph%yst(1):ph%yen(1),ph%yst(2):ph%yen(2), &
            ph%yst(3):ph%yen(3)))
       allocate(rw2b(ph%yst(1):ph%yen(1),ph%yst(2):ph%yen(2), &
            ph%yst(3):ph%yen(3)))
       allocate(cw1(sp%xst(1):sp%xen(1),sp%xst(2):sp%xen(2), &
                 sp%xst(3):sp%xen(3)))
       allocate(cw2(sp%yst(1):sp%yen(1),sp%yst(2):sp%yen(2), &
                 sp%yst(3):sp%yen(3)))
       allocate(cw22(sp%yst(1):sp%yen(1),sp%yst(2):sp%yen(2), &
            sp%yst(3):sp%yen(3)))
       allocate(cw2b(sp%yst(1):sp%yen(1),sp%yst(2):sp%yen(2), &
                 sp%yst(3):sp%yen(3)))
       allocate(cw2c(sp%yst(1):sp%yen(1),sp%yst(2):sp%yen(2), &
                 sp%yst(3):sp%yen(3)))
       allocate(kxyz(sp%yst(1):sp%yen(1),sp%yst(2):sp%yen(2), &
            sp%yst(3):sp%yen(3)))
       allocate(a(sp%yst(1):sp%yen(1),ny/2,sp%yst(3):sp%yen(3),5))
       allocate(a2(sp%yst(1):sp%yen(1),ny/2,sp%yst(3):sp%yen(3),5))
       allocate(a3(sp%yst(1):sp%yen(1),ny,sp%yst(3):sp%yen(3),5))
    else if (bcx==1 .and. bcy==1) then
       allocate(cw1(sp%xst(1):sp%xen(1),sp%xst(2):sp%xen(2), &
                 sp%xst(3):sp%xen(3)))
       allocate(cw1b(sp%xst(1):sp%xen(1),sp%xst(2):sp%xen(2), &
                 sp%xst(3):sp%xen(3)))
       allocate(cw2(sp%yst(1):sp%yen(1),sp%yst(2):sp%yen(2), &
                 sp%yst(3):sp%yen(3)))
       allocate(cw22(sp%yst(1):sp%yen(1),sp%yst(2):sp%yen(2), &
            sp%yst(3):sp%yen(3)))
       allocate(cw2b(sp%yst(1):sp%yen(1),sp%yst(2):sp%yen(2), &
                 sp%yst(3):sp%yen(3)))
       allocate(cw2c(sp%yst(1):sp%yen(1),sp%yst(2):sp%yen(2), &
                 sp%yst(3):sp%yen(3)))
       allocate(rw1(ph%xst(1):ph%xen(1),ph%xst(2):ph%xen(2), &
            ph%xst(3):ph%xen(3)))
       allocate(rw1b(ph%xst(1):ph%xen(1),ph%xst(2):ph%xen(2), &
            ph%xst(3):ph%xen(3)))
       allocate(rw2(ph%yst(1):ph%yen(1),ph%yst(2):ph%yen(2), &
            ph%yst(3):ph%yen(3)))
       allocate(rw2b(ph%yst(1):ph%yen(1),ph%yst(2):ph%yen(2), &
            ph%yst(3):ph%yen(3)))
       if (bcz==1) then
          allocate(rw3(ph%zsz(1),ph%zsz(2),ph%zsz(3)))
       end if
       allocate(kxyz(sp%xst(1):sp%xen(1),sp%xst(2):sp%xen(2), &
            sp%xst(3):sp%xen(3)))
       allocate(a(sp%yst(1):sp%yen(1),ny/2,sp%yst(3):sp%yen(3),5))
       allocate(a2(sp%yst(1):sp%yen(1),ny/2,sp%yst(3):sp%yen(3),5))
       allocate(a3(sp%yst(1):sp%yen(1),ny,sp%yst(3):sp%yen(3),5))
    end if

    call waves()

    return
  end subroutine decomp_2d_poisson_init
