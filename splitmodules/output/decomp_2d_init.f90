  subroutine decomp_2d_init(nx,ny,nz,p_row,p_col)

    implicit none

    integer, intent(IN) :: nx,ny,nz,p_row,p_col

    integer :: errorcode, ierror

#ifdef SHM_DEBUG
    character(len=80) fname
#endif

    call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierror)
    if (nproc /= p_row*p_col) then
       errorcode = 1
       call MPI_ABORT(MPI_COMM_WORLD,errorcode,ierror)
    endif
    call MPI_COMM_RANK(MPI_COMM_WORLD,nrank,ierror)

    nx_global = nx
    ny_global = ny
    nz_global = nz

    ! 2D Catersian topology
    dims(1) = p_row
    dims(2) = p_col
    periodic(1) = .false.
    periodic(2) = .false.
    call MPI_CART_CREATE(MPI_COMM_WORLD,2,dims,periodic, &
         .false., &  ! do not reorder rank
         MPI_COMM_CART, ierror)
    call MPI_CART_COORDS(MPI_COMM_CART,nrank,2,coord,ierror)

    ! these are the communicators defining sub-groups for ALLTOALL(V)
    call MPI_CART_SUB(MPI_COMM_CART,(/.true.,.false./), &
         MPI_COMM_COL,ierror)
    call MPI_CART_SUB(MPI_COMM_CART,(/.false.,.true./), &
         MPI_COMM_ROW,ierror)

    ! actually generate all 2D decomposition information
    call decomp_info_init(nx,ny,nz,main)

    ! make a copy in global variables so applications can use
    ! these to create data structures
    xstart = main%xst
    ystart = main%yst
    zstart = main%zst
    xend   = main%xen
    yend   = main%yen
    zend   = main%zen
    xsize  = main%xsz
    ysize  = main%ysz
    zsize  = main%zsz

#ifdef SHM_DEBUG
    write(fname,99) nrank
99  format('log',I2.2)
    open(10,file=fname)
    write(10,*)'I am mpi rank ', nrank, 'Total ranks ', nproc
    write(10,*)' '
    write(10,*)'Global data size:'
    write(10,*)'nx*ny*nz', nx,ny,nz
    write(10,*)' '
    write(10,*)'2D processor grid:'
    write(10,*)'p_row*p_col:', p_row, p_col
    write(10,*)' '
    write(10,*)'Portion of global data held locally:'
    write(10,*)'xsize:',xsize
    write(10,*)'ysize:',ysize
    write(10,*)'zsize:',zsize
    write(10,*)' '
    write(10,*)'How pensils are to be divided and sent in alltoallv:'
    write(10,*)'x1dist:',main%x1dist
    write(10,*)'y1dist:',main%y1dist
    write(10,*)'y2dist:',main%y2dist
    write(10,*)'z2dist:',main%z2dist
    write(10,*)' '
    write(10,*)'######Shared buffer set up after this point######'
    write(10,*)' '
    write(10,*) 'col communicator detais:'
    call print_smp_info(main%COL_INFO)
    write(10,*)' '
    write(10,*) 'row communicator detais:'
    call print_smp_info(main%ROW_INFO)
    write(10,*)' '
    write(10,*)'Buffer count and dispalcement of per-core buffers'
    write(10,*)'x1cnts:',main%x1cnts
    write(10,*)'y1cnts:',main%y1cnts
    write(10,*)'y2cnts:',main%y2cnts
    write(10,*)'z2cnts:',main%z2cnts
    write(10,*)'x1disp:',main%x1disp
    write(10,*)'y1disp:',main%y1disp
    write(10,*)'y2disp:',main%y2disp
    write(10,*)'z2disp:',main%z2disp
    write(10,*)' '
    write(10,*)'Buffer count and dispalcement of shared buffers'
    write(10,*)'x1cnts:',main%x1cnts_s
    write(10,*)'y1cnts:',main%y1cnts_s
    write(10,*)'y2cnts:',main%y2cnts_s
    write(10,*)'z2cnts:',main%z2cnts_s
    write(10,*)'x1disp:',main%x1disp_s
    write(10,*)'y1disp:',main%y1disp_s
    write(10,*)'y2disp:',main%y2disp_s
    write(10,*)'z2disp:',main%z2disp_s
    write(10,*)' '
    close(10)
#endif

    return
  end subroutine decomp_2d_init
