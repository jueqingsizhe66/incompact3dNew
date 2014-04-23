  subroutine init_neighbour

    integer :: ierror

    ! For X-pencil
    neighbour(1,1) = MPI_PROC_NULL               ! east
    neighbour(1,2) = MPI_PROC_NULL               ! west
    call MPI_CART_SHIFT(DECOMP_2D_COMM_CART_X, 0, 1, &
         neighbour(1,4), neighbour(1,3), ierror) ! north & south
    call MPI_CART_SHIFT(DECOMP_2D_COMM_CART_X, 1, 1, &
         neighbour(1,6), neighbour(1,5), ierror) ! top & bottom

    ! For Y-pencil
    call MPI_CART_SHIFT(DECOMP_2D_COMM_CART_Y, 0, 1, &
         neighbour(2,2), neighbour(2,1), ierror) ! east & west
    neighbour(2,3) = MPI_PROC_NULL               ! north
    neighbour(2,4) = MPI_PROC_NULL               ! south
    call MPI_CART_SHIFT(DECOMP_2D_COMM_CART_Y, 1, 1, &
         neighbour(2,6), neighbour(2,5), ierror) ! top & bottom

    ! For Z-pencil
    call MPI_CART_SHIFT(DECOMP_2D_COMM_CART_Z, 0, 1, &
         neighbour(3,2), neighbour(3,1), ierror) ! east & west
    call MPI_CART_SHIFT(DECOMP_2D_COMM_CART_Z, 1, 1, &
         neighbour(3,4), neighbour(3,3), ierror) ! north & south
    neighbour(3,5) = MPI_PROC_NULL               ! top
    neighbour(3,6) = MPI_PROC_NULL               ! bottom

    return
  end subroutine init_neighbour
