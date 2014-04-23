  subroutine get_smp_map(comm, nnodes, my_node, ncores, my_core, maxcor)

    use FIPC_module

    implicit none

    integer, intent(IN) :: comm
    integer, intent(OUT) :: nnodes, my_node, ncores, my_core, maxcor

    integer :: intra_comm, extra_comm
    integer :: ierror

    call FIPC_init(comm, ierror)

    ! intra_comm: communicator for processes on this shared memory node
    ! extra_comm: communicator for all rank 0 on each shared memory node
    call FIPC_ctxt_intra_comm(FIPC_ctxt_world, intra_comm, ierror)
    call FIPC_ctxt_extra_comm(FIPC_ctxt_world, extra_comm, ierror)

    call MPI_COMM_SIZE(intra_comm,  ncores, ierror)
    call MPI_COMM_RANK(intra_comm, my_core, ierror)

    ! only rank 0 on each shared memory node member of extra_comm
    ! for others extra_comm = MPI_COMM_NULL
    if (extra_comm /= MPI_COMM_NULL) then
       call MPI_COMM_SIZE(extra_comm,  nnodes, ierror)
       call MPI_COMM_RANK(extra_comm, my_node, ierror)
    end if

    ! other ranks share the same information as their leaders
    call MPI_BCAST( nnodes, 1, MPI_INTEGER, 0, intra_comm, ierror)
    call MPI_BCAST(my_node, 1, MPI_INTEGER, 0, intra_comm, ierror)

    ! maxcor
    call MPI_ALLREDUCE(ncores, maxcor, 1, MPI_INTEGER, MPI_MAX, &
         MPI_COMM_WORLD, ierror)

    call FIPC_finalize(ierror)

    return

  end subroutine get_smp_map
