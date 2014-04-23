     Module Procedure test_eq_ctxt
  End Interface
  Interface operator( /= )
     Module Procedure test_ne_ctxt
  End Interface

  ! Interfaces for communicator comparison function
  Interface operator( == )
     Module Procedure test_eq_comm
  End Interface
  Interface operator( /= )
     Module Procedure test_ne_comm
  End Interface

  !!!!> \endcond

  ! Overloaded interfaces

!!$  Interface FIPC_allreduce
!!$     Module Procedure allreduce_double
!!$  End Interface

Contains

  !> \brief Initialise FreeIPC.
  !!
  !! This routine initializes FreeIPC. It sets up a context spanning all process in the MPI communic
  !! universe_comm, i.e. FIPC_ctxt_world, and determines which processes are on the smae shared memo
  !! \param universe_comm The communicator that spans all the processes within which FreeIPC will wo
  !! \param error         On success ERROR is set to FIPC_SUCCESS. Any other value indicates error.
  !>
  Subroutine FIPC_init( universe_comm, error )

    Integer, Intent( In    ) :: universe_comm
    Integer, Intent(   Out ) :: error

    Integer :: world_comm, intra_comm, extra_comm

    ! Can only initialise once
    If( Associated( FIPC_ctxt_world%initialized ) ) Then
       error = FIPC_already_initialized
       Return
    End If

    ! Get the error values that the system V routine can return
    Call FIPC_get_errval( EACCES, EEXIST, EINVAL, ENFILE, ENOENT, ENOMEM, ENOSPC )

    ! From the universe communicator generate the 3 communicators required
    ! to specify FIPC_ctxt_world
    Call generate_base_comms( universe_comm, world_comm, intra_comm, extra_comm, error )
    If( error /= FIPC_success ) Then
       Return
    End If

    ! From those communicators generate FIPC_ctxt_world
    Call ctxt_create( world_comm, intra_comm, extra_comm, FIPC_ctxt_world, error )
    If( error /= FIPC_success ) Then
       Return
    End If

    error = FIPC_success

  End Subroutine FIPC_init
