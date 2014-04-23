  Subroutine ctxt_create( world_comm, intra_comm, extra_comm, ctxt, error )

    ! From the basic data create a context.
    !
    ! On success ERROR is set to FIPC_success and ctxt holds the new context
    ! On error error is any other value than FIPC_success

    Integer          , Intent( In    ) :: world_comm
    Integer          , Intent( In    ) :: intra_comm
    Integer          , Intent( In    ) :: extra_comm
    Type( FIPC_ctxt ), Intent(   Out ) :: ctxt
    Integer          , Intent(   Out ) :: error

    Integer( c_int ) :: semid

    Allocate( ctxt%world_comm, Stat = error )
    If( error /= 0 ) Then
       error = FIPC_allocation_failed
       Return
    End If
    Call set_comm( world_comm, ctxt%world_comm, error )
    If( error /= 0 ) Then
       Return
    End If

    Allocate( ctxt%intra_comm, Stat = error )
    If( error /= 0 ) Then
       error = FIPC_allocation_failed
       Return
    End If
    Call set_comm( intra_comm, ctxt%intra_comm, error )
    If( error /= 0 ) Then
       Return
    End If

    ! The extra communicator is only defined on processor zero
    ! of the intra processors. However it is useful to know the
    ! extra rank and size of the zero proc in the intra comm on all procs
    ! in the inter comm - this is the "node" number for this multicore
    ! node in the universe, and the number of nodes. Hence allocate the extra comm on all procs
    Allocate( ctxt%extra_comm, Stat = error )
    If( error /= 0 ) Then
       error = FIPC_allocation_failed
       Return
    End If
    If( ctxt%intra_comm%rank == 0 ) Then
       Call set_comm( extra_comm, ctxt%extra_comm, error )
       If( error /= 0 ) Then
          Return
       End If
    End If
    Call mpi_bcast( ctxt%extra_comm%rank, 1, mpi_integer, 0, ctxt%intra_comm%handle, error )
    If( error /= 0 ) Then
       Return
    End If
    Call mpi_bcast( ctxt%extra_comm%size, 1, mpi_integer, 0, ctxt%intra_comm%handle, error )
    If( error /= 0 ) Then
       Return
    End If

    ! Now set up the semaphore for this context
    Call semaphore_create( ctxt, semid, error )
    If( error /= FIPC_success ) Then
       Return
    End If
    ctxt%semid = semid

    Allocate( ctxt%initialized, Stat = error )
    If( error /= 0 ) Then
       error = FIPC_allocation_failed
       Return
    End If

    error = FIPC_success

  Contains

    Subroutine set_comm( handle, comm, error )

      Integer             , Intent( In    ) :: handle
      Type( communicator ), Intent(   Out ) :: comm
      Integer             , Intent(   Out ) :: error

      comm%handle = handle

      Call mpi_comm_size( comm%handle, comm%size, error )
      If( error /= 0 ) Then
         Return
      End If

      Call mpi_comm_rank( comm%handle, comm%rank, error )
      If( error /= 0 ) Then
         Return
      End If

      comm%initialized = .True.

    End Subroutine set_comm
