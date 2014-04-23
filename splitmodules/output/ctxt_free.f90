  Subroutine ctxt_free( ctxt, error )

    ! Free a context.
    !
    ! On success ERROR is set to FIPC_success
    ! On error error is any other value than FIPC_success

    Type( FIPC_ctxt ), Intent( InOut ) :: ctxt
    Integer          , Intent(   Out ) :: error

    Integer :: comm

    ! First free the semaphore
    Call semaphore_free( ctxt%semid, error )
    If( error /= 0 ) Then
       Return
    End If

    ! And now the communicators
    If( ctxt%intra_comm%rank == 0 ) Then
       Call mpi_comm_free( ctxt%extra_comm%handle, error )
       If( error /= 0 ) Then
          Return
       End If
       ctxt%extra_comm%initialized = .False.
       Deallocate( ctxt%extra_comm )
    End If

    ! Take a copy of the intra_node communicator so can sync at end
    Call mpi_comm_dup( ctxt%intra_comm%handle, comm, error )
    If( error /= 0 ) Then
       Return
    End If
    Call mpi_comm_free( ctxt%intra_comm%handle, error )
    If( error /= 0 ) Then
       Return
    End If
    ctxt%intra_comm%initialized = .False.
    Deallocate( ctxt%intra_comm )

    Call mpi_comm_free( ctxt%world_comm%handle, error )
    If( error /= 0 ) Then
       Return
    End If
    ctxt%world_comm%initialized = .False.
    Deallocate( ctxt%world_comm )

    ! Make sure all in sync on the node, carefully using the
    ! copy of the intra node context
    Call mpi_barrier( comm, error )
    If( error /= 0 ) Then
       Return
    End If
    Call mpi_comm_free( comm, error )
    If( error /= 0 ) Then
       Return
    End If

    Deallocate( ctxt%initialized )

    error = FIPC_success

  End Subroutine ctxt_free
