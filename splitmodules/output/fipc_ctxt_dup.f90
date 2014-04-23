  Subroutine FIPC_ctxt_dup( ctxt_1, ctxt_2, error )

    ! Duplicate the context CTXT_1, returning th new context in CTXT_2
    !
    ! On success ERROR is set to FIPC_SUCCESS. Any other value
    ! indicates error. These can be compared to the symbolic constants
    ! defined above for better diagnosis

    Type( FIPC_ctxt ), Intent( In    ) :: ctxt_1
    Type( FIPC_ctxt ), Intent(   Out ) :: ctxt_2
    Integer          , Intent(   Out ) :: error

    Integer :: world_comm_2
    Integer :: intra_comm_2
    Integer :: extra_comm_2

    ! Check things are initialized
    If( .Not. Associated( FIPC_ctxt_world%initialized ) ) Then
       error = FIPC_not_initialized
       Return
    End If

    ! Duplicate the communicators
    Call mpi_comm_dup( ctxt_1%world_comm%handle, world_comm_2, error )
    If( error /= 0 ) Then
       Return
    End If
    Call mpi_comm_dup( ctxt_1%intra_comm%handle, intra_comm_2, error )
    If( error /= 0 ) Then
       Return
    End If
    ! Extra comm only defined on rank 0 of the intra comm
    If( ctxt_1%intra_comm%rank == 0 ) Then
       Call mpi_comm_dup( ctxt_1%extra_comm%handle, extra_comm_2, error )
       If( error /= 0 ) Then
          Return
       End If
    End If

    ! Now create the new context
    Call ctxt_create( world_comm_2, intra_comm_2, extra_comm_2, ctxt_2, error )
    If( error /= 0 ) Then
       Return
    End If

    error = FIPC_success

  End Subroutine FIPC_ctxt_dup
