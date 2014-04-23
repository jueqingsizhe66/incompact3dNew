  Subroutine FIPC_ctxt_split( ctxt_1, colour, key, ctxt_2, error )

    Type( FIPC_ctxt ), Intent( In    ) :: ctxt_1
    Integer          , Intent( In    ) :: colour
    Integer          , Intent( In    ) :: key
    Type( FIPC_ctxt ), Intent(   Out ) :: ctxt_2
    Integer          , Intent(   Out ) :: error

    Integer :: world_comm_2
    Integer :: intra_comm_2
    Integer :: extra_comm_2

    Integer :: intra_rank, is_node_0

    ! Check things are initialized
    If( .Not. Associated( FIPC_ctxt_world%initialized ) ) Then
       error = FIPC_not_initialized
       Return
    End If

    ! Split the communicators
    Call mpi_comm_split( ctxt_1%world_comm%handle, colour, key, world_comm_2, error )
    If( error /= 0 ) Then
       Return
    End If
    Call mpi_comm_split( ctxt_1%intra_comm%handle, colour, key, intra_comm_2, error )
    If( error /= 0 ) Then
       Return
    End If
    ! Extra comm only defined on rank 0 of the intra comm
    ! BUT WHAT IF SPLIT WITHIN THE INTRA COMM. Need to think here ....
    Call mpi_comm_rank( intra_comm_2, intra_rank, error )
    If( error /= 0 ) Then
       Return
    End If
    is_node_0 = merge( 1, MPI_UNDEFINED, intra_rank == 0 )
    Call mpi_comm_split( world_comm_2, colour, key, extra_comm_2, error )
    If( error /= 0 ) Then
       Return
    End If

    ! Now create the new context if required
    If( colour /= FIPC_undefined ) Then
       Call ctxt_create( world_comm_2, intra_comm_2, extra_comm_2, ctxt_2, error )
    End If
    If( error /= 0 ) Then
       Return
    End If

    error = FIPC_success

  End Subroutine FIPC_ctxt_split
