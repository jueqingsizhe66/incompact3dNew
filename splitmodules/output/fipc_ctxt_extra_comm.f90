  Subroutine FIPC_ctxt_extra_comm( ctxt, comm, error )

    Type( FIPC_ctxt ), Intent( In    ) :: ctxt
    Integer          , Intent(   Out ) :: comm
    Integer          , Intent(   Out ) :: error

    ! Check things are initialized
    If( .Not. Associated( FIPC_ctxt_world%initialized ) ) Then
       error = FIPC_not_initialized
       Return
    End If

    If( ctxt%intra_comm%rank == 0 ) Then
       comm = ctxt%extra_comm%handle
    Else
       comm = FIPC_comm_null
    End If

    error = FIPC_success

  End Subroutine FIPC_ctxt_extra_comm
