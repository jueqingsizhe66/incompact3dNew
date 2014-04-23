  Subroutine FIPC_ctxt_intra_comm( ctxt, comm, error )

    Type( FIPC_ctxt ), Intent( In    ) :: ctxt
    Integer          , Intent(   Out ) :: comm
    Integer          , Intent(   Out ) :: error

    ! Check things are initialized
    If( .Not. Associated( FIPC_ctxt_world%initialized ) ) Then
       error = FIPC_not_initialized
       Return
    End If

    comm = ctxt%intra_comm%handle

    error = FIPC_success

  End Subroutine FIPC_ctxt_intra_comm
