  Subroutine FIPC_node_barrier( ctxt, error )

    Type( FIPC_ctxt ), Intent( In    ) :: ctxt
    Integer          , Intent(   Out ) :: error

    Call mpi_barrier( ctxt%intra_comm%handle, error )

  End Subroutine FIPC_node_barrier
