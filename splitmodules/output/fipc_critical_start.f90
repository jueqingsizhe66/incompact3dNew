  Subroutine FIPC_critical_start( ctxt, error )

    Type( FIPC_ctxt ), Intent( In    ) :: ctxt
    Integer          , Intent(   Out ) :: error

    ! Check things are initialized
    If( .Not. Associated( FIPC_ctxt_world%initialized ) ) Then
       error = FIPC_not_initialized
       Return
    End If

    error = fipc_crit_start( ctxt%semid )

    If( error /= 0 ) Then
       error = FIPC_critical_start_failed
       Return
    End If

    error = FIPC_success

  End Subroutine FIPC_critical_start
