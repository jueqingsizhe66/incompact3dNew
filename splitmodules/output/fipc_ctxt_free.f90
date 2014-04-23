  Subroutine FIPC_ctxt_free( ctxt, error )

    Type( FIPC_ctxt ), Intent( InOut ) :: ctxt
    Integer          , Intent(   Out ) :: error

    Type( segment_list_type ), Pointer :: p

    ! Check things are initialized
    If( .Not. Associated( FIPC_ctxt_world%initialized ) ) Then
       error = FIPC_not_initialized
       Return
    End If

    ! Check not trying to free the base context
    If( ctxt == FIPC_ctxt_world ) Then
       error = FIPC_freeing_ctxt_world
       Return
    End If

    ! Check if anybody is still using this context
    p => seg_list
    Do While( Associated( p ) )
       If( p%data%ctxt == ctxt ) Then
          error = FIPC_ctxt_in_use
          Return
       End If
       p => p%next
    End Do

    Call ctxt_free( ctxt, error )
    If( error /= FIPC_SUCCESS ) Then
       Return
    End If

    error = FIPC_success

  End Subroutine FIPC_ctxt_free
