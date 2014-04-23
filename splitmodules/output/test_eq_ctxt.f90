  Function test_eq_ctxt( a, b )

    ! Test two contexts for equality

    Logical                         :: test_eq_ctxt
    Type( FIPC_ctxt ), Intent( In ) :: a
    Type( FIPC_ctxt ), Intent( In ) :: b

    If( Associated( a%initialized ) .And. Associated( b%initialized ) ) Then
       test_eq_ctxt = a%world_comm == b%world_comm .And. &
                      a%intra_comm == b%intra_comm
       If( a%intra_comm%rank == 0 ) Then
          test_eq_ctxt = test_eq_ctxt .And. a%extra_comm == b%extra_comm
       End If
    Else
       test_eq_ctxt = .False.
    End If

  End Function test_eq_ctxt
