  Function test_ne_ctxt( a, b )

    ! Test two contexts for inequality

    Logical                    :: test_ne_ctxt
    Type( FIPC_ctxt ), Intent( In ) :: a
    Type( FIPC_ctxt ), Intent( In ) :: b

    test_ne_ctxt = .Not. a == b

  End Function test_ne_ctxt
