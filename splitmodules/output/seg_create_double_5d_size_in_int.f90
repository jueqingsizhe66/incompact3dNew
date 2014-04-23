  Subroutine seg_create_double_5d_size_in_int( ctxt, n, a, error )

    Type( FIPC_ctxt )                                    , Intent( In    ) :: ctxt
    Integer, Dimension( : )                              , Intent( In    ) :: n
    Real( c_double ), Dimension( :, :, :, :, : ), Pointer, Intent(   Out ) :: a
    Integer                                              , Intent(   Out ) :: error

    Integer, Parameter :: rank = 5
    Integer, Parameter :: type = double_vals

    Include 'FreeIPC_create.f90'

  End Subroutine seg_create_double_5d_size_in_int
