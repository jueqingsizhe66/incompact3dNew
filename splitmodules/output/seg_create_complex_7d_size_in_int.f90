  Subroutine seg_create_complex_7d_size_in_int( ctxt, n, a, error )

    Type( FIPC_ctxt )                                             , Intent( In    ) :: ctxt
    Integer, Dimension( : )                                       , Intent( In    ) :: n
    Complex( c_double ), Dimension( :, :, :, :, :, :, : ), Pointer, Intent(   Out ) :: a
    Integer                                                       , Intent(   Out ) :: error

    Integer, Parameter :: rank = 7
    Integer, Parameter :: type = complex_vals

    Include 'FreeIPC_create.f90'

  End Subroutine seg_create_complex_7d_size_in_int
