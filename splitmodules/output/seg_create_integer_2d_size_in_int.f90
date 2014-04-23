  Subroutine seg_create_integer_2d_size_in_int( ctxt, n, a, error )

    Type( FIPC_ctxt )                  , Intent( In    ) :: ctxt
    Integer, Dimension( : )            , Intent( In    ) :: n
    Integer, Dimension( :, : ), Pointer, Intent(   Out ) :: a
    Integer                            , Intent(   Out ) :: error

    Integer, Parameter :: rank = 2
    Integer, Parameter :: type = integer_vals

    Include 'FreeIPC_create.f90'

  End Subroutine seg_create_integer_2d_size_in_int
