  Subroutine seg_free_integer_3d( a, error )

    Integer, Dimension( :, :, : ), Pointer, Intent( InOut ) :: a
    Integer                               , Intent(   Out ) :: error

    Integer, Parameter :: type = integer_vals

    Integer, Dimension( :, :, : ), Pointer :: fp

    Include 'FreeIPC_free.f90'

  End Subroutine seg_free_integer_3d
