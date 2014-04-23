  Subroutine seg_free_double_4d( a, error )

    Real( c_double ), Dimension( :, :, :, : ), Pointer, Intent( InOut ) :: a
    Integer                                           , Intent(   Out ) :: error

    Integer, Parameter :: type = double_vals

    Real( c_double ), Dimension( :, :, :, : ), Pointer :: fp

    Include 'FreeIPC_free.f90'

  End Subroutine seg_free_double_4d
