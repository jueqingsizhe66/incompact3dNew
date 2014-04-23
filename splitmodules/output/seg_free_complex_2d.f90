  Subroutine seg_free_complex_2d( a, error )

    Complex( c_double ), Dimension( :, : ), Pointer, Intent( InOut ) :: a
    Integer                                        , Intent(   Out ) :: error

    Integer, Parameter :: type = complex_vals

    Complex( c_double ), Dimension( :, : ), Pointer :: fp

    Include 'FreeIPC_free.f90'

  End Subroutine seg_free_complex_2d
