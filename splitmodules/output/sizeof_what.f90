  Function sizeof_what( what ) Result( r )

    Integer( c_long ) :: r

    Integer, Intent( In ) :: what

    Select Case( what )
    Case( integer_vals )
       r = FIPC_sizeof_c_int()
    Case( double_vals )
       r = FIPC_sizeof_c_double()
    Case( complex_vals )
       r = FIPC_sizeof_c_complex()
    Case Default
       r = -1
    End Select

  End Function sizeof_what
