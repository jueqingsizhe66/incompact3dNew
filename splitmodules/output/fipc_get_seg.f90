     Function FIPC_get_seg( size, create, exclusive, perms ) Bind( C )
       Use, Intrinsic :: iso_c_binding, Only : c_int, c_long
       Implicit None
       Integer( c_int  )                      :: FIPC_get_seg
       Integer( c_long ), Value, Intent( In ) :: size
       Integer( c_int  ), Value, Intent( In ) :: create
       Integer( c_int  ), Value, Intent( In ) :: exclusive
       Integer( c_int  ), Value, Intent( In ) :: perms
     End Function FIPC_get_seg
