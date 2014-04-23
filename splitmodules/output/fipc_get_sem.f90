     Function FIPC_get_sem( create, exclusive, perms ) Bind( C )
       Use, Intrinsic :: iso_c_binding, Only : c_int
       Implicit None
       Integer( c_int )                      :: FIPC_get_sem
       Integer( c_int ), Value, Intent( In ) :: create
       Integer( c_int ), Value, Intent( In ) :: exclusive
       Integer( c_int ), Value, Intent( In ) :: perms
     End Function FIPC_get_sem
