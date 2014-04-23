     Function FIPC_remove_sem( semid ) Bind( C )
       Use, Intrinsic :: iso_c_binding, Only : c_int
       Implicit None
       Integer( c_int )                      :: FIPC_remove_sem
       Integer( c_int ), Value, Intent( In ) :: semid
     End Function FIPC_remove_sem
