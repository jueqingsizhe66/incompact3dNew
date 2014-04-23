     Function FIPC_crit_start( semid ) Bind( C )
       Use, Intrinsic :: iso_c_binding, Only : c_int
       Implicit None
       Integer( c_int )                      :: FIPC_crit_start
       Integer( c_int ), Value, Intent( In ) :: semid
     End Function FIPC_crit_start
