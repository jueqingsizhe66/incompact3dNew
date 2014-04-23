     Function FIPC_remove_seg( shmid ) Bind( C )
       Use, Intrinsic :: iso_c_binding, Only : c_int
       Implicit None
       Integer( c_int )                      :: FIPC_remove_seg
       Integer( c_int ), Value, Intent( In ) :: shmid
     End Function FIPC_remove_seg
