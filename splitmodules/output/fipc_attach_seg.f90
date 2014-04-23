     Function FIPC_attach_seg( shmid, flag ) Bind( C )
       Use, Intrinsic :: iso_c_binding, Only : c_int, c_ptr
       Implicit None
       Type( c_ptr )                          :: FIPC_attach_seg
       Integer( c_int  ), Value, Intent( In ) :: shmid
       Integer( c_int  ), Value, Intent( In ) :: flag
     End Function FIPC_attach_seg
