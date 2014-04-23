     Function FIPC_detach_seg( shmaddr ) Bind( C )
       Use, Intrinsic :: iso_c_binding, Only : c_int, c_ptr
       Implicit None
       Integer( c_int )     :: FIPC_detach_seg
       Type( c_ptr ), Value :: shmaddr ! No intent because this confuses
                                       ! the Cray compiler
     End Function FIPC_detach_seg
