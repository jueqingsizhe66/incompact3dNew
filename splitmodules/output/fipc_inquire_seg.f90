     Function FIPC_inquire_seg( shmid, n, shm_data ) Bind( C )
       Use, Intrinsic :: iso_c_binding, Only : c_int, c_long
       Implicit None
       Integer( c_int )                                   :: FIPC_inquire_seg
       Integer( c_int  ), Value         , Intent( In    ) :: shmid
       Integer( c_int  ), Value         , Intent( In    ) :: n
       Integer( c_long ), Dimension( * ), Intent(   Out ) :: shm_data
     End Function FIPC_inquire_seg
