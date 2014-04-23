  Function test_ne_comm( a, b )

    ! Test two communicators for inequality

    Logical                            :: test_ne_comm
    Type( communicator ), Intent( In ) :: a
    Type( communicator ), Intent( In ) :: b

    test_ne_comm = .Not. a == b

  End Function test_ne_comm
