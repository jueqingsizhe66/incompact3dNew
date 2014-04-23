  Function test_eq_comm( a, b )

    ! Test two communicators for equality

    Logical                            :: test_eq_comm
    Type( communicator ), Intent( In ) :: a
    Type( communicator ), Intent( In ) :: b

    If( a%initialized .And. b%initialized ) Then
       test_eq_comm = a%handle == b%handle .And. &
                      a%size   == b%size   .And. &
                      a%rank   == b%rank
    Else
       test_eq_comm = .False.
    End If

  End Function test_eq_comm
