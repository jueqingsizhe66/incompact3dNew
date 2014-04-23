  Subroutine generate_base_comms( universe_comm, world_comm, intra_comm, extra_comm, error )

    ! Try to work out the mapping of the processors onto the
    ! physical nodes, and from that generate the communicators
    ! needed to set up FIPC_ctxt_world

    Integer, Intent( In    ) :: universe_comm
    Integer, Intent(   Out ) :: world_comm
    Integer, Intent(   Out ) :: intra_comm
    Integer, Intent(   Out ) :: extra_comm
    Integer, Intent(   Out ) :: error

    ! First set the world comm.
    Call generate_base_world_comm( universe_comm, world_comm, error )
    If( error /= 0 ) Then
       Return
    End If

    ! Now the intra node communicator
    Call generate_base_intra_comm( world_comm, intra_comm, error )
    If( error /= 0 ) Then
       Return
    End If

    ! We now have an intra comm. From that it's easy to generate an extra comm.
    Call generate_base_extra_comm( world_comm, intra_comm, extra_comm, error )
    If( error /= 0 ) Then
       Return
    End If

    If( debug ) Then
       Call write_base_info( intra_comm, extra_comm, error )
       If( error /= 0 ) Then
          Return
       End If
    End If

    error = FIPC_success

  Contains

    Subroutine generate_base_world_comm( universe_comm, world_comm, error )

      ! From the universe communicator generate the base world communicator.
      ! Trivial !

      Integer, Intent( In    ) :: universe_comm
      Integer, Intent(   Out ) :: world_comm
      Integer, Intent(   Out ) :: error

      Call mpi_comm_dup( universe_comm, world_comm, error )
      If( error /= 0 ) Then
         Return
      End If

      error = 0

    End Subroutine generate_base_world_comm
