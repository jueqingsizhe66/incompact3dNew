    Subroutine generate_base_extra_comm( world_comm, intra_comm, extra_comm, error )

      ! From the world and intra communicators generate the extra communicator
      ! Note that this is defined only on processor zero of each multicore node.

      Integer, Intent( In    ) :: world_comm
      Integer, Intent( In    ) :: intra_comm
      Integer, Intent(   Out ) :: extra_comm
      Integer, Intent(   Out ) :: error

      Integer :: world_rank
      Integer :: intra_rank
      Integer :: split_val

      Call mpi_comm_rank( world_comm, world_rank, error )
      If( error /= 0 ) Then
         Return
      End If
      Call mpi_comm_rank( intra_comm, intra_rank, error )
      If( error /= 0 ) Then
         Return
      End If
      split_val = Merge( 1, MPI_UNDEFINED, intra_rank == 0 )
      Call mpi_comm_split( world_comm, split_val, world_rank, extra_comm, error )

      error = 0

    End Subroutine generate_base_extra_comm
