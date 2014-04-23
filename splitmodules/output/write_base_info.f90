    Subroutine write_base_info( intra_comm, extra_comm, error )

      Integer, Intent( In    ) :: intra_comm
      Integer, Intent( In    ) :: extra_comm
      Integer, Intent(   Out ) :: error

      Integer :: intra_size, intra_rank
      Integer :: extra_size, extra_rank
      Integer :: i

      ! Write out a bit more data
      Call mpi_comm_size( intra_comm, intra_size, error )
      If( error /= 0 ) Then
         Return
      End If
      Call mpi_comm_rank( intra_comm, intra_rank, error )
      If( error /= 0 ) Then
         Return
      End If
      If( intra_rank == 0 ) Then
         Call mpi_comm_size( extra_comm, extra_size, error )
         If( error /= 0 ) Then
            Return
         End If
         Call mpi_comm_rank( extra_comm, extra_rank, error )
         If( error /= 0 ) Then
            Return
         End If
         If( extra_rank == 0 ) Then
            Write( *, * ) 'Number of nodes found: ', extra_size
         End If
         Do i = 0 , extra_size
            Call mpi_barrier( extra_comm, error )
            If( i == extra_rank ) Then
               Write( *, * ) 'Number of processes on node ', extra_rank, ' is ', intra_size
            End If
         End Do
      End If

      error = 0

    End Subroutine write_base_info
