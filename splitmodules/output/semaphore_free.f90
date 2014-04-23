  Subroutine semaphore_free( semid, error )

    ! Free a semaphore

    Integer( c_int ), Intent( In   ) :: semid
    Integer         , Intent(  Out ) :: error

    Type( semaphore_list_type ), Pointer :: p, prev

    Integer :: retval
    Integer :: rank

    ! First find the sem
    prev => Null()
    p    => sem_list
    Do While( Associated( p ) )
       If( p%semid == semid ) Then
          Exit
       End If
       prev => p
       p    => p%next
    End Do
    If( .Not. Associated( p ) ) Then
       error = FIPC_sem_not_found
       Return
    End If

    If( debug ) Then
       ! Not really neccessary, but can make debug printing a bit clearer if
       ! things get too out of step
       Call mpi_barrier( p%intra_handle, error )
    End If

    ! And delete the semaphore
    Call mpi_comm_rank( p%intra_handle, rank, error )
    If( error /= 0 ) Then
       Return
    End if
    If( rank == 0 ) Then
       retval = fipc_remove_sem( semid )
       If( retval /= 0 ) Then
          error = FIPC_sem_remove_failed
          Return
       End If
    End If

    ! Make sure all up to date
    Call mpi_barrier( p%intra_handle, error )
    If( error /= 0 ) Then
       Return
    End if

    ! And remove the semaphore from the linked list
    ! Link up list around the one to die
    If( Associated( prev ) ) Then
       ! Not first in list
       prev%next => p%next
    Else
       ! First in list
       sem_list => p%next
    End If
    Deallocate( p )
    p => Null()

    error = FIPC_success

  End Subroutine semaphore_free
