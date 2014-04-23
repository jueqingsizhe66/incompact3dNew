  Subroutine segment_free( shmid, error )

    ! Free a segment

    Integer( c_int ), Intent( In   ) :: shmid
    Integer         , Intent(  Out ) :: error

    Type( segment_list_type ), Pointer :: p, prev

    Type( c_ptr ) :: shmaddr

    Integer :: retval

    ! First find the seg
    prev => Null()
    p    => seg_list
    Do While( Associated( p ) )
       If( p%data%shmid == shmid ) Then
          Exit
       End If
       prev => p
       p    => p%next
    End Do
    If( .Not. Associated( p ) ) Then
       error = FIPC_seg_not_found
       Return
    End If

    If( debug ) Then
       ! Not really neccessary, but can make debug printing a bit clearer if
       ! things get too out of step
       Call mpi_barrier( p%data%ctxt%intra_comm%handle, error )
    End If

    shmaddr = p%data%shmaddr
    ! Say "Buh Bye" everybody !
    retval = fipc_detach_seg( shmaddr )
    If( retval /= 0 ) Then
       error = FIPC_seg_detach_failed
       Return
    End If

    ! Make sure everybody is detached
    Call mpi_barrier( p%data%ctxt%intra_comm%handle, error )
    if( error /= 0 ) Then
       Return
    End if

    ! And delete the segment
    If( p%data%ctxt%intra_comm%rank == 0 ) Then
       retval = fipc_remove_seg( shmid )
       If( retval /= 0 ) Then
          error = FIPC_seg_remove_failed
          Return
       End If
    End If

    ! Make sure all up to date
    Call mpi_barrier( p%data%ctxt%intra_comm%handle, error )
    If( error /= 0 ) Then
       Return
    End if

    ! And remove the segment from the linked list
    ! Link up list around the one to die
    If( Associated( prev ) ) Then
       ! Not first in list
       prev%next => p%next
    Else
       ! First in list
       seg_list => p%next
    End If
    Deallocate( p%data%sizes )
    Deallocate( p )
    p => Null()

    error = FIPC_success

  End Subroutine segment_free
