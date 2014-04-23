  Subroutine FIPC_finalize( error )

    ! Finalize FIPC.
    !
    ! On success ERROR is set to FIPC_SUCCESS. Any other value
    ! indicates error. These can be compared to the symbolic constants
    ! defined above for better diagnosis

    Integer, Intent( Out ) :: error

    ! Check things are initialized
    If( .Not. Associated( FIPC_ctxt_world%initialized ) ) Then
       error = FIPC_not_initialized
       Return
    End If

    ! Free all the segments we know about
    Do While( Associated( seg_list ) )
       Call segment_free( seg_list%data%shmid, error )
       If( error /= FIPC_success ) Then
          Return
       End If
    End Do

    ! And free the context
    Call ctxt_free( FIPC_ctxt_world, error )
    If( error /= 0 ) Then
       Return
    End If

    ! Finally free any outstanding semaphores we know about
    Do While( Associated( sem_list ) )
       Call semaphore_free( sem_list%semid, error )
       If( error /= FIPC_success ) Then
          Return
       End If
    End Do

    error = FIPC_success

  End Subroutine FIPC_finalize
