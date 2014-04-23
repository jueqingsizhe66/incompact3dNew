  Subroutine semaphore_create( ctxt, semid, error )

    ! Create a semaphore of size N bytes in the context CTXT, and return
    ! it's id in SEMID. On success ERROR is equal to FIPC_SUCCESS,
    ! any other value indicates error.

    Type( FIPC_ctxt ), Intent( In    ) :: ctxt
    Integer( c_int ) , Intent(   Out ) :: semid
    Integer          , Intent(   Out ) :: error

    Type( semaphore_list_type ), Pointer :: p, this

    ! Create the sem on proc 0 in this context on this node
    If( ctxt%intra_comm%rank == 0 ) Then
       Call get_new_sem( semid )
       If( debug ) Then
          Write( *, * ) 'semid ', semid
       End If
       If( semid < 0 ) Then
          error = FIPC_sem_get_failed
          Return
       End If
    End If

    ! Tell the other procs on the node
    Call mpi_bcast( semid, 1, mpi_integer, 0, ctxt%intra_comm%handle, error )
    If( error /= 0 ) Then
       Return
    End If

    ! Add data on the sem to the end of linked list about created sems
    p => sem_list
    If( .Not. Associated( p ) ) Then
       ! First element
       Allocate( sem_list, Stat = error )
       If( error /= 0 ) Then
          error = FIPC_allocation_failed
          Return
       End If
       this => sem_list
    Else
       Do While( Associated( p%next ) )
          p => p%next
       End Do
       Allocate( p%next, Stat = error )
       If( error /= 0 ) Then
          error = FIPC_allocation_failed
          Return
       End If
       this => p%next
    End If

    this%semid        = semid
    this%intra_handle = ctxt%intra_comm%handle

    ! Make sure all up to date
    Call mpi_barrier( ctxt%intra_comm%handle, error )
    If( error /= 0 ) Then
       Return
    End if

    error = FIPC_success

  End Subroutine semaphore_create
