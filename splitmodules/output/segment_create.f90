  Subroutine segment_create( what, rank, sizes, ctxt, shmaddr, error )

    ! Create a segment of size N bytes, attach it to all procs
    ! on this node in the context CTXT, and return a c pointer
    ! to it in SHMADDR. On success ERROR is equal to FIPC_SUCCESS,
    ! any other value indicates error. On failure SHMADDR is also
    ! set to C_NULL_PTR

    Integer                          , Intent( In    ) :: what
    Integer                          , Intent( In    ) :: rank
    Integer( c_long ), Dimension( : ), Intent( In    ) :: sizes
    Type( FIPC_ctxt )                , Intent( In    ) :: ctxt
    Type( c_ptr     )                                  :: shmaddr ! No intent for Cray compiler
    Integer                          , Intent(   Out ) :: error

    Type( segment_list_type ), Pointer :: p, this

    Integer( c_int ) :: shmid

    ! Create the seg on proc 0 in this context on this node
    If( ctxt%intra_comm%rank == 0 ) Then
       Call get_new_seg( sizeof_what( what ) * Product( sizes ), shmid )
       If( shmid < 0 ) Then
          error = FIPC_seg_get_failed
          shmaddr = c_null_ptr
          Return
       End If
    End If

    ! Tell the other procs on the node
    Call mpi_bcast( shmid, 1, mpi_integer, 0, ctxt%intra_comm%handle, error )
    If( error /= 0 ) Then
       shmaddr = c_null_ptr
       Return
    End If

    ! Attach to the seg
    shmaddr = fipc_attach_seg( shmid, SEG_NOREADONLY )
    If( debug ) Then
       Write( *, * ) 'assoc status of seg in segment_create', c_associated( shmaddr )
    End If
    If( .Not.  c_associated( shmaddr ) ) Then
       error = FIPC_seg_attach_failed
       shmaddr = c_null_ptr
       Return
    End If

    ! Add data on the seg to the end of linked list about created segs
    p => seg_list
    If( .Not. Associated( p ) ) Then
       ! First element
       Allocate( seg_list, Stat = error )
       If( error /= 0 ) Then
          error = FIPC_allocation_failed
          shmaddr = c_null_ptr
          Return
       End If
       this => seg_list
    Else
       Do While( Associated( p%next ) )
          p => p%next
       End Do
       Allocate( p%next, Stat = error )
       If( error /= 0 ) Then
          error = FIPC_allocation_failed
          shmaddr = c_null_ptr
          Return
       End If
       this => p%next
    End If

    this%data%shmid   = shmid
    this%data%shmaddr = shmaddr
    this%data%ctxt    = ctxt
    this%data%type    = what
    Allocate( this%data%sizes( 1:rank ), Stat = error )
    If( error /= 0 ) Then
       error = FIPC_allocation_failed
       Return
    End If
    this%data%sizes = sizes( 1:rank )
    this%next         => Null()

    ! Make sure all up to date
    Call mpi_barrier( ctxt%intra_comm%handle, error )
    If( error /= 0 ) Then
       shmaddr = c_null_ptr
       Return
    End if

    If( debug ) Then
       Call print_seg_list( ctxt )
    End If

    error = FIPC_success

  End Subroutine segment_create
