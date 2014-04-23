  Subroutine print_seg_list( ctxt )

    Type( FIPC_ctxt ), Intent( In ) :: ctxt

    Type( segment_list_type ), Pointer :: p

    Integer( c_long ), Dimension( 1:7 ) :: buf

    Integer( c_int ) :: shmid

    Integer :: seg_n_at, rank
    Integer :: retval
    Integer :: mod_place

    Character( Len = 10 ) :: type
    Character( Len = 60 ) :: format

    If( ctxt%intra_comm%rank == 0 ) Then
       Write( *, '( "------ Shared Memory Segments --------" )' )
       Write( *, '( "shmid", t10, "type", t20, "Dimensions", t40, "Attached"  )' )
       format = '( i0, t10, a, t20, ?( i0, 1x ), t40, 1x, i0 )'
       mod_place = scan( format, '?' )
       p => seg_list
       Do While( Associated( p ) )
          shmid  = p%data%shmid
          retval = fipc_inquire_seg( shmid, Size( buf ), buf )
          seg_n_at = buf( SEG_NATTCH )
          Select Case( p%data%type )
          Case( integer_vals )
             type = 'Integer'
          Case( double_vals )
             type = 'Double'
          Case( complex_vals )
             type = 'Complex'
          Case Default
             type = 'UNKNOWN'
          End Select
          rank = Size( p%data%sizes )
          Write( format( mod_place:mod_place ), '( i1 )' ) rank
          Write( *, format ) shmid, type, p%data%sizes, seg_n_at
          p => p%next
       End Do
    End If

  End Subroutine print_seg_list
