    Subroutine generate_base_intra_comm( world_comm, intra_comm, error )

      ! Generate the intra node communicator. This is the difficult one.
      ! The basic idea is to find which other processes can attach to
      ! a shared memory segment created by a reference process.

      Integer, Intent( In    ) :: world_comm
      Integer, Intent(   Out ) :: intra_comm
      Integer, Intent(   Out ) :: error

      Type( c_ptr ) :: shmaddr
      Type( c_ptr ) :: zero_shmaddr

      Integer( c_int ), Pointer, Dimension( : ) :: seg      => Null()
      Integer( c_int ), Pointer, Dimension( : ) :: zero_seg => Null()

      Integer( c_int  ), Dimension( 1:10 ) :: zero_seg_data
      Integer( c_long ), Dimension( 1:7  ) :: shm_data
      Integer( c_long ), Dimension( 1:7  ) :: zero_shm_data

      Integer( c_int ) :: shmid, zero_shmid

      Integer :: world_size, world_rank
      Integer :: sizeof_c_int, sizeof_c_long
      Integer :: c_long_mpi_handle, c_int_mpi_handle
      Integer :: got_it, split_val
      Integer :: procs_found
      Integer :: remainder_comm, new_comm, got_comm
      Integer :: rem_rank, rem_size
      Integer :: got_rank, got_size
      Integer :: retval

      Logical :: found_seg
      Logical :: i_am_ref_proc

      Integer :: base_rank
      Call mpi_comm_rank( mpi_comm_world, base_rank, error )

      Call mpi_comm_size( world_comm, world_size, error )
      If( error /= 0 ) Then
         Return
      End If
      Call mpi_comm_rank( world_comm, world_rank, error )
      If( error /= 0 ) Then
         Return
      End If

      ! Now the difficult one, the intra node. Basic idea is to see
      ! who can attach to a shared memory seg that I create.

      ! Bit of initialisation
      shmaddr      = c_null_ptr
      zero_shmaddr = c_null_ptr

      ! Work out how to interface with the preverted C view of things
      ! Avoid mpi_sizeof because of broken implementation on Cray XT4 series.
      ! Also avoid mpi_type_create_f90_integer because of broken implementation
      ! on Cray XT4 series. ....
      sizeof_c_int  = FIPC_sizeof_c_int()
      sizeof_c_long = FIPC_sizeof_c_long()
      Call mpi_type_match_size( MPI_TYPECLASS_INTEGER, sizeof_c_int, c_int_mpi_handle, error )
      If( error /= 0 ) Then
         Return
      End If
      Call mpi_type_match_size( MPI_TYPECLASS_INTEGER, sizeof_c_long, c_long_mpi_handle, error )
      If( error /= 0 ) Then
         Return
      End If

      ! Each process creates a small seg. This is done in exclusive mode
      ! so each proc gets a different seg.
      Call get_new_seg( Int( 10 * sizeof_c_int, c_long ), shmid )
      If( shmid < 0 ) Then
         error = FIPC_seg_get_failed
         Return
      End If

      If( debug ) Then
         Write( *, * ) 'shmid ', shmid, ' on ', world_rank
      End If

      ! Try to attach to the seg
      shmaddr = fipc_attach_seg( shmid, SEG_NOREADONLY )
      If( debug ) Then
         Write( *, * ) 'Assoc status of seg ', c_associated( shmaddr ), ' on ', world_rank
      End If
      If( .Not.  c_associated( shmaddr ) ) Then
         error = FIPC_seg_attach_failed
         Return
      End If

      ! Fill the seg with some data unique to this node using Fortran
      Call c_f_pointer( shmaddr, seg, (/ 10 /) )
      Call date_and_time( values = seg( 1:8 ) )
      seg( 9  ) = shmid
!!$      seg( 10 ) = world_rank
      seg( 10 ) = base_rank

      ! Basic idea is to loop over all the procs in the world and see if
      ! the other procs can attach to the seg created by the reference proc.
      ! If they can check carefully that they are really on the same node,
      ! and if they are cut them out of list of procs to examine next time
      ! around the loop. This requires a bit of fancy communicator handling ...
      ! The most important points are
      ! a) remainder_comm spans the procs that we haven't managed yet to find a node
      !    for
      ! b) got_comm spans the list of candidate procs to be on the same node as
      !    the reference proc
      ! c) intra_comm is the final result

      ! Obviously to begin with we have all the processors to examine
      Call mpi_comm_dup( world_comm, remainder_comm, error )
      If( error /= 0 ) Then
         Return
      End If
      Call mpi_comm_size( remainder_comm, rem_size, error )
      If( error /= 0 ) Then
         Return
      End If

      ! While there are processors spanned by the remainder communicator ...
      Do While( rem_size > 0 )

         ! Proc zero in the remainer comm will be the reference proc
         Call mpi_comm_rank( remainder_comm, rem_rank, error )
         If( error /= 0 ) Then
            Return
         End If
         i_am_ref_proc = rem_rank == 0

         ! Tell all procs not yet assigned the shmid of the seg on the ref proc
         If( i_am_ref_proc ) Then
            zero_shmid = shmid
         End If
         Call mpi_bcast( zero_shmid, 1, mpi_integer, 0, remainder_comm, error )
         If( error /= 0 ) Then
            Return
         End If

         ! Now all processes except the ref proc attempt to attach to such a seg
         ! Have to be careful - if the shmid of my seg is the same as that on the
         ! ref proc I MUST be on a different node as get_new_seg creates segments
         ! in an exclusive mode, yet I'll still be able to attach. Catch this case
         If( .Not. i_am_ref_proc ) Then
            If( zero_shmid /= shmid ) Then
               zero_shmaddr = fipc_attach_seg( zero_shmid, SEG_READONLY )
               found_seg = c_associated( zero_shmaddr )
            Else
               zero_shmaddr = c_null_ptr
               found_seg = .False.
            End If
         Else
            found_seg    = .True.
            zero_shmaddr = shmaddr
         End If

         ! Get a communicator spanning the candidate procs for being
         ! on the same node as the ref proc
         If( found_seg ) Then
            split_val = 1
            If( debug ) Then
               got_it = 1
            End If
         Else
            split_val = MPI_UNDEFINED
            If( debug ) Then
               got_it = 0
            End If
         End If
         If( debug ) Then
            Call mpi_allreduce( got_it, procs_found, 1, mpi_integer, &
                 mpi_sum, remainder_comm, error )
            If( error /= 0 ) Then
               Return
            End If
            If( rem_rank == 0 ) Then
               Write( *, * ) 'procs found ', procs_found, rem_size
            End If
         End If
         Call mpi_comm_split( remainder_comm, split_val, rem_rank, new_comm, error )
         If( error /= 0 ) Then
            Return
         End If

         ! Now check carefully that the seg that these processors attached to
         ! REALLY was the seg on the ref proc - outside chance that two different nodes
         ! have a seg with the same shmid
         If( split_val == 1 ) Then

            got_comm = new_comm
            ! First check that the data describing the seg is the same as that on the ref proc
            retval = fipc_inquire_seg( zero_shmid, Size( shm_data ), shm_data )
            If( retval /= 0 ) Then
               error = FIPC_seg_inquire_failed
               Return
            End If
            If( i_am_ref_proc ) Then
               zero_shm_data = shm_data
            End If
            Call mpi_bcast( zero_shm_data, Size( zero_shm_data ), &
                 c_long_mpi_handle, 0, got_comm, error )
            If( error /= 0 ) Then
               Return
            End If
            split_val = Merge( 1, MPI_UNDEFINED, All( zero_shm_data == shm_data ) )

            ! Split off any processors that failed ...
            Call mpi_comm_split( new_comm, split_val, rem_rank, got_comm, error )
            If( error /= 0 ) Then
               Return
            End If

            ! Finished with first guess at node comm for the reference processor
            Call mpi_comm_free( new_comm, error )
            If( error /= 0 ) Then
               Return
            End If

            ! OK, turn the paranoia up to 11 ! Check that for all who survived the
            ! above test that the data in the reference seg is what we put there !!
            If( split_val == 1 ) Then
               ! Get the data on the reference processor
               If( i_am_ref_proc ) Then
                  zero_seg_data = seg
               End If
               Call mpi_bcast( zero_seg_data, Size( zero_seg_data ), &
                    c_int_mpi_handle, 0, got_comm, error )
               ! Get the Fortran pointer to the shared memory
               Call c_f_pointer( zero_shmaddr, zero_seg, (/ 10 /) )
               ! Check the data
               split_val = Merge( 1, MPI_UNDEFINED, All( zero_seg_data == zero_seg ) )
               ! Done with the Fortran pointers
               seg      => Null()
               zero_seg => Null()

               ! And get rid of the outstanding communicator
               Call mpi_comm_free( got_comm, error )
               If( error /= 0 ) Then
                  Return
               End If

            End If

         End If

         ! And at last our guess at the node comm ! Split the procs into two sets
         ! a) Those we believe are on the smae node as the reference proc
         ! b) all the others - i.e. the ones we haven't placed yet
         If( split_val == MPI_UNDEFINED ) Then
            split_val = 2
         End If
         Call mpi_comm_split( remainder_comm, split_val, rem_rank, new_comm, error )
         If( error /= 0 ) Then
            Return
         End If

         ! Done with old version of remainder comm.
         Call mpi_comm_free( remainder_comm, error )
         If( error /= 0 ) Then
            Return
         End If

         ! Store our guess at the intranode communicator and get it's details
         If( split_val == 1 ) Then
            Call mpi_comm_dup( new_comm, intra_comm, error )
            If( error /= 0 ) Then
               Return
            End If
            Call mpi_comm_size( intra_comm, got_size, error )
            If( error /= 0 ) Then
               Return
            End If
            Call mpi_comm_rank( intra_comm, got_rank, error )
            If( error /= 0 ) Then
               Return
            End If
         End If

         ! One more sanity check. The size of the intranode communicator
         ! should be the same as the number of attaches to the segment
         ! associated with the reference proc
!!$         If( i_am_ref_proc ) Then
!!$            If( got_size /= shm_data( SEG_NATTCH ) ) Then
!!$               Write( *, * ) 'sanity: ', got_size, shm_data( SEG_NATTCH )
!!$               error = FIPC_sanity_failed
!!$               Return
!!$            End If
!!$         End If

         ! Detach all procs that managed to attach to a seg with the same shmid
         ! as the ref_proc. Know we can't be detaching from our "own" seg as
         ! we caught that above
         If( found_seg ) Then
            If( .Not. i_am_ref_proc ) Then
               retval = fipc_detach_seg( zero_shmaddr )
               If( retval /= 0 ) Then
                  error = FIPC_seg_detach_failed
                  Return
               End If
            End If
         End If

         ! If we found an intracomm we're done !
         If( split_val == 1 ) Then
            Exit ! Hallelujah !!!
         End If

         ! For all left the new comm is simply the comm spanning the
         ! remaining unassigned procs
         remainder_comm = new_comm

         ! Number of procs left to deal with is simply the size of the
         ! communicator spanning them
         Call mpi_comm_size( remainder_comm, rem_size, error )
         If( error /= 0 ) Then
            Return
         End If

      End Do

      ! Clear up outstanding comm
      Call mpi_comm_free( new_comm, error )
      If( error /= 0 ) Then
         Return
      End If

      ! Be careful - avoid possible race conditions
      Call mpi_barrier( world_comm, error )
      If( error /= 0 ) Then
         Return
      End If

      ! Clear up the segs owned by me
      retval = fipc_detach_seg( shmaddr )
      If( retval /= 0 ) Then
         error = FIPC_seg_detach_failed
         Return
      End If
      retval = fipc_remove_seg( shmid )
      If( retval /= 0 ) Then
         error = FIPC_seg_remove_failed
         Return
      End If

      error = 0

    End Subroutine generate_base_intra_comm
