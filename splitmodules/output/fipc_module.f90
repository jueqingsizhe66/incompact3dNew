Module FIPC_module

  ! A library to allow Fortran programs that use MPI to easily access
  ! shared memory within multicore nodes

  ! It uses standard Fortran 2003 and a few thin C wrappers. Fortran
  ! 2003 C interoperability is used extensively, however the rest of
  ! the code is standard Fortran 95, expcept for the use of allocatble
  ! components of derived types - again a standard Fortran 2003 language
  ! feature.

  ! The shared memory facilities are provided by use of the standard System V
  ! IPC facilities. See
  !
  ! http://www.opengroup.org/onlinepubs/009695399/functions/xsh_chap02_07.html#tag_02_07
  !
  ! for details

  Use, Intrinsic :: iso_c_binding, Only : c_int, c_long, c_double, c_ptr, &
       c_f_pointer, c_null_ptr, c_associated, c_loc

  Use mpi

  Implicit None

  !
  ! Little type to handle communicators and associated data
  !
  !!!> \cond
  Type, Private :: communicator
     Private
     !> \private
     Logical :: initialized = .False.
     !> \private
     Integer :: handle
     !> \private
     Integer :: size
     !> \private
     Integer :: rank
  End Type communicator
!!!> \endcond

  !>
  !! \brief An opaque type that represents a FIPC context.
  !!
  !! This type represents a FIPC context. It is an opaque data type with no public components
  !! All operations occur within a context. It is very similar to
  !! to a mpi communicator with a bit of extra stuff held to look
  !! after the shared memory parts of the architecture.
  !!
  !! \if for_fipc_implementors
  !! It consists of a:
  !! INITIALIZED: Is this context set up ?
  !!              Note we only use the allocation status of the pointer to check this
  !! WORLD_COMM : The communicator spanning all processes in this context
  !! INTRA_COMM : A communicator spanning all process  in the context on this shared memory node
  !! EXTRA_COMM : A communicator spanning all the process zeros in all the intra_comms
  !!              on all the nodes which hold processes in WORLD_COMM
  !! SEMID      : The handle of a semaphore shared by members of the INTRA_COMM
  !> \endif
  !
  Type, Public :: FIPC_ctxt
     Private
     !> \private
     Logical             , Pointer :: initialized => Null()
     !> \private
     Type( communicator ), Pointer :: world_comm
     !> \private
     Type( communicator ), Pointer :: intra_comm
     !> \private
     Type( communicator ), Pointer :: extra_comm
     !> \private
     Integer                       :: semid
  End Type FIPC_ctxt

  !> The "default" context. Compare MPI_COMM_WORLD.
  Type( FIPC_ctxt ), Save, Public :: FIPC_ctxt_world

  ! Error flags. Note we avoid clashes with the mpi error flags.

  !> Return value to indicate succesful completion
  Integer, Parameter, Public :: FIPC_success               =  0

  !> Return value indicating an attempt was made to initialise FreeIPC when it was already initialis
  Integer, Parameter, Public :: FIPC_already_initialized   =  1 + MPI_ERR_LASTCODE

  !> Return value indicating a memory allocation failed within FreeIPC
  Integer, Parameter, Public :: FIPC_allocation_failed     =  2 + MPI_ERR_LASTCODE

  !> Return value indicating that an attempt to get a shared memory segment failed
  Integer, Parameter, Public :: FIPC_seg_get_failed        =  3 + MPI_ERR_LASTCODE

  !> Return value indicating that an attempt to attach to a shared memory segment failed
  Integer, Parameter, Public :: FIPC_seg_attach_failed     =  4 + MPI_ERR_LASTCODE

  !> Return value indicating that an attempt to inquire the properties shared memory segment failed
  Integer, Parameter, Public :: FIPC_seg_inquire_failed    =  5 + MPI_ERR_LASTCODE

  !> Return value indicating that an attempt to detach from shared memory segment failed
  Integer, Parameter, Public :: FIPC_seg_detach_failed     =  6 + MPI_ERR_LASTCODE

  !> Return value indicating that an attempt to remove a shared memory segment failed
  Integer, Parameter, Public :: FIPC_seg_remove_failed     =  7 + MPI_ERR_LASTCODE

  !> Return value indicating that a FreeIPC internal consistency check failed
  Integer, Parameter, Public :: FIPC_sanity_failed         =  8 + MPI_ERR_LASTCODE

  !> Return value indicating that an attempt was made to free FIPC_ctxt_world
  Integer, Parameter, Public :: FIPC_freeing_ctxt_world    =  9 + MPI_ERR_LASTCODE

  !> Return value indicating that FreeIPC could not identify which shared memory segment is to be fr
  Integer, Parameter, Public :: FIPC_seg_not_found         = 10 + MPI_ERR_LASTCODE

  !> Return value indicating that an attempt was made to free a NULL pointer
  Integer, Parameter, Public :: FIPC_free_null_pointer     = 11 + MPI_ERR_LASTCODE

  !> Return value indicating that an insufficent dimensions were supplied when allocating a shared m
  Integer, Parameter, Public :: FIPC_insufficient_dims     = 12 + MPI_ERR_LASTCODE

  !> Return value indicating that FreeIPC was not initalized
  Integer, Parameter, Public :: FIPC_not_initialized       = 13 + MPI_ERR_LASTCODE

  !> Return value indicating that an attempt was made to free a context that was still in use, for i
  !> a shared memory segment still exists within that context
  Integer, Parameter, Public :: FIPC_ctxt_in_use           = 14 + MPI_ERR_LASTCODE

  !> Return value indicating that an attempt to get a semaphore failed
  Integer, Parameter, Public :: FIPC_sem_get_failed        = 15 + MPI_ERR_LASTCODE

  !> Return value indicating that an attempt to remove a semaphore failed
  Integer, Parameter, Public :: FIPC_sem_remove_failed     = 16 + MPI_ERR_LASTCODE

  !> Return value indicating that FreeIPC could not identify which semaphore is to be freed
  Integer, Parameter, Public :: FIPC_sem_not_found         = 17 + MPI_ERR_LASTCODE

  !> Return value indicating that an attempt to start a critical region failed
  Integer, Parameter, Public :: FIPC_critical_start_failed = 16 + MPI_ERR_LASTCODE

  !> Return value indicating that an attempt to start a critical region failed
  Integer, Parameter, Public :: FIPC_critical_end_failed   = 17 + MPI_ERR_LASTCODE

  !> Value to indicate that this process should not be in the context that results
  !> from a FIPC_ctxt_split
  Integer, Parameter, Public :: FIPC_undefined = MPI_undefined

  !> Value to indicate that an attempt has been made to extract a non-existant communicator
  !> from a context. The most common case is trying to extract the extra communicator on a process
  !> that is not rank zero in an intra communicator
  Integer, Parameter, Public :: FIPC_comm_null = MPI_comm_null

  ! Reductions available
  !> Handle to indicate a global sum is to be performed in a reduction operation
  Integer, Parameter, Public :: FIPC_sum  = mpi_sum
  !> Handle to indicate a global product is to be performed in a reduction operation
  Integer, Parameter, Public :: FIPC_prod = mpi_prod
  !> Handle to indicate a global minimum is to be performed in a reduction operation
  Integer, Parameter, Public :: FIPC_max  = mpi_max
  !> Handle to indicate a global maximum is to be performed in a reduction operation
  Integer, Parameter, Public :: FIPC_min  = mpi_min

  ! Public interfaces
  Public :: FIPC_init
  Public :: FIPC_initialized
  Public :: FIPC_finalize
  Public :: FIPC_ctxt_dup
  Public :: FIPC_ctxt_split
  Public :: FIPC_ctxt_free
  Public :: FIPC_ctxt_intra_comm
  Public :: FIPC_ctxt_world_comm
  Public :: FIPC_ctxt_extra_comm
  Public :: FIPC_seg_create
  Public :: FIPC_seg_free
  Public :: FIPC_node_barrier
  Public :: FIPC_critical_start
  Public :: FIPC_critical_end
!!$  Public :: FIPC_allreduce

  Private


  !> \brief Create a shared memory segment.
  !!
  !! Sets up a shared memory segment on each of the shared memory nodes spanned by the context ctxt.
  !! Subroutine FIPC_seg_create(
  !!                             Type( FIPC_ctxt ),intent(in)    ctxt,\n
  !!                             Integer, Dimension( : ),intent(in) n,\n
  !!                             <choice>,Dimension( <1-7d> ),intent(out) a,\n
  !!                             Integer,intent(out) error )
  !! \param ctxt  The context within which to create the segment
  !! \param n     An array containg the dimensions of the arry which will be stored in the segment
  !! \param a     A pointer to a 1-7 dimensionsal array of type Integer( c_int ), Real( c_double ) o
  !!              Complex( complex ). On exit this points to the shared memory segment
  !! \param error On success ERROR is set to FIPC_SUCCESS. Any other value indicates error.
  !>
  Interface FIPC_seg_create
     Module Procedure seg_create_integer_1d_size_in_int
     Module Procedure seg_create_integer_2d_size_in_int
     Module Procedure seg_create_integer_3d_size_in_int
     Module Procedure seg_create_integer_4d_size_in_int
     Module Procedure seg_create_integer_5d_size_in_int
     Module Procedure seg_create_integer_6d_size_in_int
     Module Procedure seg_create_integer_7d_size_in_int
     Module Procedure seg_create_double_1d_size_in_int
     Module Procedure seg_create_double_2d_size_in_int
     Module Procedure seg_create_double_3d_size_in_int
     Module Procedure seg_create_double_4d_size_in_int
     Module Procedure seg_create_double_5d_size_in_int
     Module Procedure seg_create_double_6d_size_in_int
     Module Procedure seg_create_double_7d_size_in_int
     Module Procedure seg_create_complex_1d_size_in_int
     Module Procedure seg_create_complex_2d_size_in_int
     Module Procedure seg_create_complex_3d_size_in_int
     Module Procedure seg_create_complex_4d_size_in_int
     Module Procedure seg_create_complex_5d_size_in_int
     Module Procedure seg_create_complex_6d_size_in_int
     Module Procedure seg_create_complex_7d_size_in_int
  End Interface

  !>
  !!
  !! \brief Frees a shared memory segment.
  !!
  !! This routine frees a shared memory segment on each of the shared memory nodes spanned by the co
  !>
  Interface FIPC_seg_free
     Module Procedure seg_free_integer_1d
     Module Procedure seg_free_integer_2d
     Module Procedure seg_free_integer_3d
     Module Procedure seg_free_integer_4d
     Module Procedure seg_free_integer_5d
     Module Procedure seg_free_integer_6d
     Module Procedure seg_free_integer_7d
     Module Procedure seg_free_double_1d
     Module Procedure seg_free_double_2d
     Module Procedure seg_free_double_3d
     Module Procedure seg_free_double_4d
     Module Procedure seg_free_double_5d
     Module Procedure seg_free_double_6d
     Module Procedure seg_free_double_7d
     Module Procedure seg_free_complex_1d
     Module Procedure seg_free_complex_2d
     Module Procedure seg_free_complex_3d
     Module Procedure seg_free_complex_4d
     Module Procedure seg_free_complex_5d
     Module Procedure seg_free_complex_6d
     Module Procedure seg_free_complex_7d
  End Interface

  !!!!> \cond for_fipc_implementors

  ! Error flags that can be returned by the System V routines
  ! Need to read these from the C header files
  Integer, Private :: EACCES
  Integer, Private :: EEXIST
  Integer, Private :: EINVAL
  Integer, Private :: ENFILE
  Integer, Private :: ENOENT
  Integer, Private :: ENOMEM
  Integer, Private :: ENOSPC

  ! Flags for control of the creation of shared beasties
  Integer( c_int ), Private :: SEG_NOCREATE   = 0
  Integer( c_int ), Private :: SEG_CREATE     = 1
  Integer( c_int ), Private :: SEG_NOEXCLUDE  = 0
  Integer( c_int ), Private :: SEG_EXCLUDE    = 1
  Integer( c_int ), Private :: SEG_UREAD      = 4 * 8 * 8
  Integer( c_int ), Private :: SEG_UWRITE     = 2 * 8 * 8
  Integer( c_int ), Private :: SEG_GREAD      = 4 * 8
  Integer( c_int ), Private :: SEG_GWRITE     = 2 * 8
  Integer( c_int ), Private :: SEG_WREAD      = 4
  Integer( c_int ), Private :: SEG_WWRITE     = 2
  Integer( c_int ), Private :: SEG_NOREADONLY = 0
  Integer( c_int ), Private :: SEG_READONLY   = 1

  ! Elements of seg inquire array
  Integer( c_long ), Parameter, Private :: SEG_SEGSZ  = 1 ! size of segment in bytes
  Integer( c_long ), Parameter, Private :: SEG_LPID   = 2 ! process ID of last shared memory operati
  Integer( c_long ), Parameter, Private :: SEG_CPID   = 3 ! process ID of creator
  Integer( c_long ), Parameter, Private :: SEG_NATTCH = 4 ! number of current attaches
  Integer( c_long ), Parameter, Private :: SEG_ATIME  = 5 ! time of last shmat()
  Integer( c_long ), Parameter, Private :: SEG_DTIME  = 6 ! time of last shmdt()
  Integer( c_long ), Parameter, Private :: SEG_CTIME  = 7 ! time of last change by shmctl()

  ! For differentiating between data types
  Integer, Parameter, Private :: integer_vals = 1
  Integer, Parameter, Private :: double_vals  = 2
  Integer, Parameter, Private :: complex_vals = 3

  ! Debugging Flag. Set to false for production
  Logical, Parameter, Private :: debug = .False.

  ! Derived type for storing data about a segment
  Type, Private :: segment
     Integer( c_int  )                              :: shmid
     Type( c_ptr     )                              :: shmaddr
     Type( FIPC_ctxt )                              :: ctxt
     Integer                                        :: type
     Integer( c_long ), Dimension( : ), Allocatable :: sizes
  End Type segment

  ! Derived type for linked list for saving data about created segments
  Type, Private :: segment_list_type
     Type( segment           )          :: data
     Type( segment_list_type ), Pointer :: next => Null()
  End Type segment_list_type

  ! Linked list of data about created segments
  Type( segment_list_type ), Pointer, Private :: seg_list => Null()

  ! Derived type for storing data about a semaphore
  Type, Private :: semaphore
     Integer( c_int  ) :: semid
     Type( FIPC_ctxt ) :: ctxt
  End Type semaphore

  ! Derived type for linked list for saving data about created semaphores
  Type, Private :: semaphore_list_type
     Integer( c_int )                     :: semid
     Integer                              :: intra_handle
     Type( semaphore_list_type ), Pointer :: next => Null()
  End Type semaphore_list_type

  ! Linked list of data about created segements
  Type( semaphore_list_type ), Pointer, Private :: sem_list => Null()

  ! Largest allowed size for temporary buffers
  Integer, Parameter, Private :: max_buff_size = 1024 * 1024 / 8 ! 1 Mbyte of reals

  ! Interfaces to C wrappers.
  Interface

     Subroutine FIPC_get_errval( EACCES, EEXIST, EINVAL, ENFILE, ENOENT, ENOMEM, ENOSPC ) Bind( C )
       Use, Intrinsic :: iso_c_binding, Only : c_int
       Implicit None
       Integer( c_int ), Intent( In ) :: EACCES
       Integer( c_int ), Intent( In ) :: EEXIST
       Integer( c_int ), Intent( In ) :: EINVAL
       Integer( c_int ), Intent( In ) :: ENFILE
       Integer( c_int ), Intent( In ) :: ENOENT
       Integer( c_int ), Intent( In ) :: ENOMEM
       Integer( c_int ), Intent( In ) :: ENOSPC
     End Subroutine FIPC_get_errval
