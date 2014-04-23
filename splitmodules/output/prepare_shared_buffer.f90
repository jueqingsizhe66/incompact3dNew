  subroutine prepare_shared_buffer(C,MPI_COMM,decomp)

  implicit none

  TYPE(SMP_INFO) :: C
  INTEGER :: MPI_COMM
  TYPE(DECOMP_INFO) :: decomp

  INTEGER, ALLOCATABLE :: KTBL(:,:),NARY(:,:),KTBLALL(:,:)
  INTEGER MYSMP, MYCORE, COLOR

  integer :: ierror

  C%MPI_COMM = MPI_COMM
  CALL MPI_COMM_SIZE(MPI_COMM,C%NCPU,ierror)
  CALL MPI_COMM_RANK(MPI_COMM,C%NODE_ME,ierror)
  C%SMP_COMM  = MPI_COMM_NULL
  C%CORE_COMM = MPI_COMM_NULL
  C%SMP_ME= 0
  C%NCORE = 0
  C%CORE_ME = 0
  C%MAXCORE = 0
  C%NSMP  = 0
  C%N_SND = 0
  C%N_RCV = 0
  C%SND_P = 0
  C%RCV_P = 0
  C%SND_P_c = 0
  C%RCV_P_c = 0

  ! get the smp-node map for this communicator and set up smp communicators
  CALL GET_SMP_MAP(C%MPI_COMM, C%NSMP, MYSMP, C%NCORE, MYCORE, C%MAXCORE)
  C%SMP_ME = MYSMP + 1
  C%CORE_ME = MYCORE + 1
  ! - set up inter/intra smp-node communicators
  COLOR = MYCORE
  IF (COLOR.GT.0) COLOR = MPI_UNDEFINED
  CALL MPI_Comm_split(C%MPI_COMM, COLOR, MYSMP, C%SMP_COMM, ierror)
  CALL MPI_Comm_split(C%MPI_COMM, MYSMP, MYCORE, C%CORE_COMM, ierror)
  ! - allocate work space
  ALLOCATE(KTBL(C%MAXCORE,C%NSMP),NARY(C%NCPU,C%NCORE))
  ALLOCATE(KTBLALL(C%MAXCORE,C%NSMP))
  ! - set up smp-node/core to node_me lookup table
  KTBL = 0
  KTBL(C%CORE_ME,C%SMP_ME) = C%NODE_ME + 1
  CALL MPI_ALLREDUCE(KTBL,KTBLALL,C%NSMP*C%MAXCORE,MPI_INTEGER, &
       MPI_SUM,MPI_COMM,ierror)
  KTBL=KTBLALL
  !  IF (SUM(KTBL) /= C%NCPU*(C%NCPU+1)/2) &
  !       CALL MPI_ABORT(...

  ! compute offsets in shared SNDBUF and RCVBUF
  CALL MAPSET_SMPSHM(C, KTBL, NARY, decomp)

  DEALLOCATE(KTBL,NARY)

  return
  end subroutine prepare_shared_buffer
