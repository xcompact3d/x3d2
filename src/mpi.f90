module m_mpi
!! Single point of contact between x3d2 and MPI.
!!
!! When the build enables MPI (-DMPI, set by the WITH_MPI option, on by
!! default) this module is a thin re-export of the `mpi` module, so `use m_mpi`
!! is exactly `use mpi`.
!!
!! When MPI is disabled the module supplies serial stand-ins for the MPI
!! entities x3d2 uses. A build without MPI is a single rank, so every
!! collective reduces to the identity, ranks are 0 and communicator sizes are
!! 1. Call sites therefore stay free of preprocessor branches: they call
!! MPI_Allreduce, MPI_Comm_rank and friends unconditionally and get the right
!! serial answer. Point-to-point communication is the exception: it is only
!! reachable on more than one rank, so those stubs stop the run instead.
!!
!! Add a stub here rather than an `#ifdef MPI` at a call site when new MPI
!! calls appear. The only guards left in the tree are around communication on
!! device buffers, which cannot be typed generically here.
#ifdef MPI
  use mpi

  implicit none
#else
  use iso_fortran_env, only: int64, real32, real64, stderr => error_unit

  implicit none

  private
  public :: MPI_COMM_WORLD, MPI_COMM_SELF, MPI_COMM_NULL, MPI_SUCCESS
  public :: MPI_SUM, MPI_MAX, MPI_MIN, MPI_LAND
  public :: MPI_INTEGER, MPI_LOGICAL, MPI_REAL, MPI_DOUBLE_PRECISION, &
            MPI_COMPLEX, MPI_DOUBLE_COMPLEX
  public :: MPI_STATUS_SIZE, MPI_STATUS_IGNORE, MPI_STATUSES_IGNORE
  public :: mpi_in_place_t, MPI_IN_PLACE
  public :: MPI_Init, MPI_Finalize, MPI_Initialized, MPI_Abort
  public :: MPI_Comm_rank, MPI_Comm_size, MPI_Barrier, MPI_Wtime
  public :: MPI_Bcast, MPI_Allreduce, MPI_Reduce
  public :: MPI_Isend, MPI_Irecv, MPI_Waitall, MPI_Sendrecv

  !> Handles are opaque in MPI, so the values only have to be distinct.
  integer, parameter :: MPI_COMM_NULL = 0
  integer, parameter :: MPI_COMM_WORLD = 1
  integer, parameter :: MPI_COMM_SELF = 2

  integer, parameter :: MPI_SUCCESS = 0

  integer, parameter :: MPI_SUM = 1
  integer, parameter :: MPI_MAX = 2
  integer, parameter :: MPI_MIN = 3
  integer, parameter :: MPI_LAND = 4

  integer, parameter :: MPI_INTEGER = 1
  integer, parameter :: MPI_LOGICAL = 2
  integer, parameter :: MPI_REAL = 3
  integer, parameter :: MPI_DOUBLE_PRECISION = 4
  integer, parameter :: MPI_COMPLEX = 5
  integer, parameter :: MPI_DOUBLE_COMPLEX = 6

  integer, parameter :: MPI_STATUS_SIZE = 6
  integer, parameter :: MPI_STATUS_IGNORE(MPI_STATUS_SIZE) = 0
  integer, parameter :: MPI_STATUSES_IGNORE(MPI_STATUS_SIZE, 1) = 0

  !> MPI_IN_PLACE is a sentinel, not data. Giving it a type of its own is what
  !> lets the generic MPI_Allreduce below tell an in-place call (a no-op here)
  !> from one with separate send and receive buffers (a copy here).
  type :: mpi_in_place_t
    private
    integer :: sentinel = 0
  end type mpi_in_place_t

  type(mpi_in_place_t), parameter :: MPI_IN_PLACE = mpi_in_place_t(0)

  interface MPI_Allreduce
    module procedure allreduce_in_place
    module procedure allreduce_r32
    module procedure allreduce_r64
    module procedure allreduce_r64_1d
  end interface MPI_Allreduce

  interface MPI_Reduce
    module procedure reduce_r32
    module procedure reduce_r64
  end interface MPI_Reduce

contains

  subroutine MPI_Init(ierror)
    integer, intent(out) :: ierror

    ierror = MPI_SUCCESS
  end subroutine MPI_Init

  subroutine MPI_Finalize(ierror)
    integer, intent(out) :: ierror

    ierror = MPI_SUCCESS
  end subroutine MPI_Finalize

  subroutine MPI_Initialized(flag, ierror)
    !! Reported as initialised whether or not MPI_Init has been called. There
    !! is no MPI environment to bring up here, so callers asking whether it is
    !! safe to proceed always have their answer, and a library handed a
    !! communicator never has to wait for one.
    logical, intent(out) :: flag
    integer, intent(out) :: ierror

    flag = .true.
    ierror = MPI_SUCCESS
  end subroutine MPI_Initialized

  subroutine MPI_Abort(comm, errorcode, ierror)
    integer, intent(in) :: comm, errorcode
    integer, intent(out) :: ierror

    ierror = MPI_SUCCESS
    write (stderr, '(a,i0)') 'MPI_Abort with error code ', errorcode
    flush (stderr)
    error stop 'aborted'
  end subroutine MPI_Abort

  subroutine MPI_Comm_rank(comm, rank, ierror)
    integer, intent(in) :: comm
    integer, intent(out) :: rank, ierror

    rank = 0
    ierror = MPI_SUCCESS
  end subroutine MPI_Comm_rank

  subroutine MPI_Comm_size(comm, size, ierror)
    integer, intent(in) :: comm
    integer, intent(out) :: size, ierror

    size = 1
    ierror = MPI_SUCCESS
  end subroutine MPI_Comm_size

  subroutine MPI_Barrier(comm, ierror)
    integer, intent(in) :: comm
    integer, intent(out) :: ierror

    ierror = MPI_SUCCESS
  end subroutine MPI_Barrier

  function MPI_Wtime() result(t)
    !! Elapsed wall-clock seconds, from an arbitrary origin as MPI defines it.
    real(real64) :: t

    integer(int64) :: ticks, ticks_per_second

    call system_clock(ticks, ticks_per_second)
    t = real(ticks, real64)/real(ticks_per_second, real64)
  end function MPI_Wtime

  subroutine MPI_Bcast(buffer, count, datatype, root, comm, ierror)
    !! The only rank is the root, so its buffer is already the broadcast value.
    type(*), dimension(..), intent(inout) :: buffer
    integer, intent(in) :: count, datatype, root, comm
    integer, intent(out) :: ierror

    ierror = MPI_SUCCESS
  end subroutine MPI_Bcast

  subroutine allreduce_in_place(sendbuf, recvbuf, count, datatype, op, comm, &
                                ierror)
    !! Reducing one rank's contribution in place leaves the buffer as it is,
    !! whatever its type, shape or the operation asked for.
    type(mpi_in_place_t), intent(in) :: sendbuf
    type(*), dimension(..), intent(inout) :: recvbuf
    integer, intent(in) :: count, datatype, op, comm
    integer, intent(out) :: ierror

    ierror = MPI_SUCCESS
  end subroutine allreduce_in_place

  subroutine allreduce_r32(sendbuf, recvbuf, count, datatype, op, comm, ierror)
    real(real32), intent(in) :: sendbuf
    real(real32), intent(out) :: recvbuf
    integer, intent(in) :: count, datatype, op, comm
    integer, intent(out) :: ierror

    recvbuf = sendbuf
    ierror = MPI_SUCCESS
  end subroutine allreduce_r32

  subroutine allreduce_r64(sendbuf, recvbuf, count, datatype, op, comm, ierror)
    real(real64), intent(in) :: sendbuf
    real(real64), intent(out) :: recvbuf
    integer, intent(in) :: count, datatype, op, comm
    integer, intent(out) :: ierror

    recvbuf = sendbuf
    ierror = MPI_SUCCESS
  end subroutine allreduce_r64

  subroutine allreduce_r64_1d(sendbuf, recvbuf, count, datatype, op, comm, &
                              ierror)
    real(real64), intent(in) :: sendbuf(:)
    real(real64), intent(out) :: recvbuf(:)
    integer, intent(in) :: count, datatype, op, comm
    integer, intent(out) :: ierror

    recvbuf(1:count) = sendbuf(1:count)
    ierror = MPI_SUCCESS
  end subroutine allreduce_r64_1d

  subroutine reduce_r32(sendbuf, recvbuf, count, datatype, op, root, comm, &
                        ierror)
    real(real32), intent(in) :: sendbuf
    real(real32), intent(out) :: recvbuf
    integer, intent(in) :: count, datatype, op, root, comm
    integer, intent(out) :: ierror

    recvbuf = sendbuf
    ierror = MPI_SUCCESS
  end subroutine reduce_r32

  subroutine reduce_r64(sendbuf, recvbuf, count, datatype, op, root, comm, &
                        ierror)
    real(real64), intent(in) :: sendbuf
    real(real64), intent(out) :: recvbuf
    integer, intent(in) :: count, datatype, op, root, comm
    integer, intent(out) :: ierror

    recvbuf = sendbuf
    ierror = MPI_SUCCESS
  end subroutine reduce_r64

  ! Point-to-point communication has no serial meaning: a single rank has no
  ! partner to exchange with. Call sites reach these only after branching on a
  ! rank count above one, so the stubs exist to satisfy the compiler and stop
  ! the run if that branching is ever wrong.
  subroutine MPI_Isend(buf, count, datatype, dest, tag, comm, request, ierror)
    type(*), dimension(..), intent(in) :: buf
    integer, intent(in) :: count, datatype, dest, tag, comm
    integer, intent(out) :: request, ierror

    error stop 'MPI_Isend called, but this build was configured without MPI'
  end subroutine MPI_Isend

  subroutine MPI_Irecv(buf, count, datatype, source, tag, comm, request, &
                       ierror)
    type(*), dimension(..), intent(inout) :: buf
    integer, intent(in) :: count, datatype, source, tag, comm
    integer, intent(out) :: request, ierror

    error stop 'MPI_Irecv called, but this build was configured without MPI'
  end subroutine MPI_Irecv

  subroutine MPI_Waitall(count, array_of_requests, array_of_statuses, ierror)
    integer, intent(in) :: count
    integer, intent(inout) :: array_of_requests(*)
    integer, intent(in) :: array_of_statuses(*)
    integer, intent(out) :: ierror

    error stop 'MPI_Waitall called, but this build was configured without MPI'
  end subroutine MPI_Waitall

  subroutine MPI_Sendrecv(sendbuf, sendcount, sendtype, dest, sendtag, &
                          recvbuf, recvcount, recvtype, source, recvtag, &
                          comm, status, ierror)
    !! Reachable only where a call site has already branched on a rank count
    !! above one, as the exchanges in the CUDA Poisson solver do. It exists so
    !! that callers can name it in a `use m_mpi, only:` list that has to
    !! compile in both builds.
    type(*), dimension(..), intent(in) :: sendbuf
    type(*), dimension(..), intent(inout) :: recvbuf
    integer, intent(in) :: sendcount, sendtype, dest, sendtag
    integer, intent(in) :: recvcount, recvtype, source, recvtag, comm
    integer, intent(in) :: status(*)
    integer, intent(out) :: ierror

    error stop 'MPI_Sendrecv called, but this build was configured without MPI'
  end subroutine MPI_Sendrecv
#endif
end module m_mpi
