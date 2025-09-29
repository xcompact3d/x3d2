module m_io_factory
!! Factory module to handle allocation of I/O backends

  use m_io_base, only: io_reader_t, io_writer_t
  use m_io_dummy, only: io_dummy_reader_t, io_dummy_writer_t
#ifdef WITH_ADIOS2
  use m_io_adios2, only: io_adios2_reader_t, io_adios2_writer_t
#endif

  implicit none

  private
  public :: allocate_io_reader, allocate_io_writer
  public :: IO_BACKEND_DUMMY, IO_BACKEND_ADIOS2, IO_BACKEND_MPIIO
  public :: get_default_backend

  ! backend type constants
  integer, parameter :: IO_BACKEND_DUMMY = 0
  integer, parameter :: IO_BACKEND_ADIOS2 = 1
  integer, parameter :: IO_BACKEND_MPIIO = 2

contains

  !> Get the default I/O backend based on compile-time configuration
  function get_default_backend() result(backend)
    integer :: backend
#ifdef WITH_ADIOS2
    backend = IO_BACKEND_ADIOS2
#else
    backend = IO_BACKEND_DUMMY
#endif
  end function get_default_backend

  !> Allocation interface for readers
  subroutine allocate_io_reader(reader)
    class(io_reader_t), pointer, intent(out) :: reader

    integer :: backend

    backend = get_default_backend()

    select case (backend)
#ifdef WITH_ADIOS2
    case (IO_BACKEND_ADIOS2)
      allocate (io_adios2_reader_t :: reader)
#endif
    case (IO_BACKEND_MPIIO)
      ! TODO: implement MPI-IO backend
      ! For now, fall back to dummy
      allocate (io_dummy_reader_t :: reader)
    case default
      allocate (io_dummy_reader_t :: reader)
    end select
  end subroutine allocate_io_reader

  !> Allocation interface for writers
  subroutine allocate_io_writer(writer)
    class(io_writer_t), pointer, intent(out) :: writer

    integer :: backend

    backend = get_default_backend()

    select case (backend)
#ifdef WITH_ADIOS2
    case (IO_BACKEND_ADIOS2)
      allocate (io_adios2_writer_t :: writer)
#endif
    case (IO_BACKEND_MPIIO)
      ! TODO: implement MPI-IO backend
      ! For now, fall back to dummy
      allocate (io_dummy_writer_t :: writer)
    case default
      allocate (io_dummy_writer_t :: writer)
    end select
  end subroutine allocate_io_writer

end module m_io_factory