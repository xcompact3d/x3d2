module m_io_factory
!! Clean factory module that handles conditional compilation internally

  use m_io_base, only: io_reader_t, io_writer_t
  use m_io_dummy, only: io_dummy_reader_t, io_dummy_writer_t
#ifdef WITH_ADIOS2
  use m_io_adios2, only: io_adios2_reader_t, io_adios2_writer_t
#endif

  implicit none

  private
  public :: make_reader, make_writer, make_reader_ptr, make_writer_ptr
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

  function make_reader(backend) result(reader)
    class(io_reader_t), allocatable :: reader
    integer, intent(in), optional :: backend
    
    integer :: use_backend
    
    use_backend = get_default_backend()
    if (present(backend)) use_backend = backend
    
    select case(use_backend)
#ifdef WITH_ADIOS2
    case(IO_BACKEND_ADIOS2)
      allocate(io_adios2_reader_t :: reader)
#endif
    case(IO_BACKEND_MPIIO)
      ! TODO: allocate(io_mpiio_reader_t :: reader)
      ! for now, fall back to dummy
      allocate(io_dummy_reader_t :: reader)
    case default
      allocate(io_dummy_reader_t :: reader)
    end select
  end function make_reader

  function make_writer(backend) result(writer)
    class(io_writer_t), allocatable :: writer
    integer, intent(in), optional :: backend
    
    integer :: use_backend
    
    use_backend = get_default_backend()
    if (present(backend)) use_backend = backend
    
    select case(use_backend)
#ifdef WITH_ADIOS2
    case(IO_BACKEND_ADIOS2)
      allocate(io_adios2_writer_t :: writer)
#endif
    case(IO_BACKEND_MPIIO)
      ! TODO: allocate(io_mpiio_writer_t :: writer)
      ! for now, fall back to dummy
      allocate(io_dummy_writer_t :: writer)
    case default
      allocate(io_dummy_writer_t :: writer)
    end select
  end function make_writer

  subroutine make_reader_ptr(reader, backend)
    class(io_reader_t), pointer, intent(out) :: reader
    integer, intent(in), optional :: backend
    
    integer :: use_backend
    
    use_backend = get_default_backend()
    if (present(backend)) use_backend = backend
    
    select case(use_backend)
#ifdef WITH_ADIOS2
    case(IO_BACKEND_ADIOS2)
      allocate(io_adios2_reader_t :: reader)
#endif
    case(IO_BACKEND_MPIIO)
      ! TODO: allocate(io_mpiio_reader_t :: reader)
      ! for now, fall back to dummy
      allocate(io_dummy_reader_t :: reader)
    case default
      allocate(io_dummy_reader_t :: reader)
    end select
  end subroutine make_reader_ptr

  subroutine make_writer_ptr(writer, backend)
    class(io_writer_t), pointer, intent(out) :: writer
    integer, intent(in), optional :: backend
    
    integer :: use_backend
    
    use_backend = get_default_backend()
    if (present(backend)) use_backend = backend
    
    select case(use_backend)
#ifdef WITH_ADIOS2
    case(IO_BACKEND_ADIOS2)
      allocate(io_adios2_writer_t :: writer)
#endif
    case(IO_BACKEND_MPIIO)
      ! TODO: allocate(io_mpiio_writer_t :: writer)
      ! for now, fall back to dummy
      allocate(io_dummy_writer_t :: writer)
    case default
      allocate(io_dummy_writer_t :: writer)
    end select
  end subroutine make_writer_ptr

end module m_io_factory