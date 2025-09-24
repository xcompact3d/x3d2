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

contains

  function make_reader() result(reader)
    class(io_reader_t), allocatable :: reader
#ifdef WITH_ADIOS2
    allocate(io_adios2_reader_t :: reader)
#else
    allocate(io_dummy_reader_t :: reader)
#endif
  end function make_reader

  function make_writer() result(writer)
    class(io_writer_t), allocatable :: writer
#ifdef WITH_ADIOS2
    allocate(io_adios2_writer_t :: writer)
#else
    allocate(io_dummy_writer_t :: writer)
#endif
  end function make_writer

  subroutine make_reader_ptr(reader)
    class(io_reader_t), pointer, intent(out) :: reader
#ifdef WITH_ADIOS2
    allocate(io_adios2_reader_t :: reader)
#else
    allocate(io_dummy_reader_t :: reader)
#endif
  end subroutine make_reader_ptr

  subroutine make_writer_ptr(writer)
    class(io_writer_t), pointer, intent(out) :: writer
#ifdef WITH_ADIOS2
    allocate(io_adios2_writer_t :: writer)
#else
    allocate(io_dummy_writer_t :: writer)
#endif
  end subroutine make_writer_ptr

end module m_io_factory