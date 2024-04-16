module adios_io
  use adios2
  use mpi
  use m_allocator, only: field_t

  use, intrinsic :: iso_fortran_env, only : stdin=>input_unit, &
                                            stdout=>output_unit, &
                                            stderr=>error_unit

  implicit none

  type :: adios_io_t
    type(adios2_adios) :: adios_ctx
    integer :: irank, nproc
  contains
    procedure:: init 
    procedure:: deinit 
    procedure:: handle_fatal_error 
    procedure:: read 
    procedure:: write 
  end type

contains
  subroutine init(self, comm_in)
  class(adios_io_t), intent(inout) :: self
    integer, intent(in), optional :: comm_in

    ! TODO pass in communicator?
    ! TODO include check that MPI has been initialised
    ! TODO pass in MPI domain decomp data?

    integer :: comm
    integer :: ierr

    if(present(comm_in)) then
      comm = comm_in
    else
      comm = MPI_COMM_WORLD
    endif

    call MPI_Comm_rank(comm, self%irank, ierr)

    call adios2_init(self%adios_ctx, comm, ierr)
    if (.not.self%adios_ctx%valid) then
      call self%handle_fatal_error("Cannot initialise ADIOS context.", ierr)
    endif
  end subroutine

  subroutine deinit(self)
  class(adios_io_t), intent(inout) :: self
    integer :: ierr

    if (self%adios_ctx%valid) then
      call adios2_finalize(self%adios_ctx, ierr)
    endif
  end subroutine

  subroutine write(self, in_arr, fpath, varname, icount, ishape, istart)
  class(adios_io_t), intent(inout) :: self
  class(field_t), pointer, intent(in) :: in_arr !! Field to be outputted
    ! TODO should this be a field or a fortran array?
    character(*), intent(in) :: fpath !! Path to ouptut file
    character(*), intent(in) :: varname !! Name of variable in output file
    integer(kind=8), dimension(3), intent(in) :: icount !! Local size of in_arr
    integer(kind=8), dimension(3), intent(in) :: ishape !! Global size of in_arr
    integer(kind=8), dimension(3), intent(in) :: istart !! Local offset of in_arr

    type(adios2_io) :: io
    type(adios2_engine) :: writer
    type(adios2_variable) :: adios_var
    integer :: ierr

    call adios2_declare_io (io, self%adios_ctx, 'write', ierr)
    if (.not.io%valid) then
      call self%handle_fatal_error("Cannot create ADIOS2 IO", ierr)
    endif

    call adios2_open(writer, io, fpath, adios2_mode_write, ierr)
    if (.not.writer%valid) then
      call self%handle_fatal_error("Cannot create ADIOS2 writer", ierr)
    endif

    ! TODO tie double precision variable to array
    call adios2_define_variable(adios_var, io, varname, adios2_type_dp, &
      3, ishape, istart, icount, adios2_constant_dims, ierr)

    call adios2_begin_step(writer, adios2_step_mode_append, ierr)
    call adios2_put(writer, adios_var, in_arr%data, ierr)
    call adios2_end_step(writer, ierr)

    if (writer%valid) then
      call adios2_close(writer, ierr)
    endif

  end subroutine

  subroutine read(self, out_arr, fpath, varname, icount, ishape, istart)
  class(adios_io_t), intent(inout) :: self
  class(field_t), pointer, intent(in) :: out_arr !! Field to be read from file
    character(*), intent(in) :: fpath !! Path to input file
    character(*), intent(in) :: varname !! Name of variable in input file
    integer(kind=8), dimension(3), intent(in) :: icount !! Local size of in_arr
    integer(kind=8), dimension(3), intent(in) :: ishape !! Global size of in_arr
    integer(kind=8), dimension(3), intent(in) :: istart !! Local offset of in_arr

    integer(8), parameter :: initial_step = 0
    integer(8), parameter :: n_steps = 1

    type(adios2_io) :: io
    type(adios2_engine) :: reader
    type(adios2_variable) :: adios_var
    integer :: ierr

    call adios2_declare_io (io, self%adios_ctx, 'read', ierr)
    if (.not.io%valid) then
      call self%handle_fatal_error("Cannot create ADIOS2 IO", ierr)
    endif

    call adios2_open(reader, io, fpath, adios2_mode_read, MPI_COMM_SELF, ierr)
    if (.not.reader%valid) then
      call self%handle_fatal_error("Cannot create ADIOS2 reader", ierr)
    endif

    call adios2_begin_step(reader, adios2_step_mode_read, ierr)
    call adios2_inquire_variable(adios_var, io, varname, ierr)
    if (.not.adios_var%valid) then
      call self%handle_fatal_error("Cannot fetch ADIOS2 IO", ierr)
    endif

    call adios2_set_step_selection(adios_var, initial_step, n_steps, ierr)
    call adios2_set_selection(adios_var, 3, istart, icount, ierr)
    call adios2_get(reader, adios_var, out_arr%data, ierr)
    call adios2_end_step(reader, ierr)

    if (reader%valid) then
      call adios2_close(reader, ierr)
    endif
  end subroutine

  subroutine handle_fatal_error(self, msg, ierr)
  class(adios_io_t), intent(in) :: self
    integer, intent(in) :: ierr
    character(*), intent(in) :: msg

    write(stderr, *) "ADIOS2 Error code:", ierr
    write(stderr, *) msg
    error stop -1
  end subroutine
end module
