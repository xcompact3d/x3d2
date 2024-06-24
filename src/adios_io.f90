module m_adios_io
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
    procedure :: init 
    procedure :: deinit 
    procedure :: handle_fatal_error 
    procedure :: read 
    procedure :: write 
    procedure :: write_real
  end type

contains
  subroutine init(self, comm_in)
    class(adios_io_t), intent(inout) :: self
    integer, intent(in), optional :: comm_in

    integer :: comm
    integer :: ierr
    logical :: is_mpi_initialised

    if(present(comm_in)) then
      comm = comm_in
    else
      comm = MPI_COMM_WORLD
    endif

    call mpi_initialized(is_mpi_initialised, ierr)

    if(.not. is_mpi_initialised) then
      call self%handle_fatal_error("ADIOS must be initialised after MPI!", 0)
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

  pure function calc_local_strided_offset(offset, stride) result(local_offset)
    ! Returns the first index in the local
    ! subregion of a strided, distributed array.
    integer(8), intent(in) :: offset ! offset of subregion owned by this rank
    integer(8), intent(in) :: stride ! number of grid points between adjacent output indices
    integer(8) :: local_offset ! return value: offset into local region respecting global stride

    local_offset = mod(stride - mod(offset, stride), stride) + 1
  end function

  subroutine write_real(self, var, fpath, varname)
    class(adios_io_t), intent(inout) :: self
    real(8), intent(in) :: var
    character(*), intent(in) :: fpath !! Path to ouptut file
    character(*), intent(in) :: varname !! Name of variable in output file

    type(adios2_io) :: io
    type(adios2_engine) :: writer
    type(adios2_variable) :: adios_var
    integer :: vartype
    integer :: ierr

    call adios2_declare_io (io, self%adios_ctx, 'write', ierr)
    if (.not.io%valid) then
      call self%handle_fatal_error("Cannot create ADIOS2 IO", ierr)
    endif

    call adios2_open(writer, io, fpath, adios2_mode_write, ierr)
    if (.not.writer%valid) then
      call self%handle_fatal_error("Cannot create ADIOS2 writer", ierr)
    endif

    vartype = adios2_type_real

    call adios2_define_variable(adios_var, io, varname, vartype, ierr)
    call adios2_begin_step(writer, adios2_step_mode_append, ierr)
    call adios2_put(writer, adios_var, var, ierr)
    call adios2_end_step(writer, ierr)

    if (writer%valid) then
      call adios2_close(writer, ierr)
    endif

  end subroutine

  subroutine write(self, in_arr, fpath, varname, icount, ishape, istart, convert_to_sp_in, istride_in)
    class(adios_io_t), intent(inout) :: self
    class(field_t), pointer, intent(in) :: in_arr !! Field to be outputted
    ! TODO should this be a field or a fortran array?
    character(*), intent(in) :: fpath !! Path to ouptut file
    character(*), intent(in) :: varname !! Name of variable in output file
    integer(8), dimension(3), intent(in) :: icount !! Local size of in_arr
    integer(8), dimension(3), intent(in) :: ishape !! Global size of in_arr
    integer(8), dimension(3), intent(in) :: istart !! Local offset of in_arr

    !> if .true. input array will be converted to single precision before being
    !> outputted.
    !> Defaults to .false.
    logical, intent(in), optional :: convert_to_sp_in
    logical :: convert_to_sp = .false.
    logical :: should_stride = .false.
    !> If set, will coarsen output in each direction by only dumping every
    !> `istride` gridpoints.
    !> Defaults to (1,1,1).
    integer(8), dimension(3), intent(in), optional :: istride_in
    integer(8), dimension(3) :: istride = [1,1,1]

    real(4), allocatable :: data_sp(:,:,:)
    real(4), allocatable :: data_strided_sp(:,:,:)
    real(8), allocatable :: data_strided_dp(:,:,:)
    integer(8), dimension(3) :: istart_strided
    integer(8), dimension(3) :: icount_strided

    type(adios2_io) :: io
    type(adios2_engine) :: writer
    type(adios2_variable) :: adios_var
    integer :: vartype
    integer :: ierr

    integer :: i

    ! Set our optional inputs
    if(present(convert_to_sp_in)) then
      convert_to_sp = convert_to_sp_in
    else
      convert_to_sp = .false.
    endif
    if(present(istride_in)) istride = istride_in

    if(istride(1) < 1 .or. istride(2) < 1 .or. istride(3) < 1) then
      call self%handle_fatal_error("Output stride < 1. Cannot continue.", 0)
    endif

    call adios2_declare_io(io, self%adios_ctx, 'write', ierr)
    if (.not.io%valid) then
      call self%handle_fatal_error("Cannot create ADIOS2 IO", ierr)
    endif

    call adios2_open(writer, io, fpath, adios2_mode_write, ierr)
    if (.not.writer%valid) then
      call self%handle_fatal_error("Cannot create ADIOS2 writer", ierr)
    endif

    if(convert_to_sp) then
      print*, "converting to single precision!"
      allocate(data_sp(icount(1), icount(2), icount(3)))
      data_sp(:,:,:) = real(in_arr%data(:,:,:))
      vartype = adios2_type_real
    else
      vartype = adios2_type_dp
    endif

    should_stride = any(istride /= [1,1,1])

    if(should_stride) then
      print*, "striding!"
      do i = 1, 3
        istart_strided(i) = calc_local_strided_offset(istart(i), istride(i))
      end do
      icount_strided = icount/istride
      if(convert_to_sp) then
        allocate(data_strided_sp(&
          icount_strided(1),&
          icount_strided(2),&
          icount_strided(3)&
          ))
        data_strided_sp(:,:,:) = data_sp(&
          istart_strided(1)::istride(1),&
          istart_strided(2)::istride(2),&
          istart_strided(3)::istride(3)&
          )
      else
        allocate(data_strided_dp(&
          icount_strided(1),&
          icount_strided(2),&
          icount_strided(3)&
          ))
        data_strided_dp(:,:,:) = in_arr%data(&
          istart_strided(1)::istride(1),&
          istart_strided(2)::istride(2),&
          istart_strided(3)::istride(3)&
          )
      endif
    endif

    if(should_stride) then
      call adios2_define_variable(adios_var, io, varname, vartype, &
        3, ishape/istride, istart/istride, icount_strided, adios2_constant_dims, ierr)
    else
      call adios2_define_variable(adios_var, io, varname, vartype, &
        3, ishape, istart, icount, adios2_constant_dims, ierr)
    endif

    call adios2_begin_step(writer, adios2_step_mode_append, ierr)

    if(should_stride) then
      if(convert_to_sp) then
        call adios2_put(writer, adios_var, data_strided_sp, ierr)
      else
        call adios2_put(writer, adios_var, data_strided_dp, ierr)
      endif
    else
      if(convert_to_sp) then
        call adios2_put(writer, adios_var, data_sp, ierr)
      else
        call adios2_put(writer, adios_var, in_arr%data, ierr)
      endif
    endif

    call adios2_end_step(writer, ierr)

    if (writer%valid) then
      call adios2_close(writer, ierr)
    endif

    if(allocated(data_sp)) then
      deallocate(data_sp)
    endif

  end subroutine

  subroutine read(self, out_arr, fpath, varname, icount, ishape, istart, idump_in)
    class(adios_io_t), intent(inout) :: self
    class(field_t), pointer, intent(in) :: out_arr !! Field to be read from file
    character(*), intent(in) :: fpath !! Path to input file
    character(*), intent(in) :: varname !! Name of variable in input file
    integer(kind=8), dimension(3), intent(in) :: icount !! Local size of in_arr
    integer(kind=8), dimension(3), intent(in) :: ishape !! Global size of in_arr
    integer(kind=8), dimension(3), intent(in) :: istart !! Local offset of in_arr
    !> When multiple timesteps have been dumped into one variable, chose which
    !> dump to read.
    integer(kind=8), intent(in), optional :: idump_in
    integer(kind=8) :: idump = 0

    integer(8), parameter :: n_steps = 1 !! This version reads one dump at a time

    type(adios2_io) :: io
    type(adios2_engine) :: reader
    type(adios2_variable) :: adios_var
    integer :: ierr

    if (present(idump_in)) idump = idump_in

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
      call self%handle_fatal_error("Cannot fetch ADIOS2 variable", ierr)
    endif

    call adios2_set_step_selection(adios_var, idump, n_steps, ierr)
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
