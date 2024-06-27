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
    type(adios2_io) :: aio
    integer :: irank, nproc
    integer :: io_mode_write = adios2_mode_write
    integer :: io_mode_read = adios2_mode_read
  contains
    procedure :: init 
    procedure :: deinit 
    procedure :: handle_fatal_error 
    procedure :: read_field
    procedure :: write_field
    procedure :: write_real8
    procedure :: write_integer8
    procedure :: open_file
    procedure :: close_file
  end type

  type :: adios_file_t
    type(adios2_engine) :: engine
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
      call self%handle_fatal_error("IO: ADIOS must be initialised after MPI!", 0)
    endif

    call MPI_Comm_rank(comm, self%irank, ierr)

    call adios2_init(self%adios_ctx, comm, ierr)
    if (.not.self%adios_ctx%valid) then
      call self%handle_fatal_error("IO: Cannot initialise ADIOS context.", ierr)
    endif

    call adios2_declare_io (self%aio, self%adios_ctx, "main_io", ierr)
    if (.not.self%aio%valid) then
      call self%handle_fatal_error("IO: Cannot create ADIOS2 IO", ierr)
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

  function open_file(self, fpath, mode) result(file)
    class(adios_io_t), intent(inout) :: self
    character(*), intent(in) :: fpath !! Path to ouptut file
    integer :: mode !! Read/write mode. Should be adios2_mode_write or adios2_mode_read
    type(adios_file_t) :: file

    integer :: ierr

    call adios2_open(file%engine, self%aio, fpath, mode, ierr)
    if (.not.file%engine%valid) then
      call self%handle_fatal_error("IO: Cannot create ADIOS2 engine", ierr)
    endif

    if(mode == adios2_mode_write) then
      call adios2_begin_step(file%engine, adios2_step_mode_append, ierr)
    else if(mode == adios2_mode_read) then
      call adios2_begin_step(file%engine, adios2_step_mode_read, ierr)
    else
      error stop "Unsupported mode"
    endif

    if (ierr /= 0) then
      call self%handle_fatal_error("IO: Cannot begin step", ierr)
    endif

  end function

  subroutine close_file(self, file)
    class(adios_io_t), intent(inout) :: self
    class(adios_file_t), intent(inout) :: file

    integer :: ierr
    ! Might be abusing ADIOS here... We tend to redefined variables every time we open a file
    call adios2_remove_all_variables(self%aio, ierr) 
    call adios2_end_step(file%engine, ierr)

    if (file%engine%valid) then
      call adios2_close(file%engine, ierr)
    endif
  end subroutine

  subroutine write_real8(self, var, file, varname)
    class(adios_io_t), intent(inout) :: self
    real(8), intent(in) :: var
    type(adios_file_t), intent(inout) :: file
    character(*), intent(in) :: varname !! Name of variable in output file

    type(adios2_variable) :: adios_var
    integer :: vartype

    integer :: ierr

    vartype = adios2_type_real8

    call adios2_define_variable(adios_var, self%aio, varname, vartype, ierr)
    call adios2_put(file%engine, adios_var, var, ierr)
  end subroutine

  subroutine write_integer8(self, var, file, varname)
    class(adios_io_t), intent(inout) :: self
    integer(8), intent(in) :: var
    type(adios_file_t), intent(inout) :: file
    character(*), intent(in) :: varname !! Name of variable in output file

    type(adios2_variable) :: adios_var
    integer :: vartype

    integer :: ierr

    vartype = adios2_type_integer8

    call adios2_define_variable(adios_var, self%aio, varname, vartype, ierr)
    call adios2_put(file%engine, adios_var, var, ierr)
  end subroutine

  subroutine write_field(self, in_arr, file, varname, icount, ishape, istart, convert_to_sp_in, istride_in)
    class(adios_io_t), intent(inout) :: self
    class(field_t), pointer, intent(in) :: in_arr !! Field to be outputted
    ! TODO should this be a field or a fortran array?
    type(adios_file_t), intent(inout) :: file
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
      call adios2_define_variable(adios_var, self%aio, varname, vartype, &
        3, ishape/istride, istart/istride, icount_strided, adios2_constant_dims, ierr)
    else
      call adios2_define_variable(adios_var, self%aio, varname, vartype, &
        3, ishape, istart, icount, adios2_constant_dims, ierr)
    endif

    if(should_stride) then
      if(convert_to_sp) then
        call adios2_put(file%engine, adios_var, data_strided_sp, ierr)
      else
        call adios2_put(file%engine, adios_var, data_strided_dp, ierr)
      endif
    else
      if(convert_to_sp) then
        call adios2_put(file%engine, adios_var, data_sp, ierr)
      else
        call adios2_put(file%engine, adios_var, in_arr%data, ierr)
      endif
    endif

    if(allocated(data_sp)) then
      deallocate(data_sp)
    endif

  end subroutine

  subroutine read_field(self, out_arr, file, varname, icount, ishape, istart, idump_in)
    class(adios_io_t), intent(inout) :: self
    class(field_t), pointer, intent(in) :: out_arr !! Field to be read from file
    type(adios_file_t), intent(inout) :: file !! File already opened for reading
    character(*), intent(in) :: varname !! Name of variable in input file
    integer(kind=8), dimension(3), intent(in) :: icount !! Local size of in_arr
    integer(kind=8), dimension(3), intent(in) :: ishape !! Global size of in_arr
    integer(kind=8), dimension(3), intent(in) :: istart !! Local offset of in_arr
    !> When multiple timesteps have been dumped into one variable, chose which
    !> dump to read.
    integer(kind=8), intent(in), optional :: idump_in
    integer(kind=8) :: idump = 0

    integer(8), parameter :: n_steps = 1 !! This version reads one dump at a time

    type(adios2_variable) :: adios_var
    integer :: ierr

    if (present(idump_in)) idump = idump_in

    call adios2_inquire_variable(adios_var, self%aio, varname, ierr)
    if (.not.adios_var%valid) then
      call self%handle_fatal_error("Cannot fetch ADIOS2 variable", ierr)
    endif

    call adios2_set_step_selection(adios_var, idump, n_steps, ierr)
    if (ierr /= 0) then
      call self%handle_fatal_error("IO: Cannot set step selection", ierr)
    endif

    call adios2_set_selection(adios_var, 3, istart, icount, ierr)
    if (ierr /= 0) then
      call self%handle_fatal_error("IO: Cannot set selection", ierr)
    endif

    call adios2_get(file%engine, adios_var, out_arr%data, ierr)
    if (ierr /= 0) then
      call self%handle_fatal_error("IO: Cannot fetch data from file", ierr)
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
