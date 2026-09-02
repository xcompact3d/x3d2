program test_adios2_mixed_precision
  use adios2, only: adios2_adios, adios2_io, adios2_engine, &
                    adios2_variable, adios2_init, adios2_finalize, &
                    adios2_declare_io, adios2_set_engine, &
                    adios2_define_variable, adios2_open, adios2_close, &
                    adios2_begin_step, adios2_end_step, adios2_put, &
                    adios2_mode_write, adios2_mode_sync, &
                    adios2_step_mode_append, adios2_constant_dims, &
                    adios2_type_real, adios2_type_dp
  use iso_fortran_env, only: real64, stderr => error_unit
  use m_common, only: dp, i8, sp
  use m_io_backend, only: allocate_io_reader
  use m_io_base, only: io_reader_t, io_file_t, io_mode_read
  use mpi, only: MPI_COMM_SELF, MPI_Init, MPI_Finalize

  implicit none

  character(len=*), parameter :: filename = "test_mixed_precision_input.bp"
  integer, parameter :: nx = 3, ny = 2, nz = 2

  class(io_reader_t), allocatable :: reader
  class(io_file_t), allocatable :: file
  real(sp) :: scalar_sp
  real(real64) :: scalar_dp
  real(dp) :: scalar_read
  real(sp) :: array_sp(nx, ny, nz)
  real(real64) :: array_dp(nx, ny, nz)
  real(dp) :: array_read(nx - 1, ny, nz)
  integer(i8) :: selection_start(3), selection_count(3)
  integer :: i, j, k, ierr
  logical :: allpass

  call MPI_Init(ierr)

  scalar_sp = 1.25_sp
  scalar_dp = 1.0_real64/10.0_real64
  do k = 1, nz
    do j = 1, ny
      do i = 1, nx
        array_sp(i, j, k) = real(100*k + 10*j + i, sp)/4.0_sp
        array_dp(i, j, k) = real(100*k + 10*j + i, real64)/10.0_real64
      end do
    end do
  end do

  call write_mixed_precision_file(filename, scalar_sp, scalar_dp, &
                                  array_sp, array_dp)

  call allocate_io_reader(reader)
  call reader%init(MPI_COMM_SELF, "test_mixed_precision_read")
  file = reader%open(filename, io_mode_read, MPI_COMM_SELF)
  call file%begin_step()

  allpass = .true.
  selection_start = [1_i8, 0_i8, 0_i8]
  selection_count = int(shape(array_read), i8)

  call reader%read_data("scalar_sp", scalar_read, file)
  allpass = allpass .and. scalar_read == real(scalar_sp, dp)

  call reader%read_data("scalar_dp", scalar_read, file)
  allpass = allpass .and. scalar_read == real(scalar_dp, dp)

  call reader%read_data("array_sp", array_read, file, &
                        start_dims=selection_start, &
                        count_dims=selection_count)
  allpass = allpass .and. all(array_read == real(array_sp(2:nx, :, :), dp))

  call reader%read_data("array_dp", array_read, file, &
                        start_dims=selection_start, &
                        count_dims=selection_count)
  allpass = allpass .and. all(array_read == real(array_dp(2:nx, :, :), dp))

  call file%end_step()
  call file%close()
  call reader%finalise()
  deallocate (reader)

  call MPI_Finalize(ierr)

  if (allpass) then
    write (stderr, '(a)') 'ADIOS2 MIXED-PRECISION TEST PASSED SUCCESSFULLY.'
  else
    error stop 'ADIOS2 MIXED-PRECISION TEST FAILED.'
  end if

contains

  subroutine write_mixed_precision_file( &
    output_filename, scalar_sp, scalar_dp, array_sp, array_dp &
    )
    character(len=*), intent(in) :: output_filename
    real(sp), intent(in) :: scalar_sp
    real(real64), intent(in) :: scalar_dp
    real(sp), intent(in) :: array_sp(:, :, :)
    real(real64), intent(in) :: array_dp(:, :, :)

    type(adios2_adios) :: adios
    type(adios2_io) :: io
    type(adios2_engine) :: engine
    type(adios2_variable) :: var_scalar_sp, var_scalar_dp
    type(adios2_variable) :: var_array_sp, var_array_dp
    integer(i8) :: shape_dims(3), start_dims(3), count_dims(3)
    integer :: ierr

    shape_dims = int(shape(array_sp), i8)
    start_dims = 0_i8
    count_dims = shape_dims

    call adios2_init(adios, MPI_COMM_SELF, ierr)
    call check_adios2_error(ierr, "Failed to initialise writer")
    call adios2_declare_io(io, adios, "test_mixed_precision_write", ierr)
    call check_adios2_error(ierr, "Failed to declare writer IO")
    call adios2_set_engine(io, "BP5", ierr)
    call check_adios2_error(ierr, "Failed to set writer engine")

    call adios2_define_variable( &
      var_scalar_sp, io, "scalar_sp", adios2_type_real, ierr)
    call check_adios2_error(ierr, "Failed to define single-precision scalar")
    call adios2_define_variable( &
      var_scalar_dp, io, "scalar_dp", adios2_type_dp, ierr)
    call check_adios2_error(ierr, "Failed to define double-precision scalar")
    call adios2_define_variable( &
      var_array_sp, io, "array_sp", adios2_type_real, 3, shape_dims, &
      start_dims, count_dims, adios2_constant_dims, ierr)
    call check_adios2_error(ierr, "Failed to define single-precision array")
    call adios2_define_variable( &
      var_array_dp, io, "array_dp", adios2_type_dp, 3, shape_dims, &
      start_dims, count_dims, adios2_constant_dims, ierr)
    call check_adios2_error(ierr, "Failed to define double-precision array")

    call adios2_open( &
      engine, io, output_filename, adios2_mode_write, MPI_COMM_SELF, ierr)
    call check_adios2_error(ierr, "Failed to open writer")
    call adios2_begin_step(engine, adios2_step_mode_append, ierr)
    call check_adios2_error(ierr, "Failed to begin writer step")

    call adios2_put(engine, var_scalar_sp, scalar_sp, adios2_mode_sync, ierr)
    call check_adios2_error(ierr, "Failed to write single-precision scalar")
    call adios2_put(engine, var_scalar_dp, scalar_dp, adios2_mode_sync, ierr)
    call check_adios2_error(ierr, "Failed to write double-precision scalar")
    call adios2_put(engine, var_array_sp, array_sp, adios2_mode_sync, ierr)
    call check_adios2_error(ierr, "Failed to write single-precision array")
    call adios2_put(engine, var_array_dp, array_dp, adios2_mode_sync, ierr)
    call check_adios2_error(ierr, "Failed to write double-precision array")

    call adios2_end_step(engine, ierr)
    call check_adios2_error(ierr, "Failed to end writer step")
    call adios2_close(engine, ierr)
    call check_adios2_error(ierr, "Failed to close writer")
    call adios2_finalize(adios, ierr)
    call check_adios2_error(ierr, "Failed to finalise writer")
  end subroutine write_mixed_precision_file

  subroutine check_adios2_error(ierr, message)
    integer, intent(in) :: ierr
    character(len=*), intent(in) :: message

    if (ierr /= 0) error stop message
  end subroutine check_adios2_error

end program test_adios2_mixed_precision
