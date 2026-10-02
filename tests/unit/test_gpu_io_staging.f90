program test_gpu_io_staging
  !! Check gpu_io_staging_bytes against the device I/O staging semantics
  !! of the GPU-aware ADIOS2 write path.

  use m_common, only: i8, nbytes, is_sp
  use m_config, only: checkpoint_config_t
  use m_memory_estimate, only: gpu_io_staging_bytes

  implicit none

  integer, parameter :: dims(3) = [17, 33, 8]
  integer(i8), parameter :: n = 4488_i8 !! product(dims)
  logical :: test_pass

  test_pass = .true.

  call test_host_write_is_zero()
  call test_unit_stride_snapshot()
  call test_strided_snapshot_no_checkpoint()
  call test_checkpoint_strided()
  call test_snapshot_sp_extra()
  call test_no_padding_at_large_dims()
  call test_neither_enabled_is_zero()

  if (.not. test_pass) then
    error stop "FAIL"
  end if

contains

  subroutine test_host_write_is_zero()
    !! device_write=.false. must always give zero, even with a unit-stride
    !! snapshot enabled.
    type(checkpoint_config_t) :: cfg
    integer(i8) :: bytes

    cfg%snapshot_freq = 1
    cfg%output_stride = [1, 1, 1]
    cfg%checkpoint_freq = 0

    bytes = gpu_io_staging_bytes(dims, cfg, .false.)
    call check("host write is zero", bytes, 0_i8)
  end subroutine test_host_write_is_zero

  subroutine test_unit_stride_snapshot()
    !! Unit-stride snapshot on the device path costs one field.
    type(checkpoint_config_t) :: cfg
    integer(i8) :: bytes

    cfg%snapshot_freq = 1
    cfg%output_stride = [1, 1, 1]
    cfg%checkpoint_freq = 0
    cfg%snapshot_sp = .false.

    bytes = gpu_io_staging_bytes(dims, cfg, .true.)
    call check("unit-stride snapshot", bytes, n*int(nbytes, i8))
  end subroutine test_unit_stride_snapshot

  subroutine test_strided_snapshot_no_checkpoint()
    !! A strided snapshot falls back to the host path, so with no
    !! checkpoint enabled the device path is not taken at all.
    type(checkpoint_config_t) :: cfg
    integer(i8) :: bytes

    cfg%snapshot_freq = 1
    cfg%output_stride = [2, 2, 2]
    cfg%checkpoint_freq = 0

    bytes = gpu_io_staging_bytes(dims, cfg, .true.)
    call check("strided snapshot, no checkpoint", bytes, 0_i8)
  end subroutine test_strided_snapshot_no_checkpoint

  subroutine test_checkpoint_strided()
    !! Checkpoints never stride, so a checkpoint alone still takes the
    !! device path regardless of output_stride.
    type(checkpoint_config_t) :: cfg
    integer(i8) :: bytes

    cfg%snapshot_freq = 0
    cfg%checkpoint_freq = 5
    cfg%output_stride = [2, 2, 2]

    bytes = gpu_io_staging_bytes(dims, cfg, .true.)
    call check("checkpoint enabled, strided", bytes, n*int(nbytes, i8))
  end subroutine test_checkpoint_strided

  subroutine test_snapshot_sp_extra()
    !! snapshot_sp sizes the staging buffer to single precision instead of
    !! the build's native precision, but only when the build itself is not
    !! already single precision.
    type(checkpoint_config_t) :: cfg
    integer(i8) :: bytes, expect

    cfg%snapshot_freq = 1
    cfg%output_stride = [1, 1, 1]
    cfg%snapshot_sp = .true.

    bytes = gpu_io_staging_bytes(dims, cfg, .true.)
    if (is_sp) then
      expect = n*int(nbytes, i8)
    else
      expect = n*4_i8
    end if
    call check("single-precision snapshot copy", bytes, expect)
  end subroutine test_snapshot_sp_extra

  subroutine test_no_padding_at_large_dims()
    !! local_dims is an unpadded LOCAL field: no allocator sz-padding
    !! should show up in the result.
    integer, parameter :: big_dims(3) = [256, 256, 128]
    type(checkpoint_config_t) :: cfg
    integer(i8) :: bytes

    cfg%snapshot_freq = 1
    cfg%output_stride = [1, 1, 1]
    cfg%checkpoint_freq = 0
    cfg%snapshot_sp = .false.

    bytes = gpu_io_staging_bytes(big_dims, cfg, .true.)
    call check("no padding at 256x256x128", bytes, &
               8388608_i8*int(nbytes, i8))
  end subroutine test_no_padding_at_large_dims

  subroutine test_neither_enabled_is_zero()
    !! Neither a snapshot nor a checkpoint enabled means the device write
    !! path is never entered at all.
    type(checkpoint_config_t) :: cfg
    integer(i8) :: bytes

    cfg%snapshot_freq = 0
    cfg%checkpoint_freq = 0

    bytes = gpu_io_staging_bytes(dims, cfg, .true.)
    call check("neither snapshot nor checkpoint enabled", bytes, 0_i8)
  end subroutine test_neither_enabled_is_zero

  subroutine check(label, actual, expect)
    character(*), intent(in) :: label
    integer(i8), intent(in) :: actual, expect

    print *, label
    print *, "- Expect: ", expect
    print *, "- Got: ", actual
    if (actual /= expect) then
      test_pass = .false.
    end if
  end subroutine check

end program test_gpu_io_staging
