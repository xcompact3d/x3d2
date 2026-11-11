program test_memcheck_estimate_010
  !! Check the 010 terms of spectral_plus_mirror_bytes, the caller of the
  !! x3d2-memcheck estimate, at 2 GPUs. The stretched case adds the matrices
  !! and the three z-slab buffers and no mirror buffer, since the solver skips
  !! init_dim2_mirror there. The uniform case adds the mirror buffers only.

  use m_common, only: i8, nbytes
  use m_mesh, only: FFT_BC_010
  use m_memory_estimate, only: spectral_slab_bytes
  use m_memcheck_context, only: memcheck_ctx_t
  use m_memcheck_estimate, only: spectral_plus_mirror_bytes
  use m_test_utils, only: check_i8

  implicit none

  logical :: test_pass
  type(memcheck_ctx_t) :: ctx
  integer :: n_spec(3)
  integer(i8) :: slab, expect

  test_pass = .true.

  ctx%bc_code = FFT_BC_010
  ctx%bc_is_000 = .false.
  ctx%bc_is_010 = .true.
  ctx%bc_is_100 = .false.
  ctx%bc_is_110 = .false.
  ctx%cdims = [64, 32, 64]
  ctx%solver_cfg%lowmem_fft = .false.

  ! The spectral pencil of one GPU at 2 GPUs: (nx/2 + 1, ny/2, nz)
  n_spec = [33, 16, 64]
  slab = spectral_slab_bytes(n_spec)

  ! The matrices cover the z-slab (33, 32, 32) with ny/2 = 16 rows in each of
  ! the four arrays, doubled by the copies. The three complex z-slab buffers
  ! follow.
  ctx%domain_cfg%stretching(2) = 'top-bottom'
  expect = slab + (2_i8*4_i8*33_i8*16_i8*32_i8*5_i8 &
                   + 3_i8*33_i8*32_i8*32_i8*2_i8)*int(nbytes, i8)
  call check_i8("stretched 010, 2 GPUs", spectral_plus_mirror_bytes(ctx, 2), &
                expect, test_pass)

  ! Uniform y: the c_mirror_dev, c_slab_send_dev and two plane buffers
  ctx%domain_cfg%stretching(2) = 'uniform'
  expect = slab + (2_i8*33_i8*16_i8*64_i8 + 2_i8*33_i8*64_i8)*2_i8 &
           *int(nbytes, i8)
  call check_i8("uniform 010, 2 GPUs", spectral_plus_mirror_bytes(ctx, 2), &
                expect, test_pass)

  if (.not. test_pass) then
    error stop "FAIL"
  end if

contains

end program test_memcheck_estimate_010
