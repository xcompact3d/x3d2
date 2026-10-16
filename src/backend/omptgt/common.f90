module m_omptgt_common
  !! Constants and device-memory helpers shared by the OpenMP target backend.
  !!
  !! Device arrays are held as contiguous Fortran pointers whose target is a
  !! raw `omp_target_alloc` allocation. The host never dereferences them: they
  !! are only passed as actual arguments to the explicit-shape or assumed-size
  !! dummies of the offloaded kernels, which name those dummies in an
  !! `OMP_DEVICE_ADDR` clause to say "this array already lives on the
  !! device", so the runtime neither maps nor copies it. Referencing such a
  !! pointer directly inside a target region instead makes the runtime map
  !! it implicitly, copying the whole array into a fresh device buffer and
  !! back on every launch.
  !!
  !! `OMP_DEVICE_ADDR` is a preprocessor macro set per compiler in
  !! src/CMakeLists.txt, because no one spelling is accepted everywhere:
  !! OpenMP 5.1 says `has_device_addr`, and its compilers (GNU, Flang, Cray,
  !! ifx) allow `is_device_ptr` only on `type(c_ptr)`; NVHPC rejects
  !! `has_device_addr` but takes the OpenMP 4.5 form, `is_device_ptr` on a
  !! dummy array. The preprocessors expand macros in `!$omp` lines.
  use iso_c_binding, only: c_ptr, c_loc, c_f_pointer, c_associated, &
                           c_size_t, c_sizeof
  use omp_lib, only: omp_target_alloc, omp_target_free, omp_target_memcpy, &
                     omp_get_default_device, omp_get_initial_device

  use m_common, only: dp

  implicit none

  !> Number of points a pencil group stacks together in the leading
  !! dimension, set to the width the target hardware issues a wavefront at so
  !! that one group fills one wavefront. The last branch covers a build that
  !! selects no vendor.
#if defined(OMP_TGT_AMD)
  integer, parameter :: SZ = 64
#elif defined(OMP_TGT_NVIDIA)
  integer, parameter :: SZ = 32
#elif defined(OMP_TGT_INTEL)
  integer, parameter :: SZ = 16
#else
  integer, parameter :: SZ = 16
#endif

  interface dev_alloc
    !! Allocates a device array of the given extents.
    module procedure dev_alloc_1d, dev_alloc_2d, dev_alloc_3d
  end interface dev_alloc

  interface dev_free
    !! Frees a device array allocated with `dev_alloc`.
    module procedure dev_free_1d, dev_free_2d, dev_free_3d
  end interface dev_free

  interface dev_upload
    !! Copies a host array into a device array of the same size.
    module procedure dev_upload_1d, dev_upload_2d, dev_upload_3d
  end interface dev_upload

  interface dev_clone
    !! Allocates a device array of a host array's shape and copies it in.
    module procedure dev_clone_1d, dev_clone_2d
  end interface dev_clone

  interface dev_download
    !! Copies a device array into a host array of the same size.
    module procedure dev_download_1d, dev_download_3d
  end interface dev_download

  private :: dev_alloc_raw, dev_free_raw
  private :: dev_alloc_1d, dev_alloc_2d, dev_alloc_3d
  private :: dev_free_1d, dev_free_2d, dev_free_3d
  private :: dev_upload_1d, dev_upload_2d, dev_upload_3d, dev_upload_
  private :: dev_download_1d, dev_download_3d, dev_download_
  private :: dev_clone_1d, dev_clone_2d

contains

  type(c_ptr) function dev_alloc_raw(n) result(ptr)
    !! Allocates storage for `n` reals on the default device.
    integer, intent(in) :: n

    ptr = omp_target_alloc(max(n, 1)*c_sizeof(0._dp), &
                           omp_get_default_device())
    if (.not. c_associated(ptr)) error stop "omp_target_alloc failed"
  end function dev_alloc_raw

  subroutine dev_free_raw(ptr)
    type(c_ptr), intent(in) :: ptr

    call omp_target_free(ptr, omp_get_default_device())
  end subroutine dev_free_raw

  subroutine dev_alloc_1d(a, n)
    real(dp), pointer, contiguous, intent(out) :: a(:)
    integer, intent(in) :: n

    call c_f_pointer(dev_alloc_raw(n), a, [n])
  end subroutine dev_alloc_1d

  subroutine dev_alloc_2d(a, n1, n2)
    real(dp), pointer, contiguous, intent(out) :: a(:, :)
    integer, intent(in) :: n1, n2

    call c_f_pointer(dev_alloc_raw(n1*n2), a, [n1, n2])
  end subroutine dev_alloc_2d

  subroutine dev_alloc_3d(a, n1, n2, n3)
    real(dp), pointer, contiguous, intent(out) :: a(:, :, :)
    integer, intent(in) :: n1, n2, n3

    call c_f_pointer(dev_alloc_raw(n1*n2*n3), a, [n1, n2, n3])
  end subroutine dev_alloc_3d

  subroutine dev_free_1d(a)
    real(dp), pointer, contiguous, intent(inout) :: a(:)

    if (associated(a)) call dev_free_raw(c_loc(a))
    nullify (a)
  end subroutine dev_free_1d

  subroutine dev_free_2d(a)
    real(dp), pointer, contiguous, intent(inout) :: a(:, :)

    if (associated(a)) call dev_free_raw(c_loc(a))
    nullify (a)
  end subroutine dev_free_2d

  subroutine dev_free_3d(a)
    real(dp), pointer, contiguous, intent(inout) :: a(:, :, :)

    if (associated(a)) call dev_free_raw(c_loc(a))
    nullify (a)
  end subroutine dev_free_3d

  subroutine dev_clone_1d(a_dev, a_host)
    real(dp), pointer, contiguous, intent(out) :: a_dev(:)
    real(dp), intent(in) :: a_host(:)

    call dev_alloc_1d(a_dev, size(a_host))
    call dev_upload_1d(a_dev, a_host)
  end subroutine dev_clone_1d

  subroutine dev_clone_2d(a_dev, a_host)
    real(dp), pointer, contiguous, intent(out) :: a_dev(:, :)
    real(dp), intent(in) :: a_host(:, :)

    call dev_alloc_2d(a_dev, size(a_host, 1), size(a_host, 2))
    call dev_upload_2d(a_dev, a_host)
  end subroutine dev_clone_2d

  subroutine dev_upload_1d(a_dev, a_host)
    real(dp), pointer, contiguous, intent(in) :: a_dev(:)
    real(dp), intent(in) :: a_host(:)

    if (size(a_dev) /= size(a_host)) error stop "dev_upload: size mismatch"
    call dev_upload_(c_loc(a_dev), a_host, size(a_host))
  end subroutine dev_upload_1d

  subroutine dev_upload_2d(a_dev, a_host)
    real(dp), pointer, contiguous, intent(in) :: a_dev(:, :)
    real(dp), intent(in) :: a_host(:, :)

    if (size(a_dev) /= size(a_host)) error stop "dev_upload: size mismatch"
    call dev_upload_(c_loc(a_dev), a_host, size(a_host))
  end subroutine dev_upload_2d

  subroutine dev_upload_3d(a_dev, a_host)
    real(dp), pointer, contiguous, intent(in) :: a_dev(:, :, :)
    real(dp), intent(in) :: a_host(:, :, :)

    if (size(a_dev) /= size(a_host)) error stop "dev_upload: size mismatch"
    call dev_upload_(c_loc(a_dev), a_host, size(a_host))
  end subroutine dev_upload_3d

  subroutine dev_download_1d(a_host, a_dev)
    real(dp), intent(inout) :: a_host(:)
    real(dp), pointer, contiguous, intent(in) :: a_dev(:)

    if (size(a_dev) /= size(a_host)) error stop "dev_download: size mismatch"
    call dev_download_(a_host, c_loc(a_dev), size(a_host))
  end subroutine dev_download_1d

  subroutine dev_download_3d(a_host, a_dev)
    real(dp), intent(inout) :: a_host(:, :, :)
    real(dp), pointer, contiguous, intent(in) :: a_dev(:, :, :)

    if (size(a_dev) /= size(a_host)) error stop "dev_download: size mismatch"
    call dev_download_(a_host, c_loc(a_dev), size(a_host))
  end subroutine dev_download_3d

  subroutine dev_download_(dst, src, n)
    !! As dev_upload_, in the other direction; a non-contiguous actual is
    !! copied back out on return.
    real(dp), target, intent(inout) :: dst(*)
    type(c_ptr), intent(in) :: src
    integer, intent(in) :: n

    call memcpy_d2h(c_loc(dst(1)), src, n)
  end subroutine dev_download_

  subroutine dev_upload_(dst, src, n)
    !! The assumed-size target dummy accepts a host array of any rank and
    !! gives c_loc a valid address even when the actual argument is a
    !! non-contiguous section, which is copied in first.
    type(c_ptr), intent(in) :: dst
    integer, intent(in) :: n
    real(dp), target, intent(in) :: src(*)

    call memcpy_h2d(dst, c_loc(src(1)), n)
  end subroutine dev_upload_

  subroutine memcpy_h2d(dst, src, n)
    !! Copies `n` reals from host address `src` to device address `dst`.
    type(c_ptr), intent(in) :: dst, src
    integer, intent(in) :: n

    if (n < 1) return
    if (omp_target_memcpy(dst, src, n*c_sizeof(0._dp), 0_c_size_t, &
                          0_c_size_t, omp_get_default_device(), &
                          omp_get_initial_device()) /= 0) then
      error stop "omp_target_memcpy host to device failed"
    end if
  end subroutine memcpy_h2d

  subroutine memcpy_d2h(dst, src, n)
    !! Copies `n` reals from device address `src` to host address `dst`.
    type(c_ptr), intent(in) :: dst, src
    integer, intent(in) :: n

    if (n < 1) return
    if (omp_target_memcpy(dst, src, n*c_sizeof(0._dp), 0_c_size_t, &
                          0_c_size_t, omp_get_initial_device(), &
                          omp_get_default_device()) /= 0) then
      error stop "omp_target_memcpy device to host failed"
    end if
  end subroutine memcpy_d2h

  subroutine device_sync()
    !! Waits until the device work queued so far has completed.
    !!
    !! A target region without nowait should complete before the host
    !! carries on, but CCE was seen to return as soon as the kernel is queued
    !! (perf_penta timed 100 solves at 0.1 ms on an MI210, with correct
    !! results). OpenMP has no call to wait for the device, so this queues a
    !! region whose result the host must have: it runs after the work queued
    !! before it, and the host cannot pass it until it has.
    integer :: flag

    flag = 0
    !$omp target map(tofrom: flag)
    flag = 1
    !$omp end target
    if (flag /= 1) error stop "device_sync: target region did not run"

  end subroutine device_sync

  subroutine memcpy_d2d(dst, src, n)
    !! Copies `n` reals between two device addresses.
    type(c_ptr), intent(in) :: dst, src
    integer, intent(in) :: n

    if (n < 1) return
    if (omp_target_memcpy(dst, src, n*c_sizeof(0._dp), 0_c_size_t, &
                          0_c_size_t, omp_get_default_device(), &
                          omp_get_default_device()) /= 0) then
      error stop "omp_target_memcpy device to device failed"
    end if
  end subroutine memcpy_d2d

end module m_omptgt_common
