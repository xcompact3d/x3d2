module m_memcheck_device
  !! The backend seam of x3d2-memcheck: everything the generic estimate,
  !! report and build modules need from a compute backend, so that only one
  !! implementation (src/backend/cuda/memcheck_device.f90 for CUDA) touches
  !! backend-specific libraries.
  use m_common, only: dp, i8
  use m_mesh, only: mesh_t
  use m_allocator, only: allocator_t
  use m_base_backend, only: base_backend_t
  use m_base_case, only: base_case_t
  implicit none

  private
  public :: memcheck_device_t, busy_card_gib

  !> Card usage at start (GiB) above which another process is assumed to
  !> share the card, so the measured context size is not trusted.
  real(dp), parameter :: busy_card_gib = 0.5_dp

  type, abstract :: memcheck_device_t
  contains
    !> Select the device. ok=.false. and the line to print in msg if there
    !> is no usable one.
    procedure(init_iface), deferred :: init
    !> Total and currently free device memory.
    procedure(mem_info_iface), deferred :: mem_info
    !> Card memory used at start, in bytes: this process's device context,
    !> plus any other process already on the card.
    procedure(context_bytes_iface), deferred :: context_bytes
    !> The device-context term of the estimate, in bytes: context_bytes
    !> when the card was otherwise idle at start (at most busy_card_gib
    !> used), else the built-in A100 constant.
    procedure(context_term_bytes_iface), deferred :: context_term_bytes
    !> Pencil size (the SZ of the backend's field layout).
    procedure(sz_iface), deferred :: sz
    !> Tier 1 floor of the per-GPU overhead at ng GPUs, in bytes: the
    !> device context plus (unless bc_is_110) the distributed FFT heap.
    procedure(overhead_floor_bytes_iface), deferred :: overhead_floor_bytes
    !> Tier 2 per-GPU overhead at ng GPUs, in bytes, from a one-shot
    !> throwaway FFT plan query (run on the first call, cached after):
    !> FFT workspace/heap/data buffer plus the device context.
    procedure(overhead_query_bytes_iface), deferred :: overhead_query_bytes
    !> Whether the Tier 2 query has run and, if so, residual_gib: what it
    !> left allocated on the device (~0 if the plan was fully released).
    procedure(query_residual_gib_iface), deferred :: query_residual_gib
    !> Build the device allocator and the backend on mesh; the device owns
    !> both, so the returned pointers stay valid for as long as it does.
    procedure(make_backend_iface), deferred :: make_backend
    !> Record what the real solver of flow_case settled on (called once it
    !> is constructed), for fft_path_known and uses_distributed_fft.
    procedure(capture_solver_info_iface), deferred :: capture_solver_info
    !> Whether capture_solver_info has seen the real solver's FFT path.
    procedure(fft_path_known_iface), deferred :: fft_path_known
    !> Whether the real solver's FFT is distributed across ranks (cuFFTMp
    !> on CUDA) rather than a plain per-rank transform; only meaningful
    !> once fft_path_known.
    procedure(uses_distributed_fft_iface), deferred :: uses_distributed_fft
  end type memcheck_device_t

  abstract interface
    subroutine init_iface(self, ok, msg)
      import :: memcheck_device_t
      class(memcheck_device_t), intent(inout) :: self
      logical, intent(out) :: ok
      character(len=*), intent(out) :: msg
    end subroutine init_iface

    subroutine mem_info_iface(self, total_bytes, free_bytes)
      import :: memcheck_device_t, i8
      class(memcheck_device_t), intent(in) :: self
      integer(i8), intent(out) :: total_bytes, free_bytes
    end subroutine mem_info_iface

    function context_bytes_iface(self) result(nbytes8)
      import :: memcheck_device_t, i8
      class(memcheck_device_t), intent(in) :: self
      integer(i8) :: nbytes8
    end function context_bytes_iface

    function context_term_bytes_iface(self) result(nbytes8)
      import :: memcheck_device_t, i8
      class(memcheck_device_t), intent(in) :: self
      integer(i8) :: nbytes8
    end function context_term_bytes_iface

    function sz_iface(self) result(pencil_size)
      import :: memcheck_device_t
      class(memcheck_device_t), intent(in) :: self
      integer :: pencil_size
    end function sz_iface

    function overhead_floor_bytes_iface(self, ng, bc_is_110) result(nbytes8)
      import :: memcheck_device_t, i8
      class(memcheck_device_t), intent(in) :: self
      integer, intent(in) :: ng
      logical, intent(in) :: bc_is_110
      integer(i8) :: nbytes8
    end function overhead_floor_bytes_iface

    function overhead_query_bytes_iface(self, ng, bc_is_100, bc_is_110, &
                                        cdims, is_root) result(nbytes8)
      import :: memcheck_device_t, i8
      class(memcheck_device_t), intent(inout) :: self
      integer, intent(in) :: ng, cdims(3)
      logical, intent(in) :: bc_is_100, bc_is_110, is_root
      integer(i8) :: nbytes8
    end function overhead_query_bytes_iface

    function query_residual_gib_iface(self, residual_gib) result(ran)
      import :: memcheck_device_t, dp
      class(memcheck_device_t), intent(in) :: self
      real(dp), intent(out) :: residual_gib
      logical :: ran
    end function query_residual_gib_iface

    subroutine make_backend_iface(self, mesh, dims, allocator, backend)
      import :: memcheck_device_t, mesh_t, allocator_t, base_backend_t
      class(memcheck_device_t), target, intent(inout) :: self
      type(mesh_t), target, intent(inout) :: mesh
      integer, intent(in) :: dims(3)
      class(allocator_t), pointer, intent(out) :: allocator
      class(base_backend_t), pointer, intent(out) :: backend
    end subroutine make_backend_iface

    subroutine capture_solver_info_iface(self, flow_case)
      import :: memcheck_device_t, base_case_t
      class(memcheck_device_t), intent(inout) :: self
      class(base_case_t), intent(in) :: flow_case
    end subroutine capture_solver_info_iface

    function fft_path_known_iface(self) result(known)
      import :: memcheck_device_t
      class(memcheck_device_t), intent(in) :: self
      logical :: known
    end function fft_path_known_iface

    function uses_distributed_fft_iface(self) result(distributed)
      import :: memcheck_device_t
      class(memcheck_device_t), intent(in) :: self
      logical :: distributed
    end function uses_distributed_fft_iface
  end interface

end module m_memcheck_device
