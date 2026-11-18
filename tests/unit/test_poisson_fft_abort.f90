module m_poisson_stub
  !! A case type that overrides nothing beyond the deferred transforms, so
  !! it reaches the default periodicity reorder of poisson_fft_t.
  use m_field, only: field_t
  use m_poisson_fft, only: poisson_fft_t

  implicit none

  type, extends(poisson_fft_t) :: poisson_stub_t
  contains
    procedure :: fft_forward => stub_forward
    procedure :: fft_backward => stub_backward
    procedure :: fft_postprocess => stub_postprocess
  end type poisson_stub_t

contains

  subroutine stub_forward(self, f_in)
    class(poisson_stub_t) :: self
    class(field_t), intent(in) :: f_in
  end subroutine stub_forward

  subroutine stub_backward(self, f_out)
    class(poisson_stub_t) :: self
    class(field_t), intent(inout) :: f_out
  end subroutine stub_backward

  subroutine stub_postprocess(self)
    class(poisson_stub_t) :: self
  end subroutine stub_postprocess

end module m_poisson_stub

program test_poisson_fft_abort
  !! Reaches one of the FFT Poisson solver's aborts and must not return.
  !!
  !! ctest passes this test when the expected error message appears in the
  !! output, as test_logger_abort does. The first argument picks the abort:
  !!
  !!   - `unsupported_bc`: a mesh non-periodic in z (BC code 001), which the
  !!     backend's factory rejects before it initialises anything;
  !!   - `enforce`, `undo`: a case type that leaves enforce_periodicity or
  !!     undo_periodicity to the default of poisson_fft_t.
  use m_backend_runtime, only: backend_runtime_t, backend_is_cuda
  use m_common, only: dp
  use m_field, only: field_t
  use m_mesh, only: mesh_t
  use m_poisson_stub, only: poisson_stub_t
  use m_tdsops, only: dirps_t
  use m_test_utils, only: initialise_mpi

  implicit none

  integer :: nrank, nproc
  character(len=32) :: mode
  type(poisson_stub_t) :: stub
  type(field_t) :: f, temp
  type(mesh_t), target :: mesh
  type(backend_runtime_t), target :: runtime
  type(dirps_t) :: xdirps, ydirps, zdirps
  character(len=20) :: BC_periodic(2), BC_dirichlet(2)

  call initialise_mpi(nrank, nproc)
  mode = ''
  if (command_argument_count() >= 1) call get_command_argument(1, mode)

  select case (trim(mode))
  case ('unsupported_bc')
    BC_periodic = ['periodic', 'periodic']
    BC_dirichlet = ['dirichlet', 'dirichlet']
    mesh = mesh_t([32, 32, 33], [1, 1, nproc], [1.0_dp, 1.0_dp, 1.0_dp], &
                  BC_periodic, BC_periodic, BC_dirichlet, &
                  use_2decomp=.not. backend_is_cuda)
    call runtime%init(mesh)
    ! The factory stops on the BC code before it reads the operators, so
    ! the directions need no setting up.
    call runtime%backend%init_poisson_fft(mesh, xdirps, ydirps, zdirps)
  case ('enforce')
    call stub%enforce_periodicity(f, temp)
  case ('undo')
    call stub%undo_periodicity(f, temp)
  case default
    print *, 'unknown mode: '//trim(mode)
    error stop 1
  end select

  print *, 'the Poisson solver did not abort: this line must never be reached'

end program test_poisson_fft_abort
