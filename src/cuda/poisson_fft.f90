module m_cuda_poisson_fft
  use cudafor
  use cufft

  use m_allocator, only: field_t
  use m_common, only: dp
  use m_poisson_fft, only: poisson_fft_t
  use m_tdsops, only: dirps_t

  use m_cuda_allocator, only: cuda_field_t
  use m_cuda_complex, only: reorder_cmplx_x2y_T, reorder_cmplx_y2x_T, &
                            reorder_cmplx_y2z_T, reorder_cmplx_z2y_T, &
                            process_spectral_div_u

  implicit none

  type, extends(poisson_fft_t) :: cuda_poisson_fft_t
    !! FFT based Poisson solver

    !> Local domain sized array to store data in spectral space
    complex(dp), device, allocatable :: c_w_dev(:, :, :)
    !> Local domain sized array storing the spectral equivalence constants
    complex(dp), device, allocatable, dimension(:, :, :) :: waves_dev

    real(dp), device, allocatable, dimension(:) :: ax_dev, bx_dev, &
                                                   ay_dev, by_dev, &
                                                   az_dev, bz_dev

    real(dp), device, allocatable, dimension(:, :, :) :: f_tmp

    integer :: plan3D_fw, plan3D_bw
  contains
    procedure :: fft_forward => fft_forward_cuda
    procedure :: fft_backward => fft_backward_cuda
    procedure :: fft_postprocess => fft_postprocess_cuda
  end type cuda_poisson_fft_t

  interface cuda_poisson_fft_t
    module procedure init
  end interface cuda_poisson_fft_t

  private :: init

contains

  function init(xdirps, ydirps, zdirps) result(poisson_fft)
    implicit none

    class(dirps_t), intent(in) :: xdirps, ydirps, zdirps

    type(cuda_poisson_fft_t) :: poisson_fft

    integer :: nx, ny, nz

    integer :: ierrfft
    integer(int_ptr_kind()) :: worksize

    call poisson_fft%base_init(xdirps, ydirps, zdirps)

    nx = poisson_fft%nx; ny = poisson_fft%ny; nz = poisson_fft%nz

    allocate (poisson_fft%waves_dev(nz/2 + 1, ny, nx))
    poisson_fft%waves_dev = poisson_fft%waves

    allocate (poisson_fft%ax_dev(nx), poisson_fft%bx_dev(nx))
    allocate (poisson_fft%ay_dev(ny), poisson_fft%by_dev(ny))
    allocate (poisson_fft%az_dev(nz), poisson_fft%bz_dev(nz))
    poisson_fft%ax_dev = poisson_fft%ax; poisson_fft%bx_dev = poisson_fft%bx
    poisson_fft%ay_dev = poisson_fft%ay; poisson_fft%by_dev = poisson_fft%by
    poisson_fft%az_dev = poisson_fft%az; poisson_fft%bz_dev = poisson_fft%bz

    allocate (poisson_fft%c_w_dev(nz/2 + 1, ny, nx))

    ! We can't currently ask allocator to pass us an array with
    ! exact shape we want, so we allocate an extra array here.
    ! This will be removed when allocator is fixed.
    allocate (poisson_fft%f_tmp(nx, ny, nz))

    ! 3D plans
    ierrfft = cufftPlan3D(poisson_fft%plan3D_fw, nz, ny, nx, CUFFT_D2Z)

    ierrfft = cufftPlan3D(poisson_fft%plan3D_bw, nz, ny, nx, CUFFT_Z2D)

  end function init

  subroutine fft_forward_cuda(self, f)
    implicit none

    class(cuda_poisson_fft_t) :: self
    class(field_t), intent(in) :: f

    real(dp), device, pointer, dimension(:, :, :) :: f_dev

    type(dim3) :: blocks, threads
    integer :: ierrfft

    select type (f); type is (cuda_field_t); f_dev => f%data_d; end select

    self%f_tmp = f_dev
    ierrfft = cufftExecD2Z(self%plan3D_fw, self%f_tmp, self%c_w_dev)

  end subroutine fft_forward_cuda

  subroutine fft_backward_cuda(self, f)
    implicit none

    class(cuda_poisson_fft_t) :: self
    class(field_t), intent(inout) :: f

    real(dp), device, pointer, dimension(:, :, :) :: f_dev

    type(dim3) :: blocks, threads
    integer :: ierrfft

    select type (f); type is (cuda_field_t); f_dev => f%data_d; end select

    ierrfft = cufftExecZ2D(self%plan3D_bw, self%c_w_dev, self%f_tmp)
    f_dev = self%f_tmp

  end subroutine fft_backward_cuda

  subroutine fft_postprocess_cuda(self)
    implicit none

    class(cuda_poisson_fft_t) :: self

    complex(dp), device, dimension(:, :, :), pointer :: c_dev
    type(dim3) :: blocks, threads

    ! Postprocess
    blocks = dim3((self%nz/2 + 1), 1, 1)
    threads = dim3(self%nx, 1, 1)
    call process_spectral_div_u<<<blocks, threads>>>( & !&
      self%c_w_dev, self%waves_dev, self%nx, self%ny, self%nz, &
      self%ax_dev, self%bx_dev, self%ay_dev, self%by_dev, &
      self%az_dev, self%bz_dev &
      )

  end subroutine fft_postprocess_cuda

end module m_cuda_poisson_fft
