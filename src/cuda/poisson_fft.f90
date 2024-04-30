module m_cuda_poisson_fft
  use iso_c_binding, only: c_loc, c_ptr, c_f_pointer
  use iso_fortran_env, only: stderr => error_unit
  use cudafor
  use cufft

  use m_allocator, only: field_t
  use m_common, only: dp
  use m_poisson_fft, only: poisson_fft_t
  use m_tdsops, only: dirps_t

  use m_cuda_allocator, only: cuda_field_t
  use m_cuda_spectral, only: process_spectral_div_u

  implicit none

  type, extends(poisson_fft_t) :: cuda_poisson_fft_t
    !! FFT based Poisson solver

    !> Local domain sized array to store data in spectral space
    complex(dp), device, allocatable, dimension(:, :, :) :: c_w_dev
    !> Local domain sized array storing the spectral equivalence constants
    complex(dp), device, allocatable, dimension(:, :, :) :: waves_dev
    !> cufft requires a local domain sized storage
    complex(dp), device, allocatable, dimension(:) :: fft_worksize

    real(dp), device, allocatable, dimension(:) :: ax_dev, bx_dev, &
                                                   ay_dev, by_dev, &
                                                   az_dev, bz_dev

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

    integer :: ierr
    integer(int_ptr_kind()) :: worksize

    call poisson_fft%base_init(xdirps, ydirps, zdirps)

    nx = poisson_fft%nx; ny = poisson_fft%ny; nz = poisson_fft%nz

    allocate (poisson_fft%waves_dev(nx/2 + 1, ny, nz))
    poisson_fft%waves_dev = poisson_fft%waves

    allocate (poisson_fft%ax_dev(nx), poisson_fft%bx_dev(nx))
    allocate (poisson_fft%ay_dev(ny), poisson_fft%by_dev(ny))
    allocate (poisson_fft%az_dev(nz), poisson_fft%bz_dev(nz))
    poisson_fft%ax_dev = poisson_fft%ax; poisson_fft%bx_dev = poisson_fft%bx
    poisson_fft%ay_dev = poisson_fft%ay; poisson_fft%by_dev = poisson_fft%by
    poisson_fft%az_dev = poisson_fft%az; poisson_fft%bz_dev = poisson_fft%bz

    allocate (poisson_fft%c_w_dev(nx/2 + 1, ny, nz))
    allocate (poisson_fft%fft_worksize((nx/2 + 1)*ny*nz))

    ! 3D plans
    ierr = cufftCreate(poisson_fft%plan3D_fw)
    ierr = cufftMakePlan3D(poisson_fft%plan3D_fw, nz, ny, nx, CUFFT_D2Z, &
                           worksize)
    ierr = cufftSetWorkArea(poisson_fft%plan3D_fw, poisson_fft%fft_worksize)
    if (ierr /= 0) then
      write (stderr, *), "cuFFT Error Code: ", ierr
      error stop 'Forward 3D FFT plan generation failed'
    end if

    ierr = cufftCreate(poisson_fft%plan3D_bw)
    ierr = cufftMakePlan3D(poisson_fft%plan3D_bw, nz, ny, nx, CUFFT_Z2D, &
                           worksize)
    ierr = cufftSetWorkArea(poisson_fft%plan3D_bw, poisson_fft%fft_worksize)
    if (ierr /= 0) then
      write (stderr, *), "cuFFT Error Code: ", ierr
      error stop 'Backward 3D FFT plan generation failed'
    end if

  end function init

  subroutine fft_forward_cuda(self, f)
    implicit none

    class(cuda_poisson_fft_t) :: self
    class(field_t), intent(in) :: f

    real(dp), device, pointer, dimension(:, :, :) :: f_dev
    real(dp), device, pointer :: f_ptr
    type(c_ptr) :: f_c_ptr

    type(dim3) :: blocks, threads
    integer :: ierr

    select type (f); type is (cuda_field_t); f_dev => f%data_d; end select

    ! Using f_dev directly in cufft call causes a segfault
    ! Pointer switches below fixes the problem
    f_c_ptr = c_loc(f_dev)
    call c_f_pointer(f_c_ptr, f_ptr)

    ierr = cufftExecD2Z(self%plan3D_fw, f_ptr, self%c_w_dev)
    if (ierr /= 0) then
      write (stderr, *), "cuFFT Error Code: ", ierr
      error stop 'Forward 3D FFT execution failed'
    end if

  end subroutine fft_forward_cuda

  subroutine fft_backward_cuda(self, f)
    implicit none

    class(cuda_poisson_fft_t) :: self
    class(field_t), intent(inout) :: f

    real(dp), device, pointer, dimension(:, :, :) :: f_dev
    real(dp), device, pointer :: f_ptr
    type(c_ptr) :: f_c_ptr

    type(dim3) :: blocks, threads
    integer :: ierr

    select type (f); type is (cuda_field_t); f_dev => f%data_d; end select

    ! Using f_dev directly in cufft call causes a segfault
    ! Pointer switches below fixes the problem
    f_c_ptr = c_loc(f_dev)
    call c_f_pointer(f_c_ptr, f_ptr)

    ierr = cufftExecZ2D(self%plan3D_bw, self%c_w_dev, f_ptr)
    if (ierr /= 0) then
      write (stderr, *), "cuFFT Error Code: ", ierr
      error stop 'Backward 3D FFT execution failed'
    end if

  end subroutine fft_backward_cuda

  subroutine fft_postprocess_cuda(self)
    implicit none

    class(cuda_poisson_fft_t) :: self

    complex(dp), device, dimension(:, :, :), pointer :: c_dev
    type(dim3) :: blocks, threads
    integer :: tsize

    ! tsize is different than SZ, because here we work on a 3D Cartesian
    ! data structure, and free to specify any suitable thread/block size.
    tsize = 16
    blocks = dim3((self%ny - 1)/tsize + 1, self%nz, 1)
    threads = dim3(tsize, 1, 1)

    ! Postprocess div_u in spectral space
    call process_spectral_div_u<<<blocks, threads>>>( & !&
      self%c_w_dev, self%waves_dev, self%nx, self%ny, self%nz, &
      self%ax_dev, self%bx_dev, self%ay_dev, self%by_dev, &
      self%az_dev, self%bz_dev &
      )

  end subroutine fft_postprocess_cuda

end module m_cuda_poisson_fft
