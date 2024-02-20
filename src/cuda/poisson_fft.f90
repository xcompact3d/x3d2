module m_cuda_poisson_fft
   use cudafor
   use cufft

   use m_allocator, only: field_t
   use m_common, only: dp
   use m_poisson_fft, only: poisson_fft_t
   use m_tdsops, only: dirps_t

   use m_cuda_common, only: SZ

   implicit none

   type, extends(poisson_fft_t) :: cuda_poisson_fft_t
      !! FFT based Poisson solver
      !! It can only handle 1D decompositions along z direction.
      complex(dp), device, allocatable, dimension(:, :, :) :: &
         c_x_dev, c_y_dev, c_z_dev, waves_dev
      complex(dp), device, allocatable, dimension(:) :: ax_dev, bx_dev, &
         ay_dev, by_dev, az_dev, bz_dev
      integer :: planD2Zz, planZ2Dz, planZ2Zx, planZ2Zy
   contains
      procedure :: fft_forward => fft_forward_cuda
      procedure :: fft_backward => fft_backward_cuda
      procedure :: fft_postprocess => fft_postprocess_cuda
   end type cuda_poisson_fft_t

   interface cuda_poisson_fft_t
      module procedure init
   end interface cuda_poisson_fft_t

contains

   function init(xdirps, ydirps, zdirps) result(poisson_fft)
      implicit none

      class(dirps_t), intent(in) :: xdirps, ydirps, zdirps

      type(cuda_poisson_fft_t) :: poisson_fft

      integer :: nx, ny, nz, nx_loc, ny_loc, nz_loc

      integer :: ierrfft
      integer(int_ptr_kind()) :: worksize

      call poisson_fft%base_init(xdirps, ydirps, zdirps, SZ)

      nx = poisson_fft%nx; ny = poisson_fft%ny; nz = poisson_fft%nz
      nx_loc = nx; ny_loc = ny; nz_loc = nz

      allocate (poisson_fft%waves_dev(nx, ny, nz))
      poisson_fft%waves_dev = poisson_fft%waves

      allocate (poisson_fft%ax_dev(nx), poisson_fft%bx_dev(nx))
      allocate (poisson_fft%ay_dev(ny), poisson_fft%by_dev(ny))
      allocate (poisson_fft%az_dev(nz), poisson_fft%bz_dev(nz))
      poisson_fft%ax_dev = poisson_fft%ax; poisson_fft%bx_dev = poisson_fft%bx
      poisson_fft%ay_dev = poisson_fft%ay; poisson_fft%by_dev = poisson_fft%by
      poisson_fft%az_dev = poisson_fft%az; poisson_fft%bz_dev = poisson_fft%bz

      allocate (poisson_fft%c_x_dev(nx, SZ, (ny*(nz/2 + 1))/SZ))
      allocate (poisson_fft%c_y_dev(ny, SZ, (nx*(nz/2 + 1))/SZ))
      allocate (poisson_fft%c_z_dev(nz/2 + 1, SZ, nx*ny/SZ))

      ! plans for regular for loop executions in a single stream
      ierrfft = cufftCreate(poisson_fft%planD2Zz)
      ierrfft = cufftMakePlanMany(poisson_fft%planD2Zz, 1, nz, &
                                  nz, 1, nz, nz/2+1, 1, nz/2+1, &
                                  CUFFT_D2Z, nx*ny, worksize)
      ierrfft = cufftSetWorkArea(poisson_fft%planD2Zz, poisson_fft%c_x_dev)

      ierrfft = cufftCreate(poisson_fft%planZ2Dz)
      ierrfft = cufftMakePlanMany(poisson_fft%planZ2Dz, 1, nz, &
                                  nz/2+1, 1, nz/2+1, nz, 1, nz, &
                                  CUFFT_Z2D, nx*ny, worksize)
      ierrfft = cufftSetWorkArea(poisson_fft%planZ2Dz, poisson_fft%c_x_dev)

      ierrfft = cufftCreate(poisson_fft%planZ2Zy)
      ierrfft = cufftMakePlanMany(poisson_fft%planZ2Zy, 1, ny, &
                                  ny, 1, ny, ny, 1, ny, &
                                  CUFFT_Z2Z, nx*(nz/2 + 1), worksize)
      ierrfft = cufftSetWorkArea(poisson_fft%planZ2Zy, poisson_fft%c_x_dev)

      ierrfft = cufftCreate(poisson_fft%planZ2Zx)
      ierrfft = cufftMakePlanMany(poisson_fft%planZ2Zx, 1, nx, &
                                  nx, 1, nx, nx, 1, nx, &
                                  CUFFT_Z2Z, ny*(nz/2 + 1), worksize)
      ierrfft = cufftSetWorkArea(poisson_fft%planZ2Zx, poisson_fft%c_y_dev)

   end function init

   subroutine fft_forward_cuda(self, f_in)
      implicit none

      class(cuda_poisson_fft_t) :: self
      class(field_t), intent(in) :: f_in
   end subroutine fft_forward_cuda

   subroutine fft_backward_cuda(self, f_out)
      implicit none

      class(cuda_poisson_fft_t) :: self
      class(field_t), intent(inout) :: f_out
   end subroutine fft_backward_cuda

   subroutine fft_postprocess_cuda(self)
      implicit none

      class(cuda_poisson_fft_t) :: self
   end subroutine fft_postprocess_cuda

end module m_cuda_poisson_fft
