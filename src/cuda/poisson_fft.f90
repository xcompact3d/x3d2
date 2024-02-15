module m_cuda_poisson_fft
   use cudafor

   use m_allocator, only: field_t
   use m_common, only: dp
   use m_poisson_fft, only: poisson_fft_t
   use m_tdsops, only: dirps_t

   implicit none

   type, extends(poisson_fft_t) :: cuda_poisson_fft_t
      !! FFT based Poisson solver
      !! It can only handle 1D decompositions along z direction.
      complex(dp), device, allocatable, dimension(:, :, :) :: &
         c_x_dev, c_y_dev, c_z_dev, waves_dev
      complex(dp), device, allocatable, dimension(:) :: ax_dev, bx_dev, &
         ay_dev, by_dev, az_dev, bz_dev
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
      integer :: nx, ny, nz

      call poisson_fft%base_init(xdirps, ydirps, zdirps)

      nx = poisson_fft%nx; ny = poisson_fft%ny; nz = poisson_fft%nz

      allocate (poisson_fft%waves_dev(nx, ny, nz))
      poisson_fft%waves_dev = poisson_fft%waves

      allocate (poisson_fft%ax_dev(nx), poisson_fft%bx_dev(nx))
      allocate (poisson_fft%ay_dev(ny), poisson_fft%by_dev(ny))
      allocate (poisson_fft%az_dev(nz), poisson_fft%bz_dev(nz))
      poisson_fft%ax_dev = poisson_fft%ax; poisson_fft%bx_dev = poisson_fft%bx
      poisson_fft%ay_dev = poisson_fft%ay; poisson_fft%by_dev = poisson_fft%by
      poisson_fft%az_dev = poisson_fft%az; poisson_fft%bz_dev = poisson_fft%bz

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
