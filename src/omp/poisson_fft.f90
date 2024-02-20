module m_omp_poisson_fft
   use m_allocator, only: field_t
   use m_common, only: dp
   use m_poisson_fft, only: poisson_fft_t
   use m_tdsops, only: dirps_t

   use m_omp_common, only: SZ

   implicit none

   type, extends(poisson_fft_t) :: omp_poisson_fft_t
      !! FFT based Poisson solver
      !! It can only handle 1D decompositions along z direction.
      complex(dp), allocatable, dimension(:, :, :) :: c_x, c_y, c_z
   contains
      procedure :: fft_forward => fft_forward_omp
      procedure :: fft_backward => fft_backward_omp
      procedure :: fft_postprocess => fft_postprocess_omp
   end type omp_poisson_fft_t

   interface omp_poisson_fft_t
      module procedure init
   end interface omp_poisson_fft_t

contains

   function init(xdirps, ydirps, zdirps) result(poisson_fft)
      implicit none

      class(dirps_t), intent(in) :: xdirps, ydirps, zdirps

      type(omp_poisson_fft_t) :: poisson_fft
      integer :: nx, ny, nz

      call poisson_fft%base_init(xdirps, ydirps, zdirps, SZ)

   end function init

   subroutine fft_forward_omp(self, f_in)
      implicit none

      class(omp_poisson_fft_t) :: self
      class(field_t), intent(in) :: f_in
   end subroutine fft_forward_omp

   subroutine fft_backward_omp(self, f_out)
      implicit none

      class(omp_poisson_fft_t) :: self
      class(field_t), intent(inout) :: f_out
   end subroutine fft_backward_omp

   subroutine fft_postprocess_omp(self)
      implicit none

      class(omp_poisson_fft_t) :: self
   end subroutine fft_postprocess_omp

end module m_omp_poisson_fft
