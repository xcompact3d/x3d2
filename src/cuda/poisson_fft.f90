module m_cuda_poisson_fft
   use cudafor
   use cufft

   use m_allocator, only: field_t
   use m_common, only: dp
   use m_poisson_fft, only: poisson_fft_t
   use m_tdsops, only: dirps_t

   use m_cuda_allocator, only: cuda_field_t
   use m_cuda_common, only: SZ
   use m_cuda_complex, only: processfftdiv

   implicit none

   type, extends(poisson_fft_t) :: cuda_poisson_fft_t
      !! FFT based Poisson solver
      !! It can only handle 1D decompositions along z direction.
      complex(dp), device, pointer, dimension(:) :: c_x_dev, c_y_dev, c_z_dev
      complex(dp), device, allocatable, dimension(:, :, :) :: waves_dev

      real(dp), device, allocatable, dimension(:) :: ax_dev, bx_dev, &
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

      integer :: nx, ny, nz

      integer :: ierrfft
      integer(int_ptr_kind()) :: worksize

      call poisson_fft%base_init(xdirps, ydirps, zdirps, SZ)

      nx = poisson_fft%nx; ny = poisson_fft%ny; nz = poisson_fft%nz

      allocate (poisson_fft%waves_dev(SZ, nx, (ny*(nz/2 + 1))/SZ))
      poisson_fft%waves_dev = poisson_fft%waves

      allocate (poisson_fft%ax_dev(nx), poisson_fft%bx_dev(nx))
      allocate (poisson_fft%ay_dev(ny), poisson_fft%by_dev(ny))
      allocate (poisson_fft%az_dev(nz), poisson_fft%bz_dev(nz))
      poisson_fft%ax_dev = poisson_fft%ax; poisson_fft%bx_dev = poisson_fft%bx
      poisson_fft%ay_dev = poisson_fft%ay; poisson_fft%by_dev = poisson_fft%by
      poisson_fft%az_dev = poisson_fft%az; poisson_fft%bz_dev = poisson_fft%bz

      allocate (poisson_fft%c_x_dev(nx*ny*(nz/2 + 1)))
      allocate (poisson_fft%c_y_dev(nx*ny*(nz/2 + 1)))
      allocate (poisson_fft%c_z_dev(nx*ny*(nz/2 + 1)))

      ! set cufft plans
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

   subroutine fft_forward_cuda(self, f)
      implicit none

      class(cuda_poisson_fft_t) :: self
      class(field_t), intent(in) :: f

      real(dp), device, pointer, dimension(:, :, :) :: f_dev
      integer :: ierrfft

      select type(f); type is (cuda_field_t); f_dev => f%data_d; end select

      ! First reorder f into cartesian-like data structure

      ! Forward FFT transform in z from real to complex
      ierrfft = cufftExecD2Z(self%planD2Zz, f_dev, self%c_z_dev)

      ! Reorder from z to y

      ! In-place forward FFT in y
      ierrfft = cufftExecZ2Z(self%planZ2Zy, self%c_y_dev, self%c_y_dev, &
                             CUFFT_FORWARD)

      ! Reorder from y to x

      ! In-place forward FFT in x
      ierrfft = cufftExecZ2Z(self%planZ2Zx, self%c_x_dev, self%c_x_dev, &
                             CUFFT_FORWARD)

   end subroutine fft_forward_cuda

   subroutine fft_backward_cuda(self, f)
      implicit none

      class(cuda_poisson_fft_t) :: self
      class(field_t), intent(inout) :: f

      real(dp), device, pointer, dimension(:, :, :) :: f_dev
      integer :: ierrfft

      select type(f); type is (cuda_field_t); f_dev => f%data_d; end select

      ! In-place backward FFT in x
      ierrfft = cufftExecZ2Z(self%planZ2Zx, self%c_x_dev, self%c_x_dev, &
                             CUFFT_INVERSE)

      ! Reorder from x to y

      ! In-place backward FFT in y
      ierrfft = cufftExecZ2Z(self%planZ2Zy, self%c_y_dev, self%c_y_dev, &
                             CUFFT_INVERSE)

      ! Reorder from y to z

      ! Backward FFT transform in z from complex to real
      ierrfft = cufftExecZ2D(self%planZ2Dz, self%c_z_dev, f_dev)

      ! Finally reorder f back into our specialist data structure

   end subroutine fft_backward_cuda

   subroutine fft_postprocess_cuda(self)
      implicit none

      class(cuda_poisson_fft_t) :: self

      complex(dp), device, dimension(:, :, :), pointer :: c_dev
      type(dim3) :: blocks, threads

      ! Reshape from cartesian-like to our specialist data structure

      blocks = dim3((self%ny*(self%nz/2 + 1))/SZ, 1 , 1)
      threads = dim3(SZ, 1, 1)

      c_dev(1:SZ, 1:self%nx, 1:(self%ny*(self%nz/2 + 1))/SZ) => self%c_y_dev

      call processfftdiv<<<blocks, threads>>>( &
         c_dev, self%waves_dev, self%nx, self%ny, self%nz, &
         self%ax_dev, self%bx_dev, self%ay_dev, self%by_dev, &
         self%az_dev, self%bz_dev &
         )

      ! Reshape from our specialist data structure to cartesian-like

   end subroutine fft_postprocess_cuda

end module m_cuda_poisson_fft
