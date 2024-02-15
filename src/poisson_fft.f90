module m_poisson_fft
   use m_allocator, only: field_t
   use m_common, only: dp, pi
   use m_tdsops, only: dirps_t

   implicit none

   type, abstract :: poisson_fft_t
      !! FFT based Poisson solver
      !! It can only handle 1D decompositions along z direction.
      integer :: nx, ny, nz
      complex(dp), allocatable, dimension(:, :, :) :: waves
      complex(dp), allocatable, dimension(:) :: ax, bx, ay, by, az, bz
   contains
      procedure(fft_forward), deferred :: fft_forward
      procedure(fft_backward), deferred :: fft_backward
      procedure(fft_postprocess), deferred :: fft_postprocess
      procedure :: base_init
      procedure :: waves_set
   end type poisson_fft_t

   abstract interface
      subroutine fft_forward(self, f_in)
         import :: poisson_fft_t
         import :: field_t
         implicit none

         class(poisson_fft_t) :: self
         class(field_t), intent(in) :: f_in
      end subroutine fft_forward

      subroutine fft_backward(self, f_out)
         import :: poisson_fft_t
         import :: field_t
         implicit none

         class(poisson_fft_t) :: self
         class(field_t), intent(inout) :: f_out
      end subroutine fft_backward

      subroutine fft_postprocess(self)
         import :: poisson_fft_t
         implicit none

         class(poisson_fft_t) :: self
      end subroutine fft_postprocess
   end interface

contains

   subroutine base_init(self, xdirps, ydirps, zdirps)
      implicit none

      class(poisson_fft_t) :: self
      class(dirps_t), intent(in) :: xdirps, ydirps, zdirps

      integer :: nx, ny, nz

      self%nx = xdirps%n; self%ny = ydirps%n; self%nz = zdirps%n

      allocate (self%ax(self%nx), self%bx(self%nx))
      allocate (self%ay(self%nx), self%by(self%nx))
      allocate (self%az(self%nx), self%bz(self%nx))

      allocate (self%waves(self%nx, self%ny, self%nz))

      ! waves_set requires some of the preprocessed tdsops variables.
      call self%waves_set(xdirps, ydirps, zdirps)

   end subroutine base_init

   subroutine waves_set(self, xdirps, ydirps, zdirps)
      implicit none

      class(poisson_fft_t) :: self
      type(dirps_t), intent(in) :: xdirps, ydirps, zdirps

      complex(dp), allocatable, dimension(:) :: xkx, xk2, yky, yk2, zkz, zk2, &
                                                exs, eys, ezs

      integer :: nx, ny, nz
      real(dp) :: w, wp, rlexs, rleys, rlezs, &
                  xtt, ytt, ztt, xtt1, ytt1, ztt1, xt1, yt1, zt1
      complex(dp) :: xt2, yt2, zt2, xyzk

      integer :: i, j, k, ka, kb, ix, iy, iz


      nx = xdirps%n; ny = ydirps%n; nz = zdirps%n

      do i = 1, nx
         self%ax(i) = sin((i-1)*pi/nx)
         self%bx(i) = cos((i-1)*pi/nx)
      end do

      do i = 1, ny
         self%ay(i) = sin((i-1)*pi/ny)
         self%by(i) = cos((i-1)*pi/ny)
      end do

      do i = 1, nz
         self%az(i) = sin((i-1)*pi/nz)
         self%bz(i) = cos((i-1)*pi/nz)
      end do

      ! Now kxyz
      allocate(xkx(nx), xk2(nx), exs(nx))
      allocate(yky(ny), yk2(ny), eys(ny))
      allocate(zkz(nz), zk2(nz), ezs(nz))
      !xkx(:) = 0; xk2(:) = 0; yky(:) = 0; yk2(:) = 0; zkz(:) = 0; zk2(:) = 0

   end subroutine waves_set

end module m_poisson_fft
