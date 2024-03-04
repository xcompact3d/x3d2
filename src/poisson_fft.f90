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

   subroutine base_init(self, xdirps, ydirps, zdirps, sz)
      implicit none

      class(poisson_fft_t) :: self
      class(dirps_t), intent(in) :: xdirps, ydirps, zdirps
      integer, intent(in) :: sz

      self%nx = xdirps%n; self%ny = ydirps%n; self%nz = zdirps%n

      allocate (self%ax(self%nx), self%bx(self%nx))
      allocate (self%ay(self%nx), self%by(self%nx))
      allocate (self%az(self%nx), self%bz(self%nx))

      allocate (self%waves(sz, self%nx, (self%ny*(self%nz/2 + 1))/sz))

      ! waves_set requires some of the preprocessed tdsops variables.
      call self%waves_set(xdirps, ydirps, zdirps, sz)

   end subroutine base_init

   subroutine waves_set(self, xdirps, ydirps, zdirps, sz)
      implicit none

      class(poisson_fft_t) :: self
      type(dirps_t), intent(in) :: xdirps, ydirps, zdirps
      integer, intent(in) :: sz

      complex(dp), allocatable, dimension(:) :: xkx, xk2, yky, yk2, zkz, zk2, &
                                                exs, eys, ezs

      integer :: nx, ny, nz
      real(dp) :: w, wp, rlexs, rleys, rlezs, xtt, ytt, ztt, xt1, yt1, zt1
      complex(dp) :: xt2, yt2, zt2, xyzk

      integer :: i, j, ka, kb, ix, iy, iz


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
      xkx(:) = 0; xk2(:) = 0; yky(:) = 0; yk2(:) = 0; zkz(:) = 0; zk2(:) = 0

      ! periodic-x
      do i = 1, nx/2 + 1
         w = 2*pi*(i - 1)/nx
         wp = xdirps%stagder_v2p%a*2*xdirps%d*sin(0.5_dp*w) &
              + xdirps%stagder_v2p%b*2*xdirps%d*sin(1.5_dp*w)
         wp = wp/(1._dp + 2*xdirps%stagder_v2p%alpha*cos(w))

         xkx(i) = cmplx(1._dp, 1._dp, kind=dp)*(nx*wp/xdirps%L)
         exs(i) = cmplx(1._dp, 1._dp, kind=dp)*(nx*w/xdirps%L)
         xk2(i) = cmplx(1._dp, 1._dp, kind=dp)*(nx*wp/xdirps%L)**2
      end do
      do i = nx/2 + 2, nx
         xkx(i) = xkx(nx - i + 2)
         exs(i) = exs(nx - i + 2)
         xk2(i) = xk2(nx - i + 2)
      end do

      ! periodic-y
      do i = 1, ny/2 + 1
         w = 2*pi*(i - 1)/ny
         wp = ydirps%stagder_v2p%a*2*ydirps%d*sin(0.5_dp*w) &
              + ydirps%stagder_v2p%b*2*ydirps%d*sin(1.5_dp*w)
         wp = wp/(1._dp + 2*ydirps%stagder_v2p%alpha*cos(w))

         yky(i) = cmplx(1._dp, 1._dp, kind=dp)*(ny*wp/ydirps%L)
         eys(i) = cmplx(1._dp, 1._dp, kind=dp)*(ny*w/ydirps%L)
         yk2(i) = cmplx(1._dp, 1._dp, kind=dp)*(ny*wp/ydirps%L)**2
      end do
      do i = ny/2 + 2, ny
         yky(i) = yky(ny-i+2)
         eys(i) = eys(ny-i+2)
         yk2(i) = yk2(ny-i+2)
      end do

      ! periodic-z
      do i = 1, nz/2 + 1
         w = 2*pi*(i - 1)/nz
         wp = zdirps%stagder_v2p%a*2*zdirps%d*sin(0.5_dp*w) &
              + zdirps%stagder_v2p%b*2*zdirps%d*sin(1.5_dp*w)
         wp = wp/(1._dp + 2*zdirps%stagder_v2p%alpha*cos(w))

         zkz(i) = cmplx(1._dp, 1._dp, kind=dp)*(nz*wp/zdirps%L)
         ezs(i) = cmplx(1._dp, 1._dp, kind=dp)*(nz*w/zdirps%L)
         zk2(i) = cmplx(1._dp, 1._dp, kind=dp)*(nz*wp/zdirps%L)**2
      end do

      print*, 'set waves array'
      ! TODO: do loop ranges below are valid only for single rank runs
      do ka = 1, nz/2 + 1
         do kb = 1, ny/sz
            do j = 1, nx
               do i = 1, sz
                  ix = j; iy = (kb-1)*sz+i; iz = ka !+ xderps%nspz_st - 1
                  rlexs = real(exs(ix), kind=dp)*xdirps%d
                  rleys = real(eys(iy), kind=dp)*ydirps%d
                  rlezs = real(ezs(iz), kind=dp)*zdirps%d

                  xtt = 2*(xdirps%interpl_v2p%a*cos(rlexs*0.5_dp) &
                           + xdirps%interpl_v2p%b*cos(rlexs*1.5_dp) &
                           + xdirps%interpl_v2p%c*cos(rlexs*2.5_dp) &
                           + xdirps%interpl_v2p%d*cos(rlexs*3.5_dp))
                  ytt = 2*(ydirps%interpl_v2p%a*cos(rleys*0.5_dp) &
                           + ydirps%interpl_v2p%b*cos(rleys*1.5_dp) &
                           + ydirps%interpl_v2p%c*cos(rleys*2.5_dp) &
                           + ydirps%interpl_v2p%d*cos(rleys*3.5_dp))
                  ztt = 2*(zdirps%interpl_v2p%a*cos(rlezs*0.5_dp) &
                           + zdirps%interpl_v2p%b*cos(rlezs*1.5_dp) &
                           + zdirps%interpl_v2p%c*cos(rlezs*2.5_dp) &
                           + zdirps%interpl_v2p%d*cos(rlezs*3.5_dp))

                  xt1 = 1._dp + 2*xdirps%interpl_v2p%alpha*cos(rlexs)
                  yt1 = 1._dp + 2*ydirps%interpl_v2p%alpha*cos(rleys)
                  zt1 = 1._dp + 2*zdirps%interpl_v2p%alpha*cos(rlezs)

                  xt2 = xk2(ix)*(((ytt/yt1)*(ztt/zt1))**2)
                  yt2 = yk2(iy)*(((xtt/xt1)*(ztt/zt1))**2)
                  zt2 = zk2(iz)*(((xtt/xt1)*(ytt/yt1))**2)

                  xyzk = xt2 + yt2 + zt2
                  self%waves(i, j, ka + (kb - 1)*(nz/2 + 1)) = xyzk
               end do
            end do
         end do
      end do

   end subroutine waves_set

end module m_poisson_fft
