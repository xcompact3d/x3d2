module m_poisson_fft
  use m_allocator, only: field_t
  use m_common, only: dp, pi, VERT, CELL, DIR_X, DIR_Y, DIR_Z
  use m_tdsops, only: dirps_t
  use m_mesh, only: mesh_t, geo_t

  implicit none

  type, abstract :: poisson_fft_t
    !! FFT based Poisson solver
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

  subroutine base_init(self, mesh, xdirps, ydirps, zdirps)
    implicit none

    class(poisson_fft_t) :: self
    class(mesh_t), intent(in) :: mesh
    class(dirps_t), intent(in) :: xdirps, ydirps, zdirps

    self%nx = mesh%get_n(DIR_X, CELL)
    self%ny = mesh%get_n(DIR_Y, CELL)
    self%nz = mesh%get_n(DIR_Z, CELL)

    allocate (self%ax(self%nx), self%bx(self%nx))
    allocate (self%ay(self%ny), self%by(self%ny))
    allocate (self%az(self%nz), self%bz(self%nz))

    ! cuFFT 3D transform halves the first index.
    allocate (self%waves(self%nx/2 + 1, self%ny, self%nz))

    ! waves_set requires some of the preprocessed tdsops variables.
    call self%waves_set(mesh%geo, xdirps, ydirps, zdirps)

  end subroutine base_init

  subroutine waves_set(self, geo, xdirps, ydirps, zdirps)
    !! Spectral equivalence constants
    !!
    !! Ref. JCP 228 (2009), 5989–6015, Sec 4
    implicit none

    class(poisson_fft_t) :: self
    type(geo_t), intent(in) :: geo
    type(dirps_t), intent(in) :: xdirps, ydirps, zdirps

    complex(dp), allocatable, dimension(:) :: xkx, xk2, yky, yk2, zkz, zk2, &
                                              exs, eys, ezs

    integer :: nx, ny, nz
    real(dp) :: w, wp, rlexs, rleys, rlezs, xtt, ytt, ztt, xt1, yt1, zt1
    complex(dp) :: xt2, yt2, zt2, xyzk
    real(dp) :: d, L

    integer :: i, j, k

    nx = self%nx; ny = self%ny; nz = self%nz

    do i = 1, nx
      self%ax(i) = sin((i - 1)*pi/nx)
      self%bx(i) = cos((i - 1)*pi/nx)
    end do

    do i = 1, ny
      self%ay(i) = sin((i - 1)*pi/ny)
      self%by(i) = cos((i - 1)*pi/ny)
    end do

    do i = 1, nz
      self%az(i) = sin((i - 1)*pi/nz)
      self%bz(i) = cos((i - 1)*pi/nz)
    end do

    ! Now kxyz
    allocate (xkx(nx), xk2(nx), exs(nx))
    allocate (yky(ny), yk2(ny), eys(ny))
    allocate (zkz(nz), zk2(nz), ezs(nz))
    xkx(:) = 0; xk2(:) = 0; yky(:) = 0; yk2(:) = 0; zkz(:) = 0; zk2(:) = 0

    ! periodic-x
    d = geo%d(1)
    L = geo%L(1)
    do i = 1, nx/2 + 1
      w = 2*pi*(i - 1)/nx
      wp = xdirps%stagder_v2p%a*2*d*sin(0.5_dp*w) &
           + xdirps%stagder_v2p%b*2*d*sin(1.5_dp*w)
      wp = wp/(1._dp + 2*xdirps%stagder_v2p%alpha*cos(w))

      xkx(i) = cmplx(1._dp, 1._dp, kind=dp)*(nx*wp/L)
      exs(i) = cmplx(1._dp, 1._dp, kind=dp)*(nx*w/L)
      xk2(i) = cmplx(1._dp, 1._dp, kind=dp)*(nx*wp/L)**2
    end do
    do i = nx/2 + 2, nx
      xkx(i) = xkx(nx - i + 2)
      exs(i) = exs(nx - i + 2)
      xk2(i) = xk2(nx - i + 2)
    end do

    ! periodic-y
    d = geo%d(2)
    L = geo%L(2)
    do i = 1, ny/2 + 1
      w = 2*pi*(i - 1)/ny
      wp = ydirps%stagder_v2p%a*2*d*sin(0.5_dp*w) &
           + ydirps%stagder_v2p%b*2*d*sin(1.5_dp*w)
      wp = wp/(1._dp + 2*ydirps%stagder_v2p%alpha*cos(w))

      yky(i) = cmplx(1._dp, 1._dp, kind=dp)*(ny*wp/L)
      eys(i) = cmplx(1._dp, 1._dp, kind=dp)*(ny*w/L)
      yk2(i) = cmplx(1._dp, 1._dp, kind=dp)*(ny*wp/L)**2
    end do
    do i = ny/2 + 2, ny
      yky(i) = yky(ny - i + 2)
      eys(i) = eys(ny - i + 2)
      yk2(i) = yk2(ny - i + 2)
    end do

    ! periodic-z
    d = geo%d(3)
    L = geo%L(3)
    do i = 1, nz/2 + 1
      w = 2*pi*(i - 1)/nz
      wp = zdirps%stagder_v2p%a*2*d*sin(0.5_dp*w) &
           + zdirps%stagder_v2p%b*2*d*sin(1.5_dp*w)
      wp = wp/(1._dp + 2*zdirps%stagder_v2p%alpha*cos(w))

      zkz(i) = cmplx(1._dp, 1._dp, kind=dp)*(nz*wp/L)
      ezs(i) = cmplx(1._dp, 1._dp, kind=dp)*(nz*w/L)
      zk2(i) = cmplx(1._dp, 1._dp, kind=dp)*(nz*wp/L)**2
    end do
    do i = nz/2 + 2, nz
      zkz(i) = zkz(nz - i + 2)
      ezs(i) = ezs(nz - i + 2)
      zk2(i) = zk2(nz - i + 2)
    end do

    print *, 'waves array is correctly set only for a single rank run'
    ! TODO: do loop ranges below are valid only for single rank runs
    do i = 1, nx/2 + 1
      do j = 1, ny
        do k = 1, nz
          rlexs = real(exs(i), kind=dp)*geo%d(1)
          rleys = real(eys(j), kind=dp)*geo%d(2)
          rlezs = real(ezs(k), kind=dp)*geo%d(3)

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

          xt2 = xk2(i)*(((ytt/yt1)*(ztt/zt1))**2)
          yt2 = yk2(j)*(((xtt/xt1)*(ztt/zt1))**2)
          zt2 = zk2(k)*(((xtt/xt1)*(ytt/yt1))**2)

          xyzk = xt2 + yt2 + zt2
          self%waves(i, j, k) = xyzk
        end do
      end do
    end do

  end subroutine waves_set

end module m_poisson_fft
