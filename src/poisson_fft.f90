module m_poisson_fft
  use m_allocator, only: field_t
  use m_common, only: dp, pi, CELL
  use m_tdsops, only: dirps_t
  use m_mesh, only: mesh_t, geo_t

  implicit none

  type, abstract :: poisson_fft_t
    !! FFT based Poisson solver
    !> Global dimensions
    integer :: nx_glob, ny_glob, nz_glob
    !> Local dimensions
    integer :: nx_loc, ny_loc, nz_loc
    !> Local dimensions in the permuted slabs
    integer :: nx_perm, ny_perm, nz_perm
    !> Local dimensions in the permuted slabs in spectral space
    integer :: nx_spec, ny_spec, nz_spec
    !> Offset in y direction in the permuted slabs in spectral space
    integer :: y_sp_st
    !> Local domain sized array storing the spectral equivalence constants
    complex(dp), allocatable, dimension(:, :, :) :: waves
    !> Wave numbers in x, y, and z
    real(dp), allocatable, dimension(:) :: ax, bx, ay, by, az, bz
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
    type(dirps_t), intent(in) :: xdirps, ydirps, zdirps

    integer :: dims(3)

    dims = mesh%get_global_dims(CELL)
    self%nx_glob = dims(1); self%ny_glob = dims(2); self%nz_glob = dims(3)
    dims = mesh%get_dims(CELL)
    self%nx_loc = dims(1); self%ny_loc = dims(2); self%nz_loc = dims(3)

    ! 1D decomposition along Z in real domain, and along Y in spectral space
    if (mesh%par%nproc_dir(1) /= 1) print *, 'nproc_dir in x-dir must be 1'
    if (mesh%par%nproc_dir(2) /= 1) print *, 'nproc_dir in y-dir must be 1'
    self%nx_perm = self%nx_loc/mesh%par%nproc_dir(2)
    self%ny_perm = self%ny_loc/mesh%par%nproc_dir(3)
    self%nz_perm = self%nz_glob
    self%nx_spec = self%nx_loc/2 + 1
    self%ny_spec = self%ny_perm
    self%nz_spec = self%nz_perm

    self%y_sp_st = (self%ny_loc/mesh%par%nproc_dir(3))*mesh%par%nrank_dir(3)

    allocate (self%ax(self%nx_glob), self%bx(self%nx_glob))
    allocate (self%ay(self%ny_glob), self%by(self%ny_glob))
    allocate (self%az(self%nz_glob), self%bz(self%nz_glob))

    ! FFT 3D transform halves the first index.
    allocate (self%waves(self%nx_spec, self%ny_spec, self%nz_spec))

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

    integer :: nx, ny, nz, ix, iy, iz
    real(dp) :: w, wp, rlexs, rleys, rlezs, xtt, ytt, ztt, xt1, yt1, zt1
    complex(dp) :: xt2, yt2, zt2, xyzk
    real(dp) :: d, L

    integer :: i, j, k

    nx = self%nx_glob; ny = self%ny_glob; nz = self%nz_glob

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

    do i = 1, self%nx_spec
      do j = 1, self%ny_spec
        do k = 1, self%nz_spec
          ix = i; iy = j + self%y_sp_st; iz = k
          rlexs = real(exs(ix), kind=dp)*geo%d(1)
          rleys = real(eys(iy), kind=dp)*geo%d(2)
          rlezs = real(ezs(iz), kind=dp)*geo%d(3)

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
          self%waves(i, j, k) = xyzk
        end do
      end do
    end do

  end subroutine waves_set

end module m_poisson_fft
