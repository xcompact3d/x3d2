module m_poisson_fft
  use m_common, only: dp, pi, CELL
  use m_field, only: field_t
  use m_mesh, only: mesh_t, geo_t
  use m_tdsops, only: dirps_t

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
    !> Offset in y and z directions in the permuted slabs in spectral space
    integer :: x_sp_st, y_sp_st, z_sp_st
    !> Local domain sized array storing the spectral equivalence constants
    complex(dp), allocatable, dimension(:, :, :) :: waves
    !> Wave numbers in x, y, and z
    real(dp), allocatable, dimension(:) :: ax, bx, ay, by, az, bz
    !> Periodicity in x, y, and z
    logical :: periodic_x, periodic_y, periodic_z
    !> Procedure pointer to BC specific poisson solvers
    procedure(poisson_xxx), pointer :: poisson => null()
  contains
    procedure(fft_forward), deferred :: fft_forward
    procedure(fft_backward), deferred :: fft_backward
    procedure(fft_postprocess), deferred :: fft_postprocess_000
    procedure(fft_postprocess), deferred :: fft_postprocess_010
    procedure(field_process), deferred :: enforce_periodicity_y
    procedure(field_process), deferred :: undo_periodicity_y
    procedure :: solve_poisson
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

  abstract interface
    subroutine poisson_xxx(self, f, temp)
      import :: poisson_fft_t
      import :: field_t

      class(poisson_fft_t) :: self
      class(field_t), intent(inout) :: f, temp
    end subroutine poisson_xxx

    subroutine field_process(self, f_out, f_in)
      import :: poisson_fft_t
      import :: field_t

      class(poisson_fft_t) :: self
      class(field_t), intent(inout) :: f_out
      class(field_t), intent(in) :: f_in
    end subroutine field_process
  end interface

contains

  subroutine base_init(self, mesh, xdirps, ydirps, zdirps, n_spec, n_sp_st)
    implicit none

    class(poisson_fft_t) :: self
    type(mesh_t), intent(in) :: mesh
    type(dirps_t), intent(in) :: xdirps, ydirps, zdirps
    integer, dimension(3), intent(in) :: n_spec ! Size of the spectral pencil
    integer, dimension(3), intent(in) :: n_sp_st ! Offset of the spectral pencil

    integer :: dims(3)

    ! Decomposition is in y- and z-directions
    if (mesh%par%nproc_dir(1) /= 1) print *, 'nproc_dir in x-dir must be 1'

    dims = mesh%get_global_dims(CELL)
    self%nx_glob = dims(1); self%ny_glob = dims(2); self%nz_glob = dims(3)
    dims = mesh%get_dims(CELL)
    self%nx_loc = dims(1); self%ny_loc = dims(2); self%nz_loc = dims(3)

    self%nx_spec = n_spec(1)
    self%ny_spec = n_spec(2)
    self%nz_spec = n_spec(3)

    self%periodic_x = mesh%grid%periodic_BC(1)
    self%periodic_y = mesh%grid%periodic_BC(2)
    self%periodic_z = mesh%grid%periodic_BC(3)

    self%x_sp_st = n_sp_st(1)
    self%y_sp_st = n_sp_st(2)
    self%z_sp_st = n_sp_st(3)

    allocate (self%ax(self%nx_glob), self%bx(self%nx_glob))
    allocate (self%ay(self%ny_glob), self%by(self%ny_glob))
    allocate (self%az(self%nz_glob), self%bz(self%nz_glob))

    allocate (self%waves(self%nx_spec, self%ny_spec, self%nz_spec))

    ! waves_set requires some of the preprocessed tdsops variables.
    call self%waves_set(mesh%geo, xdirps, ydirps, zdirps)

    ! use correct procedure based on BCs
    if (self%periodic_x .and. self%periodic_y .and. self%periodic_z) then
      self%poisson => poisson_000
    else if (self%periodic_x .and. (.not. self%periodic_y) &
             .and. (self%periodic_z)) then
      self%poisson => poisson_010
    end if
  end subroutine base_init

  subroutine solve_poisson(self, f, temp)
    implicit none

    class(poisson_fft_t) :: self
    class(field_t), intent(inout) :: f, temp

    call self%poisson(f, temp)

  end subroutine solve_poisson

  subroutine poisson_000(self, f, temp)
    implicit none

    class(poisson_fft_t) :: self
    class(field_t), intent(inout) :: f, temp

    call self%fft_forward(f)
    call self%fft_postprocess_000
    call self%fft_backward(f)

  end subroutine poisson_000

  subroutine poisson_010(self, f, temp)
    implicit none

    class(poisson_fft_t) :: self
    class(field_t), intent(inout) :: f, temp

    call self%enforce_periodicity_y(temp, f)

    call self%fft_forward(temp)
    call self%fft_postprocess_010
    call self%fft_backward(temp)

    call self%undo_periodicity_y(f, temp)

  end subroutine poisson_010

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
    real(dp) :: rlexs, rleys, rlezs, xtt, ytt, ztt, xt1, yt1, zt1
    complex(dp) :: xt2, yt2, zt2, xyzk
    real(dp) :: L_x, L_y, L_z, d_x, d_y, d_z

    integer :: i, j, k

    nx = self%nx_glob; ny = self%ny_glob; nz = self%nz_glob
    L_x = geo%L(1); L_y = geo%L(2); L_z = geo%L(3)
    d_x = geo%d(1); d_y = geo%d(2); d_z = geo%d(3)

    ! Now kxyz
    allocate (xkx(nx), xk2(nx), exs(nx))
    allocate (yky(ny), yk2(ny), eys(ny))
    allocate (zkz(nz), zk2(nz), ezs(nz))
    xkx(:) = 0; xk2(:) = 0; yky(:) = 0; yk2(:) = 0; zkz(:) = 0; zk2(:) = 0

    call wave_numbers( &
      self%ax, self%bx, xkx, exs, xk2, nx, L_x, d_x, self%periodic_x, &
      xdirps%stagder_v2p%a, xdirps%stagder_v2p%b, xdirps%stagder_v2p%alpha &
      )

    call wave_numbers( &
      self%ay, self%by, yky, eys, yk2, ny, L_y, d_y, self%periodic_y, &
      ydirps%stagder_v2p%a, ydirps%stagder_v2p%b, ydirps%stagder_v2p%alpha &
      )

    call wave_numbers( &
      self%az, self%bz, zkz, ezs, zk2, nz, L_z, d_z, self%periodic_z, &
      zdirps%stagder_v2p%a, zdirps%stagder_v2p%b, zdirps%stagder_v2p%alpha &
      )

    if (self%periodic_z) then
      ! poisson 000, 100, 010, 110
      do k = 1, self%nz_spec
        do j = 1, self%ny_spec
          do i = 1, self%nx_spec
            ix = i + self%x_sp_st; iy = j + self%y_sp_st; iz = k + self%z_sp_st
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
    else if (.not. (self%periodic_x .and. self%periodic_y .and. &
                    self%periodic_z)) then
      ! poisson 111
      error stop 'No support for all non-periodic BCs yet!'
    else
      ! poisson 001, 011, 101
      error stop 'FFT Poisson solver does not support specified BCs!'
    end if

  end subroutine waves_set

  subroutine wave_numbers(a, b, k, e, k2, n, L, d, periodic, c_a, c_b, c_alpha)
    implicit none

    real(dp), dimension(:), intent(out) :: a, b
    complex(dp), dimension(:), intent(out) :: k, e, k2
    integer, intent(in) :: n
    real(dp), intent(in) :: c_a, c_b, c_alpha, L, d
    logical, intent(in) :: periodic

    real(dp) :: w, wp
    integer :: i

    do i = 1, n
      if (periodic) then
        a(i) = sin((i - 1)*pi/n)
        b(i) = cos((i - 1)*pi/n)
      else
        a(i) = sin((i - 1)*pi/2/n)
        b(i) = cos((i - 1)*pi/2/n)
      end if
    end do

    if (periodic) then
      do i = 1, n/2 + 1
        w = 2*pi*(i - 1)/n
        wp = c_a*2*d*sin(0.5_dp*w) + c_b*2*d*sin(1.5_dp*w)
        wp = wp/(1._dp + 2*c_alpha*cos(w))

        k(i) = cmplx(1._dp, 1._dp, kind=dp)*(n*wp/L)
        e(i) = cmplx(1._dp, 1._dp, kind=dp)*(n*w/L)
        k2(i) = cmplx(1._dp, 1._dp, kind=dp)*(n*wp/L)**2
      end do
      do i = n/2 + 2, n
        k(i) = k(n - i + 2)
        e(i) = e(n - i + 2)
        k2(i) = k2(n - i + 2)
      end do
    else
      do i = 1, n
        w = pi*(i - 1)/n
        wp = c_a*2*d*sin(0.5_dp*w) + c_b*2*d*sin(1.5_dp*w)
        wp = wp/(1._dp + 2*c_alpha*cos(w))

        k(i) = cmplx(1._dp, 1._dp, kind=dp)*(n*wp/L)
        e(i) = cmplx(1._dp, 1._dp, kind=dp)*(n*w/L)
        k2(i) = cmplx(1._dp, 1._dp, kind=dp)*(n*wp/L)**2
      end do
    end if

  end subroutine wave_numbers

end module m_poisson_fft
