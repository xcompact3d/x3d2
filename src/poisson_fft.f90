module m_poisson_fft
  use m_common, only: dp, pi, CELL
  use m_field, only: field_t
  use m_mesh, only: mesh_t, geo_t, FFT_BC_000, FFT_BC_010, FFT_BC_100, &
                    FFT_BC_110, FFT_BC_UNSUPPORTED
  use m_tdsops, only: dirps_t
  use m_logger, only: log_warn, abort_run

  implicit none

  private :: fill_stretched_diagonals

  !> Stretching matrix kinds of fill_stretched_diagonals
  integer, parameter, private :: STRETCH_MAT_BOTTOM = 1, STRETCH_MAT_ODD = 2, &
                                 STRETCH_MAT_EVEN = 3

  type, abstract :: poisson_fft_t
    !! FFT based Poisson solver
    type(mesh_t), pointer :: mesh => null()
    !> Global dimensions
    integer :: nx_glob, ny_glob, nz_glob
    !> Local dimensions
    integer :: nx_loc, ny_loc, nz_loc
    !> Local dimensions in the permuted slabs
    integer :: nx_perm, ny_perm, nz_perm
    !> Local dimensions in the permuted slabs in spectral space
    integer :: nx_spec, ny_spec, nz_spec
    !> Offset per spectral dimension (dim1, dim2, dim3) in spectral space
    integer :: sp_st(3)
    !> Local domain sized array storing the spectral equivalence constants
    complex(dp), allocatable, dimension(:, :, :) :: waves
    !> Wave numbers in x, y, and z
    real(dp), allocatable, dimension(:) :: ax, bx, ay, by, az, bz
    !> Wave numbers in x, y, and z
    complex(dp), allocatable, dimension(:) :: kx, ky, kz, exs, eys, ezs, &
                                              k2x, k2y, k2z
    !> Staggered grid transformation
    real(dp), allocatable, dimension(:) :: trans_x_re, trans_x_im, &
                                           trans_y_re, trans_y_im, &
                                           trans_z_re, trans_z_im
    !> Periodicity in x, y, and z
    logical :: periodic_x, periodic_y, periodic_z, &
               stretched_y = .false., stretched_y_sym
    !> Boundary condition case, see poisson_bc_code. The backend's init sets
    !> it first, ahead of base_init.
    character(len=3) :: bc_code
    !> Stretching operator matrices
    real(dp), allocatable, dimension(:, :, :, :) :: a_odd_re, a_odd_im, &
                                                    a_even_re, a_even_im, &
                                                    a_re, a_im
    !> lowmem option, only used in CUDA backend
    logical :: lowmem = .false.
  contains
    !! Each backend extends this type once per boundary condition case
    !! (000, 010, 100, 110), and each case type supplies its own transforms,
    !! spectral postprocess and periodicity reorder.
    procedure(fft_forward), deferred :: fft_forward
    procedure(fft_backward), deferred :: fft_backward
    procedure(fft_postprocess), deferred :: fft_postprocess
    procedure :: enforce_periodicity
    procedure :: undo_periodicity
    procedure :: base_init
    procedure :: solve_poisson
    procedure :: stretching_matrix
    procedure :: waves_set
    procedure :: get_km
    procedure :: get_km_re
    procedure :: get_km_im
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

  pure function poisson_bc_code(periodic) result(bc_code)
    !! Names the boundary condition case of the FFT Poisson solver, one
    !! digit per direction in x, y, z order: 0 for periodic, 1 otherwise.
    !! The 000 case is periodic everywhere, 110 is non-periodic in x and y.
    implicit none

    logical, intent(in) :: periodic(3)
    character(len=3) :: bc_code

    integer :: dir

    do dir = 1, 3
      if (periodic(dir)) then
        bc_code(dir:dir) = '0'
      else
        bc_code(dir:dir) = '1'
      end if
    end do
  end function poisson_bc_code

  pure function fft_bc_class(bc_code) result(fft_bc)
    !! The FFT_BC_* class of a poisson_bc_code case, or FFT_BC_UNSUPPORTED.
    !! The CUDA solver passes it to the helpers it shares with x3d2-memcheck,
    !! so it classifies through poisson_bc_code alone. x3d2-memcheck has no
    !! mesh and takes the class from m_mesh's classify_boundaries;
    !! test_poisson_bc_code checks that the two agree.
    implicit none

    character(len=3), intent(in) :: bc_code
    integer :: fft_bc

    select case (bc_code)
    case ('000')
      fft_bc = FFT_BC_000
    case ('010')
      fft_bc = FFT_BC_010
    case ('100')
      fft_bc = FFT_BC_100
    case ('110')
      fft_bc = FFT_BC_110
    case default
      fft_bc = FFT_BC_UNSUPPORTED
    end select
  end function fft_bc_class

  pure function spectral_dirs(bc_code) result(dirs)
    !! The physical direction held in each dimension of the spectral array.
    !! The r2c transform runs along dim1, so a case with x non-periodic
    !! permutes a periodic direction into dim1 ahead of the FFT: 100 swaps
    !! x and y, 110 brings z to the front.
    implicit none

    character(len=3), intent(in) :: bc_code
    integer :: dirs(3)

    select case (bc_code)
    case ('100')
      dirs = [2, 1, 3]
    case ('110')
      dirs = [3, 1, 2]
    case default
      dirs = [1, 2, 3]
    end select
  end function spectral_dirs

  subroutine base_init(self, mesh, xdirps, ydirps, zdirps, n_spec, n_sp_st)
    !! Sets up the state every case shares. The backend's init calls it
    !! once it has set bc_code, which the spectral extents it passes in
    !! already depend on. The supported BC codes are the cases of the
    !! backend factories, create_omp_poisson_fft and create_cuda_poisson_fft.
    implicit none
    class(poisson_fft_t) :: self
    type(mesh_t), intent(in), target :: mesh
    type(dirps_t), intent(in) :: xdirps, ydirps, zdirps
    integer, dimension(3), intent(in) :: n_spec ! Size of the spectral pencil
    integer, dimension(3), intent(in) :: n_sp_st ! Offset of the spectral pencil

    integer :: dims(3)
    self%mesh => mesh
    ! Decomposition is in y- and z-directions
    if (mesh%par%nproc_dir(1) /= 1) &
      call log_warn('nproc_dir in x-dir must be 1')

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

    self%sp_st = n_sp_st

    allocate (self%ax(self%nx_glob), self%bx(self%nx_glob))
    allocate (self%ay(self%ny_glob), self%by(self%ny_glob))
    allocate (self%az(self%nz_glob), self%bz(self%nz_glob))

    allocate (self%kx(self%nx_glob), self%k2x(self%nx_glob))
    allocate (self%ky(self%ny_glob), self%k2y(self%ny_glob))
    allocate (self%kz(self%nz_glob), self%k2z(self%nz_glob))
    allocate (self%exs(self%nx_glob))
    allocate (self%eys(self%ny_glob))
    allocate (self%ezs(self%nz_glob))

    allocate (self%trans_x_re(self%nx_spec), self%trans_x_im(self%nx_spec))
    allocate (self%trans_y_re(self%ny_spec), self%trans_y_im(self%ny_spec))
    allocate (self%trans_z_re(self%nz_spec), self%trans_z_im(self%nz_spec))

    allocate (self%waves(self%nx_spec, self%ny_spec, self%nz_spec))

    ! waves_set requires some of the preprocessed tdsops variables.
    call self%waves_set(mesh%geo, xdirps, ydirps, zdirps)

    if (mesh%geo%stretched(1) .or. mesh%geo%stretched(3)) then
      call abort_run('FFT based Poisson solver does not support stretching in &
                     &x- or z-directions!')
    end if

    ! stretching requires some coefficients matrices, which only the 010
    ! case supports
    if (mesh%geo%stretched(2)) then
      if (self%bc_code == '010') then
        self%stretched_y = .true.
        call self%stretching_matrix(mesh%geo, xdirps, ydirps, zdirps)
      else
        call abort_run('FFT based Poisson solver supports stretching in &
                       &y-direction only with periodic x and z and &
                       &non-periodic y!')
      end if
    end if
  end subroutine base_init

  subroutine solve_poisson(self, f, temp)
    !! Gathers each non-periodic line into a periodic one, so that a plain
    !! FFT can stand in for the cosine transform, solves in spectral space,
    !! and scatters the result back. The 000 case has nothing to reorder and
    !! overrides this to transform f in place.
    implicit none

    class(poisson_fft_t) :: self
    class(field_t), intent(inout) :: f, temp

    call self%enforce_periodicity(temp, f)

    call self%fft_forward(temp)
    call self%fft_postprocess
    call self%fft_backward(temp)

    call self%undo_periodicity(f, temp)

  end subroutine solve_poisson

  subroutine enforce_periodicity(self, f_out, f_in)
    !! Reorder hook of solve_poisson, overridden by every case that has a
    !! non-periodic direction.
    implicit none

    class(poisson_fft_t) :: self
    class(field_t), intent(inout) :: f_out
    class(field_t), intent(in) :: f_in

    call abort_run('enforce_periodicity is not used by this Poisson BC &
                   &case.')
  end subroutine enforce_periodicity

  subroutine undo_periodicity(self, f_out, f_in)
    !! Inverse of enforce_periodicity.
    implicit none

    class(poisson_fft_t) :: self
    class(field_t), intent(inout) :: f_out
    class(field_t), intent(in) :: f_in

    call abort_run('undo_periodicity is not used by this Poisson BC case.')
  end subroutine undo_periodicity

  subroutine stretching_matrix(self, geo, xdirps, ydirps, zdirps)
    !! Stretching necessitates a special operation in spectral space.
    !! The coefficients for the operation are stored in matrix form.
    !!
    !! Ref. JCP 228 (2009), 5989–6015, Sec 5
    implicit none

    class(poisson_fft_t) :: self
    type(geo_t), intent(in) :: geo
    type(dirps_t), intent(in) :: xdirps, ydirps, zdirps

    real(dp) :: temp, a0, a1
    integer :: i, j, k, ix, iz

    do i = 1, self%nx_spec
      temp = real(self%exs(i), kind=dp)*geo%d(1)
      self%trans_x_re(i) = 2*(xdirps%interpl_v2p%a*cos(temp*0.5_dp) &
                              + xdirps%interpl_v2p%b*cos(temp*1.5_dp) &
                              + xdirps%interpl_v2p%c*cos(temp*2.5_dp) &
                              + xdirps%interpl_v2p%d*cos(temp*3.5_dp)) &
                           /(1._dp + 2*xdirps%interpl_v2p%alpha*cos(temp))
      self%trans_x_im(i) = self%trans_x_re(i)
    end do

    do j = 1, self%ny_spec
      temp = real(self%eys(j), kind=dp)*geo%d(2)
      self%trans_y_re(j) = 2*(ydirps%interpl_v2p%a*cos(temp*0.5_dp) &
                              + ydirps%interpl_v2p%b*cos(temp*1.5_dp) &
                              + ydirps%interpl_v2p%c*cos(temp*2.5_dp) &
                              + ydirps%interpl_v2p%d*cos(temp*3.5_dp)) &
                           /(1._dp + 2*ydirps%interpl_v2p%alpha*cos(temp))
      self%trans_y_im(j) = self%trans_y_re(j)
    end do

    do k = 1, self%nz_spec
      temp = real(self%ezs(k), kind=dp)*geo%d(3)
      self%trans_z_re(k) = 2*(zdirps%interpl_v2p%a*cos(temp*0.5_dp) &
                              + zdirps%interpl_v2p%b*cos(temp*1.5_dp) &
                              + zdirps%interpl_v2p%c*cos(temp*2.5_dp) &
                              + zdirps%interpl_v2p%d*cos(temp*3.5_dp)) &
                           /(1._dp + 2*zdirps%interpl_v2p%alpha*cos(temp))
      self%trans_z_im(k) = self%trans_z_re(k)
    end do

    if (trim(geo%stretching(2)) == 'bottom') then
      self%stretched_y_sym = .false.
      allocate (self%a_re(self%nx_spec, self%ny_spec, self%nz_spec, 5))
      allocate (self%a_im(self%nx_spec, self%ny_spec, self%nz_spec, 5))

      self%a_re = 0._dp
      self%a_im = 0._dp

      a0 = (geo%alpha(2)/pi + 1._dp/(2*pi*geo%beta(2)))*geo%L(2)
      a1 = -1._dp/(4*pi*geo%beta(2))*geo%L(2)

      call fill_stretched_diagonals(self, self%a_re, self%a_im, &
                                    STRETCH_MAT_BOTTOM, a0, a1)

      ! tweak the matrix to make it not singular
      self%a_re(1, 1, 1, 3) = 1._dp; self%a_im(1, 1, 1, 3) = 1._dp
      self%a_re(1, 1, 1, 4) = 0; self%a_im(1, 1, 1, 4) = 0
      self%a_re(1, 1, 1, 5) = 0; self%a_im(1, 1, 1, 5) = 0
    else
      self%stretched_y_sym = .true.
      allocate (self%a_odd_re(self%nx_spec, self%ny_spec/2, self%nz_spec, 5))
      allocate (self%a_odd_im(self%nx_spec, self%ny_spec/2, self%nz_spec, 5))
      allocate (self%a_even_re(self%nx_spec, self%ny_spec/2, self%nz_spec, 5))
      allocate (self%a_even_im(self%nx_spec, self%ny_spec/2, self%nz_spec, 5))

      self%a_odd_re = 0._dp
      self%a_odd_im = 0._dp
      self%a_even_re = 0._dp
      self%a_even_im = 0._dp

      a0 = (geo%alpha(2)/pi + 1._dp/(2*pi*geo%beta(2)))*geo%L(2)
      select case (trim(geo%stretching(2)))
      case ('centred')
        a1 = 1._dp/(4*pi*geo%beta(2))*geo%L(2)
      case ('top-bottom')
        a1 = -1._dp/(4*pi*geo%beta(2))*geo%L(2)
      case default
        a1 = 0._dp
      end select

      call fill_stretched_diagonals(self, self%a_odd_re, self%a_odd_im, &
                                    STRETCH_MAT_ODD, a0, a1)
      call fill_stretched_diagonals(self, self%a_even_re, self%a_even_im, &
                                    STRETCH_MAT_EVEN, a0, a1)

      do k = 1, self%nz_spec
        do i = 1, self%nx_spec
          ix = i + self%sp_st(1); iz = k + self%sp_st(3)
          if (get_real(self%k2x(ix)) < 1e-15 &
              .and. get_real(self%k2z(iz)) < 1e-15) then
            self%a_odd_re(i, 1, k, 3) = 1._dp
            self%a_odd_im(i, 1, k, 3) = 1._dp
            self%a_odd_re(i, 1, k, 4) = 0
            self%a_odd_im(i, 1, k, 4) = 0
            self%a_odd_re(i, 1, k, 5) = 0
            self%a_odd_im(i, 1, k, 5) = 0
          end if
        end do
      end do
    end if

  end subroutine stretching_matrix

  subroutine fill_stretched_diagonals(self, a_re, a_im, matrix_kind, a0, a1)
    !! Fill the five diagonals of one stretching matrix. Row j of the matrix
    !! couples the cosine modes iy_scale*iy + iy_off + (-2:2)*shift, with
    !! iy = j + sp_st(2). The caller zeroes a_re and a_im beforehand.
    !! The matrix kind sets nrows, iy_scale, iy_off and shift: 'bottom' takes
    !! all the modes, 'odd' and 'even' take every second mode.
    !! The kind also selects the first and last rows: 'bottom' and 'odd' double
    !! the off-diagonals of row 1, as the ky=0 cosine mode has half the weight
    !! of the other modes. 'even' has its own edge coefficients and no doubling.
    !! No kind has a diagonal + 1 on its last row, as there is no mode above.
    !! Leave it zero, so the loop never reads ky past ny_glob.
    implicit none

    class(poisson_fft_t), intent(in) :: self
    real(dp), intent(inout) :: a_re(:, :, :, :), a_im(:, :, :, :)
    integer, intent(in) :: matrix_kind
    real(dp), intent(in) :: a0, a1

    real(dp) :: c1, c2, c_first
    complex(dp) :: km_a1
    integer :: i, j, k, ix, iy, iz, r, nup1
    integer :: nrows, iy_scale, iy_off, shift
    logical :: is_even

    select case (matrix_kind)
    case (STRETCH_MAT_BOTTOM)
      nrows = self%ny_spec; iy_scale = 1; iy_off = 0; shift = 1
    case (STRETCH_MAT_ODD)
      nrows = self%ny_spec/2; iy_scale = 2; iy_off = -1; shift = 2
    case (STRETCH_MAT_EVEN)
      nrows = self%ny_spec/2; iy_scale = 2; iy_off = 0; shift = 2
    case default
      call abort_run('fill_stretched_diagonals: unknown matrix kind')
    end select

    is_even = matrix_kind == STRETCH_MAT_EVEN
    ! First-row coefficient of the even diagonal + 1. Compute it outside the loop
    ! to keep the rounding of the old code (no fused multiply-add in the loop).
    c_first = a0*a1 - a1*a1
    nup1 = nrows - 1

    do k = 1, self%nz_spec
      do j = 1, nrows
        do i = 1, self%nx_spec
          ix = i + self%sp_st(1); iy = j + self%sp_st(2); iz = k + self%sp_st(3)
          r = iy_scale*iy + iy_off

          ! diagonal
          c1 = a0*a0
          c2 = a1*a1
          if (iy == 1) then
            if (is_even) c1 = a0*a0 - a1*a1
            km_a1 = self%get_km(ix, r + shift, iz)
          else if (iy == nrows) then
            if (is_even) c1 = (a0 + a1)*(a0 + a1)
            km_a1 = self%get_km(ix, r - shift, iz)
          else
            km_a1 = self%get_km(ix, r - shift, iz) &
                    + self%get_km(ix, r + shift, iz)
          end if
          a_re(i, j, k, 3) = &
            -(get_real(self%kx(ix))*self%trans_y_re(r) &
              *self%trans_z_re(iz))**2 &
            - (get_real(self%kz(iz))*self%trans_y_re(r) &
               *self%trans_x_re(ix))**2 &
            - c1*self%get_km_re(ix, r, iz)**2 &
            - c2*self%get_km_re(ix, r, iz) &
            *get_real(km_a1)
          a_im(i, j, k, 3) = &
            -(get_imag(self%kx(ix))*self%trans_y_im(r) &
              *self%trans_z_im(iz))**2 &
            - (get_imag(self%kz(iz))*self%trans_y_im(r) &
               *self%trans_x_im(ix))**2 &
            - c1*self%get_km_im(ix, r, iz)**2 &
            - c2*self%get_km_im(ix, r, iz) &
            *get_imag(km_a1)

          ! diagonal + 1
          if (j <= nup1) then
            c1 = a0*a1
            c2 = a0*a1
            if (iy == 1) then
              if (is_even) then
                c1 = c_first
              else
                c1 = 2*a0*a1
                c2 = 2*a0*a1
              end if
            else if (is_even .and. iy == nrows - 1) then
              c2 = (a0 + a1)*a1
            end if
            a_re(i, j, k, 4) = &
              c1*(self%get_km_re(ix, r, iz) &
                  *self%get_km_re(ix, r + shift, iz)) &
              + c2*self%get_km_re(ix, r + shift, iz)**2
            a_im(i, j, k, 4) = &
              c1*(self%get_km_im(ix, r, iz) &
                  *self%get_km_im(ix, r + shift, iz)) &
              + c2*self%get_km_im(ix, r + shift, iz)**2
          end if

          ! diagonal + 2
          if (j <= nrows - 2) then
            c1 = a1*a1
            if (iy == 1 .and. .not. is_even) c1 = 2*a1*a1
            a_re(i, j, k, 5) = &
              -(c1*self%get_km_re(ix, r + shift, iz) &
                *self%get_km_re(ix, r + 2*shift, iz))
            a_im(i, j, k, 5) = &
              -(c1*self%get_km_im(ix, r + shift, iz) &
                *self%get_km_im(ix, r + 2*shift, iz))
          end if

          ! diagonal - 1
          if (j >= 2) then
            c1 = a0*a1
            c2 = a0*a1
            if (is_even .and. iy == 2) then
              c2 = (a0 + a1)*a1
            else if (is_even .and. iy == nrows) then
              c1 = (a0 + a1)*a1
            end if
            a_re(i, j, k, 2) = &
              c1*(self%get_km_re(ix, r, iz) &
                  *self%get_km_re(ix, r - shift, iz)) &
              + c2*self%get_km_re(ix, r - shift, iz)**2
            a_im(i, j, k, 2) = &
              c1*(self%get_km_im(ix, r, iz) &
                  *self%get_km_im(ix, r - shift, iz)) &
              + c2*self%get_km_im(ix, r - shift, iz)**2
          end if

          ! diagonal - 2
          if (j >= 3) then
            a_re(i, j, k, 1) = &
              -(a1*a1*self%get_km_re(ix, r - shift, iz) &
                *self%get_km_re(ix, r - 2*shift, iz))
            a_im(i, j, k, 1) = &
              -(a1*a1*self%get_km_im(ix, r - shift, iz) &
                *self%get_km_im(ix, r - 2*shift, iz))
          end if
        end do
      end do
    end do
  end subroutine fill_stretched_diagonals

  subroutine waves_set(self, geo, xdirps, ydirps, zdirps)
  !! Spectral equivalence constants
  !!
  !! Ref. JCP 228 (2009), 5989–6015, Sec 4
    implicit none

    class(poisson_fft_t) :: self
    type(geo_t), intent(in) :: geo
    type(dirps_t), intent(in) :: xdirps, ydirps, zdirps

    integer :: i, j, k, ix, iy, iz
    real(dp) :: rlexs, rleys, rlezs, xtt, ytt, ztt, xt1, yt1, zt1
    complex(dp) :: xt2, yt2, zt2, xyzk

    call wave_numbers( &
      self%ax, self%bx, self%kx, self%exs, self%k2x, &
      self%nx_glob, geo%L(1), geo%d(1), self%periodic_x, &
      xdirps%stagder_v2p%a, xdirps%stagder_v2p%b, xdirps%stagder_v2p%alpha &
      )

    call wave_numbers( &
      self%ay, self%by, self%ky, self%eys, self%k2y, &
      self%ny_glob, geo%L(2), geo%d(2), self%periodic_y, &
      ydirps%stagder_v2p%a, ydirps%stagder_v2p%b, ydirps%stagder_v2p%alpha &
      )

    call wave_numbers( &
      self%az, self%bz, self%kz, self%ezs, self%k2z, &
      self%nz_glob, geo%L(3), geo%d(3), self%periodic_z, &
      zdirps%stagder_v2p%a, zdirps%stagder_v2p%b, zdirps%stagder_v2p%alpha &
      )

    ! The spectral array is laid out as spectral_dirs gives, e.g. for 110
    ! it is (nz/2 + 1, nx, ny): dim1 indexes the z r2c modes, dim2 the x
    ! modes and dim3 the y modes. Each layout keeps its own loop nest with
    ! the mode numbers written out. Mapping the indices through
    ! spectral_dirs in one shared loop stops the compiler from vectorising
    ! the inner loop and hoisting the other two directions out of it, which
    ! changes waves by a few ulp under -ffast-math.
    select case (self%bc_code)
    case ('110')
      do k = 1, self%nz_spec
        do j = 1, self%ny_spec
          do i = 1, self%nx_spec
            iz = i + self%sp_st(1)
            ix = j + self%sp_st(2)
            iy = k + self%sp_st(3)

            rlexs = real(self%exs(ix), kind=dp)*geo%d(1)
            rleys = real(self%eys(iy), kind=dp)*geo%d(2)
            rlezs = real(self%ezs(iz), kind=dp)*geo%d(3)

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

            xt2 = self%k2x(ix)*((ytt/yt1)*(ztt/zt1))**2
            yt2 = self%k2y(iy)*((xtt/xt1)*(ztt/zt1))**2
            zt2 = self%k2z(iz)*((xtt/xt1)*(ytt/yt1))**2

            xyzk = xt2 + yt2 + zt2
            self%waves(i, j, k) = xyzk
          end do
        end do
      end do
    case ('100')
      do k = 1, self%nz_spec
        do j = 1, self%ny_spec
          do i = 1, self%nx_spec
            iy = i + self%sp_st(1)
            ix = j + self%sp_st(2)
            iz = k + self%sp_st(3)

            rlexs = real(self%exs(ix), kind=dp)*geo%d(1)
            rleys = real(self%eys(iy), kind=dp)*geo%d(2)
            rlezs = real(self%ezs(iz), kind=dp)*geo%d(3)

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

            xt2 = self%k2x(ix)*((ytt/yt1)*(ztt/zt1))**2
            yt2 = self%k2y(iy)*((xtt/xt1)*(ztt/zt1))**2
            zt2 = self%k2z(iz)*((xtt/xt1)*(ytt/yt1))**2

            xyzk = xt2 + yt2 + zt2
            self%waves(i, j, k) = xyzk
          end do
        end do
      end do
    case default
      do k = 1, self%nz_spec
        do j = 1, self%ny_spec
          do i = 1, self%nx_spec
            ix = i + self%sp_st(1)
            iy = j + self%sp_st(2)
            iz = k + self%sp_st(3)

            rlexs = real(self%exs(ix), kind=dp)*geo%d(1)
            rleys = real(self%eys(iy), kind=dp)*geo%d(2)
            rlezs = real(self%ezs(iz), kind=dp)*geo%d(3)

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

            xt2 = self%k2x(ix)*((ytt/yt1)*(ztt/zt1))**2
            yt2 = self%k2y(iy)*((xtt/xt1)*(ztt/zt1))**2
            zt2 = self%k2z(iz)*((xtt/xt1)*(ytt/yt1))**2

            xyzk = xt2 + yt2 + zt2
            self%waves(i, j, k) = xyzk
          end do
        end do
      end do
    end select

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

  real(dp) function get_km_re(self, i, j, k) result(re)
    implicit none

    class(poisson_fft_t) :: self
    integer, intent(in) :: i, j, k

    re = get_real(self%get_km(i, j, k))
  end function get_km_re

  real(dp) function get_km_im(self, i, j, k) result(re)
    implicit none

    class(poisson_fft_t) :: self
    integer, intent(in) :: i, j, k

    re = get_imag(self%get_km(i, j, k))
  end function get_km_im

  complex(dp) function get_km(self, i, j, k) result(km)
    implicit none

    class(poisson_fft_t) :: self
    integer, intent(in) :: i, j, k

    km = cmplx(self%trans_x_re(i)*get_real(self%ky(j))*self%trans_z_re(k), &
               self%trans_x_im(i)*get_imag(self%ky(j))*self%trans_z_im(k), &
               kind=dp)
  end function get_km

  real(dp) function get_real(complx) result(re)
    implicit none

    complex(dp), intent(in) :: complx

    re = real(complx, kind=dp)
  end function get_real

  real(dp) function get_imag(complx) result(im)
    implicit none

    complex(dp), intent(in) :: complx

    im = aimag(complx)
  end function get_imag

end module m_poisson_fft
