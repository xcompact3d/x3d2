module m_poisson_fft
  use m_common, only: dp, pi, CELL
  use m_field, only: field_t
  use m_mesh, only: mesh_t, geo_t, FFT_BC_000, FFT_BC_010, FFT_BC_100, &
                    FFT_BC_110, FFT_BC_111, FFT_BC_UNSUPPORTED
  use m_tdsops, only: dirps_t, tdsops_t
  use m_logger, only: log_warn, abort_run

  implicit none

  private :: fill_stretched_diagonals, interp_transfer

  !> The 1D factors of waves_set along one spectral dimension, see
  !! spectral_line. The _im factors differ from the _re ones only along a
  !! non-periodic r2c direction.
  type, private :: spectral_line_t
    real(dp), allocatable, dimension(:) :: k2_re, k2_im, t_re, t_im
  end type spectral_line_t

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
    !! (000, 010, 100, and 11x for both 110 and 111), and each case type
    !! supplies its own transforms, spectral postprocess and periodicity
    !! reorder.
    procedure(fft_forward), deferred :: fft_forward
    procedure(fft_backward), deferred :: fft_backward
    procedure(fft_postprocess), deferred :: fft_postprocess
    procedure :: enforce_periodicity
    procedure :: undo_periodicity
    procedure :: base_init
    procedure :: solve_poisson
    procedure :: stretching_matrix
    procedure :: waves_set
    procedure, private :: spectral_line
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
    case ('111')
      fft_bc = FFT_BC_111
    case default
      fft_bc = FFT_BC_UNSUPPORTED
    end select
  end function fft_bc_class

  pure function spectral_dirs(bc_code) result(dirs)
    !! The physical direction held in each dimension of the spectral array.
    !! The r2c transform runs along dim1, so a case with x non-periodic
    !! permutes a periodic direction into dim1 ahead of the FFT: 100 swaps
    !! x and y, 110 brings z to the front. 111 has no periodic direction
    !! and keeps the 110 layout, with z folded so that the r2c transform
    !! stands in for a cosine transform (see waves_set).
    implicit none

    character(len=3), intent(in) :: bc_code
    integer :: dirs(3)

    select case (bc_code)
    case ('100')
      dirs = [2, 1, 3]
    case ('110', '111')
      dirs = [3, 1, 2]
    case default
      dirs = [1, 2, 3]
    end select
  end function spectral_dirs

  subroutine base_init(self, mesh, xdirps, ydirps, zdirps, n_spec, n_sp_st, &
                       n_blk, blk_st)
    !! Sets up the state every case shares. The backend's init calls it
    !! once it has set bc_code, which the spectral extents it passes in
    !! already depend on. The supported BC codes are the cases of the
    !! backend factories, create_omp_poisson_fft and create_cuda_poisson_fft.
    !! n_blk and blk_st are passed on to stretching_matrix.
    implicit none
    class(poisson_fft_t) :: self
    type(mesh_t), intent(in), target :: mesh
    type(dirps_t), intent(in) :: xdirps, ydirps, zdirps
    integer, dimension(3), intent(in) :: n_spec ! Size of the spectral pencil
    integer, dimension(3), intent(in) :: n_sp_st ! Offset of the spectral pencil
    !> Size and offset of the block the stretching matrices cover
    integer, dimension(3), optional, intent(in) :: n_blk, blk_st

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
        call self%stretching_matrix(mesh%geo, xdirps, ydirps, zdirps, &
                                    n_blk, blk_st)
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

  subroutine stretching_matrix(self, geo, xdirps, ydirps, zdirps, n_blk, &
                               blk_st)
    !! Stretching necessitates a special operation in spectral space.
    !! The coefficients for the operation are stored in matrix form.
    !! The matrices and trans_* cover the spectral block of size n_blk at
    !! offset blk_st, the local spectral pencil by default. trans_* hold the
    !! modes of the block and take the local index.
    !!
    !! Ref. JCP 228 (2009), 5989–6015, Sec 5
    implicit none

    class(poisson_fft_t) :: self
    type(geo_t), intent(in) :: geo
    type(dirps_t), intent(in) :: xdirps, ydirps, zdirps
    integer, dimension(3), optional, intent(in) :: n_blk, blk_st

    real(dp) :: a0, a1
    integer :: i, k, ix, iz
    integer :: n(3), st(3)

    n = [self%nx_spec, self%ny_spec, self%nz_spec]
    st = self%sp_st
    if (present(n_blk)) n = n_blk
    if (present(blk_st)) st = blk_st

    allocate (self%trans_x_re(n(1)), self%trans_x_im(n(1)))
    allocate (self%trans_y_re(n(2)), self%trans_y_im(n(2)))
    allocate (self%trans_z_re(n(3)), self%trans_z_im(n(3)))

    self%trans_x_re = interp_transfer( &
                      xdirps%interpl_v2p, &
                      real(self%exs(st(1) + 1:st(1) + n(1)), kind=dp)*geo%d(1))
    self%trans_y_re = interp_transfer( &
                      ydirps%interpl_v2p, &
                      real(self%eys(st(2) + 1:st(2) + n(2)), kind=dp)*geo%d(2))
    self%trans_z_re = interp_transfer( &
                      zdirps%interpl_v2p, &
                      real(self%ezs(st(3) + 1:st(3) + n(3)), kind=dp)*geo%d(3))
    self%trans_x_im = self%trans_x_re
    self%trans_y_im = self%trans_y_re
    self%trans_z_im = self%trans_z_re

    if (trim(geo%stretching(2)) == 'bottom') then
      self%stretched_y_sym = .false.
      allocate (self%a_re(n(1), n(2), n(3), 5))
      allocate (self%a_im(n(1), n(2), n(3), 5))

      self%a_re = 0._dp
      self%a_im = 0._dp

      a0 = (geo%alpha(2)/pi + 1._dp/(2*pi*geo%beta(2)))*geo%L(2)
      a1 = -1._dp/(4*pi*geo%beta(2))*geo%L(2)

      call fill_stretched_diagonals(self, self%a_re, self%a_im, &
                                    STRETCH_MAT_BOTTOM, a0, a1, n, st)

      ! tweak the matrix to make it not singular
      self%a_re(1, 1, 1, 3) = 1._dp; self%a_im(1, 1, 1, 3) = 1._dp
      self%a_re(1, 1, 1, 4) = 0; self%a_im(1, 1, 1, 4) = 0
      self%a_re(1, 1, 1, 5) = 0; self%a_im(1, 1, 1, 5) = 0
    else
      self%stretched_y_sym = .true.
      allocate (self%a_odd_re(n(1), n(2)/2, n(3), 5))
      allocate (self%a_odd_im(n(1), n(2)/2, n(3), 5))
      allocate (self%a_even_re(n(1), n(2)/2, n(3), 5))
      allocate (self%a_even_im(n(1), n(2)/2, n(3), 5))

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
                                    STRETCH_MAT_ODD, a0, a1, n, st)
      call fill_stretched_diagonals(self, self%a_even_re, self%a_even_im, &
                                    STRETCH_MAT_EVEN, a0, a1, n, st)

      do k = 1, n(3)
        do i = 1, n(1)
          ix = i + st(1); iz = k + st(3)
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

  subroutine fill_stretched_diagonals(self, a_re, a_im, matrix_kind, a0, a1, &
                                      n, st)
    !! Fill the five diagonals of one stretching matrix on the spectral block
    !! of size n at offset st. Row j of the matrix couples the cosine modes
    !! iy_scale*iy + iy_off + (-2:2)*shift, with iy = j + st(2). trans_*
    !! take the local index of the block. The caller zeroes a_re and a_im
    !! beforehand.
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
    integer, intent(in) :: n(3), st(3)

    real(dp) :: c1, c2, c_first
    complex(dp) :: km_a1
    integer :: i, j, k, ix, iy, iz, r, nup1
    integer :: nrows, iy_scale, iy_off, shift
    logical :: is_even

    select case (matrix_kind)
    case (STRETCH_MAT_BOTTOM)
      nrows = n(2); iy_scale = 1; iy_off = 0; shift = 1
    case (STRETCH_MAT_ODD)
      nrows = n(2)/2; iy_scale = 2; iy_off = -1; shift = 2
    case (STRETCH_MAT_EVEN)
      nrows = n(2)/2; iy_scale = 2; iy_off = 0; shift = 2
    case default
      call abort_run('fill_stretched_diagonals: unknown matrix kind')
      nrows = 0; iy_scale = 0; iy_off = 0; shift = 0
    end select

    is_even = matrix_kind == STRETCH_MAT_EVEN
    ! First-row coefficient of the even diagonal + 1. Compute it outside the loop
    ! to keep the rounding of the old code (no fused multiply-add in the loop).
    c_first = a0*a1 - a1*a1
    nup1 = nrows - 1

    do k = 1, n(3)
      do j = 1, nrows
        do i = 1, n(1)
          ix = i + st(1); iy = j + st(2); iz = k + st(3)
          r = iy_scale*iy + iy_off

          ! diagonal
          c1 = a0*a0
          c2 = a1*a1
          if (iy == 1) then
            if (is_even) c1 = a0*a0 - a1*a1
            km_a1 = self%get_km(i, r + shift, k)
          else if (iy == nrows) then
            if (is_even) c1 = (a0 + a1)*(a0 + a1)
            km_a1 = self%get_km(i, r - shift, k)
          else
            km_a1 = self%get_km(i, r - shift, k) &
                    + self%get_km(i, r + shift, k)
          end if
          a_re(i, j, k, 3) = &
            -(get_real(self%kx(ix))*self%trans_y_re(r - st(2)) &
              *self%trans_z_re(k))**2 &
            - (get_real(self%kz(iz))*self%trans_y_re(r - st(2)) &
               *self%trans_x_re(i))**2 &
            - c1*self%get_km_re(i, r, k)**2 &
            - c2*self%get_km_re(i, r, k) &
            *get_real(km_a1)
          a_im(i, j, k, 3) = &
            -(get_imag(self%kx(ix))*self%trans_y_im(r - st(2)) &
              *self%trans_z_im(k))**2 &
            - (get_imag(self%kz(iz))*self%trans_y_im(r - st(2)) &
               *self%trans_x_im(i))**2 &
            - c1*self%get_km_im(i, r, k)**2 &
            - c2*self%get_km_im(i, r, k) &
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
              c1*(self%get_km_re(i, r, k) &
                  *self%get_km_re(i, r + shift, k)) &
              + c2*self%get_km_re(i, r + shift, k)**2
            a_im(i, j, k, 4) = &
              c1*(self%get_km_im(i, r, k) &
                  *self%get_km_im(i, r + shift, k)) &
              + c2*self%get_km_im(i, r + shift, k)**2
          end if

          ! diagonal + 2
          if (j <= nrows - 2) then
            c1 = a1*a1
            if (iy == 1 .and. .not. is_even) c1 = 2*a1*a1
            a_re(i, j, k, 5) = &
              -(c1*self%get_km_re(i, r + shift, k) &
                *self%get_km_re(i, r + 2*shift, k))
            a_im(i, j, k, 5) = &
              -(c1*self%get_km_im(i, r + shift, k) &
                *self%get_km_im(i, r + 2*shift, k))
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
              c1*(self%get_km_re(i, r, k) &
                  *self%get_km_re(i, r - shift, k)) &
              + c2*self%get_km_re(i, r - shift, k)**2
            a_im(i, j, k, 2) = &
              c1*(self%get_km_im(i, r, k) &
                  *self%get_km_im(i, r - shift, k)) &
              + c2*self%get_km_im(i, r - shift, k)**2
          end if

          ! diagonal - 2
          if (j >= 3) then
            a_re(i, j, k, 1) = &
              -(a1*a1*self%get_km_re(i, r - shift, k) &
                *self%get_km_re(i, r - 2*shift, k))
            a_im(i, j, k, 1) = &
              -(a1*a1*self%get_km_im(i, r - shift, k) &
                *self%get_km_im(i, r - 2*shift, k))
          end if
        end do
      end do
    end do
  end subroutine fill_stretched_diagonals

  subroutine waves_set(self, geo, xdirps, ydirps, zdirps)
    !! Spectral equivalence constants: the modified wavenumber of the
    !! discrete div(grad) on each spectral entry,
    !!   k2x*(ty*tz)**2 + k2y*(tx*tz)**2 + k2z*(tx*ty)**2,
    !! where t is the transfer function of the staggered interpolation. The
    !! expression is symmetric in the three directions, so it is built from
    !! the 1D factors of each spectral dimension (spectral_line), whatever
    !! direction spectral_dirs puts there.
    !!
    !! The real and imaginary parts differ only when the r2c direction (dim1)
    !! is non-periodic, as in the 111 case. The r2c transform of the folded
    !! line, with the half-shift twiddle of the postprocess, leaves the cosine
    !! mode k in the real part of entry k and the mode N - k in its imaginary
    !! part, so the imaginary part takes the factors of the mirrored mode.
    !! Incompact3d's waves does the same for bcz == 1.
    !!
    !! Ref. JCP 228 (2009), 5989–6015, Sec 4
    implicit none

    class(poisson_fft_t) :: self
    type(geo_t), intent(in) :: geo
    type(dirps_t), intent(in) :: xdirps, ydirps, zdirps

    type(spectral_line_t) :: l1, l2, l3
    integer :: i, j, k, dirs(3)
    real(dp) :: w_re, w_im

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

    dirs = spectral_dirs(self%bc_code)
    call self%spectral_line(dirs(1), geo, xdirps, ydirps, zdirps, &
                            self%nx_spec, self%sp_st(1), .true., l1)
    call self%spectral_line(dirs(2), geo, xdirps, ydirps, zdirps, &
                            self%ny_spec, self%sp_st(2), .false., l2)
    call self%spectral_line(dirs(3), geo, xdirps, ydirps, zdirps, &
                            self%nz_spec, self%sp_st(3), .false., l3)

    do k = 1, self%nz_spec
      do j = 1, self%ny_spec
        do i = 1, self%nx_spec
          w_re = l1%k2_re(i)*(l2%t_re(j)*l3%t_re(k))**2 &
                 + l2%k2_re(j)*(l1%t_re(i)*l3%t_re(k))**2 &
                 + l3%k2_re(k)*(l1%t_re(i)*l2%t_re(j))**2
          w_im = l1%k2_im(i)*(l2%t_im(j)*l3%t_im(k))**2 &
                 + l2%k2_im(j)*(l1%t_im(i)*l3%t_im(k))**2 &
                 + l3%k2_im(k)*(l1%t_im(i)*l2%t_im(j))**2
          self%waves(i, j, k) = cmplx(w_re, w_im, kind=dp)
        end do
      end do
    end do

  end subroutine waves_set

  subroutine spectral_line(self, dir, geo, xdirps, ydirps, zdirps, n, st, &
                           is_r2c, line)
    !! The 1D factors of waves_set along the spectral dimension that holds
    !! direction dir: the squared modified wavenumber k2 and the
    !! interpolation transfer function t of the n local modes from offset st.
    !! Along the r2c dimension (is_r2c) of a non-periodic direction, the
    !! imaginary factors are those of the mirrored mode N - k (see
    !! waves_set). The mean mode has no mirror; its imaginary part carries
    !! no data, so it keeps its own factors.
    implicit none

    class(poisson_fft_t), intent(in) :: self
    integer, intent(in) :: dir, n, st
    type(geo_t), intent(in) :: geo
    type(dirps_t), intent(in) :: xdirps, ydirps, zdirps
    logical, intent(in) :: is_r2c
    type(spectral_line_t), intent(out) :: line

    complex(dp), allocatable, dimension(:) :: k2, e
    real(dp), allocatable, dimension(:) :: t
    integer :: i, m, n_glob
    logical :: periodic

    select case (dir)
    case (1)
      k2 = self%k2x; e = self%exs; periodic = self%periodic_x
      t = interp_transfer(xdirps%interpl_v2p, real(e, kind=dp)*geo%d(1))
    case (2)
      k2 = self%k2y; e = self%eys; periodic = self%periodic_y
      t = interp_transfer(ydirps%interpl_v2p, real(e, kind=dp)*geo%d(2))
    case default
      k2 = self%k2z; e = self%ezs; periodic = self%periodic_z
      t = interp_transfer(zdirps%interpl_v2p, real(e, kind=dp)*geo%d(3))
    end select
    n_glob = size(k2)

    allocate (line%k2_re(n), line%k2_im(n), line%t_re(n), line%t_im(n))
    do i = 1, n
      line%k2_re(i) = real(k2(st + i), kind=dp)
      line%t_re(i) = t(st + i)
      m = st + i
      if (is_r2c .and. (.not. periodic) .and. m > 1) m = n_glob - m + 2
      line%k2_im(i) = real(k2(m), kind=dp)
      line%t_im(i) = t(m)
    end do

  end subroutine spectral_line

  pure elemental real(dp) function interp_transfer(interp, w) result(t)
    !! Transfer function of the staggered (v2p) interpolation at the scaled
    !! wavenumber w = k*d, the factor waves_set and stretching_matrix apply
    !! to the wavenumbers of the other two directions.
    implicit none

    type(tdsops_t), intent(in) :: interp
    real(dp), intent(in) :: w

    t = 2*(interp%a*cos(w*0.5_dp) + interp%b*cos(w*1.5_dp) &
           + interp%c*cos(w*2.5_dp) + interp%d*cos(w*3.5_dp)) &
        /(1._dp + 2*interp%alpha*cos(w))

  end function interp_transfer

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
