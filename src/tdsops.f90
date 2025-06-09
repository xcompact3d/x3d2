module m_tdsops
  use iso_fortran_env, only: stderr => error_unit

  use m_common, only: dp, pi, VERT, CELL, &
                      BC_PERIODIC, BC_NEUMANN, BC_DIRICHLET

  implicit none

  type :: tdsops_t
      !! Tridiagonal Solver Operators class.
      !!
      !! Operator arrays are preprocessed in this class based on the arguments
      !! provided. dist_fw and dist_bw are used in the first phase of the
      !! distributed tridiagonal solver algorithm. dist_sa and dist_sc are used
      !! in the final substitution phase. See the kernels_dist.f90 files in the
      !! relevant backend folders.
      !! coeff arrays define the specific rules of building the RHS
      !! corresponding to the tridiagonal system to be solved, and used only in
      !! the first phase of the distributed algorithm when building the RHS.
      !! If a boundary condition is defined then coeffs_s and coeffs_e differ
      !! from coeffs array and define the RHS rule for the first and last 4
      !! entries in the tridiagonal system (n_halo = 4).
      !!
      !! This class does not know about the current rank or its relative
      !! location among other ranks. All the operator arrays here are used when
      !! executing a distributed tridiagonal solver phase one or two.
    real(dp), allocatable, dimension(:) :: dist_fw, dist_bw, & !! fw/bw phase
                                           dist_sa, dist_sc, & !! back subs.
                                           dist_af !! the auxiliary factors
    real(dp), allocatable, dimension(:) :: thom_f, thom_s, thom_w, thom_p
    real(dp), allocatable :: stretch(:), stretch_correct(:)
    real(dp), allocatable :: coeffs(:), coeffs_s(:, :), coeffs_e(:, :)
    real(dp) :: alpha, a, b, c = 0._dp, d = 0._dp !! Compact scheme coeffs
    logical :: periodic
    integer :: n_tds !! Tridiagonal system size
    integer :: n_rhs !! Right-hand-side builder size
    integer :: move = 0 !! move between vertices and cell centres
    integer :: n_halo !! number of halo points
  contains
    procedure :: deriv_1st, deriv_2nd, interpl_mid, stagder_1st
    procedure :: preprocess_dist, preprocess_thom
  end type tdsops_t

  interface tdsops_t
    module procedure tdsops_init
  end interface tdsops_t

  type :: dirps_t
    !! Directional tridiagonal solver container.
    !!
    !! This class contains the preprocessed tridiagonal solvers for operating
    !! in each coordinate direction.
    class(tdsops_t), allocatable :: der1st, der1st_sym, der2nd, der2nd_sym, &
      stagder_v2p, stagder_p2v, interpl_v2p, interpl_p2v
    integer :: dir
  end type dirps_t

contains

  function tdsops_init( &
    n_tds, delta, operation, scheme, bc_start, bc_end, &
    stretch, stretch_correct, n_halo, from_to, sym, c_nu, nu0_nu &
    ) result(tdsops)
    !! Constructor function for the tdsops_t class.
    !!
    !! 'n_tds', 'delta', 'operation', 'scheme', 'bc_start', and 'bc_end' are
    !! necessary arguments. The remaining arguments are optional.
    !!
    !! 'stretch' is for obtaining the correct derivations in a stretched mesh
    !! 'stretch_correct' is for correcting the second derivative with the first
    !!
    !! 'from_to' is necessary for interpolation and staggared derivative, and
    !! it can be 'v2p' or 'p2v'.
    !! If the specific region the instance is operating is not a boundary
    !! region, then 'bc_start' and 'bc_end' are BC_HALO.
    !!
    !! 'sym' is relevant when the BC is free-slip. If sym is .true. then it
    !! means the field we operate on is assumed to be an even function
    !! (symmetric, cos type) accross the boundary. If it is .false. it means
    !! the field is assumed to be an odd function (anti-symmetric, sin type).
    !!
    !! 'c_nu', 'nu0_nu' are relevant when operation is second order
    !! derivative and scheme is compact6-hyperviscous.
    implicit none

    type(tdsops_t) :: tdsops !! return value of the function

    integer, intent(in) :: n_tds !! Tridiagonal system size
    real(dp), intent(in) :: delta !! Grid spacing
    character(*), intent(in) :: operation, scheme
    integer, intent(in) :: bc_start, bc_end !! Boundary Cond.
    real(dp), optional, intent(in) :: stretch(:) !! Stretching coefficients
    real(dp), optional, intent(in) :: stretch_correct(:) !! Stretch correction
    integer, optional, intent(in) :: n_halo !! Number of halo cells
    character(*), optional, intent(in) :: from_to !! 'v2p' or 'p2v'
    logical, optional, intent(in) :: sym !! (==npaire), only for Neumann BCs
    real(dp), optional, intent(in) :: c_nu, nu0_nu !! params for hypervisc.

    integer :: n, n_stencil

    tdsops%n_tds = n_tds

    ! we need special treatment in the right-hand-side build stage for
    ! the very last point in the domain if output length is smaller than
    ! the input length
    if (present(from_to)) then
      if ((bc_end == BC_NEUMANN .or. bc_end == BC_DIRICHLET) &
          .and. from_to == 'v2p') then
        tdsops%n_rhs = n_tds + 1
      else
        tdsops%n_rhs = n_tds
      end if
    else
      tdsops%n_rhs = n_tds
    end if

    if (present(n_halo)) then
      tdsops%n_halo = n_halo
      if (n_halo /= 4) then
        write (stderr, '("Warning: n_halo is set to ", i2, "be careful! &
                          &The default is 4 and there are quite a few &
                          &places where things are hardcoded assuming &
                          &n_halo is 4.")') n_halo
      end if
    else
      tdsops%n_halo = 4
    end if

    ! n_rhs >= n_tds, n is used when its better to allocate a larger size
    n = tdsops%n_rhs

    ! preprocessed coefficient arrays for the distributed algorithm
    allocate (tdsops%dist_fw(n), tdsops%dist_bw(n))
    allocate (tdsops%dist_sa(n), tdsops%dist_sc(n))
    allocate (tdsops%dist_af(n))

    ! preprocessed coefficient arrays for the Thomas algorithm
    allocate (tdsops%thom_f(n), tdsops%thom_s(n))
    allocate (tdsops%thom_w(n), tdsops%thom_p(n))

    ! RHS coefficient arrays
    n_stencil = 2*tdsops%n_halo + 1
    allocate (tdsops%coeffs(n_stencil))
    allocate (tdsops%coeffs_s(n_stencil, tdsops%n_halo))
    allocate (tdsops%coeffs_e(n_stencil, tdsops%n_halo))

    allocate (tdsops%stretch(n_tds))
    if (present(stretch)) then
      tdsops%stretch(:) = stretch(:)
    else
      tdsops%stretch(:) = 1._dp
    end if

    allocate (tdsops%stretch_correct(n_tds))
    if (present(stretch_correct)) then
      tdsops%stretch_correct(:) = stretch_correct(:)
    else
      tdsops%stretch_correct(:) = 0._dp
    end if

    tdsops%periodic = bc_start == BC_PERIODIC .and. bc_end == BC_PERIODIC

    if (operation == 'first-deriv') then
      call tdsops%deriv_1st(delta, scheme, bc_start, bc_end, sym)
    else if (operation == 'second-deriv') then
      call tdsops%deriv_2nd(delta, scheme, bc_start, bc_end, sym, &
                            c_nu, nu0_nu)
    else if (operation == 'interpolate') then
      call tdsops%interpl_mid(scheme, from_to, bc_start, bc_end, sym)
    else if (operation == 'stag-deriv') then
      call tdsops%stagder_1st(delta, scheme, from_to, bc_start, bc_end, sym)
    else
      error stop 'operation is not defined'
    end if

    select case (from_to)
    case ('v2p')
      tdsops%move = 1
    case ('p2v')
      tdsops%move = -1
    case default
      tdsops%move = 0
    end select

    if (tdsops%dist_sa(n_tds) > 1d-16) then
      print *, 'There are ', n_tds, 'points in a subdomain, it may be too few!'
      print *, 'The entry distributed solver disregards in "' &
        //operation//'" operation is:', tdsops%dist_sa(n_tds)
      print *, 'It may result in numerical errors with the distributed solver!'
    end if

  end function tdsops_init

  subroutine deriv_1st(self, delta, scheme, bc_start, bc_end, sym)
    implicit none

    class(tdsops_t), intent(inout) :: self
    real(dp), intent(in) :: delta
    character(*), intent(in) :: scheme
    integer, intent(in) :: bc_start, bc_end
    logical, optional, intent(in) :: sym

    real(dp), allocatable :: dist_b(:)
    real(dp) :: alpha, afi, bfi
    integer :: i, n, n_halo
    logical :: symmetry

    if (self%n_halo < 2) error stop 'First derivative require n_halo >= 2'

    if (present(sym)) then
      symmetry = sym
    else
      symmetry = .false.
    end if

    ! alpha is alfa

    select case (scheme)
    case ('compact6')
      alpha = 1._dp/3._dp
      afi = 7._dp/9._dp/delta
      bfi = 1._dp/36._dp/delta
    case default
      error stop 'scheme is not defined'
    end select

    self%alpha = alpha
    self%a = afi; self%b = bfi

    self%coeffs(:) = [0._dp, 0._dp, -bfi, -afi, &
                      0._dp, &
                      afi, bfi, 0._dp, 0._dp]

    do i = 1, self%n_halo
      self%coeffs_s(:, i) = self%coeffs(:)
      self%coeffs_e(:, i) = self%coeffs(:)
    end do

    self%dist_sa(:) = alpha; self%dist_sc(:) = alpha

    n = self%n_tds
    n_halo = self%n_halo

    allocate (dist_b(self%n_rhs))
    dist_b(:) = 1._dp

    select case (bc_start)
    case (BC_NEUMANN)
      if (symmetry) then
        ! sym == .true.; d(uu)/dx, dv/dx, dw/dx
        !                d(vv)/dy, du/dy, dw/dy
        !                d(ww)/dz, du/dz, dv/dz
        self%dist_sa(1) = 0._dp
        self%dist_sc(1) = 0._dp
        self%coeffs_s(:, 1) = [0._dp, 0._dp, 0._dp, 0._dp, &
                               0._dp, &
                               0._dp, 0._dp, 0._dp, 0._dp]
        self%coeffs_s(:, 2) = [0._dp, 0._dp, 0._dp, -afi, &
                               -bfi, &
                               afi, bfi, 0._dp, 0._dp]
      else
        ! sym == .false.; d(uv)/dx, d(uw)/dx, du/dx
        !                 d(vu)/dy, d(vw)/dy, dv/dy
        !                 d(wu)/dz, d(wv)/dz, dw/dz
        self%dist_sa(1) = 0._dp
        self%dist_sc(1) = 2*alpha
        self%coeffs_s(:, 1) = [0._dp, 0._dp, 0._dp, 0._dp, &
                               0._dp, &
                               2*afi, 2*bfi, 0._dp, 0._dp]
        self%coeffs_s(:, 2) = [0._dp, 0._dp, 0._dp, -afi, &
                               bfi, &
                               afi, bfi, 0._dp, 0._dp]
      end if
    case (BC_DIRICHLET)
      ! first line
      self%dist_sa(1) = 0._dp
      self%dist_sc(1) = 2._dp
      self%coeffs_s(:, 1) = [0._dp, 0._dp, 0._dp, 0._dp, &
                             -2.5_dp, &
                             2._dp, 0.5_dp, 0._dp, 0._dp]
      self%coeffs_s(:, 1) = self%coeffs_s(:, 1)/delta
      ! second line
      self%dist_sa(2) = 0.25_dp
      self%dist_sc(2) = 0.25_dp
      self%coeffs_s(:, 2) = [0._dp, 0._dp, 0._dp, -0.75_dp, &
                             0._dp, &
                             0.75_dp, 0._dp, 0._dp, 0._dp]
      self%coeffs_s(:, 2) = self%coeffs_s(:, 2)/delta
    end select

    select case (bc_end)
    case (BC_NEUMANN)
      if (symmetry) then
        ! sym == .true.; d(uu)/dx, dv/dx, dw/dx
        !                d(vv)/dy, du/dy, dw/dy
        !                d(ww)/dz, du/dz, dv/dz
        self%dist_sa(n) = 0._dp
        self%dist_sc(n) = 0._dp
        self%coeffs_e(:, n_halo) = [0._dp, 0._dp, 0._dp, 0._dp, &
                                    0._dp, &
                                    0._dp, 0._dp, 0._dp, 0._dp]
        self%coeffs_e(:, n_halo - 1) = [0._dp, 0._dp, -bfi, -afi, &
                                        bfi, &
                                        afi, 0._dp, 0._dp, 0._dp]
      else
        ! sym == .false.; d(uv)/dx, d(uw)/dx, du/dx
        !                 d(vu)/dy, d(vw)/dy, dv/dy
        !                 d(wu)/dz, d(wv)/dz, dw/dz
        self%dist_sa(n) = 2*alpha
        self%dist_sc(n) = 0._dp
        self%coeffs_e(:, n_halo) = [0._dp, 0._dp, -2*bfi, -2*afi, &
                                    0._dp, &
                                    0._dp, 0._dp, 0._dp, 0._dp]
        self%coeffs_e(:, n_halo - 1) = [0._dp, 0._dp, -bfi, -afi, &
                                        -bfi, &
                                        afi, 0._dp, 0._dp, 0._dp]
      end if
    case (BC_DIRICHLET)
      ! last line
      self%dist_sa(n) = 2._dp
      self%dist_sc(n) = 0._dp
      self%coeffs_e(:, n_halo) = [0._dp, 0._dp, -0.5_dp, -2._dp, &
                                  2.5_dp, &
                                  0._dp, 0._dp, 0._dp, 0._dp]
      self%coeffs_e(:, n_halo) = self%coeffs_e(:, n_halo)/delta
      ! second last line
      self%dist_sa(n - 1) = 0.25_dp
      self%dist_sc(n - 1) = 0.25_dp
      self%coeffs_e(:, n_halo - 1) = [0._dp, 0._dp, 0._dp, -0.75_dp, &
                                      0._dp, &
                                      0.75_dp, 0._dp, 0._dp, 0._dp]
      self%coeffs_e(:, n_halo - 1) = self%coeffs_e(:, n_halo - 1)/delta
    end select

    call self%preprocess_thom(dist_b)
    call self%preprocess_dist(dist_b)

  end subroutine deriv_1st

  subroutine deriv_2nd(self, delta, scheme, bc_start, bc_end, sym, &
                       c_nu, nu0_nu)
    implicit none

    class(tdsops_t), intent(inout) :: self
    real(dp), intent(in) :: delta
    character(*), intent(in) :: scheme
    integer, intent(in) :: bc_start, bc_end
    logical, optional, intent(in) :: sym
    real(dp), optional, intent(in) :: c_nu, nu0_nu

    real(dp), allocatable :: dist_b(:)
    real(dp) :: alpha, asi, bsi, csi, dsi
    real(dp) :: dpis3, xnpi2, xmpi2, den, d2, temp1, temp2
    integer :: i, n, n_halo
    logical :: symmetry

    if (self%n_halo < 4) error stop 'Second derivative require n_halo >= 4'

    if (present(sym)) then
      symmetry = sym
    else
      symmetry = .false.
    end if

    d2 = delta*delta

    ! alpha is alsa

    select case (scheme)
    case ('compact6')
      alpha = 2._dp/11._dp
      asi = 12._dp/11._dp/d2
      bsi = 3._dp/44._dp/d2
      csi = 0._dp
      dsi = 0._dp
    case ('compact6-hyperviscous')
      if (present(c_nu) .and. present(nu0_nu)) then
        dpis3 = 2._dp*pi/3._dp
        xnpi2 = pi*pi*(1._dp + nu0_nu)
        xmpi2 = dpis3*dpis3*(1._dp + c_nu*nu0_nu)
        den = 405._dp*xnpi2 - 640._dp*xmpi2 + 144._dp
        alpha = 0.5_dp - (320._dp*xmpi2 - 1296._dp)/den
        asi = -(4329._dp*xnpi2/8._dp - 32._dp*xmpi2 &
                - 140._dp*xnpi2*xmpi2 + 286._dp)/den/d2
        bsi = (2115._dp*xnpi2 - 1792._dp*xmpi2 &
               - 280._dp*xnpi2*xmpi2 + 1328._dp)/den/(4._dp*d2)
        csi = -(7695._dp*xnpi2/8._dp + 288._dp*xmpi2 &
                - 180._dp*xnpi2*xmpi2 - 2574._dp)/den/(9._dp*d2)
        dsi = (198._dp*xnpi2 + 128._dp*xmpi2 &
               - 40._dp*xnpi2*xmpi2 - 736._dp)/den/(16._dp*d2)
      else
        error stop 'compact6-hyperviscous requires c_nu and nu0_nu'
      end if
    case default
      error stop 'scheme is not defined'
    end select

    self%alpha = alpha
    self%a = asi; self%b = bsi; self%c = csi; self%d = dsi

    self%coeffs(:) = [dsi, csi, bsi, asi, &
                      -2._dp*(asi + bsi + csi + dsi), &
                      asi, bsi, csi, dsi]

    do i = 1, self%n_halo
      self%coeffs_s(:, i) = self%coeffs(:)
      self%coeffs_e(:, i) = self%coeffs(:)
    end do

    self%dist_sa(:) = alpha; self%dist_sc(:) = alpha

    n = self%n_tds
    n_halo = self%n_halo

    allocate (dist_b(self%n_rhs))
    dist_b(:) = 1._dp

    select case (bc_start)
    case (BC_NEUMANN)
      if (symmetry) then
        ! sym == .true.; d2v/dx2, d2w/dx2
        !                d2u/dy2, d2w/dy2
        !                d2u/dz2, d2v/dz2
        self%dist_sa(1) = 0._dp
        self%dist_sc(1) = 2*alpha
        self%coeffs_s(:, 1) = [0._dp, 0._dp, 0._dp, 0._dp, &
                               -2*asi - 2*bsi - 2*csi - 2*dsi, &
                               2*asi, 2*bsi, 2*csi, 2*dsi]
        self%coeffs_s(:, 2) = [0._dp, 0._dp, 0._dp, asi, &
                               -2*asi - bsi - 2*csi - 2*dsi, &
                               asi + csi, bsi + dsi, csi, dsi]
        self%coeffs_s(:, 3) = [0._dp, 0._dp, bsi, asi + csi, &
                               -2*asi - 2*bsi - 2*csi - dsi, &
                               asi, bsi, csi, dsi]
        self%coeffs_s(:, 4) = [0._dp, csi, bsi + dsi, asi, &
                               -2*asi - 2*bsi - 2*csi - 2*dsi, &
                               asi, bsi, csi, dsi]
      else
        ! sym == .false.; d2u/dx2
        !                 d2v/dy2
        !                 d2w/dz2
        self%dist_sa(1) = 0._dp
        self%dist_sc(1) = 0._dp
        self%coeffs_s(:, 1) = [0._dp, 0._dp, 0._dp, 0._dp, &
                               0._dp, &
                               0._dp, 0._dp, 0._dp, 0._dp]
        self%coeffs_s(:, 2) = [0._dp, 0._dp, 0._dp, asi, &
                               -2*asi - 3*bsi - 2*csi - 2*dsi, &
                               asi - csi, bsi - dsi, csi, dsi]
        self%coeffs_s(:, 3) = [0._dp, 0._dp, bsi, asi - csi, &
                               -2*asi - 2*bsi - 2*csi - 3*dsi, &
                               asi, bsi, csi, dsi]
        self%coeffs_s(:, 4) = [0._dp, -csi, bsi - dsi, asi, &
                               -2*asi - 2*bsi - 2*csi - 2*dsi, &
                               asi, bsi, csi, dsi]
      end if
    case (BC_DIRICHLET)
      ! first line
      self%dist_sa(1) = 0._dp
      self%dist_sc(1) = 11._dp
      self%coeffs_s(:, 1) = [0._dp, 0._dp, 0._dp, 0._dp, &
                             13._dp/d2, &
                             -27._dp/d2, 15._dp/d2, -1._dp/d2, 0._dp]
      ! second line
      self%dist_sa(2) = 0.1_dp
      self%dist_sc(2) = 0.1_dp
      self%coeffs_s(:, 2) = [0._dp, 0._dp, 0._dp, 1.2_dp/d2, &
                             -2.4_dp/d2, &
                             1.2_dp/d2, 0._dp, 0._dp, 0._dp]
      ! third line
      self%dist_sa(3) = 2._dp/11._dp
      self%dist_sc(3) = 2._dp/11._dp
      temp1 = 3._dp/44._dp/d2; temp2 = 12._dp/11._dp/d2
      self%coeffs_s(:, 3) = [0._dp, 0._dp, temp1, temp2, &
                             -2._dp*(temp1 + temp2), &
                             temp2, temp1, 0._dp, 0._dp]
      ! fourth line is same as third
      self%dist_sa(4) = 2._dp/11._dp
      self%dist_sc(4) = 2._dp/11._dp
      self%coeffs_s(:, 4) = self%coeffs_s(:, 3)
    end select

    select case (bc_end)
    case (BC_NEUMANN)
      if (symmetry) then
        ! sym == .true.; d2v/dx2, d2w/dx2
        !                d2u/dy2, d2w/dy2
        !                d2u/dz2, d2v/dz2
        self%dist_sa(n) = 2*alpha
        self%dist_sc(n) = 0._dp
        self%coeffs_e(:, 4) = [dsi, csi, bsi, asi, &
                               -2*asi - 2*bsi - 2*csi - 2*dsi, &
                               0._dp, 0._dp, 0._dp, 0._dp]
        self%coeffs_e(:, 3) = [dsi, csi, bsi + dsi, asi + csi, &
                               -2*asi - bsi - 2*csi - 2*dsi, &
                               asi, 0._dp, 0._dp, 0._dp]
        self%coeffs_e(:, 2) = [dsi, csi, bsi, asi, &
                               -2*asi - 2*bsi - 2*csi - dsi, &
                               asi + csi, bsi, 0._dp, 0._dp]
        self%coeffs_e(:, 1) = [dsi, csi, bsi, asi, &
                               -2*asi - 2*bsi - 2*csi - 2*dsi, &
                               asi, bsi + dsi, csi, 0._dp]
      else
        ! sym == .false.; d2u/dx2
        !                 d2v/dy2
        !                 d2w/dz2
        self%dist_sa(n) = 0._dp
        self%dist_sc(n) = 0._dp
        self%coeffs_e(:, 4) = [0._dp, 0._dp, 0._dp, 0._dp, &
                               0._dp, &
                               0._dp, 0._dp, 0._dp, 0._dp]
        self%coeffs_e(:, 3) = [dsi, csi, bsi - dsi, asi - csi, &
                               -2*asi - 3*bsi - 2*csi - 2*dsi, &
                               asi, 0._dp, 0._dp, 0._dp]
        self%coeffs_e(:, 2) = [dsi, csi, bsi, asi, &
                               -2*asi - 2*bsi - 2*csi - 3*dsi, &
                               asi - csi, bsi, 0._dp, 0._dp]
        self%coeffs_e(:, 1) = [dsi, csi, bsi, asi, &
                               -2*asi - 2*bsi - 2*csi - 2*dsi, &
                               asi, bsi - dsi, -csi, 0._dp]
      end if
    case (BC_DIRICHLET)
      ! last line
      self%dist_sa(n) = 11._dp
      self%dist_sc(n) = 0._dp
      self%coeffs_e(:, 4) = [0._dp, -1._dp/d2, 15._dp/d2, -27._dp/d2, &
                             13._dp/d2, &
                             0._dp, 0._dp, 0._dp, 0._dp]
      ! second last line
      self%dist_sa(n - 1) = 0.1_dp
      self%dist_sc(n - 1) = 0.1_dp
      self%coeffs_e(:, 3) = [0._dp, 0._dp, 0._dp, 1.2_dp/d2, &
                             -2.4_dp/d2, &
                             1.2_dp/d2, 0._dp, 0._dp, 0._dp]
      ! third last line
      self%dist_sa(n - 2) = 2._dp/11._dp
      self%dist_sc(n - 2) = 2._dp/11._dp
      temp1 = 3._dp/44._dp/d2; temp2 = 12._dp/11._dp/d2
      self%coeffs_e(:, 2) = [0._dp, 0._dp, temp1, temp2, &
                             -2._dp*(temp1 + temp2), &
                             temp2, temp1, 0._dp, 0._dp]
      ! fourth last line is same as third last
      self%dist_sa(n - 3) = 2._dp/11._dp
      self%dist_sc(n - 3) = 2._dp/11._dp
      self%coeffs_e(:, 1) = self%coeffs_e(:, 2)
    end select

    call self%preprocess_thom(dist_b)
    call self%preprocess_dist(dist_b)

  end subroutine deriv_2nd

  subroutine interpl_mid(self, scheme, from_to, bc_start, bc_end, sym)
    implicit none

    class(tdsops_t), intent(inout) :: self
    character(*), intent(in) :: scheme, from_to
    integer, intent(in) :: bc_start, bc_end
    logical, optional, intent(in) :: sym

    real(dp), allocatable :: dist_b(:)
    real(dp) :: alpha, aici, bici, cici, dici
    integer :: i, n, n_halo

    if (self%n_halo < 4) error stop 'Interpolation require n_halo >= 4'

    ! alpha is ailcai

    select case (scheme)
    case ('classic')
      alpha = 0.3_dp
      aici = 0.75_dp
      bici = 0.05_dp
      cici = 0._dp
      dici = 0._dp
    case ('optimised')
      alpha = 0.461658_dp
      dici = 0.00146508_dp
      aici = (75._dp + 70._dp*alpha - 640._dp*dici)/128._dp
      bici = (-25._dp + 126._dp*alpha + 2304._dp*dici)/256._dp
      cici = (3._dp - 10._dp*alpha - 1280._dp*dici)/256._dp
    case ('aggressive')
      alpha = 0.49_dp
      aici = (75._dp + 70._dp*alpha)/128._dp
      bici = (-25._dp + 126._dp*alpha)/256._dp
      cici = (3._dp - 10._dp*alpha)/256._dp
      dici = 0._dp
    case default
      error stop 'scheme is not defined'
    end select

    self%alpha = alpha
    self%a = aici; self%b = bici; self%c = cici; self%d = dici

    select case (from_to)
    case ('v2p')
      self%coeffs(:) = [0._dp, dici, cici, bici, &
                        aici, &
                        aici, bici, cici, dici]
    case ('p2v')
      self%coeffs(:) = [dici, cici, bici, aici, &
                        aici, &
                        bici, cici, dici, 0._dp]
    end select

    do i = 1, self%n_halo
      self%coeffs_s(:, i) = self%coeffs(:)
      self%coeffs_e(:, i) = self%coeffs(:)
    end do

    self%dist_sa(:) = alpha; self%dist_sc(:) = alpha

    n = self%n_tds
    n_halo = self%n_halo

    allocate (dist_b(self%n_rhs))
    dist_b(:) = 1._dp

    if ((bc_start == BC_DIRICHLET) .or. (bc_start == BC_NEUMANN)) then
      self%dist_sa(1) = 0._dp

      select case (from_to)
      case ('v2p')
        ! sym is always .true.
        dist_b(1) = 1._dp + alpha
        self%coeffs_s(:, 1) = [0._dp, 0._dp, 0._dp, 0._dp, &
                               aici, &
                               aici + bici, bici + cici, cici + dici, dici]
        self%coeffs_s(:, 2) = [0._dp, 0._dp, 0._dp, bici, &
                               aici + cici, &
                               aici + dici, bici, cici, dici]
        self%coeffs_s(:, 3) = [0._dp, 0._dp, cici, bici + dici, &
                               aici, &
                               aici, bici, cici, dici]
      case ('p2v')
        ! sym is always .true.
        self%dist_sc(1) = 2*alpha
        self%coeffs_s(:, 1) = [0._dp, 0._dp, 0._dp, 0._dp, &
                               2*aici, &
                               2*bici, 2*cici, 2*dici, 0._dp]
        self%coeffs_s(:, 2) = [0._dp, 0._dp, 0._dp, aici + bici, &
                               aici + cici, &
                               bici + dici, cici, dici, 0._dp]
        self%coeffs_s(:, 3) = [0._dp, 0._dp, bici + cici, aici + dici, &
                               aici, &
                               bici, cici, dici, 0._dp]
        self%coeffs_s(:, 4) = [0._dp, cici + dici, bici, aici, &
                               aici, &
                               bici, cici, dici, 0._dp]
      end select
    end if

    if ((bc_end == BC_DIRICHLET) .or. (bc_end == BC_NEUMANN)) then
      self%dist_sc(n) = 0._dp

      select case (from_to)
      case ('v2p')
        ! sym is always .true.
        dist_b(n) = 1._dp + alpha
        self%coeffs_e(:, 4) = 0._dp
        self%coeffs_e(:, 3) = [0._dp, dici, cici + dici, bici + cici, &
                               aici + bici, &
                               aici, 0._dp, 0._dp, 0._dp]
        self%coeffs_e(:, 2) = [0._dp, dici, cici, bici, &
                               aici + dici, &
                               aici + cici, bici, 0._dp, 0._dp]
        self%coeffs_e(:, 1) = [0._dp, dici, cici, bici, &
                               aici, &
                               aici, bici + dici, cici, 0._dp]
      case ('p2v')
        ! sym is always .true.
        self%dist_sa(n) = 2*alpha
        self%coeffs_e(:, 4) = [2*dici, 2*cici, 2*bici, 2*aici, &
                               0._dp, &
                               0._dp, 0._dp, 0._dp, 0._dp]
        self%coeffs_e(:, 3) = [dici, cici, bici + dici, aici + cici, &
                               aici + bici, &
                               0._dp, 0._dp, 0._dp, 0._dp]
        self%coeffs_e(:, 2) = [dici, cici, bici, aici, &
                               aici + dici, &
                               bici + cici, 0._dp, 0._dp, 0._dp]
        self%coeffs_e(:, 1) = [dici, cici, bici, aici, &
                               aici, &
                               bici, cici + dici, 0._dp, 0._dp]
      end select
    end if

    call self%preprocess_thom(dist_b)
    call self%preprocess_dist(dist_b)

  end subroutine interpl_mid

  subroutine stagder_1st(self, delta, scheme, from_to, bc_start, bc_end, sym)
    implicit none

    class(tdsops_t), intent(inout) :: self
    real(dp), intent(in) :: delta
    character(*), intent(in) :: scheme, from_to
    integer, intent(in) :: bc_start, bc_end
    logical, optional, intent(in) :: sym

    real(dp), allocatable :: dist_b(:)
    real(dp) :: alpha, aci, bci
    integer :: i, n, n_halo

    if (self%n_halo < 2) error stop 'Staggared deriv require n_halo >= 2'

    ! alpha is alcai

    select case (scheme)
    case ('compact6')
      alpha = 9._dp/62._dp
      aci = 63._dp/62._dp/delta
      bci = 17._dp/62._dp/3._dp/delta
    case default
      error stop 'scheme is not defined'
    end select

    self%alpha = alpha
    self%a = aci; self%b = bci

    select case (from_to)
    case ('v2p')
      self%coeffs(:) = [0._dp, 0._dp, 0._dp, -bci, &
                        -aci, &
                        aci, bci, 0._dp, 0._dp]
    case ('p2v')
      self%coeffs(:) = [0._dp, 0._dp, -bci, -aci, &
                        aci, &
                        bci, 0._dp, 0._dp, 0._dp]
    end select

    do i = 1, self%n_halo
      self%coeffs_s(:, i) = self%coeffs(:)
      self%coeffs_e(:, i) = self%coeffs(:)
    end do

    self%dist_sa(:) = alpha; self%dist_sc(:) = alpha

    n = self%n_tds
    n_halo = self%n_halo

    allocate (dist_b(self%n_rhs))
    dist_b(:) = 1._dp

    if ((bc_start == BC_DIRICHLET) .or. (bc_start == BC_NEUMANN)) then
      self%dist_sa(1) = 0._dp

      select case (from_to)
      case ('v2p')
        ! sym is always .false.
        dist_b(1) = 1._dp + alpha
        self%coeffs_s(:, 1) = [0._dp, 0._dp, 0._dp, 0._dp, &
                               -aci - 2*bci, &
                               aci + bci, bci, 0._dp, 0._dp]
        self%coeffs_s(:, 2) = [0._dp, 0._dp, 0._dp, -bci, &
                               -aci, &
                               aci, bci, 0._dp, 0._dp]
      case ('p2v')
        ! sym is always .true.
        self%dist_sc(1) = 0._dp
        self%coeffs_s(:, 1) = 0._dp
        self%coeffs_s(:, 2) = [0._dp, 0._dp, 0._dp, -aci - bci, &
                               aci, &
                               bci, 0._dp, 0._dp, 0._dp]
      end select
    end if

    if ((bc_end == BC_DIRICHLET) .or. (bc_end == BC_NEUMANN)) then
      self%dist_sc(n) = 0._dp

      select case (from_to)
      case ('v2p')
        ! sym is always .false.
        dist_b(n) = 1._dp + alpha
        self%coeffs_e(:, n_halo) = 0._dp
        self%coeffs_e(:, n_halo - 1) = [0._dp, 0._dp, 0._dp, -bci, &
                                        -aci - bci, &
                                        aci + 2*bci, 0._dp, 0._dp, 0._dp]
      case ('p2v')
        ! sym is always .true.
        self%dist_sa(n) = 0._dp
        self%coeffs_e(:, n_halo) = 0._dp
        self%coeffs_e(:, n_halo - 1) = [0._dp, 0._dp, -bci, -aci, &
                                        aci + bci, &
                                        0._dp, 0._dp, 0._dp, 0._dp]
      end select
    end if

    call self%preprocess_thom(dist_b)
    call self%preprocess_dist(dist_b)

  end subroutine stagder_1st

  subroutine preprocess_dist(self, dist_b)
    implicit none

    class(tdsops_t), intent(inout) :: self

    real(dp), dimension(:), intent(in) :: dist_b

    integer :: i

    ! Ref DOI: 10.1109/MCSE.2021.3130544
    ! Algorithm 3 in page 4
    ! First two lines first
    do i = 1, 2
      self%dist_sa(i) = self%dist_sa(i)/dist_b(i)
      self%dist_sc(i) = self%dist_sc(i)/dist_b(i)
      self%dist_bw(i) = self%dist_sc(i)
      self%dist_af(i) = 1._dp/dist_b(i)
    end do

    ! Then the remaining in the forward pass
    do i = 3, self%n_tds
      ! Algorithm 3 in ref obtains 'r' coeffs on the fly in line 7.
      ! As we have to solve many RHSs with the same tridiagonal system,
      ! it is better to do a preprocessing first.
      ! So lets store 'r' coeff in dist_fw array.
      self%dist_fw(i) = 1._dp/(dist_b(i) &
                               - self%dist_sa(i)*self%dist_sc(i - 1))
      ! dist_af is 'a_i' in line 7 of Algorithm 3 in ref.
      self%dist_af(i) = self%dist_sa(i)
      ! We store a_i^* and c_i^* in dist_sa and dist_sc because
      ! we need them later in the substitution phase.
      self%dist_sa(i) = -self%dist_fw(i)*self%dist_sa(i) &
                        *self%dist_sa(i - 1)
      self%dist_sc(i) = self%dist_fw(i)*self%dist_sc(i)
    end do

    ! backward pass starting in line 12 of Algorithm 3.
    do i = self%n_tds - 2, 2, -1
      self%dist_sa(i) = self%dist_sa(i) &
                        - self%dist_sc(i)*self%dist_sa(i + 1)
      self%dist_bw(i) = self%dist_sc(i)
      self%dist_sc(i) = -self%dist_sc(i)*self%dist_sc(i + 1)
    end do

    ! Line 17 and 18 are tricky
    ! First we have a new 'r', we need it.
    ! And for 'r' we need c_0^*...
    ! Now examine closely, c_0^* is set in line 4 and never overwritten!
    ! So we can use dist_sc(1) as is in place of c_0^*.
    ! We need to store this new 'r' somewhere ...
    ! dist_fw(1) is never used, so store this extra 'r' factor here instead
    self%dist_fw(1) = 1._dp/(1._dp - self%dist_sc(1)*self%dist_sa(2))

    ! Finally Line 19 and 20 in Algorithm 3 in ref.
    self%dist_sa(1) = self%dist_fw(1)*self%dist_sa(1)
    self%dist_sc(1) = -self%dist_fw(1)*self%dist_sc(1)*self%dist_sc(2)

  end subroutine preprocess_dist

  subroutine preprocess_thom(self, b)
    implicit none

    class(tdsops_t), intent(inout) :: self
    real(dp), dimension(:), intent(in) :: b

    integer :: i, n

    n = self%n_tds

    self%thom_w = b
    self%thom_f = self%dist_sc
    if (self%periodic) then
      self%thom_w(1) = 2._dp
      self%thom_w(n) = 1._dp + self%alpha*self%alpha
    end if

    self%thom_s(1) = 0._dp
    do i = 2, n
      self%thom_s(i) = self%dist_sa(i)/self%thom_w(i - 1)
      self%thom_w(i) = self%thom_w(i) - self%thom_f(i - 1)*self%thom_s(i)
    end do
    do i = 1, n
      self%thom_w(i) = 1._dp/self%thom_w(i)
    end do

    self%thom_p = [-1._dp, (0._dp, i=2, n - 1), self%alpha]
    do i = 2, n
      self%thom_p(i) = self%thom_p(i) - self%thom_p(i - 1)*self%thom_s(i)
    end do
    self%thom_p(n) = self%thom_p(n)*self%thom_w(n)
    do i = n - 1, 1, -1
      self%thom_p(i) = self%thom_w(i)*(self%thom_p(i) &
                                       - self%thom_f(i)*self%thom_p(i + 1))
    end do

  end subroutine preprocess_thom

end module m_tdsops

