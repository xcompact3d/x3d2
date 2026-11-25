module m_abl
  use iso_fortran_env, only: int64
!! Reusable, backend-agnostic Atmospheric Boundary Layer driver.
!!
!! Ported from Incompact3d's `Case-ABL.f90` (itype=10). The driver owns the
!! ABL physics and is owned by a case; it never uses `m_solver`/`m_base_case`,
!! taking fields/mesh/coords in (same pattern as `m_ibm`).
!!
!! The neutral rough-wall model is configured here and evaluated as part of
!! the LES SGS divergence, where it can replace the bottom-plane SGS flux
!! without overwriting other momentum terms.
  use m_allocator, only: allocator_t
  use m_base_backend, only: base_backend_t
  use m_common, only: dp, pi, DIR_C, DIR_X, VERT
  use m_config, only: abl_config_t
  use m_field, only: field_t
  use m_logger, only: log_info, abort_run, str
  use m_mesh, only: mesh_t
  use m_mpi, only: MPI_COMM_WORLD, MPI_INTEGER, MPI_IN_PLACE, MPI_SUM, &
                   MPI_Allreduce, MPI_Bcast
  use m_les, only: les_t, abl_wall_sample_location

  implicit none

  private
  public :: abl_t

  type :: abl_t
    class(base_backend_t), pointer :: backend => null()
    class(mesh_t), pointer :: mesh => null()
    type(allocator_t), pointer :: host_allocator => null()
    type(abl_config_t) :: cfg
    ! Cached Rayleigh damping coefficient field (coeff*lambda(y)), built once.
    class(field_t), pointer :: damp_coeff => null()
    ! Fringe region(s), set up by configure_fringe: the source band starts
    ! at x vertex fringe_start, with weights fringe_lambda. The shifted
    ! periodic band is blended shift_offset vertices further on, the
    ! concurrent precursor band concurrent_offset vertices further on.
    integer :: fringe_start = 0
    integer :: shift_offset = 0
    integer :: concurrent_offset = 0
    integer :: fringe_nshift = 0
    real(dp), allocatable :: fringe_lambda(:)
    !> Rank holding each global z plane
    integer, allocatable :: z_plane_rank(:)
    !> Number of global x vertices in the precursor region
    integer :: n_precursor = 0
  contains
    procedure :: initialise
    procedure :: apply_forcing
    procedure :: configure_fringe
    procedure :: apply_fringe
    procedure :: configure_wall_boundary_correction
    procedure :: apply_damping
    procedure :: build_damping
  end type abl_t

  interface abl_t
    module procedure init
  end interface abl_t

contains

  function init(backend, mesh, host_allocator, cfg) result(abl)
    implicit none

    class(base_backend_t), target, intent(inout) :: backend
    type(mesh_t), target, intent(inout) :: mesh
    type(allocator_t), target, intent(inout) :: host_allocator
    type(abl_config_t), intent(in) :: cfg
    type(abl_t) :: abl

    abl%backend => backend
    abl%mesh => mesh
    abl%host_allocator => host_allocator
    abl%cfg = cfg

  end function init

  subroutine initialise(self, u, v, w)
    !! Sets the t=0 velocity profile: log-law when the drive is pressure-
    !! gradient / mass-conserve, otherwise the geostrophic wind, plus noise
    !! (Incompact3d init_abl, 68-92).
    implicit none

    class(abl_t) :: self
    class(field_t), intent(inout) :: u, v, w

    class(field_t), pointer :: hu, hv, hw
    integer :: i, j, k, dims(3), seed_size, date_values(8), base_seed, ierr
    integer, allocatable :: seed(:)
    integer(int64) :: clock, state
    real(dp) :: coords(3), y, prof, noise(3)
    logical :: log_law

    dims = self%mesh%get_dims(VERT)

    hu => self%host_allocator%get_block(DIR_C)
    hv => self%host_allocator%get_block(DIR_C)
    hw => self%host_allocator%get_block(DIR_C)

    ! As Incompact3d init_abl, each MPI rank draws its own noise stream, in
    ! u, v, w order. With seed = 0 the base seed comes from rank 0's clock,
    ! as in Incompact3d; it is reported so the run can be repeated. A fixed
    ! seed repeats the noise for the same compiler and rank count only, as
    ! random_number is processor dependent.
    base_seed = self%cfg%seed
    if (base_seed == 0) then
      if (self%mesh%par%is_root()) then
        call system_clock(count=clock)
        if (clock <= 0) then
          ! No usable clock: the standard gives -huge, some runtimes give 0
          call date_and_time(values=date_values)
          clock = ((((int(date_values(1), int64)*12 + date_values(2))*31 &
                     + date_values(3))*24 + date_values(5))*60 &
                   + date_values(6))*60000_int64 &
                  + date_values(7)*1000_int64 + date_values(8)
        end if
        base_seed = 1 + int(modulo(clock, 2147483646_int64))
      end if
      call MPI_Bcast(base_seed, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    end if
    call log_info('ABL noise seed: '//str(base_seed))

    ! Spread the base seed and rank over the processor's seed with a MINSTD
    ! generator. The integer arithmetic is done in int64 and every word lies
    ! in [1, 2**31 - 2], so nothing overflows a default integer, however many
    ! words the seed has or ranks the run uses, and it is never all zero.
    call random_seed(size=seed_size)
    allocate (seed(seed_size))
    state = 1 + modulo(int(base_seed, int64) &
                       + 63946_int64*(self%mesh%par%nrank + 1), &
                       2147483646_int64)
    do i = 1, seed_size
      state = modulo(48271_int64*state, 2147483647_int64)
      seed(i) = int(state)
    end do
    call random_seed(put=seed)

    call random_number(hu%data(1:dims(1), 1:dims(2), 1:dims(3)))
    call random_number(hv%data(1:dims(1), 1:dims(2), 1:dims(3)))
    call random_number(hw%data(1:dims(1), 1:dims(2), 1:dims(3)))

    noise = self%cfg%init_noise
    log_law = self%cfg%pressure_gradient .or. self%cfg%mass_conserve

    do k = 1, dims(3)
      do j = 1, dims(2)
        do i = 1, dims(1)
          coords = self%mesh%get_coordinates(i, j, k)
          y = coords(2)
          if (log_law) then
            prof = self%cfg%u_star/self%cfg%kappa &
                   *log((y + self%cfg%z0)/self%cfg%z0)
          else
            prof = self%cfg%u_geo(1)
          end if
          hu%data(i, j, k) = prof* &
                             (1._dp + noise(1)* &
                              (2._dp*hu%data(i, j, k) - 1._dp))
          hv%data(i, j, k) = noise(2)*(2._dp*hv%data(i, j, k) - 1._dp)
          hw%data(i, j, k) = noise(3)*(2._dp*hw%data(i, j, k) - 1._dp)
        end do
      end do
    end do

    call self%backend%set_field_data(u, hu%data)
    call self%backend%set_field_data(v, hv%data)
    call self%backend%set_field_data(w, hw%data)

    call self%host_allocator%release_block(hu)
    call self%host_allocator%release_block(hv)
    call self%host_allocator%release_block(hw)

    call u%set_data_loc(VERT)
    call v%set_data_loc(VERT)
    call w%set_data_loc(VERT)

  end subroutine initialise

  subroutine apply_forcing(self, du, dv, dw, u, v, w)
    !! Adds the ABL driving forces to the momentum RHS each substep
    !! (Incompact3d momentum_forcing_abl, 296-320).
    implicit none

    class(abl_t) :: self
    class(field_t), intent(inout) :: du, dv, dw
    class(field_t), intent(in) :: u, v, w

    real(dp) :: f

    f = self%cfg%coriolis_freq

    ! Driving force: pressure gradient, or the geostrophic-balance term when
    ! the pressure gradient is off. With a concurrent precursor these stay
    ! uniform in x, unlike Incompact3d, which zeroes the pressure gradient
    ! for x >= pdl: a force that varies in x only differs from its x mean
    ! by the gradient of a periodic potential, which the pressure absorbs,
    ! so the mask would only weaken the drive by pdl/Lx everywhere.
    if (self%cfg%pressure_gradient) then
      call self%backend%field_shift(du, self%cfg%u_star**2/self%cfg%delta)
    else if (self%cfg%coriolis) then
      call self%backend%field_shift(du, -f*self%cfg%u_geo(3))
      call self%backend%field_shift(dw, f*self%cfg%u_geo(1))
    end if

    ! Coriolis on the resolved field.
    if (self%cfg%coriolis) then
      call self%backend%vecadd(f, w, 1._dp, du)
      call self%backend%vecadd(-f, u, 1._dp, dw)
    end if

    ! Rayleigh damping layer near the domain top.
    if (self%cfg%damping) then
      call self%apply_damping(du, dv, dw, u, v, w)
    end if

  end subroutine apply_forcing

  subroutine apply_damping(self, du, dv, dw, u, v, w)
    !! Rayleigh sponge relaxing the flow toward the reference profile over the
    !! top damping layer (Incompact3d damping_zone, neutral branch).
    !!
    !! As in Incompact3d, the target is the same whatever drives the flow:
    !! u relaxes to the log-law value at delta, and v, w to u_geo(2), u_geo(3).
    !! That suits the pressure-gradient and mass-conserve drives, whose
    !! profile is the log law. With only the Coriolis drive the free stream
    !! is u_geo, so the sponge pulls u toward a different value than the
    !! drive does. The layer is centred on delta and 0.1*delta thick, so it
    !! acts only when the domain extends above delta - 0.05*delta.
    implicit none

    class(abl_t) :: self
    class(field_t), intent(inout) :: du, dv, dw
    class(field_t), intent(in) :: u, v, w

    class(field_t), pointer :: tmp
    real(dp) :: u_ref

    if (.not. associated(self%damp_coeff)) call self%build_damping()

    u_ref = self%cfg%u_star/self%cfg%kappa*log(self%cfg%delta/self%cfg%z0)

    tmp => self%backend%allocator%get_block(DIR_X, VERT)

    call self%backend%veccopy(tmp, u)
    call self%backend%field_shift(tmp, -u_ref)
    call self%backend%vecmult(tmp, self%damp_coeff)
    call self%backend%vecadd(-1._dp, tmp, 1._dp, du)

    call self%backend%veccopy(tmp, v)
    call self%backend%field_shift(tmp, -self%cfg%u_geo(2))
    call self%backend%vecmult(tmp, self%damp_coeff)
    call self%backend%vecadd(-1._dp, tmp, 1._dp, dv)

    call self%backend%veccopy(tmp, w)
    call self%backend%field_shift(tmp, -self%cfg%u_geo(3))
    call self%backend%vecmult(tmp, self%damp_coeff)
    call self%backend%vecadd(-1._dp, tmp, 1._dp, dw)

    call self%backend%allocator%release_block(tmp)

  end subroutine apply_damping

  subroutine build_damping(self)
    !! Precomputes the damping coefficient coeff*lambda(y), constant in x/z,
    !! or zero in the successor region of a concurrent precursor, as in
    !! Incompact3d. The damping depends on y and on the flow, so unlike an
    !! x-only body force its mask does act locally.
    implicit none

    class(abl_t) :: self

    class(field_t), pointer :: h
    integer :: i, j, k, dims(3)
    real(dp) :: coords(3), x, y, coeff, dheight, ylo, yhi, lambda
    ! Nondimensional relaxation rate, as hard-coded in Incompact3d's
    ! damping_zone (wvar = fifteen). It is empirical: there it replaced a
    ! commented-out formula that is singular when the domain height equals
    ! delta. The sponge rate is then 15*u_star/delta, i.e. a timescale of
    ! delta/(15*u_star), about 150 s for the Deskos (2019) case.
    real(dp), parameter :: wvar = 15._dp

    self%damp_coeff => self%backend%allocator%get_block(DIR_X, VERT)
    call self%damp_coeff%set_data_loc(VERT)

    h => self%host_allocator%get_block(DIR_C)
    dims = self%mesh%get_dims(VERT)

    dheight = 0.1_dp*self%cfg%delta
    coeff = wvar*self%cfg%u_star/self%cfg%delta
    ylo = self%cfg%delta - 0.5_dp*dheight
    yhi = self%cfg%delta + 0.5_dp*dheight

    do k = 1, dims(3)
      do j = 1, dims(2)
        do i = 1, dims(1)
          coords = self%mesh%get_coordinates(i, j, k)
          y = coords(2)
          if (y >= yhi) then
            lambda = 1._dp
          else if (y >= ylo) then
            lambda = 0.5_dp*(1._dp - cos(pi*(y - ylo)/dheight))
          else
            lambda = 0._dp
          end if
          ! Incompact3d damping_zone: no damping where x >= pdl
          x = real(self%mesh%par%n_offset(1) + i - 1, dp)*self%mesh%geo%d(1)
          if (self%cfg%concurrent_precursor .and. &
              x >= self%cfg%precursor_length) lambda = 0._dp
          h%data(i, j, k) = coeff*lambda
        end do
      end do
    end do

    call self%backend%set_field_data(self%damp_coeff, h%data)
    call self%host_allocator%release_block(h)

  end subroutine build_damping

  subroutine configure_wall_boundary_correction(self, les)
    !! Configure the neutral ABL wall model in the LES closure.
    !!
    !! The floor is no-slip, so the velocity driving the drag law is sampled
    !! at dsampling*dy above it, matching Incompact3d's wall_sgs_noslip.
    implicit none

    class(abl_t), intent(in) :: self
    type(les_t), intent(inout) :: les
    real(dp) :: sampling_height, sample_weight
    integer :: dims(3), sample_plane

    ! The wall stress is written at local y-plane 2 and sampled at a local
    ! plane, so the wall must sit on the first plane of every rank.
    if (self%mesh%par%nproc_dir(2) /= 1) &
      call abort_run('ABL wall model does not support a decomposed y &
                     &direction.')

    dims = self%mesh%get_dims(VERT)
    if (dims(2) < 3) &
      call abort_run('ABL wall model requires at least 3 y vertices.')

    call abl_wall_sample_location(self%mesh, self%cfg%dsampling, &
                                  sample_plane, sample_weight, &
                                  sampling_height)
    if (sampling_height <= self%cfg%z0) &
      call abort_run('ABL wall-model sampling height must exceed z0.')

    call les%configure_abl_wall_boundary( &
      self%cfg%kappa, self%cfg%z0, sampling_height, sample_plane, &
      sample_weight)
  end subroutine configure_wall_boundary_correction

  subroutine configure_fringe(self)
    !! Set up Incompact3d's fringe_region for the shifted periodic
    !! (ishiftedper) and concurrent precursor (iconcprec) options, alone or
    !! together. The ramp and band sizes follow the legacy code exactly.
    !!
    !! Shifted periodic alone: over the last third of the periodic-x domain
    !! the flow is blended towards the flow a third of the domain upstream,
    !! shifted spanwise by Lz/8, so streaks cannot lock onto the periodic box.
    !!
    !! Concurrent precursor: the domain holds a precursor region
    !! x < precursor_length and a successor region behind it. The end of the
    !! successor region is blended back into the end of the precursor
    !! region, so the precursor's inlet (the domain outlet, through the
    !! periodic boundary) never sees the successor's flow. With shifted
    !! periodic as well, the precursor region has its own shifted periodic
    !! fringe, and the successor band takes the shifted flow too.
    !!
    !! The flow rate through every x plane is the same, so drag in the
    !! successor region (turbines) would slow the precursor under a fixed
    !! force. mass_conserve holds the flow rate instead, which keeps the
    !! precursor region equivalent to a standalone periodic ABL; the config
    !! requires it with a concurrent precursor.
    implicit none

    class(abl_t), intent(inout) :: self

    real(dp), allocatable :: lambda(:)
    real(dp) :: lx, dx, x, pdl, fre, frs, frd
    integer :: global_dims(3), dims(3), nx, nz, nfe, npe, i, i_first
    integer :: z_offset, ierr
    logical :: shifted, concurrent

    shifted = self%cfg%shifted_periodic
    concurrent = self%cfg%concurrent_precursor
    if (.not. (shifted .or. concurrent)) return

    if (self%mesh%par%nproc_dir(1) /= 1) &
      call abort_run('ABL shifted_periodic and concurrent_precursor require &
                     &an undecomposed x direction.')

    global_dims = self%mesh%get_global_dims(VERT)
    nx = global_dims(1)
    nz = global_dims(3)
    lx = self%mesh%geo%L(1)
    dx = self%mesh%geo%d(1)
    pdl = self%cfg%precursor_length

    if (shifted) then
      self%fringe_nshift = int(self%mesh%geo%L(3)/8._dp/self%mesh%geo%d(3))
      if (self%fringe_nshift < 1) &
        call abort_run('ABL shifted_periodic needs at least 8 z points.')
    else
      self%fringe_nshift = 0
    end if

    if (concurrent) then
      if (pdl >= lx) &
        call abort_run('ABL precursor_length must be less than the domain &
                       &length in x.')
      ! Vertices with x = (i - 1)*dx < precursor_length
      self%n_precursor = 0
      do i = 1, nx
        if (real(i - 1, dp)*dx < pdl) self%n_precursor = i
      end do
    end if

    ! Legacy edl = 2/3 and frl = 1/6
    if (shifted .and. .not. concurrent) then
      fre = lx*2._dp/3._dp
      npe = nx
    else if (shifted .and. concurrent) then
      fre = pdl*2._dp/3._dp
      npe = int(nx*pdl/lx)
    else
      fre = pdl
      npe = 0
    end if
    frs = fre - fre/6._dp
    frd = fre - frs
    nfe = int(nx*fre/lx)
    allocate (lambda(nfe))
    do i = 1, nfe
      x = real(i - 1, dp)*dx
      if (x < frs) then
        lambda(i) = 0._dp
      else if (x < fre - frd/4._dp) then
        lambda(i) = 0.5_dp*(1._dp - cos(4._dp*pi/3._dp*(x - frs)/frd))
      else if (x < fre) then
        lambda(i) = 1._dp
      else
        lambda(i) = 0._dp
      end if
    end do
    ! A loop rather than findloc, which nvfortran lacks for logical arrays
    i_first = 0
    do i = nfe, 1, -1
      if (lambda(i) > 0._dp) i_first = i
    end do
    if (i_first == 0) &
      call abort_run('ABL fringe region is narrower than one grid spacing.')
    self%fringe_start = i_first
    self%fringe_lambda = lambda(i_first:nfe)

    ! The blended bands must lie downstream of the source band, which is
    ! read before it is blended, and the concurrent precursor band in the
    ! successor region
    if (shifted) then
      self%shift_offset = npe - nfe
      if (i_first + self%shift_offset <= nfe) &
        call abort_run('ABL shifted periodic fringe overlaps its source.')
    end if
    if (concurrent) then
      self%concurrent_offset = nx - nfe
      if (i_first + self%concurrent_offset <= self%n_precursor) &
        call abort_run('ABL concurrent precursor fringe reaches into the &
                       &precursor region: shorten precursor_length.')
    end if

    ! x and y are whole on every rank, so each rank holds one z block
    dims = self%mesh%get_dims(VERT)
    allocate (self%z_plane_rank(nz))
    self%z_plane_rank = 0
    z_offset = self%mesh%par%n_offset(3)
    self%z_plane_rank(z_offset + 1:z_offset + dims(3)) = self%mesh%par%nrank
    call MPI_Allreduce(MPI_IN_PLACE, self%z_plane_rank, nz, MPI_INTEGER, &
                       MPI_SUM, MPI_COMM_WORLD, ierr)

    if (shifted) &
      call log_info('ABL shifted periodic: spanwise shift of ' &
                    //str(self%fringe_nshift)//' planes, x vertices ' &
                    //str(i_first)//'-'//str(nfe)//' blended into ' &
                    //str(i_first + self%shift_offset)//'-' &
                    //str(nfe + self%shift_offset))

    if (concurrent) &
      call log_info('ABL concurrent precursor: x vertices 1-' &
                    //str(self%n_precursor)//' damped, x vertices ' &
                    //str(i_first)//'-'//str(nfe)//' blended into ' &
                    //str(i_first + self%concurrent_offset)//'-'//str(nx))
  end subroutine configure_fringe

  subroutine apply_fringe(self, u, v, w)
    !! Apply the fringe region(s) to the velocity, at the start of each
    !! stage as Incompact3d does from boundary_conditions_abl.
    !!
    !! Incompact3d blends both bands from a copy of the field taken before
    !! either blend. Here the bands are blended one after the other, which is
    !! the same because configure_fringe checks that neither blended band
    !! overlaps the source band.
    implicit none

    class(abl_t) :: self
    class(field_t), intent(inout) :: u, v, w

    call blend(u)
    call blend(v)
    call blend(w)

  contains

    subroutine blend(f)
      class(field_t), intent(inout) :: f

      if (self%cfg%shifted_periodic) &
        call self%backend%field_shift_blend_x_band( &
        f, self%fringe_lambda, self%fringe_start, self%shift_offset, &
        self%fringe_nshift, self%z_plane_rank)
      if (self%cfg%concurrent_precursor) &
        call self%backend%field_shift_blend_x_band( &
        f, self%fringe_lambda, self%fringe_start, self%concurrent_offset, &
        self%fringe_nshift, self%z_plane_rank)
    end subroutine blend

  end subroutine apply_fringe

end module m_abl
