module m_config
  !! Contains all the namelist configurations used in x3d2
  use iso_fortran_env, only: iostat_end
  use m_common
  use m_logger, only: log_warn, abort_run

  implicit none

  integer, parameter :: n_species_max = 99

  !! Maximum number of additional snapshot output fields accepted from the
  !! input file. Some entries, such as ``species``, can expand to multiple
  !! written snapshot fields.
  integer, parameter :: MAX_OUTPUT_FIELDS = 10

  type, abstract :: base_config_t
    !! All config types have a method read to initialise their data
  contains
    procedure(read), deferred :: read
  end type base_config_t

  type, extends(base_config_t) :: domain_config_t
    character(len=30) :: flow_case_name
    real(dp) :: L_global(3)
    integer :: dims_global(3), nproc_dir(3)
    character(len=20) :: BC_x(2), BC_y(2), BC_z(2)
    character(len=20) :: stretching(3)
    real(dp) :: beta(3)
  contains
    procedure :: read => read_domain_nml
  end type domain_config_t

  type, extends(base_config_t) :: solver_config_t
    real(dp) :: Re, dt
    logical :: ibm_on
    real(dp), dimension(:), allocatable :: pr_species
    integer :: n_iters, n_output, n_species
    logical :: lowmem_transeq, lowmem_fft
    logical :: spatial_filter          !! explicit low-pass filter on velocity
    real(dp) :: filter_alpha           !! filter parameter, -0.5 < alpha < 0.5
    character(3) :: poisson_solver_type, time_intg
    character(30) :: der1st_scheme, der2nd_scheme, &
                     interpl_scheme, stagder_scheme
    real(dp) :: c_nu, nu0_nu !! hyperviscous (iSVV) der2nd parameters
  contains
    procedure :: read => read_solver_nml
  end type solver_config_t

  type, extends(base_config_t) :: les_config_t
    !! Configuration for explicit sub-grid-scale modelling.
    character(len=20) :: model = 'none'
    real(dp) :: smagorinsky_constant = 0.14_dp
    logical :: wall_damping = .false.
    real(dp) :: wall_damping_n = 3._dp
  contains
    procedure :: read => read_les_nml
  end type les_config_t

  type, extends(base_config_t) :: channel_config_t
    real(dp) :: omega_rot
    real(dp) :: init_noise(3)
    real(dp) :: inlet_noise(3)
    logical :: rotation
    integer :: n_rotate
  contains
    procedure :: read => read_channel_nml
  end type channel_config_t

  type, extends(base_config_t) :: cylinder_config_t
    real(dp) :: init_noise(3)
    real(dp) :: inlet_noise(3)
  contains
    procedure :: read => read_cylinder_nml
  end type cylinder_config_t

  type, extends(base_config_t) :: wind_turbine_config_t
    !! Configuration for the wind-turbine carrier case. Inflow/outflow and
    !! initial noise mirror the cylinder case; the turbine block selects and
    !! drives the turbine forcing model.
    real(dp) :: init_noise(3) = 0._dp
    real(dp) :: inlet_noise(3) = 0._dp  !! perturbation applied at inlet each substep
    real(dp) :: bc_start_u = 1._dp   !! uniform inflow u (Dirichlet at i=1)
    real(dp) :: bc_start_v = 0._dp
    real(dp) :: bc_start_w = 0._dp
    integer :: iturbine = 0          !! 0: none, 1: actuator line, 2: disc
    integer :: iturboutput = 1       !! turbine diagnostics output frequency
    character(len=256) :: adm_coords = '' !! actuator-disc (.ad) coords file
    real(dp) :: rho_air = 1._dp      !! air density used in thrust/power
    real(dp) :: T_relax = -1._dp     !! ADM velocity-filter time; < 0 disables
    !> When false, a missing wind_turbine_nml keeps the defaults (no
    !> turbines), for a case that carries turbines optionally (the ABL)
    logical :: required = .true.
  contains
    procedure :: read => read_wind_turbine_nml
  end type wind_turbine_config_t

  type, extends(base_config_t) :: abl_config_t
    real(dp) :: z0 = 0.1_dp            !! aerodynamic roughness length
    real(dp) :: u_star = 0._dp         !! friction velocity
    real(dp) :: delta = 1._dp          !! boundary-layer depth
    real(dp) :: kappa = 0.41_dp        !! von Karman constant
    real(dp) :: dsampling = 3._dp      !! wall-model sampling height, in dy
    real(dp) :: u_geo(3) = 0._dp       !! geostrophic wind UG
    real(dp) :: coriolis_freq = 0._dp  !! Coriolis frequency f
    real(dp) :: init_noise(3) = 0._dp
    !> Seed for the initial noise; 0 draws one from the clock
    integer :: seed = 0
    real(dp) :: u_bulk = 0._dp         !! target bulk velocity (mass_conserve)
    integer :: profile_start_iter = -1 !! first profile sample (-1 = disabled)
    character(len=256) :: profile_file = 'abl_profile.csv'
    logical :: pressure_gradient = .false.
    logical :: coriolis = .false.
    logical :: mass_conserve = .false.
    logical :: damping = .false.
    !> Blend the outlet towards a spanwise-shifted copy of the flow
    !> (Incompact3d ishiftedper); periodic x only
    logical :: shifted_periodic = .false.
  contains
    procedure :: read => read_abl_nml
  end type abl_config_t

  type, extends(base_config_t) :: inflow_config_t
    !! Precursor inflow planes (docs/source/user/input_file.rst)
    logical :: write_planes = .false.  !! sample and write y-z planes
    !> Global x vertex index of the sampling plane; 0 selects the last one
    integer :: plane_x_index = 0
    integer :: ntimesteps = 100        !! planes per file
    integer :: sample_every = 1        !! iterations between samples
    integer :: start_iter = 0          !! first sample at start_iter + sample_every
    character(len=256) :: inflow_prefix = 'inflow'
    !> Successor: impose planes read from inflow_read_prefix at the x inlet
    logical :: read_planes = .false.
    integer :: ninflows = 1            !! files cycled through by a successor
    character(len=256) :: inflow_read_prefix = 'inflow'
  contains
    procedure :: read => read_inflow_nml
  end type inflow_config_t

  type, extends(base_config_t) :: stats_config_t
    integer :: initstat = 0          !! iteration to start accumulating (0 = disabled)
    integer :: istatfreq = 1         !! accumulate every N steps
    integer :: istatout = 0          !! write stats every N steps (0 = disabled)
    character(len=256) :: stats_prefix = "statistics"
  contains
    procedure :: read => read_stats_nml
  end type stats_config_t

  type, extends(base_config_t) :: checkpoint_config_t
    integer :: checkpoint_freq = 0                         !! Frequency of checkpointing (0 = off)
    integer :: snapshot_freq = 0                           !! Frequency of snapshots (0 = off)
    logical :: keep_checkpoint = .true.                    !! If false, only keep latest checkpoint
    character(len=256) :: checkpoint_prefix = "checkpoint"
    character(len=256) :: snapshot_prefix = "snapshot"
    logical :: restart_from_checkpoint = .false.
    character(len=256) :: restart_file = ""
    integer, dimension(3) :: output_stride = [2, 2, 2]     !! Spatial stride for snapshot output
    logical :: snapshot_sp = .false.                       !! if true, snapshot in single precision
    character(len=32) :: output_fields(MAX_OUTPUT_FIELDS) = '' !! additional snapshot output fields
  contains
    procedure :: read => read_checkpoint_nml
  end type checkpoint_config_t

  abstract interface
    subroutine read(self, nml_file, nml_string) !&
      !! Assigns the member variables either from a file or text source.
      !!
      !! nml_file can be an absolute or relative path
      !! nml_string is a character string that contains the namelist.
      !! For example, nml_string="&foobar_nml foo=0, bar='this'/"
      import :: base_config_t

      class(base_config_t) :: self
      character(*), optional, intent(in) :: nml_file
      character(*), optional, intent(in) :: nml_string
    end subroutine read
  end interface

contains

  subroutine read_domain_nml(self, nml_file, nml_string)
    implicit none

    class(domain_config_t) :: self
    character(*), optional, intent(in) :: nml_file
    character(*), optional, intent(in) :: nml_string

    integer :: unit

    character(len=20) :: flow_case_name
    real(dp), dimension(3) :: L_global
    integer, dimension(3) :: dims_global
    integer, dimension(3) :: nproc_dir
    character(len=20) :: BC_x(2), BC_y(2), BC_z(2)
    character(len=20) :: stretching(3) = ['uniform', 'uniform', 'uniform']
    real(dp), dimension(3) :: beta

    namelist /domain_settings/ flow_case_name, L_global, dims_global, &
      nproc_dir, BC_x, BC_y, BC_z, stretching, beta

    if (present(nml_file) .and. present(nml_string)) then
      call abort_run('Reading domain config failed! Provide only a file name &
                     &or source, not both.', collective=.true.)
    else if (present(nml_file)) then
      open (newunit=unit, file=nml_file)
      read (unit, nml=domain_settings)
      close (unit)
    else if (present(nml_string)) then
      read (nml_string, nml=domain_settings)
    else
      call abort_run('Reading domain config failed! Provide at least one of &
                     &the following: file name or source', collective=.true.)
    end if

    self%flow_case_name = flow_case_name
    self%L_global = L_global
    self%dims_global = dims_global
    self%nproc_dir = nproc_dir
    self%BC_x = BC_x
    self%BC_y = BC_y
    self%BC_z = BC_z
    self%stretching = stretching
    self%beta = beta

  end subroutine read_domain_nml

  subroutine read_solver_nml(self, nml_file, nml_string)
    implicit none

    class(solver_config_t) :: self
    character(*), optional, intent(in) :: nml_file
    character(*), optional, intent(in) :: nml_string

    integer :: unit

    real(dp) :: Re, dt
    logical :: ibm_on = .false.
    real(dp), dimension(n_species_max) :: pr_species = 1._dp
    integer :: n_iters, n_output, n_species = 0
    !> triggers the low memory implementations
    logical :: lowmem_transeq = .false., lowmem_fft = .false.
    !! Explicit low-pass filter on the velocity, off by default. Needed by
    !! wall-modelled ABL runs, where the 2*dx mode the compact schemes cannot
    !! dissipate otherwise breaks momentum conservation in the convection.
    logical :: spatial_filter = .false.
    !! Incompact3d's value. Smaller alpha is far more dissipative: at 0.40 the
    !! filter damps every resolved wavenumber about ten times harder per
    !! application than at 0.49, which over a long run scrubs out the smallest
    !! resolved eddies. The filter's system is only marginally diagonally
    !! dominant here, so it is solved with Thomas rather than the distributed
    !! algorithm wherever the direction is not decomposed.
    real(dp) :: filter_alpha = 0.49_dp
    character(3) :: time_intg
    character(3) :: poisson_solver_type = 'FFT'
    character(30) :: der1st_scheme = 'compact6', der2nd_scheme = 'compact6', &
                     interpl_scheme = 'classic', stagder_scheme = 'compact6'
    !> hyperviscous (iSVV) der2nd parameters; only used when der2nd_scheme is
    !> 'compact6-hyperviscous' (Incompact3d cnu, nu0nu). Harmless defaults otherwise.
    real(dp) :: c_nu = 0.44_dp, nu0_nu = 0._dp

    namelist /solver_params/ Re, dt, n_iters, n_output, poisson_solver_type, &
      spatial_filter, filter_alpha, &
      n_species, pr_species, lowmem_transeq, lowmem_fft, &
      time_intg, der1st_scheme, der2nd_scheme, interpl_scheme, &
      stagder_scheme, ibm_on, c_nu, nu0_nu

    if (present(nml_file) .and. present(nml_string)) then
      call abort_run('Reading solver config failed! Provide only a file name &
                     &or source, not both.', collective=.true.)
    else if (present(nml_file)) then
      open (newunit=unit, file=nml_file)
      read (unit, nml=solver_params)
      close (unit)
    else if (present(nml_string)) then
      read (nml_string, nml=solver_params)
    else
      call abort_run('Reading solver config failed! Provide at least one of &
                     &the following: file name or source', collective=.true.)
    end if

    self%Re = Re
    self%dt = dt
    self%n_iters = n_iters
    self%n_output = n_output
    self%spatial_filter = spatial_filter
    self%filter_alpha = filter_alpha
    self%ibm_on = ibm_on
    self%n_species = n_species
    if (n_species > 0) self%pr_species = pr_species(1:n_species)
    self%lowmem_transeq = lowmem_transeq
    self%lowmem_fft = lowmem_fft
    self%poisson_solver_type = poisson_solver_type
    self%time_intg = time_intg
    self%der1st_scheme = der1st_scheme
    self%der2nd_scheme = der2nd_scheme
    self%interpl_scheme = interpl_scheme
    self%stagder_scheme = stagder_scheme
    self%c_nu = c_nu
    self%nu0_nu = nu0_nu

  end subroutine read_solver_nml

  subroutine read_les_nml(self, nml_file, nml_string)
    class(les_config_t) :: self
    character(*), optional, intent(in) :: nml_file
    character(*), optional, intent(in) :: nml_string

    integer :: unit, ierr
    character(len=20) :: model
    real(dp) :: smagorinsky_constant, wall_damping_n
    logical :: wall_damping

    namelist /les_params/ model, smagorinsky_constant, wall_damping, &
      wall_damping_n

    model = self%model
    smagorinsky_constant = self%smagorinsky_constant
    wall_damping = self%wall_damping
    wall_damping_n = self%wall_damping_n

    if (present(nml_file) .and. present(nml_string)) then
      call abort_run('Reading LES config failed!  Provide only a file name or &
                     &source.', collective=.true.)
    else if (present(nml_file)) then
      open (newunit=unit, file=nml_file, iostat=ierr)
      if (ierr /= 0) &
        call abort_run('Opening LES config file failed.', collective=.true.)
      read (unit, nml=les_params, iostat=ierr)
      close (unit)
      ! Existing input files may omit this optional namelist. End-of-file
      ! therefore selects the defaults, while a malformed block remains fatal.
      if (ierr > 0) &
        call abort_run('Reading LES config failed.', collective=.true.)
    else if (present(nml_string)) then
      read (nml_string, nml=les_params)
    else
      call abort_run('Reading LES config failed! Provide a file name or &
                     &source.', collective=.true.)
    end if

    select case (trim(model))
    case ('none', 'smagorinsky')
    case default
      call abort_run('Unknown LES model. Use "none" or "smagorinsky".', &
                     collective=.true.)
    end select
    if (smagorinsky_constant <= 0._dp) &
      call abort_run('smagorinsky_constant must be positive.', &
                     collective=.true.)
    if (wall_damping_n <= 0._dp) &
      call abort_run('wall_damping_n must be positive.', collective=.true.)

    self%model = trim(model)
    self%smagorinsky_constant = smagorinsky_constant
    self%wall_damping = wall_damping
    self%wall_damping_n = wall_damping_n
  end subroutine read_les_nml

  subroutine read_channel_nml(self, nml_file, nml_string)
    implicit none

    class(channel_config_t) :: self
    character(*), optional, intent(in) :: nml_file
    character(*), optional, intent(in) :: nml_string

    integer :: unit

    real(dp) :: init_noise(3)
    real(dp) :: inlet_noise(3)
    real(dp) :: omega_rot
    logical :: rotation
    integer :: n_rotate

    namelist /channel_nml/ init_noise, inlet_noise, &
      rotation, omega_rot, n_rotate

    ! Default to no noise if the namelist omits these entries.
    init_noise = 0._dp
    inlet_noise = 0._dp

    if (present(nml_file) .and. present(nml_string)) then
      call abort_run('Reading channel config failed! Provide only a file name &
                     &or source, not both.', collective=.true.)
    else if (present(nml_file)) then
      open (newunit=unit, file=nml_file)
      read (unit, nml=channel_nml)
      close (unit)
    else if (present(nml_string)) then
      read (nml_string, nml=channel_nml)
    else
      call abort_run('Reading channel config failed! Provide at least one of &
                     &the following: file name or source', collective=.true.)
    end if

    self%init_noise = init_noise
    self%inlet_noise = inlet_noise
    self%rotation = rotation
    self%omega_rot = omega_rot
    self%n_rotate = n_rotate

  end subroutine read_channel_nml

  subroutine read_cylinder_nml(self, nml_file, nml_string)
    implicit none

    class(cylinder_config_t) :: self
    character(*), optional, intent(in) :: nml_file
    character(*), optional, intent(in) :: nml_string

    integer :: unit

    real(dp) :: init_noise(3)
    real(dp) :: inlet_noise(3)

    namelist /cylinder_nml/ init_noise, inlet_noise

    ! Default to no noise if the namelist omits these entries.
    init_noise = 0._dp
    inlet_noise = 0._dp

    if (present(nml_file) .and. present(nml_string)) then
      call abort_run('Reading cylinder config failed! Provide only a file &
                     &name or source, not both.', collective=.true.)
    else if (present(nml_file)) then
      open (newunit=unit, file=nml_file)
      read (unit, nml=cylinder_nml)
      close (unit)
    else if (present(nml_string)) then
      read (nml_string, nml=cylinder_nml)
    else
      call abort_run('Reading cylinder config failed! Provide at least one of &
                     &the following: file name or source', collective=.true.)
    end if

    self%init_noise = init_noise
    self%inlet_noise = inlet_noise

  end subroutine read_cylinder_nml

  subroutine read_wind_turbine_nml(self, nml_file, nml_string)
    implicit none

    class(wind_turbine_config_t) :: self
    character(*), optional, intent(in) :: nml_file
    character(*), optional, intent(in) :: nml_string

    integer :: unit, ierr
    character(len=256) :: msg

    real(dp) :: init_noise(3) = 0._dp
    real(dp) :: inlet_noise(3) = 0._dp
    real(dp) :: bc_start_u = 1._dp
    real(dp) :: bc_start_v = 0._dp
    real(dp) :: bc_start_w = 0._dp
    integer :: iturbine = 0
    integer :: iturboutput = 1
    character(len=256) :: adm_coords = ''
    real(dp) :: rho_air = 1._dp
    real(dp) :: T_relax = -1._dp

    namelist /wind_turbine_nml/ init_noise, inlet_noise, bc_start_u, &
      bc_start_v, bc_start_w, iturbine, iturboutput, adm_coords, &
      rho_air, T_relax

    if (present(nml_file) .and. present(nml_string)) then
      call abort_run('Reading wind_turbine config failed! Provide only a file &
                     &name or source, not both.', collective=.true.)
    else if (present(nml_file)) then
      open (newunit=unit, file=nml_file)
      read (unit, nml=wind_turbine_nml, iostat=ierr, iomsg=msg)
      close (unit)
      ! An end of file means the namelist is absent, unless its group is in
      ! the file and only the closing / is missing
      if (ierr == iostat_end .and. .not. self%required) then
        if (nml_group_in_file(nml_file, 'wind_turbine_nml')) &
          call abort_run('Reading wind_turbine_nml failed: the file ends &
                         &before its closing /.', collective=.true.)
      else if (ierr /= 0) then
        call abort_run('Reading wind_turbine_nml failed: '//trim(msg), &
                       collective=.true.)
      end if
    else if (present(nml_string)) then
      read (nml_string, nml=wind_turbine_nml)
    else
      call abort_run('Reading wind_turbine config failed! Provide at least &
                     &one of the following: file name or source', &
                     collective=.true.)
    end if

    self%init_noise = init_noise
    self%inlet_noise = inlet_noise
    self%bc_start_u = bc_start_u
    self%bc_start_v = bc_start_v
    self%bc_start_w = bc_start_w
    self%iturbine = iturbine
    self%iturboutput = iturboutput
    self%adm_coords = adm_coords
    self%rho_air = rho_air
    self%T_relax = T_relax

    if (self%iturboutput < 0) &
      call abort_run('wind_turbine_nml: iturboutput must be non-negative.', &
                     collective=.true.)
    if (self%iturbine == 2 .and. len_trim(self%adm_coords) == 0) &
      call abort_run('wind_turbine_nml: iturbine = 2 needs adm_coords.', &
                     collective=.true.)

  end subroutine read_wind_turbine_nml

  subroutine read_abl_nml(self, nml_file, nml_string)
    implicit none

    class(abl_config_t) :: self
    character(*), optional, intent(in) :: nml_file
    character(*), optional, intent(in) :: nml_string

    integer :: unit

    real(dp) :: z0, u_star, delta, kappa, dsampling
    real(dp) :: u_geo(3), coriolis_freq, init_noise(3), u_bulk
    integer :: profile_start_iter, seed
    character(len=256) :: profile_file
    logical :: pressure_gradient, coriolis, mass_conserve, damping
    logical :: shifted_periodic

    namelist /abl_nml/ z0, u_star, delta, kappa, dsampling, &
      u_geo, coriolis_freq, &
      init_noise, seed, u_bulk, profile_start_iter, profile_file, &
      pressure_gradient, coriolis, mass_conserve, damping, shifted_periodic

    ! Defaults
    z0 = 0.1_dp
    u_star = 0._dp
    delta = 1._dp
    kappa = 0.41_dp
    dsampling = 3._dp
    u_geo = 0._dp
    coriolis_freq = 0._dp
    init_noise = 0._dp
    seed = 0
    u_bulk = 0._dp
    profile_start_iter = -1
    profile_file = 'abl_profile.csv'
    pressure_gradient = .false.
    coriolis = .false.
    mass_conserve = .false.
    damping = .false.
    shifted_periodic = .false.

    if (present(nml_file) .and. present(nml_string)) then
      call abort_run('Reading ABL config failed! Provide only a file name or &
                     &source, not both.', collective=.true.)
    else if (present(nml_file)) then
      open (newunit=unit, file=nml_file)
      read (unit, nml=abl_nml)
      close (unit)
    else if (present(nml_string)) then
      read (nml_string, nml=abl_nml)
    else
      call abort_run('Reading ABL config failed! Provide at least one of the &
                     &following: file name or source', collective=.true.)
    end if

    self%z0 = z0
    self%u_star = u_star
    self%delta = delta
    self%kappa = kappa
    self%dsampling = dsampling
    self%u_geo = u_geo
    self%coriolis_freq = coriolis_freq
    self%init_noise = init_noise
    self%seed = seed
    self%u_bulk = u_bulk
    self%profile_start_iter = profile_start_iter
    self%profile_file = trim(profile_file)
    self%pressure_gradient = pressure_gradient
    self%coriolis = coriolis
    self%mass_conserve = mass_conserve
    self%damping = damping
    self%shifted_periodic = shifted_periodic

    call validate_abl_config(self)

  end subroutine read_abl_nml

  subroutine validate_abl_config(config)
    !! Validate relationships between ABL namelist values in one place.
    type(abl_config_t), intent(in) :: config

    ! The driving mechanism is checked by the case, which may instead be
    ! driven by a precursor inflow.

    if (config%z0 <= 0._dp .or. config%delta <= config%z0 .or. &
        config%kappa <= 0._dp) &
      call abort_run('ABL config error: require delta > z0 > 0 and &
                     &kappa > 0.', collective=.true.)

    if (any(config%init_noise < 0._dp)) &
      call abort_run('ABL config error: init_noise must not be negative.', &
                     collective=.true.)
    if (config%seed < 0) &
      call abort_run('ABL config error: seed must not be negative.', &
                     collective=.true.)

    if (config%profile_start_iter < -1) &
      call abort_run('ABL config error: profile_start_iter must be -1 or &
                     &greater.', collective=.true.)
    if (config%profile_start_iter >= 0 .and. &
        len_trim(config%profile_file) == 0) &
      call abort_run('ABL config error: profile_file must not be empty.', &
                     collective=.true.)
  end subroutine validate_abl_config

  subroutine read_checkpoint_nml(self, nml_file, nml_string)
    implicit none

    class(checkpoint_config_t) :: self
    character(*), optional, intent(in) :: nml_file
    character(*), optional, intent(in) :: nml_string

    integer :: unit, ierr

    integer :: checkpoint_freq = 0
    integer :: snapshot_freq = 0
    logical :: keep_checkpoint = .false.
    character(len=256) :: checkpoint_prefix = "checkpoint"
    character(len=256) :: snapshot_prefix = "snapshot"
    logical :: restart_from_checkpoint = .false.
    character(len=256) :: restart_file = ""
    integer, dimension(3) :: output_stride = [1, 1, 1]
    logical :: snapshot_sp = .false.
    character(len=32) :: output_fields(MAX_OUTPUT_FIELDS) = ''

    namelist /checkpoint_params/ checkpoint_freq, snapshot_freq, &
      keep_checkpoint, checkpoint_prefix, snapshot_prefix, &
      restart_from_checkpoint, restart_file, output_stride, snapshot_sp, &
      output_fields
    if (present(nml_file) .and. present(nml_string)) then
      call abort_run('Reading checkpoint config failed! Provide only a file &
                     &name or source, not both.', collective=.true.)
    else if (present(nml_file)) then
      open (newunit=unit, file=nml_file, iostat=ierr)
      if (ierr == 0) then
        read (unit, nml=checkpoint_params, iostat=ierr)

        if (ierr /= 0 .and. ierr /= -1) &
          call log_warn('Error in checkpoint_params namelist, using defaults')
      end if
      close (unit)
    else if (present(nml_string)) then
      read (nml_string, nml=checkpoint_params)
    else
      call abort_run('Reading checkpoint config failed! Provide at least one &
                     &of the following: file name or source', &
                     collective=.true.)
    end if

    self%checkpoint_freq = checkpoint_freq
    self%snapshot_freq = snapshot_freq
    self%keep_checkpoint = keep_checkpoint
    self%checkpoint_prefix = checkpoint_prefix
    self%snapshot_prefix = snapshot_prefix
    self%restart_from_checkpoint = restart_from_checkpoint
    self%restart_file = restart_file
    self%output_stride = output_stride
    self%snapshot_sp = snapshot_sp
    self%output_fields = output_fields
  end subroutine read_checkpoint_nml

  subroutine read_inflow_nml(self, nml_file, nml_string)
    !! The namelist is optional: a missing inflow_nml keeps the defaults.
    implicit none

    class(inflow_config_t) :: self
    character(*), optional, intent(in) :: nml_file
    character(*), optional, intent(in) :: nml_string

    integer :: unit, ierr

    logical :: write_planes, read_planes
    integer :: plane_x_index, ntimesteps, sample_every, start_iter, ninflows
    character(len=256) :: inflow_prefix, inflow_read_prefix

    namelist /inflow_nml/ write_planes, plane_x_index, ntimesteps, &
      sample_every, start_iter, inflow_prefix, read_planes, ninflows, &
      inflow_read_prefix

    write_planes = .false.
    plane_x_index = 0
    ntimesteps = 100
    sample_every = 1
    start_iter = 0
    inflow_prefix = 'inflow'
    read_planes = .false.
    ninflows = 1
    inflow_read_prefix = 'inflow'

    if (present(nml_file) .and. present(nml_string)) then
      call abort_run('Reading inflow config failed! Provide only a file name &
                     &or source, not both.', collective=.true.)
    else if (present(nml_file)) then
      open (newunit=unit, file=nml_file, iostat=ierr)
      if (ierr == 0) then
        read (unit, nml=inflow_nml, iostat=ierr)
        ! -1 is end of file: the namelist is absent
        if (ierr /= 0 .and. ierr /= -1) &
          call abort_run('Reading inflow config failed! Error in &
                         &inflow_nml.', collective=.true.)
      end if
      close (unit)
    else if (present(nml_string)) then
      read (nml_string, nml=inflow_nml)
    else
      call abort_run('Reading inflow config failed! Provide at least one of &
                     &the following: file name or source', collective=.true.)
    end if

    self%write_planes = write_planes
    self%plane_x_index = plane_x_index
    self%ntimesteps = ntimesteps
    self%sample_every = sample_every
    self%start_iter = start_iter
    self%inflow_prefix = trim(inflow_prefix)
    self%read_planes = read_planes
    self%ninflows = ninflows
    self%inflow_read_prefix = trim(inflow_read_prefix)

    if (self%ntimesteps < 1) &
      call abort_run('inflow_nml: ntimesteps must be at least 1.', &
                     collective=.true.)
    if (self%sample_every < 1) &
      call abort_run('inflow_nml: sample_every must be at least 1.', &
                     collective=.true.)
    if (self%start_iter < 0) &
      call abort_run('inflow_nml: start_iter must be non-negative.', &
                     collective=.true.)
    if (self%plane_x_index < 0) &
      call abort_run('inflow_nml: plane_x_index must be non-negative.', &
                     collective=.true.)
    if (self%ninflows < 1) &
      call abort_run('inflow_nml: ninflows must be at least 1.', &
                     collective=.true.)
    if (self%write_planes .and. self%read_planes .and. &
        self%inflow_prefix == self%inflow_read_prefix) &
      call abort_run('inflow_nml: inflow_prefix and inflow_read_prefix must &
                     &differ, or the run overwrites the planes it reads.', &
                     collective=.true.)
  end subroutine read_inflow_nml

  subroutine read_stats_nml(self, nml_file, nml_string)
    implicit none

    class(stats_config_t) :: self
    character(*), optional, intent(in) :: nml_file
    character(*), optional, intent(in) :: nml_string

    integer :: unit, ierr

    integer :: initstat = 0
    integer :: istatfreq = 1
    integer :: istatout = 0
    character(len=256) :: stats_prefix = "statistics"

    namelist /stats_params/ initstat, istatfreq, istatout, stats_prefix

    if (present(nml_file) .and. present(nml_string)) then
      call abort_run('Reading stats config failed! Provide only a file name &
                     &or source, not both.', collective=.true.)
    else if (present(nml_file)) then
      open (newunit=unit, file=nml_file, iostat=ierr)
      if (ierr == 0) then
        read (unit, nml=stats_params, iostat=ierr)
        if (ierr /= 0 .and. ierr /= -1) &
          call log_warn('Error in stats_params namelist, using defaults')
      end if
      close (unit)
    else if (present(nml_string)) then
      read (nml_string, nml=stats_params)
    else
      call abort_run('Reading stats config failed! Provide at least one of &
                     &the following: file name or source', collective=.true.)
    end if

    self%initstat = initstat
    self%istatfreq = istatfreq
    self%istatout = istatout
    self%stats_prefix = stats_prefix
  end subroutine read_stats_nml

  pure logical function has_output_field(config, name)
    !! Check whether a field name is present in the output_fields list.
    type(checkpoint_config_t), intent(in) :: config
    character(*), intent(in) :: name

    has_output_field = any(config%output_fields == name)
  end function has_output_field

  pure logical function output_field_active(checkpoint_cfg, name) &
    result(active)
    !! Whether a snapshot output field is actually driven at runtime - a
    !! field only triggers its associated allocations if BOTH it is listed
    !! in output_fields AND snapshotting is enabled. Mirrors the pattern
    !! used for keep_pressure (src/case/base_case.f90:120-122), shared by
    !! run() and x3d2-memcheck.
    type(checkpoint_config_t), intent(in) :: checkpoint_cfg
    character(*), intent(in) :: name

    active = has_output_field(checkpoint_cfg, name) .and. &
             checkpoint_cfg%snapshot_freq > 0
  end function output_field_active

  logical function nml_group_in_file(nml_file, group) result(found)
    !! Whether a line of the file opens the namelist group (&group, in any
    !! case), so an optional namelist that is present but unterminated can be
    !! told apart from an absent one.
    character(*), intent(in) :: nml_file
    character(*), intent(in) :: group

    character(len=1024) :: line
    integer :: unit, ios, i, c, n

    found = .false.
    n = len(group) + 1
    open (newunit=unit, file=nml_file, status='old', action='read', &
          iostat=ios)
    if (ios /= 0) return
    do
      read (unit, '(a)', iostat=ios) line
      if (ios /= 0) exit
      line = adjustl(line)
      do i = 1, n
        c = iachar(line(i:i))
        if (c >= iachar('A') .and. c <= iachar('Z')) &
          line(i:i) = achar(c + 32)
      end do
      ! A blank after the name, so a longer group name does not match
      found = line(1:n) == '&'//group .and. line(n + 1:n + 1) == ' '
      if (found) exit
    end do
    close (unit)
  end function nml_group_in_file

end module m_config
