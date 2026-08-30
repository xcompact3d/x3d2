module m_config
  !! Contains all the namelist configurations used in x3d2
  use iso_fortran_env, only: stderr => error_unit

  use m_common

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
  contains
    procedure :: read => read_solver_nml
  end type solver_config_t

  type, extends(base_config_t) :: les_config_t
    !! Configuration for explicit sub-grid-scale modelling.
    character(len=20) :: model = 'none'
    real(dp) :: smagorinsky_constant = 0.14_dp
    logical :: wall_damping = .false.
    real(dp) :: wall_damping_n = 3._dp
    real(dp) :: von_karman_constant = 0.4_dp
    real(dp) :: roughness_length = 0._dp
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

  type, extends(base_config_t) :: abl_config_t
    real(dp) :: z0 = 0.1_dp            !! aerodynamic roughness length
    real(dp) :: u_star = 0._dp         !! friction velocity
    real(dp) :: delta = 1._dp          !! boundary-layer depth
    real(dp) :: kappa = 0.41_dp        !! von Karman constant
    real(dp) :: dsampling = 3._dp      !! wall-model sampling height, in dy
    real(dp) :: u_geo(3) = 0._dp       !! geostrophic wind UG
    real(dp) :: coriolis_freq = 0._dp  !! Coriolis frequency f
    real(dp) :: init_noise(3) = 0._dp
    real(dp) :: u_bulk = 0._dp         !! target bulk velocity (mass_conserve)
    integer :: profile_start_iter = -1 !! first profile sample (-1 = disabled)
    character(len=256) :: profile_file = 'abl_profile.csv'
    logical :: pressure_gradient = .false.
    logical :: coriolis = .false.
    logical :: mass_conserve = .false.
    logical :: damping = .false.
  contains
    procedure :: read => read_abl_nml
  end type abl_config_t

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
      error stop 'Reading domain config failed! &
                 &Provide only a file name or source, not both.'
    else if (present(nml_file)) then
      open (newunit=unit, file=nml_file)
      read (unit, nml=domain_settings)
      close (unit)
    else if (present(nml_string)) then
      read (nml_string, nml=domain_settings)
    else
      error stop 'Reading domain config failed! &
                 &Provide at least one of the following: file name or source'
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
    !! Incompact3d uses 0.49, but it solves the filter with a full Thomas
    !! algorithm. x3d2's distributed solver truncates, and the filter system
    !! is only marginally diagonally dominant as alpha -> 0.5, so 0.49 leaves
    !! a 1e-3 error at zero wavenumber on 32 points per rank. 0.4 removes the
    !! 2*dx mode just as exactly while staying well within the solver.
    real(dp) :: filter_alpha = 0.4_dp
    character(3) :: time_intg
    character(3) :: poisson_solver_type = 'FFT'
    character(30) :: der1st_scheme = 'compact6', der2nd_scheme = 'compact6', &
                     interpl_scheme = 'classic', stagder_scheme = 'compact6'

    namelist /solver_params/ Re, dt, n_iters, n_output, poisson_solver_type, &
      spatial_filter, filter_alpha, &
      n_species, pr_species, lowmem_transeq, lowmem_fft, &
      time_intg, der1st_scheme, der2nd_scheme, interpl_scheme, &
      stagder_scheme, ibm_on

    if (present(nml_file) .and. present(nml_string)) then
      error stop 'Reading solver config failed! &
                 &Provide only a file name or source, not both.'
    else if (present(nml_file)) then
      open (newunit=unit, file=nml_file)
      read (unit, nml=solver_params)
      close (unit)
    else if (present(nml_string)) then
      read (nml_string, nml=solver_params)
    else
      error stop 'Reading solver config failed! &
                 &Provide at least one of the following: file name or source'
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

  end subroutine read_solver_nml

  subroutine read_les_nml(self, nml_file, nml_string)
    class(les_config_t) :: self
    character(*), optional, intent(in) :: nml_file
    character(*), optional, intent(in) :: nml_string

    integer :: unit, ierr
    character(len=20) :: model
    real(dp) :: smagorinsky_constant, wall_damping_n, von_karman_constant
    real(dp) :: roughness_length
    logical :: wall_damping

    namelist /les_params/ model, smagorinsky_constant, wall_damping, &
      wall_damping_n, von_karman_constant, roughness_length

    model = self%model
    smagorinsky_constant = self%smagorinsky_constant
    wall_damping = self%wall_damping
    wall_damping_n = self%wall_damping_n
    von_karman_constant = self%von_karman_constant
    roughness_length = self%roughness_length

    if (present(nml_file) .and. present(nml_string)) then
      error stop 'Reading LES config failed! &
                 & Provide only a file name or source.'
    else if (present(nml_file)) then
      open (newunit=unit, file=nml_file, iostat=ierr)
      if (ierr /= 0) error stop 'Opening LES config file failed.'
      read (unit, nml=les_params, iostat=ierr)
      close (unit)
      ! Existing input files may omit this optional namelist. End-of-file
      ! therefore selects the defaults, while a malformed block remains fatal.
      if (ierr > 0) error stop 'Reading LES config failed.'
    else if (present(nml_string)) then
      read (nml_string, nml=les_params)
    else
      error stop 'Reading LES config failed! Provide a file name or source.'
    end if

    select case (trim(model))
    case ('none', 'smagorinsky')
    case default
      error stop 'Unknown LES model. Use "none" or "smagorinsky".'
    end select
    if (smagorinsky_constant <= 0._dp) &
      error stop 'smagorinsky_constant must be positive.'
    if (von_karman_constant <= 0._dp .or. wall_damping_n <= 0._dp) &
      error stop 'LES wall-damping constants must be positive.'
    if (roughness_length < 0._dp) &
      error stop 'roughness_length must not be negative.'

    self%model = trim(model)
    self%smagorinsky_constant = smagorinsky_constant
    self%wall_damping = wall_damping
    self%wall_damping_n = wall_damping_n
    self%von_karman_constant = von_karman_constant
    self%roughness_length = roughness_length
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
      error stop 'Reading channel config failed! &
                 &Provide only a file name or source, not both.'
    else if (present(nml_file)) then
      open (newunit=unit, file=nml_file)
      read (unit, nml=channel_nml)
      close (unit)
    else if (present(nml_string)) then
      read (nml_string, nml=channel_nml)
    else
      error stop 'Reading channel config failed! &
                 &Provide at least one of the following: file name or source'
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
      error stop 'Reading cylinder config failed! &
                 &Provide only a file name or source, not both.'
    else if (present(nml_file)) then
      open (newunit=unit, file=nml_file)
      read (unit, nml=cylinder_nml)
      close (unit)
    else if (present(nml_string)) then
      read (nml_string, nml=cylinder_nml)
    else
      error stop 'Reading cylinder config failed! &
                 &Provide at least one of the following: file name or source'
    end if

    self%init_noise = init_noise
    self%inlet_noise = inlet_noise

  end subroutine read_cylinder_nml

  subroutine read_abl_nml(self, nml_file, nml_string)
    implicit none

    class(abl_config_t) :: self
    character(*), optional, intent(in) :: nml_file
    character(*), optional, intent(in) :: nml_string

    integer :: unit

    real(dp) :: z0, u_star, delta, kappa, dsampling
    real(dp) :: u_geo(3), coriolis_freq, init_noise(3), u_bulk
    integer :: profile_start_iter
    character(len=256) :: profile_file
    logical :: pressure_gradient, coriolis, mass_conserve, damping

    namelist /abl_nml/ z0, u_star, delta, kappa, dsampling, &
      u_geo, coriolis_freq, &
      init_noise, u_bulk, profile_start_iter, profile_file, &
      pressure_gradient, coriolis, mass_conserve, damping

    ! Defaults
    z0 = 0.1_dp
    u_star = 0._dp
    delta = 1._dp
    kappa = 0.41_dp
    dsampling = 3._dp
    u_geo = 0._dp
    coriolis_freq = 0._dp
    init_noise = 0._dp
    u_bulk = 0._dp
    profile_start_iter = -1
    profile_file = 'abl_profile.csv'
    pressure_gradient = .false.
    coriolis = .false.
    mass_conserve = .false.
    damping = .false.

    if (present(nml_file) .and. present(nml_string)) then
      error stop 'Reading ABL config failed! &
                 &Provide only a file name or source, not both.'
    else if (present(nml_file)) then
      open (newunit=unit, file=nml_file)
      read (unit, nml=abl_nml)
      close (unit)
    else if (present(nml_string)) then
      read (nml_string, nml=abl_nml)
    else
      error stop 'Reading ABL config failed! &
                 &Provide at least one of the following: file name or source'
    end if

    self%z0 = z0
    self%u_star = u_star
    self%delta = delta
    self%kappa = kappa
    self%dsampling = dsampling
    self%u_geo = u_geo
    self%coriolis_freq = coriolis_freq
    self%init_noise = init_noise
    self%u_bulk = u_bulk
    self%profile_start_iter = profile_start_iter
    self%profile_file = trim(profile_file)
    self%pressure_gradient = pressure_gradient
    self%coriolis = coriolis
    self%mass_conserve = mass_conserve
    self%damping = damping

    call validate_abl_config(self)

  end subroutine read_abl_nml

  subroutine validate_abl_config(config)
    !! Validate relationships between ABL namelist values in one place.
    type(abl_config_t), intent(in) :: config

    ! Damping is a sponge rather than a driving mechanism.
    if (.not. (config%pressure_gradient .or. config%coriolis .or. &
               config%mass_conserve)) &
      error stop 'ABL config error: enable pressure_gradient, coriolis, &
                 &or mass_conserve.'

    if (config%z0 <= 0._dp .or. config%delta <= config%z0 .or. &
        config%kappa <= 0._dp) &
      error stop 'ABL config error: require delta > z0 > 0 and kappa > 0.'

    if (any(config%init_noise < 0._dp)) &
      error stop 'ABL config error: init_noise must not be negative.'

    if (config%profile_start_iter < -1) &
      error stop 'ABL config error: profile_start_iter must be -1 or greater.'
    if (config%profile_start_iter >= 0 .and. &
        len_trim(config%profile_file) == 0) &
      error stop 'ABL config error: profile_file must not be empty.'
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
      error stop 'Reading checkpoint config failed! &
                 &Provide only a file name or source, not both.'
    else if (present(nml_file)) then
      open (newunit=unit, file=nml_file, iostat=ierr)
      if (ierr == 0) then
        read (unit, nml=checkpoint_params, iostat=ierr)

        if (ierr /= 0 .and. ierr /= -1) &
          print *, 'WARNING: Error in checkpoint_params namelist, &
          & using defaults'
      end if
      close (unit)
    else if (present(nml_string)) then
      read (nml_string, nml=checkpoint_params)
    else
      error stop 'Reading checkpoint config failed! &
                 &Provide at least one of the following: file name or source'
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
      error stop 'Reading stats config failed! &
                 &Provide only a file name or source, not both.'
    else if (present(nml_file)) then
      open (newunit=unit, file=nml_file, iostat=ierr)
      if (ierr == 0) then
        read (unit, nml=stats_params, iostat=ierr)
        if (ierr /= 0 .and. ierr /= -1) &
          print *, 'WARNING: Error in stats_params namelist, &
          & using defaults'
      end if
      close (unit)
    else if (present(nml_string)) then
      read (nml_string, nml=stats_params)
    else
      error stop 'Reading stats config failed! &
                 &Provide at least one of the following: file name or source'
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

end module m_config
