module m_config
  !! Contains all the namelist configurations used in x3d2
  use iso_fortran_env, only: stderr => error_unit

  use m_common

  implicit none

  integer, parameter :: n_species_max = 99

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
    character(3) :: poisson_solver_type, time_intg
    character(30) :: der1st_scheme, der2nd_scheme, &
                     interpl_scheme, stagder_scheme
  contains
    procedure :: read => read_solver_nml
  end type solver_config_t

  type, extends(base_config_t) :: channel_config_t
    real(dp) :: noise, omega_rot
    logical :: rotation
    integer :: n_rotate
  contains
    procedure :: read => read_channel_nml
  end type channel_config_t

  type, extends(base_config_t) :: checkpoint_config_t
    integer :: checkpoint_freq = 0                         !! Frequency of checkpointing (0 = off)
    integer :: snapshot_freq = 0                           !! Frequency of snapshots (0 = off)
    logical :: keep_checkpoint = .true.                    !! If false, only keep latest checkpoint
    character(len=256) :: checkpoint_prefix = "checkpoint"
    character(len=256) :: snapshot_prefix = "snapshot"
    logical :: restart_from_checkpoint = .false.
    character(len=256) :: restart_file = ""
    integer, dimension(3) :: output_stride = [2, 2, 2]     !! Spatial stride for snapshot output
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
    character(3) :: time_intg
    character(3) :: poisson_solver_type = 'FFT'
    character(30) :: der1st_scheme = 'compact6', der2nd_scheme = 'compact6', &
                     interpl_scheme = 'classic', stagder_scheme = 'compact6'

    namelist /solver_params/ Re, dt, n_iters, n_output, poisson_solver_type, &
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

  subroutine read_channel_nml(self, nml_file, nml_string)
    implicit none

    class(channel_config_t) :: self
    character(*), optional, intent(in) :: nml_file
    character(*), optional, intent(in) :: nml_string

    integer :: unit

    real(dp) :: noise, omega_rot
    logical :: rotation
    integer :: n_rotate

    namelist /channel_nml/ noise, rotation, omega_rot, n_rotate

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

    self%noise = noise
    self%rotation = rotation
    self%omega_rot = omega_rot
    self%n_rotate = n_rotate

  end subroutine read_channel_nml

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

    namelist /checkpoint_params/ checkpoint_freq, snapshot_freq, &
      keep_checkpoint, checkpoint_prefix, snapshot_prefix, &
      restart_from_checkpoint, restart_file, output_stride

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
  end subroutine read_checkpoint_nml

end module m_config

