module m_config
  !! Contains all the namelist configurations used in x3d2
  use iso_fortran_env, only: stderr => error_unit

  use m_common

  implicit none

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
  contains
    procedure :: read => read_domain_nml
  end type domain_config_t

  type, extends(base_config_t) :: solver_config_t
    real(dp) :: Re, dt
    integer :: n_iters, n_output
    character(3) :: poisson_solver_type, time_intg
    character(30) :: der1st_scheme, der2nd_scheme, &
                     interpl_scheme, stagder_scheme
  contains
    procedure :: read => read_solver_nml
  end type solver_config_t

  abstract interface
    subroutine read(self, file_name, src) !&
      !! Assigns the member variables either from a file or text source.
      !!
      !! file_name can be an absolute or relative path
      !! src is a character string that contains the namelist definition.
      !! For example, src="&foobar_nml foo=0, bar='this'/"
      import :: base_config_t

      class(base_config_t) :: self
      character(*), optional, intent(in) :: file_name
      character(*), optional, intent(in) :: src
    end subroutine read
  end interface

contains

  subroutine read_domain_nml(self, file_name, src)
    implicit none

    class(domain_config_t) :: self
    character(*), optional, intent(in) :: file_name
    character(*), optional, intent(in) :: src

    integer :: unit

    character(len=20) :: flow_case_name
    real(dp), dimension(3) :: L_global
    integer, dimension(3) :: dims_global
    integer, dimension(3) :: nproc_dir
    character(len=20) :: BC_x(2), BC_y(2), BC_z(2)

    namelist /domain_settings/ flow_case_name, L_global, dims_global, &
      nproc_dir, BC_x, BC_y, BC_z

    if (present(file_name) .and. present(src)) then
      error stop 'Reading domain config failed! &
                 &Provide only a file name or source, not both.'
    else if (present(file_name)) then
      open (newunit=unit, file=file_name)
      read (unit, nml=domain_settings)
      close (unit)
    else if (present(src)) then
      read (src, nml=domain_settings)
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

  end subroutine read_domain_nml

  subroutine read_solver_nml(self, file_name, src)
    implicit none

    class(solver_config_t) :: self
    character(*), optional, intent(in) :: file_name
    character(*), optional, intent(in) :: src

    integer :: unit

    real(dp) :: Re, dt
    integer :: n_iters, n_output
    character(3) :: poisson_solver_type = 'FFT', time_intg
    character(30) :: der1st_scheme = 'compact6', der2nd_scheme = 'compact6', &
                     interpl_scheme = 'classic', stagder_scheme = 'compact6'

    namelist /solver_params/ Re, dt, n_iters, n_output, poisson_solver_type, &
      time_intg, der1st_scheme, der2nd_scheme, interpl_scheme, stagder_scheme

    if (present(file_name) .and. present(src)) then
      error stop 'Reading solver config failed! &
                 &Provide only a file name or source, not both.'
    else if (present(file_name)) then
      open (newunit=unit, file=file_name)
      read (unit, nml=solver_params)
      close (unit)
    else if (present(src)) then
      read (src, nml=solver_params)
    else
      error stop 'Reading solver config failed! &
                 &Provide at least one of the following: file name or source'
    end if

    self%Re = Re
    self%dt = dt
    self%n_iters = n_iters
    self%n_output = n_output
    self%poisson_solver_type = poisson_solver_type
    self%time_intg = time_intg
    self%der1st_scheme = der1st_scheme
    self%der2nd_scheme = der2nd_scheme
    self%interpl_scheme = interpl_scheme
    self%stagder_scheme = stagder_scheme

  end subroutine read_solver_nml

end module m_config

