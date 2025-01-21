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

  abstract interface
    subroutine read(self, file_name, src) !&
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
      open (100, file=file_name)
      read (100, nml=domain_settings)
      close (100)
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

end module m_config

