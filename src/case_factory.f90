module m_case_factory
  !! Flow case dispatch shared by xcompact and x3d2-memcheck.
  use m_allocator, only: allocator_t
  use m_base_backend, only: base_backend_t
  use m_base_case, only: base_case_t
  use m_logger, only: abort_run
  use m_mesh, only: mesh_t
  use m_case_abl, only: case_abl_t
  use m_case_channel, only: case_channel_t
  use m_case_cylinder, only: case_cylinder_t
  use m_case_generic, only: case_generic_t
  use m_case_tgv, only: case_tgv_t
  use m_case_wind_turbine, only: case_wind_turbine_t

  implicit none

  private
  public :: allocate_flow_case, is_supported_flow_case

contains

  subroutine allocate_flow_case(name, backend, mesh, host_allocator, flow_case)
    !! Build the flow case with the given name.
    character(len=*), intent(in) :: name
    class(base_backend_t), pointer, intent(in) :: backend
    type(mesh_t), target, intent(inout) :: mesh
    type(allocator_t), target, intent(inout) :: host_allocator
    class(base_case_t), allocatable, intent(out) :: flow_case

    select case (trim(name))
    case ('abl')
      allocate (case_abl_t :: flow_case)
      flow_case = case_abl_t(backend, mesh, host_allocator)
    case ('channel')
      allocate (case_channel_t :: flow_case)
      flow_case = case_channel_t(backend, mesh, host_allocator)
    case ('cylinder')
      allocate (case_cylinder_t :: flow_case)
      flow_case = case_cylinder_t(backend, mesh, host_allocator)
    case ('generic')
      allocate (case_generic_t :: flow_case)
      flow_case = case_generic_t(backend, mesh, host_allocator)
    case ('tgv')
      allocate (case_tgv_t :: flow_case)
      flow_case = case_tgv_t(backend, mesh, host_allocator)
    case ('wind_turbine')
      allocate (case_wind_turbine_t :: flow_case)
      flow_case = case_wind_turbine_t(backend, mesh, host_allocator)
    case default
      call abort_run('Undefined flow_case: '//trim(name), collective=.true.)
    end select
  end subroutine allocate_flow_case

  pure logical function is_supported_flow_case(name)
    !! True when allocate_flow_case has a branch for the name.
    character(len=*), intent(in) :: name

    select case (trim(name))
    case ('abl', 'channel', 'cylinder', 'generic', 'tgv', 'wind_turbine')
      is_supported_flow_case = .true.
    case default
      is_supported_flow_case = .false.
    end select
  end function is_supported_flow_case

end module m_case_factory
