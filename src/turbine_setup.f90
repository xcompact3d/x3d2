module m_turbine_setup
  !! Turbine set-up and tear-down shared by the cases that carry turbines
  !! (wind_turbine and abl). This is the only place the `iturbine` switch
  !! appears (at allocation time); downstream a case calls the polymorphic
  !! turbine object with no branching.
  use m_adm, only: adm_t
  use m_alm, only: alm_t
  use m_config, only: wind_turbine_config_t
  use m_io_manager, only: io_manager_t
  use m_logger, only: abort_run
  use m_solver, only: solver_t
  use m_turbine_dummy, only: turbine_dummy_t
  use m_turbine_model, only: turbine_model_t

  implicit none

  private
  public :: init_turbine_model, finalise_turbine_model

contains

  subroutine init_turbine_model(turbine, cfg, solver, io_mgr)
    !! Allocate the turbine model wind_turbine_nml selects, initialise it and
    !! register it for checkpointing
    implicit none

    class(turbine_model_t), pointer, intent(inout) :: turbine
    type(wind_turbine_config_t), intent(in) :: cfg
    class(solver_t), intent(in) :: solver
    type(io_manager_t), intent(inout) :: io_mgr

    select case (cfg%iturbine)
    case (1)
      allocate (alm_t :: turbine)
      select type (turbine)
      type is (alm_t)
        turbine%turbine_files = cfg%alm_turbines
        turbine%rho_air = cfg%rho_air
        turbine%eps_factor = cfg%eps_factor
        turbine%nu = solver%nu
        turbine%u_inflow = cfg%bc_start_u
        turbine%stats_start = max(1, io_mgr%stats_mgr%config%initstat)
        turbine%stats_freq = max(1, io_mgr%stats_mgr%config%istatfreq)
      end select
    case (2)
      allocate (adm_t :: turbine)
      select type (turbine)
      type is (adm_t)
        turbine%coords_file = cfg%adm_coords
        turbine%rho_air = cfg%rho_air
        turbine%T_relax = cfg%T_relax
        turbine%stats_start = max(1, io_mgr%stats_mgr%config%initstat)
        turbine%stats_freq = max(1, io_mgr%stats_mgr%config%istatfreq)
      end select
    case (0)
      allocate (turbine_dummy_t :: turbine)
    case default
      call abort_run('Unsupported wind turbine model selected.', &
                     collective=.true.)
    end select

    call turbine%init(solver%backend, solver%mesh, solver%host_allocator, &
                      solver%dt)
    call io_mgr%register_checkpoint_state(turbine)
    call turbine%setup_output(solver%mesh%par%is_root(), io_mgr%is_restart())
  end subroutine init_turbine_model

  subroutine finalise_turbine_model(turbine, io_mgr)
    !! Release a turbine model set up by init_turbine_model; a no-op when
    !! the case has none
    implicit none

    class(turbine_model_t), pointer, intent(inout) :: turbine
    type(io_manager_t), intent(inout) :: io_mgr

    if (.not. associated(turbine)) return
    call turbine%finalise()
    call io_mgr%unregister_checkpoint_state()
    deallocate (turbine)
  end subroutine finalise_turbine_model

end module m_turbine_setup
