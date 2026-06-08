module m_turbine_model
  !! Abstract interface for turbine forcing models (actuator disc, actuator
  !! line)
  !!
  !! A turbine model is a backend-agnostic physics object: it samples the flow
  !! at actuator points and accumulates a momentum source back onto the RHS
  !! fields (du, dv, dw). It is owned and driven entirely by the flow case
  !! (see m_case_wind_turbine) and therefore knows nothing about the solver,
  !! the case, or the `iturbine` switch - the case decides which concrete model
  !! to allocate and when to call it.
  !!
  !! Concrete implementations:
  !!   - turbine_dummy_t (m_turbine_dummy): no-op / iturbine=0, the fallback
  !!   - alm_t           (m_alm)          : actuator line model
  !!   - adm_t           (m_adm)          : actuator disc model
  use m_allocator, only: allocator_t
  use m_base_backend, only: base_backend_t
  use m_common, only: dp
  use m_field, only: field_t
  use m_mesh, only: mesh_t

  implicit none

  type, abstract :: turbine_model_t
  contains
    procedure(init_iface), deferred :: init
    procedure(update_iface), deferred :: update
    procedure(project_forces_iface), deferred :: project_forces
    procedure(write_output_iface), deferred :: write_output
  end type turbine_model_t

  abstract interface
    subroutine init_iface(self, backend, mesh, host_allocator, dt)
      !! Set up turbine geometry/state (read config files, place actuator
      !! points, build smearing kernels, ...). Called once after the solver
      !! exists. The model caches whatever context it needs (backend, mesh,
      !! host_allocator) so the per-step hooks below take only fields.
      import :: turbine_model_t, base_backend_t, mesh_t, allocator_t, dp
      implicit none
      class(turbine_model_t), intent(inout) :: self
      class(base_backend_t), target, intent(in) :: backend
      type(mesh_t), target, intent(in) :: mesh
      type(allocator_t), target, intent(in) :: host_allocator
      real(dp), intent(in) :: dt
    end subroutine init_iface

    subroutine update_iface(self, t, dt)
      !! Advance the turbine internal state (blade kinematics, rotor speed,
      !! velocity filtering) by one step. Called once per substep before
      !! project_forces.
      import :: turbine_model_t, dp
      implicit none
      class(turbine_model_t), intent(inout) :: self
      real(dp), intent(in) :: t, dt
    end subroutine update_iface

    subroutine project_forces_iface(self, du, dv, dw, u, v, w)
      !! Sample the velocity field at the actuator points and accumulate the
      !! turbine momentum source into the RHS fields du, dv, dw. Uses the
      !! backend/mesh cached at init.
      import :: turbine_model_t, field_t
      implicit none
      class(turbine_model_t), intent(inout) :: self
      class(field_t), intent(inout) :: du, dv, dw
      class(field_t), intent(in) :: u, v, w
    end subroutine project_forces_iface

    subroutine write_output_iface(self, iter, is_root)
      !! Write per-turbine diagnostics (power, thrust, loads).
      import :: turbine_model_t
      implicit none
      class(turbine_model_t), intent(inout) :: self
      integer, intent(in) :: iter
      logical, intent(in) :: is_root
    end subroutine write_output_iface
  end interface

end module m_turbine_model
