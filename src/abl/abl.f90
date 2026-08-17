module m_abl
!! Reusable, backend-agnostic Atmospheric Boundary Layer driver.
!!
!! Ported from Incompact3d's `Case-ABL.f90` (itype=10). The driver owns the
!! ABL physics and is owned by a case; it never uses `m_solver`/`m_base_case`,
!! taking fields/mesh/coords in (same pattern as `m_ibm`).
!!
!! This skeleton provides the type, `init`, and the three physics entry points.
!! `initialise`/`apply_forcing` are filled in issue #317 commit 2;
!! `apply_wall_model` is blocked on the LES `nut` field (issue #321) and lands
!! in commit 3.
  use m_allocator, only: allocator_t
  use m_base_backend, only: base_backend_t
  use m_common, only: dp, VERT
  use m_config, only: abl_config_t
  use m_field, only: field_t
  use m_mesh, only: mesh_t

  implicit none

  private
  public :: abl_t

  type :: abl_t
    class(base_backend_t), pointer :: backend => null()
    class(mesh_t), pointer :: mesh => null()
    type(allocator_t), pointer :: host_allocator => null()
    type(abl_config_t) :: cfg
    real(dp) :: dt = 0._dp
  contains
    procedure :: initialise
    procedure :: apply_forcing
    procedure :: apply_wall_model
  end type abl_t

  interface abl_t
    module procedure init
  end interface abl_t

contains

  function init(backend, mesh, host_allocator, cfg, dt) result(abl)
    implicit none

    class(base_backend_t), target, intent(inout) :: backend
    type(mesh_t), target, intent(inout) :: mesh
    type(allocator_t), target, intent(inout) :: host_allocator
    type(abl_config_t), intent(in) :: cfg
    real(dp), intent(in) :: dt
    type(abl_t) :: abl

    abl%backend => backend
    abl%mesh => mesh
    abl%host_allocator => host_allocator
    abl%cfg = cfg
    abl%dt = dt

  end function init

  subroutine initialise(self, u, v, w)
    !! Sets the t=0 velocity profile.
    implicit none

    class(abl_t) :: self
    class(field_t), intent(inout) :: u, v, w

    ! Placeholder zero field; the log-law/geostrophic profile lands in
    ! issue #317 commit 2.
    call u%fill(0._dp)
    call v%fill(0._dp)
    call w%fill(0._dp)

    call u%set_data_loc(VERT)
    call v%set_data_loc(VERT)
    call w%set_data_loc(VERT)

  end subroutine initialise

  subroutine apply_forcing(self, du, dv, dw, u, v, w)
    !! Adds the ABL driving forces to the momentum RHS each substep.
    implicit none

    class(abl_t) :: self
    class(field_t), intent(inout) :: du, dv, dw
    class(field_t), intent(in) :: u, v, w

    ! Pressure-gradient / Coriolis / damping forcing lands in issue #317
    ! commit 2.

  end subroutine apply_forcing

  subroutine apply_wall_model(self, u, v, w, nut)
    !! Imposes the neutral log-law wall stress at the bottom face, blended
    !! with the resolved SGS stress `nut`.
    implicit none

    class(abl_t) :: self
    class(field_t), intent(inout) :: u, v, w
    class(field_t), intent(in) :: nut

    ! Wall model is blocked on the LES `nut` field (issue #321) and lands in
    ! issue #317 commit 3.

  end subroutine apply_wall_model

end module m_abl
