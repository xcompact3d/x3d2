module m_abl
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
  use m_mesh, only: mesh_t
  use m_les, only: les_t

  implicit none

  private
  public :: abl_t

  type :: abl_t
    class(base_backend_t), pointer :: backend => null()
    class(mesh_t), pointer :: mesh => null()
    type(allocator_t), pointer :: host_allocator => null()
    type(abl_config_t) :: cfg
    real(dp) :: dt = 0._dp
    ! Cached Rayleigh damping coefficient field (coeff*lambda(y)), built once.
    class(field_t), pointer :: damp_coeff => null()
  contains
    procedure :: initialise
    procedure :: apply_forcing
    procedure :: configure_wall_boundary_correction
    procedure :: apply_damping
    procedure :: build_damping
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
    !! Sets the t=0 velocity profile: log-law when the drive is pressure-
    !! gradient / mass-conserve, otherwise the geostrophic wind, plus noise
    !! (Incompact3d init_abl, 68-92).
    implicit none

    class(abl_t) :: self
    class(field_t), intent(inout) :: u, v, w

    class(field_t), pointer :: hu, hv, hw
    integer :: i, j, k, dims(3)
    real(dp) :: coords(3), y, prof, noise(3)
    logical :: log_law

    dims = self%mesh%get_dims(VERT)

    hu => self%host_allocator%get_block(DIR_C)
    hv => self%host_allocator%get_block(DIR_C)
    hw => self%host_allocator%get_block(DIR_C)

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
    ! the pressure gradient is off.
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
    !! Precomputes the damping coefficient coeff*lambda(y), constant in x/z.
    implicit none

    class(abl_t) :: self

    class(field_t), pointer :: h
    integer :: i, j, k, dims(3)
    real(dp) :: coords(3), y, coeff, dheight, ylo, yhi, lambda
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
    real(dp) :: dy, sampling_height
    integer :: dims(3), sample_plane

    dims = self%mesh%get_dims(VERT)
    if (dims(2) < 3) &
      error stop 'ABL wall model requires at least 3 y vertices.'

    dy = abs(self%mesh%geo%vert_coords(2, 2) &
             - self%mesh%geo%vert_coords(1, 2))
    ! Vertex 1 is the wall, so a sampling height of n*dy is vertex n+1.
    sample_plane = nint(self%cfg%dsampling) + 1
    if (sample_plane < 2 .or. sample_plane > dims(2)) &
      error stop 'ABL dsampling puts the wall-model sample outside the domain.'
    sampling_height = self%mesh%geo%vert_coords(sample_plane, 2) &
                      - self%mesh%geo%vert_coords(1, 2)
    if (sampling_height <= self%cfg%z0) &
      error stop 'ABL wall-model sampling height must exceed z0.'

    call les%configure_abl_wall_boundary( &
      self%cfg%kappa, self%cfg%z0, sampling_height, sample_plane)
  end subroutine configure_wall_boundary_correction

end module m_abl
