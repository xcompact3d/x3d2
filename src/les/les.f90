module m_les
  !! Shared physics for explicit eddy-viscosity LES models.
  !!
  !! Derivatives and data movement are orchestrated here, while pointwise
  !! operations are dispatched to the selected computational backend.
  use m_base_backend, only: base_backend_t
  use m_common, only: dp, DIR_X, DIR_Y, DIR_Z, DIR_C, VERT, &
                      RDR_X2Y, RDR_X2Z, RDR_Y2X, RDR_Z2X
  use m_config, only: les_config_t
  use m_field, only: field_t
  use m_mesh, only: mesh_t
  use m_tdsops, only: dirps_t

  implicit none

  private
  public :: les_t, smagorinsky_nut, strain_rate_magnitude, &
            filter_width, wall_damped_mixing_length

  type :: les_t
    character(len=20) :: model = 'none'
    real(dp) :: smagorinsky_constant = 0.14_dp
    logical :: wall_damping = .false.
    real(dp) :: wall_damping_n = 3._dp
    real(dp) :: von_karman_constant = 0.4_dp
    real(dp) :: roughness_length = 0._dp
    logical :: neutral_wall_model = .false.
    real(dp) :: wall_model_kappa = 0.4_dp
    real(dp) :: wall_model_roughness = 0._dp
    class(field_t), pointer :: nut => null()
    class(field_t), pointer :: mixing_length_sq => null()
  contains
    procedure :: mixing_length
    procedure :: nut_from_gradient
    procedure :: compute_nut
    procedure :: apply_sgs_stress
    procedure :: configure_neutral_wall
    procedure :: finalise
  end type les_t

  interface les_t
    module procedure les_init
  end interface les_t

contains

  pure function les_init(config) result(les)
    type(les_config_t), intent(in) :: config
    type(les_t) :: les

    les%model = config%model
    les%smagorinsky_constant = config%smagorinsky_constant
    les%wall_damping = config%wall_damping
    les%wall_damping_n = config%wall_damping_n
    les%von_karman_constant = config%von_karman_constant
    les%roughness_length = config%roughness_length
  end function les_init

  pure real(dp) function filter_width(spacing) result(delta)
    !! Volume-equivalent grid-filter width, Delta=(dx*dy*dz)^(1/3).
    real(dp), intent(in) :: spacing(3)

    if (any(spacing <= 0._dp)) then
      delta = 0._dp
    else
      delta = product(spacing)**(1._dp/3._dp)
    end if
  end function filter_width

  pure real(dp) function strain_rate_magnitude(velocity_gradient) result(smag)
    !! |S|=sqrt(2 S_ij S_ij), where S_ij=(du_i/dx_j+du_j/dx_i)/2.
    real(dp), intent(in) :: velocity_gradient(3, 3)
    real(dp) :: sij_sq

    sij_sq = velocity_gradient(1, 1)**2 + &
             velocity_gradient(2, 2)**2 + &
             velocity_gradient(3, 3)**2 + &
             0.5_dp*(velocity_gradient(1, 2) + &
                     velocity_gradient(2, 1))**2 + &
             0.5_dp*(velocity_gradient(1, 3) + &
                     velocity_gradient(3, 1))**2 + &
             0.5_dp*(velocity_gradient(2, 3) + &
                     velocity_gradient(3, 2))**2
    smag = sqrt(2._dp*sij_sq)
  end function strain_rate_magnitude

  pure real(dp) function wall_damped_mixing_length( &
    delta, wall_distance, smagorinsky_constant, von_karman_constant, &
    exponent, roughness_length) result(length)
    !! Mason-Thomson blend used by Incompact3d's ABL Smagorinsky path.
    real(dp), intent(in) :: delta, wall_distance, smagorinsky_constant
    real(dp), intent(in) :: von_karman_constant, exponent, roughness_length
    real(dp) :: wall_scale

    if (delta <= 0._dp .or. wall_distance + roughness_length <= 0._dp) then
      length = 0._dp
      return
    end if

    wall_scale = von_karman_constant*(wall_distance + roughness_length)/delta
    length = delta*(smagorinsky_constant**(-exponent) + &
                    wall_scale**(-exponent))**(-1._dp/exponent)
  end function wall_damped_mixing_length

  pure real(dp) function smagorinsky_nut( &
    velocity_gradient, length) result(nut)
    real(dp), intent(in) :: velocity_gradient(3, 3)
    real(dp), intent(in) :: length

    nut = length**2*strain_rate_magnitude(velocity_gradient)
  end function smagorinsky_nut

  pure real(dp) function mixing_length( &
    self, spacing, wall_distance) result(length)
    !! Smagorinsky mixing length, optionally Mason-Thomson wall-damped.
    !!
    !! Returns zero when wall damping is requested without a wall distance,
    !! which in turn zeroes the eddy viscosity built from it.
    class(les_t), intent(in) :: self
    real(dp), intent(in) :: spacing(3)
    real(dp), optional, intent(in) :: wall_distance
    real(dp) :: delta

    delta = filter_width(spacing)
    length = self%smagorinsky_constant*delta
    if (self%wall_damping) then
      if (.not. present(wall_distance)) then
        length = 0._dp
        return
      end if
      length = wall_damped_mixing_length( &
               delta, wall_distance, &
               self%smagorinsky_constant, self%von_karman_constant, &
               self%wall_damping_n, self%roughness_length &
               )
    end if
  end function mixing_length

  pure real(dp) function nut_from_gradient( &
    self, velocity_gradient, spacing, wall_distance) result(nut)
    class(les_t), intent(in) :: self
    real(dp), intent(in) :: velocity_gradient(3, 3), spacing(3)
    real(dp), optional, intent(in) :: wall_distance

    if (trim(self%model) == 'none') then
      nut = 0._dp
      return
    end if

    nut = smagorinsky_nut(velocity_gradient, &
                          self%mixing_length(spacing, wall_distance))
  end function nut_from_gradient

  subroutine configure_neutral_wall(self, kappa, roughness_length)
    class(les_t), intent(inout) :: self
    real(dp), intent(in) :: kappa, roughness_length

    if (trim(self%model) /= 'smagorinsky') &
      error stop 'The neutral ABL wall model requires Smagorinsky LES.'
    if (.not. self%wall_damping) &
      error stop 'The neutral ABL wall model requires LES wall damping.'
    if (kappa <= 0._dp .or. roughness_length <= 0._dp) &
      error stop 'ABL wall-model constants must be positive.'
    if (abs(self%von_karman_constant - kappa) > &
        100._dp*epsilon(1._dp)*max(1._dp, abs(kappa))) &
      error stop 'ABL and LES von Karman constants must match.'
    if (abs(self%roughness_length - roughness_length) > &
        100._dp*epsilon(1._dp)*max(1._dp, abs(roughness_length))) &
      error stop 'ABL and LES roughness lengths must match.'

    self%neutral_wall_model = .true.
    self%wall_model_kappa = kappa
    self%wall_model_roughness = roughness_length
  end subroutine configure_neutral_wall

  subroutine compute_nut(self, backend, mesh, u, v, w, &
                         xdirps, ydirps, zdirps)
    !! Compute the Smagorinsky eddy-viscosity field in DIR_X layout.
    class(les_t), intent(inout) :: self
    class(base_backend_t), intent(inout) :: backend
    type(mesh_t), intent(in) :: mesh
    class(field_t), intent(in) :: u, v, w
    type(dirps_t), intent(in) :: xdirps, ydirps, zdirps

    class(field_t), pointer :: dudx, dudy, dudz
    class(field_t), pointer :: dvdx, dvdy, dvdz
    class(field_t), pointer :: dwdx, dwdy, dwdz
    class(field_t), pointer :: sgs_du, sgs_dv, sgs_dw
    real(dp) :: sampling_height

    if (u%dir /= DIR_X .or. v%dir /= DIR_X .or. w%dir /= DIR_X) then
      error stop 'LES velocity fields must use DIR_X layout.'
    end if

    if (.not. associated(self%nut)) then
      self%nut => backend%allocator%get_block(DIR_X, VERT)
    end if

    select case (trim(self%model))
    case ('none')
      call self%nut%fill(0._dp)
      return
    case ('smagorinsky')
    case default
      error stop 'Unsupported LES model in compute_nut.'
    end select

    if (.not. associated(self%mixing_length_sq)) then
      self%mixing_length_sq => backend%allocator%get_block(DIR_X, VERT)
      call initialise_mixing_length(self, backend, mesh)
    end if

    call compute_velocity_gradients( &
      backend, u, v, w, xdirps, ydirps, zdirps, &
      dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz)

    call backend%compute_smagorinsky_nut( &
      self%nut, self%mixing_length_sq, &
      dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz)

    call release_velocity_gradients( &
      backend, dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz)
  end subroutine compute_nut

  subroutine apply_sgs_stress(self, backend, mesh, du, dv, dw, u, v, w, &
                              xdirps, ydirps, zdirps)
    !! Add div(2*nut*S_ij) to the momentum right-hand-side fields.
    class(les_t), intent(inout) :: self
    class(base_backend_t), intent(inout) :: backend
    type(mesh_t), intent(in) :: mesh
    class(field_t), intent(inout) :: du, dv, dw
    class(field_t), intent(in) :: u, v, w
    type(dirps_t), intent(in) :: xdirps, ydirps, zdirps

    class(field_t), pointer :: dudx, dudy, dudz
    class(field_t), pointer :: dvdx, dvdy, dvdz
    class(field_t), pointer :: dwdx, dwdy, dwdz

    if (trim(self%model) == 'none') return
    if (trim(self%model) /= 'smagorinsky') &
      error stop 'Unsupported LES model in apply_sgs_stress.'
    if (du%dir /= DIR_X .or. dv%dir /= DIR_X .or. dw%dir /= DIR_X) then
      error stop 'LES momentum RHS fields must use DIR_X layout.'
    end if
    if (u%dir /= DIR_X .or. v%dir /= DIR_X .or. w%dir /= DIR_X) then
      error stop 'LES velocity fields must use DIR_X layout.'
    end if

    if (.not. associated(self%nut)) &
      self%nut => backend%allocator%get_block(DIR_X, VERT)
    if (.not. associated(self%mixing_length_sq)) then
      self%mixing_length_sq => backend%allocator%get_block(DIR_X, VERT)
      call initialise_mixing_length(self, backend, mesh)
    end if

    call compute_velocity_gradients( &
      backend, u, v, w, xdirps, ydirps, zdirps, &
      dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz)
    call backend%compute_smagorinsky_nut( &
      self%nut, self%mixing_length_sq, &
      dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz)

    if (self%neutral_wall_model) then
      ! Keep the SGS contribution separate until the bottom boundary has been
      ! replaced. This prevents the wall model from overwriting convection or
      ! molecular diffusion already present in du/dv/dw.
      sgs_du => backend%allocator%get_block(DIR_X, VERT)
      sgs_dv => backend%allocator%get_block(DIR_X, VERT)
      sgs_dw => backend%allocator%get_block(DIR_X, VERT)
      call sgs_du%fill(0._dp)
      call sgs_dv%fill(0._dp)
      call sgs_dw%fill(0._dp)

      call add_sgs_terms(backend, sgs_du, sgs_dv, sgs_dw, self%nut, &
                         dudx, dudy, dudz, dvdx, dvdy, dvdz, &
                         dwdx, dwdy, dwdz, xdirps, ydirps, zdirps)

      sampling_height = 0.5_dp*abs(mesh%geo%vert_coords(2, 2) - &
                                   mesh%geo%vert_coords(1, 2))
      call backend%apply_neutral_wall_flux( &
        sgs_du, sgs_dv, sgs_dw, u, w, self%nut, &
        dudy, dvdx, dwdy, dvdz, self%wall_model_kappa, &
        self%wall_model_roughness, sampling_height)

      call backend%vecadd(1._dp, sgs_du, 1._dp, du)
      call backend%vecadd(1._dp, sgs_dv, 1._dp, dv)
      call backend%vecadd(1._dp, sgs_dw, 1._dp, dw)
      call backend%allocator%release_block(sgs_du)
      call backend%allocator%release_block(sgs_dv)
      call backend%allocator%release_block(sgs_dw)
    else
      call add_sgs_terms(backend, du, dv, dw, self%nut, &
                         dudx, dudy, dudz, dvdx, dvdy, dvdz, &
                         dwdx, dwdy, dwdz, xdirps, ydirps, zdirps)
    end if

    call release_velocity_gradients( &
      backend, dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz)
  end subroutine apply_sgs_stress

  subroutine add_sgs_terms(backend, du, dv, dw, nut, &
                           dudx, dudy, dudz, dvdx, dvdy, dvdz, &
                           dwdx, dwdy, dwdz, xdirps, ydirps, zdirps)
    class(base_backend_t), intent(inout) :: backend
    class(field_t), intent(inout) :: du, dv, dw
    class(field_t), intent(in) :: nut
    class(field_t), intent(in) :: dudx, dudy, dudz
    class(field_t), intent(in) :: dvdx, dvdy, dvdz
    class(field_t), intent(in) :: dwdx, dwdy, dwdz
    type(dirps_t), intent(in) :: xdirps, ydirps, zdirps

    call add_normal_stress(backend, du, nut, dudx, xdirps)
    call add_normal_stress(backend, dv, nut, dvdy, ydirps)
    call add_normal_stress(backend, dw, nut, dwdz, zdirps)
    call add_shear_stress(backend, du, ydirps, dv, xdirps, nut, dudy, dvdx)
    call add_shear_stress(backend, du, zdirps, dw, xdirps, nut, dudz, dwdx)
    call add_shear_stress(backend, dv, zdirps, dw, ydirps, nut, dvdz, dwdy)
  end subroutine add_sgs_terms

  subroutine compute_velocity_gradients( &
    backend, u, v, w, xdirps, ydirps, zdirps, &
    dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz)
    class(base_backend_t), intent(inout) :: backend
    class(field_t), intent(in) :: u, v, w
    type(dirps_t), intent(in) :: xdirps, ydirps, zdirps
    class(field_t), pointer, intent(out) :: dudx, dudy, dudz
    class(field_t), pointer, intent(out) :: dvdx, dvdy, dvdz
    class(field_t), pointer, intent(out) :: dwdx, dwdy, dwdz

    call derivative_to_x(backend, dudx, u, xdirps)
    call derivative_to_x(backend, dvdx, v, xdirps)
    call derivative_to_x(backend, dwdx, w, xdirps)
    call derivative_to_x(backend, dudy, u, ydirps)
    call derivative_to_x(backend, dvdy, v, ydirps)
    call derivative_to_x(backend, dwdy, w, ydirps)
    call derivative_to_x(backend, dudz, u, zdirps)
    call derivative_to_x(backend, dvdz, v, zdirps)
    call derivative_to_x(backend, dwdz, w, zdirps)
  end subroutine compute_velocity_gradients

  subroutine release_velocity_gradients( &
    backend, dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz)
    class(base_backend_t), intent(inout) :: backend
    class(field_t), pointer, intent(inout) :: dudx, dudy, dudz
    class(field_t), pointer, intent(inout) :: dvdx, dvdy, dvdz
    class(field_t), pointer, intent(inout) :: dwdx, dwdy, dwdz

    call backend%allocator%release_block(dudx)
    call backend%allocator%release_block(dudy)
    call backend%allocator%release_block(dudz)
    call backend%allocator%release_block(dvdx)
    call backend%allocator%release_block(dvdy)
    call backend%allocator%release_block(dvdz)
    call backend%allocator%release_block(dwdx)
    call backend%allocator%release_block(dwdy)
    call backend%allocator%release_block(dwdz)
  end subroutine release_velocity_gradients

  subroutine add_normal_stress(backend, rhs, nut, gradient, direction)
    class(base_backend_t), intent(inout) :: backend
    class(field_t), intent(inout) :: rhs
    class(field_t), intent(in) :: nut, gradient
    type(dirps_t), intent(in) :: direction

    class(field_t), pointer :: stress

    stress => backend%allocator%get_block(DIR_X, VERT)
    call backend%compute_sgs_stress( &
      stress, nut, gradient, gradient, 2._dp, 0._dp)
    call add_stress_derivative(backend, rhs, stress, direction)
    call backend%allocator%release_block(stress)
  end subroutine add_normal_stress

  subroutine add_shear_stress( &
    backend, rhs_a, direction_a, rhs_b, &
    direction_b, nut, gradient_a, gradient_b &
    )
    class(base_backend_t), intent(inout) :: backend
    class(field_t), intent(inout) :: rhs_a, rhs_b
    type(dirps_t), intent(in) :: direction_a, direction_b
    class(field_t), intent(in) :: nut, gradient_a, gradient_b

    class(field_t), pointer :: stress

    stress => backend%allocator%get_block(DIR_X, VERT)
    call backend%compute_sgs_stress( &
      stress, nut, gradient_a, gradient_b, 1._dp, 1._dp)
    call add_stress_derivative(backend, rhs_a, stress, direction_a)
    call add_stress_derivative(backend, rhs_b, stress, direction_b)
    call backend%allocator%release_block(stress)
  end subroutine add_shear_stress

  subroutine add_stress_derivative(backend, rhs, stress, direction)
    class(base_backend_t), intent(inout) :: backend
    class(field_t), intent(inout) :: rhs
    class(field_t), intent(in) :: stress
    type(dirps_t), intent(in) :: direction

    class(field_t), pointer :: derivative

    call derivative_to_x(backend, derivative, stress, direction)
    call backend%vecadd(1._dp, derivative, 1._dp, rhs)
    call backend%allocator%release_block(derivative)
  end subroutine add_stress_derivative

  subroutine derivative_to_x(backend, derivative, velocity, dirps)
    class(base_backend_t), intent(inout) :: backend
    class(field_t), pointer, intent(out) :: derivative
    class(field_t), intent(in) :: velocity
    type(dirps_t), intent(in) :: dirps

    class(field_t), pointer :: velocity_dir, derivative_dir

    derivative => backend%allocator%get_block(DIR_X)
    select case (dirps%dir)
    case (DIR_X)
      call backend%tds_solve(derivative, velocity, dirps%der1st)
    case (DIR_Y)
      velocity_dir => backend%allocator%get_block(DIR_Y)
      derivative_dir => backend%allocator%get_block(DIR_Y)
      call backend%reorder(velocity_dir, velocity, RDR_X2Y)
      call backend%tds_solve(derivative_dir, velocity_dir, dirps%der1st)
      call backend%reorder(derivative, derivative_dir, RDR_Y2X)
      call backend%allocator%release_block(velocity_dir)
      call backend%allocator%release_block(derivative_dir)
    case (DIR_Z)
      velocity_dir => backend%allocator%get_block(DIR_Z)
      derivative_dir => backend%allocator%get_block(DIR_Z)
      call backend%reorder(velocity_dir, velocity, RDR_X2Z)
      call backend%tds_solve(derivative_dir, velocity_dir, dirps%der1st)
      call backend%reorder(derivative, derivative_dir, RDR_Z2X)
      call backend%allocator%release_block(velocity_dir)
      call backend%allocator%release_block(derivative_dir)
    case default
      error stop 'Invalid derivative direction in LES.'
    end select
  end subroutine derivative_to_x

  subroutine initialise_mixing_length(self, backend, mesh)
    class(les_t), intent(in) :: self
    class(base_backend_t), intent(inout) :: backend
    type(mesh_t), intent(in) :: mesh

    real(dp), allocatable :: mixing_data(:, :, :)
    real(dp) :: spacing(3), wall_distance, y_lower
    integer :: dims(3), dims_padded(3), i, j, k

    dims = mesh%get_dims(VERT)
    dims_padded = backend%allocator%get_padded_dims(DIR_C)
    allocate (mixing_data(dims_padded(1), dims_padded(2), dims_padded(3)))
    mixing_data = 0._dp

    y_lower = 0._dp
    if (trim(mesh%geo%stretching(2)) == 'centred') &
      y_lower = -0.5_dp*mesh%geo%L(2)

    do k = 1, dims(3)
      spacing(3) = spacing_at_vertex(mesh, k, 3, dims(3))
      do j = 1, dims(2)
        spacing(2) = spacing_at_vertex(mesh, j, 2, dims(2))
        wall_distance = max(mesh%geo%vert_coords(j, 2) - y_lower, 0._dp)
        do i = 1, dims(1)
          spacing(1) = spacing_at_vertex(mesh, i, 1, dims(1))
          mixing_data(i, j, k) = self%mixing_length(spacing, wall_distance)**2
        end do
      end do
    end do

    call backend%set_field_data(self%mixing_length_sq, mixing_data)
    deallocate (mixing_data)
  end subroutine initialise_mixing_length

  pure real(dp) function spacing_at_vertex(mesh, index, direction, n) &
    result(spacing)
    type(mesh_t), intent(in) :: mesh
    integer, intent(in) :: index, direction, n

    if (.not. mesh%geo%stretched(direction) .or. n == 1) then
      spacing = mesh%geo%d(direction)
    else if (index < n) then
      spacing = abs(mesh%geo%vert_coords(index + 1, direction) - &
                    mesh%geo%vert_coords(index, direction))
    else
      spacing = abs(mesh%geo%vert_coords(index, direction) - &
                    mesh%geo%vert_coords(index - 1, direction))
    end if
  end function spacing_at_vertex

  subroutine finalise(self, backend)
    class(les_t), intent(inout) :: self
    class(base_backend_t), intent(inout) :: backend

    if (associated(self%nut)) then
      call backend%allocator%release_block(self%nut)
      nullify (self%nut)
    end if
    if (associated(self%mixing_length_sq)) then
      call backend%allocator%release_block(self%mixing_length_sq)
      nullify (self%mixing_length_sq)
    end if
  end subroutine finalise

end module m_les
