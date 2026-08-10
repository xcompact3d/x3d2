module m_les
  !! Shared physics for explicit eddy-viscosity LES models.
  !!
  !! The field-level implementation will use backend derivative and pointwise
  !! operations.  Keeping these scalar relations here gives both backends one
  !! definition of the model and makes the closure independently testable.
  use m_common, only: dp
  use m_config, only: les_config_t

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
  contains
    procedure :: nut_from_gradient
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
    delta, wall_distance, smagorinsky_constant, von_karman_constant, exponent, &
    roughness_length) result(length)
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
    velocity_gradient, mixing_length) result(nut)
    real(dp), intent(in) :: velocity_gradient(3, 3)
    real(dp), intent(in) :: mixing_length

    nut = mixing_length**2*strain_rate_magnitude(velocity_gradient)
  end function smagorinsky_nut

  pure real(dp) function nut_from_gradient( &
    self, velocity_gradient, spacing, wall_distance) result(nut)
    class(les_t), intent(in) :: self
    real(dp), intent(in) :: velocity_gradient(3, 3), spacing(3)
    real(dp), optional, intent(in) :: wall_distance
    real(dp) :: delta, length

    if (trim(self%model) == 'none') then
      nut = 0._dp
      return
    end if

    delta = filter_width(spacing)
    length = self%smagorinsky_constant*delta
    if (self%wall_damping) then
      if (.not. present(wall_distance)) then
        nut = 0._dp
        return
      end if
      length = wall_damped_mixing_length(delta, wall_distance, &
        self%smagorinsky_constant, self%von_karman_constant, &
        self%wall_damping_n, self%roughness_length)
    end if
    nut = smagorinsky_nut(velocity_gradient, length)
  end function nut_from_gradient

end module m_les
