module m_velocity
  use m_vector3d, only: vector3d
  implicit none


  type :: velocity_field
     class(vector3d) :: x_slab, y_slab, z_slab
  end type velocity_field

contains

  pure function F(self)
    class(velocity_field), intent(in) :: self
    type(velocity_field) :: tranported_field
    call async_rotate_x_to_z(x_slab_x)
    x_dir_contrib = x_slab%transport()
    call from_comms_buffer(z_slab%get("u"))
    call from_comms_buffer(z_slab%get("v"))
    call from_comms_buffer(z_slab%get("w"))
    z_dir_contrib = z_slab%transport(diffeng)
    call async_rotate_z_to_x(z_dir_contrib, x_slab_z)
    call transpose(x_slab, y_slab)
    transported_field%x_slab = &
         sum_contributions(x_slab, y_slab, x_slab_z)
  end function F

  subroutine setup(dims, ybcs, zbcs)
    m_y_direction_buffer = get_slab(dims(2), ybcs)
    m_z_direction_buffer = get_slab(dims(3), zbcs)
  end subroutine setup
end module m_velocity
