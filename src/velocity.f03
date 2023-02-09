module m_velocity_field
  use m_vector3d, only: vector3d
  implicit none


  type :: velocity_field
     class(vector3d) :: x_slab, y_slab, z_slab
  end type velocity_field

contains

  pure function F(self)
    class(velocity_field), intent(in) :: self
    type(velocity_field) :: tranported_field

    class(vector3d), allocatable :: x_dir_contrib, z_dir_contrib

    call async_rotate_x_to_z(x_slab_x)

    x_dir_contrib = x_slab%transport()

    call from_comms_buffer(z_slab%u(1))
    call from_comms_buffer(z_slab%u(2))
    call from_comms_buffer(z_slab%u(3))
    z_dir_contrib = z_slab%transport()

    call async_rotate_z_to_x(z_dir_contrib, x_slab_z)
    call rotate_x_to_y(x_slab, y_slab)

    transported_field%x_slab = &
         sum_contributions(x_slab, y_slab, x_slab_z)
  end function F
end module m_velocity_field

