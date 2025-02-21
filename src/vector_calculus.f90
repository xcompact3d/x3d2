module m_vector_calculus
  use iso_fortran_env, only: stderr => error_unit

  use m_allocator, only: allocator_t
  use m_base_backend, only: base_backend_t
  use m_common, only: dp, DIR_X, DIR_Y, DIR_Z, &
                      RDR_X2Y, RDR_X2Z, RDR_Y2X, RDR_Y2Z, RDR_Z2X, RDR_Z2Y
  use m_field, only: field_t
  use m_tdsops, only: tdsops_t

  implicit none

  type :: vector_calculus_t
    !! Defines vector calculus operators
    class(base_backend_t), pointer :: backend
  contains
    procedure :: curl
    procedure :: divergence_v2c
    procedure :: gradient_c2v
    procedure :: laplacian
  end type vector_calculus_t

  interface vector_calculus_t
    module procedure init
  end interface vector_calculus_t

contains

  function init(backend) result(vector_calculus)
    implicit none

    class(base_backend_t), target, intent(inout) :: backend
    type(vector_calculus_t) :: vector_calculus

    vector_calculus%backend => backend

  end function init

  subroutine curl(self, o_i_hat, o_j_hat, o_k_hat, u, v, w, &
                  x_der1st, y_der1st, z_der1st)
    !! Curl of a vector field (u, v, w).
    !!
    !! Evaluated at the data_loc defined by u, v, w fields.
    !!
    !! All the input and output fields are in DIR_X layout.
    implicit none

    class(vector_calculus_t) :: self
    !> Vector components of the output vector field Omega
    class(field_t), intent(inout) :: o_i_hat, o_j_hat, o_k_hat
    class(field_t), intent(in) :: u, v, w
    class(tdsops_t), intent(in) :: x_der1st, y_der1st, z_der1st

    class(field_t), pointer :: u_y, u_z, v_z, w_y, dwdy_y, dvdz_z, dvdz_x, &
      dudz_z, dudz_x, dudy_y, dudy_x

    if (o_i_hat%dir /= DIR_X .or. o_j_hat%dir /= DIR_X &
        .or. o_k_hat%dir /= DIR_X .or. u%dir /= DIR_X .or. v%dir /= DIR_X &
        .or. w%dir /= DIR_X) then
      error stop 'Error in curl input/output field %dirs: &
                  &outputs and inputs must be in DIR_X layout.'
    end if

    ! omega_i_hat = dw/dy - dv/dz
    ! omega_j_hat = du/dz - dw/dx
    ! omega_k_hat = dv/dx - du/dy

    ! omega_i_hat
    ! dw/dy
    w_y => self%backend%allocator%get_block(DIR_Y)
    dwdy_y => self%backend%allocator%get_block(DIR_Y)
    call self%backend%reorder(w_y, w, RDR_X2Y)
    call self%backend%tds_solve(dwdy_y, w_y, y_der1st)

    call self%backend%reorder(o_i_hat, dwdy_y, RDR_Y2X)

    call self%backend%allocator%release_block(w_y)
    call self%backend%allocator%release_block(dwdy_y)

    ! dv/dz
    v_z => self%backend%allocator%get_block(DIR_Z)
    dvdz_z => self%backend%allocator%get_block(DIR_Z)
    call self%backend%reorder(v_z, v, RDR_X2Z)
    call self%backend%tds_solve(dvdz_z, v_z, z_der1st)

    dvdz_x => self%backend%allocator%get_block(DIR_X)
    call self%backend%reorder(dvdz_x, dvdz_z, RDR_Z2X)

    call self%backend%allocator%release_block(v_z)
    call self%backend%allocator%release_block(dvdz_z)

    ! omega_i_hat = dw/dy - dv/dz
    call self%backend%vecadd(-1._dp, dvdz_x, 1._dp, o_i_hat)

    call self%backend%allocator%release_block(dvdz_x)

    ! omega_j_hat
    ! du/dz
    u_z => self%backend%allocator%get_block(DIR_Z)
    dudz_z => self%backend%allocator%get_block(DIR_Z)
    call self%backend%reorder(u_z, u, RDR_X2Z)
    call self%backend%tds_solve(dudz_z, u_z, z_der1st)

    dudz_x => self%backend%allocator%get_block(DIR_X)
    call self%backend%reorder(dudz_x, dudz_z, RDR_Z2X)

    call self%backend%allocator%release_block(u_z)
    call self%backend%allocator%release_block(dudz_z)

    ! dw/dx
    call self%backend%tds_solve(o_j_hat, w, x_der1st)

    ! omega_j_hat = du/dz - dw/dx
    call self%backend%vecadd(1._dp, dudz_x, -1._dp, o_j_hat)

    call self%backend%allocator%release_block(dudz_x)

    ! omega_k_hat
    ! dv/dx
    call self%backend%tds_solve(o_k_hat, v, x_der1st)

    ! du/dy
    u_y => self%backend%allocator%get_block(DIR_Y)
    dudy_y => self%backend%allocator%get_block(DIR_Y)
    call self%backend%reorder(u_y, u, RDR_X2Y)
    call self%backend%tds_solve(dudy_y, u_y, y_der1st)

    dudy_x => self%backend%allocator%get_block(DIR_X)
    call self%backend%reorder(dudy_x, dudy_y, RDR_Y2X)

    call self%backend%allocator%release_block(u_y)
    call self%backend%allocator%release_block(dudy_y)

    ! omega_k_hat = dv/dx - du/dy
    call self%backend%vecadd(-1._dp, dudy_x, 1._dp, o_k_hat)

    call self%backend%allocator%release_block(dudy_x)

  end subroutine curl

  subroutine divergence_v2c(self, div_u, u, v, w, &
                            x_stagder_v2c, x_interpl_v2c, &
                            y_stagder_v2c, y_interpl_v2c, &
                            z_stagder_v2c, z_interpl_v2c)
    !! Divergence of a vector field (u, v, w).
    !!
    !! Evaluated at the cell centers (data_loc=CELL)
    !! Input fields are at vertices (data_loc=VERT)
    !!
    !! Input fields are in DIR_X data layout.
    !! Output field is in DIR_Z data layout.
    implicit none

    class(vector_calculus_t) :: self
    class(field_t), intent(inout) :: div_u
    class(field_t), intent(in) :: u, v, w
    class(tdsops_t), intent(in) :: x_stagder_v2c, x_interpl_v2c, &
      y_stagder_v2c, y_interpl_v2c, &
      z_stagder_v2c, z_interpl_v2c

    class(field_t), pointer :: du_x, dv_x, dw_x, &
      u_y, v_y, w_y, du_y, dv_y, dw_y, &
      u_z, w_z, dw_z

    if (div_u%dir /= DIR_Z .or. u%dir /= DIR_X .or. v%dir /= DIR_X &
        .or. w%dir /= DIR_X) then
      error stop 'Error in divergence_v2c input/output field dirs: &
                  &output must be in DIR_Z, inputs must be in DIR_X layout.'
    end if

    du_x => self%backend%allocator%get_block(DIR_X)
    dv_x => self%backend%allocator%get_block(DIR_X)
    dw_x => self%backend%allocator%get_block(DIR_X)

    ! Staggared der for u field in x
    ! Interpolation for v field in x
    ! Interpolation for w field in x
    call self%backend%tds_solve(du_x, u, x_stagder_v2c)
    call self%backend%tds_solve(dv_x, v, x_interpl_v2c)
    call self%backend%tds_solve(dw_x, w, x_interpl_v2c)

    ! request fields from the allocator
    u_y => self%backend%allocator%get_block(DIR_Y)
    v_y => self%backend%allocator%get_block(DIR_Y)
    w_y => self%backend%allocator%get_block(DIR_Y)

    ! reorder data from x orientation to y orientation
    call self%backend%reorder(u_y, du_x, RDR_X2Y)
    call self%backend%reorder(v_y, dv_x, RDR_X2Y)
    call self%backend%reorder(w_y, dw_x, RDR_X2Y)

    call self%backend%allocator%release_block(du_x)
    call self%backend%allocator%release_block(dv_x)
    call self%backend%allocator%release_block(dw_x)

    du_y => self%backend%allocator%get_block(DIR_Y)
    dv_y => self%backend%allocator%get_block(DIR_Y)
    dw_y => self%backend%allocator%get_block(DIR_Y)

    ! similar to the x direction, obtain derivatives in y.
    call self%backend%tds_solve(du_y, u_y, y_interpl_v2c)
    call self%backend%tds_solve(dv_y, v_y, y_stagder_v2c)
    call self%backend%tds_solve(dw_y, w_y, y_interpl_v2c)

    ! we don't need the velocities in y orientation any more, so release
    ! them to open up space.
    ! It is important that this doesn't actually deallocate any memory,
    ! it just makes the corresponding memory space available for use.
    call self%backend%allocator%release_block(u_y)
    call self%backend%allocator%release_block(v_y)
    call self%backend%allocator%release_block(w_y)

    ! just like in y direction, get some fields for the z derivatives.
    u_z => self%backend%allocator%get_block(DIR_Z)
    w_z => self%backend%allocator%get_block(DIR_Z)

    ! du_y = dv_y + du_y
    call self%backend%vecadd(1._dp, dv_y, 1._dp, du_y)

    ! reorder from y to z
    call self%backend%reorder(u_z, du_y, RDR_Y2Z)
    call self%backend%reorder(w_z, dw_y, RDR_Y2Z)

    ! release all the unnecessary blocks.
    call self%backend%allocator%release_block(du_y)
    call self%backend%allocator%release_block(dv_y)
    call self%backend%allocator%release_block(dw_y)

    dw_z => self%backend%allocator%get_block(DIR_Z)

    ! get the derivatives in z
    call self%backend%tds_solve(div_u, u_z, z_interpl_v2c)
    call self%backend%tds_solve(dw_z, w_z, z_stagder_v2c)

    ! div_u = div_u + dw_z
    call self%backend%vecadd(1._dp, dw_z, 1._dp, div_u)

    ! div_u array is in z orientation

    ! there is no need to keep velocities in z orientation around, so release
    call self%backend%allocator%release_block(u_z)
    call self%backend%allocator%release_block(w_z)
    call self%backend%allocator%release_block(dw_z)

  end subroutine divergence_v2c

  subroutine gradient_c2v(self, dpdx, dpdy, dpdz, p, &
                          x_stagder_c2v, x_interpl_c2v, &
                          y_stagder_c2v, y_interpl_c2v, &
                          z_stagder_c2v, z_interpl_c2v)
    !! Gradient of a scalar field 'p'.
    !!
    !! Evaluated at the vertices (data_loc=VERT)
    !! Input field is at cell centers (data_loc=CELL)
    !!
    !! Input field is in DIR_Z data layout.
    !! Output fields (dpdx, dpdy, dpdz) are in DIR_X data layout.
    implicit none

    class(vector_calculus_t) :: self
    class(field_t), intent(inout) :: dpdx, dpdy, dpdz
    class(field_t), intent(in) :: p
    class(tdsops_t), intent(in) :: x_stagder_c2v, x_interpl_c2v, &
      y_stagder_c2v, y_interpl_c2v, &
      z_stagder_c2v, z_interpl_c2v

    class(field_t), pointer :: p_sxy_z, dpdz_sxy_z, &
      p_sxy_y, dpdz_sxy_y, &
      p_sx_y, dpdy_sx_y, dpdz_sx_y, &
      p_sx_x, dpdy_sx_x, dpdz_sx_x

    if (dpdx%dir /= DIR_X .or. dpdy%dir /= DIR_X .or. dpdz%dir /= DIR_X &
        .or. p%dir /= DIR_Z) then
      error stop 'Error in gradient_c2v input/output field dirs: &
                  &outputs must be in DIR_X, input must be in DIR_Z layout.'
    end if

    p_sxy_z => self%backend%allocator%get_block(DIR_Z)
    dpdz_sxy_z => self%backend%allocator%get_block(DIR_Z)

    ! Staggared der for p field in z
    ! Interpolation for p field in z
    call self%backend%tds_solve(p_sxy_z, p, z_interpl_c2v)
    call self%backend%tds_solve(dpdz_sxy_z, p, z_stagder_c2v)

    ! request fields from the allocator
    p_sxy_y => self%backend%allocator%get_block(DIR_Y)
    dpdz_sxy_y => self%backend%allocator%get_block(DIR_Y)

    ! reorder data from z orientation to y orientation
    call self%backend%reorder(p_sxy_y, p_sxy_z, RDR_Z2Y)
    call self%backend%reorder(dpdz_sxy_y, dpdz_sxy_z, RDR_Z2Y)

    call self%backend%allocator%release_block(p_sxy_z)
    call self%backend%allocator%release_block(dpdz_sxy_z)

    p_sx_y => self%backend%allocator%get_block(DIR_Y)
    dpdy_sx_y => self%backend%allocator%get_block(DIR_Y)
    dpdz_sx_y => self%backend%allocator%get_block(DIR_Y)

    ! similar to the z direction, obtain derivatives in y.
    call self%backend%tds_solve(p_sx_y, p_sxy_y, y_interpl_c2v)
    call self%backend%tds_solve(dpdy_sx_y, p_sxy_y, y_stagder_c2v)
    call self%backend%tds_solve(dpdz_sx_y, dpdz_sxy_y, y_interpl_c2v)

    ! release memory
    call self%backend%allocator%release_block(p_sxy_y)
    call self%backend%allocator%release_block(dpdz_sxy_y)

    ! just like in y direction, get some fields for the x derivatives.
    p_sx_x => self%backend%allocator%get_block(DIR_X)
    dpdy_sx_x => self%backend%allocator%get_block(DIR_X)
    dpdz_sx_x => self%backend%allocator%get_block(DIR_X)

    ! reorder from y to x
    call self%backend%reorder(p_sx_x, p_sx_y, RDR_Y2X)
    call self%backend%reorder(dpdy_sx_x, dpdy_sx_y, RDR_Y2X)
    call self%backend%reorder(dpdz_sx_x, dpdz_sx_y, RDR_Y2X)

    ! release all the y directional fields.
    call self%backend%allocator%release_block(p_sx_y)
    call self%backend%allocator%release_block(dpdy_sx_y)
    call self%backend%allocator%release_block(dpdz_sx_y)

    ! get the derivatives in x
    call self%backend%tds_solve(dpdx, p_sx_x, x_stagder_c2v)
    call self%backend%tds_solve(dpdy, dpdy_sx_x, x_interpl_c2v)
    call self%backend%tds_solve(dpdz, dpdz_sx_x, x_interpl_c2v)

    ! release temporary x fields
    call self%backend%allocator%release_block(p_sx_x)
    call self%backend%allocator%release_block(dpdy_sx_x)
    call self%backend%allocator%release_block(dpdz_sx_x)

  end subroutine gradient_c2v

  subroutine laplacian(self, lapl_u, u, x_der2nd, y_der2nd, z_der2nd)
    !! Laplacian of a scalar field 'u'.
    !!
    !! Evaluated at the data_loc defined by the input u field
    !!
    !! Input and output fields are in DIR_X layout.
    implicit none

    class(vector_calculus_t) :: self
    class(field_t), intent(inout) :: lapl_u
    class(field_t), intent(in) :: u

    class(tdsops_t), intent(in) :: x_der2nd, y_der2nd, z_der2nd

    class(field_t), pointer :: u_y, d2u_y, u_z, d2u_z

    if (u%dir /= DIR_X .or. lapl_u%dir /= DIR_X) then
      error stop 'Error in laplacian input/output field %dirs: &
                  &outputs and inputs must be in DIR_X layout.'
    end if

    ! d2u/dx2
    call self%backend%tds_solve(lapl_u, u, x_der2nd)

    ! y directional temporary fields
    u_y => self%backend%allocator%get_block(DIR_Y)
    d2u_y => self%backend%allocator%get_block(DIR_Y)

    call self%backend%reorder(u_y, u, RDR_X2Y)

    ! d2u/dy2
    call self%backend%tds_solve(d2u_y, u_y, y_der2nd)

    ! add the y derivative component into the result field
    call self%backend%sum_yintox(lapl_u, d2u_y)

    ! release y directional fields
    call self%backend%allocator%release_block(u_y)
    call self%backend%allocator%release_block(d2u_y)

    ! z directional temporary fields
    u_z => self%backend%allocator%get_block(DIR_Z)
    d2u_z => self%backend%allocator%get_block(DIR_Z)

    call self%backend%reorder(u_z, u, RDR_X2Z)

    ! d2u/dz2
    call self%backend%tds_solve(d2u_z, u_z, z_der2nd)

    ! add the z derivative component into the result field
    call self%backend%sum_zintox(lapl_u, d2u_z)

    ! release z directional fields
    call self%backend%allocator%release_block(u_z)
    call self%backend%allocator%release_block(d2u_z)

  end subroutine laplacian

end module m_vector_calculus
