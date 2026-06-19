module m_postprocess
  !! Computation of derived fields for snapshot output.
  !!
  !! Provides routines to compute quantities that are not part of
  !! the core time-stepping algorithm but are useful for analysis
  !! and visualisation (e.g. pressure on the vertex grid, vorticity
  !! magnitude, Q-criterion).

  use m_common, only: dp, DIR_X, DIR_Y, DIR_Z, DIR_C, VERT, &
                      RDR_X2Y, RDR_X2Z, RDR_Y2X, RDR_Z2X
  use m_field, only: field_t
  use m_solver, only: solver_t

  implicit none

  private
  public :: compute_derived_fields, compute_pressure_vert

contains

  subroutine compute_derived_fields(solver, output_vorticity, &
                                    output_qcriterion)
    !! Compute derived fields from the velocity gradient tensor.
    !!
    !! All 9 components of the velocity gradient tensor are computed
    !! once, then used to evaluate whichever quantities are enabled
    !! (vorticity magnitude, Q-criterion).
    implicit none

    class(solver_t), intent(inout) :: solver
    logical, intent(in) :: output_vorticity, output_qcriterion

    class(field_t), pointer :: &
      dudx, dudy_x, dudz_x, &
      dvdx, dvdy_x, dvdz_x, &
      dwdx, dwdy_x, dwdz_x
    class(field_t), pointer :: &
      u_y, dudy_y, u_z, dudz_z, &
      v_y, dvdy_y, v_z, dvdz_z, &
      w_y, dwdy_y, w_z, dwdz_z
    class(field_t), pointer :: &
      dudx_h, dudy_h, dudz_h, &
      dvdx_h, dvdy_h, dvdz_h, &
      dwdx_h, dwdy_h, dwdz_h, &
      vort_h, qcrit_h

    ! Allocate output fields on first use
    if (output_vorticity .and. .not. associated(solver%vort)) &
      solver%vort => solver%backend%allocator%get_block(DIR_X, VERT)
    if (output_qcriterion .and. .not. associated(solver%qcrit)) &
      solver%qcrit => solver%backend%allocator%get_block(DIR_X, VERT)

    ! Compute all 9 velocity gradient components
    ! du/dx (x-derivative in x-pencil, no reorder needed)
    dudx => solver%backend%allocator%get_block(DIR_X)
    call solver%backend%tds_solve(dudx, solver%u, &
                                  solver%xdirps%der1st)

    ! dv/dx
    dvdx => solver%backend%allocator%get_block(DIR_X)
    call solver%backend%tds_solve(dvdx, solver%v, &
                                  solver%xdirps%der1st)

    ! dw/dx
    dwdx => solver%backend%allocator%get_block(DIR_X)
    call solver%backend%tds_solve(dwdx, solver%w, &
                                  solver%xdirps%der1st)

    ! du/dy
    u_y => solver%backend%allocator%get_block(DIR_Y)
    dudy_y => solver%backend%allocator%get_block(DIR_Y)
    call solver%backend%reorder(u_y, solver%u, RDR_X2Y)
    call solver%backend%tds_solve(dudy_y, u_y, &
                                  solver%ydirps%der1st)
    dudy_x => solver%backend%allocator%get_block(DIR_X)
    call solver%backend%reorder(dudy_x, dudy_y, RDR_Y2X)
    call solver%backend%allocator%release_block(u_y)
    call solver%backend%allocator%release_block(dudy_y)

    ! dv/dy
    v_y => solver%backend%allocator%get_block(DIR_Y)
    dvdy_y => solver%backend%allocator%get_block(DIR_Y)
    call solver%backend%reorder(v_y, solver%v, RDR_X2Y)
    call solver%backend%tds_solve(dvdy_y, v_y, &
                                  solver%ydirps%der1st)
    dvdy_x => solver%backend%allocator%get_block(DIR_X)
    call solver%backend%reorder(dvdy_x, dvdy_y, RDR_Y2X)
    call solver%backend%allocator%release_block(v_y)
    call solver%backend%allocator%release_block(dvdy_y)

    ! dw/dy
    w_y => solver%backend%allocator%get_block(DIR_Y)
    dwdy_y => solver%backend%allocator%get_block(DIR_Y)
    call solver%backend%reorder(w_y, solver%w, RDR_X2Y)
    call solver%backend%tds_solve(dwdy_y, w_y, &
                                  solver%ydirps%der1st)
    dwdy_x => solver%backend%allocator%get_block(DIR_X)
    call solver%backend%reorder(dwdy_x, dwdy_y, RDR_Y2X)
    call solver%backend%allocator%release_block(w_y)
    call solver%backend%allocator%release_block(dwdy_y)

    ! du/dz
    u_z => solver%backend%allocator%get_block(DIR_Z)
    dudz_z => solver%backend%allocator%get_block(DIR_Z)
    call solver%backend%reorder(u_z, solver%u, RDR_X2Z)
    call solver%backend%tds_solve(dudz_z, u_z, &
                                  solver%zdirps%der1st)
    dudz_x => solver%backend%allocator%get_block(DIR_X)
    call solver%backend%reorder(dudz_x, dudz_z, RDR_Z2X)
    call solver%backend%allocator%release_block(u_z)
    call solver%backend%allocator%release_block(dudz_z)

    ! dv/dz
    v_z => solver%backend%allocator%get_block(DIR_Z)
    dvdz_z => solver%backend%allocator%get_block(DIR_Z)
    call solver%backend%reorder(v_z, solver%v, RDR_X2Z)
    call solver%backend%tds_solve(dvdz_z, v_z, &
                                  solver%zdirps%der1st)
    dvdz_x => solver%backend%allocator%get_block(DIR_X)
    call solver%backend%reorder(dvdz_x, dvdz_z, RDR_Z2X)
    call solver%backend%allocator%release_block(v_z)
    call solver%backend%allocator%release_block(dvdz_z)

    ! dw/dz
    w_z => solver%backend%allocator%get_block(DIR_Z)
    dwdz_z => solver%backend%allocator%get_block(DIR_Z)
    call solver%backend%reorder(w_z, solver%w, RDR_X2Z)
    call solver%backend%tds_solve(dwdz_z, w_z, &
                                  solver%zdirps%der1st)
    dwdz_x => solver%backend%allocator%get_block(DIR_X)
    call solver%backend%reorder(dwdz_x, dwdz_z, RDR_Z2X)
    call solver%backend%allocator%release_block(w_z)
    call solver%backend%allocator%release_block(dwdz_z)

    ! All 9 derivatives now in x-pencil layout:
    ! dudx, dudy_x, dudz_x
    ! dvdx, dvdy_x, dvdz_x
    ! dwdx, dwdy_x, dwdz_x

    dudx_h => solver%host_allocator%get_block(DIR_C, dudx%data_loc)
    dudy_h => solver%host_allocator%get_block(DIR_C, dudy_x%data_loc)
    dudz_h => solver%host_allocator%get_block(DIR_C, dudz_x%data_loc)
    dvdx_h => solver%host_allocator%get_block(DIR_C, dvdx%data_loc)
    dvdy_h => solver%host_allocator%get_block(DIR_C, dvdy_x%data_loc)
    dvdz_h => solver%host_allocator%get_block(DIR_C, dvdz_x%data_loc)
    dwdx_h => solver%host_allocator%get_block(DIR_C, dwdx%data_loc)
    dwdy_h => solver%host_allocator%get_block(DIR_C, dwdy_x%data_loc)
    dwdz_h => solver%host_allocator%get_block(DIR_C, dwdz_x%data_loc)

    call solver%backend%get_field_data(dudx_h%data, dudx, DIR_C)
    call solver%backend%get_field_data(dudy_h%data, dudy_x, DIR_C)
    call solver%backend%get_field_data(dudz_h%data, dudz_x, DIR_C)
    call solver%backend%get_field_data(dvdx_h%data, dvdx, DIR_C)
    call solver%backend%get_field_data(dvdy_h%data, dvdy_x, DIR_C)
    call solver%backend%get_field_data(dvdz_h%data, dvdz_x, DIR_C)
    call solver%backend%get_field_data(dwdx_h%data, dwdx, DIR_C)
    call solver%backend%get_field_data(dwdy_h%data, dwdy_x, DIR_C)
    call solver%backend%get_field_data(dwdz_h%data, dwdz_x, DIR_C)

    !! Vorticity magnitude:
    !! \( |\boldsymbol{\omega}| = \sqrt{(\partial w/\partial y - \partial v/\partial z)^2
    !! + (\partial u/\partial z - \partial w/\partial x)^2
    !! + (\partial v/\partial x - \partial u/\partial y)^2} \)
    if (output_vorticity) then
      vort_h => solver%host_allocator%get_block(DIR_C, solver%vort%data_loc)
      vort_h%data = sqrt( &
                     (dwdy_h%data - dvdz_h%data)**2 + &
                     (dudz_h%data - dwdx_h%data)**2 + &
                     (dvdx_h%data - dudy_h%data)**2 &
                     )
      call solver%backend%set_field_data(solver%vort, vort_h%data, DIR_C)
      call solver%host_allocator%release_block(vort_h)
    end if

    !! Q-criterion:
    !! \( Q = -\frac{1}{2}(u_{x}^2 + v_{y}^2 + w_{z}^2)
    !! - u_{y}v_{x} - u_{z}w_{x} - v_{z}w_{y} \)
    if (output_qcriterion) then
      qcrit_h => solver%host_allocator%get_block(DIR_C, solver%qcrit%data_loc)
      qcrit_h%data = &
        -0.5_dp*(dudx_h%data**2 &
                 + dvdy_h%data**2 &
                 + dwdz_h%data**2) &
        - dudy_h%data*dvdx_h%data &
        - dudz_h%data*dwdx_h%data &
        - dvdz_h%data*dwdy_h%data
      call solver%backend%set_field_data(solver%qcrit, qcrit_h%data, DIR_C)
      call solver%host_allocator%release_block(qcrit_h)
    end if

    call solver%host_allocator%release_block(dudx_h)
    call solver%host_allocator%release_block(dudy_h)
    call solver%host_allocator%release_block(dudz_h)
    call solver%host_allocator%release_block(dvdx_h)
    call solver%host_allocator%release_block(dvdy_h)
    call solver%host_allocator%release_block(dvdz_h)
    call solver%host_allocator%release_block(dwdx_h)
    call solver%host_allocator%release_block(dwdy_h)
    call solver%host_allocator%release_block(dwdz_h)

    ! Release all gradient fields
    call solver%backend%allocator%release_block(dudx)
    call solver%backend%allocator%release_block(dvdx)
    call solver%backend%allocator%release_block(dwdx)
    call solver%backend%allocator%release_block(dudy_x)
    call solver%backend%allocator%release_block(dvdy_x)
    call solver%backend%allocator%release_block(dwdy_x)
    call solver%backend%allocator%release_block(dudz_x)
    call solver%backend%allocator%release_block(dvdz_x)
    call solver%backend%allocator%release_block(dwdz_x)

  end subroutine compute_derived_fields

  subroutine compute_pressure_vert(solver)
    !! Interpolates the pressure field from CELL (DIR_Z) to VERT
    !! (DIR_X) for snapshot output and rescales from pseudo-pressure
    !! to physical (kinematic) pressure.
    implicit none

    class(solver_t), intent(inout) :: solver

    if (.not. associated(solver%pressure)) then
      error stop 'compute_pressure_vert: pressure not yet computed'
    end if

    ! Allocate pressure_vert on first use
    if (.not. associated(solver%pressure_vert)) then
      solver%pressure_vert => &
        solver%backend%allocator%get_block(DIR_X, VERT)
    end if

    call solver%vector_calculus%interpl_c2v( &
      solver%pressure_vert, solver%pressure, &
      solver%xdirps%interpl_p2v, &
      solver%ydirps%interpl_p2v, &
      solver%zdirps%interpl_p2v &
      )

    ! Rescale pseudo-pressure to physical pressure: p'/dt -> p
    call solver%backend%vecadd( &
      1._dp/solver%dt, solver%pressure_vert, &
      0._dp, solver%pressure_vert &
      )

  end subroutine compute_pressure_vert

end module m_postprocess
