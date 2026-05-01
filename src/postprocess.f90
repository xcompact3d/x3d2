module m_postprocess
  !! Computation of derived fields for snapshot output.
  !!
  !! Provides routines to compute quantities that are not part of
  !! the core time-stepping algorithm but are useful for analysis
  !! and visualisation (e.g. pressure on the vertex grid, vorticity
  !! magnitude, Q-criterion).

  use m_common, only: dp, DIR_X, DIR_Y, DIR_Z, VERT, &
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

    ! Lazily allocate output fields
    if (output_vorticity .and. .not. associated(solver%vort)) &
      solver%vort => solver%backend%allocator%get_block(DIR_X, VERT)
    if (output_qcriterion .and. .not. associated(solver%qcrit)) &
      solver%qcrit => solver%backend%allocator%get_block(DIR_X, VERT)

    ! --- Compute all 9 velocity gradient components ---

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

    ! Vorticity magnitude:
    ! |w| = sqrt((dw/dy-dv/dz)^2+(du/dz-dw/dx)^2+(dv/dx-du/dy)^2)
    if (output_vorticity) then
      solver%vort%data = sqrt( &
                         (dwdy_x%data - dvdz_x%data)**2 + &
                         (dudz_x%data - dwdx%data)**2 + &
                         (dvdx%data - dudy_x%data)**2 &
                         )
    end if

    ! Q-criterion:
    ! Q = -0.5*(dudx^2+dvdy^2+dwdz^2)
    !     - dudy*dvdx - dudz*dwdx - dvdz*dwdy
    if (output_qcriterion) then
      solver%qcrit%data = &
        -0.5_dp*(dudx%data**2 &
                 + dvdy_x%data**2 &
                 + dwdz_x%data**2) &
        - dudy_x%data*dvdx%data &
        - dudz_x%data*dwdx%data &
        - dvdz_x%data*dwdy_x%data
    end if

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

    ! Lazy allocate pressure_vert on first call
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
