module m_tdsops_setup
  use m_common, only: dp, VERT, CELL, BC_NEUMANN, BC_DIRICHLET
  use m_base_backend, only: base_backend_t
  use m_mesh, only: mesh_t
  use m_tdsops, only: dirps_t

  implicit none

contains

  subroutine allocate_tdsops(dirps, backend, mesh, der1st_scheme, &
                             der2nd_scheme, interpl_scheme, stagder_scheme)
    type(dirps_t), intent(inout) :: dirps
    class(base_backend_t), intent(in) :: backend
    type(mesh_t), intent(in) :: mesh
    character(*), intent(in) :: der1st_scheme, der2nd_scheme, &
                                interpl_scheme, stagder_scheme

    integer :: dir, bc_start, bc_end, bc_mp_start, bc_mp_end, n_vert, n_cell, i
    real(dp) :: d

    dir = dirps%dir
    d = mesh%geo%d(dir)

    bc_start = mesh%grid%BCs(dir, 1)
    bc_end = mesh%grid%BCs(dir, 2)

    ! For the FFT based poisson solver, the BC for the pressure has to be
    ! Neumann. This is not strictly compatible with Navier-Stokes equations,
    ! but it does not affect the quality of the simulation.
    ! Thus, if the BC is Dirichlet, we enforce Neumann for the midpoint &
    ! operators. (It could be BC_HALO too, if so we just keep it as is)
    if (bc_start == BC_DIRICHLET) then
      bc_mp_start = BC_NEUMANN
    else
      bc_mp_start = bc_start
    end if
    if (bc_end == BC_DIRICHLET) then
      bc_mp_end = BC_NEUMANN
    else
      bc_mp_end = bc_end
    end if

    n_vert = mesh%get_n(dir, VERT)
    n_cell = mesh%get_n(dir, CELL)

    call backend%alloc_tdsops( &
      dirps%der1st, n_vert, d, 'first-deriv', der1st_scheme, &
      bc_start, bc_end, stretch=mesh%geo%vert_ds(1:n_vert, dir) &
      )
    call backend%alloc_tdsops( &
      dirps%der1st_sym, n_vert, d, 'first-deriv', der1st_scheme, &
      bc_start, bc_end, stretch=mesh%geo%vert_ds(1:n_vert, dir), &
      sym=.true. &
      )
    call backend%alloc_tdsops( &
      dirps%der2nd, n_vert, d, 'second-deriv', der2nd_scheme, &
      bc_start, bc_end, stretch=mesh%geo%vert_ds2(1:n_vert, dir), &
      stretch_correct=mesh%geo%vert_d2s(1:n_vert, dir) &
      )
    call backend%alloc_tdsops( &
      dirps%der2nd_sym, n_vert, d, 'second-deriv', der2nd_scheme, &
      bc_start, bc_end, stretch=mesh%geo%vert_ds2(1:n_vert, dir), &
      stretch_correct=mesh%geo%vert_d2s(1:n_vert, dir), &
      sym=.true. &
      )
    call backend%alloc_tdsops( &
      dirps%stagder_v2p, n_cell, d, 'stag-deriv', stagder_scheme, &
      bc_mp_start, bc_mp_end, from_to='v2p', &
      stretch=mesh%geo%midp_ds(1:n_cell, dir) &
      )
    call backend%alloc_tdsops( &
      dirps%stagder_p2v, n_vert, d, 'stag-deriv', stagder_scheme, &
      bc_mp_start, bc_mp_end, from_to='p2v', &
      stretch=mesh%geo%vert_ds(1:n_vert, dir) &
      )
    call backend%alloc_tdsops( &
      dirps%interpl_v2p, n_cell, d, 'interpolate', interpl_scheme, &
      bc_mp_start, bc_mp_end, from_to='v2p', stretch=[(1._dp, i=1, n_cell)] &
      )
    call backend%alloc_tdsops( &
      dirps%interpl_p2v, n_vert, d, 'interpolate', interpl_scheme, &
      bc_mp_start, bc_mp_end, from_to='p2v', stretch=[(1._dp, i=1, n_vert)] &
      )

  end subroutine

end module m_tdsops_setup
