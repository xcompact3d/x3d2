module m_tdsops_setup
  use m_common, only: dp, VERT, CELL, BC_NEUMANN, BC_DIRICHLET
  use m_base_backend, only: base_backend_t
  use m_mesh, only: mesh_t
  use m_tdsops, only: dirps_t

  implicit none

contains

  subroutine allocate_tdsops(dirps, backend, mesh, der1st_scheme, &
                             der2nd_scheme, interpl_scheme, stagder_scheme, &
                             c_nu, nu0_nu, filter_alpha)
    type(dirps_t), intent(inout) :: dirps
    class(base_backend_t), intent(in) :: backend
    type(mesh_t), intent(in) :: mesh
    character(*), intent(in) :: der1st_scheme, der2nd_scheme, &
                                interpl_scheme, stagder_scheme
    !! hyperviscous (iSVV) der2nd params; only required when der2nd_scheme is
    !! 'compact6-hyperviscous'.
    real(dp), optional, intent(in) :: c_nu, nu0_nu
    !! When present, also build the low-pass filter operators for this
    !! direction. Absent means the case does not filter.
    real(dp), optional, intent(in) :: filter_alpha

    integer :: dir, bc_start, bc_end, bc_mp_start, bc_mp_end, n_vert, n_cell, i
    real(dp) :: d

    dir = dirps%dir
    d = mesh%geo%d(dir)

    bc_start = mesh%grid%BCs(dir, 1)
    bc_end = mesh%grid%BCs(dir, 2)

    ! The FFT based Poisson solver needs a Neumann BC for the pressure.
    ! This does not match the Navier-Stokes equations exactly.
    ! It does not affect the quality of the simulation.
    ! The midpoint operators use Neumann where the BC is Dirichlet.
    ! They keep the BC as it is when it is `BC_HALO`.
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
      stretch_correct=mesh%geo%vert_d2s(1:n_vert, dir), &
      c_nu=c_nu, nu0_nu=nu0_nu &
      )
    call backend%alloc_tdsops( &
      dirps%der2nd_sym, n_vert, d, 'second-deriv', der2nd_scheme, &
      bc_start, bc_end, stretch=mesh%geo%vert_ds2(1:n_vert, dir), &
      stretch_correct=mesh%geo%vert_d2s(1:n_vert, dir), &
      sym=.true., c_nu=c_nu, nu0_nu=nu0_nu &
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

    if (present(filter_alpha)) then
      ! Same parity split as the first derivatives: the component along this
      ! direction is odd across a free-slip boundary, the other two even.
      call backend%alloc_tdsops( &
        dirps%lowpass, n_vert, d, 'filter', der1st_scheme, &
        bc_start, bc_end, filter_alpha=filter_alpha &
        )
      call backend%alloc_tdsops( &
        dirps%lowpass_sym, n_vert, d, 'filter', der1st_scheme, &
        bc_start, bc_end, sym=.true., filter_alpha=filter_alpha &
        )
      ! The filter's system is only marginally diagonally dominant as
      ! filter_alpha approaches 0.5, so the distributed solver's truncation
      ! is no longer negligible: at Incompact3d's 0.49 it leaves a 1e-3 error
      ! at zero wavenumber on 64 points. Solve it exactly with Thomas where
      ! the direction is not decomposed; Thomas is serial along the line.
      if (mesh%par%nproc_dir(dir) == 1) then
        dirps%lowpass%prefer_thomas = .true.
        dirps%lowpass_sym%prefer_thomas = .true.
      else if (mesh%par%is_root()) then
        print *, 'WARNING: the low-pass filter is decomposed in direction ', &
          dir, '; it does not exactly preserve the mean and its result &
          &depends on the decomposition.'
      end if
    end if

  end subroutine

end module m_tdsops_setup
