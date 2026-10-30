program test_alm
  !! Tests the actuator line model pieces: airfoil polar interpolation, the
  !! rotation helper, the NREL 5 MW controller, the backend point sampling
  !! (trilinear, across ranks and periodic boundaries) and Gaussian source
  !! projection, and a full ALM step whose projected momentum source must
  !! balance the element forces.
  use m_mpi, only: MPI_COMM_WORLD, MPI_IN_PLACE, MPI_SUM, MPI_Init, &
                   MPI_Finalize, MPI_Comm_rank, MPI_Comm_size, MPI_Barrier, &
                   MPI_Allreduce

  use m_alm, only: alm_t
  use m_alm_airfoil, only: airfoil_t, read_airfoil, eval_static_coeffs
  use m_alm_controller, only: controller_t, init_controller, &
                              operate_controller
  use m_alm_element, only: quat_rot
  use m_allocator, only: allocator_t
  use m_backend_runtime, only: backend_runtime_t
  use m_base_backend, only: base_backend_t
  use m_common, only: dp, pi, MPI_X3D2_DP, DIR_C, DIR_X, VERT
  use m_field, only: field_t
  use m_mesh, only: mesh_t

  implicit none

  class(base_backend_t), pointer :: backend
  class(allocator_t), pointer :: allocator
  class(allocator_t), pointer :: host_allocator
  type(mesh_t), target :: mesh
  type(backend_runtime_t), target :: runtime
  type(airfoil_t) :: airfoil
  type(alm_t) :: alm
  type(controller_t) :: ctrl
  class(field_t), pointer :: u, v, w, du, dv, dw, host_data
  character(len=9) :: dirichlet_bcs(2), periodic_bcs(2)
  integer :: ierr, nrank, nproc, i, j, k, dims(3), gdims(3), s, t, ib, ie
  real(dp) :: cl, cd, cm, tol, proj_tol, coords(3), d(3), cell_vol, eps
  real(dp) :: expected
  real(dp) :: points(3, 3), vals(3, 3), coeffs(3, 1), src(3, 1), f_tot(3)
  real(dp) :: sums(3)
  logical :: allpass, ok

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, nrank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
  allpass = .true.
  ! Scaled with the working precision, so the test also holds for a
  ! single-precision build
  tol = 1000._dp*epsilon(1._dp)
  proj_tol = max(1.e-8_dp, 100._dp*epsilon(1._dp))
  ok = .true.

  ! Polars: two Re tables, CL linear in alpha, CD and CM constant per table
  if (nrank == 0) then
    open (unit=10, file='test_alm.air', status='replace', action='write')
    write (10, '(A)') 'Title: test foil'
    write (10, '(A)') 'Thickness to Chord Ratio: 0.2'
    write (10, '(A)') 'Zero Lift AOA (deg): 0.0'
    write (10, '(A)') 'Reverse Camber Direction: 0'
    write (10, '(A)') ''
    write (10, '(A)') 'Reynolds Number: 1e4'
    write (10, '(A)') 'AOA (deg) CL CD Cm25'
    write (10, '(A)') '-180.0 -2.0 0.1 0.0'
    write (10, '(A)') '0.0 0.0 0.1 0.0'
    write (10, '(A)') '180.0 2.0 0.1 0.0'
    write (10, '(A)') ''
    write (10, '(A)') 'Reynolds Number: 2e4'
    write (10, '(A)') 'AOA (deg) CL CD Cm25'
    write (10, '(A)') '-180.0 -4.0 0.3 0.2'
    write (10, '(A)') '0.0 0.0 0.3 0.2'
    write (10, '(A)') '180.0 4.0 0.3 0.2'
    close (10)
  end if
  call MPI_Barrier(MPI_COMM_WORLD, ierr)

  call read_airfoil(airfoil, 'test_alm.air')
  ok = ok .and. airfoil%n_re == 2 &
       .and. close_enough(airfoil%tc, 0.2_dp, tol)
  ! Midway in Re, a quarter of the way up in alpha
  call eval_static_coeffs(airfoil, 1.5e4_dp, 45._dp, cl, cd, cm)
  ok = ok .and. close_enough(cl, 0.75_dp, tol) &
       .and. close_enough(cd, 0.2_dp, tol) &
       .and. close_enough(cm, 0.1_dp, tol)
  ! Outside the tabulated Re range the nearest table is used
  call eval_static_coeffs(airfoil, 1.e3_dp, -90._dp, cl, cd, cm)
  ok = ok .and. close_enough(cl, -1._dp, tol)
  call eval_static_coeffs(airfoil, 1.e6_dp, 90._dp, cl, cd, cm)
  ok = ok .and. close_enough(cl, 2._dp, tol) &
       .and. close_enough(cd, 0.3_dp, tol)

  ! Right-handed rotation: z rotated a quarter turn about x goes to -y
  ok = ok .and. &
       all(abs(quat_rot([1._dp, 2._dp, 3._dp], pi/2._dp, &
                        [2._dp, 0._dp, 0._dp], [1._dp, 2._dp, 2._dp]) &
               - [1._dp, 1._dp, 2._dp]) < tol)
  call report('airfoil and rotation')

  ! NREL 5 MW controller at dt = vs_dt = pc_dt. In region 2 the first call
  ! already commands the torque law (k*w^2), held at constant speed
  call init_nrel_controller(ctrl)
  do s = 1, 3
    call operate_controller(ctrl, s*ctrl%vs_dt, 3, 100._dp/ctrl%gear_box_ratio)
    ctrl%status = ctrl%status + 1
    ok = ok .and. close_enough(ctrl%gen_trq, ctrl%vs_rgn2_k*100._dp**2, tol) &
         .and. all(ctrl%pit_com == 0._dp)
    if (.not. ok .and. nrank == 0) print *, 'step', s, 'torque', ctrl%gen_trq
  end do
  ! A speed jump is then rate limited over one interval, not two
  expected = ctrl%gen_trq + ctrl%vs_max_rat*ctrl%vs_dt
  call operate_controller(ctrl, 4*ctrl%vs_dt, 3, 130._dp/ctrl%gear_box_ratio)
  ok = ok .and. close_enough(ctrl%gen_trq, expected, tol)
  if (.not. ok .and. nrank == 0) print *, 'rate limited torque', ctrl%gen_trq

  ! Starting above rated speed, the first call commands region 3 torque and
  ! pitches at the maximum rate
  call init_nrel_controller(ctrl)
  call operate_controller(ctrl, ctrl%vs_dt, 3, 130._dp/ctrl%gear_box_ratio)
  ok = ok .and. close_enough(ctrl%gen_trq, ctrl%vs_rt_pwr/130._dp, tol) &
       .and. all(abs(ctrl%pit_com - ctrl%pc_max_rat*ctrl%pc_dt) &
                 < tol*ctrl%pc_max_rat*ctrl%pc_dt)
  if (.not. ok .and. nrank == 0) print *, 'region 3 first call', &
    ctrl%gen_trq, ctrl%pit_com
  call report('controller')

  ! Inflow/outflow in x, periodic y and z, decomposed in z so sampling and
  ! sources straddle rank boundaries
  dirichlet_bcs = ['dirichlet', 'dirichlet']
  periodic_bcs = ['periodic ', 'periodic ']
  mesh = mesh_t([33, 32, 32], [1, 1, nproc], [1._dp, 1._dp, 1._dp], &
                dirichlet_bcs, periodic_bcs, periodic_bcs)
  call runtime%init(mesh, separate_host_allocator=.true.)
  backend => runtime%backend
  allocator => runtime%allocator
  host_allocator => runtime%host_allocator

  u => allocator%get_block(DIR_X, VERT)
  v => allocator%get_block(DIR_X, VERT)
  w => allocator%get_block(DIR_X, VERT)
  du => allocator%get_block(DIR_X, VERT)
  dv => allocator%get_block(DIR_X, VERT)
  dw => allocator%get_block(DIR_X, VERT)
  host_data => host_allocator%get_block(DIR_C, VERT)
  dims = mesh%get_dims(VERT)
  gdims = mesh%get_global_dims(VERT)
  d = mesh%geo%d
  cell_vol = product(d)

  ! Trilinear sampling is exact for a trilinear field; w is periodic in y
  do k = 1, dims(3)
    do j = 1, dims(2)
      do i = 1, dims(1)
        coords = mesh%get_coordinates(i, j, k)
        host_data%data(i, j, k) = 1._dp + 2._dp*coords(1) &
                                  - 3._dp*coords(2) + 0.5_dp*coords(3) &
                                  + coords(1)*coords(2)*coords(3)
      end do
    end do
  end do
  call backend%set_field_data(u, host_data%data)
  host_data%data = 2._dp
  call backend%set_field_data(v, host_data%data)
  do k = 1, dims(3)
    do j = 1, dims(2)
      do i = 1, dims(1)
        coords = mesh%get_coordinates(i, j, k)
        host_data%data(i, j, k) = cos(2._dp*pi*coords(2))
      end do
    end do
  end do
  call backend%set_field_data(w, host_data%data)

  points(:, 1) = [0.31_dp, 0.47_dp, 0.52_dp]
  points(:, 2) = [0.5_dp, 0.5_dp, 0.5_dp]           ! on a node
  points(:, 3) = [0.2_dp, 1._dp - 0.25_dp*d(2), 0.6_dp] ! past the last y node
  call backend%field_sample_points(vals, points, u, v, w)
  call MPI_Allreduce(MPI_IN_PLACE, vals, size(vals), MPI_X3D2_DP, MPI_SUM, &
                     MPI_COMM_WORLD, ierr)
  do s = 1, 2
    expected = 1._dp + 2._dp*points(1, s) - 3._dp*points(2, s) &
               + 0.5_dp*points(3, s) + product(points(:, s))
    ok = ok .and. close_enough(vals(1, s), expected, tol) &
         .and. close_enough(vals(2, s), 2._dp, tol)
  end do
  ! Between the last y node and the first (wrapped) one
  expected = 0.25_dp*cos(2._dp*pi*(1._dp - d(2))) + 0.75_dp
  ok = ok .and. close_enough(vals(3, 3), expected, tol)
  call report('sampling')

  ! Gaussian source: the kernel sums to eps^3 pi^(3/2) over the grid
  host_data%data = 0._dp
  call backend%set_field_data(du, host_data%data)
  call backend%set_field_data(dv, host_data%data)
  call backend%set_field_data(dw, host_data%data)
  eps = 2._dp*d(1)
  src(:, 1) = [0.5_dp, 0.45_dp, 0.55_dp]
  coeffs(:, 1) = [1._dp, -2._dp, 3._dp]
  call backend%field_add_gaussian_sources(du, dv, dw, src, coeffs, eps, &
                                          6._dp*eps)
  sums = [backend%field_volume_integral(du), &
          backend%field_volume_integral(dv), &
          backend%field_volume_integral(dw)]*cell_vol/(eps**3*pi**1.5_dp)
  ok = ok .and. all(abs(sums - coeffs(:, 1)) < proj_tol*abs(coeffs(:, 1)))
  if (.not. ok .and. nrank == 0) print *, 'projected integrals', sums

  ! A source straddling the periodic y and z boundaries keeps its whole
  ! kernel: the part beyond y = 0 and z = 0 wraps to the far side, in z
  ! onto another rank
  host_data%data = 0._dp
  call backend%set_field_data(du, host_data%data)
  call backend%set_field_data(dv, host_data%data)
  call backend%set_field_data(dw, host_data%data)
  src(:, 1) = [0.5_dp, 0.01_dp, 0.02_dp]
  call backend%field_add_gaussian_sources(du, dv, dw, src, coeffs, eps, &
                                          6._dp*eps)
  sums = [backend%field_volume_integral(du), &
          backend%field_volume_integral(dv), &
          backend%field_volume_integral(dw)]*cell_vol/(eps**3*pi**1.5_dp)
  ok = ok .and. all(abs(sums - coeffs(:, 1)) < proj_tol*abs(coeffs(:, 1)))
  if (.not. ok .and. nrank == 0) print *, 'wrapped integrals', sums

  ! Next to the non-periodic x inflow the part beyond x = 0 is cut, as in
  ! Incompact3d: only the in-domain share of the x kernel remains
  host_data%data = 0._dp
  call backend%set_field_data(du, host_data%data)
  call backend%set_field_data(dv, host_data%data)
  call backend%set_field_data(dw, host_data%data)
  src(:, 1) = [0.02_dp, 0.45_dp, 0.55_dp]
  call backend%field_add_gaussian_sources(du, dv, dw, src, coeffs, eps, &
                                          6._dp*eps)
  expected = 0._dp
  do i = 0, gdims(1) - 1
    expected = expected + exp(-((i*d(1) - src(1, 1))/eps)**2)
  end do
  expected = expected*d(1)/(eps*sqrt(pi))
  sums = [backend%field_volume_integral(du), &
          backend%field_volume_integral(dv), &
          backend%field_volume_integral(dw)]*cell_vol/(eps**3*pi**1.5_dp)
  ok = ok .and. expected < 0.9_dp &
       .and. all(abs(sums - expected*coeffs(:, 1)) &
                 < proj_tol*abs(coeffs(:, 1)))
  if (.not. ok .and. nrank == 0) print *, 'cut integrals', sums, expected
  call report('projection')

  ! Full ALM step: one non-rotating 2-bladed rotor in uniform flow. The
  ! projected momentum source must integrate to -sum(F)/rho.
  if (nrank == 0) then
    open (unit=10, file='test_alm.al', status='replace', action='write')
    write (10, '(A)') 'R : 0.2'
    write (10, '(A)') 'Spanwise : 0.0 0.0 1.0'
    write (10, '(A)') 'NStations : 6'
    do s = 0, 5
      write (10, '(4(ES24.16,1X))') 0.2_dp + 0.16_dp*s, 0.2_dp, &
        10._dp - s, 0.2_dp
    end do
    close (10)
    open (unit=10, file='test_alm.turb', status='replace', action='write')
    write (10, '(A)') '&TurbineSpecs'
    write (10, '(A)') "name = 'test_rotor'"
    write (10, '(A)') 'origin = 0.5, 0.5, 0.5'
    write (10, '(A)') 'numblades = 2'
    write (10, '(A)') "blade_geom = 'test_alm.al'"
    write (10, '(A)') 'numfoil = 1'
    write (10, '(A)') "afname = 'test_alm.air'"
    write (10, '(A)') 'TypeFlag = 1'
    write (10, '(A)') 'OperFlag = 1'
    write (10, '(A)') 'uref = 1.0'
    write (10, '(A)') 'tsr = 2.0'
    write (10, '(A)') 'RotFlag = 1'
    write (10, '(A)') '/'
    close (10)
  end if
  call MPI_Barrier(MPI_COMM_WORLD, ierr)

  host_data%data = 1._dp
  call backend%set_field_data(u, host_data%data)
  host_data%data = 0._dp
  call backend%set_field_data(v, host_data%data)
  call backend%set_field_data(w, host_data%data)
  call backend%set_field_data(du, host_data%data)
  call backend%set_field_data(dv, host_data%data)
  call backend%set_field_data(dw, host_data%data)

  alm%turbine_files = ['test_alm.turb']
  alm%rho_air = 1.2_dp
  alm%eps_factor = 1.5_dp
  alm%nu = 1.e-5_dp
  call alm%init(backend, mesh, host_allocator, 0.01_dp)
  call alm%setup_output(nrank == 0, .false.)
  call alm%update(0.01_dp, 0.01_dp)
  call alm%project_forces(du, dv, dw, u, v, w)

  f_tot = 0._dp
  do t = 1, alm%n_turb
    do ib = 1, alm%turbine(t)%n_blades
      do ie = 1, alm%turbine(t)%blade(ib)%n_elem
        f_tot = f_tot + alm%turbine(t)%blade(ib)%ef(:, ie)
      end do
    end do
  end do
  sums = [backend%field_volume_integral(du), &
          backend%field_volume_integral(dv), &
          backend%field_volume_integral(dw)]*cell_vol
  ! The kernels are cut off at 10 eps but the rotor lies 6 eps from the
  ! nearest boundary, so a small share is lost there
  ok = ok &
       .and. all(abs(sums + f_tot/alm%rho_air) &
                 < 1.e-3_dp*norm2(f_tot)/alm%rho_air) &
       .and. norm2(f_tot) > 0._dp &
       .and. close_enough(alm%turbine(1)%angular_vel, 10._dp, tol) &
       .and. close_enough(alm%turbine(1)%azim_angle, 0.1_dp, tol)

  ! A repeated substage adds the same source again without recomputing it
  call alm%update(0.01_dp, 0.01_dp)
  call alm%project_forces(du, dv, dw, u, v, w)
  sums = [backend%field_volume_integral(du), &
          backend%field_volume_integral(dv), &
          backend%field_volume_integral(dw)]*cell_vol
  ok = ok &
       .and. all(abs(sums + 2._dp*f_tot/alm%rho_air) &
                 < 2.e-3_dp*norm2(f_tot)/alm%rho_air) &
       .and. close_enough(alm%turbine(1)%azim_angle, 0.1_dp, tol)
  if (.not. ok .and. nrank == 0) print *, 'projected and expected', &
    sums, -2._dp*f_tot/alm%rho_air
  call report('ALM step')

  call alm%write_output(1, nrank == 0)
  call alm%finalise()

  if (nrank == 0) then
    call delete_file('test_alm.air')
    call delete_file('test_alm.al')
    call delete_file('test_alm.turb')
    call delete_file('test_rotor.perf')
    call delete_file('test_rotor_blade_1.load')
    call delete_file('test_rotor_blade_2.load')
    if (allpass) then
      print *, 'PASS'
    else
      error stop 'ALM test failed'
    end if
  end if

  call allocator%release_block(u)
  call allocator%release_block(v)
  call allocator%release_block(w)
  call allocator%release_block(du)
  call allocator%release_block(dv)
  call allocator%release_block(dw)
  call host_allocator%release_block(host_data)

  call MPI_Finalize(ierr)

contains

  subroutine report(section)
    !! Fold the checks of one section into allpass, naming it on failure.
    character(len=*), intent(in) :: section

    if (.not. ok .and. nrank == 0) print *, trim(section)//' checks failed'
    allpass = allpass .and. ok
    ok = .true.
  end subroutine report

  subroutine init_nrel_controller(ctrl)
    !! A fresh controller with the NREL 5 MW parameters of the
    !! Incompact3d example.
    type(controller_t), intent(out) :: ctrl

    call init_controller(ctrl, generator_inertia=534.116_dp, &
                         gb_ratio=97._dp, rated_gen_speed=121.6805_dp, &
                         rated_limit_gen_torque=15000._dp, &
                         cut_in_gen_speed=70.16224_dp, &
                         region2_start_gen_speed=91.21091_dp, kgen=2.33_dp, &
                         rated_power=5296610._dp, maximum_torque=47402.91_dp)
  end subroutine init_nrel_controller

  logical function close_enough(actual, expected, relative_tol)
    real(dp), intent(in) :: actual, expected, relative_tol

    close_enough = abs(actual - expected) &
                   <= relative_tol*max(1._dp, abs(expected))
  end function close_enough

  subroutine delete_file(filename)
    character(len=*), intent(in) :: filename
    integer :: unit, ios

    open (newunit=unit, file=filename, status='old', iostat=ios)
    if (ios == 0) close (unit, status='delete')
  end subroutine delete_file

end program test_alm
