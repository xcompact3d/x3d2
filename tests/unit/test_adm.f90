program test_adm
  !! Tests ADM velocity sampling, filtering, force projection, running means,
  !! and scalar-output columns independently of the checkpoint backend, then
  !! checks that two discs sharing the kernel fields each sample only their
  !! own inflow.
  use m_mpi, only: MPI_COMM_WORLD, MPI_Init, MPI_Finalize, MPI_Comm_rank, &
                   MPI_Barrier

  use m_adm, only: adm_t
  use m_allocator, only: allocator_t
  use m_backend_runtime, only: backend_runtime_t
  use m_base_backend, only: base_backend_t
  use m_common, only: dp, pi, DIR_C, DIR_X, VERT
  use m_field, only: field_t
  use m_mesh, only: mesh_t

  implicit none

  class(base_backend_t), pointer :: backend
  class(allocator_t), pointer :: allocator
  class(allocator_t), pointer :: host_allocator
  type(mesh_t), target :: mesh
  type(backend_runtime_t), target :: runtime
  type(adm_t) :: adm, adm_pair
  class(field_t), pointer :: u, v, w, du, dv, dw, host_data
  character(len=8) :: periodic_bcs(2)
  character(len=512) :: output_header
  integer :: ierr, nrank, i, j, k, dims(3)
  real(dp) :: expected_speed, expected_thrust, sum_du, sum_dv, sum_dw, tol
  real(dp) :: thrust_coeff, coords(3)
  logical :: allpass

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, nrank, ierr)

  periodic_bcs = ['periodic', 'periodic']
  mesh = mesh_t([17, 17, 17], [1, 1, 1], &
                [1._dp, 1._dp, 1._dp], &
                periodic_bcs, periodic_bcs, periodic_bcs)

  call runtime%init(mesh, separate_host_allocator=.true.)
  backend => runtime%backend
  allocator => runtime%allocator
  host_allocator => runtime%host_allocator

  if (nrank == 0) then
    open (unit=10, file='test_adm.ad', status='replace', action='write')
    write (10, '(A)') 'CoR_x CoR_y CoR_z Yaw Tilt D C_T alpha'
    write (10, '(*(ES24.16,1X))') 0.5_dp, 0.5_dp, 0.5_dp, &
      0._dp, 0._dp, 0.5_dp, 0.75_dp, 0.25_dp
    close (10)
  end if
  call MPI_Barrier(MPI_COMM_WORLD, ierr)

  adm%coords_file = 'test_adm.ad'
  adm%rho_air = 1.225_dp
  adm%T_relax = 1._dp
  call adm%init(backend, mesh, host_allocator, 0.01_dp)
  call adm%setup_output(nrank == 0, .false.)

  u => allocator%get_block(DIR_X, VERT)
  v => allocator%get_block(DIR_X, VERT)
  w => allocator%get_block(DIR_X, VERT)
  du => allocator%get_block(DIR_X, VERT)
  dv => allocator%get_block(DIR_X, VERT)
  dw => allocator%get_block(DIR_X, VERT)
  host_data => host_allocator%get_block(DIR_C, VERT)

  host_data%data = 10._dp
  call backend%set_field_data(u, host_data%data)
  host_data%data = 0._dp
  call backend%set_field_data(v, host_data%data)
  call backend%set_field_data(w, host_data%data)

  call backend%vecadd(0._dp, adm%gamma, 0._dp, du)
  call backend%vecadd(0._dp, adm%gamma, 0._dp, dv)
  call backend%vecadd(0._dp, adm%gamma, 0._dp, dw)
  call adm%update(0.01_dp, 0.01_dp)
  call adm%project_forces(du, dv, dw, u, v, w)

  expected_thrust = 0.5_dp*adm%rho_air &
                    *(0.75_dp/(1._dp - 0.25_dp)**2) &
                    *10._dp**2*pi*0.5_dp**2/4._dp
  sum_du = backend%field_volume_integral(du)*adm%cell_vol
  sum_dv = backend%field_volume_integral(dv)*adm%cell_vol
  sum_dw = backend%field_volume_integral(dw)*adm%cell_vol
  tol = 1000._dp*epsilon(1._dp)

  allpass = close_enough(adm%disc(1)%U_disc, 10._dp, tol) &
            .and. close_enough(adm%disc(1)%U_disc_filt, 10._dp, tol) &
            .and. close_enough(adm%disc(1)%thrust, expected_thrust, tol) &
           .and. close_enough(adm%disc(1)%power, 10._dp*expected_thrust, tol) &
            .and. close_enough(sum_du, -expected_thrust/adm%rho_air, tol) &
            .and. close_enough(sum_dv, 0._dp, tol) &
            .and. close_enough(sum_dw, 0._dp, tol)

  ! A repeated RK stage at the same physical time: the filter is not advanced,
  ! but U_disc and the thrust track the current substage velocity.
  host_data%data = 20._dp
  call backend%set_field_data(u, host_data%data)
  call backend%vecadd(0._dp, adm%gamma, 0._dp, du)
  call backend%vecadd(0._dp, adm%gamma, 0._dp, dv)
  call backend%vecadd(0._dp, adm%gamma, 0._dp, dw)
  call adm%update(0.01_dp, 0.01_dp)
  call adm%project_forces(du, dv, dw, u, v, w)
  expected_thrust = 0.5_dp*adm%rho_air &
                    *(0.75_dp/(1._dp - 0.25_dp)**2) &
                    *20._dp**2*pi*0.5_dp**2/4._dp
  sum_du = backend%field_volume_integral(du)*adm%cell_vol
  allpass = allpass &
            .and. close_enough(adm%disc(1)%U_disc, 20._dp, tol) &
            .and. close_enough(sum_du, -expected_thrust/adm%rho_air, tol)

  ! The next physical timestep exercises the T_relax filter.
  call backend%vecadd(0._dp, adm%gamma, 0._dp, du)
  call backend%vecadd(0._dp, adm%gamma, 0._dp, dv)
  call backend%vecadd(0._dp, adm%gamma, 0._dp, dw)
  call adm%update(0.02_dp, 0.01_dp)
  call adm%project_forces(du, dv, dw, u, v, w)
  expected_speed = (0.01_dp/1.01_dp)*20._dp + (1._dp/1.01_dp)*10._dp
  expected_thrust = 0.5_dp*adm%rho_air &
                    *(0.75_dp/(1._dp - 0.25_dp)**2) &
                    *expected_speed**2*pi*0.5_dp**2/4._dp
  sum_du = backend%field_volume_integral(du)*adm%cell_vol
  allpass = allpass &
            .and. close_enough(adm%disc(1)%U_disc, 20._dp, tol) &
            .and. close_enough(adm%disc(1)%U_disc_filt, expected_speed, tol) &
            .and. adm%sample_count == 2 &
            .and. close_enough(adm%disc(1)%U_disc_mean, &
                               0.5_dp*(10._dp + expected_speed), tol) &
            .and. close_enough(adm%disc(1)%thrust, expected_thrust, tol) &
            .and. close_enough(adm%disc(1)%power, &
                               expected_speed*expected_thrust, tol) &
            .and. close_enough(sum_du, -expected_thrust/adm%rho_air, tol)

  call adm%write_output(2, nrank == 0)
  call adm%finalise()

  ! Two discs, one in each half of the domain in x, with u = 10 upstream of
  ! x = 0.5 and u = 20 downstream.
  if (nrank == 0) then
    open (unit=10, file='test_adm_pair.ad', status='replace', action='write')
    write (10, '(A)') 'CoR_x CoR_y CoR_z Yaw Tilt D C_T alpha'
    write (10, '(*(ES24.16,1X))') 0.25_dp, 0.5_dp, 0.5_dp, &
      0._dp, 0._dp, 0.4_dp, 0.75_dp, 0.25_dp
    write (10, '(*(ES24.16,1X))') 0.75_dp, 0.5_dp, 0.5_dp, &
      0._dp, 0._dp, 0.4_dp, 0.75_dp, 0.25_dp
    close (10)
  end if
  call MPI_Barrier(MPI_COMM_WORLD, ierr)

  adm_pair%coords_file = 'test_adm_pair.ad'
  adm_pair%rho_air = 1.225_dp
  call adm_pair%init(backend, mesh, host_allocator, 0.01_dp)

  dims = mesh%get_dims(VERT)
  host_data%data = 0._dp
  do k = 1, dims(3)
    do j = 1, dims(2)
      do i = 1, dims(1)
        coords = mesh%get_coordinates(i, j, k)
        if (coords(1) < 0.5_dp) then
          host_data%data(i, j, k) = 10._dp
        else
          host_data%data(i, j, k) = 20._dp
        end if
      end do
    end do
  end do
  call backend%set_field_data(u, host_data%data)
  host_data%data = 0._dp
  call backend%set_field_data(v, host_data%data)
  call backend%set_field_data(w, host_data%data)
  call backend%set_field_data(du, host_data%data)
  call backend%set_field_data(dv, host_data%data)
  call backend%set_field_data(dw, host_data%data)

  call adm_pair%update(0.01_dp, 0.01_dp)
  call adm_pair%project_forces(du, dv, dw, u, v, w)

  thrust_coeff = 0.5_dp*adm_pair%rho_air*(0.75_dp/(1._dp - 0.25_dp)**2) &
                 *pi*0.4_dp**2/4._dp
  expected_thrust = thrust_coeff*(10._dp**2 + 20._dp**2)
  sum_du = backend%field_volume_integral(du)*adm_pair%cell_vol
  sum_dv = backend%field_volume_integral(dv)*adm_pair%cell_vol
  allpass = allpass &
            .and. close_enough(adm_pair%disc(1)%U_disc, 10._dp, tol) &
            .and. close_enough(adm_pair%disc(2)%U_disc, 20._dp, tol) &
            .and. close_enough(sum_du, -expected_thrust/adm_pair%rho_air, tol) &
            .and. close_enough(sum_dv, 0._dp, tol)
  call adm_pair%finalise()

  if (nrank == 0) then
    open (unit=11, file='disc1.adm', status='old', action='read')
    read (11, '(A)') output_header
    close (11, status='delete')
    allpass = allpass &
              .and. index(output_header, 'U_disc_filt') > 0 &
              .and. index(output_header, 'thrust_mean') > 0
    open (unit=10, file='test_adm.ad', status='old')
    close (10, status='delete')
    open (unit=10, file='test_adm_pair.ad', status='old')
    close (10, status='delete')
    if (allpass) then
      print *, 'PASS'
    else
      error stop 'ADM test failed'
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

  logical function close_enough(actual, expected, relative_tol)
    real(dp), intent(in) :: actual, expected, relative_tol
    close_enough = abs(actual - expected) &
                   <= relative_tol*max(1._dp, abs(expected))
  end function close_enough

end program test_adm
