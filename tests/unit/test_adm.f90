program test_adm
  use mpi

  use m_adm, only: adm_t
  use m_allocator, only: allocator_t
  use m_base_backend, only: base_backend_t
  use m_common, only: dp, pi, DIR_C, DIR_X, VERT
  use m_field, only: field_t
  use m_mesh, only: mesh_t
#ifdef CUDA
  use m_cuda_allocator, only: cuda_allocator_t
  use m_cuda_backend, only: cuda_backend_t
  use m_cuda_common, only: SZ
#else
  use m_omp_backend, only: omp_backend_t
  use m_omp_common, only: SZ
#endif

  implicit none

  class(base_backend_t), pointer :: backend
  class(allocator_t), pointer :: allocator
  type(allocator_t), target :: host_allocator
  class(mesh_t), allocatable, target :: mesh
#ifdef CUDA
  type(cuda_allocator_t), target :: cuda_allocator
  type(cuda_backend_t), target :: cuda_backend
#else
  type(allocator_t), target :: omp_allocator
  type(omp_backend_t), target :: omp_backend
#endif
  type(adm_t) :: adm
  class(field_t), pointer :: u, v, w, du, dv, dw, host_data
  character(len=8) :: periodic_bcs(2)
  integer :: ierr, nrank
  real(dp) :: expected_speed, expected_thrust, sum_du, sum_dv, sum_dw, tol
  logical :: allpass

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, nrank, ierr)

  periodic_bcs = ['periodic', 'periodic']
  mesh = mesh_t([17, 17, 17], [1, 1, 1], &
                [1._dp, 1._dp, 1._dp], &
                periodic_bcs, periodic_bcs, periodic_bcs)

#ifdef CUDA
  cuda_allocator = cuda_allocator_t(mesh%get_dims(VERT), SZ)
  allocator => cuda_allocator
  host_allocator = allocator_t(mesh%get_dims(VERT), SZ)
  cuda_backend = cuda_backend_t(mesh, allocator)
  backend => cuda_backend
#else
  omp_allocator = allocator_t(mesh%get_dims(VERT), SZ)
  allocator => omp_allocator
  host_allocator = allocator_t(mesh%get_dims(VERT), SZ)
  omp_backend = omp_backend_t(mesh, allocator)
  backend => omp_backend
#endif

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

  call backend%vecadd(0._dp, adm%disc(1)%gamma_disc, 0._dp, du)
  call backend%vecadd(0._dp, adm%disc(1)%gamma_disc, 0._dp, dv)
  call backend%vecadd(0._dp, adm%disc(1)%gamma_disc, 0._dp, dw)
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

  ! A repeated RK stage at the same physical time reuses the force.
  host_data%data = 20._dp
  call backend%set_field_data(u, host_data%data)
  call backend%vecadd(0._dp, adm%disc(1)%gamma_disc, 0._dp, du)
  call backend%vecadd(0._dp, adm%disc(1)%gamma_disc, 0._dp, dv)
  call backend%vecadd(0._dp, adm%disc(1)%gamma_disc, 0._dp, dw)
  call adm%update(0.01_dp, 0.01_dp)
  call adm%project_forces(du, dv, dw, u, v, w)
  sum_du = backend%field_volume_integral(du)*adm%cell_vol
  allpass = allpass &
            .and. close_enough(adm%disc(1)%U_disc, 10._dp, tol) &
            .and. close_enough(sum_du, -expected_thrust/adm%rho_air, tol)

  ! The next physical timestep exercises the T_relax filter.
  call backend%vecadd(0._dp, adm%disc(1)%gamma_disc, 0._dp, du)
  call backend%vecadd(0._dp, adm%disc(1)%gamma_disc, 0._dp, dv)
  call backend%vecadd(0._dp, adm%disc(1)%gamma_disc, 0._dp, dw)
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
            .and. close_enough(adm%disc(1)%thrust, expected_thrust, tol) &
            .and. close_enough(adm%disc(1)%power, &
                               expected_speed*expected_thrust, tol) &
            .and. close_enough(sum_du, -expected_thrust/adm%rho_air, tol)

  if (nrank == 0) then
    open (unit=10, file='test_adm.ad', status='old')
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
