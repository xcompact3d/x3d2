program test_transeq
  use iso_fortran_env, only: stderr => error_unit
  use mpi

  use m_common, only: dp, pi, BC_PERIODIC
  use m_omp_common, only: SZ
  use m_omp_exec_dist, only: exec_dist_transeq_compact
  use m_omp_sendrecv, only: sendrecv_fields
  use m_tdsops, only: tdsops_t

  implicit none

  logical :: allpass = .true.
  real(dp), allocatable, dimension(:, :, :) :: u, v, r_u
  real(dp), allocatable, dimension(:, :, :) :: du, dud, d2u ! intermediate solution arrays
  real(dp), allocatable, dimension(:, :, :) :: &
    du_recv_s, du_recv_e, du_send_s, du_send_e, &
    dud_recv_s, dud_recv_e, dud_send_s, dud_send_e, &
    d2u_recv_s, d2u_recv_e, d2u_send_s, d2u_send_e

  real(dp), allocatable, dimension(:, :, :) :: &
    u_send_s, u_send_e, u_recv_s, u_recv_e, &
    v_send_s, v_send_e, v_recv_s, v_recv_e

  type(tdsops_t) :: der1st, der2nd

  integer :: n, n_block, i, j, k, n_halo, n_iters
  integer :: n_glob
  integer :: nrank, nproc, pprev, pnext
  integer :: ierr

  real(dp) :: dx, dx_per, nu, norm_du, tol = 1d-8

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, nrank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)

  if (nrank == 0) print *, 'Parallel run with', nproc, 'ranks'

  pnext = modulo(nrank - nproc + 1, nproc)
  pprev = modulo(nrank - 1, nproc)

  n_glob = 32*4
  n = n_glob/nproc
  n_block = 32*32/SZ
  n_iters = 1

  nu = 1._dp

  allocate (u(SZ, n, n_block), v(SZ, n, n_block), r_u(SZ, n, n_block))

  ! main input fields
  ! field for storing the result
  ! intermediate solution fields
  allocate (du(SZ, n, n_block))
  allocate (dud(SZ, n, n_block))
  allocate (d2u(SZ, n, n_block))

  dx_per = 2*pi/n_glob
  dx = 2*pi/(n_glob - 1)

  do k = 1, n_block
    do j = 1, n
      do i = 1, SZ
        u(i, j, k) = sin((j - 1 + nrank*n)*dx_per)
        v(i, j, k) = cos((j - 1 + nrank*n)*dx_per)
      end do
    end do
  end do

  n_halo = 4

  ! arrays for exchanging data between ranks
  allocate (u_send_s(SZ, n_halo, n_block))
  allocate (u_send_e(SZ, n_halo, n_block))
  allocate (u_recv_s(SZ, n_halo, n_block))
  allocate (u_recv_e(SZ, n_halo, n_block))
  allocate (v_send_s(SZ, n_halo, n_block))
  allocate (v_send_e(SZ, n_halo, n_block))
  allocate (v_recv_s(SZ, n_halo, n_block))
  allocate (v_recv_e(SZ, n_halo, n_block))

  allocate (du_send_s(SZ, 1, n_block), du_send_e(SZ, 1, n_block))
  allocate (du_recv_s(SZ, 1, n_block), du_recv_e(SZ, 1, n_block))
  allocate (dud_send_s(SZ, 1, n_block), dud_send_e(SZ, 1, n_block))
  allocate (dud_recv_s(SZ, 1, n_block), dud_recv_e(SZ, 1, n_block))
  allocate (d2u_send_s(SZ, 1, n_block), d2u_send_e(SZ, 1, n_block))
  allocate (d2u_recv_s(SZ, 1, n_block), d2u_recv_e(SZ, 1, n_block))

  ! preprocess the operator and coefficient arrays
  der1st = tdsops_t(n, dx_per, operation='first-deriv', scheme='compact6', &
                    bc_start=BC_PERIODIC, bc_end=BC_PERIODIC)
  der2nd = tdsops_t(n, dx_per, operation='second-deriv', scheme='compact6', &
                    bc_start=BC_PERIODIC, bc_end=BC_PERIODIC)

  u_send_s(:, :, :) = u(:, 1:4, :)
  u_send_e(:, :, :) = u(:, n - n_halo + 1:n, :)
  v_send_s(:, :, :) = v(:, 1:4, :)
  v_send_e(:, :, :) = v(:, n - n_halo + 1:n, :)

  ! halo exchange
  call sendrecv_fields(u_recv_s, u_recv_e, &
                       u_send_s, u_send_e, &
                       SZ*4*n_block, nproc, pprev, pnext)

  call sendrecv_fields(v_recv_s, v_recv_e, &
                       v_send_s, v_send_e, &
                       SZ*4*n_block, nproc, pprev, pnext)

  call exec_dist_transeq_compact( &
    r_u, &
    du, dud, d2u, &
    du_send_s, du_send_e, du_recv_s, du_recv_e, &
    dud_send_s, dud_send_e, dud_recv_s, dud_recv_e, &
    d2u_send_s, d2u_send_e, d2u_recv_s, d2u_recv_e, &
    u, u_recv_s, u_recv_e, &
    v, v_recv_s, v_recv_e, &
    der1st, der1st, der2nd, nu, nproc, pprev, pnext, n_block &
    )

  ! check error
  r_u = r_u - (-v*v + 0.5_dp*u*u - nu*u)
  norm_du = norm2(r_u)
  norm_du = norm_du*norm_du/n_glob/n_block/SZ
  call MPI_Allreduce(MPI_IN_PLACE, norm_du, 1, MPI_DOUBLE_PRECISION, &
                     MPI_SUM, MPI_COMM_WORLD, ierr)
  norm_du = sqrt(norm_du)

  if (nrank == 0) print *, 'error norm', norm_du

  if (nrank == 0) then
    if (norm_du > tol) then
      allpass = .false.
      write (stderr, '(a)') 'Check second derivatives... failed'
    else
      write (stderr, '(a)') 'Check second derivatives... passed'
    end if
  end if

  if (allpass) then
    if (nrank == 0) write (stderr, '(a)') 'ALL TESTS PASSED SUCCESSFULLY.'
  else
    error stop 'SOME TESTS FAILED.'
  end if

  call MPI_Finalize(ierr)

end program test_transeq

