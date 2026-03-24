program test_cuda_transeq
  use iso_fortran_env, only: stderr => error_unit
  use cudafor
  use mpi

  use m_common, only: dp, pi, MPI_X3D2_DP, BC_PERIODIC
  use m_cuda_common, only: SZ
  use m_cuda_exec_dist, only: exec_dist_transeq_3fused
  use m_cuda_sendrecv, only: sendrecv_fields
  use m_cuda_tdsops, only: cuda_tdsops_t

  implicit none

  logical :: allpass = .true.
  real(dp), allocatable, dimension(:, :, :) :: u, v, r_u
  real(dp), device, allocatable, dimension(:, :, :) :: &
    u_dev, v_dev, r_u_dev, & ! main fields u, v and result r_u
    dud_dev, d2u_dev ! intermediate solution arrays
  real(dp), device, allocatable, dimension(:, :, :) :: &
    du_recv_s_dev, du_recv_e_dev, du_send_s_dev, du_send_e_dev, &
    dud_recv_s_dev, dud_recv_e_dev, dud_send_s_dev, dud_send_e_dev, &
    d2u_recv_s_dev, d2u_recv_e_dev, d2u_send_s_dev, d2u_send_e_dev

  real(dp), device, allocatable, dimension(:, :, :) :: &
    u_send_s_dev, u_send_e_dev, u_recv_s_dev, u_recv_e_dev, &
    v_send_s_dev, v_send_e_dev, v_recv_s_dev, v_recv_e_dev

  type(cuda_tdsops_t) :: der1st, der2nd

  integer :: n, n_block, i, j, k, n_halo
  integer :: n_glob
  integer :: nrank, nproc, pprev, pnext
  integer :: ierr, ndevs, devnum

  type(dim3) :: blocks, threads
  real(dp), parameter :: residual_tol = 1.0e-8_dp
  real(dp) :: dx, dx_per, nu, norm_du

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, nrank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)

  if (nrank == 0) print *, 'Parallel run with', nproc, 'ranks'

  ierr = cudaGetDeviceCount(ndevs)
  ierr = cudaSetDevice(mod(nrank, ndevs)) ! round-robin
  ierr = cudaGetDevice(devnum)

  !print*, 'I am rank', nrank, 'I am running on device', devnum
  pnext = modulo(nrank - nproc + 1, nproc)
  pprev = modulo(nrank - 1, nproc)

  n_glob = 512
  n = n_glob/nproc
  n_block = 512*512/SZ

  nu = 1._dp

  allocate (u(SZ, n, n_block), v(SZ, n, n_block), r_u(SZ, n, n_block))

  ! main input fields
  allocate (u_dev(SZ, n, n_block), v_dev(SZ, n, n_block))
  ! field for storing the result
  allocate (r_u_dev(SZ, n, n_block))
  ! intermediate solution fields
  allocate (dud_dev(SZ, n, n_block))
  allocate (d2u_dev(SZ, n, n_block))

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

  ! move data to device
  u_dev = u
  v_dev = v

  n_halo = 4

  ! arrays for exchanging data between ranks
  allocate (u_send_s_dev(SZ, n_halo, n_block))
  allocate (u_send_e_dev(SZ, n_halo, n_block))
  allocate (u_recv_s_dev(SZ, n_halo, n_block))
  allocate (u_recv_e_dev(SZ, n_halo, n_block))
  allocate (v_send_s_dev(SZ, n_halo, n_block))
  allocate (v_send_e_dev(SZ, n_halo, n_block))
  allocate (v_recv_s_dev(SZ, n_halo, n_block))
  allocate (v_recv_e_dev(SZ, n_halo, n_block))

  allocate (du_send_s_dev(SZ, 1, n_block), du_send_e_dev(SZ, 1, n_block))
  allocate (du_recv_s_dev(SZ, 1, n_block), du_recv_e_dev(SZ, 1, n_block))
  allocate (dud_send_s_dev(SZ, 1, n_block), dud_send_e_dev(SZ, 1, n_block))
  allocate (dud_recv_s_dev(SZ, 1, n_block), dud_recv_e_dev(SZ, 1, n_block))
  allocate (d2u_send_s_dev(SZ, 1, n_block), d2u_send_e_dev(SZ, 1, n_block))
  allocate (d2u_recv_s_dev(SZ, 1, n_block), d2u_recv_e_dev(SZ, 1, n_block))

  ! preprocess the operator and coefficient arrays
  der1st = cuda_tdsops_t(n, dx_per, operation='first-deriv', &
                         scheme='compact6', &
                         bc_start=BC_PERIODIC, bc_end=BC_PERIODIC)
  der2nd = cuda_tdsops_t(n, dx_per, operation='second-deriv', &
                         scheme='compact6', &
                         bc_start=BC_PERIODIC, bc_end=BC_PERIODIC)

  blocks = dim3(n_block, 1, 1)
  threads = dim3(SZ, 1, 1)

  call run_kernel()

  ! check error
  r_u = r_u_dev
  r_u = r_u - (-v*v + 0.5_dp*u*u - nu*u)
  norm_du = norm2(r_u)
  norm_du = norm_du*norm_du/real(n_glob*n_block*SZ, dp)
  call MPI_Allreduce(MPI_IN_PLACE, norm_du, 1, MPI_X3D2_DP, &
                     MPI_SUM, MPI_COMM_WORLD, ierr)
  norm_du = sqrt(norm_du)

  if (nrank == 0) print *, 'error norm', norm_du

  if (nrank == 0) then
    if (norm_du > residual_tol) then
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

contains

  subroutine run_kernel()
    u_send_s_dev(:, :, :) = u_dev(:, 1:n_halo, :)
    u_send_e_dev(:, :, :) = u_dev(:, n - n_halo + 1:n, :)
    v_send_s_dev(:, :, :) = v_dev(:, 1:n_halo, :)
    v_send_e_dev(:, :, :) = v_dev(:, n - n_halo + 1:n, :)

    call sendrecv_fields(u_recv_s_dev, u_recv_e_dev, &
                         u_send_s_dev, u_send_e_dev, &
                         SZ*n_halo*n_block, nproc, pprev, pnext)

    call sendrecv_fields(v_recv_s_dev, v_recv_e_dev, &
                         v_send_s_dev, v_send_e_dev, &
                         SZ*n_halo*n_block, nproc, pprev, pnext)

    call exec_dist_transeq_3fused(r_u_dev, &
                                  u_dev, u_recv_s_dev, u_recv_e_dev, &
                                  v_dev, v_recv_s_dev, v_recv_e_dev, &
                                  dud_dev, d2u_dev, &
                                  du_send_s_dev, du_send_e_dev, &
                                  du_recv_s_dev, du_recv_e_dev, &
                                  dud_send_s_dev, dud_send_e_dev, &
                                  dud_recv_s_dev, dud_recv_e_dev, &
                                  d2u_send_s_dev, d2u_send_e_dev, &
                                  d2u_recv_s_dev, d2u_recv_e_dev, &
                                  der1st, der1st, der2nd, nu, nproc, pprev, pnext, &
                                  blocks, threads)
  end subroutine run_kernel

end program test_cuda_transeq
