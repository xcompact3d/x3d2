program test_cuda_transeq
  use iso_fortran_env, only: stderr => error_unit
  use cudafor

  use m_common, only: dp, pi, BC_PERIODIC
  use m_cuda_common, only: SZ
  use m_cuda_exec_dist, only: exec_dist_transeq_3fused
  use m_cuda_sendrecv, only: sendrecv_fields
  use m_cuda_tdsops, only: cuda_tdsops_t
  use m_backend_runtime, only: select_device
  use m_test_utils, only: initialise_mpi, finalise_test, relative_l2_error

  implicit none

  logical :: allpass = .true.
  real(dp), allocatable, dimension(:, :, :) :: u, v, r_u, ref
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

  integer :: n, n_block, n_halo, n_glob
  integer :: nrank, nproc, pprev, pnext

  type(dim3) :: blocks, threads
  ! The tolerance bounds the relative L2 error of the transport equation
  ! right hand side. The diffusion term's roundoff floor is ~eps/dx^2: at
  ! n=512 this is ~8e-4 in single precision, so the tolerance must account
  ! for it. Largest value measured on the CUDA backend:
  !   double precision: 4.0869e-12
  !   single precision: 1.1078e-03
  ! The tolerance is the smallest {1,2,5}x10^k value at least 10x these, a
  ! margin of about 12 in double precision and 18 in single precision.
#ifdef SINGLE_PREC
  real(dp), parameter :: residual_tol = 2.0e-2_dp
#else
  real(dp), parameter :: residual_tol = 5.0e-11_dp
#endif
  real(dp) :: dx_per, nu, norm_du

  call initialise_mpi(nrank, nproc, pprev, pnext)
  if (nrank == 0) print *, 'Parallel run with', nproc, 'ranks'
  call select_device(nrank)
  call setup_geometry()
  call allocate_fields()
  call initialise_input()
  call setup_backend()
  call run_kernel()
  call check_result()
  call finalise_test(allpass, nrank)

contains

  subroutine setup_geometry()
    n_glob = 512
    n = n_glob/nproc
    n_block = 512*512/SZ
    n_halo = 4
    nu = 1._dp
    dx_per = 2*pi/n_glob
  end subroutine setup_geometry

  subroutine allocate_fields()
    allocate (u(SZ, n, n_block), v(SZ, n, n_block), r_u(SZ, n, n_block), &
              ref(SZ, n, n_block))

    ! main input fields
    allocate (u_dev(SZ, n, n_block), v_dev(SZ, n, n_block))
    ! field for storing the result
    allocate (r_u_dev(SZ, n, n_block))
    ! intermediate solution fields
    allocate (dud_dev(SZ, n, n_block))
    allocate (d2u_dev(SZ, n, n_block))

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
  end subroutine allocate_fields

  subroutine initialise_input()
    integer :: i, j, k

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
  end subroutine initialise_input

  subroutine setup_backend()
    ! preprocess the operator and  coefficient arrays
    der1st = cuda_tdsops_t(n, dx_per, operation='first-deriv', &
                           scheme='compact6', &
                           bc_start=BC_PERIODIC, bc_end=BC_PERIODIC)
    der2nd = cuda_tdsops_t(n, dx_per, operation='second-deriv', &
                           scheme='compact6', &
                           bc_start=BC_PERIODIC, bc_end=BC_PERIODIC)

    blocks = dim3(n_block, 1, 1)
    threads = dim3(SZ, 1, 1)
  end subroutine setup_backend

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

  subroutine check_result()
    ! check error
    r_u = r_u_dev
    ref = -v*v + 0.5_dp*u*u - nu*u
    r_u = r_u - ref
    norm_du = relative_l2_error(r_u, ref)

    if (nrank == 0) print *, 'error norm', norm_du

    if (nrank == 0) then
      if (norm_du > residual_tol) then
        allpass = .false.
        write (stderr, '(a)') 'Check second derivatives... failed'
      else
        write (stderr, '(a)') 'Check second derivatives... passed'
      end if
    end if

  end subroutine check_result

end program test_cuda_transeq
