program test_transeq
  use iso_fortran_env, only: stderr => error_unit

  use m_common, only: dp, pi, BC_PERIODIC
  use m_omp_common, only: SZ
  use m_omp_exec_dist, only: exec_dist_transeq_compact
  use m_omp_sendrecv, only: sendrecv_fields
  use m_tdsops, only: tdsops_t
  use m_test_utils, only: initialise_mpi, finalise_test, relative_l2_error

  implicit none

  logical :: allpass = .true.
  real(dp), allocatable, dimension(:, :, :) :: u, v, r_u, ref
  real(dp), allocatable, dimension(:, :, :) :: dud, d2u ! intermediate solution arrays
  real(dp), allocatable, dimension(:, :, :) :: &
    du_recv_s, du_recv_e, du_send_s, du_send_e, &
    dud_recv_s, dud_recv_e, dud_send_s, dud_send_e, &
    d2u_recv_s, d2u_recv_e, d2u_send_s, d2u_send_e

  real(dp), allocatable, dimension(:, :, :) :: &
    u_send_s, u_send_e, u_recv_s, u_recv_e, &
    v_send_s, v_send_e, v_recv_s, v_recv_e

  type(tdsops_t) :: der1st, der2nd

  integer :: n, n_block, n_halo, n_glob
  integer :: nrank, nproc, pprev, pnext

  real(dp) :: dx_per, nu, norm_du
  ! The tolerance bounds the relative L2 error of the transport equation
  ! right hand side. Roundoff floor of the second-derivative term is
  ! ~eps/dx^2, which at 96^3 cells is ~1e-4 in single precision. Largest
  ! value measured, over the OpenMP backend binary of the OpenMP and CUDA
  ! builds:
  !   double precision: 1.6562e-10
  !   single precision: 1.1935e-04
  ! The tolerance is the smallest {1,2,5}x10^k value at least 10x these, a
  ! margin of about 12 in double precision and 17 in single precision.
#ifdef SINGLE_PREC
  real(dp) :: tol = 2e-3
#else
  real(dp) :: tol = 2d-9
#endif

  call initialise_mpi(nrank, nproc, pprev, pnext)
  if (nrank == 0) print *, 'Parallel run with', nproc, 'ranks'
  call setup_geometry()
  call allocate_fields()
  call initialise_input()
  call setup_tdsops()
  call run_kernel()
  call check_result()
  call finalise_test(allpass, nrank)

contains


  subroutine setup_geometry()
    n_glob = 32*4
    n = n_glob/nproc
    n_block = 32*32/SZ
    n_halo = 4
    nu = 1._dp
    dx_per = 2*pi/n_glob
  end subroutine setup_geometry

  subroutine allocate_fields()
    allocate (u(SZ, n, n_block), v(SZ, n, n_block), r_u(SZ, n, n_block), &
              ref(SZ, n, n_block))

    ! main input fields
    ! field for storing the result
    ! intermediate solution fields
    allocate (dud(SZ, n, n_block))
    allocate (d2u(SZ, n, n_block))

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
  end subroutine initialise_input

  subroutine setup_tdsops()
    ! preprocess the operator and coefficient array
    der1st = tdsops_t(n, dx_per, operation='first-deriv', scheme='compact6', &
                      bc_start=BC_PERIODIC, bc_end=BC_PERIODIC)
    der2nd = tdsops_t(n, dx_per, operation='second-deriv', scheme='compact6', &
                      bc_start=BC_PERIODIC, bc_end=BC_PERIODIC)
  end subroutine setup_tdsops

  subroutine run_kernel()
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
      r_u, dud, d2u, &
      du_send_s, du_send_e, du_recv_s, du_recv_e, &
      dud_send_s, dud_send_e, dud_recv_s, dud_recv_e, &
      d2u_send_s, d2u_send_e, d2u_recv_s, d2u_recv_e, &
      u, u_recv_s, u_recv_e, &
      v, v_recv_s, v_recv_e, &
      der1st, der1st, der2nd, nu, nproc, pprev, pnext, n_block &
      )
  end subroutine run_kernel

  subroutine check_result()
    ! check error
    ref = -v*v + 0.5_dp*u*u - nu*u
    r_u = r_u - ref
    norm_du = relative_l2_error(r_u, ref)

    if (nrank == 0) print *, 'error norm', norm_du

    if (nrank == 0) then
      if (norm_du > tol) then
        allpass = .false.
        write (stderr, '(a)') 'Check second derivatives... failed'
      else
        write (stderr, '(a)') 'Check second derivatives... passed'
      end if
    end if

  end subroutine check_result

end program test_transeq
