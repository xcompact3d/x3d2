program test_relative_l2_error
  !! Checks m_test_utils::relative_l2_error against values known in advance.

  use m_common, only: dp
  use m_test_utils, only: finalise_test, initialise_mpi, relative_l2_error

  implicit none

  ! dp is single in SINGLE_PREC builds, so the tolerance scales with epsilon.
  real(dp), parameter :: tol = 100*epsilon(1.0_dp)

  integer :: nrank, nproc, r, n_r
  logical :: allpass
  real(dp), allocatable :: ref(:, :, :), err(:, :, :)
  real(dp) :: res, expected, num, den, a_r, c_r

  call initialise_mpi(nrank, nproc)
  allpass = .true.

  ! 1. Exact agreement gives exactly zero.
  allocate (ref(3, 4, 5), err(3, 4, 5))
  ref = 1.0_dp
  err = 0.0_dp
  res = relative_l2_error(err, ref)
  if (res /= 0.0_dp) then
    print *, 'exact agreement: expected 0, got ', res
    allpass = .false.
  end if
  deallocate (ref, err)

  ! 2. A known relative error of 0.1 on a non-constant reference.
  allocate (ref(3, 4, 5), err(3, 4, 5))
  call random_number(ref)
  ref = ref + 0.5_dp + real(nrank, dp)
  err = 0.1_dp*ref
  res = relative_l2_error(err, ref)
  if (abs(res - 0.1_dp) > tol) then
    print *, 'known error: expected 0.1, got ', res
    allpass = .false.
  end if
  deallocate (ref, err)

  ! 3. Uneven partitions. Rank r holds n_r = 2 + r points along the first
  ! dimension, ref = a_r = 1 + r and err = c_r*ref with c_r = 0.1*(r + 1).
  ! The global error is sqrt(sum n_r (c_r a_r)^2 / sum n_r a_r^2), computed
  ! here from nproc without MPI. For nproc >= 2 it differs from the mean of
  ! c_r, so the test also fails for a per-rank average of the local errors.
  n_r = 2 + nrank
  a_r = 1.0_dp + real(nrank, dp)
  c_r = 0.1_dp*real(nrank + 1, dp)
  allocate (ref(n_r, 4, 5), err(n_r, 4, 5))
  ref = a_r
  err = c_r*ref
  res = relative_l2_error(err, ref)

  num = 0.0_dp
  den = 0.0_dp
  do r = 0, nproc - 1
    num = num + real(2 + r, dp)*(0.1_dp*real(r + 1, dp)*(1.0_dp + real(r, dp)))**2
    den = den + real(2 + r, dp)*(1.0_dp + real(r, dp))**2
  end do
  expected = sqrt(num/den)
  if (abs(res - expected) > tol) then
    print *, 'uneven partitions: expected ', expected, ' got ', res
    allpass = .false.
  end if
  if (nproc >= 2) then
    ! The mean of c_r over the ranks is 0.05*(nproc + 1).
    if (abs(expected - 0.05_dp*real(nproc + 1, dp)) <= 1.0e-3_dp) then
      print *, 'uneven partitions: global value equals the mean of c_r'
      allpass = .false.
    end if
  end if

  call finalise_test(allpass, nrank)

end program test_relative_l2_error
