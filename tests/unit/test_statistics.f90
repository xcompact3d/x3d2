program test_statistics
  !! Unit tests for the accumulate_mean subroutine in m_stats.

  use m_common, only: dp
  use m_stats, only: accumulate_mean
  use m_test_utils, only: initialise_mpi, finalise_test

  implicit none

  integer :: nrank, nproc
  logical :: allpass = .true.

  call initialise_mpi(nrank, nproc)

  call test_constant_field(allpass)
  call test_known_mean(allpass)
  call test_rms_nonzero(allpass)
  call test_reynolds_stress(allpass)

  call finalise_test(allpass, nrank)

contains

  subroutine test_constant_field(allpass)
    !! Constant field c=1: umean==1, uprime==0.
    logical, intent(inout) :: allpass
    real(dp), parameter :: c = 1.0_dp
    integer, parameter :: n_samples = 100
    real(dp) :: umean(1, 1, 1), uumean(1, 1, 1), val(1, 1, 1)
    real(dp) :: uprime, stat_inc
    integer :: n
    real(dp), parameter :: tol = 1.0e4_dp*epsilon(1.0_dp)

    umean = 0.0_dp; uumean = 0.0_dp

    do n = 1, n_samples
      stat_inc = 1.0_dp/n
      val(1, 1, 1) = c
      call accumulate_mean(umean, val, stat_inc)
      val(1, 1, 1) = c*c
      call accumulate_mean(uumean, val, stat_inc)
    end do

    uprime = sqrt(max(0.0_dp, uumean(1, 1, 1) - umean(1, 1, 1)**2))

    if (abs(umean(1, 1, 1) - c) > tol) then
      print *, 'FAIL test_constant_field: umean =', umean(1, 1, 1), &
        'expected', c
      allpass = .false.
    end if
    if (abs(uprime) > tol) then
      print *, 'FAIL test_constant_field: uprime =', uprime, 'expected 0'
      allpass = .false.
    end if
  end subroutine test_constant_field

  subroutine test_known_mean(allpass)
    !! Values 1..N: mean must equal (N+1)/2.
    logical, intent(inout) :: allpass
    integer, parameter :: n_samples = 50
    real(dp) :: umean(1, 1, 1), val(1, 1, 1)
    real(dp) :: stat_inc, expected
    integer :: n
    real(dp), parameter :: tol = 1.0e4_dp*epsilon(1.0_dp)

    umean = 0.0_dp
    expected = (n_samples + 1)/2.0_dp

    do n = 1, n_samples
      stat_inc = 1.0_dp/n
      val(1, 1, 1) = real(n, dp)
      call accumulate_mean(umean, val, stat_inc)
    end do

    if (abs(umean(1, 1, 1) - expected) > tol) then
      print *, 'FAIL test_known_mean: umean =', umean(1, 1, 1), &
        'expected', expected
      allpass = .false.
    end if
  end subroutine test_known_mean

  subroutine test_rms_nonzero(allpass)
    !! Alternating -1/+1: mean==0, uprime==1.
    logical, intent(inout) :: allpass
    integer, parameter :: n_samples = 200
    real(dp) :: umean(1, 1, 1), uumean(1, 1, 1), val(1, 1, 1)
    real(dp) :: uprime, stat_inc
    integer :: n
    real(dp), parameter :: tol = 1.0e4_dp*epsilon(1.0_dp)

    umean = 0.0_dp; uumean = 0.0_dp

    do n = 1, n_samples
      val(1, 1, 1) = merge(1.0_dp, -1.0_dp, mod(n, 2) == 0)
      stat_inc = 1.0_dp/n
      call accumulate_mean(umean, val, stat_inc)
      val(1, 1, 1) = val(1, 1, 1)**2
      call accumulate_mean(uumean, val, stat_inc)
    end do

    uprime = sqrt(max(0.0_dp, uumean(1, 1, 1) - umean(1, 1, 1)**2))

    if (abs(umean(1, 1, 1)) > tol) then
      print *, 'FAIL test_rms_nonzero: umean =', umean(1, 1, 1), &
        'expected 0'
      allpass = .false.
    end if
    if (abs(uprime - 1.0_dp) > tol) then
      print *, 'FAIL test_rms_nonzero: uprime =', uprime, 'expected 1'
      allpass = .false.
    end if
  end subroutine test_rms_nonzero

  subroutine test_reynolds_stress(allpass)
    !! u=v=1..N: Reynolds stress <u'v'> = <uv> - <u><v> must equal var(u).
    !! u=-v: Reynolds stress must equal -var(u).
    logical, intent(inout) :: allpass
    integer, parameter :: n_samples = 100
    real(dp) :: umean(1, 1, 1), vmean(1, 1, 1)
    real(dp) :: uumean(1, 1, 1), uvmean(1, 1, 1)
    real(dp) :: val(1, 1, 1), stat_inc
    real(dp) :: uprime2, reynolds_stress
    integer :: n
    real(dp), parameter :: tol = 1.0e4_dp*epsilon(1.0_dp)

    ! Perfectly correlated: u == v
    umean = 0.0_dp; vmean = 0.0_dp
    uumean = 0.0_dp; uvmean = 0.0_dp

    do n = 1, n_samples
      stat_inc = 1.0_dp/n
      val(1, 1, 1) = real(n, dp)
      call accumulate_mean(umean, val, stat_inc)
      call accumulate_mean(vmean, val, stat_inc)
      val(1, 1, 1) = real(n, dp)**2
      call accumulate_mean(uumean, val, stat_inc)
      call accumulate_mean(uvmean, val, stat_inc)
    end do

    uprime2 = uumean(1, 1, 1) - umean(1, 1, 1)**2
    reynolds_stress = uvmean(1, 1, 1) - umean(1, 1, 1)*vmean(1, 1, 1)

    if (abs(reynolds_stress - uprime2) > tol*abs(uprime2)) then
      print *, 'FAIL test_reynolds_stress: <u''v''> /= var(u) for u==v'
      allpass = .false.
    end if

    ! Anticorrelated: u = -v
    umean = 0.0_dp; vmean = 0.0_dp
    uumean = 0.0_dp; uvmean = 0.0_dp

    do n = 1, n_samples
      stat_inc = 1.0_dp/n
      val(1, 1, 1) = real(n, dp)
      call accumulate_mean(umean, val, stat_inc)
      val(1, 1, 1) = -real(n, dp)
      call accumulate_mean(vmean, val, stat_inc)
      val(1, 1, 1) = real(n, dp)**2
      call accumulate_mean(uumean, val, stat_inc)
      val(1, 1, 1) = -real(n, dp)**2
      call accumulate_mean(uvmean, val, stat_inc)
    end do

    uprime2 = uumean(1, 1, 1) - umean(1, 1, 1)**2
    reynolds_stress = uvmean(1, 1, 1) - umean(1, 1, 1)*vmean(1, 1, 1)

    if (abs(reynolds_stress + uprime2) > tol*abs(uprime2)) then
      print *, 'FAIL test_reynolds_stress: <u''v''> /= -var(u) for u==-v'
      allpass = .false.
    end if
  end subroutine test_reynolds_stress

end program test_statistics
