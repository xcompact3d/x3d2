program test_statistics
  !! Unit tests for the accumulate_mean subroutine in m_stats.

  use m_common, only: dp
  use m_stats, only: accumulate_mean
  use mpi

  implicit none

  integer :: ierr

  call MPI_Init(ierr)

  call test_constant_field()
  call test_known_mean()
  call test_rms_nonzero()
  call test_reynolds_stress()

  print *, 'All statistics accumulation tests passed'

  call MPI_Finalize(ierr)

contains

  subroutine test_constant_field()
    !! Constant field c=1: umean==1, uprime==0.
    real(dp), parameter :: c = 1.0_dp
    integer, parameter :: n_samples = 100
    real(dp) :: umean(1, 1, 1), uumean(1, 1, 1), val(1, 1, 1)
    real(dp) :: uprime, stat_inc
    integer :: n
    real(dp), parameter :: tol = 1.0e-12_dp

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
      error stop 'test_constant_field: umean /= constant'
    end if
    if (abs(uprime) > tol) then
      print *, 'FAIL test_constant_field: uprime =', uprime, 'expected 0'
      error stop 'test_constant_field: uprime /= 0 for constant field'
    end if
  end subroutine test_constant_field

  subroutine test_known_mean()
    !! Values 1..N: mean must equal (N+1)/2.
    integer, parameter :: n_samples = 50
    real(dp) :: umean(1, 1, 1), val(1, 1, 1)
    real(dp) :: stat_inc, expected
    integer :: n
    real(dp), parameter :: tol = 1.0e-12_dp

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
      error stop 'test_known_mean: mean /= analytical mean'
    end if
  end subroutine test_known_mean

  subroutine test_rms_nonzero()
    !! Alternating -1/+1: mean==0, uprime==1.
    integer, parameter :: n_samples = 200
    real(dp) :: umean(1, 1, 1), uumean(1, 1, 1), val(1, 1, 1)
    real(dp) :: uprime, stat_inc
    integer :: n
    real(dp), parameter :: tol = 1.0e-10_dp

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
      error stop 'test_rms_nonzero: mean of alternating series /= 0'
    end if
    if (abs(uprime - 1.0_dp) > tol) then
      print *, 'FAIL test_rms_nonzero: uprime =', uprime, 'expected 1'
      error stop 'test_rms_nonzero: rms of alternating series /= 1'
    end if
  end subroutine test_rms_nonzero

  subroutine test_reynolds_stress()
    !! u=v=1..N: Reynolds stress <u'v'> = <uv> - <u><v> must equal var(u).
    !! u=-v: Reynolds stress must equal -var(u).
    integer, parameter :: n_samples = 100
    real(dp) :: umean(1, 1, 1), vmean(1, 1, 1)
    real(dp) :: uumean(1, 1, 1), uvmean(1, 1, 1)
    real(dp) :: val(1, 1, 1), stat_inc
    real(dp) :: uprime2, reynolds_stress
    integer :: n
    real(dp), parameter :: tol = 1.0e-10_dp

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
      error stop 'test_reynolds_stress: correlated case failed'
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
      error stop 'test_reynolds_stress: anticorrelated case failed'
    end if
  end subroutine test_reynolds_stress

end program test_statistics
