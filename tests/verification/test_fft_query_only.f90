program test_fft_query_only
  !! Check the query-only mode of the x3d2-memcheck FFT workspace query. The
  !! predicted heap is the plan memory plus the workspace rounded up to 1 GiB.
  !! Below 1 GiB of workspace, NVSHMEM fits the workspace into spare heap,
  !! so the prediction is an upper bound there: it is not below the measured
  !! heap and at most one rounding step (1 GiB) plus 0.05 GiB above it. Small
  !! grids always fit, so the probe runs for them and the predicted heap is
  !! not used. From 512^3 the prediction must match the measured heap within
  !! 0.1 GiB.
  !! The query-only mode allocates no xtdesc buffer. The rounding helper is
  !! also checked on its own.

  use m_common, only: i8, is_sp
  use m_memory_estimate, only: predicted_fft_heap_bytes, round_up_bytes
  use m_cuda_memcheck_device, only: fft_workspace_bytes_query
  use m_test_utils, only: initialise_mpi, finalise_test

  implicit none

  integer, parameter :: grids(3) = [256, 384, 512]
  !> The bound from 512^3. A dropped rounding of the workspace moves the
  !> prediction by up to 1 GiB.
  integer(i8), parameter :: max_diff_bytes = 107374182_i8
  !> The grids below this size only need an upper bound.
  integer, parameter :: exact_from = 512
  integer(i8), parameter :: gib = 1024_i8**3
  !> The upper bound allows one rounding step plus 0.05 GiB.
  integer(i8), parameter :: max_above_bytes = gib + 53687091_i8

  integer :: nrank, nproc, stat_n, stat_q, ig, cdims(3)
  integer(i8) :: ws_n, heap_n, xt_n, ws_q, heap_q, xt_q
  integer(i8) :: predicted, measured, diff
  logical :: used_n, used_q, allpass, in_bound

  call initialise_mpi(nrank, nproc)
  allpass = .true.

  ! The rounding helper: exact multiples stay, anything above goes up.
  if (round_up_bytes(gib, gib) /= gib .or. &
      round_up_bytes(gib + 1_i8, gib) /= 2_i8*gib .or. &
      round_up_bytes(gib - 1_i8, gib) /= gib .or. &
      round_up_bytes(0_i8, gib) /= 0_i8) then
    if (nrank == 0) print *, 'FAIL: round_up_bytes'
    allpass = .false.
  end if

  do ig = 1, size(grids)
    cdims = grids(ig)
    call fft_workspace_bytes_query(.false., .false., cdims, .true., &
                                   nrank == 0, ws_n, heap_n, xt_n, used_n, &
                                   stat_n)
    call fft_workspace_bytes_query(.false., .false., cdims, .true., &
                                   nrank == 0, ws_q, heap_q, xt_q, used_q, &
                                   stat_q, query_only=.true.)

    if (stat_n /= 0 .or. stat_q /= 0) then
      if (nrank == 0) print *, 'FAIL: query status ', stat_n, stat_q
      allpass = .false.
      cycle
    end if
    if (.not. (used_n .and. used_q)) then
      if (nrank == 0) print *, 'FAIL: the query did not use cuFFTMp'
      allpass = .false.
      cycle
    end if
    if (xt_q /= 0_i8) then
      if (nrank == 0) print *, 'FAIL: query-only mode reported an xtdesc size'
      allpass = .false.
    end if
    predicted = predicted_fft_heap_bytes(heap_q, ws_q)
    measured = heap_n
    if (nrank == 0) print '(a,i0,2(a,f8.3),a)', 'grid ', grids(ig), &
      ': predicted ', real(predicted)/1024.**3, ' GiB, measured ', &
      real(measured)/1024.**3, ' GiB'
    diff = predicted - measured
    ! In single precision the 512^3 workspace is under 1 GiB, so the
    ! prediction is an upper bound there.
    if (grids(ig) >= exact_from .and. .not. is_sp) then
      in_bound = abs(diff) <= max_diff_bytes
    else
      in_bound = diff >= 0_i8 .and. diff <= max_above_bytes
    end if
    if (.not. in_bound) then
      if (nrank == 0) print *, 'FAIL: the predicted heap is outside the bound'
      allpass = .false.
    end if
  end do

  call finalise_test(allpass, nrank)
end program test_fft_query_only
