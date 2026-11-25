module m_cuda_dist_reduced
  !! Exact solve of the reduced system of the distributed solver on the
  !! device (see m_dist_reduced for the system). Each rank packs its first
  !! and last values after the first phase into one buffer, and a single
  !! gather along the direction collects every rank's straight from device
  !! memory, as the halo exchange of m_cuda_sendrecv already relies on a
  !! GPU-aware MPI. The system is then solved with one thread per line, with
  !! the factorisation made when the operator was set up. Nothing goes
  !! through the host.
  use cudafor
  use m_mpi, only: MPI_Allgather

  use m_common, only: dp, MPI_X3D2_DP
  use m_logger, only: abort_run, check_mpi
  use m_mesh_content, only: par_t
  use m_dist_reduced, only: direction_comm
  use m_cuda_errors, only: check_cuda
  use m_cuda_kernels_dist, only: pack_reduced_rows, solve_reduced_lines
  use m_cuda_tdsops, only: cuda_tdsops_t

  implicit none

  private
  public :: solve_reduced_exact_cuda

  !> This rank's packed rows and the gathered rows of every rank, kept
  !! between calls so that the solve every step does not allocate
  real(dp), device, allocatable, save :: send_rows(:), all_rows(:)

contains

  subroutine solve_reduced_exact_cuda(x_s, x_e, d_s, d_e, tdsops, par, &
                                      dir, blocks, threads)
    !! x_s and x_e get this rank's first and last unknowns of each line,
    !! from d_s and d_e, the first and last values after the first phase.
    !! Line l must be the same physical line on every rank along dir.
    implicit none

    real(dp), device, dimension(:, :, :), intent(out) :: x_s, x_e
    real(dp), device, dimension(:, :, :), intent(in) :: d_s, d_e
    type(cuda_tdsops_t), intent(in) :: tdsops
    type(par_t), intent(in) :: par
    integer, intent(in) :: dir
    !> The launch configuration of the first phase: one thread per line
    type(dim3), intent(in) :: blocks, threads

    integer :: np, nl, comm, ierr

    np = par%nproc_dir(dir)
    nl = threads%x*blocks%x
    if (tdsops%reduced%m /= 2*np &
        .or. tdsops%reduced%me /= par%nrank_dir(dir)) &
      call abort_run('solve_reduced_exact_cuda: the reduced system was '// &
                     'not factorised for this rank layout.')

    if (allocated(send_rows)) then
      if (size(send_rows) /= 2*nl .or. size(all_rows) /= 2*nl*np) &
        deallocate (send_rows, all_rows)
    end if
    if (.not. allocated(send_rows)) &
      allocate (send_rows(2*nl), all_rows(2*nl*np))

    comm = direction_comm(par, dir)

    call pack_reduced_rows<<<blocks, threads>>>(send_rows, d_s, d_e) !&

    ! The packing must be complete before MPI reads device memory
    ierr = cudaStreamSynchronize(0)
    call check_cuda(ierr, 'cudaStreamSynchronize in solve_reduced_exact_cuda')

#ifdef MPI
    call MPI_Allgather(send_rows, 2*nl, MPI_X3D2_DP, all_rows, 2*nl, &
                       MPI_X3D2_DP, comm, ierr)
    call check_mpi(ierr, 'MPI_Allgather in solve_reduced_exact_cuda')
#endif

    call solve_reduced_lines<<<blocks, threads>>>( & !&
      x_s, x_e, all_rows, &
      tdsops%reduced_lo_dev, tdsops%reduced_den_dev, &
      tdsops%reduced_cp_dev, tdsops%reduced_z_dev, &
      tdsops%reduced%corner_g, tdsops%reduced%zfac, np, &
      tdsops%reduced%me &
      )

  end subroutine solve_reduced_exact_cuda

end module m_cuda_dist_reduced
