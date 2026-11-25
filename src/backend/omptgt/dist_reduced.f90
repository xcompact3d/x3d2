module m_omptgt_dist_reduced
  !! Exact solve of the reduced system of the distributed solver on the
  !! device for the OpenMP target backend, following
  !! cuda/dist_reduced.f90 (see m_dist_reduced for the system). Each rank
  !! packs its first and last values after the first phase into one buffer,
  !! and a single gather along the direction collects every rank's straight
  !! from device memory, as the halo exchange of m_omptgt_sendrecv already
  !! relies on a GPU-aware MPI. The system is then solved one line per
  !! iteration, with the factorisation made when the operator was set up.
  !! Nothing goes through the host.
  use m_mpi, only: MPI_Allgather

  use m_common, only: dp, MPI_X3D2_DP
  use m_logger, only: abort_run, check_mpi
  use m_mesh_content, only: par_t
  use m_dist_reduced, only: direction_comm
  use m_omptgt_common, only: SZ, dev_alloc, dev_free, device_sync
  use m_omptgt_kernels_dist, only: pack_reduced_rows, solve_reduced_lines
  use m_omptgt_tdsops, only: omptgt_tdsops_t

  implicit none

  private
  public :: solve_reduced_exact_omptgt

  !> This rank's packed rows and the gathered rows of every rank, kept
  !! between calls so that the solve every step does not allocate
  real(dp), pointer, contiguous, save :: send_rows(:) => null(), &
                                         all_rows(:) => null()

contains

  subroutine solve_reduced_exact_omptgt(x_s, x_e, d_s, d_e, tdsops, par, &
                                        dir, n_groups)
    !! x_s and x_e get this rank's first and last unknowns of the first
    !! n_groups groups of lines, from d_s and d_e, the first and last values
    !! after the first phase. Line l must be the same physical line on every
    !! rank along dir.
    implicit none

    real(dp), pointer, contiguous, dimension(:, :, :), intent(in) :: &
      x_s, x_e, d_s, d_e
    type(omptgt_tdsops_t), intent(in) :: tdsops
    type(par_t), intent(in) :: par
    integer, intent(in) :: dir, n_groups

    integer :: np, nl, comm, ierr

    np = par%nproc_dir(dir)
    nl = SZ*n_groups
    if (tdsops%reduced%m /= 2*np &
        .or. tdsops%reduced%me /= par%nrank_dir(dir)) &
      call abort_run('solve_reduced_exact_omptgt: the reduced system was '// &
                     'not factorised for this rank layout.')

    if (associated(send_rows)) then
      if (size(send_rows) /= 2*nl .or. size(all_rows) /= 2*nl*np) then
        call dev_free(send_rows)
        call dev_free(all_rows)
      end if
    end if
    if (.not. associated(send_rows)) then
      call dev_alloc(send_rows, 2*nl)
      call dev_alloc(all_rows, 2*nl*np)
    end if

    comm = direction_comm(par, dir)

    call pack_reduced_rows(send_rows, d_s, d_e, n_groups)

    ! The packing must be complete before MPI reads device memory
    call device_sync()

#ifdef MPI
    call MPI_Allgather(send_rows, 2*nl, MPI_X3D2_DP, all_rows, 2*nl, &
                       MPI_X3D2_DP, comm, ierr)
    call check_mpi(ierr, 'MPI_Allgather in solve_reduced_exact_omptgt')
#endif

    call solve_reduced_lines( &
      x_s, x_e, all_rows, tdsops%reduced_lo_dev, tdsops%reduced_den_dev, &
      tdsops%reduced_cp_dev, tdsops%reduced_z_dev, tdsops%reduced%corner_g, &
      tdsops%reduced%zfac, np, tdsops%reduced%me, n_groups &
      )

  end subroutine solve_reduced_exact_omptgt

end module m_omptgt_dist_reduced
