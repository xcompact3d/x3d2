program x3d2_memcheck
  !! Standalone GPU memory estimator for an x3d2 .x3d input file.
  !!   --static   The formula estimate only.
  !!   (no flag)  The formula estimate, then an FFT plan probe on the whole
  !!              grid.
  !!   --build    The above, then a real case build and one substep on 1 GPU,
  !!              measured with cudaMemGetInfo, and a CHECK line against the
  !!              estimate.
  !! The real build runs in the same process, in a scratch directory under
  !! the invoking directory. The probe gives the CUDA context and the
  !! NVSHMEM heap of the GPU in use. The GPU-aware IO staging buffer is
  !! added by formula.
  !! Usage: x3d2-memcheck <input.x3d> [--static | --build [--extensive <n>]]
  use mpi
  use m_common, only: i8
  use m_memory_estimate, only: output_field_active, peak_fields_lookup
  use m_memcheck_context, only: memcheck_ctx_t, parse_args, read_config, &
                                classify_bc, resolve_gpu_io_mode
  use m_memcheck_estimate, only: to_gib
  use m_memcheck_report, only: report
  use m_memcheck_build, only: run_build_mode

  implicit none

  type(memcheck_ctx_t) :: ctx
  integer :: ierr, nproc
  logical :: device_ok
  character(len=128) :: device_msg

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, ctx%irank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
  if (nproc /= 1) error stop 'x3d2-memcheck: run single-rank (mpirun -n 1).'

  call ctx%device%init(device_ok, device_msg)
  if (.not. device_ok) then
    print '(a)', trim(device_msg)
    call no_device_exit()
  end if

  call parse_args(ctx)
  call read_config(ctx)
  call resolve_gpu_io_mode(ctx)
  call classify_bc(ctx)
  call query_card_gib()

  ctx%peak_fields = peak_fields_lookup( &
                    trim(ctx%domain_cfg%flow_case_name), &
                    ctx%solver_cfg%time_intg, &
                    ctx%solver_cfg%n_species, &
                    trim(ctx%les_cfg%model) /= 'none', &
                    ctx%solver_cfg%ibm_on, &
                    output_field_active(ctx%checkpoint_cfg, 'vorticity'), &
                    output_field_active(ctx%checkpoint_cfg, 'qcriterion'), &
                    ctx%solver_cfg%lowmem_transeq)

  ctx%multi_gpu_supported = .not. (ctx%bc_is_010 .or. ctx%bc_is_110)

  call report(ctx)
  if (trim(ctx%run_mode) == 'BUILD') call run_build_mode(ctx)

  call MPI_Finalize(ierr)

  ! Exit code: 0 FITS, 1 BORDERLINE, 2 DOES_NOT_FIT or UNSUPPORTED.
  select case (trim(ctx%final_verdict))
  case ('FITS')
    call exit(0)
  case ('BORDERLINE')
    call exit(1)
  case default
    call exit(2)
  end select

contains

  subroutine no_device_exit()
    !! Shuts MPI down, then stops.
    integer :: ierr_finalize

    call MPI_Finalize(ierr_finalize)
    error stop 'x3d2-memcheck: no usable CUDA device'
  end subroutine no_device_exit

  subroutine query_card_gib()
    integer(i8) :: free_b, total_b

    call ctx%device%mem_info(total_b, free_b)
    ctx%card_gib = to_gib(total_b)
    ctx%card_free_gib = to_gib(free_b)
  end subroutine query_card_gib

end program x3d2_memcheck
