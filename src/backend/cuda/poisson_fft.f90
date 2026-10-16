module m_cuda_poisson_fft
  use iso_c_binding, only: c_loc, c_ptr, c_f_pointer, c_int, c_float, &
                           c_double_complex, c_float_complex
  use cudafor
  use cufftXt
  use cufft
  use m_mpi, only: MPI_COMM_WORLD, MPI_COMPLEX, MPI_DOUBLE_COMPLEX, &
                   MPI_STATUS_IGNORE, MPI_Sendrecv, &
                   MPI_Alltoall

  use m_common, only: dp, CELL, is_sp, MPI_X3D2_DP
  use m_field, only: field_t
  use m_mesh
  use m_poisson_fft, only: poisson_fft_t, poisson_bc_code, spectral_dirs
  use m_tdsops, only: dirps_t

  use m_cuda_allocator, only: cuda_field_t
  use m_cuda_sendrecv, only: default_stream
  use m_cuda_spectral, only: memcpy3D, memcpy3D_with_transpose, &
                             memcpy3D_with_transpose_back, &
                             transpose_xyz_to_zxy, transpose_zxy_to_xyz, &
                             process_spectral_000, process_spectral_010, &
                             process_spectral_110_norm_z, &
                             process_spectral_110_x_pair_fw, &
                             process_spectral_110_y_pair_fw, &
                             process_spectral_110_poisson, &
                             process_spectral_110_y_pair_bw, &
                             process_spectral_110_x_pair_bw, &
                             process_spectral_110_x_pair_fw_mirror, &
                             process_spectral_110_x_pair_bw_mirror, &
                             process_spectral_110_z_bw, &
                             enforce_periodicity_x, undo_periodicity_x, &
                             enforce_periodicity_y, undo_periodicity_y, &
                             enforce_periodicity_xy, undo_periodicity_xy, &
                             process_spectral_010_fw, &
                             process_spectral_010_poisson, &
                             process_spectral_010_bw, &
                             process_spectral_100_fw, &
                             process_spectral_100_pair_fw, &
                             process_spectral_100_pair_bw, &
                             pack_spectral_plane, &
                             pack_110_xyz_to_yslab, unpack_110_yslab, &
                             pack_110_yslab, unpack_110_yslab_to_xyz
  use m_logger, only: log_info, log_warn, abort_run, check_mpi, str
  use m_cuda_errors, only: check_cuda, check_cufft

  implicit none

  type, abstract, extends(poisson_fft_t) :: cuda_poisson_fft_t
    !! FFT based Poisson solver, the state and the plain transform that
    !! every boundary condition case shares. Each case extends this type
    !! with its own buffers, see cuda_poisson_000_t to cuda_poisson_110_t.

    !> Local domain sized array storing the spectral equivalence constants
    complex(dp), device, allocatable, dimension(:, :, :) :: waves_dev
    !> Wave numbers in x, y, and z
    real(dp), device, allocatable, dimension(:) :: ax_dev, bx_dev, &
                                                   ay_dev, by_dev, &
                                                   az_dev, bz_dev
    !> Forward and backward FFT transform plans
    integer :: plan3D_fw, plan3D_bw

    !> Flag to indicate whether cuFFTMp is used
    logical :: use_cufftmp = .true.

    !> cuFFTMp object manages decomposition and data storage
    type(cudaLibXtDesc), pointer :: xtdesc

    !> Standard cuFFT storage
    complex(dp), device, allocatable, dimension(:, :, :) :: c_dev
    !> cuFFT real workspace (input/output of R2C/C2R)
    real(dp), device, allocatable, dimension(:, :, :) :: r_dev
    !> dim2 mirror of the local spectral slab, used by the 100 and 110 cases
    !> to reach the paired split partner that lives on another rank
    complex(dp), device, allocatable, dimension(:, :, :) :: c_mirror_dev
    !> Staging copy of the local slab, sent in place of c_dev itself.
    !> In the cuFFTMp path c_dev is a c_f_pointer onto descriptor%data(1),
    !> which is cufftXtMalloc memory and therefore lives in the NVSHMEM
    !> symmetric heap. UCX does not reliably classify such a pointer as
    !> device memory: when it falls back to the host CMA transport it calls
    !> process_vm_readv on it and the job aborts with EFAULT. Copying into
    !> ordinary device memory first makes the exchange independent of which
    !> transport UCX happens to select.
    complex(dp), device, allocatable, dimension(:, :, :) :: c_slab_send_dev
    !> The single dim2 plane that falls outside the mirrored slab
    complex(dp), device, allocatable, dimension(:, :) :: c_plane_send_dev, &
                                                         c_plane_recv_dev
    !> Ranks holding the mirror slab and the extra mirror plane
    integer :: mirror_slab_rank = 0, mirror_plane_rank = 0
  contains
    procedure :: init => init_cuda
    procedure :: check_decomposition
    procedure :: spectral_extents
    procedure :: init_case
    procedure :: create_plans
    procedure :: init_dim2_mirror
    procedure :: fft_forward => fft_forward_cuda
    procedure :: fft_backward => fft_backward_cuda
    procedure :: exchange_dim2_mirror
  end type cuda_poisson_fft_t

  type, extends(cuda_poisson_fft_t) :: cuda_poisson_000_t
    !! Periodic in x, y and z
  contains
    procedure :: solve_poisson => solve_poisson_000_cuda
    procedure :: fft_postprocess => fft_postprocess_000_cuda
  end type cuda_poisson_000_t

  type, extends(cuda_poisson_fft_t) :: cuda_poisson_010_t
    !! Periodic in x and z, non-periodic in y, the only case that supports
    !! stretching in y
    !> Stretching operator matrices stores
    real(dp), device, allocatable, dimension(:, :, :, :) :: &
      store_a_odd_re_dev, store_a_odd_im_dev, &
      store_a_even_re_dev, store_a_even_im_dev, &
      store_a_re_dev, store_a_im_dev
    !> Stretching operator matrices
    real(dp), device, allocatable, dimension(:, :, :, :) :: &
      a_odd_re_dev, a_odd_im_dev, a_even_re_dev, a_even_im_dev, &
      a_re_dev, a_im_dev
  contains
    procedure :: check_decomposition => check_decomposition_010
    procedure :: init_case => init_case_010
    procedure :: fft_postprocess => fft_postprocess_010_cuda
    procedure :: enforce_periodicity => enforce_periodicity_y_cuda
    procedure :: undo_periodicity => undo_periodicity_y_cuda
  end type cuda_poisson_010_t

  type, extends(cuda_poisson_fft_t) :: cuda_poisson_100_t
    !! Non-periodic in x, periodic in y and z
  contains
    procedure :: check_decomposition => check_decomposition_100
    procedure :: init_case => init_case_100
    procedure :: fft_forward => fft_forward_100_cuda
    procedure :: fft_backward => fft_backward_100_cuda
    procedure :: fft_postprocess => fft_postprocess_100_cuda
    procedure :: enforce_periodicity => enforce_periodicity_x_cuda
    procedure :: undo_periodicity => undo_periodicity_x_cuda
  end type cuda_poisson_100_t

  type, extends(cuda_poisson_fft_t) :: cuda_poisson_110_t
    !! Non-periodic in x and y, periodic in z
    !> Transposed real workspace (nz, nx, ny), single rank only
    real(dp), device, allocatable, dimension(:, :, :) :: r_dev_110
    !> All-to-all buffers redistributing the field between the physical
    !> Z-slab (nx, ny, nz_l) and the cuFFTMp Y-slab (2*(nz/2+1), nx, ny_l),
    !> sized (nz_l, nx, ny_l, nproc). Ordinary device memory, for the same
    !> reason as c_slab_send_dev.
    real(dp), device, allocatable, dimension(:, :, :, :) :: r_a2a_send_dev, &
                                                            r_a2a_recv_dev
  contains
    procedure :: check_decomposition => check_decomposition_110
    procedure :: spectral_extents => spectral_extents_110
    procedure :: init_case => init_case_110
    procedure :: create_plans => create_plans_110
    procedure :: fft_forward => fft_forward_110_cuda
    procedure :: fft_backward => fft_backward_110_cuda
    procedure :: fft_postprocess => fft_postprocess_110_cuda
    procedure :: enforce_periodicity => enforce_periodicity_xy_cuda
    procedure :: undo_periodicity => undo_periodicity_xy_cuda
    procedure :: redistribute_110_to_yslab
    procedure :: redistribute_110_to_zslab
  end type cuda_poisson_110_t

  ! Explicit C interfaces for cuFFT functions that nvfortran has trouble with
  interface
    integer(c_int) function cufftExecR2C_C(plan, idata, odata) &
      bind(C, name='cufftExecR2C')
      use iso_c_binding
      integer(c_int), value :: plan
      type(c_ptr), value :: idata
      type(c_ptr), value :: odata
    end function cufftExecR2C_C

    integer(c_int) function cufftExecC2R_C(plan, idata, odata) &
      bind(C, name='cufftExecC2R')
      use iso_c_binding
      integer(c_int), value :: plan
      type(c_ptr), value :: idata
      type(c_ptr), value :: odata
    end function cufftExecC2R_C
  end interface

  private :: create_fft_plan, sync_or_abort, yslab_view_110, &
             redistribution_launch_config, redistribution_alltoall, &
             init_cuda, require_z_decomposition

contains

  subroutine sync_or_abort(context)
    !! Waits for the default stream so that MPI can read the device buffers.
    !! `context` names the caller in the message printed on failure.
    implicit none

    character(len=*), intent(in) :: context

    integer :: ierr

    ierr = cudaStreamSynchronize(default_stream)
    call check_cuda(ierr, 'CUDA synchronisation before the '//context)

  end subroutine sync_or_abort

  subroutine create_fft_plan(plan, use_cufftmp, nx, ny, nz, &
                             plan_type, plan_name)
    !! Helper subroutine to create FFT plan with automatic cuFFTMp fallback
    implicit none

    integer, intent(inout) :: plan
    logical, intent(inout) :: use_cufftmp
    integer, intent(in) :: nx, ny, nz, plan_type
    character(*), intent(in) :: plan_name

    integer :: ierr
    integer(int_ptr_kind()) :: worksize
    logical :: cufftmp_failed

    cufftmp_failed = .false.
    ierr = cufftCreate(plan)

#ifndef MPI
    ! cuFFTMp distributes a transform over MPI ranks, so a build without MPI
    ! links plain cuFFT instead (see src/CMakeLists.txt) and has no cuFFTMp
    ! entry points to call.  A serial build is one rank, which plain cuFFT
    ! handles on its own.
    use_cufftmp = .false.
#endif

#ifdef MPI
    if (use_cufftmp) then
      ! Try to attach MPI communicator for cuFFTMp
      ierr = cufftMpAttachComm(plan, CUFFT_COMM_MPI, MPI_COMM_WORLD)
      if (ierr /= 0) then
        call log_warn('cuFFTMp: MPI attach failed (error '//str(ierr)//')')
        cufftmp_failed = .true.
      else
        ! MPI attach succeeded, create the plan
        ierr = cufftMakePlan3D(plan, nz, ny, nx, plan_type, worksize)
        if (ierr /= 0) then
          call log_warn('cuFFTMp: Plan creation failed (error ' &
                        //str(ierr)//')')
          cufftmp_failed = .true.
        end if
      end if

      ! if cuFFTMp failed at any stage, fall back to cuFFT
      if (cufftmp_failed) then
        call log_warn('Falling back to cuFFT')
        use_cufftmp = .false.
        ierr = cufftDestroy(plan)
        ierr = cufftCreate(plan)
      end if
    end if
#endif

    ! create plan with cuFFT
    if (.not. use_cufftmp) then
      ierr = cufftMakePlan3D(plan, nz, ny, nx, plan_type, worksize)
      call check_cufft(ierr, trim(plan_name)//' 3D FFT plan generation')
    end if

  end subroutine create_fft_plan

  subroutine create_cuda_poisson_fft(poisson_fft, mesh, xdirps, ydirps, &
                                     zdirps, lowmem)
    !! Allocates the solver type of the mesh's boundary condition case and
    !! initialises it.
    implicit none

    class(poisson_fft_t), pointer, intent(out) :: poisson_fft
    type(mesh_t), target, intent(in) :: mesh
    type(dirps_t), intent(in) :: xdirps, ydirps, zdirps
    logical, optional, intent(in) :: lowmem

    select case (poisson_bc_code(mesh%grid%periodic_BC))
    case ('000')
      allocate (cuda_poisson_000_t :: poisson_fft)
    case ('010')
      allocate (cuda_poisson_010_t :: poisson_fft)
    case ('100')
      allocate (cuda_poisson_100_t :: poisson_fft)
    case ('110')
      allocate (cuda_poisson_110_t :: poisson_fft)
    case default
      call abort_run("Unsupported periodicity combination!!")
    end select

    select type (poisson_fft)
    class is (cuda_poisson_fft_t)
      call poisson_fft%init(mesh, xdirps, ydirps, zdirps, lowmem)
    end select

  end subroutine create_cuda_poisson_fft

  subroutine init_cuda(self, mesh, xdirps, ydirps, zdirps, lowmem)
    implicit none

    class(cuda_poisson_fft_t) :: self
    type(mesh_t), target, intent(in) :: mesh
    type(dirps_t), intent(in) :: xdirps, ydirps, zdirps
    logical, optional, intent(in) :: lowmem

    integer :: nx, ny, nz

    integer :: dims_glob(3), dims_loc(3), n_spec(3), n_sp_st(3)

    integer :: fw_plan_type, bw_plan_type

    ! 1D decomposition along Z in real domain, and along Y in spectral space
    if (mesh%par%nproc_dir(2) /= 1) &
      call log_warn('nproc_dir in y-dir must be 1')
    call self%check_decomposition(mesh)

    ! Work out the spectral dimensions in the permuted state
    dims_glob = mesh%get_global_dims(CELL)
    dims_loc = mesh%get_dims(CELL)
    call self%spectral_extents(mesh, dims_glob, dims_loc, n_spec, n_sp_st)

    call self%base_init(mesh, xdirps, ydirps, zdirps, n_spec, n_sp_st)

    nx = self%nx_glob
    ny = self%ny_glob
    nz = self%nz_glob

    if (is_sp) then
      fw_plan_type = CUFFT_R2C
      bw_plan_type = CUFFT_C2R
    else
      fw_plan_type = CUFFT_D2Z
      bw_plan_type = CUFFT_Z2D
    end if

    allocate (self%waves_dev(self%nx_spec, self%ny_spec, self%nz_spec))
    self%waves_dev = self%waves

    allocate (self%ax_dev(nx), self%bx_dev(nx))
    allocate (self%ay_dev(ny), self%by_dev(ny))
    allocate (self%az_dev(nz), self%bz_dev(nz))
    self%ax_dev = self%ax; self%bx_dev = self%bx
    self%ay_dev = self%ay; self%by_dev = self%by
    self%az_dev = self%az; self%bz_dev = self%bz

    ! will store the a matrix coefficients in GPU memory if (.not. lowmem)
    ! and do a device-to-device copy at each iter. Otherwise copy from host.
    ! lowmem is .false. by default
    if (present(lowmem)) self%lowmem = lowmem

    ! Try cuFFTMp first, with automatic fallback to cuFFT if not supported
    self%use_cufftmp = .true.

    call self%init_case(mesh)

    call self%create_plans(mesh, dims_loc, fw_plan_type, bw_plan_type)

    ! use_cufftmp is only settled once create_fft_plan has run, because plain
    ! cuFFT is a fallback taken inside it, so unlike the decomposition checks
    ! above this one cannot move ahead of the plans. The dim2 mirror exists
    ! exactly on the multi-rank paths of the 100 and 110 cases.
    if (allocated(self%c_mirror_dev) .and. (.not. self%use_cufftmp)) then
      call abort_run('The 100 and 110 cases on multiple ranks need cuFFTMp, &
                     &plain cuFFT cannot decompose the transform.')
    end if

  end subroutine init_cuda

  subroutine check_decomposition(self, mesh)
    !! Stops on a decomposition the case cannot run on. The 000 case runs on
    !! any 1D z decomposition.
    implicit none

    class(cuda_poisson_fft_t) :: self
    type(mesh_t), intent(in) :: mesh
  end subroutine check_decomposition

  subroutine require_z_decomposition(mesh)
    !! The 100 and 110 cases on multiple ranks run cuFFTMp on a 1D z
    !! decomposition.
    implicit none

    type(mesh_t), intent(in) :: mesh

    if (mesh%par%nproc > 1) then
      if (mesh%par%nproc_dir(1) /= 1 .or. mesh%par%nproc_dir(2) /= 1) then
        call abort_run('The 100 and 110 cases on multiple ranks need a 1D z &
                       &decomposition, nproc_dir must be [1, 1, nproc].')
      end if
    end if

  end subroutine require_z_decomposition

  subroutine check_decomposition_010(self, mesh)
    implicit none

    class(cuda_poisson_010_t) :: self
    type(mesh_t), intent(in) :: mesh

    if (mesh%par%nproc > 1) then
      call abort_run('Multiple ranks are not yet supported for the 010 case &
                     &in the CUDA backend!')
    end if

  end subroutine check_decomposition_010

  subroutine check_decomposition_100(self, mesh)
    implicit none

    class(cuda_poisson_100_t) :: self
    type(mesh_t), intent(in) :: mesh

    integer :: dims_loc(3)

    call require_z_decomposition(mesh)

    dims_loc = mesh%get_dims(CELL)
    if (mesh%par%nproc > 1) then
      if (mod(dims_loc(1), mesh%par%nproc_dir(3)) /= 0) then
        call abort_run('The 100 case on multiple ranks needs the number of x &
                       &cells to divide by the number of ranks.')
      end if
    end if

  end subroutine check_decomposition_100

  subroutine check_decomposition_110(self, mesh)
    implicit none

    class(cuda_poisson_110_t) :: self
    type(mesh_t), intent(in) :: mesh

    integer :: dims_loc(3)

    call require_z_decomposition(mesh)

    dims_loc = mesh%get_dims(CELL)
    if (mesh%par%nproc > 1) then
      if (mod(dims_loc(1), mesh%par%nproc_dir(3)) /= 0 .or. &
          mod(dims_loc(2), mesh%par%nproc_dir(3)) /= 0) then
        call abort_run('The 110 case on multiple ranks needs the number of x &
                       &and y cells to divide by the number of ranks.')
      end if
    end if

  end subroutine check_decomposition_110

  subroutine spectral_extents(self, mesh, dims_glob, dims_loc, n_spec, &
                              n_sp_st)
    !! Local extents and offsets of the spectral slab, whose dimensions are
    !! laid out as spectral_dirs gives. dim1 holds the r2c modes and dim2 is
    !! split across ranks, as cuFFTMp leaves it, while dim3 is z, whole.
    !! This serves the 000, 010 and 100 cases.
    implicit none

    class(cuda_poisson_fft_t) :: self
    type(mesh_t), intent(in) :: mesh
    integer, intent(in) :: dims_glob(3), dims_loc(3)
    integer, intent(out) :: n_spec(3), n_sp_st(3)

    integer :: dirs(3)

    dirs = spectral_dirs(poisson_bc_code(mesh%grid%periodic_BC))

    n_spec(1) = dims_loc(dirs(1))/2 + 1
    n_spec(2) = dims_loc(dirs(2))/mesh%par%nproc_dir(3)
    n_spec(3) = dims_glob(3)

    n_sp_st(1) = 0
    n_sp_st(2) = dims_loc(dirs(2))/mesh%par%nproc_dir(3) &
                 *mesh%par%nrank_dir(3)
    n_sp_st(3) = 0

  end subroutine spectral_extents

  subroutine spectral_extents_110(self, mesh, dims_glob, dims_loc, n_spec, &
                                  n_sp_st)
    !! R2C Z-transpose: spectral array is (nz/2+1, nx, ny) on 1 rank.
    !! dim1 = Z R2C modes, dim2 = X spectrum, dim3 = Y full spectrum.
    !! On P ranks the cuFFTMp Y-slab transform shuffles the spectral
    !! output to split dim2 (X) instead of leaving it whole, same as the
    !! 100 case's dim2 split.
    implicit none

    class(cuda_poisson_110_t) :: self
    type(mesh_t), intent(in) :: mesh
    integer, intent(in) :: dims_glob(3), dims_loc(3)
    integer, intent(out) :: n_spec(3), n_sp_st(3)

    n_spec(1) = dims_glob(3)/2 + 1
    if (mesh%par%nproc > 1) then
      n_spec(2) = dims_loc(1)/mesh%par%nproc_dir(3)
      n_sp_st(2) = dims_loc(1)/mesh%par%nproc_dir(3)*mesh%par%nrank_dir(3)
    else
      n_spec(2) = dims_loc(1)
      n_sp_st(2) = 0
    end if
    n_spec(3) = dims_loc(2)

    n_sp_st(1) = 0
    n_sp_st(3) = 0

  end subroutine spectral_extents_110

  subroutine init_case(self, mesh)
    !! Sets up the device state of the case ahead of the FFT plans. The 000
    !! case needs nothing beyond the shared state.
    implicit none

    class(cuda_poisson_fft_t) :: self
    type(mesh_t), intent(in) :: mesh
  end subroutine init_case

  subroutine init_case_010(self, mesh)
    implicit none

    class(cuda_poisson_010_t) :: self
    type(mesh_t), intent(in) :: mesh

    ! if stretching in y is 'centred' or 'top-bottom'
    if (self%stretched_y .and. self%stretched_y_sym) then
      self%a_odd_re_dev = self%a_odd_re
      self%a_odd_im_dev = self%a_odd_im
      self%a_even_re_dev = self%a_even_re
      self%a_even_im_dev = self%a_even_im
      if (.not. self%lowmem) then
        self%store_a_odd_re_dev = self%a_odd_re
        self%store_a_odd_im_dev = self%a_odd_im
        self%store_a_even_re_dev = self%a_even_re
        self%store_a_even_im_dev = self%a_even_im
      end if
    !! if stretching in y is 'bottom'
    else if (self%stretched_y .and. (.not. self%stretched_y_sym)) then
      self%a_re_dev = self%a_re
      self%a_im_dev = self%a_im
      if (.not. self%lowmem) then
        self%store_a_re_dev = self%a_re
        self%store_a_im_dev = self%a_im
      end if
    end if

  end subroutine init_case_010

  subroutine init_case_100(self, mesh)
    implicit none

    class(cuda_poisson_100_t) :: self
    type(mesh_t), intent(in) :: mesh

    call self%init_dim2_mirror(mesh)

  end subroutine init_case_100

  subroutine init_case_110(self, mesh)
    implicit none

    class(cuda_poisson_110_t) :: self
    type(mesh_t), intent(in) :: mesh

    call self%init_dim2_mirror(mesh)

  end subroutine init_case_110

  subroutine init_dim2_mirror(self, mesh)
    !! The 100 and 110 cases pair each dim2 spectral mode with its mirror,
    !! which cuFFTMp may place on another rank. Slab index r of P owns global
    !! dim2 modes [r*m + 1, (r + 1)*m] with m = ny_spec, and the mirror of
    !! that range is slab P - 1 - r shifted up by one plane, so the extra
    !! plane comes from slab P - r. Both partners are symmetric. A single
    !! rank holds every partner itself and allocates nothing.
    implicit none

    class(cuda_poisson_fft_t) :: self
    type(mesh_t), intent(in) :: mesh

    if (mesh%par%nproc == 1) return

    self%mirror_slab_rank = mesh%par%nproc_dir(3) - 1 - mesh%par%nrank_dir(3)
    self%mirror_plane_rank = mod(mesh%par%nproc_dir(3) &
                                 - mesh%par%nrank_dir(3), &
                                 mesh%par%nproc_dir(3))

    allocate (self%c_mirror_dev(self%nx_spec, self%ny_spec, self%nz_spec))
    allocate (self%c_slab_send_dev(self%nx_spec, self%ny_spec, self%nz_spec))
    allocate (self%c_plane_send_dev(self%nx_spec, self%nz_spec))
    allocate (self%c_plane_recv_dev(self%nx_spec, self%nz_spec))

  end subroutine init_dim2_mirror

  subroutine create_plans(self, mesh, dims_loc, fw_plan_type, bw_plan_type)
    !! Standard R2C/C2R plans on the grid permuted into the spectral layout,
    !! with cuFFTMp fallback. This serves the 000, 010 and 100 cases.
    implicit none

    class(cuda_poisson_fft_t) :: self
    type(mesh_t), intent(in) :: mesh
    integer, intent(in) :: dims_loc(3), fw_plan_type, bw_plan_type

    integer :: ierr, dirs(3), n_glob(3)
    logical :: fw_was_cufftmp

    integer :: fft_n1, fft_n2, fft_n3, fft_n1_loc, fft_n2_loc

    ! 000, 010: FFT on (nx, ny, nz); 100: transpose X<->Y, FFT on (ny, nx, nz)
    dirs = spectral_dirs(self%bc_code)
    n_glob = [self%nx_glob, self%ny_glob, self%nz_glob]
    fft_n1 = n_glob(dirs(1)); fft_n2 = n_glob(dirs(2))
    fft_n3 = n_glob(dirs(3))
    fft_n1_loc = dims_loc(dirs(1)); fft_n2_loc = dims_loc(dirs(2))

    call create_fft_plan(self%plan3D_fw, self%use_cufftmp, &
                         fft_n1, fft_n2, fft_n3, fw_plan_type, 'Forward')
    fw_was_cufftmp = self%use_cufftmp

    ! Create backward FFT plan with automatic cuFFTMp detection/fallback
    call create_fft_plan(self%plan3D_bw, self%use_cufftmp, &
                         fft_n1, fft_n2, fft_n3, bw_plan_type, 'Backward')

    ! If backward plan forced fallback, rebuild forward plan too
    if (fw_was_cufftmp .and. (.not. self%use_cufftmp)) then
      ierr = cufftDestroy(self%plan3D_fw)
      call create_fft_plan(self%plan3D_fw, self%use_cufftmp, &
                           fft_n1, fft_n2, fft_n3, fw_plan_type, 'Forward')
    end if

    ! Allocate storage - cuFFTMp uses xtdesc, cuFFT uses c_dev
    if (self%use_cufftmp) then
      ierr = cufftXtMalloc(self%plan3D_fw, self%xtdesc, &
                           CUFFT_XT_FORMAT_INPLACE)
      call check_cufft(ierr, 'cufftXtMalloc')
    else
      ! For R2C: nx_spec is nx/2+1 (or ny/2+1 for 100 case)
      allocate (self%c_dev(self%nx_spec, self%ny_spec, self%nz_spec))
      allocate (self%r_dev(fft_n1_loc, fft_n2_loc, self%nz_loc))
    end if

    if (self%use_cufftmp) then
      call log_info('Using cuFFTMp for FFT')
    else
      call log_info('Using cuFFT for FFT')
    end if

  end subroutine create_plans

  subroutine create_plans_110(self, mesh, dims_loc, fw_plan_type, &
                              bw_plan_type)
    !! R2C with Z-transpose, on (nz, nx, ny). A single rank keeps the plain
    !! cuFFT path; on multiple ranks cuFFTMp distributes the transform, and
    !! the physical Z-slab is redistributed into its Y-slab input through an
    !! MPI_Alltoall (see redistribute_110_to_yslab).
    implicit none

    class(cuda_poisson_110_t) :: self
    type(mesh_t), intent(in) :: mesh
    integer, intent(in) :: dims_loc(3), fw_plan_type, bw_plan_type

    integer :: ierr, nx, ny, nz

    nx = self%nx_glob
    ny = self%ny_glob
    nz = self%nz_glob

    self%use_cufftmp = mesh%par%nproc > 1

    call create_fft_plan(self%plan3D_fw, self%use_cufftmp, &
                         nz, nx, ny, fw_plan_type, 'Forward 110')
    call create_fft_plan(self%plan3D_bw, self%use_cufftmp, &
                         nz, nx, ny, bw_plan_type, 'Backward 110')

    if (self%use_cufftmp) then
      ierr = cufftXtMalloc(self%plan3D_fw, self%xtdesc, &
                           CUFFT_XT_FORMAT_INPLACE)
      call check_cufft(ierr, 'cufftXtMalloc')

      ! All-to-all buffers between the physical Z-slab (nx, ny, nz_l) and
      ! the cuFFTMp Y-slab (nz, nx, ny_l), sized (nz_l, nx, ny_l, nproc).
      allocate (self%r_a2a_send_dev(self%nz_loc, nx, &
                                    ny/mesh%par%nproc_dir(3), &
                                    mesh%par%nproc_dir(3)))
      allocate (self%r_a2a_recv_dev(self%nz_loc, nx, &
                                    ny/mesh%par%nproc_dir(3), &
                                    mesh%par%nproc_dir(3)))

      call log_info('Using cuFFTMp R2C (Z-transpose) for FFT (110 case)')
    else
      ! Transposed real workspace (nz, nx, ny)
      allocate (self%r_dev_110(nz, nx, ny))
      ! Spectral complex workspace (nz/2+1, nx, ny)
      allocate (self%c_dev(self%nx_spec, self%ny_spec, self%nz_spec))

      call log_info('Using cuFFT R2C (Z-transpose) for FFT (110 case)')
    end if

  end subroutine create_plans_110

  subroutine solve_poisson_000_cuda(self, f, temp)
    !! The 000 case has nothing to reorder, so it transforms f in place.
    implicit none

    class(cuda_poisson_000_t) :: self
    class(field_t), intent(inout) :: f, temp

    call self%fft_forward(f)
    call self%fft_postprocess
    call self%fft_backward(f)

  end subroutine solve_poisson_000_cuda

  subroutine yslab_view_110(self, d_dev)
    !! Points d_dev at the cuFFTMp Y-slab real view (2*(nz/2+1), nx, ny_l)
    !! of the descriptor allocation.
    implicit none

    class(cuda_poisson_110_t), intent(in) :: self
    real(dp), device, pointer, dimension(:, :, :), intent(out) :: d_dev

    type(cudaXtDesc), pointer :: descriptor
    integer :: ny_l

    ny_l = self%ny_glob/self%mesh%par%nproc_dir(3)
    call c_f_pointer(self%xtdesc%descriptor, descriptor)
    call c_f_pointer(descriptor%data(1), d_dev, &
                     [2*(self%nz_glob/2 + 1), self%nx_glob, ny_l])

  end subroutine yslab_view_110

  subroutine fft_forward_110_cuda(self, f)
    !! Forward FFT for 110 case: transpose (nx,ny,nz)->(nz,nx,ny) then R2C
    implicit none

    class(cuda_poisson_110_t) :: self
    class(field_t), intent(in) :: f

    real(dp), device, pointer :: padded_dev(:, :, :), d_dev(:, :, :)
    integer :: ierr, tpb
    type(dim3) :: blocks, threads

    select type (f)
    type is (cuda_field_t)
      padded_dev => f%data_d
    end select

    if (self%use_cufftmp) then
      ! Redistribute the physical Z-slab (nx, ny, nz_l) into the cuFFTMp
      ! Y-slab (2*(nz/2+1), nx, ny_l) through an MPI_Alltoall, in place of
      ! the local transpose the single rank path uses.
      call yslab_view_110(self, d_dev)

      call self%redistribute_110_to_yslab(padded_dev, d_dev)

      ierr = cufftXtExecDescriptor(self%plan3D_fw, self%xtdesc, self%xtdesc, &
                                   CUFFT_FORWARD)
    else
      ! Transpose (nx, ny, nz) -> (nz, nx, ny)
      tpb = min(self%ny_loc, 256)
      blocks = dim3(self%nz_loc, (self%ny_loc - 1)/tpb + 1, 1)
      threads = dim3(tpb, 1, 1)

      call transpose_xyz_to_zxy<<<blocks, threads>>>( & !&
        self%r_dev_110, padded_dev, self%nx_loc, self%ny_loc, self%nz_loc &
        )

      ! R2C FFT on transposed data
#ifdef SINGLE_PREC
      ierr = cufftExecR2C_C(self%plan3D_fw, c_loc(self%r_dev_110), &
                            c_loc(self%c_dev))
#else
      ierr = cufftExecD2Z(self%plan3D_fw, self%r_dev_110, self%c_dev)
#endif
    end if

    call check_cufft(ierr, 'Forward R2C FFT (110 Z-transpose case)')

  end subroutine fft_forward_110_cuda

  subroutine fft_backward_110_cuda(self, f)
    !! Backward FFT for 110 case: C2R then transpose (nz,nx,ny)->(nx,ny,nz)
    implicit none

    class(cuda_poisson_110_t) :: self
    class(field_t), intent(inout) :: f

    real(dp), device, pointer :: padded_dev(:, :, :), d_dev(:, :, :)
    integer :: ierr, tpb
    type(dim3) :: blocks, threads

    select type (f)
    type is (cuda_field_t)
      padded_dev => f%data_d
    end select

    if (self%use_cufftmp) then
      ierr = cufftXtExecDescriptor(self%plan3D_bw, self%xtdesc, self%xtdesc, &
                                   CUFFT_INVERSE)
    else
      ! C2R FFT
#ifdef SINGLE_PREC
      ierr = cufftExecC2R_C(self%plan3D_bw, c_loc(self%c_dev), &
                            c_loc(self%r_dev_110))
#else
      ierr = cufftExecZ2D(self%plan3D_bw, self%c_dev, self%r_dev_110)
#endif
    end if

    call check_cufft(ierr, 'Backward C2R FFT (110 Z-transpose case)')

    ! ensure untouched padded/halo entries are deterministic
    padded_dev = 0._dp

    if (self%use_cufftmp) then
      ! Redistribute the cuFFTMp Y-slab back into the physical Z-slab
      ! (nx, ny, nz_l) through an MPI_Alltoall, in place of the local
      ! transpose the single rank path uses.
      call yslab_view_110(self, d_dev)

      call self%redistribute_110_to_zslab(d_dev, padded_dev)
    else
      ! Transpose (nz, nx, ny) -> (nx, ny, nz)
      tpb = min(self%ny_loc, 256)
      blocks = dim3(self%nz_loc, (self%ny_loc - 1)/tpb + 1, 1)
      threads = dim3(tpb, 1, 1)

      call transpose_zxy_to_xyz<<<blocks, threads>>>( & !&
        padded_dev, self%r_dev_110, self%nx_loc, self%ny_loc, self%nz_loc &
        )
    end if
  end subroutine fft_backward_110_cuda

  subroutine redistribution_launch_config(ny_l, nz_l, blocks, threads)
    !! Launch configuration shared by the four pack and unpack kernels of the
    !! 110 slab redistribution.
    implicit none

    integer, intent(in) :: ny_l, nz_l
    type(dim3), intent(out) :: blocks, threads

    integer :: tpb

    tpb = min(ny_l, 256)
    blocks = dim3(nz_l, (ny_l - 1)/tpb + 1, 1)
    threads = dim3(tpb, 1, 1)
  end subroutine redistribution_launch_config

  subroutine redistribution_alltoall(self, ny_l, slab)
    !! Syncs the device, then swaps the packed 110 slab buffers between all
    !! ranks. The slab text ('Y' or 'Z') names the target in the failure
    !! message.
    implicit none

    class(cuda_poisson_110_t) :: self
    integer, intent(in) :: ny_l
    character(len=*), intent(in) :: slab

    integer :: nx, nz_l, count, ierr_mpi

    nx = self%nx_glob
    nz_l = self%nz_loc

    ! Ordinary device memory, but still owed a sync before MPI reads it, for
    ! the same reason as the dim2 mirror exchange (see exchange_dim2_mirror).
    call sync_or_abort('110 slab redistribution')

    count = nz_l*nx*ny_l

#ifdef MPI
    call MPI_Alltoall(self%r_a2a_send_dev, count, MPI_X3D2_DP, &
                      self%r_a2a_recv_dev, count, MPI_X3D2_DP, &
                      MPI_COMM_WORLD, ierr_mpi)
#else
    call abort_run('The 110 slab redistribution needs more than one rank, but &
                   &this build was configured without MPI')
#endif
    call check_mpi(ierr_mpi, 'The 110 slab redistribution to the ' &
                   //slab//'-slab')
  end subroutine redistribution_alltoall

  subroutine redistribute_110_to_yslab(self, zslab, yslab)
    !! Redistributes the 110 field from the physical Z-slab (nx, ny, nz_l)
    !! to the cuFFTMp Y-slab real view (2*(nz/2+1), nx, ny_l), ahead of the
    !! forward transform. See the comment above pack_110_xyz_to_yslab.
    implicit none

    class(cuda_poisson_110_t) :: self
    real(dp), device, dimension(:, :, :), intent(in) :: zslab
    real(dp), device, dimension(:, :, :), intent(inout) :: yslab

    type(dim3) :: blocks, threads
    integer :: nx, ny_l, nz_l, nproc

    nx = self%nx_glob
    nz_l = self%nz_loc
    nproc = self%mesh%par%nproc_dir(3)
    ny_l = size(yslab, 3)

    call redistribution_launch_config(ny_l, nz_l, blocks, threads)

    call pack_110_xyz_to_yslab<<<blocks, threads>>>( & !&
      self%r_a2a_send_dev, zslab, nx, ny_l, nz_l, nproc &
      )

    call redistribution_alltoall(self, ny_l, 'Y')

    call unpack_110_yslab<<<blocks, threads>>>( & !&
      yslab, self%r_a2a_recv_dev, nx, ny_l, nz_l, nproc &
      )

  end subroutine redistribute_110_to_yslab

  subroutine redistribute_110_to_zslab(self, yslab, zslab)
    !! Inverse of redistribute_110_to_yslab, after the backward transform.
    implicit none

    class(cuda_poisson_110_t) :: self
    real(dp), device, dimension(:, :, :), intent(in) :: yslab
    real(dp), device, dimension(:, :, :), intent(inout) :: zslab

    type(dim3) :: blocks, threads
    integer :: nx, ny_l, nz_l, nproc

    nx = self%nx_glob
    nz_l = self%nz_loc
    nproc = self%mesh%par%nproc_dir(3)
    ny_l = size(yslab, 3)

    call redistribution_launch_config(ny_l, nz_l, blocks, threads)

    call pack_110_yslab<<<blocks, threads>>>( & !&
      self%r_a2a_send_dev, yslab, nx, ny_l, nz_l, nproc &
      )

    call redistribution_alltoall(self, ny_l, 'Z')

    call unpack_110_yslab_to_xyz<<<blocks, threads>>>( & !&
      zslab, self%r_a2a_recv_dev, nx, ny_l, nz_l, nproc &
      )

  end subroutine redistribute_110_to_zslab

  subroutine fft_forward_100_cuda(self, f)
  !! Forward FFT for non-periodic-X case
  !! We transpose X<->Y so that the non-periodic direction becomes the
  !! "Y" direction in the transposed space, then use the same FFT approach as 010
    implicit none

    class(cuda_poisson_100_t) :: self
    class(field_t), intent(in) :: f

    real(dp), device, pointer :: padded_dev(:, :, :), d_dev(:, :, :)
    type(cudaXtDesc), pointer :: descriptor

    integer :: tsize, ierr
    type(dim3) :: blocks, threads

    select type (f)
    type is (cuda_field_t)
      padded_dev => f%data_d
    end select

    if (self%use_cufftmp) then
      call c_f_pointer(self%xtdesc%descriptor, descriptor)

      ! For 100 case with transposed FFT plan (ny, nx, nz):
      ! The descriptor has shape (2*(ny/2 +1), nx, nz) for R2C output (see: issue #279)
      call c_f_pointer(descriptor%data(1), d_dev, &
                       [2*(self%ny_loc/2 + 1), self%nx_loc, self%nz_loc])

      tsize = 16
      blocks = dim3((self%nx_loc - 1)/tsize + 1, self%nz_loc, 1)
      threads = dim3(tsize, 1, 1)

      ! Copy with transpose: src(nx,ny,nz) -> dst(ny,nx,nz)
      call memcpy3D_with_transpose<<<blocks, threads>>>( & !&
        d_dev, padded_dev, self%nx_loc, self%ny_loc, self%nz_loc &
        )

      ierr = cufftXtExecDescriptor(self%plan3D_fw, self%xtdesc, self%xtdesc, &
                                   CUFFT_FORWARD)
    else
      ! using standard cuFFT with explicit real workspace
      call c_f_pointer(c_loc(self%r_dev), d_dev, &
                       [self%ny_loc, self%nx_loc, self%nz_loc])

      tsize = 16
      blocks = dim3((self%nx_loc - 1)/tsize + 1, self%nz_loc, 1)
      threads = dim3(tsize, 1, 1)

      ! copy with transpose: src(nx,ny,nz) -> dst(ny,nx,nz)
      call memcpy3D_with_transpose<<<blocks, threads>>>( & !&
        d_dev, padded_dev, self%nx_loc, self%ny_loc, self%nz_loc &
        )

#ifdef SINGLE_PREC
      ierr = cufftExecR2C_C(self%plan3D_fw, c_loc(self%r_dev), &
                            c_loc(self%c_dev))
#else
      ierr = cufftExecD2Z(self%plan3D_fw, self%r_dev, self%c_dev)
#endif
    end if

    call check_cufft(ierr, 'Forward 3D FFT execution (100 case)')

  end subroutine fft_forward_100_cuda

  subroutine fft_backward_100_cuda(self, f)
    implicit none

    class(cuda_poisson_100_t) :: self
    class(field_t), intent(inout) :: f

    real(dp), device, pointer :: padded_dev(:, :, :), d_dev(:, :, :)
    type(cudaXtDesc), pointer :: descriptor

    integer :: tsize, ierr
    type(dim3) :: blocks, threads

    select type (f)
    type is (cuda_field_t)
      padded_dev => f%data_d
    end select

    if (self%use_cufftmp) then
      ierr = cufftXtExecDescriptor(self%plan3D_bw, self%xtdesc, self%xtdesc, &
                                   CUFFT_INVERSE)
    else
      ! using standard cuFFT with explicit real workspace
#ifdef SINGLE_PREC
      ierr = cufftExecC2R_C(self%plan3D_bw, c_loc(self%c_dev), &
                            c_loc(self%r_dev))
#else
      ierr = cufftExecZ2D(self%plan3D_bw, self%c_dev, self%r_dev)
#endif
    end if

    call check_cufft(ierr, 'Backward 3D FFT execution (100 case)')

    if (self%use_cufftmp) then
      call c_f_pointer(self%xtdesc%descriptor, descriptor)
      call c_f_pointer(descriptor%data(1), d_dev, &
                       [2*(self%ny_loc/2 + 1), self%nx_loc, self%nz_loc])

      tsize = 16
      blocks = dim3((self%nx_loc - 1)/tsize + 1, self%nz_loc, 1)
      threads = dim3(tsize, 1, 1)

      ! Copy with transpose back: src(ny,nx,nz) -> dst(nx,ny,nz)
      call memcpy3D_with_transpose_back<<<blocks, threads>>>( & !&
        padded_dev, d_dev, self%nx_loc, self%ny_loc, self%nz_loc &
        )
    else
      ! ensure untouched padded/halo entries are deterministic
      padded_dev = 0._dp

      ! cast r_dev as (ny, nx, nz) real buffer
      call c_f_pointer(c_loc(self%r_dev), d_dev, &
                       [self%ny_loc, self%nx_loc, self%nz_loc])

      tsize = 16
      blocks = dim3((self%nx_loc - 1)/tsize + 1, self%nz_loc, 1)
      threads = dim3(tsize, 1, 1)

      ! copy with transpose back: src(ny,nx,nz) -> dst(nx,ny,nz)
      call memcpy3D_with_transpose_back<<<blocks, threads>>>( & !&
        padded_dev, d_dev, self%nx_loc, self%ny_loc, self%nz_loc &
        )
    end if

  end subroutine fft_backward_100_cuda

  subroutine fft_forward_cuda(self, f)
    implicit none

    class(cuda_poisson_fft_t) :: self
    class(field_t), intent(in) :: f

    real(dp), device, pointer :: padded_dev(:, :, :), d_dev(:, :, :)

    type(cudaXtDesc), pointer :: descriptor

    integer :: tsize, ierr
    type(dim3) :: blocks, threads

    select type (f)
    type is (cuda_field_t)
      padded_dev => f%data_d
    end select

    if (self%use_cufftmp) then
      ! using cuFFTMp
      ! tsize is different than SZ, because here we work on a 3D Cartesian
      ! data structure, and free to specify any suitable thread/block size.
      tsize = 16
      blocks = dim3((self%ny_loc - 1)/tsize + 1, self%nz_loc, 1)
      threads = dim3(tsize, 1, 1)

      call c_f_pointer(self%xtdesc%descriptor, descriptor)
      call c_f_pointer(descriptor%data(1), d_dev, &
                       [2*(self%nx_loc/2 + 1), self%ny_loc, self%nz_loc])

      call memcpy3D<<<blocks, threads>>>( & !&
        d_dev, padded_dev, self%nx_loc, self%ny_loc, self%nz_loc &
        )

      ierr = cufftXtExecDescriptor(self%plan3D_fw, self%xtdesc, self%xtdesc, &
                                   CUFFT_FORWARD)
    else
      ! using standard cuFFT with explicit real workspace
      call c_f_pointer(c_loc(self%r_dev), d_dev, &
                       [self%nx_loc, self%ny_loc, self%nz_loc])

      tsize = 16
      blocks = dim3((self%ny_loc - 1)/tsize + 1, self%nz_loc, 1)
      threads = dim3(tsize, 1, 1)

      call memcpy3D<<<blocks, threads>>>( & !&
        d_dev, padded_dev, self%nx_loc, self%ny_loc, self%nz_loc &
        )

#ifdef SINGLE_PREC
      ierr = cufftExecR2C_C(self%plan3D_fw, c_loc(self%r_dev), &
                            c_loc(self%c_dev))
#else
      ierr = cufftExecD2Z(self%plan3D_fw, self%r_dev, self%c_dev)
#endif
    end if

    call check_cufft(ierr, 'Forward 3D FFT execution')

  end subroutine fft_forward_cuda

  subroutine fft_backward_cuda(self, f)
    implicit none

    class(cuda_poisson_fft_t) :: self
    class(field_t), intent(inout) :: f

    real(dp), device, pointer :: padded_dev(:, :, :), d_dev(:, :, :)

    type(cudaXtDesc), pointer :: descriptor

    integer :: tsize, ierr
    type(dim3) :: blocks, threads

    select type (f)
    type is (cuda_field_t)
      padded_dev => f%data_d
    end select

    if (self%use_cufftmp) then
      ! using cuFFTMp
      ierr = cufftXtExecDescriptor(self%plan3D_bw, self%xtdesc, self%xtdesc, &
                                   CUFFT_INVERSE)
    else
      ! using standard cuFFT with explicit real workspace
#ifdef SINGLE_PREC
      ierr = cufftExecC2R_C(self%plan3D_bw, c_loc(self%c_dev), &
                            c_loc(self%r_dev))
#else
      ierr = cufftExecZ2D(self%plan3D_bw, self%c_dev, self%r_dev)
#endif
    end if

    call check_cufft(ierr, 'Backward 3D FFT execution')

    if (self%use_cufftmp) then
      ! cuFFTMp path: cast xtdesc as (nx+2, ny, nz) real buffer
      call c_f_pointer(self%xtdesc%descriptor, descriptor)
      call c_f_pointer(descriptor%data(1), d_dev, &
                       [2*(self%nx_loc/2 + 1), self%ny_loc, self%nz_loc])
    else
      ! ensure untouched padded/halo entries are deterministic
      padded_dev = 0._dp

      ! cuFFT path: cast r_dev as (nx, ny, nz) real buffer
      call c_f_pointer(c_loc(self%r_dev), d_dev, &
                       [self%nx_loc, self%ny_loc, self%nz_loc])
    end if

    tsize = 16
    blocks = dim3((self%ny_loc - 1)/tsize + 1, self%nz_loc, 1)
    threads = dim3(tsize, 1, 1)
    call memcpy3D<<<blocks, threads>>>( & !&
      padded_dev, d_dev, self%nx_loc, self%ny_loc, self%nz_loc &
      )

  end subroutine fft_backward_cuda

  subroutine fft_postprocess_000_cuda(self)
    implicit none

    class(cuda_poisson_000_t) :: self

    type(cudaXtDesc), pointer :: descriptor

    complex(dp), device, dimension(:, :, :), pointer :: c_dev
    type(dim3) :: blocks, threads
    integer :: tsize

    ! tsize is different than SZ, because here we work on a 3D Cartesian
    ! data structure, and free to specify any suitable thread/block size.
    tsize = 16
    blocks = dim3((self%ny_spec - 1)/tsize + 1, self%nz_spec, 1)
    threads = dim3(tsize, 1, 1)

    ! Get pointer to the appropriate FFT data storage
    if (self%use_cufftmp) then
      ! cuFFTMp path: get pointer from xtdesc
      call c_f_pointer(self%xtdesc%descriptor, descriptor)
      call c_f_pointer(descriptor%data(1), c_dev, &
                       [self%nx_spec, self%ny_spec, self%nz_spec])
    else
      ! cuFFT path: get pointer to self%c_dev
      call c_f_pointer(c_loc(self%c_dev), c_dev, &
                       [self%nx_spec, self%ny_spec, self%nz_spec])
    end if

    ! Postprocess div_u in spectral space
    call process_spectral_000<<<blocks, threads>>>( & !&
      c_dev, self%waves_dev, self%nx_spec, self%ny_spec, self%sp_st(2), &
      self%nx_glob, self%ny_glob, self%nz_glob, &
      self%ax_dev, self%bx_dev, self%ay_dev, self%by_dev, &
      self%az_dev, self%bz_dev &
      )

  end subroutine fft_postprocess_000_cuda

  subroutine fft_postprocess_100_cuda(self)
    implicit none

    class(cuda_poisson_100_t) :: self

    type(cudaXtDesc), pointer :: descriptor

    complex(dp), device, dimension(:, :, :), pointer :: c_dev
    type(dim3) :: blocks, threads
    integer :: tsize

    ! Get pointer to the appropriate FFT data storage
    if (self%use_cufftmp) then
      call c_f_pointer(self%xtdesc%descriptor, descriptor)
      call c_f_pointer(descriptor%data(1), c_dev, &
                       [self%nx_spec, self%ny_spec, self%nz_spec])
    else
      call c_f_pointer(c_loc(self%c_dev), c_dev, &
                       [self%nx_spec, self%ny_spec, self%nz_spec])
    end if

    ! tsize is different than SZ, because here we work on a 3D Cartesian
    ! data structure, and free to specify any suitable thread/block size.
    tsize = 16
    blocks = dim3((self%nx_spec - 1)/tsize + 1, self%nz_spec, 1)
    threads = dim3(tsize, 1, 1)

    ! Note the spectral component names describe the 010 layout, so in
    ! the 100 case nx_spec counts Y (R2C) modes, ny_spec counts the local
    ! X modes, and sp_st(2) is the X mode offset. On a single rank the X
    ! paired split reads its partner straight out of div_u, which is what
    ! the 010 kernel already does, so it is reused with swapped arguments:
    ! - Swap nx <-> ny for grid sizes
    ! - Swap ax,bx <-> ay,by for wave coefficients
    ! - Use sp_st(2) for the dim2 offset (X modes after transpose)
    if (self%mesh%par%nproc == 1) then
      call process_spectral_010<<<blocks, threads>>>( & !&
        c_dev, self%waves_dev, &
        self%nx_spec, self%ny_spec, self%sp_st(2), &
        self%ny_glob, self%nx_glob, self%nz_glob, &  ! SWAP nx <-> ny
        self%ay_dev, self%by_dev, &  ! dim1, periodic Y
        self%ax_dev, self%bx_dev, &  ! dim2, non-periodic X
        self%az_dev, self%bz_dev &   ! dim3, Z
        )
      return
    end if

    ! On multiple ranks the X paired split reaches its partner through a
    ! mirror exchange, so the same work is staged across three launches.
    ! The same argument swap as above applies.

    ! Stage 1: normalisation, then the Y (R2C) and Z postprocess
    call process_spectral_100_fw<<<blocks, threads>>>( & !&
      c_dev, self%nx_spec, self%ny_spec, &
      self%ny_glob, self%nx_glob, self%nz_glob, &  ! SWAP nx <-> ny
      self%ay_dev, self%by_dev, &  ! dim1 uses Y coefficients
      self%az_dev, self%bz_dev &   ! dim3 uses Z coefficients
      )

    ! Fetch the X pairing partners that live on other ranks
    call self%exchange_dim2_mirror(c_dev)

    ! Stage 2: X paired split (forward), then the Poisson solve
    call process_spectral_100_pair_fw<<<blocks, threads>>>( & !&
      c_dev, self%c_mirror_dev, self%c_plane_recv_dev, self%waves_dev, &
      self%nx_spec, self%ny_spec, self%sp_st(2), &  ! dims and dim2 offset
      self%ny_glob, self%nz_glob, &  ! SWAP nx <-> ny
      self%ax_dev, self%bx_dev &  ! dim2 uses X coefficients
      )

    ! Stage 2 rewrote every mode, so the mirror is stale and is refreshed
    call self%exchange_dim2_mirror(c_dev)

    ! Stage 3: X paired recombination (backward), then the Y and Z
    ! postprocess
    call process_spectral_100_pair_bw<<<blocks, threads>>>( & !&
      c_dev, self%c_mirror_dev, self%c_plane_recv_dev, &
      self%nx_spec, self%ny_spec, self%sp_st(2), &  ! dims and dim2 offset
      self%ny_glob, self%nz_glob, &  ! SWAP nx <-> ny
      self%ay_dev, self%by_dev, &  ! dim1 uses Y coefficients
      self%ax_dev, self%bx_dev, &  ! dim2 uses X coefficients
      self%az_dev, self%bz_dev &   ! dim3 uses Z coefficients
      )
  end subroutine fft_postprocess_100_cuda

  subroutine exchange_dim2_mirror(self, c_dev)
    !! Fills c_mirror_dev and c_plane_recv_dev with the dim2 mirror of the
    !! local spectral slab, so that the paired split can read its
    !! partner mode without going looking for it in another rank's memory.
    !!
    !! Slab index r of P owns global dim2 modes [r*m + 1, (r + 1)*m], with
    !! m = ny_spec. Local mode j pairs with global mode
    !! (P - r)*m - j + 2, which is local mode m - j + 2 of slab P - 1 - r
    !! for j >= 2, and the first mode of slab P - r for j = 1. One slab
    !! sized exchange plus one plane sized exchange therefore covers the
    !! whole pairing, and both partners are symmetric so the rank that
    !! sends is the rank that receives.
    implicit none

    class(cuda_poisson_fft_t) :: self
    complex(dp), device, dimension(:, :, :), intent(in) :: c_dev

    type(dim3) :: blocks, threads
    integer :: tsize, n_slab, n_plane, mpi_cplx, nrank, ierr_mpi
    integer :: tag_slab = 2345, tag_plane = 2346

    ! MPI cannot take the first dim2 plane strided out of c_dev, so pack it
    tsize = 16
    blocks = dim3((self%nx_spec - 1)/tsize + 1, self%nz_spec, 1)
    threads = dim3(tsize, 1, 1)
    call pack_spectral_plane<<<blocks, threads>>>( & !&
  self%c_plane_send_dev, c_dev, self%nx_spec, self%nz_spec &
  )

    nrank = self%mesh%par%nrank

    ! Copy the slab out of the descriptor allocation before MPI sees it, see
    ! the comment on c_slab_send_dev. Done ahead of the synchronisation below
    ! so that one sync covers the packing kernel and this copy together.
    if (self%mirror_slab_rank /= nrank) then
      self%c_slab_send_dev = c_dev
    end if

    if (is_sp) then
      mpi_cplx = MPI_COMPLEX
    else
      mpi_cplx = MPI_DOUBLE_COMPLEX
    end if

    n_slab = self%nx_spec*self%ny_spec*self%nz_spec
    n_plane = self%nx_spec*self%nz_spec

    ! Make sure the staged kernel and the packing above have landed before
    ! MPI reads the device buffers. This is only reached on multiple ranks.
    ! fft_postprocess_100_cuda returns early for the single rank path, and
    ! fft_postprocess_110_cuda only gets here inside its use_cufftmp branch.
    call sync_or_abort('dim2 mirror exchange')

    ! On the middle rank when P is odd the mirror is the local slab itself.
    ! It is still copied out rather than read in place, so that the pairing
    ! kernels always see pre-update values.
    if (self%mirror_slab_rank == nrank) then
      self%c_mirror_dev = c_dev
    else
#ifdef MPI
      call MPI_Sendrecv(self%c_slab_send_dev, n_slab, mpi_cplx, &
                        self%mirror_slab_rank, tag_slab, &
                        self%c_mirror_dev, n_slab, mpi_cplx, &
                        self%mirror_slab_rank, tag_slab, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr_mpi)
#else
      call abort_run('The dim2 mirror slab exchange needs more than one rank, &
                     &but this build was configured without MPI')
#endif
      call check_mpi(ierr_mpi, 'The dim2 mirror slab exchange')
    end if

    if (self%mirror_plane_rank == nrank) then
      self%c_plane_recv_dev = self%c_plane_send_dev
    else
#ifdef MPI
      call MPI_Sendrecv(self%c_plane_send_dev, n_plane, mpi_cplx, &
                        self%mirror_plane_rank, tag_plane, &
                        self%c_plane_recv_dev, n_plane, mpi_cplx, &
                        self%mirror_plane_rank, tag_plane, &
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr_mpi)
#else
      call abort_run('The dim2 mirror plane exchange needs more than one &
                     &rank, but this build was configured without MPI')
#endif
      call check_mpi(ierr_mpi, 'The dim2 mirror plane exchange')
    end if

  end subroutine exchange_dim2_mirror

  subroutine fft_postprocess_010_cuda(self)
    implicit none

    class(cuda_poisson_010_t) :: self

    type(cudaXtDesc), pointer :: descriptor

    complex(dp), device, dimension(:, :, :), pointer :: c_dev
    type(dim3) :: blocks, threads
    integer :: tsize, off, inc

    ! tsize is different than SZ, because here we work on a 3D Cartesian
    ! data structure, and free to specify any suitable thread/block size.
    tsize = 16
    blocks = dim3((self%nx_spec - 1)/tsize + 1, self%nz_spec, 1)
    threads = dim3(tsize, 1, 1)

    ! Get pointer to the appropriate FFT data storage
    if (self%use_cufftmp) then
      ! cuFFTMp path: get pointer from xtdesc
      call c_f_pointer(self%xtdesc%descriptor, descriptor)
      call c_f_pointer(descriptor%data(1), c_dev, &
                       [self%nx_spec, self%ny_spec, self%nz_spec])
    else
      ! cuFFT path: get pointer to self%c_dev
      call c_f_pointer(c_loc(self%c_dev), c_dev, &
                       [self%nx_spec, self%ny_spec, self%nz_spec])
    end if

    ! Postprocess div_u in spectral space
    if (.not. self%stretched_y) then
      call process_spectral_010<<<blocks, threads>>>( & !&
        c_dev, self%waves_dev, self%nx_spec, self%ny_spec, self%sp_st(2), &
        self%nx_glob, self%ny_glob, self%nz_glob, &
        self%ax_dev, self%bx_dev, self%ay_dev, self%by_dev, &
        self%az_dev, self%bz_dev &
        )
    else
      call process_spectral_010_fw<<<blocks, threads>>>( & !&
        c_dev, self%nx_spec, self%ny_spec, self%sp_st(2), &
        self%nx_glob, self%ny_glob, self%nz_glob, &
        self%ax_dev, self%bx_dev, self%ay_dev, self%by_dev, &
        self%az_dev, self%bz_dev &
        )

      ! if stretching in y is 'centred' or 'top-bottom'
      if (self%stretched_y_sym) then
        ! copy from host to device if lowmem else from device stores
        if (self%lowmem) then
          self%a_odd_re_dev = self%a_odd_re
          self%a_odd_im_dev = self%a_odd_im
          self%a_even_re_dev = self%a_even_re
          self%a_even_im_dev = self%a_even_im
        else
          self%a_odd_re_dev = self%store_a_odd_re_dev
          self%a_odd_im_dev = self%store_a_odd_im_dev
          self%a_even_re_dev = self%store_a_even_re_dev
          self%a_even_im_dev = self%store_a_even_im_dev
        end if
        ! start from the first odd entry and continue with odd ones
        off = 0
        inc = 2
        call process_spectral_010_poisson<<<blocks, threads>>>( & !&
          c_dev, self%a_odd_re_dev, self%a_odd_im_dev, off, inc, &
          self%nx_spec, self%ny_spec/2, &
          self%nx_glob, self%ny_glob, self%nz_glob &
          )
        ! start from the first even entry and continue with even ones
        off = 1
        call process_spectral_010_poisson<<<blocks, threads>>>( & !&
          c_dev, self%a_even_re_dev, self%a_even_im_dev, off, inc, &
          self%nx_spec, self%ny_spec/2, &
          self%nx_glob, self%ny_glob, self%nz_glob &
          )
      !! if stretching in y is 'bottom'
      else
        ! copy from host to device if lowmem else from device stores
        if (self%lowmem) then
          self%a_re_dev = self%a_re
          self%a_im_dev = self%a_im
        else
          self%a_re_dev = self%store_a_re_dev
          self%a_im_dev = self%store_a_im_dev
        end if
        ! start from the first entry and increment 1
        off = 0
        inc = 1
        call process_spectral_010_poisson<<<blocks, threads>>>( & !&
          c_dev, self%a_re_dev, self%a_im_dev, off, inc, &
          self%nx_spec, self%ny_spec, &
          self%nx_glob, self%ny_glob, self%nz_glob &
          )
      end if

      call process_spectral_010_bw<<<blocks, threads>>>( & !&
        c_dev, self%nx_spec, self%ny_spec, self%sp_st(2), &
        self%nx_glob, self%ny_glob, self%nz_glob, &
        self%ax_dev, self%bx_dev, self%ay_dev, self%by_dev, &
        self%az_dev, self%bz_dev &
        )
    end if

  end subroutine fft_postprocess_010_cuda

  subroutine fft_postprocess_110_cuda(self)
    !! Spectral post-processing for 110 case (R2C Z-transpose).
    !! Spectral array is (nz/2+1, nx, ny) = (nx_spec, ny_spec, nz_spec).
    !! Thread i->X (ny_spec), blockIdx%y k->Y (nz_spec), serial j->Z (nx_spec).
    !! On a single rank X and Y are both whole, so 7 kernel launches avoid
    !! cross-block race conditions. On multiple ranks cuFFTMp splits X
    !! across ranks, so the X pairing reaches its partner through the
    !! dim2 mirror exchange (see exchange_dim2_mirror) instead, adding two
    !! more launches.
    implicit none

    class(cuda_poisson_110_t) :: self

    type(cudaXtDesc), pointer :: descriptor
    complex(dp), device, dimension(:, :, :), pointer :: c_dev
    type(dim3) :: blocks, threads
    integer :: tsize, nz_h

    ! Get pointer to the appropriate FFT data storage
    if (self%use_cufftmp) then
      call c_f_pointer(self%xtdesc%descriptor, descriptor)
      call c_f_pointer(descriptor%data(1), c_dev, &
                       [self%nx_spec, self%ny_spec, self%nz_spec])
    else
      call c_f_pointer(c_loc(self%c_dev), c_dev, &
                       [self%nx_spec, self%ny_spec, self%nz_spec])
    end if

    nz_h = self%nx_spec  ! = nz/2+1 (1st dim of spectral array)

    tsize = 16
    ! threads over X (ny_spec = local nx), blocks over Y (nz_spec = ny)
    blocks = dim3((self%ny_spec - 1)/tsize + 1, self%nz_spec, 1)
    threads = dim3(tsize, 1, 1)

    ! Step 1: normalise + Z periodic forward
    call process_spectral_110_norm_z<<<blocks, threads>>>( & !&
      c_dev, nz_h, &
      self%ny_spec, self%nz_spec, self%nz_glob, &
      self%nx_glob, self%ny_glob, &
      self%az_dev, self%bz_dev &
      )

    ! On multiple ranks X is split across ranks by cuFFTMp, reached through
    ! the dim2 mirror. Y stays whole (self%sp_st(3) is always 0), so its
    ! pair kernels are the same on one rank and on several.
    if (.not. self%use_cufftmp) then
      ! Step 2: X paired split (forward)
      call process_spectral_110_x_pair_fw<<<blocks, threads>>>( & !&
        c_dev, nz_h, &
        self%ny_spec, self%nz_spec, self%sp_st(2), &
        self%ax_dev, self%bx_dev &
        )
    else
      ! Fetch the X pairing partners that live on other ranks
      call self%exchange_dim2_mirror(c_dev)

      ! Step 2: X paired split (forward), mirror-aware
      call process_spectral_110_x_pair_fw_mirror<<<blocks, threads>>>( & !&
        c_dev, self%c_mirror_dev, self%c_plane_recv_dev, nz_h, &
        self%ny_spec, self%nz_spec, self%sp_st(2), &
        self%ax_dev, self%bx_dev &
        )
    end if

    ! Step 3: Y paired split (forward)
    call process_spectral_110_y_pair_fw<<<blocks, threads>>>( & !&
      c_dev, nz_h, &
      self%ny_spec, self%nz_spec, self%sp_st(3), &
      self%ay_dev, self%by_dev &
      )

    ! Step 4: Poisson solve
    call process_spectral_110_poisson<<<blocks, threads>>>( & !&
      c_dev, self%waves_dev, nz_h, &
      self%ny_spec, self%nz_spec, self%nz_glob, self%sp_st(2), &
      self%nx_glob &
      )

    ! Step 5: Y paired recombine (backward)
    call process_spectral_110_y_pair_bw<<<blocks, threads>>>( & !&
      c_dev, nz_h, &
      self%ny_spec, self%nz_spec, self%sp_st(3), &
      self%ay_dev, self%by_dev &
      )

    if (.not. self%use_cufftmp) then
      ! Step 6: X paired recombine (backward)
      call process_spectral_110_x_pair_bw<<<blocks, threads>>>( & !&
        c_dev, nz_h, &
        self%ny_spec, self%nz_spec, self%sp_st(2), &
        self%ax_dev, self%bx_dev &
        )
    else
      ! Steps 2 to 5 rewrote every mode, so the mirror is stale
      call self%exchange_dim2_mirror(c_dev)

      ! Step 6: X paired recombine (backward), mirror-aware
      call process_spectral_110_x_pair_bw_mirror<<<blocks, threads>>>( & !&
        c_dev, self%c_mirror_dev, self%c_plane_recv_dev, nz_h, &
        self%ny_spec, self%nz_spec, self%sp_st(2), &
        self%ax_dev, self%bx_dev &
        )
    end if

    ! Step 7: Z periodic undo (backward)
    call process_spectral_110_z_bw<<<blocks, threads>>>( & !&
      c_dev, nz_h, &
      self%ny_spec, self%nz_spec, &
      self%az_dev, self%bz_dev &
      )

  end subroutine fft_postprocess_110_cuda

  subroutine enforce_periodicity_x_cuda(self, f_out, f_in)
    implicit none

    class(cuda_poisson_100_t) :: self
    class(field_t), intent(inout) :: f_out
    class(field_t), intent(in) :: f_in

    real(dp), device, pointer, dimension(:, :, :) :: f_out_dev, f_in_dev
    type(dim3) :: blocks, threads

    select type (f_out)
    type is (cuda_field_t)
      f_out_dev => f_out%data_d
    end select
    select type (f_in)
    type is (cuda_field_t)
      f_in_dev => f_in%data_d
    end select

    blocks = dim3(self%nz_loc, 1, 1)
    threads = dim3(self%ny_loc, 1, 1)

    ! cuFFT path can consume stale padded/halo entries.
    if (.not. self%use_cufftmp) f_out_dev = 0._dp

    call enforce_periodicity_x<<<blocks, threads>>>( & !&
      f_out_dev, f_in_dev, self%nx_glob &
      )

  end subroutine enforce_periodicity_x_cuda

  subroutine undo_periodicity_x_cuda(self, f_out, f_in)
    implicit none

    class(cuda_poisson_100_t) :: self
    class(field_t), intent(inout) :: f_out
    class(field_t), intent(in) :: f_in

    real(dp), device, pointer, dimension(:, :, :) :: f_out_dev, f_in_dev
    type(dim3) :: blocks, threads

    select type (f_out)
    type is (cuda_field_t)
      f_out_dev => f_out%data_d
    end select
    select type (f_in)
    type is (cuda_field_t)
      f_in_dev => f_in%data_d
    end select

    blocks = dim3(self%nz_loc, 1, 1)
    threads = dim3(self%ny_loc, 1, 1)

    ! cuFFT path can consume stale padded/halo entries.
    if (.not. self%use_cufftmp) f_out_dev = 0._dp

    call undo_periodicity_x<<<blocks, threads>>>( & !&
      f_out_dev, f_in_dev, self%nx_glob &
      )

  end subroutine undo_periodicity_x_cuda

  subroutine enforce_periodicity_y_cuda(self, f_out, f_in)
    implicit none

    class(cuda_poisson_010_t) :: self
    class(field_t), intent(inout) :: f_out
    class(field_t), intent(in) :: f_in

    real(dp), device, pointer, dimension(:, :, :) :: f_out_dev, f_in_dev
    type(dim3) :: blocks, threads

    select type (f_out)
    type is (cuda_field_t)
      f_out_dev => f_out%data_d
    end select
    select type (f_in)
    type is (cuda_field_t)
      f_in_dev => f_in%data_d
    end select

    blocks = dim3(self%nz_loc, 1, 1)
    threads = dim3(self%nx_loc, 1, 1)

    ! cuFFT path can consume stale padded/halo entries.
    if (.not. self%use_cufftmp) f_out_dev = 0._dp

    call enforce_periodicity_y<<<blocks, threads>>>( & !&
      f_out_dev, f_in_dev, self%ny_glob &
      )

  end subroutine enforce_periodicity_y_cuda

  subroutine undo_periodicity_y_cuda(self, f_out, f_in)
    implicit none

    class(cuda_poisson_010_t) :: self
    class(field_t), intent(inout) :: f_out
    class(field_t), intent(in) :: f_in

    real(dp), device, pointer, dimension(:, :, :) :: f_out_dev, f_in_dev
    type(dim3) :: blocks, threads

    select type (f_out)
    type is (cuda_field_t)
      f_out_dev => f_out%data_d
    end select
    select type (f_in)
    type is (cuda_field_t)
      f_in_dev => f_in%data_d
    end select

    blocks = dim3(self%nz_loc, 1, 1)
    threads = dim3(self%nx_loc, 1, 1)

    ! cuFFT path can consume stale padded/halo entries.
    if (.not. self%use_cufftmp) f_out_dev = 0._dp

    call undo_periodicity_y<<<blocks, threads>>>( & !&
      f_out_dev, f_in_dev, self%ny_glob &
      )

  end subroutine undo_periodicity_y_cuda

  subroutine enforce_periodicity_xy_cuda(self, f_out, f_in)
    implicit none

    class(cuda_poisson_110_t) :: self
    class(field_t), intent(inout) :: f_out
    class(field_t), intent(in) :: f_in

    real(dp), device, pointer, dimension(:, :, :) :: f_out_dev, f_in_dev
    type(dim3) :: blocks, threads
    integer :: tpb

    select type (f_out)
    type is (cuda_field_t)
      f_out_dev => f_out%data_d
    end select
    select type (f_in)
    type is (cuda_field_t)
      f_in_dev => f_in%data_d
    end select

    f_out_dev = 0._dp

    tpb = min(self%ny_glob, 256)
    blocks = dim3(self%nz_loc, (self%ny_glob - 1)/tpb + 1, 1)
    threads = dim3(tpb, 1, 1)

    call enforce_periodicity_xy<<<blocks, threads>>>( & !&
      f_out_dev, f_in_dev, self%nx_glob, self%ny_glob &
      )

  end subroutine enforce_periodicity_xy_cuda

  subroutine undo_periodicity_xy_cuda(self, f_out, f_in)
    implicit none

    class(cuda_poisson_110_t) :: self
    class(field_t), intent(inout) :: f_out
    class(field_t), intent(in) :: f_in

    real(dp), device, pointer, dimension(:, :, :) :: f_out_dev, f_in_dev
    type(dim3) :: blocks, threads
    integer :: tpb

    select type (f_out)
    type is (cuda_field_t)
      f_out_dev => f_out%data_d
    end select
    select type (f_in)
    type is (cuda_field_t)
      f_in_dev => f_in%data_d
    end select

    f_out_dev = 0._dp

    tpb = min(self%ny_glob, 256)
    blocks = dim3(self%nz_loc, (self%ny_glob - 1)/tpb + 1, 1)
    threads = dim3(tpb, 1, 1)

    call undo_periodicity_xy<<<blocks, threads>>>( & !&
      f_out_dev, f_in_dev, self%nx_glob, self%ny_glob &
      )

  end subroutine undo_periodicity_xy_cuda

end module m_cuda_poisson_fft
