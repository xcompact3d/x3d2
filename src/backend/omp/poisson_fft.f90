module m_omp_poisson_fft

  use decomp_2d_constants, only: PHYSICAL_IN_X
  use decomp_2d, only: decomp_info, decomp_info_init, get_decomp_dims, &
                       transpose_x_to_y, transpose_y_to_x, &
                       transpose_z_to_y, transpose_y_to_z
  use decomp_2d_fft, only: decomp_2d_fft_init, decomp_2d_fft_3d, &
                           decomp_2d_fft_get_ph, decomp_2d_fft_get_sp

  use m_common, only: dp, CELL
  use m_field, only: field_t
  use m_mesh, only: mesh_t
  use m_poisson_fft, only: poisson_fft_t, poisson_bc_code, spectral_dirs
  use m_tdsops, only: dirps_t

  use m_omp_spectral, only: process_spectral_000, process_spectral_010, &
                            process_spectral_010_fw, &
                            process_spectral_010_poisson, &
                            process_spectral_010_bw, &
                            process_spectral_100, process_spectral_110, &
                            process_spectral_110_fw, &
                            process_spectral_110_solve, &
                            process_spectral_110_bw
  use m_logger, only: log_info, abort_run

  implicit none

  type, abstract, extends(poisson_fft_t) :: omp_poisson_fft_t
    !! FFT based Poisson solver, the state and the plain transform that
    !! every boundary condition case shares. Each case extends this type
    !! with its own buffers, see omp_poisson_000_t to omp_poisson_110_t.
    complex(dp), allocatable, dimension(:, :, :) :: c_x
    !> 2decomp processor grid. p_row splits y in the x-pencil, p_col splits z.
    integer :: p_row = 1, p_col = 1
    !> Spectral decomposition, and its y-pencil staging buffers. The FFT
    !! leaves the spectral slab in a z-pencil, which p_col splits along its
    !! dim2. In the 100 case that is the x modes, in the 010 case the y
    !! modes. The paired split in the postprocess couples a mode with its
    !! mirror in that same direction, so it needs that dimension whole.
    !! The y-pencil of the very same decomposition has dim2 complete for any
    !! processor grid, so 2decomp's own transpose delivers exactly the
    !! layout the pairing wants. In the 100 case and in the 010 case with a
    !! z split (p_col > 1), c_pair and waves_pair hold that y-pencil, and
    !! waves is redistributed once, at init.
    !! sp is associated for every case, since init also reads the spectral
    !! slab extents from it, but c_pair and waves_pair exist only when the
    !! pairing actually needs the hop. The 110 case also keeps a c_pair for
    !! its own hop, without a waves_pair.
    type(decomp_info), pointer :: sp => null()
    complex(dp), allocatable, dimension(:, :, :) :: c_pair, waves_pair
  contains
    procedure :: init => init_omp
    procedure :: init_buffers
    procedure :: init_pair_buffers
    procedure :: fft_forward => fft_forward_omp
    procedure :: fft_backward => fft_backward_omp
  end type omp_poisson_fft_t

  type, extends(omp_poisson_fft_t) :: omp_poisson_000_t
    !! Periodic in x, y and z
  contains
    procedure :: solve_poisson => solve_poisson_000_omp
    procedure :: fft_postprocess => fft_postprocess_000_omp
  end type omp_poisson_000_t

  type, extends(omp_poisson_fft_t) :: omp_poisson_010_t
    !! Periodic in x and z, non-periodic in y
    !> Working copies of the stretching matrices. The pentadiagonal
    !! elimination overwrites its matrices, so each solve starts from a copy
    !! of the pristine matrices of the base class.
    real(dp), allocatable, dimension(:, :, :, :) :: a_work_re, a_work_im
  contains
    procedure :: init_buffers => init_buffers_010
    procedure :: fft_postprocess => fft_postprocess_010_omp
    procedure :: enforce_periodicity => enforce_periodicity_y_omp
    procedure :: undo_periodicity => undo_periodicity_y_omp
  end type omp_poisson_010_t

  type, extends(omp_poisson_fft_t) :: omp_poisson_100_t
    !! Non-periodic in x, periodic in y and z
    !> Physical space decomposition owned by the FFT. The 100 case builds
    !! the transform on the swapped cell grid (ny, nx, nz), so ph%xsz is
    !! exactly the layout the transform reads, and ph%ysz is what a local
    !! dim1/dim2 swap of the untransposed x-pencil lands on. decomp_main
    !! cannot serve here because it is built on vertex dims.
    type(decomp_info), pointer :: ph => null()
    !> Y-pencil staging buffer of ph, only needed when y is decomposed.
    real(dp), allocatable, dimension(:, :, :) :: r_yp
    !> Exact sized, x-y transposed real work buffer.
    !! The FFT is initialised with the x and y dimensions swapped so that
    !! the r2c transform runs along the periodic y direction and the
    !! non-periodic x direction keeps its full mode range in dim2, which
    !! is what the paired split in the spectral postprocess needs. The
    !! buffer is also exactly (ny, nx, nz) rather than padded up to a
    !! multiple of SZ, because 2decomp plans its FFTW transforms over the
    !! flat array and a padded one would be read with the wrong stride.
    real(dp), allocatable, dimension(:, :, :) :: r_tr
  contains
    procedure :: init_buffers => init_buffers_100
    procedure :: fft_forward => fft_forward_100_omp
    procedure :: fft_backward => fft_backward_100_omp
    procedure :: fft_postprocess => fft_postprocess_100_omp
    procedure :: enforce_periodicity => enforce_periodicity_x_omp
    procedure :: undo_periodicity => undo_periodicity_x_omp
  end type omp_poisson_100_t

  type, extends(omp_poisson_fft_t) :: omp_poisson_110_t
    !! Non-periodic in x and y, periodic in z
    !> Physical space decomposition owned by the FFT, whose ph%xsz is the
    !! (nz, nx_loc, ny_loc) layout the transform reads.
    type(decomp_info), pointer :: ph => null()
    !> Transposed real work buffer (nz, nx_loc, ny_loc), the x-pencil of ph.
    real(dp), allocatable, dimension(:, :, :) :: r_tr
    !> Cell-dims decomposition, built separately because decomp_main is on
    !! vertex dims. Its x-pencil is the solver field layout (x whole, y
    !! split by p_row, z split by p_col); its z-pencil is, after a local
    !! (i,j,k)->(k,i,j) transpose, the x-pencil of ph.
    !! r_xp/r_yp/r_yp2/r_zp are its x-, y-, second y- and z-pencil staging
    !! buffers, only allocated on multiple ranks. r_yp2 is only needed when
    !! y is split (p_row > 1), since that is when the xy-periodicity fold
    !! needs a second y-pencil buffer alongside r_yp.
    type(decomp_info) :: dc
    real(dp), allocatable, dimension(:, :, :) :: r_xp, r_yp, r_yp2, r_zp
  contains
    procedure :: init_buffers => init_buffers_110
    procedure :: fft_forward => fft_forward_110_omp
    procedure :: fft_backward => fft_backward_110_omp
    procedure :: fft_postprocess => fft_postprocess_110_omp
    procedure :: enforce_periodicity => enforce_periodicity_xy_omp
    procedure :: undo_periodicity => undo_periodicity_xy_omp
  end type omp_poisson_110_t

  private :: init_omp, init_buffers, init_buffers_010, init_buffers_100, &
             init_buffers_110, fold_index, fold_destinations

contains

  subroutine create_omp_poisson_fft(poisson_fft, mesh, xdirps, ydirps, &
                                    zdirps, lowmem)
    !! Allocates the solver type of the mesh's boundary condition case and
    !! initialises it.
    implicit none

    class(poisson_fft_t), pointer, intent(out) :: poisson_fft
    type(mesh_t), target, intent(in) :: mesh
    class(dirps_t), intent(in) :: xdirps, ydirps, zdirps
    logical, optional, intent(in) :: lowmem

    select case (poisson_bc_code(mesh%grid%periodic_BC))
    case ('000')
      allocate (omp_poisson_000_t :: poisson_fft)
    case ('010')
      allocate (omp_poisson_010_t :: poisson_fft)
    case ('100')
      allocate (omp_poisson_100_t :: poisson_fft)
    case ('110')
      allocate (omp_poisson_110_t :: poisson_fft)
    case default
      call abort_run('Requested BCs are not supported in FFT-based Poisson &
                     &solver!')
    end select

    select type (poisson_fft)
    class is (omp_poisson_fft_t)
      call poisson_fft%init(mesh, xdirps, ydirps, zdirps, lowmem)
    end select

  end subroutine create_omp_poisson_fft

  subroutine init_omp(self, mesh, xdirps, ydirps, zdirps, lowmem)
    implicit none

    class(omp_poisson_fft_t) :: self
    type(mesh_t), target, intent(in) :: mesh
    class(dirps_t), intent(in) :: xdirps, ydirps, zdirps
    logical, optional, intent(in) :: lowmem
    integer :: dims(3), grid_dims(2), fft_dims(3)

    call log_info("Initialising 2decomp&fft")

    self%bc_code = poisson_bc_code(mesh%grid%periodic_BC)

    if (present(lowmem)) then
      call log_info('lowmem_fft has no impact in the OpenMP backend.')
    end if

    ! Get global cell dims
    dims = mesh%get_global_dims(CELL)

    ! Read the 2decomp processor grid (p_row, p_col), which the 100, 110 and
    ! 010 cases all use.
    grid_dims = get_decomp_dims()
    self%p_row = grid_dims(1)
    self%p_col = grid_dims(2)

    ! Build the transform on the cell grid permuted into the spectral
    ! layout of the case, so that the spectral slab comes back in the
    ! layout waves_set builds, e.g. (ny/2 + 1, nx, nz) for the 100 case.
    fft_dims = dims(spectral_dirs(self%bc_code))
    call decomp_2d_fft_init(PHYSICAL_IN_X, fft_dims(1), fft_dims(2), &
                            fft_dims(3))
    ! For PHYSICAL_IN_X the spectral slab is a z-pencil of sp, which is all
    ! decomp_2d_fft_get_size hands back. Read sp directly, subtracting one
    ! to convert the start position into an offset.
    self%sp => decomp_2d_fft_get_sp()
    if (mesh%geo%stretched(2) .and. self%bc_code == '010' &
        .and. self%p_col > 1) then
      ! The stretched solve runs on the y-pencil, where y is whole, so the
      ! stretching matrices must cover that block.
      call self%base_init(mesh, xdirps, ydirps, zdirps, &
                          self%sp%zsz, self%sp%zst - 1, &
                          n_blk=self%sp%ysz, blk_st=self%sp%yst - 1)
    else
      call self%base_init(mesh, xdirps, ydirps, zdirps, &
                          self%sp%zsz, self%sp%zst - 1)
    end if

    ! Only the top-bottom stretching on one rank is solved
    if (mesh%geo%stretched(2) .and. &
        (trim(mesh%geo%stretching(2)) /= 'top-bottom' &
         .or. self%p_row*self%p_col > 1)) then
      call abort_run('OpenMP backends FFT based Poisson solver supports only &
                     &top-bottom y stretching on one rank!')
    end if

    allocate (self%c_x(self%nx_spec, self%ny_spec, self%nz_spec))

    call self%init_buffers(mesh, dims)

  end subroutine init_omp

  subroutine init_buffers(self, mesh, dims)
    !! Allocates the work buffers of the case, called at the end of init
    !! with the global cell dims. The 000 case needs none beyond c_x.
    implicit none

    class(omp_poisson_fft_t) :: self
    type(mesh_t), intent(in) :: mesh
    integer, intent(in) :: dims(3)
  end subroutine init_buffers

  subroutine init_pair_buffers(self, with_waves)
    !! Allocates c_pair, the y-pencil of sp, and when with_waves is set
    !! also waves_pair, redistributing waves into it once.
    implicit none

    class(omp_poisson_fft_t) :: self
    logical, intent(in) :: with_waves

    allocate (self%c_pair(self%sp%ysz(1), self%sp%ysz(2), self%sp%ysz(3)))
    if (with_waves) then
      allocate (self%waves_pair(self%sp%ysz(1), self%sp%ysz(2), &
                                self%sp%ysz(3)))
      call transpose_z_to_y(self%waves, self%waves_pair, self%sp)
    end if

  end subroutine init_pair_buffers

  subroutine init_buffers_010(self, mesh, dims)
    implicit none

    class(omp_poisson_010_t) :: self
    type(mesh_t), intent(in) :: mesh
    integer, intent(in) :: dims(3)

    if (self%p_row > 1) then
      ! The periodicity reorder of the 010 case needs the y line whole.
      call abort_run('The OpenMP 010 Poisson solver does not support a y &
                     &decomposition, split along z instead.')
    end if

    ! Pairing a mode with its mirror in dim2 (the y modes) needs that
    ! dimension whole, which the y-pencil of sp provides.
    if (self%p_col > 1) call self%init_pair_buffers(with_waves=.true.)

    ! The top-bottom stretching matrices are symmetric, odd and even alike
    if (self%stretched_y) then
      allocate (self%a_work_re, mold=self%a_odd_re)
      allocate (self%a_work_im, mold=self%a_odd_im)
    end if

  end subroutine init_buffers_010

  subroutine init_buffers_100(self, mesh, dims)
    implicit none

    class(omp_poisson_100_t) :: self
    type(mesh_t), intent(in) :: mesh
    integer, intent(in) :: dims(3)

    ! Pairing a mode with its mirror in dim2 (the x modes) needs that
    ! dimension whole, which the y-pencil of sp provides.
    if (self%p_col > 1) call self%init_pair_buffers(with_waves=.true.)

    self%ph => decomp_2d_fft_get_ph()
    if (self%p_row > 1) then
      ! Reaching the transposed x-pencil takes a real redistribution once
      ! y is split: local dim1/dim2 swap -> y-pencil of ph -> x-pencil.
      allocate (self%r_yp(self%ph%ysz(1), self%ph%ysz(2), self%ph%ysz(3)))
      allocate (self%r_tr(self%ph%xsz(1), self%ph%xsz(2), self%ph%xsz(3)))
    else
      ! p_row == 1 makes the y-pencil identical to the x-pencil, so the
      ! swap alone reaches the transposed x-pencil with no communication.
      allocate (self%r_tr(self%ny_loc, self%nx_loc, self%nz_loc))
    end if

  end subroutine init_buffers_100

  subroutine init_buffers_110(self, mesh, dims)
    implicit none

    class(omp_poisson_110_t) :: self
    type(mesh_t), intent(in) :: mesh
    integer, intent(in) :: dims(3)

    ! ph%xsz is exactly (nz, nx_loc, ny_loc) for the 110 fft plan, on one
    ! rank as much as many.
    self%ph => decomp_2d_fft_get_ph()
    allocate (self%r_tr(self%ph%xsz(1), self%ph%xsz(2), self%ph%xsz(3)))

    if (mesh%par%nproc > 1) then
      ! Cell-dims decomposition of the solver field (x whole, y split by
      ! p_row, z split by p_col). decomp_main cannot serve here because
      ! it is built on vertex dims. Its z-pencil, after a local
      ! (i,j,k)->(k,i,j) transpose, is exactly the x-pencil of ph.
      call decomp_info_init(dims(1), dims(2), dims(3), self%dc)
      if (.not. all(self%dc%zsz == [self%ph%xsz(2), self%ph%xsz(3), &
                                    self%ph%xsz(1)])) then
        call abort_run('The 110 case cell-dims decomposition does not match &
                       &the FFT physical decomposition.')
      end if
      ! The mesh splits the vertices, and the last rank of a non-periodic
      ! direction has one cell fewer than vertices, while dc splits the
      ! cells directly. The two can hand a rank different y ranges (257
      ! cells on 2 ranks are 129 + 128 for the mesh but 128 + 129 for dc),
      ! and the copies between the solver field and r_xp assume they
      ! agree. A cell count divisible by p_row always agrees.
      if (.not. all(self%dc%xsz == [self%nx_loc, self%ny_loc, &
                                    self%nz_loc])) then
        call abort_run('The OpenMP 110 Poisson solver needs the y cell &
                       &count to be divisible by the number of ranks &
                       &along y, or split along z instead.')
      end if

      allocate (self%r_xp(self%dc%xsz(1), self%dc%xsz(2), self%dc%xsz(3)))
      allocate (self%r_yp(self%dc%ysz(1), self%dc%ysz(2), self%dc%ysz(3)))
      allocate (self%r_zp(self%dc%zsz(1), self%dc%zsz(2), self%dc%zsz(3)))

      if (self%p_row > 1) then
        allocate (self%r_yp2(self%dc%ysz(1), self%dc%ysz(2), &
                             self%dc%ysz(3)))
      end if

      if (self%p_col > 1) call self%init_pair_buffers(with_waves=.false.)
    end if

  end subroutine init_buffers_110

  subroutine solve_poisson_000_omp(self, f, temp)
    !! The 000 case has nothing to reorder, so it transforms f in place.
    implicit none

    class(omp_poisson_000_t) :: self
    class(field_t), intent(inout) :: f, temp

    call self%fft_forward(f)
    call self%fft_postprocess
    call self%fft_backward(f)

  end subroutine solve_poisson_000_omp

  subroutine fft_forward_omp(self, f_in)
    implicit none

    class(omp_poisson_fft_t) :: self
    class(field_t), intent(in) :: f_in

    call decomp_2d_fft_3d( &
      f_in%data(1:self%nx_loc, 1:self%ny_loc, 1:self%nz_loc), self%c_x &
      )

  end subroutine fft_forward_omp

  subroutine fft_forward_100_omp(self, f_in)
    !! Forward FFT for the non-periodic x case.
    !! Transposes x and y so that the non-periodic direction sits in dim2,
    !! where it keeps its full mode range, then transforms the (ny, nx, nz)
    !! buffer with the swapped plan set up in init.
    implicit none

    class(omp_poisson_100_t) :: self
    class(field_t), intent(in) :: f_in

    integer :: i, j, k, n1, n2, n3

    if (self%p_row > 1) then
      ! Strip the SZ padding and swap dim1 and dim2 in one pass, which
      ! lands on the y-pencil of ph, then redistribute to its x-pencil.
      n1 = self%ph%ysz(1); n2 = self%ph%ysz(2)
      n3 = self%ph%ysz(3)
      !$omp parallel do collapse(2)
      do k = 1, n3
        do i = 1, n2
          do j = 1, n1
            self%r_yp(j, i, k) = f_in%data(i, j, k)
          end do
        end do
      end do
      !$omp end parallel do

      call transpose_y_to_x(self%r_yp, self%r_tr, self%ph)
    else
      !$omp parallel do collapse(2)
      do k = 1, self%nz_loc
        do i = 1, self%nx_loc
          do j = 1, self%ny_loc
            self%r_tr(j, i, k) = f_in%data(i, j, k)
          end do
        end do
      end do
      !$omp end parallel do
    end if

    call decomp_2d_fft_3d(self%r_tr, self%c_x)

  end subroutine fft_forward_100_omp

  subroutine fft_forward_110_omp(self, f_in)
    !! Forward FFT for the non-periodic x and y case. Transposes the
    !! physical field (nx, ny, nz) to (nz, nx, ny) so that the r2c
    !! transform runs along the periodic z direction, then transforms
    !! with the plan set up in init.
    !!
    !! On multiple ranks the solver field's x-pencil is redistributed
    !! through the cell-dims decomposition dc to reach the FFT's x-pencil.
    !! When y is split (p_row > 1), enforce_periodicity_xy_omp could only
    !! fold x locally, so the y fold happens here instead, in the
    !! y-pencil of dc where the whole y range is local.
    implicit none

    class(omp_poisson_110_t) :: self
    class(field_t), intent(in) :: f_in

    integer :: i, j, k
    integer, allocatable :: dst_j(:)

    if (self%mesh%par%nproc == 1) then
      !$omp parallel do collapse(2)
      do k = 1, self%nz_loc
        do j = 1, self%ny_loc
          do i = 1, self%nx_loc
            self%r_tr(k, i, j) = f_in%data(i, j, k)
          end do
        end do
      end do
      !$omp end parallel do
    else
      !$omp parallel do collapse(2)
      do k = 1, self%nz_loc
        do j = 1, self%ny_loc
          do i = 1, self%nx_loc
            self%r_xp(i, j, k) = f_in%data(i, j, k)
          end do
        end do
      end do
      !$omp end parallel do

      call transpose_x_to_y(self%r_xp, self%r_yp, self%dc)

      if (self%p_row > 1) then
        ! Scatter, as in enforce_periodicity_y_omp
        dst_j = fold_destinations(self%ny_glob)

        !$omp parallel do collapse(2)
        do k = 1, self%dc%ysz(3)
          do j = 1, self%ny_glob
            do i = 1, self%dc%ysz(1)
              self%r_yp2(i, dst_j(j), k) = self%r_yp(i, j, k)
            end do
          end do
        end do
        !$omp end parallel do

        call transpose_y_to_z(self%r_yp2, self%r_zp, self%dc)
      else
        call transpose_y_to_z(self%r_yp, self%r_zp, self%dc)
      end if

      !$omp parallel do collapse(2)
      do k = 1, self%dc%zsz(3)
        do j = 1, self%dc%zsz(2)
          do i = 1, self%dc%zsz(1)
            self%r_tr(k, i, j) = self%r_zp(i, j, k)
          end do
        end do
      end do
      !$omp end parallel do
    end if

    call decomp_2d_fft_3d(self%r_tr, self%c_x)

  end subroutine fft_forward_110_omp

  subroutine fft_backward_omp(self, f_out)
    implicit none

    class(omp_poisson_fft_t) :: self
    class(field_t), intent(inout) :: f_out

    call decomp_2d_fft_3d( &
      self%c_x, f_out%data(1:self%nx_loc, 1:self%ny_loc, 1:self%nz_loc) &
      )

  end subroutine fft_backward_omp

  subroutine fft_backward_100_omp(self, f_out)
    !! Backward FFT for the non-periodic x case, undoing the x-y transpose
    !! applied by fft_forward_100_omp. Only the unpadded extents are
    !! written, which is all undo_periodicity_x_omp goes on to read.
    implicit none

    class(omp_poisson_100_t) :: self
    class(field_t), intent(inout) :: f_out

    integer :: i, j, k, n1, n2, n3

    call decomp_2d_fft_3d(self%c_x, self%r_tr)

    if (self%p_row > 1) then
      call transpose_x_to_y(self%r_tr, self%r_yp, self%ph)

      n1 = self%ph%ysz(1); n2 = self%ph%ysz(2)
      n3 = self%ph%ysz(3)
      !$omp parallel do collapse(2)
      do k = 1, n3
        do j = 1, n1
          do i = 1, n2
            f_out%data(i, j, k) = self%r_yp(j, i, k)
          end do
        end do
      end do
      !$omp end parallel do
    else
      !$omp parallel do collapse(2)
      do k = 1, self%nz_loc
        do j = 1, self%ny_loc
          do i = 1, self%nx_loc
            f_out%data(i, j, k) = self%r_tr(j, i, k)
          end do
        end do
      end do
      !$omp end parallel do
    end if

  end subroutine fft_backward_100_omp

  subroutine fft_backward_110_omp(self, f_out)
    !! Backward FFT for the non-periodic x and y case, undoing the
    !! transpose applied by fft_forward_110_omp. Only the unpadded
    !! extents are written, which is all undo_periodicity_xy_omp goes on
    !! to read.
    !!
    !! On multiple ranks this is the exact reverse of fft_forward_110_omp,
    !! undoing the y fold (when p_row > 1) in the y-pencil of dc before
    !! redistributing back to the solver field's x-pencil.
    implicit none

    class(omp_poisson_110_t) :: self
    class(field_t), intent(inout) :: f_out

    integer :: i, j, k
    integer, allocatable :: dst_j(:)

    call decomp_2d_fft_3d(self%c_x, self%r_tr)

    if (self%mesh%par%nproc == 1) then
      !$omp parallel do collapse(2)
      do k = 1, self%nz_loc
        do j = 1, self%ny_loc
          do i = 1, self%nx_loc
            f_out%data(i, j, k) = self%r_tr(k, i, j)
          end do
        end do
      end do
      !$omp end parallel do
    else
      !$omp parallel do collapse(2)
      do k = 1, self%dc%zsz(3)
        do j = 1, self%dc%zsz(2)
          do i = 1, self%dc%zsz(1)
            self%r_zp(i, j, k) = self%r_tr(k, i, j)
          end do
        end do
      end do
      !$omp end parallel do

      call transpose_z_to_y(self%r_zp, self%r_yp, self%dc)

      if (self%p_row > 1) then
        dst_j = fold_index([(j, j=1, self%ny_glob)], self%ny_glob)

        !$omp parallel do collapse(2)
        do k = 1, self%dc%ysz(3)
          do j = 1, self%ny_glob
            do i = 1, self%dc%ysz(1)
              self%r_yp2(i, dst_j(j), k) = self%r_yp(i, j, k)
            end do
          end do
        end do
        !$omp end parallel do

        call transpose_y_to_x(self%r_yp2, self%r_xp, self%dc)
      else
        call transpose_y_to_x(self%r_yp, self%r_xp, self%dc)
      end if

      f_out%data(1:self%nx_loc, 1:self%ny_loc, 1:self%nz_loc) = self%r_xp

    end if

  end subroutine fft_backward_110_omp

  subroutine fft_postprocess_000_omp(self)
    implicit none

    class(omp_poisson_000_t) :: self

    call process_spectral_000( &
      div_u=self%c_x, waves=self%waves, &
      nx_spec=self%nx_spec, ny_spec=self%ny_spec, nz_spec=self%nz_spec, &
      x_sp_st=self%sp_st(1), y_sp_st=self%sp_st(2), &
      z_sp_st=self%sp_st(3), &
      nx=self%nx_glob, ny=self%ny_glob, nz=self%nz_glob, &
      ax=self%ax, bx=self%bx, ay=self%ay, by=self%by, &
      az=self%az, bz=self%bz &
      )

  end subroutine fft_postprocess_000_omp

  subroutine fft_postprocess_010_omp(self)
    implicit none

    class(omp_poisson_010_t) :: self

    if (self%p_col > 1) then
      ! The FFT leaves y split across p_col in the z-pencil. Hop to the
      ! y-pencil, where y is whole, pair there, and hop back.
      call transpose_z_to_y(self%c_x, self%c_pair, self%sp)
      call process_spectral_010( &
        div_u=self%c_pair, waves=self%waves_pair, &
        nx_spec=self%sp%ysz(1), ny_spec=self%sp%ysz(2), &
        nz_spec=self%sp%ysz(3), &
        x_sp_st=self%sp%yst(1) - 1, y_sp_st=self%sp%yst(2) - 1, &
        z_sp_st=self%sp%yst(3) - 1, &
        nx=self%nx_glob, ny=self%ny_glob, nz=self%nz_glob, &
        ax=self%ax, bx=self%bx, ay=self%ay, by=self%by, &
        az=self%az, bz=self%bz &
        )
      call transpose_y_to_z(self%c_pair, self%c_x, self%sp)
    else if (self%stretched_y) then
      call process_spectral_010_fw( &
        div_u=self%c_x, &
        nx_spec=self%nx_spec, ny_spec=self%ny_spec, nz_spec=self%nz_spec, &
        x_sp_st=self%sp_st(1), y_sp_st=self%sp_st(2), &
        z_sp_st=self%sp_st(3), &
        nx=self%nx_glob, ny=self%ny_glob, nz=self%nz_glob, &
        ax=self%ax, bx=self%bx, ay=self%ay, by=self%by, &
        az=self%az, bz=self%bz &
        )

      ! The odd rows, then the even rows, each with its own matrices
      self%a_work_re = self%a_odd_re
      self%a_work_im = self%a_odd_im
      call process_spectral_010_poisson( &
        div_u=self%c_x, a_re=self%a_work_re, a_im=self%a_work_im, &
        off=0, inc=2, nx_spec=self%nx_spec, n=self%ny_spec/2, &
        nz_spec=self%nz_spec, x_sp_st=self%sp_st(1), &
        z_sp_st=self%sp_st(3), &
        nx=self%nx_glob, ny=self%ny_glob, nz=self%nz_glob &
        )

      self%a_work_re = self%a_even_re
      self%a_work_im = self%a_even_im
      call process_spectral_010_poisson( &
        div_u=self%c_x, a_re=self%a_work_re, a_im=self%a_work_im, &
        off=1, inc=2, nx_spec=self%nx_spec, n=self%ny_spec/2, &
        nz_spec=self%nz_spec, x_sp_st=self%sp_st(1), &
        z_sp_st=self%sp_st(3), &
        nx=self%nx_glob, ny=self%ny_glob, nz=self%nz_glob &
        )

      call process_spectral_010_bw( &
        div_u=self%c_x, &
        nx_spec=self%nx_spec, ny_spec=self%ny_spec, nz_spec=self%nz_spec, &
        x_sp_st=self%sp_st(1), y_sp_st=self%sp_st(2), &
        z_sp_st=self%sp_st(3), &
        nx=self%nx_glob, ny=self%ny_glob, nz=self%nz_glob, &
        ax=self%ax, bx=self%bx, ay=self%ay, by=self%by, &
        az=self%az, bz=self%bz &
        )
    else
      call process_spectral_010( &
        div_u=self%c_x, waves=self%waves, &
        nx_spec=self%nx_spec, ny_spec=self%ny_spec, nz_spec=self%nz_spec, &
        x_sp_st=self%sp_st(1), y_sp_st=self%sp_st(2), &
        z_sp_st=self%sp_st(3), &
        nx=self%nx_glob, ny=self%ny_glob, nz=self%nz_glob, &
        ax=self%ax, bx=self%bx, ay=self%ay, by=self%by, &
        az=self%az, bz=self%bz &
        )
    end if

  end subroutine fft_postprocess_010_omp

  subroutine fft_postprocess_100_omp(self)
    !! Post-process div U* in spectral space for the non-periodic x case.
    !!
    !! After the x-y transpose applied by fft_forward_100_omp, dim1 holds
    !! the periodic y r2c modes and dim2 holds the non-periodic x modes.
    !! process_spectral_100 handles that arrangement directly.
    !!
    !! The CUDA backend stages the equivalent work across three kernels so
    !! that a mode can fetch its pairing partner from another GPU. On a
    !! single rank every partner is local, so the staging collapses back
    !! into this one fused pass.
    implicit none

    class(omp_poisson_100_t) :: self

    if (self%p_col > 1) then
      ! Hop to the y-pencil, where dim2 is whole, pair there, and hop back.
      call transpose_z_to_y(self%c_x, self%c_pair, self%sp)
      call process_spectral_100( &
        div_u=self%c_pair, waves=self%waves_pair, &
        nx_spec=self%sp%ysz(1), ny_spec=self%sp%ysz(2), &
        nz_spec=self%sp%ysz(3), &
        x_sp_st=self%sp%yst(1) - 1, y_sp_st=self%sp%yst(2) - 1, &
        z_sp_st=self%sp%yst(3) - 1, &
        nx=self%nx_glob, ny=self%ny_glob, nz=self%nz_glob, &
        ax=self%ax, bx=self%bx, ay=self%ay, by=self%by, &
        az=self%az, bz=self%bz &
        )
      call transpose_y_to_z(self%c_pair, self%c_x, self%sp)
    else
      call process_spectral_100( &
        div_u=self%c_x, waves=self%waves, &
        nx_spec=self%nx_spec, ny_spec=self%ny_spec, nz_spec=self%nz_spec, &
        x_sp_st=self%sp_st(1), y_sp_st=self%sp_st(2), &
        z_sp_st=self%sp_st(3), &
        nx=self%nx_glob, ny=self%ny_glob, nz=self%nz_glob, &
        ax=self%ax, bx=self%bx, ay=self%ay, by=self%by, &
        az=self%az, bz=self%bz &
        )
    end if

  end subroutine fft_postprocess_100_omp

  subroutine fft_postprocess_110_omp(self)
    !! Post-process div U* in spectral space for the non-periodic x and y
    !! case. After the (nx, ny, nz) -> (nz, nx, ny) transpose applied by
    !! fft_forward_110_omp, dim1 holds the periodic z r2c modes, dim2
    !! holds the non-periodic x modes and dim3 holds the non-periodic y
    !! modes.
    !!
    !! The z-pencil of sp (which c_x is) splits dim2 (x modes) across
    !! p_col, so the X paired split and recombine, which each couple a
    !! mode with its mirror in dim2, hop to the y-pencil of sp (where
    !! dim2 is whole) and back, the same way fft_postprocess_100_omp does.
    !! The Y paired split/solve/recombine run directly on the z-pencil,
    !! where dim3 (Y modes) is always whole.
    implicit none

    class(omp_poisson_110_t) :: self

    if (self%p_col == 1) then
      call process_spectral_110( &
        div_u=self%c_x, waves=self%waves, &
        n1=self%nx_spec, n2=self%ny_spec, n3=self%nz_spec, &
        st1=self%sp_st(1), st2=self%sp_st(2), st3=self%sp_st(3), &
        nx=self%nx_glob, ny=self%ny_glob, nz=self%nz_glob, &
        ax=self%ax, bx=self%bx, ay=self%ay, by=self%by, &
        az=self%az, bz=self%bz &
        )
    else
      call transpose_z_to_y(self%c_x, self%c_pair, self%sp)
      call process_spectral_110_fw( &
        div_u=self%c_pair, n1=self%sp%ysz(1), n2=self%sp%ysz(2), &
        n3=self%sp%ysz(3), st1=self%sp%yst(1) - 1, st2=self%sp%yst(2) - 1, &
        nx=self%nx_glob, ny=self%ny_glob, nz=self%nz_glob, &
        ax=self%ax, bx=self%bx, az=self%az, bz=self%bz &
        )
      call transpose_y_to_z(self%c_pair, self%c_x, self%sp)

      call process_spectral_110_solve( &
        div_u=self%c_x, waves=self%waves, &
        n1=self%nx_spec, n2=self%ny_spec, n3=self%nz_spec, &
        st3=self%sp_st(3), ay=self%ay, by=self%by &
        )

      call transpose_z_to_y(self%c_x, self%c_pair, self%sp)
      call process_spectral_110_bw( &
        div_u=self%c_pair, n1=self%sp%ysz(1), n2=self%sp%ysz(2), &
        n3=self%sp%ysz(3), st1=self%sp%yst(1) - 1, st2=self%sp%yst(2) - 1, &
        nx=self%nx_glob, ny=self%ny_glob, nz=self%nz_glob, &
        ax=self%ax, bx=self%bx, az=self%az, bz=self%bz &
        )
      call transpose_y_to_z(self%c_pair, self%c_x, self%sp)
    end if

  end subroutine fft_postprocess_110_omp

  subroutine enforce_periodicity_x_omp(self, f_out, f_in)
    !! Gathers the non-periodic x line into a periodic one so that a plain
    !! FFT can stand in for the cosine transform.
    implicit none

    class(omp_poisson_100_t) :: self
    class(field_t), intent(inout) :: f_out
    class(field_t), intent(in) :: f_in

    integer :: i, j, k
    integer, allocatable :: src_i(:)

    src_i = fold_index([(i, i=1, self%nx_glob)], self%nx_glob)

    !$omp parallel do collapse(2)
    do k = 1, self%nz_loc
      do j = 1, self%ny_loc
        do i = 1, self%nx_glob
          f_out%data(i, j, k) = f_in%data(src_i(i), j, k)
        end do
      end do
    end do
    !$omp end parallel do

  end subroutine enforce_periodicity_x_omp

  subroutine undo_periodicity_x_omp(self, f_out, f_in)
    !! Scatters the gathered x line back to its original ordering.
    implicit none

    class(omp_poisson_100_t) :: self
    class(field_t), intent(inout) :: f_out
    class(field_t), intent(in) :: f_in

    integer :: i, j, k
    integer, allocatable :: dst_i(:)

    ! Scatter through the gather map rather than gather through its
    ! inverse: the reads stay contiguous, which measured faster.
    dst_i = fold_index([(i, i=1, self%nx_glob)], self%nx_glob)

    ! The gathered ordering only covers 1..nx_glob, but the block is padded
    ! up to a multiple of SZ and reorder copies the whole padded extent back
    ! out to the pressure field, so the untouched entries have to be cleared
    ! rather than left holding whatever the recycled block came with.
    !$omp parallel workshare
    f_out%data = 0._dp
    !$omp end parallel workshare

    !$omp parallel do collapse(2)
    do k = 1, self%nz_loc
      do j = 1, self%ny_loc
        do i = 1, self%nx_glob
          f_out%data(dst_i(i), j, k) = f_in%data(i, j, k)
        end do
      end do
    end do
    !$omp end parallel do

  end subroutine undo_periodicity_x_omp

  subroutine enforce_periodicity_y_omp(self, f_out, f_in)
    implicit none

    class(omp_poisson_010_t) :: self
    class(field_t), intent(inout) :: f_out
    class(field_t), intent(in) :: f_in

    integer :: i, j, k
    integer, allocatable :: dst_j(:)

    ! Scatter through the inverse map rather than gather through the
    ! fold: the reads stay contiguous, which measured faster.
    dst_j = fold_destinations(self%ny_glob)

    !$omp parallel do collapse(2)
    do k = 1, self%nz_loc
      do j = 1, self%ny_glob
        do i = 1, self%nx_loc
          f_out%data(i, dst_j(j), k) = f_in%data(i, j, k)
        end do
      end do
    end do
    !$omp end parallel do

  end subroutine enforce_periodicity_y_omp

  subroutine undo_periodicity_y_omp(self, f_out, f_in)
    implicit none

    class(omp_poisson_010_t) :: self
    class(field_t), intent(inout) :: f_out
    class(field_t), intent(in) :: f_in

    integer :: i, j, k
    integer, allocatable :: dst_j(:)

    ! Scatter through the gather map, as in undo_periodicity_x_omp.
    dst_j = fold_index([(j, j=1, self%ny_glob)], self%ny_glob)

    !$omp parallel do collapse(2)
    do k = 1, self%nz_loc
      do j = 1, self%ny_glob
        do i = 1, self%nx_loc
          f_out%data(i, dst_j(j), k) = f_in%data(i, j, k)
        end do
      end do
    end do
    !$omp end parallel do

  end subroutine undo_periodicity_y_omp

  subroutine enforce_periodicity_xy_omp(self, f_out, f_in)
    !! Gathers the non-periodic x and, when it is not split across ranks,
    !! y lines into periodic ones so that a plain FFT can stand in for
    !! the cosine transform in both directions.
    !!
    !! When p_row > 1 the y direction is split across ranks, so only the
    !! x fold can be done locally here; fft_forward_110_omp folds y once
    !! the field has been redistributed to a y-pencil where the whole y
    !! range is local.
    implicit none

    class(omp_poisson_110_t) :: self
    class(field_t), intent(inout) :: f_out
    class(field_t), intent(in) :: f_in

    integer :: i, j, k
    integer, allocatable :: src_i(:), dst_i(:), dst_j(:)

    ! The maps stay in global terms: they are a property of the global
    ! line, and nx_loc == nx_glob always for the 110 case's x-pencil
    ! field, and ny_loc == ny_glob when p_row == 1 since y is not split.
    if (self%p_row == 1) then
      ! Scatter through the inverse maps, as in enforce_periodicity_y_omp
      dst_i = fold_destinations(self%nx_glob)
      dst_j = fold_destinations(self%ny_glob)

      !$omp parallel do collapse(2)
      do k = 1, self%nz_loc
        do j = 1, self%ny_loc
          do i = 1, self%nx_loc
            f_out%data(dst_i(i), dst_j(j), k) = f_in%data(i, j, k)
          end do
        end do
      end do
      !$omp end parallel do
    else
      ! A fold in x alone is faster as a gather, as in
      ! enforce_periodicity_x_omp
      src_i = fold_index([(i, i=1, self%nx_glob)], self%nx_glob)

      !$omp parallel do collapse(2)
      do k = 1, self%nz_loc
        do j = 1, self%ny_loc
          do i = 1, self%nx_loc
            f_out%data(i, j, k) = f_in%data(src_i(i), j, k)
          end do
        end do
      end do
      !$omp end parallel do
    end if

  end subroutine enforce_periodicity_xy_omp

  subroutine undo_periodicity_xy_omp(self, f_out, f_in)
    !! Scatters the gathered x and, when it is not split across ranks, y
    !! lines back to their original ordering.
    !!
    !! Mirrors enforce_periodicity_xy_omp: the combined undo runs locally
    !! when p_row == 1, otherwise only x is undone here, since
    !! fft_backward_110_omp already undoes y while the whole y range is
    !! still local, in the y-pencil.
    implicit none

    class(omp_poisson_110_t) :: self
    class(field_t), intent(inout) :: f_out
    class(field_t), intent(in) :: f_in

    integer :: i, j, k
    integer, allocatable :: dst_i(:), dst_j(:)

    ! Scatter through the gather maps, as in undo_periodicity_x_omp. They
    ! stay in global terms, as in enforce_periodicity_xy_omp.
    dst_i = fold_index([(i, i=1, self%nx_glob)], self%nx_glob)

    ! The gathered ordering only covers the unpadded extent, but the
    ! block is padded up to a multiple of SZ and reorder copies the
    ! whole padded extent back out to the pressure field, so the
    ! untouched entries have to be cleared rather than left holding
    ! whatever the recycled block came with.
    !$omp parallel workshare
    f_out%data = 0._dp
    !$omp end parallel workshare

    if (self%p_row == 1) then
      dst_j = fold_index([(j, j=1, self%ny_glob)], self%ny_glob)

      !$omp parallel do collapse(2)
      do k = 1, self%nz_loc
        do j = 1, self%ny_loc
          do i = 1, self%nx_loc
            f_out%data(dst_i(i), dst_j(j), k) = f_in%data(i, j, k)
          end do
        end do
      end do
      !$omp end parallel do
    else
      !$omp parallel do collapse(2)
      do k = 1, self%nz_loc
        do j = 1, self%ny_loc
          do i = 1, self%nx_loc
            f_out%data(dst_i(i), j, k) = f_in%data(i, j, k)
          end do
        end do
      end do
      !$omp end parallel do
    end if

  end subroutine undo_periodicity_xy_omp

  pure elemental integer function fold_index(i, n) result(src)
    !! Source index of entry i of a non-periodic line of n points gathered
    !! into a periodic one: the odd entries ascending, then the even
    !! entries descending, so that a plain FFT can stand in for the cosine
    !! transform. For an odd n the last entry lands at the centre, n/2 + 1.
    !! The single definition of the fold: the undo routines scatter
    !! through it, and the enforce routines gather through it (x) or
    !! scatter through its inverse, fold_destinations (y, xy). Each is the
    !! form that measured fastest.
    implicit none

    integer, intent(in) :: i, n

    if (i <= n/2) then
      src = 2*i - 1
    else if (mod(n, 2) == 1 .and. i == n/2 + 1) then
      src = n
    else
      src = 2*n - 2*i + 2
    end if

  end function fold_index

  pure function fold_destinations(n) result(dst)
    !! Where each entry of a non-periodic line of n points lands when it is
    !! gathered by fold_index: the inverse permutation, built from
    !! fold_index itself so the fold keeps a single definition.
    implicit none

    integer, intent(in) :: n
    integer, allocatable :: dst(:)

    integer :: i

    allocate (dst(n))
    dst(fold_index([(i, i=1, n)], n)) = [(i, i=1, n)]

  end function fold_destinations
end module m_omp_poisson_fft
