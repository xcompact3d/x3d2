module m_stats
  !! Online accumulation of flow field statistics.
  !!
  !! Running means are maintained using the incremental update:
  !! \[ \bar{x}_n = \bar{x}_{n-1} + \frac{x_n - \bar{x}_{n-1}}{n} \]
  !! applied to \(u, v, w, u^2, v^2, w^2, uv, uw, vw\), plus pressure
  !! mean and scalar statistics (\(\phi, \phi^2\)) when present.
  !! At write time the fluctuation quantities are derived:
  !! \[ u' = \sqrt{\max(0,\, \overline{u^2} - \bar{u}^2)}, \quad
  !!    \langle u'v' \rangle = \overline{uv} - \bar{u}\bar{v} \]
  !!
  !! Statistics are written at `istatout` frequency. Accumulation starts
  !! at `initstat` and samples every `istatfreq` steps.
  !!
  !! Note: velocity fields must be at `VERT` data location when `update`
  !! is called (i.e. after `pressure_correction`).
  use mpi, only: MPI_COMM_WORLD, MPI_Comm_rank
  use m_common, only: dp, i8, DIR_C, DIR_X, VERT, get_argument
  use m_config, only: stats_config_t
  use m_field, only: field_t
  use m_solver, only: solver_t
  use m_io_session, only: writer_session_t, reader_session_t

  implicit none

  private
  public :: stats_manager_t, accumulate_mean

  type :: stats_manager_t
    type(stats_config_t) :: config
    integer :: sample_count = 0
    logical :: is_active = .false.
    !! Running means of first moments
    real(dp), allocatable :: umean(:, :, :)
    real(dp), allocatable :: vmean(:, :, :)
    real(dp), allocatable :: wmean(:, :, :)
    !! Running means of second moments
    real(dp), allocatable :: uumean(:, :, :)  !! \(\overline{u^2}\)
    real(dp), allocatable :: vvmean(:, :, :)  !! \(\overline{v^2}\)
    real(dp), allocatable :: wwmean(:, :, :)  !! \(\overline{w^2}\)
    real(dp), allocatable :: uvmean(:, :, :)  !! \(\overline{uv}\)
    real(dp), allocatable :: uwmean(:, :, :)  !! \(\overline{uw}\)
    real(dp), allocatable :: vwmean(:, :, :)  !! \(\overline{vw}\)
    !! Pressure mean
    real(dp), allocatable :: pmean(:, :, :)
    !! Scalar (species) statistics
    integer :: nspecies = 0
    real(dp), allocatable :: phimean(:, :, :, :)     !! \(\overline{\phi}\)
    real(dp), allocatable :: phiphimean(:, :, :, :)  !! \(\overline{\phi^2}\)
  contains
    procedure :: init
    procedure :: update
    procedure :: write_stats
    procedure :: write_checkpoint
    procedure :: read_checkpoint
    procedure :: finalise
  end type stats_manager_t

contains

  pure subroutine accumulate_mean(mean, val, stat_inc)
    !! Incremental running-mean update:
    !! \( \bar{x}_n = \bar{x}_{n-1} + (x_n - \bar{x}_{n-1})/n \)
    !! where `stat_inc = 1/n`.
    real(dp), intent(inout) :: mean(:, :, :)
    real(dp), intent(in) :: val(:, :, :)
    real(dp), intent(in) :: stat_inc

    mean = mean + (val - mean)*stat_inc
  end subroutine accumulate_mean

  subroutine init(self, solver, comm)
    !! Initialise the statistics manager: read config and allocate accumulators.
    class(stats_manager_t), intent(inout) :: self
    class(solver_t), intent(in) :: solver
    integer, intent(in) :: comm

    integer :: dims(3), myrank, ierr

    self%config = stats_config_t()
    call self%config%read(nml_file=get_argument(1))

    if (self%config%initstat <= 0) return

    self%is_active = .true.
    dims = solver%mesh%get_dims(VERT)

    allocate (self%umean(dims(1), dims(2), dims(3)), source=0.0_dp)
    allocate (self%vmean(dims(1), dims(2), dims(3)), source=0.0_dp)
    allocate (self%wmean(dims(1), dims(2), dims(3)), source=0.0_dp)
    allocate (self%uumean(dims(1), dims(2), dims(3)), source=0.0_dp)
    allocate (self%vvmean(dims(1), dims(2), dims(3)), source=0.0_dp)
    allocate (self%wwmean(dims(1), dims(2), dims(3)), source=0.0_dp)
    allocate (self%uvmean(dims(1), dims(2), dims(3)), source=0.0_dp)
    allocate (self%uwmean(dims(1), dims(2), dims(3)), source=0.0_dp)
    allocate (self%vwmean(dims(1), dims(2), dims(3)), source=0.0_dp)

    if (solver%keep_pressure) then
      allocate (self%pmean(dims(1), dims(2), dims(3)), source=0.0_dp)
    end if

    self%nspecies = solver%nspecies
    if (self%nspecies > 0) then
      allocate (self%phimean(dims(1), dims(2), dims(3), self%nspecies), &
                source=0.0_dp)
      allocate (self%phiphimean(dims(1), dims(2), dims(3), self%nspecies), &
                source=0.0_dp)
    end if

    call MPI_Comm_rank(comm, myrank, ierr)
    if (myrank == 0) then
      print *, 'Statistics: initstat=', self%config%initstat, &
        ' istatfreq=', self%config%istatfreq, &
        ' istatout=', self%config%istatout
    end if
  end subroutine init

  subroutine update(self, solver, iter)
    !! Accumulate running means for the current iteration.
    !! Velocity fields must be at `VERT` data location.
    class(stats_manager_t), intent(inout) :: self
    class(solver_t), intent(in) :: solver
    integer, intent(in) :: iter

    class(field_t), pointer :: u_tmp, v_tmp, w_tmp, p_tmp, phi_tmp
    integer :: dims(3), is
    real(dp) :: stat_inc

    if (.not. self%is_active) return
    if (iter < self%config%initstat) return
    if (mod(iter - self%config%initstat, self%config%istatfreq) /= 0) return

    dims = solver%mesh%get_dims(VERT)

    u_tmp => solver%host_allocator%get_block(DIR_C, VERT)
    v_tmp => solver%host_allocator%get_block(DIR_C, VERT)
    w_tmp => solver%host_allocator%get_block(DIR_C, VERT)

    call solver%backend%get_field_data(u_tmp%data, solver%u)
    call solver%backend%get_field_data(v_tmp%data, solver%v)
    call solver%backend%get_field_data(w_tmp%data, solver%w)

    self%sample_count = self%sample_count + 1
    stat_inc = 1.0_dp/real(self%sample_count, dp)

    associate ( &
      u => u_tmp%data(1:dims(1), 1:dims(2), 1:dims(3)), &
      v => v_tmp%data(1:dims(1), 1:dims(2), 1:dims(3)), &
      w => w_tmp%data(1:dims(1), 1:dims(2), 1:dims(3)) &
      )
      call accumulate_mean(self%umean, u, stat_inc)
      call accumulate_mean(self%vmean, v, stat_inc)
      call accumulate_mean(self%wmean, w, stat_inc)
      call accumulate_mean(self%uumean, u*u, stat_inc)
      call accumulate_mean(self%vvmean, v*v, stat_inc)
      call accumulate_mean(self%wwmean, w*w, stat_inc)
      call accumulate_mean(self%uvmean, u*v, stat_inc)
      call accumulate_mean(self%uwmean, u*w, stat_inc)
      call accumulate_mean(self%vwmean, v*w, stat_inc)
    end associate

    ! Pressure mean (from pressure_vert, already on VERT grid)
    if (allocated(self%pmean) .and. associated(solver%pressure_vert)) then
      p_tmp => solver%host_allocator%get_block(DIR_C, VERT)
      call solver%backend%get_field_data(p_tmp%data, solver%pressure_vert)
      associate (p => p_tmp%data(1:dims(1), 1:dims(2), 1:dims(3)))
        call accumulate_mean(self%pmean, p, stat_inc)
      end associate
      call solver%host_allocator%release_block(p_tmp)
    end if

    ! Scalar (species) statistics
    do is = 1, self%nspecies
      phi_tmp => solver%host_allocator%get_block(DIR_C, VERT)
      call solver%backend%get_field_data(phi_tmp%data, &
                                         solver%species(is)%ptr)
      associate (phi => phi_tmp%data(1:dims(1), 1:dims(2), 1:dims(3)))
        call accumulate_mean(self%phimean(:, :, :, is), phi, stat_inc)
        call accumulate_mean(self%phiphimean(:, :, :, is), phi*phi, stat_inc)
      end associate
      call solver%host_allocator%release_block(phi_tmp)
    end do

    call solver%host_allocator%release_block(u_tmp)
    call solver%host_allocator%release_block(v_tmp)
    call solver%host_allocator%release_block(w_tmp)
  end subroutine update

  subroutine write_stats(self, solver, timestep, comm)
    !! Write statistics to file. Output fields are mean velocities
    !! (`umean`, `vmean`, `wmean`), RMS fluctuations (`uprime`, `vprime`,
    !! `wprime`), and Reynolds stresses (`uvmean`, `uwmean`, `vwmean`).
    class(stats_manager_t), intent(inout) :: self
    class(solver_t), intent(in) :: solver
    integer, intent(in) :: timestep
    integer, intent(in) :: comm

    type(writer_session_t) :: writer_session
    character(len=256) :: filename
    integer(i8), dimension(3) :: shape_dims, start_dims, count_dims
    integer :: dims(3), myrank, ierr, is
    real(dp), allocatable :: uprime(:, :, :)
    real(dp), allocatable :: vprime(:, :, :)
    real(dp), allocatable :: wprime(:, :, :)
    real(dp), allocatable :: uv(:, :, :)
    real(dp), allocatable :: uw(:, :, :)
    real(dp), allocatable :: vw(:, :, :)
    real(dp), allocatable :: phiprime(:, :, :)
    character(len=64) :: field_name

    if (.not. self%is_active) return
    if (self%config%istatout <= 0) return
    if (mod(timestep, self%config%istatout) /= 0) return

    call MPI_Comm_rank(comm, myrank, ierr)

    write (filename, '(A,A,I0.6,A)') &
      trim(self%config%stats_prefix), '_', timestep, '.bp'

    dims = solver%mesh%get_dims(VERT)
    shape_dims = int(solver%mesh%get_global_dims(VERT), i8)
    start_dims = int(solver%mesh%par%n_offset, i8)
    count_dims = int(dims, i8)

    allocate (uprime(dims(1), dims(2), dims(3)))
    allocate (vprime(dims(1), dims(2), dims(3)))
    allocate (wprime(dims(1), dims(2), dims(3)))
    allocate (uv(dims(1), dims(2), dims(3)))
    allocate (uw(dims(1), dims(2), dims(3)))
    allocate (vw(dims(1), dims(2), dims(3)))

    uprime = sqrt(max(0.0_dp, self%uumean - self%umean**2))
    vprime = sqrt(max(0.0_dp, self%vvmean - self%vmean**2))
    wprime = sqrt(max(0.0_dp, self%wwmean - self%wmean**2))
    uv = self%uvmean - self%umean*self%vmean
    uw = self%uwmean - self%umean*self%wmean
    vw = self%vwmean - self%vmean*self%wmean

    call writer_session%open(filename, comm)
    if (writer_session%is_session_functional() .and. myrank == 0) then
      print *, 'Writing statistics: ', trim(filename), &
        ' (samples=', self%sample_count, ')'
    end if

    call writer_session%write_data("sample_count", self%sample_count)
    call writer_session%write_data("umean", self%umean, &
                                   shape_dims, start_dims, count_dims)
    call writer_session%write_data("vmean", self%vmean, &
                                   shape_dims, start_dims, count_dims)
    call writer_session%write_data("wmean", self%wmean, &
                                   shape_dims, start_dims, count_dims)
    call writer_session%write_data("uprime", uprime, &
                                   shape_dims, start_dims, count_dims)
    call writer_session%write_data("vprime", vprime, &
                                   shape_dims, start_dims, count_dims)
    call writer_session%write_data("wprime", wprime, &
                                   shape_dims, start_dims, count_dims)
    call writer_session%write_data("uvmean", uv, &
                                   shape_dims, start_dims, count_dims)
    call writer_session%write_data("uwmean", uw, &
                                   shape_dims, start_dims, count_dims)
    call writer_session%write_data("vwmean", vw, &
                                   shape_dims, start_dims, count_dims)

    ! Pressure mean (only when output_pressure is enabled)
    if (allocated(self%pmean)) then
      call writer_session%write_data("pmean", self%pmean, &
                                     shape_dims, start_dims, count_dims)
    end if

    ! Scalar statistics
    if (self%nspecies > 0) then
      allocate (phiprime(dims(1), dims(2), dims(3)))
      do is = 1, self%nspecies
        write (field_name, '(A,I0)') 'phimean_', is
        call writer_session%write_data( &
          trim(field_name), self%phimean(:, :, :, is), &
          shape_dims, start_dims, count_dims &
          )

        phiprime = sqrt(max(0.0_dp, self%phiphimean(:, :, :, is) &
                            - self%phimean(:, :, :, is)**2))
        write (field_name, '(A,I0)') 'phiprime_', is
        call writer_session%write_data( &
          trim(field_name), phiprime, shape_dims, start_dims, count_dims &
          )
      end do
      deallocate (phiprime)
    end if

    call writer_session%close()

    deallocate (uprime, vprime, wprime, uv, uw, vw)
  end subroutine write_stats

  subroutine write_checkpoint(self, solver, writer_session)
    !! Write running means into an already-open checkpoint session.
    !! Saves all accumulator arrays and sample_count so that
    !! statistics can be resumed exactly after a restart.
    class(stats_manager_t), intent(inout) :: self
    class(solver_t), intent(in) :: solver
    type(writer_session_t), intent(inout) :: writer_session

    integer(i8), dimension(3) :: shape_dims, start_dims, count_dims
    integer :: dims(3), is
    character(len=64) :: field_name

    if (.not. self%is_active) return

    dims = solver%mesh%get_dims(VERT)
    shape_dims = int(solver%mesh%get_global_dims(VERT), i8)
    start_dims = int(solver%mesh%par%n_offset, i8)
    count_dims = int(dims, i8)

    call writer_session%write_data("stats_sample_count", self%sample_count)

    ! First moments
    call writer_session%write_data("stats_umean", self%umean, &
                                   shape_dims, start_dims, count_dims)
    call writer_session%write_data("stats_vmean", self%vmean, &
                                   shape_dims, start_dims, count_dims)
    call writer_session%write_data("stats_wmean", self%wmean, &
                                   shape_dims, start_dims, count_dims)

    ! Second moments
    call writer_session%write_data("stats_uumean", self%uumean, &
                                   shape_dims, start_dims, count_dims)
    call writer_session%write_data("stats_vvmean", self%vvmean, &
                                   shape_dims, start_dims, count_dims)
    call writer_session%write_data("stats_wwmean", self%wwmean, &
                                   shape_dims, start_dims, count_dims)
    call writer_session%write_data("stats_uvmean", self%uvmean, &
                                   shape_dims, start_dims, count_dims)
    call writer_session%write_data("stats_uwmean", self%uwmean, &
                                   shape_dims, start_dims, count_dims)
    call writer_session%write_data("stats_vwmean", self%vwmean, &
                                   shape_dims, start_dims, count_dims)

    ! Pressure mean
    if (allocated(self%pmean)) then
      call writer_session%write_data("stats_pmean", self%pmean, &
                                     shape_dims, start_dims, count_dims)
    end if

    ! Scalar statistics
    do is = 1, self%nspecies
      write (field_name, '(A,I0)') 'stats_phimean_', is
      call writer_session%write_data( &
        trim(field_name), self%phimean(:, :, :, is), &
        shape_dims, start_dims, count_dims &
        )
      write (field_name, '(A,I0)') 'stats_phiphimean_', is
      call writer_session%write_data( &
        trim(field_name), self%phiphimean(:, :, :, is), &
        shape_dims, start_dims, count_dims &
        )
    end do
  end subroutine write_checkpoint

  subroutine read_checkpoint(self, solver, reader_session)
    !! Restore running means from an already-open checkpoint session.
    !! If stats were not present in the checkpoint, statistics start fresh.
    class(stats_manager_t), intent(inout) :: self
    class(solver_t), intent(in) :: solver
    type(reader_session_t), intent(inout) :: reader_session

    integer(i8), dimension(3) :: start_dims, count_dims
    integer :: dims(3), myrank, ierr, is
    character(len=64) :: field_name

    if (.not. self%is_active) return

    call MPI_Comm_rank(MPI_COMM_WORLD, myrank, ierr)

    dims = solver%mesh%get_dims(VERT)
    start_dims = int(solver%mesh%par%n_offset, i8)
    count_dims = int(dims, i8)

    call reader_session%read_data("stats_sample_count", self%sample_count)

    ! First moments
    call reader_session%read_data("stats_umean", self%umean, &
                                  start_dims=start_dims, &
                                  count_dims=count_dims)
    call reader_session%read_data("stats_vmean", self%vmean, &
                                  start_dims=start_dims, &
                                  count_dims=count_dims)
    call reader_session%read_data("stats_wmean", self%wmean, &
                                  start_dims=start_dims, &
                                  count_dims=count_dims)

    ! Second moments
    call reader_session%read_data("stats_uumean", self%uumean, &
                                  start_dims=start_dims, &
                                  count_dims=count_dims)
    call reader_session%read_data("stats_vvmean", self%vvmean, &
                                  start_dims=start_dims, &
                                  count_dims=count_dims)
    call reader_session%read_data("stats_wwmean", self%wwmean, &
                                  start_dims=start_dims, &
                                  count_dims=count_dims)
    call reader_session%read_data("stats_uvmean", self%uvmean, &
                                  start_dims=start_dims, &
                                  count_dims=count_dims)
    call reader_session%read_data("stats_uwmean", self%uwmean, &
                                  start_dims=start_dims, &
                                  count_dims=count_dims)
    call reader_session%read_data("stats_vwmean", self%vwmean, &
                                  start_dims=start_dims, &
                                  count_dims=count_dims)

    ! Pressure mean
    if (allocated(self%pmean)) then
      call reader_session%read_data("stats_pmean", self%pmean, &
                                    start_dims=start_dims, &
                                    count_dims=count_dims)
    end if

    ! Scalar statistics
    do is = 1, self%nspecies
      write (field_name, '(A,I0)') 'stats_phimean_', is
      call reader_session%read_data( &
        trim(field_name), self%phimean(:, :, :, is), &
        start_dims=start_dims, count_dims=count_dims &
        )
      write (field_name, '(A,I0)') 'stats_phiphimean_', is
      call reader_session%read_data( &
        trim(field_name), self%phiphimean(:, :, :, is), &
        start_dims=start_dims, count_dims=count_dims &
        )
    end do

    if (myrank == 0) then
      print *, 'Stats checkpoint restored: sample_count=', self%sample_count
    end if
  end subroutine read_checkpoint

  subroutine finalise(self)
    !! Deallocate all accumulator arrays.
    class(stats_manager_t), intent(inout) :: self

    if (allocated(self%umean)) deallocate (self%umean)
    if (allocated(self%vmean)) deallocate (self%vmean)
    if (allocated(self%wmean)) deallocate (self%wmean)
    if (allocated(self%uumean)) deallocate (self%uumean)
    if (allocated(self%vvmean)) deallocate (self%vvmean)
    if (allocated(self%wwmean)) deallocate (self%wwmean)
    if (allocated(self%uvmean)) deallocate (self%uvmean)
    if (allocated(self%uwmean)) deallocate (self%uwmean)
    if (allocated(self%vwmean)) deallocate (self%vwmean)
    if (allocated(self%pmean)) deallocate (self%pmean)
    if (allocated(self%phimean)) deallocate (self%phimean)
    if (allocated(self%phiphimean)) deallocate (self%phiphimean)
  end subroutine finalise

end module m_stats
