module m_stats
  !! Online accumulation of flow field statistics.
  !!
  !! Nine running means are maintained using the incremental update:
  !! \[ \bar{x}_n = \bar{x}_{n-1} + \frac{x_n - \bar{x}_{n-1}}{n} \]
  !! applied to \(u, v, w, u^2, v^2, w^2, uv, uw, vw\).
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
  use m_common, only: dp, i8, DIR_C, VERT, get_argument
  use m_config, only: stats_config_t
  use m_field, only: field_t
  use m_solver, only: solver_t
  use m_io_session, only: writer_session_t

  implicit none

  private
  public :: stats_manager_t

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
  contains
    procedure :: init
    procedure :: update
    procedure :: write_stats
    procedure :: handle_restart
    procedure :: finalise
  end type stats_manager_t

contains

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

    class(field_t), pointer :: u_tmp, v_tmp, w_tmp
    integer :: dims(3)
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
      self%umean  = self%umean  + (u   - self%umean)  * stat_inc
      self%vmean  = self%vmean  + (v   - self%vmean)  * stat_inc
      self%wmean  = self%wmean  + (w   - self%wmean)  * stat_inc
      self%uumean = self%uumean + (u*u - self%uumean) * stat_inc
      self%vvmean = self%vvmean + (v*v - self%vvmean) * stat_inc
      self%wwmean = self%wwmean + (w*w - self%wwmean) * stat_inc
      self%uvmean = self%uvmean + (u*v - self%uvmean) * stat_inc
      self%uwmean = self%uwmean + (u*w - self%uwmean) * stat_inc
      self%vwmean = self%vwmean + (v*w - self%vwmean) * stat_inc
    end associate

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
    integer :: dims(3), myrank, ierr
    real(dp), allocatable :: uprime(:, :, :)
    real(dp), allocatable :: vprime(:, :, :)
    real(dp), allocatable :: wprime(:, :, :)
    real(dp), allocatable :: uv(:, :, :)
    real(dp), allocatable :: uw(:, :, :)
    real(dp), allocatable :: vw(:, :, :)

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
    call writer_session%write_data("umean",  self%umean, &
                                   shape_dims, start_dims, count_dims)
    call writer_session%write_data("vmean",  self%vmean, &
                                   shape_dims, start_dims, count_dims)
    call writer_session%write_data("wmean",  self%wmean, &
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

    call writer_session%close()

    deallocate (uprime, vprime, wprime, uv, uw, vw)
  end subroutine write_stats

  subroutine handle_restart(self, solver, comm)
    !! Zero-initialise all accumulators on restart, so statistics accumulate
    !! from scratch rather than resuming from a previous run.
    class(stats_manager_t), intent(inout) :: self
    class(solver_t), intent(in) :: solver
    integer, intent(in) :: comm

    integer :: myrank, ierr

    if (.not. self%is_active) return

    call MPI_Comm_rank(comm, myrank, ierr)
    if (myrank == 0) then
      print *, 'Statistics: re-initialising accumulators for restart'
    end if

    self%sample_count = 0
    self%umean  = 0.0_dp
    self%vmean  = 0.0_dp
    self%wmean  = 0.0_dp
    self%uumean = 0.0_dp
    self%vvmean = 0.0_dp
    self%wwmean = 0.0_dp
    self%uvmean = 0.0_dp
    self%uwmean = 0.0_dp
    self%vwmean = 0.0_dp
  end subroutine handle_restart

  subroutine finalise(self)
    !! Deallocate all accumulator arrays.
    class(stats_manager_t), intent(inout) :: self

    if (allocated(self%umean))  deallocate (self%umean)
    if (allocated(self%vmean))  deallocate (self%vmean)
    if (allocated(self%wmean))  deallocate (self%wmean)
    if (allocated(self%uumean)) deallocate (self%uumean)
    if (allocated(self%vvmean)) deallocate (self%vvmean)
    if (allocated(self%wwmean)) deallocate (self%wwmean)
    if (allocated(self%uvmean)) deallocate (self%uvmean)
    if (allocated(self%uwmean)) deallocate (self%uwmean)
    if (allocated(self%vwmean)) deallocate (self%vwmean)
  end subroutine finalise

end module m_stats
