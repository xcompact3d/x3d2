module m_alm
  !! Actuator Line Model (ALM) turbine forcing (iturbine = 1).
  !!
  !! Ported from Incompact3d's actuator line model (G. Deskos, 2015-2017;
  !! `acl_*.f90`, `airfoils.f90`). Each blade (and the optional tower) is a
  !! line of actuator elements; the blade-element loads computed from the
  !! airfoil polars (m_alm_turbine, m_alm_element) are projected back onto
  !! the flow as a momentum source.
  !!
  !! Once per time step, on the first substage:
  !!   1. the fluid velocity at every element centre \(\mathbf{x}_s\) is
  !!      interpolated trilinearly from the grid,
  !!   2. the element forces \(\mathbf{F}_s\) follow from the polars,
  !!   3. they are smeared onto the grid with an isotropic Gaussian kernel,
  !!      \[
  !!      \mathbf{f}(\mathbf{x}) = -\frac{1}{\rho}\sum_s \mathbf{F}_s\,
  !!      \frac{1}{\epsilon^3\pi^{3/2}}
  !!      \exp\left(-\frac{|\mathbf{x}-\mathbf{x}_s|^2}{\epsilon^2}\right),
  !!      \qquad \epsilon = \epsilon_f\,(\Delta x\Delta y\Delta z)^{1/3},
  !!      \]
  !!      truncated at \(|\mathbf{x}-\mathbf{x}_s| = 10\epsilon\),
  !!   4. the rotors are advanced to the next time step.
  !! The same sources are projected onto the momentum RHS at every substage
  !! of the step, as Incompact3d adds its force field. The sampling
  !! and the projection run on the backend; the blade bookkeeping is
  !! replicated on every rank.
  !!
  !! Turbines are described by Incompact3d `.turb` files (see m_alm_turbine).
  !! The standalone actuator lines of Incompact3d (`NActuatorlines`) are not
  !! ported. Only uniform grids are supported.
  use m_mpi, only: MPI_COMM_WORLD, MPI_IN_PLACE, MPI_SUM, MPI_Allreduce, &
                   MPI_Bcast

  use m_allocator, only: allocator_t
  use m_base_backend, only: base_backend_t
  use m_common, only: dp, pi, MPI_X3D2_DP, VERT
  use m_field, only: field_t
  use m_io_session, only: reader_session_t, writer_session_t
  use m_logger, only: log_info, abort_run, str
  use m_mesh, only: mesh_t
  use m_scalar_series, only: scalar_series_t
  use m_turbine_model, only: turbine_model_t
  use m_alm_controller, only: max_ctrl_blades
  use m_alm_element, only: deg2rad
  use m_alm_output, only: init_loads_file, write_line_loads
  use m_alm_turbine, only: turbine_t, read_turbine_file, &
                           set_turbine_geometry, compute_turbine_forces, &
                           advance_turbine, upstream_probe, OPER_CONSTANT

  implicit none

  private
  public :: alm_t

  !! Kernel truncation radius in units of epsilon (Incompact3d's
  !! dist2_max = 100 epsilon^2)
  real(dp), parameter :: r_cut_eps = 10._dp

  type, extends(turbine_model_t) :: alm_t
    ! Set by the case before init
    character(len=256), allocatable :: turbine_files(:)
    real(dp) :: rho_air = 1._dp
    real(dp) :: eps_factor = 2._dp   !! kernel width in cell sizes
    real(dp) :: nu = 0._dp           !! kinematic viscosity (element Re)
    real(dp) :: u_inflow = 0._dp     !! initial reference speed (NREL mode)
    integer :: stats_start = 1
    integer :: stats_freq = 1

    integer :: n_turb = 0
    type(turbine_t), allocatable :: turbine(:)
    integer :: n_src = 0              !! actuator points over all turbines
    integer :: n_tower_elem = 0
    real(dp) :: eps, r_cut, prefactor
    real(dp) :: dt
    real(dp) :: ctime = 0._dp         !! rotor time, as Incompact3d's ctime
    real(dp) :: last_update_time = -huge(1._dp)
    integer :: current_iter = 0
    logical :: new_step = .true.
    integer :: sample_count = 0
    real(dp), allocatable :: thrust_mean(:), torque_mean(:), power_mean(:)
    !! Sources of the current step: positions and acceleration coefficients
    real(dp), allocatable :: src_points(:, :), src_coeffs(:, :)
    type(scalar_series_t), allocatable :: perf(:)
    class(base_backend_t), pointer :: backend => null()
    type(mesh_t), pointer :: mesh => null()
    type(allocator_t), pointer :: host_allocator => null()
  contains
    procedure :: init => init_alm
    procedure :: update => update_alm
    procedure :: project_forces => project_forces_alm
    procedure :: write_output => write_output_alm
    procedure :: setup_output => setup_output_alm
    procedure :: write_checkpoint => write_checkpoint_alm
    procedure :: read_checkpoint => read_checkpoint_alm
    procedure :: finalise => finalise_alm
    procedure :: gather_points
    procedure :: sample_velocity
  end type alm_t

contains

  subroutine init_alm(self, backend, mesh, host_allocator, dt)
    implicit none
    class(alm_t), intent(inout) :: self
    class(base_backend_t), target, intent(in) :: backend
    type(mesh_t), target, intent(in) :: mesh
    type(allocator_t), target, intent(in) :: host_allocator
    real(dp), intent(in) :: dt

    integer :: t

    self%backend => backend
    self%mesh => mesh
    self%host_allocator => host_allocator
    self%dt = dt

    if (.not. allocated(self%turbine_files)) &
      call abort_run('ALM: no turbine files given (alm_turbines).', &
                     collective=.true.)
    if (size(self%turbine_files) < 1) &
      call abort_run('ALM: no turbine files given (alm_turbines).', &
                     collective=.true.)
    if (self%rho_air <= 0._dp) &
      call abort_run('ALM: rho_air must be positive.', collective=.true.)
    if (self%eps_factor <= 0._dp) &
      call abort_run('ALM: eps_factor must be positive.', collective=.true.)
    if (self%nu <= 0._dp) &
      call abort_run('ALM: the viscosity must be positive.', &
                     collective=.true.)
    if (any(mesh%geo%stretching /= 'uniform')) &
      call abort_run('ALM: stretched grids are not supported.', &
                     collective=.true.)

    call log_info('Actuator line model enabled')
    self%n_turb = size(self%turbine_files)
    call log_info('Number of turbines: '//str(self%n_turb))
    allocate (self%turbine(self%n_turb))
    do t = 1, self%n_turb
      call read_turbine_file(self%turbine(t), self%turbine_files(t), &
                             self%u_inflow)
      call set_turbine_geometry(self%turbine(t), restart=.false.)
    end do

    self%n_src = 0
    self%n_tower_elem = 0
    do t = 1, self%n_turb
      associate (tb => self%turbine(t))
        self%n_src = self%n_src + tb%n_blades*tb%blade(1)%n_elem
        if (tb%has_tower) &
          self%n_tower_elem = self%n_tower_elem + tb%tower%n_elem
      end associate
    end do
    self%n_src = self%n_src + self%n_tower_elem

    self%eps = self%eps_factor*product(mesh%geo%d)**(1._dp/3._dp)
    self%r_cut = r_cut_eps*self%eps
    self%prefactor = 1._dp/(self%eps**3*pi**1.5_dp)
    call log_info('ALM: '//str(self%n_src)//' actuator points, epsilon = ' &
                  //str(self%eps, '(es12.5)'))

    allocate (self%thrust_mean(self%n_turb), self%torque_mean(self%n_turb), &
              self%power_mean(self%n_turb))
    self%thrust_mean = 0._dp
    self%torque_mean = 0._dp
    self%power_mean = 0._dp

    allocate (self%src_points(3, self%n_src), self%src_coeffs(3, self%n_src))
  end subroutine init_alm

  subroutine update_alm(self, t, dt)
    implicit none
    class(alm_t), intent(inout) :: self
    real(dp), intent(in) :: t, dt

    ! The forces are rebuilt on the first substage of each time step only.
    self%new_step = t /= self%last_update_time
    self%last_update_time = t
    self%dt = dt
    self%current_iter = nint(t/dt)
  end subroutine update_alm

  subroutine project_forces_alm(self, du, dv, dw, u, v, w)
    implicit none
    class(alm_t), intent(inout) :: self
    class(field_t), intent(inout) :: du, dv, dw
    class(field_t), intent(in) :: u, v, w

    real(dp) :: vel(3, self%n_src), tower_noise(self%n_tower_elem)
    real(dp) :: probe(3, 1), probe_vel(3, 1), stat_inc
    integer :: t, ib, n, m, ierr
    logical :: accumulate_stats

    if (self%new_step) then
      ! 1. Sample the flow at the actuator points
      call self%gather_points(self%src_points)
      call check_points_in_domain(self%mesh, self%src_points)
      call self%sample_velocity(self%src_points, vel, u, v, w)

      ! Tower lift noise, drawn on root so every rank applies the same force
      if (self%n_tower_elem > 0) then
        if (self%mesh%par%is_root()) call random_number(tower_noise)
        call MPI_Bcast(tower_noise, self%n_tower_elem, MPI_X3D2_DP, 0, &
                       MPI_COMM_WORLD, ierr)
      end if

      ! 2. Element and rotor loads
      n = 0
      m = 0
      do t = 1, self%n_turb
        associate (tb => self%turbine(t))
          do ib = 1, tb%n_blades
            tb%blade(ib)%evel = vel(:, n + 1:n + tb%blade(ib)%n_elem)
            n = n + tb%blade(ib)%n_elem
          end do
          if (tb%has_tower) then
            tb%tower%evel = vel(:, n + 1:n + tb%tower%n_elem)
            n = n + tb%tower%n_elem
          end if
          if (tb%has_tower) then
            call compute_turbine_forces(tb, self%rho_air, self%nu, &
                                        self%ctime, &
                                        tower_noise(m + 1:m + tb%tower%n_elem))
            m = m + tb%tower%n_elem
          else
            call compute_turbine_forces(tb, self%rho_air, self%nu, &
                                        self%ctime, tower_noise(1:0))
          end if
        end associate
      end do

      accumulate_stats = self%current_iter >= self%stats_start &
                         .and. mod(self%current_iter - self%stats_start, &
                                   self%stats_freq) == 0
      if (accumulate_stats) then
        self%sample_count = self%sample_count + 1
        stat_inc = 1._dp/real(self%sample_count, dp)
        do t = 1, self%n_turb
          self%thrust_mean(t) = self%thrust_mean(t) + stat_inc &
                                *(self%turbine(t)%thrust - self%thrust_mean(t))
          self%torque_mean(t) = self%torque_mean(t) + stat_inc &
                                *(self%turbine(t)%torque - self%torque_mean(t))
          self%power_mean(t) = self%power_mean(t) + stat_inc &
                               *(self%turbine(t)%power - self%power_mean(t))
        end do
      end if

      ! 3. Kernel coefficients of the acceleration, held for the step
      call gather_forces(self, self%src_coeffs)
      self%src_coeffs = -self%src_coeffs*self%prefactor/self%rho_air

      ! 4. Advance the rotors to the next time step
      self%ctime = self%last_update_time
      do t = 1, self%n_turb
        if (self%turbine(t)%oper /= OPER_CONSTANT) then
          probe(:, 1) = nearest_grid_point(self%mesh, &
                                           upstream_probe(self%turbine(t)))
          call self%sample_velocity(probe, probe_vel, u, v, w)
          self%turbine(t)%u_upstream = probe_vel(:, 1)
        end if
        call advance_turbine(self%turbine(t), self%ctime, self%dt)
      end do
    end if

    ! Smear the step's sources onto the RHS. Projecting at every substage
    ! touches only the grid points near the rotors, so it is cheaper than
    ! keeping a full force field and adding it each substage.
    call self%backend%field_add_gaussian_sources( &
      du, dv, dw, self%src_points, self%src_coeffs, self%eps, self%r_cut)
    self%new_step = .false.
  end subroutine project_forces_alm

  subroutine gather_points(self, points)
    !! Element centres of all blades and towers, in the source order used
    !! throughout: turbine by turbine, blades then tower.
    implicit none
    class(alm_t), intent(in) :: self
    real(dp), intent(out) :: points(:, :)

    integer :: t, ib, n, ne

    n = 0
    do t = 1, self%n_turb
      associate (tb => self%turbine(t))
        do ib = 1, tb%n_blades
          ne = tb%blade(ib)%n_elem
          points(:, n + 1:n + ne) = tb%blade(ib)%pe
          n = n + ne
        end do
        if (tb%has_tower) then
          ne = tb%tower%n_elem
          points(:, n + 1:n + ne) = tb%tower%pe
          n = n + ne
        end if
      end associate
    end do
  end subroutine gather_points

  subroutine gather_forces(self, forces)
    implicit none
    class(alm_t), intent(in) :: self
    real(dp), intent(out) :: forces(:, :)

    integer :: t, ib, n, ne

    n = 0
    do t = 1, self%n_turb
      associate (tb => self%turbine(t))
        do ib = 1, tb%n_blades
          ne = tb%blade(ib)%n_elem
          forces(:, n + 1:n + ne) = tb%blade(ib)%ef
          n = n + ne
        end do
        if (tb%has_tower) then
          ne = tb%tower%n_elem
          forces(:, n + 1:n + ne) = tb%tower%ef
          n = n + ne
        end if
      end associate
    end do
  end subroutine gather_forces

  subroutine sample_velocity(self, points, vel, u, v, w)
    !! Trilinear interpolation of (u, v, w) at `points` over all ranks.
    implicit none
    class(alm_t), intent(inout) :: self
    real(dp), intent(in) :: points(:, :)
    real(dp), intent(out) :: vel(:, :)
    class(field_t), intent(in) :: u, v, w

    integer :: ierr

    call self%backend%field_sample_points(vel, points, u, v, w)
    call MPI_Allreduce(MPI_IN_PLACE, vel, size(vel), MPI_X3D2_DP, MPI_SUM, &
                       MPI_COMM_WORLD, ierr)
  end subroutine sample_velocity

  subroutine check_points_in_domain(mesh, points)
    type(mesh_t), intent(in) :: mesh
    real(dp), intent(in) :: points(:, :)

    integer :: s

    do s = 1, size(points, 2)
      if (any(points(:, s) < 0._dp) .or. any(points(:, s) > mesh%geo%L)) &
        call abort_run('ALM: actuator point ('//str(points(1, s), '(es12.5)') &
                       //', '//str(points(2, s), '(es12.5)')//', ' &
                       //str(points(3, s), '(es12.5)') &
                       //') lies outside the domain.', collective=.true.)
    end do
  end subroutine check_points_in_domain

  function nearest_grid_point(mesh, point) result(node)
    !! The VERT grid point closest to `point` on a uniform grid.
    type(mesh_t), intent(in) :: mesh
    real(dp), intent(in) :: point(3)
    real(dp) :: node(3)

    integer :: gdims(3), dir, g

    gdims = mesh%get_global_dims(VERT)
    do dir = 1, 3
      g = nint(point(dir)/mesh%geo%d(dir))
      if (mesh%grid%periodic_BC(dir)) then
        g = modulo(g, gdims(dir))
      else
        g = min(max(g, 0), gdims(dir) - 1)
      end if
      node(dir) = g*mesh%geo%d(dir)
    end do
  end function nearest_grid_point

  subroutine write_output_alm(self, iter, is_root)
    implicit none
    class(alm_t), intent(inout) :: self
    integer, intent(in) :: iter
    logical, intent(in) :: is_root

    real(dp) :: time
    integer :: t, ib, n, ne

    if (.not. is_root) return
    time = real(iter, dp)*self%dt
    n = 0
    do t = 1, self%n_turb
      associate (tb => self%turbine(t))
        call log_info('ALM turbine '//trim(tb%name)//' iter '//str(iter) &
                      //': CT = '//str(tb%ct, '(es12.5)') &
                      //'  CP = '//str(tb%cp, '(es12.5)') &
                      //'  thrust = '//str(tb%thrust, '(es12.5)') &
                      //'  torque = '//str(tb%torque, '(es12.5)') &
                      //'  power = '//str(tb%power, '(es12.5)'))
        if (allocated(self%perf)) then
          call self%perf(t)%write_step(time, [ &
                                       tb%azim_angle/(2._dp*pi), &
                                       tb%angular_vel, tb%delta_omega, &
                                       tb%controller%gen_speed, &
                                       tb%controller%gen_trq, &
                                       -tb%cbp/deg2rad, &
                                       tb%u_upstream, tb%uref, &
                                       tb%thrust, tb%torque, tb%power, &
                                       tb%ct, tb%cp, &
                                       self%thrust_mean(t), &
                                       self%torque_mean(t), &
                                       self%power_mean(t)])
        end if
        ! src_points holds the positions sampled this step, in source order
        do ib = 1, tb%n_blades
          ne = tb%blade(ib)%n_elem
          call write_line_loads(tb%blade(ib), &
                                self%src_points(:, n + 1:n + ne), &
                                iter, time, is_root)
          n = n + ne
        end do
        if (tb%has_tower) then
          ne = tb%tower%n_elem
          call write_line_loads(tb%tower, self%src_points(:, n + 1:n + ne), &
                                iter, time, is_root)
          n = n + ne
        end if
      end associate
    end do
  end subroutine write_output_alm

  subroutine setup_output_alm(self, is_root, append)
    class(alm_t), intent(inout) :: self
    logical, intent(in) :: is_root, append

    character(len=16), parameter :: columns(18) = &
                                    ['revolutions     ', &
                                     'omega           ', &
                                     'delta_omega     ', &
                                     'gen_speed       ', &
                                     'gen_torque      ', &
                                     'pitch_deg       ', &
                                     'u_upstream      ', &
                                     'v_upstream      ', &
                                     'w_upstream      ', &
                                     'u_ref           ', &
                                     'thrust          ', &
                                     'torque          ', &
                                     'power           ', &
                                     'CT              ', &
                                     'CP              ', &
                                     'thrust_mean     ', &
                                     'torque_mean     ', &
                                     'power_mean      ']
    integer :: t, ib

    allocate (self%perf(self%n_turb))
    do t = 1, self%n_turb
      associate (tb => self%turbine(t))
        call self%perf(t)%init(trim(tb%name)//'.perf', columns, is_root, &
                               append)
        do ib = 1, tb%n_blades
          call init_loads_file(tb%blade(ib), is_root, append)
        end do
        if (tb%has_tower) call init_loads_file(tb%tower, is_root, append)
      end associate
    end do
  end subroutine setup_output_alm

  subroutine write_checkpoint_alm(self, writer)
    class(alm_t), intent(inout) :: self
    type(writer_session_t), intent(inout) :: writer

    integer :: t, k, ib, ie

    call writer%write_data('alm_n_turb', self%n_turb)
    call writer%write_data('alm_sample_count', self%sample_count)
    call writer%write_data('alm_ctime', self%ctime)
    do t = 1, self%n_turb
      associate (tb => self%turbine(t), c => self%turbine(t)%controller)
        call put('azim_angle', tb%azim_angle)
        call put('angular_vel', tb%angular_vel)
        call put('delta_omega', tb%delta_omega)
        call put('cbp', tb%cbp)
        call put('cbp_old', tb%cbp_old)
        call put('uref', tb%uref)
        call put('thrust_mean', self%thrust_mean(t))
        call put('torque_mean', self%torque_mean(t))
        call put('power_mean', self%power_mean(t))
        call put('ctrl_gen_speed', c%gen_speed)
        call put('ctrl_gen_speed_f', c%gen_speed_f)
        call put('ctrl_gen_trq', c%gen_trq)
        call put('ctrl_gk', c%gk)
        call put('ctrl_int_spd_err', c%int_spd_err)
        call put('ctrl_last_gen_trq', c%last_gen_trq)
        call put('ctrl_last_time', c%last_time)
        call put('ctrl_last_time_pc', c%last_time_pc)
        call put('ctrl_last_time_vs', c%last_time_vs)
        do k = 1, max_ctrl_blades
          call put('ctrl_bl_pitch_'//str(k), c%bl_pitch(k))
          call put('ctrl_pit_com_'//str(k), c%pit_com(k))
        end do
        call writer%write_data('alm_ctrl_status_'//str(t), c%status)
        ! The end corrections use the previous step's angles of attack
        do ib = 1, tb%n_blades
          do ie = 1, tb%blade(ib)%n_elem
            call put('eaoa_'//str(ib)//'_'//str(ie), tb%blade(ib)%eaoa(ie))
          end do
        end do
      end associate
    end do

  contains

    subroutine put(name, val)
      character(len=*), intent(in) :: name
      real(dp), intent(in) :: val
      call writer%write_data('alm_'//name//'_'//str(t), val)
    end subroutine put
  end subroutine write_checkpoint_alm

  subroutine read_checkpoint_alm(self, reader)
    class(alm_t), intent(inout) :: self
    type(reader_session_t), intent(inout) :: reader

    integer :: checkpoint_n_turb, t, k, ib, ie

    call reader%read_data('alm_n_turb', checkpoint_n_turb)
    if (checkpoint_n_turb /= self%n_turb) &
      call abort_run('ALM: checkpoint turbine count does not match &
                     &configuration.')
    call reader%read_data('alm_sample_count', self%sample_count)
    call reader%read_data('alm_ctime', self%ctime)
    do t = 1, self%n_turb
      associate (tb => self%turbine(t), c => self%turbine(t)%controller)
        call get('azim_angle', tb%azim_angle)
        call get('angular_vel', tb%angular_vel)
        call get('delta_omega', tb%delta_omega)
        call get('cbp', tb%cbp)
        call get('cbp_old', tb%cbp_old)
        call get('uref', tb%uref)
        call get('thrust_mean', self%thrust_mean(t))
        call get('torque_mean', self%torque_mean(t))
        call get('power_mean', self%power_mean(t))
        call get('ctrl_gen_speed', c%gen_speed)
        call get('ctrl_gen_speed_f', c%gen_speed_f)
        call get('ctrl_gen_trq', c%gen_trq)
        call get('ctrl_gk', c%gk)
        call get('ctrl_int_spd_err', c%int_spd_err)
        call get('ctrl_last_gen_trq', c%last_gen_trq)
        call get('ctrl_last_time', c%last_time)
        call get('ctrl_last_time_pc', c%last_time_pc)
        call get('ctrl_last_time_vs', c%last_time_vs)
        do k = 1, max_ctrl_blades
          call get('ctrl_bl_pitch_'//str(k), c%bl_pitch(k))
          call get('ctrl_pit_com_'//str(k), c%pit_com(k))
        end do
        call reader%read_data('alm_ctrl_status_'//str(t), c%status)

        ! Rebuild the blades at the restored azimuth and pitch
        call set_turbine_geometry(tb, restart=.true.)
        do ib = 1, tb%n_blades
          do ie = 1, tb%blade(ib)%n_elem
            call get('eaoa_'//str(ib)//'_'//str(ie), tb%blade(ib)%eaoa(ie))
          end do
        end do
      end associate
    end do

  contains

    subroutine get(name, val)
      character(len=*), intent(in) :: name
      real(dp), intent(inout) :: val
      call reader%read_data('alm_'//name//'_'//str(t), val)
    end subroutine get
  end subroutine read_checkpoint_alm

  subroutine finalise_alm(self)
    class(alm_t), intent(inout) :: self

    integer :: t

    if (allocated(self%perf)) then
      do t = 1, size(self%perf)
        call self%perf(t)%finalise()
      end do
      deallocate (self%perf)
    end if
    if (allocated(self%src_points)) deallocate (self%src_points)
    if (allocated(self%src_coeffs)) deallocate (self%src_coeffs)
    if (allocated(self%turbine)) deallocate (self%turbine)
  end subroutine finalise_alm

end module m_alm
