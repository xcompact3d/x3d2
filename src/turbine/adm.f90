module m_adm
  !! Actuator Disc Model (ADM) turbine forcing (iturbine = 2).
  !!
  !! Each turbine is a disc that applies a thrust force computed from a thrust
  !! coefficient and the disc-averaged inflow velocity. The force is
  !! regularised and projected onto the Eulerian grid with a per-disc kernel
  !! \(\gamma\),
  !! built once at init and normalised to sum to one. This operation is also
  !! commonly called force spreading or smearing.
  !!
  !! The kernels of all discs share two fields, so device memory does not grow
  !! with the number of discs: `gamma` holds the kernel values and `disc_id`
  !! labels each grid point with the disc whose kernel is largest there (zero
  !! where every kernel vanishes). Kernels of separate discs must therefore
  !! not overlap; init stops if a significant part of any kernel lies where
  !! another disc's kernel is larger. The disc-averaged velocities of all
  !! discs come from one labelled weighted sum and one MPI reduction, and the
  !! projected forces of all discs are added to the momentum RHS in one
  !! labelled axpy.
  !!
  !! For yaw angle \(\psi\) and tilt angle \(\theta\), the rotor-normal vector
  !! used by the model is
  !! \[
  !! \mathbf{n}
  !! = \left(\cos\psi\cos\theta,\ \sin\theta,\ -\sin\psi\right).
  !! \]
  !! For a grid point \(\mathbf{x}\) and disc centre \(\mathbf{x}_c\), the
  !! normal and radial distances are
  !! \[
  !! \delta_N = (\mathbf{x}-\mathbf{x}_c)\cdot\mathbf{n},
  !! \qquad
  !! \delta_R =
  !! \left\lVert
  !! (\mathbf{x}-\mathbf{x}_c)-\delta_N\mathbf{n}
  !! \right\rVert_2.
  !! \]
  !! The unnormalised anisotropic regularisation kernel is Gaussian in the
  !! disc-normal direction and eighth-order super-Gaussian in the radial
  !! direction:
  !! \[
  !! \widetilde{\gamma}(\mathbf{x})
  !! = \exp\left[
  !! -\left(
  !! \frac{\delta_N^2}{(h/2)^2}
  !! + \frac{\delta_R^8}{(D/2)^8}
  !! \right)
  !! \right],
  !! \qquad
  !! h = \max\left(\frac{D}{8},\,1.5\Delta_n\right),
  !! \]
  !! where
  !! \[
  !! \Delta_n =
  !! \sqrt{(\Delta x\,n_x)^2+(\Delta y\,n_y)^2+(\Delta z\,n_z)^2}.
  !! \]
  !! The discrete kernel is normalised globally over all MPI ranks:
  !! \[
  !! \gamma_i =
  !! \frac{\widetilde{\gamma}_i}
  !! {\sum_j\widetilde{\gamma}_j},
  !! \qquad
  !! \sum_i\gamma_i=1.
  !! \]
  !! Consequently, the instantaneous disc-normal velocity is the weighted
  !! global average
  !! \[
  !! U_d = \sum_i \gamma_i\,\mathbf{u}_i\cdot\mathbf{n}.
  !! \]
  !! an optional first-order low-pass filter is applied:
  !! \[
  !! \beta = \frac{\Delta t/T_{\mathrm{relax}}}
  !! {1+\Delta t/T_{\mathrm{relax}}},
  !! \qquad
  !! U_F^{n} = \beta U_d^{n} + (1-\beta)U_F^{n-1}.
  !! \]
  !! On the first update, or when \(T_{\mathrm{relax}}<0\), \(U_F=U_d\).
  !! The local thrust coefficient, thrust, and power are
  !! \[
  !! C_T' = \frac{C_T}{(1-\alpha)^2},
  !! \qquad
  !! T = \frac{1}{2}\rho C_T' A U_F^2,
  !! \qquad
  !! P = T U_F,
  !! \qquad
  !! A = \frac{\pi D^2}{4}.
  !! \]
  !! Finally, the acceleration added to the momentum right-hand side at grid
  !! point \(i\) is
  !! \[
  !! \mathbf{f}_i =
  !! -\frac{T}{\rho\,\Delta V}\,\gamma_i\,\mathbf{n}.
  !! \]
  !!
  use iso_fortran_env, only: stderr => error_unit
  use m_mpi, only: MPI_COMM_WORLD, MPI_IN_PLACE, MPI_SUM, MPI_Allreduce

  use m_allocator, only: allocator_t
  use m_base_backend, only: base_backend_t
  use m_common, only: dp, pi, MPI_X3D2_DP, DIR_X, DIR_C, VERT
  use m_field, only: field_t
  use m_io_session, only: reader_session_t, writer_session_t
  use m_mesh, only: mesh_t
  use m_scalar_series, only: scalar_series_t
  use m_turbine_model, only: turbine_model_t

  implicit none

  type :: actuator_disc_t
    real(dp) :: coords(3)         !! disc centre
    real(dp) :: yaw, tilt         !! rotor yaw/tilt angles [rad]
    real(dp) :: D, area           !! diameter and swept area
    real(dp) :: thick             !! kernel thickness h along rot_N
    real(dp) :: C_T, alpha        !! thrust and induction coefficients
    real(dp) :: rot_N(3)          !! unit rotor-axis vector
    real(dp) :: U_disc = 0._dp       !! instantaneous disc-averaged speed
    real(dp) :: U_disc_filt = 0._dp  !! filtered disc-averaged speed
    real(dp) :: thrust = 0._dp    !! instantaneous thrust
    real(dp) :: power = 0._dp     !! instantaneous power
    real(dp) :: U_disc_mean = 0._dp
    real(dp) :: thrust_mean = 0._dp
    real(dp) :: power_mean = 0._dp
  end type actuator_disc_t

  type, extends(turbine_model_t) :: adm_t
    integer :: n_turb = 0
    real(dp) :: rho_air = 1._dp
    real(dp) :: T_relax = -1._dp
    real(dp) :: cell_vol
    real(dp) :: dt
    real(dp) :: last_update_time = -huge(1._dp)
    integer :: current_iter = 0
    integer :: stats_start = 1
    integer :: stats_freq = 1
    integer :: sample_count = 0
    logical :: filter_initialized = .false.
    logical :: recompute_forces = .true.
    character(len=256) :: coords_file = ''
    type(actuator_disc_t), allocatable :: disc(:)
    class(field_t), pointer :: gamma => null()   !! smearing kernels (DIR_X)
    class(field_t), pointer :: disc_id => null() !! disc label (DIR_X)
    type(scalar_series_t), allocatable :: output(:)
    class(base_backend_t), pointer :: backend => null()
    type(mesh_t), pointer :: mesh => null()
    type(allocator_t), pointer :: host_allocator => null()
  contains
    procedure :: init => init_adm
    procedure :: update => update_adm
    procedure :: project_forces => project_forces_adm
    procedure :: write_output => write_output_adm
    procedure :: setup_output => setup_output_adm
    procedure :: write_checkpoint => write_checkpoint_adm
    procedure :: read_checkpoint => read_checkpoint_adm
    procedure :: finalise => finalise_adm
  end type adm_t

contains

  subroutine init_adm(self, backend, mesh, host_allocator, dt)
    implicit none
    class(adm_t), intent(inout) :: self
    class(base_backend_t), target, intent(in) :: backend
    type(mesh_t), target, intent(in) :: mesh
    type(allocator_t), target, intent(in) :: host_allocator
    real(dp), intent(in) :: dt

    real(dp), parameter :: overlap_tol = 1.e-6_dp
    real(dp), allocatable :: disc_params(:, :), kernel_sums(:, :)
    class(field_t), pointer :: gamma_host, disc_id_host
    integer :: t, t_max, i, j, k, dims(3), ierr
    real(dp) :: coords(3), delta, gamma_val, gamma_max
    real(dp) :: yaw, tilt, D, rot_N(3)

    self%backend => backend
    self%mesh => mesh
    self%host_allocator => host_allocator
    self%dt = dt
    self%cell_vol = product(mesh%geo%d)

    if (self%rho_air <= 0._dp) then
      error stop 'ADM: rho_air must be positive.'
    end if
    if (self%T_relax == 0._dp) then
      error stop 'ADM: T_relax must be positive or negative to disable.'
    end if
    if (any(mesh%geo%stretching /= 'uniform')) then
      error stop 'ADM: stretched grids are not supported.'
    end if

    call read_ad_file(self%coords_file, disc_params)
    self%n_turb = size(disc_params, dim=1)
    allocate (self%disc(self%n_turb))

    if (mesh%par%is_root()) then
      print *, 'Actuator disc model enabled'
      print *, 'Number of actuator discs:', self%n_turb
    end if

    do t = 1, self%n_turb
      yaw = disc_params(t, 4)*pi/180._dp
      tilt = disc_params(t, 5)*pi/180._dp
      D = disc_params(t, 6)
      if (D <= 0._dp) error stop 'ADM: rotor diameter must be positive.'
      if (abs(1._dp - disc_params(t, 8)) <= epsilon(1._dp)) then
        error stop 'ADM: induction coefficient alpha must not equal one.'
      end if
      rot_N(1) = cos(yaw)*cos(tilt)
      rot_N(2) = sin(tilt)
      rot_N(3) = sin(yaw)

      self%disc(t)%coords(:) = disc_params(t, 1:3)
      self%disc(t)%yaw = yaw
      self%disc(t)%tilt = tilt
      self%disc(t)%D = D
      self%disc(t)%C_T = disc_params(t, 7)
      self%disc(t)%alpha = disc_params(t, 8)
      self%disc(t)%rot_N(:) = rot_N(:)
      self%disc(t)%area = pi*(D**2)/4._dp
      delta = sqrt((mesh%geo%d(1)*rot_N(1))**2 &
                   + (mesh%geo%d(2)*rot_N(2))**2 &
                   + (mesh%geo%d(3)*rot_N(3))**2)
      self%disc(t)%thick = max(D/8._dp, delta*1.5_dp)
    end do

    ! Label every grid point with the disc whose kernel is largest there and
    ! keep that kernel value. kernel_sums(1, t) sums the full kernel of disc t
    ! and kernel_sums(2, t) only the part kept where disc t is the largest.
    dims = mesh%get_dims(VERT)
    gamma_host => self%host_allocator%get_block(DIR_C)
    disc_id_host => self%host_allocator%get_block(DIR_C)
    gamma_host%data(:, :, :) = 0._dp
    disc_id_host%data(:, :, :) = 0._dp
    allocate (kernel_sums(2, self%n_turb))
    kernel_sums = 0._dp
    do k = 1, dims(3)
      do j = 1, dims(2)
        do i = 1, dims(1)
          coords = mesh%get_coordinates(i, j, k)
          gamma_max = 0._dp
          t_max = 0
          do t = 1, self%n_turb
            gamma_val = disc_kernel(self%disc(t), coords)
            kernel_sums(1, t) = kernel_sums(1, t) + gamma_val
            if (gamma_val > gamma_max) then
              gamma_max = gamma_val
              t_max = t
            end if
          end do
          gamma_host%data(i, j, k) = gamma_max
          disc_id_host%data(i, j, k) = real(t_max, dp)
          if (t_max > 0) then
            kernel_sums(2, t_max) = kernel_sums(2, t_max) + gamma_max
          end if
        end do
      end do
    end do

    call MPI_Allreduce(MPI_IN_PLACE, kernel_sums, size(kernel_sums), &
                       MPI_X3D2_DP, MPI_SUM, MPI_COMM_WORLD, ierr)
    do t = 1, self%n_turb
      if (kernel_sums(1, t) <= tiny(1._dp)) then
        error stop 'ADM: disc kernel does not overlap the computational mesh.'
      end if
      if (kernel_sums(1, t) - kernel_sums(2, t) &
          > overlap_tol*kernel_sums(1, t)) then
        if (mesh%par%is_root()) then
          write (stderr, '(A,I0,A,ES10.3,A)') ' ADM: disc ', t, &
            ' loses a fraction ', 1._dp - kernel_sums(2, t)/kernel_sums(1, t), &
            ' of its kernel to neighbouring discs.'
        end if
        error stop 'ADM: disc kernels overlap; move the discs further apart.'
      end if
    end do

    ! Normalise each disc's kernel to sum to one across all ranks. The
    ! disc-averaged velocity is then the gamma-weighted sum of u over the
    ! points labelled with that disc.
    do k = 1, dims(3)
      do j = 1, dims(2)
        do i = 1, dims(1)
          t = nint(disc_id_host%data(i, j, k))
          if (t > 0) then
            gamma_host%data(i, j, k) = gamma_host%data(i, j, k) &
                                       /kernel_sums(2, t)
          end if
        end do
      end do
    end do

    self%gamma => self%backend%allocator%get_block(DIR_X, VERT)
    self%disc_id => self%backend%allocator%get_block(DIR_X, VERT)
    call self%backend%set_field_data(self%gamma, gamma_host%data)
    call self%backend%set_field_data(self%disc_id, disc_id_host%data)

    call self%host_allocator%release_block(gamma_host)
    call self%host_allocator%release_block(disc_id_host)
  end subroutine init_adm

  subroutine update_adm(self, t, dt)
    implicit none
    class(adm_t), intent(inout) :: self
    real(dp), intent(in) :: t, dt

    ! Advance the ADM time filter once per physical timestep (first substage).
    ! On repeated substages the filter is frozen, but U_disc and the thrust are
    ! still recomputed from the current substage velocity.
    self%recompute_forces = t /= self%last_update_time
    self%last_update_time = t
    self%dt = dt
    self%current_iter = nint(t/dt)
  end subroutine update_adm

  subroutine project_forces_adm(self, du, dv, dw, u, v, w)
    implicit none
    class(adm_t), intent(inout) :: self
    class(field_t), intent(inout) :: du, dv, dw
    class(field_t), intent(in) :: u, v, w

    integer :: t, ierr
    real(dp) :: disc_avgs(3, self%n_turb), coeffs(3, self%n_turb)
    real(dp) :: rot_N(3), C_T_prime, filter_weight, stat_inc, U_thrust
    logical :: accumulate_stats

    accumulate_stats = self%recompute_forces &
                       .and. self%current_iter >= self%stats_start &
                       .and. mod(self%current_iter - self%stats_start, &
                                 self%stats_freq) == 0
    if (accumulate_stats) then
      self%sample_count = self%sample_count + 1
      stat_inc = 1._dp/real(self%sample_count, dp)
    end if

    ! Always recompute U_disc from the current stage velocity so the force
    ! is consistent with each RK substage. The disc-weighted u, v and w of
    ! every disc are summed locally in one pass, then reduced together in a
    ! single MPI call.
    call self%backend%field_labelled_sums(disc_avgs, self%gamma, &
                                          self%disc_id, u, v, w)
    call MPI_Allreduce(MPI_IN_PLACE, disc_avgs, size(disc_avgs), &
                       MPI_X3D2_DP, MPI_SUM, MPI_COMM_WORLD, ierr)

    do t = 1, self%n_turb
      rot_N(:) = self%disc(t)%rot_N(:)
      self%disc(t)%U_disc = disc_avgs(1, t)*rot_N(1) &
                            + disc_avgs(2, t)*rot_N(2) &
                            - disc_avgs(3, t)*rot_N(3)

      if (self%recompute_forces) then
        ! First substage: advance the time filter.
        if (.not. self%filter_initialized .or. self%T_relax < 0._dp) then
          self%disc(t)%U_disc_filt = self%disc(t)%U_disc
        else
          filter_weight = (self%dt/self%T_relax) &
                          /(1._dp + self%dt/self%T_relax)
          self%disc(t)%U_disc_filt = filter_weight*self%disc(t)%U_disc &
                                     + (1._dp - filter_weight) &
                                     *self%disc(t)%U_disc_filt
        end if
        U_thrust = self%disc(t)%U_disc_filt
      else
        ! Subsequent RK substages: filter not advanced; use instantaneous
        ! U_disc so the thrust tracks the current stage velocity.
        U_thrust = self%disc(t)%U_disc
      end if

      ! Thrust from the local (Calaf-style) thrust coefficient.
      C_T_prime = self%disc(t)%C_T/(1._dp - self%disc(t)%alpha)**2
      self%disc(t)%thrust = 0.5_dp*self%rho_air*C_T_prime &
                            *U_thrust**2*self%disc(t)%area
      self%disc(t)%power = self%disc(t)%thrust*U_thrust

      if (accumulate_stats) then
        self%disc(t)%U_disc_mean = self%disc(t)%U_disc_mean &
                                   + stat_inc*(self%disc(t)%U_disc_filt &
                                               - self%disc(t)%U_disc_mean)
        self%disc(t)%thrust_mean = self%disc(t)%thrust_mean &
                                   + stat_inc*(self%disc(t)%thrust &
                                               - self%disc(t)%thrust_mean)
        self%disc(t)%power_mean = self%disc(t)%power_mean &
                                  + stat_inc*(self%disc(t)%power &
                                              - self%disc(t)%power_mean)
      end if

      ! Thrust acceleration (opposing the inflow) for the momentum RHS.
      ! Density remains in the thrust/power diagnostics and cancels here.
      coeffs(1, t) = -self%disc(t)%thrust*rot_N(1) &
                     /(self%rho_air*self%cell_vol)
      coeffs(2, t) = -self%disc(t)%thrust*rot_N(2) &
                     /(self%rho_air*self%cell_vol)
      coeffs(3, t) = self%disc(t)%thrust*rot_N(3) &
                     /(self%rho_air*self%cell_vol)
    end do
    call self%backend%field_labelled_axpy(coeffs, self%gamma, self%disc_id, &
                                          du, dv, dw)
    if (self%recompute_forces) self%filter_initialized = .true.
    self%recompute_forces = .false.
  end subroutine project_forces_adm

  subroutine write_output_adm(self, iter, is_root)
    implicit none
    class(adm_t), intent(inout) :: self
    integer, intent(in) :: iter
    logical, intent(in) :: is_root
    integer :: t
    real(dp) :: values(7)

    if (.not. is_root) return
    do t = 1, self%n_turb
      print '(A,I0,A,I0,3(A,ES12.5))', ' ADM disc ', t, ' iter ', iter, &
        ': U_disc_filt = ', self%disc(t)%U_disc_filt, &
        '  thrust = ', self%disc(t)%thrust, &
        '  power = ', self%disc(t)%power
      if (allocated(self%output)) then
        values = [self%disc(t)%U_disc, self%disc(t)%U_disc_filt, &
                  self%disc(t)%power, self%disc(t)%thrust, &
                  self%disc(t)%U_disc_mean, self%disc(t)%power_mean, &
                  self%disc(t)%thrust_mean]
        call self%output(t)%write_step(real(iter, dp)*self%dt, values)
      end if
    end do
  end subroutine write_output_adm

  subroutine setup_output_adm(self, is_root, append)
    class(adm_t), intent(inout) :: self
    logical, intent(in) :: is_root, append

    character(len=32) :: filename
    character(len=16), parameter :: columns(7) = &
                                    ['U_disc          ', &
                                     'U_disc_filt     ', &
                                     'power           ', &
                                     'thrust          ', &
                                     'U_disc_mean     ', &
                                     'power_mean      ', &
                                     'thrust_mean     ']
    integer :: t

    allocate (self%output(self%n_turb))
    do t = 1, self%n_turb
      write (filename, '("disc",I0,".adm")') t
      call self%output(t)%init(filename, columns, is_root, append)
    end do
  end subroutine setup_output_adm

  subroutine write_checkpoint_adm(self, writer)
    class(adm_t), intent(inout) :: self
    type(writer_session_t), intent(inout) :: writer

    character(len=64) :: name
    integer :: t

    call writer%write_data('adm_n_turb', self%n_turb)
    call writer%write_data('adm_sample_count', self%sample_count)
    do t = 1, self%n_turb
      write (name, '("adm_U_disc_",I0)') t
      call writer%write_data(trim(name), self%disc(t)%U_disc)
      write (name, '("adm_U_disc_filt_",I0)') t
      call writer%write_data(trim(name), self%disc(t)%U_disc_filt)
      write (name, '("adm_power_",I0)') t
      call writer%write_data(trim(name), self%disc(t)%power)
      write (name, '("adm_thrust_",I0)') t
      call writer%write_data(trim(name), self%disc(t)%thrust)
      write (name, '("adm_U_disc_mean_",I0)') t
      call writer%write_data(trim(name), self%disc(t)%U_disc_mean)
      write (name, '("adm_power_mean_",I0)') t
      call writer%write_data(trim(name), self%disc(t)%power_mean)
      write (name, '("adm_thrust_mean_",I0)') t
      call writer%write_data(trim(name), self%disc(t)%thrust_mean)
    end do
  end subroutine write_checkpoint_adm

  subroutine read_checkpoint_adm(self, reader)
    class(adm_t), intent(inout) :: self
    type(reader_session_t), intent(inout) :: reader

    character(len=64) :: name
    integer :: checkpoint_n_turb, t

    call reader%read_data('adm_n_turb', checkpoint_n_turb)
    if (checkpoint_n_turb /= self%n_turb) then
      error stop 'ADM: checkpoint turbine count does not match configuration.'
    end if
    call reader%read_data('adm_sample_count', self%sample_count)
    do t = 1, self%n_turb
      write (name, '("adm_U_disc_",I0)') t
      call reader%read_data(trim(name), self%disc(t)%U_disc)
      write (name, '("adm_U_disc_filt_",I0)') t
      call reader%read_data(trim(name), self%disc(t)%U_disc_filt)
      write (name, '("adm_power_",I0)') t
      call reader%read_data(trim(name), self%disc(t)%power)
      write (name, '("adm_thrust_",I0)') t
      call reader%read_data(trim(name), self%disc(t)%thrust)
      write (name, '("adm_U_disc_mean_",I0)') t
      call reader%read_data(trim(name), self%disc(t)%U_disc_mean)
      write (name, '("adm_power_mean_",I0)') t
      call reader%read_data(trim(name), self%disc(t)%power_mean)
      write (name, '("adm_thrust_mean_",I0)') t
      call reader%read_data(trim(name), self%disc(t)%thrust_mean)
    end do
    self%filter_initialized = .true.
  end subroutine read_checkpoint_adm

  subroutine finalise_adm(self)
    class(adm_t), intent(inout) :: self

    integer :: t

    if (allocated(self%output)) then
      do t = 1, size(self%output)
        call self%output(t)%finalise()
      end do
      deallocate (self%output)
    end if
    if (associated(self%gamma)) then
      call self%backend%allocator%release_block(self%gamma)
      nullify (self%gamma)
    end if
    if (associated(self%disc_id)) then
      call self%backend%allocator%release_block(self%disc_id)
      nullify (self%disc_id)
    end if
    if (allocated(self%disc)) deallocate (self%disc)
  end subroutine finalise_adm

  pure real(dp) function disc_kernel(disc, coords) result(gamma_val)
    !! Unnormalised smearing kernel of `disc` at `coords`: Gaussian in the
    !! rotor-normal direction (delta_N), super-Gaussian (^8) in the radial
    !! direction (delta_R) to give a flat-topped disc.
    implicit none
    type(actuator_disc_t), intent(in) :: disc
    real(dp), intent(in) :: coords(3)

    real(dp) :: dx, dy, dz, delta_N, delta_R

    dx = coords(1) - disc%coords(1)
    dy = coords(2) - disc%coords(2)
    dz = coords(3) - disc%coords(3)
    delta_N = dx*disc%rot_N(1) + dy*disc%rot_N(2) - dz*disc%rot_N(3)

    dx = coords(1) - delta_N*disc%rot_N(1)
    dy = coords(2) - delta_N*disc%rot_N(2)
    dz = coords(3) + delta_N*disc%rot_N(3)
    delta_R = sqrt((dx - disc%coords(1))**2 &
                   + (dy - disc%coords(2))**2 &
                   + (dz - disc%coords(3))**2)

    gamma_val = exp(-((delta_N/(disc%thick/2._dp))**2 &
                      + (delta_R/(disc%D/2._dp))**8))
  end function disc_kernel

  subroutine read_ad_file(filename, disc_params)
    !! Read turbine disc parameters from a `.ad` file. Format: one header line
    !! followed by one row per disc of
    !!   CoR(x) CoR(y) CoR(z) Yaw[deg] Tilt[deg] RotorDiam C_T alpha
    implicit none
    character(*), intent(in) :: filename
    real(dp), allocatable, intent(out) :: disc_params(:, :)

    integer :: unit, ios, n_turb, t
    character(len=512) :: line

    if (len_trim(filename) == 0) then
      error stop 'ADM: no disc coordinates file specified (adm_coords).'
    end if

    open (newunit=unit, file=trim(filename), status='old', action='read', &
          iostat=ios)
    if (ios /= 0) then
      write (stderr, *) 'ADM: cannot open disc coordinates file: ', &
        trim(filename)
      error stop
    end if

    ! Count data rows (skip the header line).
    read (unit, '(A)', iostat=ios) line  ! header
    n_turb = 0
    do
      read (unit, '(A)', iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) == 0) cycle
      n_turb = n_turb + 1
    end do

    if (n_turb == 0) error stop 'ADM: no discs found in coordinates file.'

    allocate (disc_params(n_turb, 8))
    rewind (unit)
    read (unit, '(A)', iostat=ios) line  ! header again
    do t = 1, n_turb
      read (unit, *, iostat=ios) disc_params(t, :)
      if (ios /= 0) error stop 'ADM: error reading disc coordinates row.'
    end do
    close (unit)
  end subroutine read_ad_file

end module m_adm
