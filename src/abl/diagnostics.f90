module m_abl_diagnostics
  !! Horizontal/time profile diagnostics for the neutral ABL validation.
  !!
  !! This is deliberately separate from `abl_t`: the ABL driver owns physics,
  !! while this case-owned object owns validation state and CSV output.
  use mpi, only: MPI_COMM_WORLD, MPI_Allreduce, MPI_IN_PLACE, MPI_INTEGER, &
                 MPI_SUM

  use m_allocator, only: allocator_t
  use m_base_backend, only: base_backend_t
  use m_common, only: dp, DIR_C, MPI_X3D2_DP, VERT
  use m_config, only: abl_config_t
  use m_field, only: field_t
  use m_mesh, only: mesh_t

  implicit none

  private
  public :: abl_diagnostics_t, neutral_log_law, friction_velocity

  type :: abl_diagnostics_t
    class(base_backend_t), pointer :: backend => null()
    type(mesh_t), pointer :: mesh => null()
    type(allocator_t), pointer :: host_allocator => null()
    type(abl_config_t) :: cfg
    real(dp), allocatable :: profile_sum(:, :)
    real(dp) :: wall_stress_sum(2) = 0._dp
    real(dp) :: first_time = 0._dp
    integer :: sample_count = 0
    integer :: last_iter = -1
  contains
    procedure :: sample
    procedure, private :: write_profile
  end type abl_diagnostics_t

  interface abl_diagnostics_t
    module procedure init
  end interface abl_diagnostics_t

contains

  function init(backend, mesh, host_allocator, cfg) result(diagnostics)
    class(base_backend_t), target, intent(inout) :: backend
    type(mesh_t), target, intent(inout) :: mesh
    type(allocator_t), target, intent(inout) :: host_allocator
    type(abl_config_t), intent(in) :: cfg
    type(abl_diagnostics_t) :: diagnostics
    integer :: dims(3)

    diagnostics%backend => backend
    diagnostics%mesh => mesh
    diagnostics%host_allocator => host_allocator
    diagnostics%cfg = cfg

    if (cfg%profile_start_iter < 0) return
    if (mesh%par%nproc_dir(2) /= 1) &
      error stop 'ABL profile diagnostics require an undecomposed y direction.'
    dims = mesh%get_global_dims(VERT)
    allocate (diagnostics%profile_sum(3, dims(2)), source=0._dp)
  end function init

  pure real(dp) function neutral_log_law(y, u_star, kappa, z0) result(u_log)
    real(dp), intent(in) :: y, u_star, kappa, z0

    if (y <= 0._dp .or. z0 <= 0._dp .or. kappa <= 0._dp) then
      u_log = 0._dp
    else
      u_log = u_star/kappa*log(y/z0)
    end if
  end function neutral_log_law

  pure real(dp) function friction_velocity(tau_x, tau_z) result(u_star)
    !! Friction velocity from the magnitude of a mean wall-stress vector.
    real(dp), intent(in) :: tau_x, tau_z

    u_star = sqrt(sqrt(tau_x**2 + tau_z**2))
  end function friction_velocity

  subroutine sample(self, iter, time, u, v, w)
    class(abl_diagnostics_t), intent(inout) :: self
    integer, intent(in) :: iter
    real(dp), intent(in) :: time
    class(field_t), intent(in) :: u, v, w

    class(field_t), pointer :: hu, hv, hw
    real(dp), allocatable :: profile(:, :)
    real(dp) :: wall_stress(2), drag_coeff
    real(dp) :: u_sample, w_sample, speed
    integer :: dims(3), plane_count, ierr
    integer :: i, j, k

    if (self%cfg%profile_start_iter < 0) return
    if (iter < self%cfg%profile_start_iter .or. iter == self%last_iter) return

    dims = self%mesh%get_dims(VERT)
    allocate (profile(3, dims(2)), source=0._dp)

    hu => self%host_allocator%get_block(DIR_C, VERT)
    hv => self%host_allocator%get_block(DIR_C, VERT)
    hw => self%host_allocator%get_block(DIR_C, VERT)
    call self%backend%get_field_data(hu%data, u)
    call self%backend%get_field_data(hv%data, v)
    call self%backend%get_field_data(hw%data, w)

    do k = 1, dims(3)
      do j = 1, dims(2)
        do i = 1, dims(1)
          profile(1, j) = profile(1, j) + hu%data(i, j, k)
          profile(2, j) = profile(2, j) + hv%data(i, j, k)
          profile(3, j) = profile(3, j) + hw%data(i, j, k)
        end do
      end do
    end do

    plane_count = dims(1)*dims(3)
    call MPI_Allreduce(MPI_IN_PLACE, profile, size(profile), MPI_X3D2_DP, &
                       MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_Allreduce(MPI_IN_PLACE, plane_count, 1, MPI_INTEGER, &
                       MPI_SUM, MPI_COMM_WORLD, ierr)
    profile = profile/real(plane_count, dp)

    drag_coeff = (self%cfg%kappa/log( &
                  0.5_dp*self%mesh%geo%d(2)/self%cfg%z0))**2
    u_sample = 0.5_dp*(profile(1, 1) + profile(1, 2))
    w_sample = 0.5_dp*(profile(3, 1) + profile(3, 2))
    speed = sqrt(u_sample**2 + w_sample**2)
    wall_stress(1) = -drag_coeff*u_sample*speed
    wall_stress(2) = -drag_coeff*w_sample*speed

    self%sample_count = self%sample_count + 1
    if (self%sample_count == 1) self%first_time = time
    self%profile_sum = self%profile_sum + profile
    self%wall_stress_sum = self%wall_stress_sum + wall_stress
    self%last_iter = iter

    call self%write_profile(time)

    call self%host_allocator%release_block(hu)
    call self%host_allocator%release_block(hv)
    call self%host_allocator%release_block(hw)
    deallocate (profile)
  end subroutine sample

  subroutine write_profile(self, time)
    class(abl_diagnostics_t), intent(in) :: self
    real(dp), intent(in) :: time

    real(dp) :: profile(3, size(self%profile_sum, 2))
    real(dp) :: tau(2), diagnosed_u_star, y, u_log
    integer :: j, unit

    if (.not. self%mesh%par%is_root()) return

    profile = self%profile_sum/real(self%sample_count, dp)
    tau = self%wall_stress_sum/real(self%sample_count, dp)
    diagnosed_u_star = friction_velocity(tau(1), tau(2))

    open (newunit=unit, file=trim(self%cfg%profile_file), &
          status='replace', action='write')
    write (unit, '(A,I0)') '# sample_count = ', self%sample_count
    write (unit, '(A,ES20.12)') '# time_start = ', self%first_time
    write (unit, '(A,ES20.12)') '# time_end = ', time
    write (unit, '(A,ES20.12)') '# imposed_u_star = ', self%cfg%u_star
    write (unit, '(A,ES20.12)') '# diagnosed_u_star = ', diagnosed_u_star
    write (unit, '(A,ES20.12)') '# mean_tau_x = ', tau(1)
    write (unit, '(A,ES20.12)') '# mean_tau_z = ', tau(2)
    write (unit, '(A)') '# y,u_mean,v_mean,w_mean,u_log'
    do j = 1, size(profile, 2)
      y = self%mesh%geo%vert_coords(j, 2)
      u_log = neutral_log_law( &
              y, self%cfg%u_star, self%cfg%kappa, self%cfg%z0)
      write (unit, '(ES20.12,4(",",ES20.12))') &
        y, profile(1, j), profile(2, j), profile(3, j), u_log
    end do
    close (unit)
  end subroutine write_profile

end module m_abl_diagnostics
