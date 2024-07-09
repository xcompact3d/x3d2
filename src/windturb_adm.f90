module m_windturb_adm
  use iso_fortran_env, only: stderr => error_unit
  use mpi

  use m_allocator, only: allocator_t, field_t
  use m_base_backend, only: base_backend_t
  use m_common, only: dp, pi, DIR_X, DIR_C, VERT
  use m_field, only: field_t
  use m_mesh, only: mesh_t

  implicit none

  type :: actuator_disc_t
    integer :: id
    real(dp) :: coords(3) ! center of rotation
    real(dp) :: yaw ! rotor yaw angle, in radians w.r.t. y-axis
    real(dp) :: tilt ! rotor tilt angle, in radians w.r.t. z-axis
    real(dp) :: D, area ! actuator disk diameter and area
    real(dp) :: C_T ! thrust coefficient
    real(dp) :: alpha ! induction coefficient
    real(dp) :: rot_N(3) ! axis of rotation
    real(dp) :: U_disc, U_disc_filt ! disk-averaged speed and filtered speed
    real(dp) :: power, thrust ! instantaneous quantities
    real(dp) :: U_disc_tavg, power_tavg, thrust_tavg ! time averaged quantities
    class(field_t), pointer :: gamma_disc
  end type actuator_disc_t

  type :: windturb_adm_t
    integer :: n_turb
    real(dp) :: rho_air, cell_vol
    type(actuator_disc_t), allocatable :: disc(:)
    class(base_backend_t), pointer :: backend
    class(mesh_t), pointer :: mesh
    type(allocator_t), pointer :: host_allocator
  contains
    procedure :: compute_sources
  end type windturb_adm_t

  interface windturb_adm_t
    module procedure init
  end interface windturb_adm_t

contains

  function init(backend, mesh, host_allocator, disc_params) &
    result(windturb_adm)
    implicit none

    class(base_backend_t), target, intent(inout) :: backend
    type(mesh_t), target, intent(inout) :: mesh
    type(allocator_t), target, intent(inout) :: host_allocator
    real(dp), dimension(:, :), intent(in) :: disc_params
    type(windturb_adm_t) :: windturb_adm

    class(field_t), pointer :: gamma_host
    real(dp) :: disc_coords(3), yaw, tilt, D, rot_N(3)
    integer :: t, i, j, k, dims(3), ierr
    real(dp) :: coords(3), delta, dx, dy, dz, delta_N, delta_R, disc_thick, &
                gamma_val, gamma_tot

    windturb_adm%backend => backend
    windturb_adm%mesh => mesh
    windturb_adm%host_allocator => host_allocator

    windturb_adm%n_turb = size(disc_params, dim=1)
    windturb_adm%rho_air = 1._dp
    windturb_adm%cell_vol = product(mesh%geo%d)

    allocate (windturb_adm%disc(windturb_adm%n_turb))

    dims = mesh%get_dims(VERT)

    gamma_host => windturb_adm%host_allocator%get_block(DIR_C)

    do t = 1, windturb_adm%n_turb
      disc_coords(:) = disc_params(t, 1:3)
      yaw = disc_params(t, 4)*pi/180._dp
      tilt = disc_params(t, 5)*pi/180._dp
      D = disc_params(t, 6)
      rot_N(1) = cos(yaw)*cos(tilt)
      rot_N(2) = sin(tilt)
      rot_N(3) = sin(yaw)
      windturb_adm%disc(t)%coords(:) = disc_coords(:)
      windturb_adm%disc(t)%yaw = yaw
      windturb_adm%disc(t)%tilt = tilt
      windturb_adm%disc(t)%D = D
      windturb_adm%disc(t)%C_T = disc_params(t, 7)
      windturb_adm%disc(t)%alpha = disc_params(t, 8)
      windturb_adm%disc(t)%rot_N(:) = rot_N(:)
      windturb_adm%disc(t)%area = pi*(D**2)/4._dp

      gamma_host%data(:, :, :) = 0._dp
      gamma_tot = 0._dp

      delta = sqrt((mesh%geo%d(1)*rot_N(1))**2 &
                   + (mesh%geo%d(2)*rot_N(2))**2 &
                   + (mesh%geo%d(3)*rot_N(3))**2)
      disc_thick = max(D/8._dp, delta*1.5_dp)
      do k = 1, dims(3)
        do j = 1, dims(2)
          do i = 1, dims(1)
            coords = mesh%get_coordinates(i, j, k)
            dx = coords(1) - disc_coords(1)
            dy = coords(2) - disc_coords(2)
            dz = coords(3) - disc_coords(3)
            delta_N = dx*rot_N(1) + dy*rot_N(2) - dz*rot_N(3)

            dx = coords(1) - delta_N*rot_N(1)
            dy = coords(2) - delta_N*rot_N(2)
            dz = coords(3) + delta_N*rot_N(3)
            delta_R = sqrt((dx - disc_coords(1))**2 &
                           + (dy - disc_coords(2))**2 &
                           + (dz - disc_coords(3))**2)

            gamma_val = exp(-((delta_N/(disc_thick/2._dp))**2 &
                              + (delta_R/(D/2._dp))**8))
            gamma_host%data(i, j, k) = gamma_val
            gamma_tot = gamma_tot + gamma_val
          end do
        end do
      end do

      ! normalise the gamma field
      call MPI_Allreduce(MPI_IN_PLACE, gamma_tot, 1, MPI_DOUBLE_PRECISION, &
                         MPI_SUM, MPI_COMM_WORLD, ierr)
      gamma_host%data(:, :, :) = gamma_host%data(:, :, :)/gamma_tot

      windturb_adm%disc(t)%gamma_disc &
        => windturb_adm%backend%allocator%get_block(DIR_X)
      call windturb_adm%backend%set_field_data( &
        windturb_adm%disc(t)%gamma_disc, gamma_host%data &
        )
    end do

    call windturb_adm%host_allocator%release_block(gamma_host)

  end function init

  subroutine compute_sources(self, Fx, Fy, Fz, u, v, w)
    implicit none

    class(windturb_adm_t) :: self
    class(field_t), intent(inout) :: Fx, Fy, Fz
    class(field_t), intent(in) :: u, v, w

    real(dp) :: u_avg, v_avg, w_avg, rot_N(3)
    real(dp) :: coeff_x, coeff_y, coeff_z, C_T_prime
    real(dp), allocatable, dimension(:) :: vel_avg
    integer :: t, ierr

    allocate (vel_avg(self%n_turb))

    do t = 1, self%n_turb
      rot_N(:) = self%disc(t)%rot_N(:)
      u_avg = self%backend%scalar_product(self%disc(t)%gamma_disc, u)
      v_avg = self%backend%scalar_product(self%disc(t)%gamma_disc, v)
      w_avg = self%backend%scalar_product(self%disc(t)%gamma_disc, w)
      vel_avg(t) = u_avg*rot_N(1) + v_avg*rot_N(2) - w_avg*rot_N(3)
    end do

    ! obtain the average velocity combining all the subdomains
    call MPI_Allreduce(vel_avg, self%disc%U_disc, self%n_turb, &
                       MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

    deallocate (vel_avg)

    ! no filtering for the moment
    self%disc(:)%U_disc_filt = self%disc(:)%U_disc

    ! use a hack to zeroise the force fields
    call self%backend%vecadd(0._dp, Fx, 0._dp, Fx)
    call self%backend%vecadd(0._dp, Fy, 0._dp, Fy)
    call self%backend%vecadd(0._dp, Fz, 0._dp, Fz)

    do t = 1, self%n_turb
      C_T_prime = self%disc(t)%C_T/(1._dp - self%disc(t)%alpha)**2
      self%disc(t)%thrust = 0.5_dp*self%rho_air*C_T_prime &
                            *self%disc(t)%U_disc_filt**2*self%disc(t)%area
      self%disc(t)%power = self%disc(t)%thrust*self%disc(t)%U_disc_filt
      coeff_x = -self%disc(t)%thrust*self%disc(t)%rot_N(1)/self%cell_vol
      coeff_y = -self%disc(t)%thrust*self%disc(t)%rot_N(2)/self%cell_vol
      coeff_z = self%disc(t)%thrust*self%disc(t)%rot_N(3)/self%cell_vol
      call self%backend%vecadd(coeff_x, self%disc(t)%gamma_disc, 1._dp, Fx)
      call self%backend%vecadd(coeff_y, self%disc(t)%gamma_disc, 1._dp, Fy)
      call self%backend%vecadd(coeff_z, self%disc(t)%gamma_disc, 1._dp, Fz)
    end do

  end subroutine compute_sources

end module m_windturb_adm
