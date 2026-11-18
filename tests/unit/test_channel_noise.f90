program test_channel_noise
  !! The channel initial noise must repeat exactly, must scale each velocity
  !! component with its own init_noise, must ignore inlet_noise, and must
  !! differ between two ranks in z.
   use m_mpi, only: MPI_Init, MPI_Finalize, MPI_Comm_size, MPI_Allreduce, &
                    MPI_COMM_WORLD, MPI_IN_PLACE, MPI_LOGICAL, MPI_LAND, &
                    MPI_MIN, MPI_MAX

   use m_allocator, only: allocator_t
   use m_backend_runtime, only: backend_runtime_t
   use m_base_backend, only: base_backend_t
   use m_case_channel, only: case_channel_t
   use m_common, only: dp, DIR_X, DIR_C, VERT, MPI_X3D2_DP
   use m_mesh, only: mesh_t
   use m_solver, only: solver_t

   implicit none

   type(backend_runtime_t), target :: runtime
   type(mesh_t), target :: mesh
   class(allocator_t), pointer :: allocator
   class(base_backend_t), pointer :: backend
   type(case_channel_t) :: flow_case

   integer, parameter :: dims_global(3) = [16, 17, 8]
   real(dp), parameter :: lengths(3) = [4._dp, 2._dp, 2._dp]
   real(dp), parameter :: tol = 100._dp*epsilon(1._dp)
   real(dp), parameter :: noise_a(3) = [0.1_dp, 0.2_dp, 0.3_dp]
   real(dp), parameter :: noise_b(3) = [0.05_dp, 0.05_dp, 0.05_dp]
   real(dp), parameter :: noise_off(3) = [0._dp, 0._dp, 0._dp]
   real(dp), parameter :: ratio(3) = [2._dp, 4._dp, 6._dp]
   real(dp), parameter :: inlet_on(3) = [0.5_dp, 0.5_dp, 0.5_dp]
   real(dp), parameter :: inlet_off(3) = [0._dp, 0._dp, 0._dp]
   character(len=9), parameter :: bc_periodic(2) = ['periodic ', 'periodic ']
   character(len=9), parameter :: bc_wall(2) = ['dirichlet', 'dirichlet']

   real(dp), allocatable :: a(:, :, :, :), b(:, :, :, :)
   real(dp), allocatable :: noise_a_part(:, :, :, :), noise_b_part(:, :, :, :)
   real(dp) :: checksum, checksum_min, checksum_max
   integer :: dims(3), dims_padded(3), nproc, ierr, i
   logical :: all_pass

   call MPI_Init(ierr)
   call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)
   all_pass = .true.

   ! All ranks are in z, thus each rank has the same x and y coordinates
   mesh = mesh_t(dims_global, [1, 1, nproc], lengths, &
                 bc_periodic, bc_wall, bc_periodic)
   dims = mesh%get_dims(VERT)

   call runtime%init(mesh)
   allocator => runtime%allocator
   backend => runtime%backend

   ! The initial condition uses only these members of the solver
   allocate (solver_t :: flow_case%solver)
   flow_case%solver%mesh => mesh
   flow_case%solver%backend => backend
   flow_case%solver%host_allocator => runtime%host_allocator
   flow_case%solver%u => allocator%get_block(DIR_X, VERT)
   flow_case%solver%v => allocator%get_block(DIR_X, VERT)
   flow_case%solver%w => allocator%get_block(DIR_X, VERT)

   dims_padded = allocator%get_padded_dims(DIR_C)
   allocate (a(dims_padded(1), dims_padded(2), dims_padded(3), 3))
   allocate (b, mold=a)

   call initial_field(noise_a, inlet_on, a)
   call initial_field(noise_a, inlet_on, b)
   call check('a second initialisation repeats the field', &
              maxval(abs(active(a) - active(b))) == 0._dp, all_pass)

   call initial_field(noise_a, inlet_off, b)
   call check('inlet_noise does not change the field', &
              maxval(abs(active(a) - active(b))) == 0._dp, all_pass)

   ! Without noise, u is the laminar profile and v, w are zero. Subtract
   ! this field to get the noise of each component.
   call initial_field(noise_off, inlet_on, b)
   noise_a_part = active(a) - active(b)
   call check('each component has noise', &
              minval(maxval(maxval(maxval(abs(noise_a_part), 1), 1), 1)) &
              > 0._dp, all_pass)

   ! The draws do not depend on init_noise, thus the noise scales with it
   call initial_field(noise_b, inlet_on, a)
   noise_b_part = active(a) - active(b)
   do i = 1, 3
      call check('component '//achar(48 + i)//' scales with its init_noise', &
                 maxval(abs(noise_a_part(:, :, :, i) &
                            - ratio(i)*noise_b_part(:, :, :, i))) <= tol, &
                 all_pass)
   end do

   ! v is noise only. Ranks with the same seed have the same sum.
   if (nproc > 1) then
      checksum = sum(noise_a_part(:, :, :, 2))
      call MPI_Allreduce(checksum, checksum_min, 1, MPI_X3D2_DP, &
                         MPI_MIN, MPI_COMM_WORLD, ierr)
      call MPI_Allreduce(checksum, checksum_max, 1, MPI_X3D2_DP, &
                         MPI_MAX, MPI_COMM_WORLD, ierr)
      call check('the ranks have different noise', &
                 checksum_max > checksum_min, all_pass)
   end if

   call allocator%release_block(flow_case%solver%u)
   call allocator%release_block(flow_case%solver%v)
   call allocator%release_block(flow_case%solver%w)

   call MPI_Allreduce(MPI_IN_PLACE, all_pass, 1, MPI_LOGICAL, MPI_LAND, &
                      MPI_COMM_WORLD, ierr)
   if (.not. all_pass) error stop 'FAIL'
   print *, 'PASS'
   call MPI_Finalize(ierr)

contains

   subroutine initial_field(init_noise, inlet_noise, fields)
    !! Set the channel initial condition and return u, v, w.
      real(dp), intent(in) :: init_noise(3), inlet_noise(3)
      real(dp), intent(out) :: fields(:, :, :, :)

      flow_case%channel_cfg%init_noise = init_noise
      flow_case%channel_cfg%inlet_noise = inlet_noise
      call flow_case%initial_conditions()
      call backend%get_field_data(fields(:, :, :, 1), flow_case%solver%u)
      call backend%get_field_data(fields(:, :, :, 2), flow_case%solver%v)
      call backend%get_field_data(fields(:, :, :, 3), flow_case%solver%w)
   end subroutine initial_field

   function active(fields) result(inner)
    !! The unpadded part of the fields
      real(dp), intent(in) :: fields(:, :, :, :)
      real(dp), allocatable :: inner(:, :, :, :)

      inner = fields(1:dims(1), 1:dims(2), 1:dims(3), :)
   end function active

   subroutine check(label, condition, pass)
      character(*), intent(in) :: label
      logical, intent(in) :: condition
      logical, intent(inout) :: pass

      if (condition) then
         print *, 'PASS: ', label
      else
         print *, 'FAIL: ', label
         pass = .false.
      end if
   end subroutine check

end program test_channel_noise
