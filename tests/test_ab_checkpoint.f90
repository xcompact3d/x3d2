program test_ab_checkpoint
  use iso_fortran_env, only: stderr => error_unit
  use mpi
  use m_common, only: dp, DIR_X, VERT
  use m_mesh, only: mesh_t
  use m_allocator, only: allocator_t
  use m_base_backend, only: base_backend_t
  use m_field, only: flist_t
  use m_time_integrator, only: time_intg_t, init
  use m_checkpoint_manager, only: checkpoint_manager_t
  use m_solver, only: solver_t
  use m_omp_backend, only: omp_backend_t
  use m_omp_common, only: SZ

  implicit none

  type old_history_t
    real(dp), allocatable :: data(:, :, :)
  end type old_history_t

  type(mesh_t), allocatable :: mesh
  type(allocator_t), target :: host_allocator
  class(allocator_t), pointer :: allocator
  type(omp_backend_t), target :: omp_backend
  class(base_backend_t), pointer :: backend
  type(solver_t) :: solver_cont, solver_restart
  type(checkpoint_manager_t) :: chk_mgr_write, chk_mgr_restart
  type(flist_t), allocatable :: curr(:)
  type(flist_t), allocatable :: deriv(:)
  integer :: ierr, irank, nproc
  integer :: dims_global(3), nproc_dir(3)
  real(dp) :: L_global(3)
  character(len=20) :: BC_x(2), BC_y(2), BC_z(2)
  type(old_history_t), allocatable :: olds_ref(:)
  integer :: iter, idx, checkpoint_iter, nolds_total
  real(dp), parameter :: dt = 1.0e-3_dp
  real(dp), parameter :: tol = 1.0e-15_dp
  character(len=*), parameter :: ckpt_prefix = 'ab3_test_checkpoint'
  logical :: allpass
  real(dp) :: diff

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, irank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, nproc, ierr)

  allpass = .true.

  dims_global = [1, 1, 1]
  nproc_dir = [1, 1, 1]
  L_global = [1.0_dp, 1.0_dp, 1.0_dp]
  BC_x = ['periodic', 'periodic']
  BC_y = ['periodic', 'periodic']
  BC_z = ['periodic', 'periodic']

  mesh = mesh_t(dims_global, nproc_dir, L_global, BC_x, BC_y, BC_z)

  host_allocator = allocator_t(mesh%get_dims(VERT), SZ)
  allocator => host_allocator
  omp_backend = omp_backend_t(mesh, allocator)
  backend => omp_backend

  call init_solver(solver_cont, backend, mesh, allocator, dt)
  call init_checkpoint_config(chk_mgr_write, ckpt_prefix, 5)

  call setup_curr(curr, solver_cont)
  call allocate_derivs(allocator, deriv)

  checkpoint_iter = 5
  do iter = 1, checkpoint_iter
    call fill_derivs(deriv, iter)
    call solver_cont%time_integrator%step(curr, deriv, solver_cont%dt)
    solver_cont%current_iter = iter
    if (iter == checkpoint_iter) then
      call chk_mgr_write%handle_checkpoint_step(solver_cont, iter, &
                                                MPI_COMM_WORLD)
    end if
  end do

  call capture_old_history(solver_cont, olds_ref)

  call release_derivs(allocator, deriv)
  deallocate (curr)

  call solver_cont%time_integrator%finalize()
  call release_solver_fields(solver_cont)

  call init_solver(solver_restart, backend, mesh, allocator, dt)
  call init_checkpoint_config(chk_mgr_restart, ckpt_prefix, 5)
  chk_mgr_restart%config%restart_from_checkpoint = .true.
  chk_mgr_restart%config%restart_file = trim(ckpt_prefix)//'_000005.bp'

  call chk_mgr_restart%handle_restart(solver_restart, MPI_COMM_WORLD)

  nolds_total = solver_restart%time_integrator%nolds*solver_restart%nvars
  idx = 0
  do iter = 1, solver_restart%nvars
    do checkpoint_iter = 1, solver_restart%time_integrator%nolds
      idx = idx + 1
      ! Each olds(i,j) must match bit-for-bit - this is the regression guard
      diff = maxval(abs( &
                    solver_restart%time_integrator%olds( &
                    iter, checkpoint_iter)%ptr%data - olds_ref(idx)%data))
      if (diff > tol .and. irank == 0) then
        write (stderr, '(a,i0,a,i0,a,es12.5)') 'Mismatch in olds(', iter, &
          ',', checkpoint_iter, ') diff=', diff
        allpass = .false.
      end if
    end do
  end do

  call solver_restart%time_integrator%finalize()
  call release_solver_fields(solver_restart)

  call cleanup_checkpoint_files(ckpt_prefix, irank)
  if (allocated(olds_ref)) then
    do idx = 1, size(olds_ref)
      if (allocated(olds_ref(idx)%data)) deallocate (olds_ref(idx)%data)
    end do
    deallocate (olds_ref)
  end if

  call MPI_Finalize(ierr)

  if (irank == 0) then
    if (allpass) then
      write (stderr, '(a)') 'AB checkpoint test passed.'
    else
      error stop 'AB checkpoint test failed.'
    end if
  end if

contains

  subroutine init_solver(solver, backend, mesh, allocator, dt)
    type(solver_t), intent(inout) :: solver
    class(base_backend_t), pointer, intent(inout) :: backend
    type(mesh_t), target, intent(inout) :: mesh
    class(allocator_t), pointer, intent(inout) :: allocator
    real(dp), intent(in) :: dt

    solver%backend => backend
    solver%mesh => mesh
    solver%host_allocator => allocator
    solver%nvars = 3
    solver%nspecies = 0
    solver%dt = dt
    solver%n_iters = 0
    solver%n_output = 0
    solver%current_iter = 0
    solver%ngrid = product(mesh%get_global_dims(VERT))
    solver%u => allocator%get_block(DIR_X)
    solver%v => allocator%get_block(DIR_X)
    solver%w => allocator%get_block(DIR_X)
    call solver%u%set_data_loc(VERT)
    call solver%v%set_data_loc(VERT)
    call solver%w%set_data_loc(VERT)
    call solver%u%fill(0.0_dp)
    call solver%v%fill(0.0_dp)
    call solver%w%fill(0.0_dp)
    solver%time_integrator = init(backend, allocator, 'AB3', solver%nvars)
  end subroutine init_solver

  subroutine release_solver_fields(solver)
    type(solver_t), intent(inout) :: solver
    if (associated(solver%backend)) then
      if (associated(solver%u)) &
        call solver%backend%allocator%release_block(solver%u)
      if (associated(solver%v)) &
        call solver%backend%allocator%release_block(solver%v)
      if (associated(solver%w)) &
        call solver%backend%allocator%release_block(solver%w)
    end if
    nullify (solver%u)
    nullify (solver%v)
    nullify (solver%w)
  end subroutine release_solver_fields

  subroutine init_checkpoint_config(manager, prefix, freq)
    type(checkpoint_manager_t), intent(inout) :: manager
    character(len=*), intent(in) :: prefix
    integer, intent(in) :: freq

    manager%config%checkpoint_freq = freq
    manager%config%snapshot_freq = 0
    manager%config%keep_checkpoint = .true.
    manager%config%checkpoint_prefix = prefix
    manager%config%snapshot_prefix = 'snapshot'
    manager%config%restart_from_checkpoint = .false.
    manager%config%output_stride = [1, 1, 1]
    manager%last_checkpoint_step = -1
    manager%full_resolution = [1, 1, 1]
  end subroutine init_checkpoint_config

  subroutine setup_curr(list, solver)
    type(flist_t), allocatable, intent(out) :: list(:)
    type(solver_t), intent(inout) :: solver

    allocate (list(3))
    list(1)%ptr => solver%u
    list(2)%ptr => solver%v
    list(3)%ptr => solver%w
  end subroutine setup_curr

  subroutine allocate_derivs(allocator, deriv)
    class(allocator_t), pointer, intent(inout) :: allocator
    type(flist_t), allocatable, intent(out) :: deriv(:)
    integer :: i

    allocate (deriv(3))
    do i = 1, 3
      deriv(i)%ptr => allocator%get_block(DIR_X)
      call deriv(i)%ptr%fill(0.0_dp)
    end do
  end subroutine allocate_derivs

  subroutine release_derivs(allocator, deriv)
    class(allocator_t), pointer, intent(inout) :: allocator
    type(flist_t), allocatable, intent(inout) :: deriv(:)
    integer :: i

    if (.not. allocated(deriv)) return
    do i = 1, size(deriv)
      if (associated(deriv(i)%ptr)) call allocator%release_block(deriv(i)%ptr)
    end do
    deallocate (deriv)
  end subroutine release_derivs

  subroutine fill_derivs(deriv, istep)
    type(flist_t), intent(inout) :: deriv(:)
    integer, intent(in) :: istep
    integer :: i
    real(dp) :: base
    do i = 1, size(deriv)
      base = real(istep, dp) + 0.1_dp*real(i - 1, dp)
      deriv(i)%ptr%data = base
    end do
  end subroutine fill_derivs

  subroutine capture_old_history(solver, olds_ref)
    type(solver_t), intent(in) :: solver
    type(old_history_t), allocatable, intent(out) :: olds_ref(:)
    integer :: nolds_total, i, j, idx
    integer :: old_dims(3)

    nolds_total = solver%time_integrator%nolds*solver%nvars
    allocate (olds_ref(nolds_total))
    idx = 0
    do i = 1, solver%nvars
      do j = 1, solver%time_integrator%nolds
        idx = idx + 1
        old_dims = shape(solver%time_integrator%olds(i, j)%ptr%data)
        allocate (olds_ref(idx)%data(old_dims(1), old_dims(2), old_dims(3)))
        olds_ref(idx)%data = solver%time_integrator%olds(i, j)%ptr%data
      end do
    end do
  end subroutine capture_old_history

  subroutine cleanup_checkpoint_files(prefix, irank)
    character(len=*), intent(in) :: prefix
    integer, intent(in) :: irank
    integer :: ierr
    if (irank /= 0) return
    call execute_command_line('rm -rf '//trim(prefix)//'_*.bp', exitstat=ierr)
    call execute_command_line('rm -f '//trim(prefix)//'_temp.bp', &
                              exitstat=ierr)
  end subroutine cleanup_checkpoint_files

end program test_ab_checkpoint
