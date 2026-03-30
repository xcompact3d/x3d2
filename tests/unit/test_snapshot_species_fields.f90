program test_snapshot_species_fields
  use iso_fortran_env, only: stderr => error_unit
  use mpi

  use m_allocator, only: allocator_t
  use m_base_backend, only: base_backend_t
  use m_common, only: dp, DIR_X, VERT
  use m_io_field_utils, only: field_ptr_t, setup_field_arrays, &
                              cleanup_field_arrays
  use m_mesh, only: mesh_t
  use m_omp_backend, only: omp_backend_t
  use m_omp_common, only: SZ
  use m_solver, only: solver_t

  implicit none

  integer, parameter :: nspecies = 2
  real(dp), parameter :: tol = 1.0e-12_dp

  type(mesh_t), target :: mesh
  type(allocator_t), target :: omp_allocator
  class(allocator_t), pointer :: allocator
  type(omp_backend_t), target :: omp_backend
  class(base_backend_t), pointer :: backend
  type(solver_t) :: solver
  type(field_ptr_t), allocatable :: field_ptrs(:), host_fields(:)
  character(len=32) :: field_names(3)

  integer :: ierr, irank
  integer :: dims(3)
  logical :: allpass

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, irank, ierr)

  allpass = .true.

  mesh = mesh_t([4, 4, 4], [1, 1, 1], [1.0_dp, 1.0_dp, 1.0_dp], &
                ['periodic', 'periodic'], &
                ['periodic', 'periodic'], &
                ['periodic', 'periodic'])

  omp_allocator = allocator_t(mesh%get_dims(VERT), SZ)
  allocator => omp_allocator
  omp_backend = omp_backend_t(mesh, allocator)
  backend => omp_backend

  call init_solver_with_species(solver, backend, mesh, allocator)

  field_names = [character(len=32) :: 'u', 'phi_1', 'phi_2']
  call setup_field_arrays(solver, field_names, field_ptrs, host_fields)

  dims = mesh%get_dims(VERT)

  if (maxval(abs(host_fields(1)%ptr%data(1:dims(1), 1:dims(2), 1:dims(3)) &
                 - 11.0_dp)) > tol) then
    allpass = .false.
    if (irank == 0) print *, 'Velocity field was not copied correctly'
  end if

  if (maxval(abs(host_fields(2)%ptr%data(1:dims(1), 1:dims(2), 1:dims(3)) &
                 - 101.0_dp)) > tol) then
    allpass = .false.
    if (irank == 0) print *, 'phi_1 was not mapped to species 1'
  end if

  if (maxval(abs(host_fields(3)%ptr%data(1:dims(1), 1:dims(2), 1:dims(3)) &
                 - 202.0_dp)) > tol) then
    allpass = .false.
    if (irank == 0) print *, 'phi_2 was not mapped to species 2'
  end if

  call cleanup_field_arrays(solver, field_ptrs, host_fields)
  call release_solver_fields(solver)

  call MPI_Finalize(ierr)

  if (irank == 0) then
    if (allpass) then
      write (stderr, '(a)') 'Snapshot species field test passed.'
    else
      error stop 'Snapshot species field test failed.'
    end if
  end if

contains

  subroutine init_solver_with_species(solver, backend, mesh, allocator)
    type(solver_t), intent(inout) :: solver
    class(base_backend_t), pointer, intent(inout) :: backend
    type(mesh_t), target, intent(inout) :: mesh
    class(allocator_t), pointer, intent(inout) :: allocator

    integer :: is

    solver%backend => backend
    solver%mesh => mesh
    solver%host_allocator => allocator
    solver%nvars = 3 + nspecies
    solver%nspecies = nspecies
    solver%ngrid = product(mesh%get_global_dims(VERT))

    solver%u => allocator%get_block(DIR_X)
    solver%v => allocator%get_block(DIR_X)
    solver%w => allocator%get_block(DIR_X)
    call solver%u%set_data_loc(VERT)
    call solver%v%set_data_loc(VERT)
    call solver%w%set_data_loc(VERT)
    call solver%u%fill(11.0_dp)
    call solver%v%fill(22.0_dp)
    call solver%w%fill(33.0_dp)

    allocate (solver%species(nspecies))
    do is = 1, nspecies
      solver%species(is)%ptr => allocator%get_block(DIR_X)
      call solver%species(is)%ptr%set_data_loc(VERT)
      call solver%species(is)%ptr%fill(real(101*is, dp))
    end do
  end subroutine init_solver_with_species

  subroutine release_solver_fields(solver)
    type(solver_t), intent(inout) :: solver

    integer :: is

    if (associated(solver%backend)) then
      if (associated(solver%u)) &
        call solver%backend%allocator%release_block(solver%u)
      if (associated(solver%v)) &
        call solver%backend%allocator%release_block(solver%v)
      if (associated(solver%w)) &
        call solver%backend%allocator%release_block(solver%w)

      if (associated(solver%species)) then
        do is = 1, size(solver%species)
          if (associated(solver%species(is)%ptr)) then
            call solver%backend%allocator%release_block(solver%species(is)%ptr)
          end if
        end do
        deallocate (solver%species)
      end if
    end if

    nullify (solver%u)
    nullify (solver%v)
    nullify (solver%w)
    nullify (solver%species)
  end subroutine release_solver_fields

end program test_snapshot_species_fields
