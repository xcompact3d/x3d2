program test_precorrect_walls
  !! precorrect_walls adds the kept tangential correction on each Dirichlet
  !! face. Where a Dirichlet x-face meets a Dirichlet y-face (an inflow on a
  !! wall-bounded domain) the edge must get it once, from the y-face, as in
  !! Incompact3d's pre_correc. Adding it from both faces is a feedback the
  !! projection cannot undo, and the edge blows up.
  use m_common, only: dp, DIR_X, DIR_C, VERT
  use m_field, only: field_t
  use m_mesh, only: mesh_t
  use m_allocator, only: allocator_t
  use m_base_backend, only: base_backend_t
  use m_time_integrator, only: init
  use m_solver, only: solver_t
  use m_backend_runtime, only: backend_runtime_t
  use m_test_utils, only: initialise_mpi, finalise_test, global_all

  implicit none

  type(mesh_t), allocatable, target :: mesh
  type(backend_runtime_t), target :: runtime
  class(allocator_t), pointer :: allocator
  class(base_backend_t), pointer :: backend
  type(solver_t) :: solver

  real(dp), parameter :: dt = 1e-3_dp, gx = 1._dp, gy = 2._dp, gz = 3._dp
  character(len=9), parameter :: bc_x(2) = ['dirichlet', 'dirichlet'], &
                                 bc_y(2) = ['dirichlet', 'neumann  '], &
                                 bc_z(2) = ['periodic ', 'periodic ']

  real(dp), allocatable :: data(:, :, :), want(:, :, :)
  integer :: dims(3), padded(3), irank, nproc
  logical :: allpass

  call initialise_mpi(irank, nproc)
  allpass = .true.

  ! ny spans several y-blocks on both backends
  mesh = mesh_t([8, 37, 4], [1, 1, 1], [1._dp, 1._dp, 1._dp], &
                bc_x, bc_y, bc_z)
  call runtime%init(mesh)
  allocator => runtime%allocator
  backend => runtime%backend
  dims = mesh%get_dims(VERT)
  padded = allocator%get_padded_dims(DIR_C)
  allocate (data(padded(1), padded(2), padded(3)), want(padded(1), &
                                                        padded(2), padded(3)))

  solver%backend => backend
  solver%mesh => mesh
  solver%u => new_field(0._dp)
  solver%v => new_field(0._dp)
  solver%w => new_field(0._dp)
  solver%time_integrator = init(backend, allocator, 'AB3', 3)
  solver%time_integrator%stage_dt = dt
  solver%dp_last_dt = dt
  solver%dpdx_last => new_field(gx)
  solver%dpdy_last => new_field(gy)
  solver%dpdz_last => new_field(gz)

  call solver%precorrect_walls(solver%u, solver%v, solver%w)

  ! u: tangential only to the y-faces, Dirichlet only at the floor
  want = 0._dp
  want(1:dims(1), 1, 1:dims(3)) = gx
  call check('u', solver%u)

  ! v: tangential only to the x-faces, but not on the floor edge
  want = 0._dp
  want(1, 2:dims(2), 1:dims(3)) = gy
  want(dims(1), 2:dims(2), 1:dims(3)) = gy
  call check('v', solver%v)

  ! w: tangential to both; the floor edges get it once
  want = 0._dp
  want(1:dims(1), 1, 1:dims(3)) = gz
  want(1, 2:dims(2), 1:dims(3)) = gz
  want(dims(1), 2:dims(2), 1:dims(3)) = gz
  call check('w', solver%w)

  call global_all(allpass)
  call finalise_test(allpass, irank)

contains

  function new_field(value) result(f)
    real(dp), intent(in) :: value
    class(field_t), pointer :: f

    f => allocator%get_block(DIR_X)
    call f%set_data_loc(VERT)
    call f%fill(value)
  end function new_field

  subroutine check(label, f)
    character(*), intent(in) :: label
    class(field_t), intent(in) :: f

    real(dp) :: err

    call backend%get_field_data(data, f)
    err = maxval(abs(data(1:dims(1), 1:dims(2), 1:dims(3)) &
                     - want(1:dims(1), 1:dims(2), 1:dims(3))))
    if (err > 0._dp) then
      allpass = .false.
      print *, 'FAIL: ', label, ' max error ', err, ' floor-edge values ', &
        data(1, 1, 1), data(dims(1), 1, 1)
    end if
  end subroutine check

end program test_precorrect_walls
