module m_base_poisson_cg

  use m_allocator, only: field_t
  use m_base_backend, only: base_backend_t
  use m_common, only: dp, &
                      RDR_X2Z, RDR_Y2Z, RDR_C2Z, &
                      RDR_Z2X, RDR_Z2Y, RDR_Z2C, &
                      DIR_X, DIR_Y, DIR_Z, DIR_C, &
                      CELL, VERT, &
                      BC_PERIODIC
  use m_mesh, only: mesh_t
  use m_tdsops, only: tdsops_t, dirps_t
  use m_vector_calculus, only: vector_calculus_t

  implicit none

  private

  type, public :: laplace_operator_t
    !! Operator that computes the Laplacian of a field.
    private
    type(dirps_t) :: xdirps, ydirps, zdirps
    type(vector_calculus_t) :: vector_calculus
  contains
    procedure :: apply => poissmult
  end type laplace_operator_t

  interface laplace_operator_t
    !! Public constructor for the laplace_operator_t type.
    procedure init_lapl
  end interface laplace_operator_t

  type, abstract, public :: poisson_precon_t
    !! Wrapper definition of the Poisson preconditioner.
    !! User code should use this class which will instantiate the
    !! backend-specific code (determined at compile time).
    private
  contains
    private
    ! Applies the preconditioner to compute b=Px, mostly for testing purposes.
    procedure(apply_precon), public, deferred :: apply
  end type poisson_precon_t

  abstract interface
    subroutine apply_precon(self, p, b, backend)
      ! Applies the preconditioner to compute b=Px, mostly for testing purposes.
      import poisson_precon_t
      import field_t
      import base_backend_t
      implicit none
      class(poisson_precon_t) :: self
      class(field_t), intent(in) :: p    ! Pressure solution
      class(field_t), intent(inout) :: b ! The evaluated matrix-vector product
      class(base_backend_t), intent(in) :: backend
    end subroutine apply_precon
  end interface

  type, abstract, public :: poisson_solver_t
    !! Base definition of the backend-specific implementation of the iterative
    !! Poisson solver.
  contains
    private
    ! Solves the Poisson problem
    procedure(solve), public, deferred :: solve
  end type poisson_solver_t

  abstract interface
    subroutine solve(self, p, f, backend)
      import poisson_solver_t
      import field_t
      import base_backend_t
      implicit none
      class(poisson_solver_t) :: self
      class(field_t), intent(inout) :: p
      class(field_t), intent(in) :: f
      class(base_backend_t), intent(in) :: backend
    end subroutine solve
  end interface

contains

  function init_lapl(backend, mesh) result(lapl)
    !! Public constructor for the laplace_operator_t type.
    class(base_backend_t), target, intent(in) :: backend
    type(mesh_t), intent(in) :: mesh
    type(laplace_operator_t) :: lapl

    ! TODO: read these from config
    character(len=*), parameter :: stagder_scheme = "compact6"
    character(len=*), parameter :: interpl_scheme = "classic"

    lapl%vector_calculus = vector_calculus_t(backend)

    lapl%xdirps%dir = DIR_X; lapl%ydirps%dir = DIR_Y; lapl%zdirps%dir = DIR_Z

    call allocate_lapl_tdsops(lapl%xdirps, backend, mesh, stagder_scheme, interpl_scheme)
    call allocate_lapl_tdsops(lapl%ydirps, backend, mesh, stagder_scheme, interpl_scheme)
    call allocate_lapl_tdsops(lapl%zdirps, backend, mesh, stagder_scheme, interpl_scheme)

  end function init_lapl

  subroutine allocate_lapl_tdsops(dirps, backend, mesh, stagder_scheme, interpl_scheme)
    type(dirps_t), intent(inout) :: dirps
    class(base_backend_t), intent(in) :: backend
    type(mesh_t), intent(in) :: mesh
    character(len=*), intent(in) :: stagder_scheme, interpl_scheme

    integer :: dir, bc_start, bc_end, n_vert, n_cell
    real(dp) :: d

    dir = dirps%dir
    bc_start = mesh%grid%BCs(dir, 1)
    bc_end = mesh%grid%BCs(dir, 2)
    d = mesh%geo%d(dir)

    n_vert = mesh%get_n(dir, VERT)
    n_cell = mesh%get_n(dir, CELL)

    call backend%alloc_tdsops(dirps%interpl_v2p, n_cell, d, 'interpolate', &
                              interpl_scheme, bc_start, bc_end, from_to='v2p')
    call backend%alloc_tdsops(dirps%interpl_p2v, n_vert, d, 'interpolate', &
                              interpl_scheme, bc_start, bc_end, from_to='p2v')
    call backend%alloc_tdsops(dirps%stagder_v2p, n_cell, d, 'stag-deriv', &
                              stagder_scheme, bc_start, bc_end, from_to='v2p')
    call backend%alloc_tdsops(dirps%stagder_p2v, n_vert, d, 'stag-deriv', &
                              stagder_scheme, bc_start, bc_end, from_to='p2v')
    
  end subroutine allocate_lapl_tdsops

  subroutine poissmult(self, f, p)
    !! Computes the action of the Laplace operator, i.e. `f = Ax` where `A` is
    !! the discrete Laplacian.
    class(laplace_operator_t) :: self
    class(field_t), intent(inout) :: f ! The output field
    class(field_t), intent(in) :: p    ! The input field

    class(field_t), pointer :: f_z, p_z
    
    integer :: reorder_op, reorder_op2z

    if (p%dir == DIR_Z .and. f%dir == DIR_Z) then
      call poissmult_dirz(self, f, p)
    else
      if (p%dir /= f%dir) then
        error stop "Currently orientations of P and F must match"
      end if
      if (f%dir == DIR_X) then
        reorder_op2z = RDR_X2Z
        reorder_op = RDR_Z2X
      else if (f%dir == DIR_Y) then
        reorder_op2z = RDR_Y2Z
        reorder_op = RDR_Z2Y
      else if (f%dir == DIR_C) then
        reorder_op2z = RDR_C2Z
        reorder_op = RDR_Z2C
      else
        error stop "Unsupported Poisson orientation"
      end if

      f_z => self%vector_calculus%backend%allocator%get_block(DIR_Z, CELL)
      p_z => self%vector_calculus%backend%allocator%get_block(DIR_Z, CELL)

      call self%vector_calculus%backend%reorder(p_z, p, reorder_op2z)
      
      call poissmult_dirz(self, f_z, p_z)

      call self%vector_calculus%backend%reorder(f, f_z, reorder_op)

      call self%vector_calculus%backend%allocator%release_block(f_z)
      call self%vector_calculus%backend%allocator%release_block(p_z)
    end if
    
  end subroutine poissmult

  subroutine poissmult_dirz(lapl, f, p)
    class(laplace_operator_t), intent(in) :: lapl
    class(field_t), intent(inout) :: f ! The output field
    class(field_t), intent(in) :: p    ! The input field

    if (p%dir /= DIR_Z .or. f%dir /= DIR_Z) then
      error stop "Currently orientations of P and F must be in Z"
    end if
    if (p%data_loc /= CELL .or. f%data_loc /= CELL) then
      error stop "The pressure Poisson equation must be evaluated at cell centres"
    end if

    ! call lapl%vector_calculus%divgrad(f, p, &
    !   lapl%xdirps%stagder_p2v, lapl%xdirps%interpl_p2v, &
    !   lapl%ydirps%stagder_p2v, lapl%ydirps%interpl_p2v, &
    !   lapl%zdirps%stagder_p2v, lapl%zdirps%interpl_p2v, &
    !   lapl%xdirps%stagder_v2p, lapl%xdirps%interpl_v2p, &
    !   lapl%ydirps%stagder_v2p, lapl%ydirps%interpl_v2p, &
    !   lapl%zdirps%stagder_v2p, lapl%zdirps%interpl_v2p)
    call lapl%vector_calculus%divgrad_stag(f, p, &
      lapl%xdirps%stagder_p2v, lapl%ydirps%stagder_p2v, lapl%zdirps%stagder_p2v, &
      lapl%xdirps%stagder_v2p, lapl%ydirps%stagder_v2p, lapl%zdirps%stagder_v2p)

  end subroutine poissmult_dirz
  

end module m_base_poisson_cg
