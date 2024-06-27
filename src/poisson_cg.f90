module m_poisson_cg
  !! Module defining a Poisson solver based on the (preconditioned) Conjugate
  !! Gradient method.

  use m_mesh, only: mesh_t
  use m_common, only: RDR_X2Y, RDR_X2Z, DIR_Y, DIR_Z
  use m_allocator, only: allocator_t, field_t
  use m_base_backend, only: base_backend_t
  use m_tdsops, only: dirps_t

  implicit none

  private
  public :: init_lapl

  type, public :: laplace_operator_t
    !! Operator that computes the Laplacian of a field.

    ! Prevent default access to components of type.
    private
  contains
    procedure :: apply => poissmult
  end type laplace_operator_t

  type, abstract, public :: poisson_cg_t
    !! Conjugate Gradient based Poisson solver.
    !! Supports any decomposition that is also supported by the underlying
    !! finite difference schemes.

    ! Prevent default access to components of type.
    private

    type(laplace_operator_t), pointer :: lapl
  contains
  end type poisson_cg_t

  interface
    module function init_cg(mesh, backend) result(poisson_cg)
      class(mesh_t), intent(in) :: mesh
      class(poisson_cg_t), allocatable :: poisson_cg
      class(base_backend_t), pointer, intent(in) :: backend
    end function init_cg
  end interface

  interface poisson_cg_t
    !! Public constructor for the poisson_cg_t type.
    module procedure init_cg
  end interface poisson_cg_t

contains

  subroutine init_lapl(lapl)

    type(laplace_operator_t), pointer :: lapl

    print *, "Initialise LAPL"
    
    allocate(lapl)
    
  end subroutine init_lapl
  
  subroutine poissmult(self, f, p, backend)
    !! Computes the action of the Laplace operator, i.e. `f = Ax` where `A` is
    !! the discrete Laplacian.
    class(laplace_operator_t) :: self
    class(field_t), intent(inout) :: f ! The output field
    class(field_t), intent(in) :: p    ! The input field
    class(base_backend_t), intent(in) :: backend

    ! Compute d2pdx2
    call compute_der2nd(f, p, backend, backend%xdirps)
    
    ! Compute d2pdy2, d2pdz2 and accumulate
    call compute_and_acc_der2nd(f, p, backend, backend%ydirps, RDR_X2Y)
    call compute_and_acc_der2nd(f, p, backend, backend%zdirps, RDR_X2Z)

  end subroutine poissmult

  subroutine compute_and_acc_der2nd(f, p, backend, dirps, reorder_op)

    class(field_t), intent(inout) :: f
    class(field_t), intent(in) :: p
    class(base_backend_t), intent(in) :: backend
    class(dirps_t), intent(in) :: dirps
    integer, intent(in) :: reorder_op

    class(field_t), pointer :: p_i ! P in operation order
    class(field_t), pointer :: f_i ! F in operation order

    integer :: DIR

    if (reorder_op == RDR_X2Y) then
      DIR = DIR_Y
    else if (reorder_op == RDR_X2Z) then
      DIR = DIR_Z
    else
      error stop "Unsupported reordering operation"
    end if
    
    p_i => backend%allocator%get_block(DIR)
    f_i => backend%allocator%get_block(DIR)
    call backend%reorder(p_i, p, reorder_op)

    call compute_der2nd(f_i, p_i, backend, dirps)
    if (reorder_op == RDR_X2Y) then
      call backend%sum_yintox(f, f_i)
    else if (reorder_op == RDR_X2Z) then
      call backend%sum_zintox(f, f_i)
    else
      error stop "Unsupported reordering operation"
    end if

    call backend%allocator%release_block(p_i)
    call backend%allocator%release_block(f_i)
    
  end subroutine compute_and_acc_der2nd

  subroutine compute_der2nd(d2fdx2, f, backend, dirps)
    !! Computes the 2nd derivative of a field
    class(field_t), intent(inout) :: d2fdx2
    class(field_t), intent(in) :: f
    class(base_backend_t), intent(in) :: backend
    class(dirps_t), intent(in) :: dirps

    call backend%tds_solve(d2fdx2, f, dirps, dirps%der2nd)
    
  end subroutine compute_der2nd

end module m_poisson_cg
