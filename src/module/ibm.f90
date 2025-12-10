module m_ibm
!! This module implements the IBM capabilities.
!!
!! When iibm = 0, the IBM object is never used.
!!
!! When iibm = 1, the basic IBM capability is used.
!! It only requires ep1, a 3D field, as input.
!! This field should be one (zero) in the fluid (solid)
!! domain.
  use iso_fortran_env, only: stderr => error_unit
  use mpi

  use m_io_session, only: reader_session_t
  use m_allocator, only: allocator_t, field_t
  use m_base_backend, only: base_backend_t
  use m_common, only: dp, i8, pi, DIR_X, DIR_C, VERT
  use m_field, only: field_t
  use m_mesh, only: mesh_t

  implicit none

  private
  public :: ibm_t, iibm_basic

  integer, parameter :: iibm_basic = 1

  type :: ibm_t
    class(base_backend_t), pointer :: backend => null()
    class(mesh_t), pointer :: mesh => null()
    type(allocator_t), pointer :: host_allocator => null()
    integer :: iibm = 0
    class(field_t), pointer :: ep1 => null()
  contains
    procedure :: body
  end type ibm_t

  interface ibm_t
    module procedure init
  end interface ibm_t

contains

  function init(backend, mesh, host_allocator) result(ibm)
    !! Initialize the basic IBM
    implicit none

    class(base_backend_t), target, intent(inout) :: backend
    type(mesh_t), target, intent(inout) :: mesh
    type(allocator_t), target, intent(inout) :: host_allocator
    type(ibm_t) :: ibm

    integer :: i, j, k
    integer :: dims(3)
    real(dp), allocatable :: field_data(:, :, :)
    class(field_t), pointer :: ep1
    type(reader_session_t) :: reader_session
    character(len=*), parameter :: ibm_file = "ibm.bp"
    integer(i8) :: start_dims(3), count_dims(3), iibm_i8

    ibm%backend => backend
    ibm%mesh => mesh
    ibm%host_allocator => host_allocator

    ! Open a session to read the IBM configuration file
    call reader_session%open(ibm_file, MPI_COMM_WORLD)

    ! Read the iibm parameter
    call reader_session%read_data("iibm", iibm_i8)
    ibm%iibm = int(iibm_i8, kind=4)

    ! Basic IBM only needs ep1 on the vertices
    if (ibm%iibm == iibm_basic) then

      ! Read the vertex mask ep1
      !
      ! The mask was written in python in C order
      ! start_dims and count_dims are thus reversed
      ! The resulting output is in reversed order
      dims = mesh%get_dims(VERT)
      start_dims = int(ibm%mesh%par%n_offset(3:1:-1), i8)
      count_dims = int(dims(3:1:-1), i8)

      ! Allocate field_data with the expected (reversed) dimensions from the read
      allocate (field_data(count_dims(1), count_dims(2), count_dims(3)))
      call reader_session%read_data("ep1", field_data, start_dims, count_dims)

      ! Get and fill a block on the host
      ! The order of the data is corrected in the loop below
      ep1 => ibm%host_allocator%get_block(DIR_C)
      call ep1%fill(1.0_dp)
      do i = 1, dims(1)
        do j = 1, dims(2)
          do k = 1, dims(3)
            ep1%data(i, j, k) = field_data(k, j, i)
          end do
        end do
      end do

      ! Get a block on the device and move the data
      ibm%ep1 => ibm%backend%allocator%get_block(DIR_X)
      call ibm%backend%set_field_data(ibm%ep1, ep1%data)

      ! Free memory
      call ibm%host_allocator%release_block(ep1)
      deallocate (field_data)

    else
      ! IBM disabled or unsupported type
      ibm%ep1 => null()
    end if

    ! Closing the session handles all file and reader finalisation
    call reader_session%close()

  end function init

  subroutine body(self, u, v, w)
    !! Apply basic IBM before the pressure solver
    implicit none

    class(ibm_t) :: self
    class(field_t), intent(inout) :: u, v, w

    if (self%iibm == iibm_basic) then

      ! vel = vel * ep1
      !
      ! FIXME : currently velocity is zero in the solid.
      !         It should be dt * grad(p^n).
      !         After reconstruction, it should be
      !   dt * grad(p^n) - dt * grad(p^n+1)
      !         Currently grad(p^n) is not available
      call self%backend%vecmult(u, self%ep1)
      call self%backend%vecmult(v, self%ep1)
      call self%backend%vecmult(w, self%ep1)

    end if

  end subroutine body

end module m_ibm
