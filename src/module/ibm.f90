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

  use m_adios2_io, only : adios2_reader_t, adios2_file_t, adios2_mode_read
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
    type(field_t), pointer :: ep1 => null()
  contains
    procedure :: body
  end type ibm_t

  interface ibm_t
    module procedure init
  end interface ibm_t

contains

  function init(backend, mesh, host_allocator, ibm_type) result(ibm)
    !! Initialize the basic IBM
    implicit none

    class(base_backend_t), target, intent(inout) :: backend
    type(mesh_t), target, intent(inout) :: mesh
    type(allocator_t), target, intent(inout) :: host_allocator
    integer, intent(in) :: ibm_type
    type(ibm_t) :: ibm

    integer :: dims(3)
    real(dp), allocatable :: field_data(:, :, :)
    class(field_t), pointer :: ep1
    type(adios2_reader_t) :: reader
    type(adios2_file_t) :: file
    integer(i8), dimension(3) :: start_dims, count_dims
    ! The variables below will be removed when the IBM preprocessor is ready
    integer :: i, j, k
    real(dp) :: coords(3), dx, dy, dz, dR, x0, y0, R0

    ibm%backend => backend
    ibm%mesh => mesh
    ibm%host_allocator => host_allocator

    ! Open the IBM file and read the iibm parameter
    if (.false.) then
    file = reader%open("filename", adios2_mode_read, MPI_COMM_WORLD)
    call reader%begin_step(file)
    call reader%read_data("iibm", ibm%iibm, file)
    else
       ! This will be removed when the IBM preprocessor is ready
       ibm%iibm = ibm_type
    end if

    ! Basic IBM only needs ep1 in the vertices
    if (ibm%iibm == iibm_basic) then

      ! Get a field on the host
      ep1 => ibm%host_allocator%get_block(DIR_C)

      ! Read the file ep1
      if (.false.) then
      dims = mesh%get_dims(VERT)
      start_dims = int(ibm%mesh%par%n_offset, i8)
      count_dims = int(dims, i8)
      call reader%read_data("ep1", field_data, file, start_dims, count_dims)
      ep1%data(1:dims(1), 1:dims(2), 1:dims(3)) = field_data
      else

        ! This will be removed when the IBM preprocessor is ready
        ! Compute analytical ep1
        x0 = 1.d0
        y0 = 3.d0
        R0 = 0.25d0
        dims = mesh%get_dims(VERT)
        do k = 1, dims(3)
          do j = 1, dims(2)
            do i = 1, dims(1)
              coords = mesh%get_coordinates(i, j, k)
              dx = coords(1)
              dy = coords(2)
              dz = coords(3)

              dR = sqrt((dx - x0)**2 + (dy - y0)**2)
              if (dR <= R0) then
                ep1%data(i, j, k) = 0.d0
              else
                ep1%data(i, j, k) = 1.d0
              end if

            end do
          end do
        end do

      end if

      ! Get a block on the device
      ibm%ep1 => ibm%backend%allocator%get_block(DIR_X)

      ! Move the local host array ep1 to ibm%ep1
      call ibm%backend%set_field_data(ibm%ep1, ep1%data)

      ! Free memory
      call ibm%host_allocator%release_block(ep1)
      deallocate (field_data)

    end if

    ! Close the IBM file
    if (.false.) then
    call reader%close(file)
    end if

  end function init

  subroutine body(self, u, v, w)
    !! Apply basic IBM before the pressure solver
    implicit none

    class(ibm_t) :: self
    class(field_t), intent(inout) :: u, v, w

    ! vel = vel * ep1
    call self%backend%vecmult(u, self%ep1)
    call self%backend%vecmult(v, self%ep1)
    call self%backend%vecmult(w, self%ep1)

  end subroutine body

end module m_ibm
