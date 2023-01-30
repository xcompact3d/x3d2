module m_vector3d_simd
  use m_vector3d, only: vector3d
  type, extends(vector3d) :: vector3d_simd
     real, allocatable :: u(:,:,:)
     real, allocatable :: v(:,:,:)
     real, allocatable :: w(:,:,:)
   contains
     procedure, public :: transport
     procedure, public :: div
  end type vector3d_simd

  interface vector3d_simd
     module procedure construct
  end interface vector3d_simd

  real, allocatable :: du(:, :, :)
  real, allocatable :: d2u(:, :, :)
  real, allocatable :: u2(:, :, :)
  real, allocatable :: du2(:, :, :)

contains

  subroutine ensure_work_arrays(dims)
    integer, intent(in) :: dims(3)

    if (.not. allocated(du)) then
       allocate(du(dims(1), dims(2), dims(3)))
       allocate(d2u(dims(1), dims(2), dims(3)))
       allocate(u2(dims(1), dims(2), dims(3)))
       allocate(du2(dims(1), dims(2), dims(3)))
    else
       if ( &
            & dims(1) > dim(du, 1) .or. &
            & dims(2) > dim(du, 2) .or. &
            & dims(3) > dim(du, 3) &
            & ) then
          du = reshape(du, dims)
          d2u = reshape(d2u, dims)
          u2 = reshape(u2, dims)
          du2 = reshape(du2, dims)
       end if
    end if
  end subroutine ensure_work_arrays

  type(vector3d_simd) function construct(name, dims)
    character(*), intent(in) :: name
    integer, intent(in) :: dims(3)

    call ensure_work_arrays(dims)

    call mpi_comm_rank(MPI_COMM_WORLD, construct%rankid, errcode)
    call mpi_comm_size(MPI_COMM_WORLD, construct%nranks, errcode)

    allocate(construct%u(dims(1), dims(2), dims(3)))
    allocate(construct%v(dims(1), dims(2), dims(3)))
    allocate(construct%w(dims(1), dims(2), dims(3)))

    construct%name = name
  end function construct

  subroutine transport(self, rslt)
    class(vector3d_simd), intent(in) :: self
    class(vector3d), intent(inout) :: rslt

    select type (rslt)
    type is (vector3d_simd)
       rslt%u = self%u + 1.
       rslt%v = self%v + 1.
       rslt%w = self%w + 1.
    class default
       error stop
    end select
  end subroutine transport

  subroutine div(self, rslt)
    class(vector3d_simd), intent(in) :: self
    class(vector3d), intent(inout) :: rslt

    select type (rslt)
    type is (vector3d_simd)
       rslt%u = self%u + 1.
       rslt%v = self%v + 1.
       rslt%w = self%w + 1.
    class default
       error stop
    end select
  end subroutine div
end module m_vector3d_simd
