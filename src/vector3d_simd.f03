module m_vector3d_simd
  use m_allocator, only: allocator_t, memblock_t
  use m_vector3d, only: vector3d
  use m_diffengine, only: diffengine

  type, extends(vector3d) :: vector3d_simd
     integer :: dims(3)
     type(allocator_t), pointer :: allocator
     type(memblock_t), pointer :: u1, u2, u3
   contains
     procedure, public :: transport
     procedure, public :: div
     procedure, public :: u => get_component_ptr
  end type vector3d_simd

  interface vector3d_simd
     module procedure construct
  end interface vector3d_simd

contains

  function make_vector3d_simd(allocator, diffeng, diffeng2) result(self)
    type(allocator_t), pointer, intent(in) :: allocator
    type(diffengine), intent(in) :: diffeng, diffeng2
    type(vector3d_simd) :: self

    call mpi_comm_rank(MPI_COMM_WORLD, self%rankid, errcode)
    call mpi_comm_size(MPI_COMM_WORLD, self%nranks, errcode)

    self%allocator => allocator
    self%diffeng = diffeng
    self%diffeng2 = diffeng2

    blockptr1 => self%allocator%get_block()
    blockptr2 => self%allocator%get_block()
    blockptr3 => self%allocator%get_block()
  end function make_vector3d_simd

  function get_component_ptr(self, i) result(ptr)
    class(vector3d_simd), intent(in) :: self
    integer, intent(in) :: i
    real, pointer :: ptr(:, :, :)

    select case(i)
    case (1)
       ptr => u1%data
    case (2)
       ptr => u2%data
    case (3)
       ptr => u3%data
    end select
  end function get_component_ptr

  pure function transport(self)
    class(vector3d_simd), intent(in) :: self
    class(vector3d), allocatable :: transport
    allocate(vector3d_simd::transport)
    transport = vector3d_simd(&
         & self%allocator, self%diffeng, self%diffeng2 &
         & )
    tranport%set_storage()
    select type (tranport)
    type is (vector3d_simd)
       call self%transport_dir(self%u(1), transport%u(1))
       call self%transport_dir(self%u(2), transport%u(2))
       call self%transport_dir(self%u(3), transport%u(3))
    class default
       error stop
    end select
  end function transport

  pure subroutine transport_dir(self, u, u_dir, rslt)
    class(vector3d_simd), intent(in) :: self
    real, intent(in) :: u(:, :, :)
    real, intent(in) :: u_dir(:, :, :)
    real, intent(out) :: rslt(:, :, :)

    du => self%allocator%get_block()
    d2u => self%allocator%get_block()
    u2 => self%allocator%get_block()
    du2 => self%allocator%get_block()

    layers: do k = 1, size(u, 3)
       call diffeng%diff(u(:, :, k), du%data)
       call diffeng2%diff(u(:, :, k), d2u%data)
       do j = 1, size(u, 2)
          !$omp simd
          do i = 1, size(u, 1)
             u2%data(i, j) = u(i, j, k) * u_dir(i, j, k)
          end do
          !$omp end simd
       end do
       call diffeng%diff(du2%data, u2%data)
       reshape(transport_dir, shape(u))
       do j = 1, size(u, 2)
          !$omp simd
          do i = 1, size(u, 1)
             rslt(i, j, k) = -0.5 * &
                  (u(i, j, k) * du%data(i, j) + du2%data(i, j)) &
                  & + xnu * d2u%data(i, j)
          !$omp end simd
          end do
       end do
    end do layers

    call self%allocator%release_block(du)
    call self%allocator%release_block(d2u)
    call self%allocator%release_block(u2)
    call self%allocator%release_block(du2)
  end function transport_dir

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
