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

    self%u1 => self%allocator%get_block()
    self%u2 => self%allocator%get_block()
    self%u3 => self%allocator%get_block()
  end function make_vector3d_simd

  function transport(self)
    class(vector3d_simd), intent(in) :: self
    class(vector3d), allocatable :: transport
    allocate(vector3d_simd::transport)
    transport = vector3d_simd( &
         & self%allocator, self%diffeng, self%diffeng2 &
         & )
    select type (transport)
    type is (vector3d_simd)
       call self%transport_dir(1, transport%u(1))
       call self%transport_dir(2, transport%u(2))
       call self%transport_dir(3, transport%u(3))
    class default
       error stop
    end select
  end function transport

  subroutine transport_dir(self, dim, rslt)
    class(vector3d_simd), intent(in) :: self
    integer, intent(in) :: dim
    real, intent(out) :: rslt(:, :, :)
    integer :: i, j, k, SZ, n

    real, allocatable, dimension(:, :) :: du, d2u, usq, dusq

    SZ = size(rslt, 1)
    n = size(rslt, 2)

    allocate(du(SZ, n), d2u(SZ, n))
    allocate(usq(SZ, n), dusq(SZ, n))

    associate(u => self%u(dim), u_dir => self%u(1))
      layers: do k = 1, size(u, 3)
         call self%diffeng%diff(u(:, :, k), du)
         call self%diffeng2%diff(u(:, :, k), d2u)
         do j = 1, size(u, 2)
            !$omp simd
            do i = 1, size(u, 1)
               usq(i, j) = u(i, j, k) * u_dir(i, j, k)
            end do
            !$omp end simd
         end do
         call self%diffeng%diff(dusq, usq)
         do j = 1, size(u, 2)
            !$omp simd
            do i = 1, size(u, 1)
               rslt(i, j, k) = -0.5 * &
                    (u(i, j, k) * du(i, j) + dusq(i, j)) &
                    & + self%xnu * d2u(i, j)
               !$omp end simd
            end do
         end do
      end do layers
    end associate
  end subroutine transport_dir
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
