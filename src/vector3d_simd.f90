module m_slab_cpu

   use mpi

   use m_allocator, only: allocator_t, memblock_t
   use m_slab, only: slab
   use m_diffengine, only: diffengine_t

   implicit none

   type, extends(slab) :: slab_cpu_t
      private
      type(diffengine_t) :: diffeng, diffeng2
   contains
      procedure, public :: transport
      procedure, public :: transport_dir
      procedure, public :: div
   end type slab_cpu_t

   interface slab_cpu_t
      module procedure make_slab_cpu_t
   end interface slab_cpu_t

contains

   function make_slab_cpu(allocator, diffeng, diffeng2) result(self)
      type(allocator_t), pointer, intent(in) :: allocator
      type(diffengine_t), intent(in) :: diffeng, diffeng2
      type(slab_cpu_t) :: self
      integer :: SZ, n

      integer :: errcode

      call mpi_comm_rank(MPI_COMM_WORLD, self%rankid, errcode)
      call mpi_comm_size(MPI_COMM_WORLD, self%nranks, errcode)

      self%allocator => allocator
      self%diffeng = diffeng
      self%diffeng2 = diffeng2

      self%u1 => self%allocator%get_block()
      self%u2 => self%allocator%get_block()
      self%u3 => self%allocator%get_block()
   end function make_slab_cpu

   function transport(self)
      class(slab_cpu_t), intent(in) :: self
      class(vector3d), allocatable :: transport
      transport = slab_cpu_t(&
           & self%allocator, self%diffeng, self%diffeng2 &
           & )
      call self%transport_dir(1, transport%u(1))
      call self%transport_dir(2, transport%u(2))
      call self%transport_dir(3, transport%u(3))
   end function transport

   subroutine transport_dir(self, dim, rslt)
      class(slab_cpu_t), intent(in) :: self
      integer, intent(in) :: dim
      real, intent(out) :: rslt(:, :, :)
      integer :: i, j, k, SZ, n

      real, allocatable, dimension(:, :) :: du, d2u, usq, dusq

      SZ = size(rslt, 1)
      n = size(rslt, 2)

      allocate (du(SZ, n), d2u(SZ, n))
      allocate (usq(SZ, n), dusq(SZ, n))

      associate (u => self%u(dim), u_dir => self%u(1))
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

   function div(self) result(rslt)
      class(slab_cpu_t), intent(in) :: self
      class(vector3d), allocatable :: rslt
      allocate (slab_cpu_t :: rslt)
      rslt = slab_cpu_t( &
           & self%allocator, self%diffeng, self%diffeng2 &
           & )
      rslt%u(1) = self%u(1) + 1.
      rslt%u(2) = self%u(2) + 1.
      rslt%u(3) = self%u(3) + 1.
   end function div
end module m_slab_cpu
