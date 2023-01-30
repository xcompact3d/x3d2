program main
  use mpi
  use m_vector3d_simd, only: vector3d_simd
  implicit none

  integer :: rankid, nranks, errcode
  type(vector3d_simd) :: a

  call mpi_init(errcode)
  call mpi_comm_rank(MPI_COMM_WORLD, rankid, errcode)
  call mpi_comm_size(MPI_COMM_WORLD, nranks, errcode)

  a = vector3d_simd('velocity', [8, 8, 4])

  call mpi_finalize(errcode)

end program main
