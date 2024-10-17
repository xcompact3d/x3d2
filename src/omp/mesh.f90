module m_omp_mesh
    use mpi
    use decomp_2d, only: decomp_2d_init, DECOMP_2D_COMM_CART_X, xsize, xstart
    use m_mesh, only: mesh_t
    use m_common, only: dp, DIR_X, DIR_Y, DIR_Z, DIR_C, &
                        CELL, VERT, none, X_FACE, Y_FACE, Z_FACE, &
                        X_EDGE, Y_EDGE, Z_EDGE
    use m_field, only: field_t


    type, extends(mesh_t) :: omp_mesh_t
    contains
      procedure :: domain_decomposition => domain_decomposition_2decompfft
    end type omp_mesh_t

  interface omp_mesh_t
    module procedure init
  end interface omp_mesh_t

  private init

    contains

  function init(dims_global, nproc_dir, L_global, &
                     periodic_BC) result(mesh)
    integer, dimension(3), intent(in) :: dims_global
    integer, dimension(3), intent(in) :: nproc_dir ! Number of proc in each direction
    real(dp), dimension(3), intent(in) :: L_global
    logical, dimension(3), optional, intent(in) :: periodic_BC

    type(omp_mesh_t), allocatable :: mesh
    allocate(omp_mesh_t :: mesh) 

    mesh%mesh_t = mesh_t(dims_global, nproc_dir, L_global, &
                     periodic_BC)

    call domain_decomposition_2decompfft(mesh%grid, mesh%par)

  end function

end module m_omp_mesh