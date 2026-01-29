module m_omp_poisson_fft
  !! FFT-based Poisson solver for OMP backend.
  !!
  !! Solves $\nabla^2 \phi = f$ using spectral methods with 2DECOMP&FFT library.
  !! Transforms to Fourier space, solves diagonal system in spectral space,
  !! then transforms back to physical space.
  !!
  !! **Algorithm:**
  !! 1. Forward FFT: physical $\rightarrow$ spectral space
  !! 2. Spectral solve: $\phi_k = f_k / k^2$ (with modifications for boundary conditions)
  !! 3. Backward FFT: spectral $\rightarrow$ physical space
  !!
  !! **Boundary conditions:**
  !! - (0,0,0): Periodic in all directions
  !! - (0,1,0): Dirichlet in Y, periodic in X/Z (uses symmetry transform)
  !!
  !! **Parallelisation:** MPI via 2DECOMP&FFT pencil decomposition
  !!
  !! **Limitation:** Does not support Y-direction grid stretching

  use decomp_2d_constants, only: PHYSICAL_IN_X
  use decomp_2d_fft, only: decomp_2d_fft_init, decomp_2d_fft_3d, &
                           decomp_2d_fft_get_size

  use m_common, only: dp, CELL
  use m_field, only: field_t
  use m_mesh, only: mesh_t
  use m_poisson_fft, only: poisson_fft_t
  use m_tdsops, only: dirps_t

  use m_omp_spectral, only: process_spectral_000, process_spectral_010

  implicit none

  type, extends(poisson_fft_t) :: omp_poisson_fft_t
      !! FFT based Poisson solver
    complex(dp), allocatable, dimension(:, :, :) :: c_x  !! Spectral space buffer (X-pencil oriented)
  contains
    procedure :: fft_forward => fft_forward_omp           !! Transform to spectral space
    procedure :: fft_backward => fft_backward_omp         !! Transform to physical space
    procedure :: fft_postprocess_000 => fft_postprocess_000_omp  !! Spectral solve for (0,0,0) BCs
    procedure :: fft_postprocess_010 => fft_postprocess_010_omp  !! Spectral solve for (0,1,0) BCs
    procedure :: enforce_periodicity_y => enforce_periodicity_y_omp  !! Symmetry transform for Y Dirichlet
    procedure :: undo_periodicity_y => undo_periodicity_y_omp        !! Inverse symmetry transform
  end type omp_poisson_fft_t

  interface omp_poisson_fft_t
    module procedure init
  end interface omp_poisson_fft_t

  private :: init

contains

  function init(mesh, xdirps, ydirps, zdirps, lowmem) result(poisson_fft)
    !! Initialise FFT-based Poisson solver.
    !!
    !! Sets up 2DECOMP&FFT library and allocates spectral space buffers.
    !! Computes wavenumbers and coefficients for spectral solve.
    !!
    !! **Error checking:** Fails if Y-direction grid stretching requested
    !! (not supported by FFT method).
    implicit none

    type(mesh_t), intent(in) :: mesh                 !! Mesh with grid spacing
    class(dirps_t), intent(in) :: xdirps, ydirps, zdirps  !! Spectral operators
    logical, optional, intent(in) :: lowmem          !! Low-memory flag (ignored for OMP)
    integer, dimension(3) :: istart, iend, isize     !! Local spectral dimensions
    integer :: dims(3)                               !! Global grid dimensions

    type(omp_poisson_fft_t) :: poisson_fft           !! Initialised solver

    if (mesh%par%is_root()) then
      print *, "Initialising 2decomp&fft"
    end if

    if (present(lowmem)) then
      print *, 'lowmem_fft has no impact in the OpenMP backend.'
    end if

    ! Get global cell dims
    dims = mesh%get_global_dims(CELL)

    ! Work out the spectral dimensions in the permuted state
    call decomp_2d_fft_init(PHYSICAL_IN_X, dims(1), dims(2), dims(3))
    call decomp_2d_fft_get_size(istart, iend, isize)
    ! Converts a start position into an offset
    istart(:) = istart(:) - 1

    call poisson_fft%base_init(mesh, xdirps, ydirps, zdirps, isize, istart)

    if (mesh%geo%stretched(2)) then
      error stop 'OpenMP backends FFT based Poisson solver does not support&
                  & stretching in y-direction yet!'
    end if

    allocate (poisson_fft%c_x(poisson_fft%nx_spec, poisson_fft%ny_spec, &
                              poisson_fft%nz_spec))

  end function init

  subroutine fft_forward_omp(self, f_in)
    !! Forward FFT: physical space to spectral space.
    !!
    !! Transforms input field from physical (real) to spectral (complex)
    !! representation using 2DECOMP&FFT. Result stored in `self%c_x`.
    implicit none

    class(omp_poisson_fft_t) :: self      !! Solver instance
    class(field_t), intent(in) :: f_in    !! Physical space field (RHS)

    call decomp_2d_fft_3d(f_in%data, self%c_x)

  end subroutine fft_forward_omp

  subroutine fft_backward_omp(self, f_out)
    !! Backward FFT: spectral space to physical space.
    !!
    !! Transforms spectral solution back to physical (real) space using
    !! inverse FFT. Reads from `self%c_x`, writes to output field.
    implicit none

    class(omp_poisson_fft_t) :: self         !! Solver instance
    class(field_t), intent(inout) :: f_out   !! Physical space solution

    call decomp_2d_fft_3d(self%c_x, f_out%data)

  end subroutine fft_backward_omp

  subroutine fft_postprocess_000_omp(self)
    !! Spectral solve for (0,0,0) boundary conditions.
    !!
    !! Solves Poisson equation in spectral space for fully periodic domain.
    !! Divides each Fourier mode by its corresponding $k^2$ eigenvalue.
    !!
    !! **Formula:** $\hat{\phi}_k = \hat{f}_k / (k_x^2 + k_y^2 + k_z^2)$
    implicit none

    class(omp_poisson_fft_t) :: self  !! Solver instance

    call process_spectral_000( &
      self%c_x, self%waves, self%nx_spec, self%ny_spec, self%nz_spec, &
      self%x_sp_st, self%y_sp_st, self%z_sp_st, &
      self%nx_glob, self%ny_glob, self%nz_glob, &
      self%ax, self%bx, self%ay, self%by, self%az, self%bz &
      )

  end subroutine fft_postprocess_000_omp

  subroutine fft_postprocess_010_omp(self)
    !! Spectral solve for (0,1,0) boundary conditions.
    !!
    !! Solves Poisson equation with Dirichlet BCs in Y-direction,
    !! periodic in X and Z. Uses modified wavenumbers accounting for
    !! symmetry transformation (sine series in Y).
    !!
    !! **Formula:** Modified $k_y$ for sine series representation
    implicit none

    class(omp_poisson_fft_t) :: self  !! Solver instance

    call process_spectral_010( &
      self%c_x, self%waves, self%nx_spec, self%ny_spec, self%nz_spec, &
      self%x_sp_st, self%y_sp_st, self%z_sp_st, &
      self%nx_glob, self%ny_glob, self%nz_glob, &
      self%ax, self%bx, self%ay, self%by, self%az, self%bz &
      )

  end subroutine fft_postprocess_010_omp

  subroutine enforce_periodicity_y_omp(self, f_out, f_in)
    !! Apply symmetry transform for Y Dirichlet boundary conditions.
    !!
    !! Converts physical field to symmetric/antisymmetric representation
    !! suitable for sine series FFT. Used before forward FFT when Y-direction
    !! has Dirichlet (non-periodic) BCs.
    !!
    !! **Transformation:** Maps domain to symmetric extension for sine basis.
    implicit none

    class(omp_poisson_fft_t) :: self       !! Solver instance
    class(field_t), intent(inout) :: f_out  !! Transformed field
    class(field_t), intent(in) :: f_in      !! Original field

    integer :: i, j, k

    !$omp parallel do
    do k = 1, self%nz_loc
      do j = 1, self%ny_glob/2
        do i = 1, self%nx_loc
          f_out%data(i, j, k) = f_in%data(i, 2*(j - 1) + 1, k)
        end do
      end do
      do j = self%ny_glob/2 + 1, self%ny_glob
        do i = 1, self%nx_loc
          f_out%data(i, j, k) = f_in%data(i, 2*self%ny_glob - 2*j + 2, k)
        end do
      end do
    end do
    !$omp end parallel do

  end subroutine enforce_periodicity_y_omp

  subroutine undo_periodicity_y_omp(self, f_out, f_in)
    !! Inverse symmetry transform for Y Dirichlet boundary conditions.
    !!
    !! Converts symmetric/antisymmetric representation back to physical
    !! field. Used after backward FFT when Y-direction has Dirichlet BCs.
    !!
    !! **Transformation:** Extracts physical domain from symmetric extension.
    implicit none

    class(omp_poisson_fft_t) :: self       !! Solver instance
    class(field_t), intent(inout) :: f_out  !! Physical field
    class(field_t), intent(in) :: f_in      !! Transformed field

    integer :: i, j, k

    !$omp parallel do
    do k = 1, self%nz_loc
      do i = 1, self%nx_loc
        do j = 1, self%ny_glob/2
          f_out%data(i, 2*j - 1, k) = f_in%data(i, j, k)
        end do
        do j = 1, self%ny_glob/2
          f_out%data(i, 2*j, k) = f_in%data(i, self%ny_glob - j + 1, k)
        end do
      end do
    end do
    !$omp end parallel do

  end subroutine undo_periodicity_y_omp

end module m_omp_poisson_fft
