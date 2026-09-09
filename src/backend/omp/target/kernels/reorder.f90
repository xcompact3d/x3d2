!!! src/backend/omp/target/kernels/reorder.f90
!!
!! Transposing reorder kernels for the OpenMP target offload backend.
!!
!! A directional field is stored as (SZ, n, n_groups), where the leading
!! dimension holds SZ consecutive points along one cartesian direction. Which
!! direction that is depends on the storage: as laid out in m_ordering, DIR_C,
!! DIR_Y and DIR_Z all run their leading dimension along x, and only DIR_X runs
!! it along y.
!!
!! A reorder that changes which direction the leading dimension runs along is
!! therefore a transpose, and a kernel that simply walks the input has either
!! its reads or its writes strided by SZ. The kernels here avoid that by
!! staging an SZ x SZ tile: a first worksharing loop fills the tile, reading
!! the input contiguously, and a second drains it transposed, writing the
!! output contiguously.
!!
!! The tile is declared private to the `teams` construct, so all the threads of
!! a team work on one copy of it, and the implicit barrier that ends the first
!! worksharing loop makes the fill complete before the drain reads it.
!!
!! `target`, `teams` and `distribute` are written as three separate directives
!! rather than as one combined construct, which matters here: a clause on a
!! combined construct is applied to every constituent construct that accepts
!! it, and all three accept `private`. That would leave the number of copies of
!! the tile up to the compiler - one per kernel launch, one per team, or one
!! per thread - and only one per team is correct. Splitting the directives leaves
!! `teams` as the only construct the clause can attach to. Whether that copy
!! then lands in the device's shared memory is left to the compiler; OpenMP 5.0
!! can ask for it explicitly with `allocate(omp_pteam_mem_alloc: tile)`, which
!! not every compiler this backend builds with accepts yet.
!!
!! Only the transposing reorders are here, i.e. those between DIR_X and any
!! other direction. C2Y, C2Z, Y2Z and their inverses keep the leading dimension
!! along x on both sides, so they are contiguous whichever way they are walked,
!! and are left to the generic index mapping in the backend.

module m_omptgt_kernels_reorder

  use m_common, only: dp
  use m_omp_common, only: SZ

  implicit none

  private
  public :: reorder_omptgt_c2x, reorder_omptgt_x2c, &
            reorder_omptgt_x2y, reorder_omptgt_y2x, &
            reorder_omptgt_x2z, reorder_omptgt_z2x

contains

  subroutine reorder_omptgt_c2x(u_x, u_c, nx, ny, nz)
    real(dp), dimension(:, :, :), intent(inout) :: u_x
    real(dp), dimension(:, :, :), intent(in) :: u_c
    integer, intent(in) :: nx, ny, nz

    real(dp) :: tile(SZ, SZ)
    integer :: i, j, b_x, b_y, k

    !$omp target has_device_addr(u_x, u_c)
    !$omp teams private(tile)
    !$omp distribute collapse(3)
    do k = 1, nz
      do b_y = 1, ny/SZ
        do b_x = 1, nx/SZ
          !$omp parallel
          !$omp do collapse(2)
          do j = 1, SZ
            do i = 1, SZ
              tile(i, j) = u_c(i + (b_x - 1)*SZ, j + (b_y - 1)*SZ, k)
            end do
          end do
          !$omp end do
          !$omp do collapse(2)
          do j = 1, SZ
            do i = 1, SZ
              u_x(i, j + (b_x - 1)*SZ, b_y + (ny/SZ)*(k - 1)) = tile(j, i)
            end do
          end do
          !$omp end do
          !$omp end parallel
        end do
      end do
    end do
    !$omp end distribute
    !$omp end teams
    !$omp end target

  end subroutine reorder_omptgt_c2x

  subroutine reorder_omptgt_x2c(u_c, u_x, nx, ny, nz)
    real(dp), dimension(:, :, :), intent(inout) :: u_c
    real(dp), dimension(:, :, :), intent(in) :: u_x
    integer, intent(in) :: nx, ny, nz

    real(dp) :: tile(SZ, SZ)
    integer :: i, j, b_x, b_y, k

    !$omp target has_device_addr(u_c, u_x)
    !$omp teams private(tile)
    !$omp distribute collapse(3)
    do k = 1, nz
      do b_y = 1, ny/SZ
        do b_x = 1, nx/SZ
          !$omp parallel
          !$omp do collapse(2)
          do j = 1, SZ
            do i = 1, SZ
              tile(i, j) = u_x(i, j + (b_x - 1)*SZ, b_y + (ny/SZ)*(k - 1))
            end do
          end do
          !$omp end do
          !$omp do collapse(2)
          do j = 1, SZ
            do i = 1, SZ
              u_c(i + (b_x - 1)*SZ, j + (b_y - 1)*SZ, k) = tile(j, i)
            end do
          end do
          !$omp end do
          !$omp end parallel
        end do
      end do
    end do
    !$omp end distribute
    !$omp end teams
    !$omp end target

  end subroutine reorder_omptgt_x2c

  subroutine reorder_omptgt_x2y(u_y, u_x, nx, ny, nz)
    real(dp), dimension(:, :, :), intent(inout) :: u_y
    real(dp), dimension(:, :, :), intent(in) :: u_x
    integer, intent(in) :: nx, ny, nz

    real(dp) :: tile(SZ, SZ)
    integer :: i, j, b_x, b_y, k

    !$omp target has_device_addr(u_y, u_x)
    !$omp teams private(tile)
    !$omp distribute collapse(3)
    do k = 1, nz
      do b_y = 1, ny/SZ
        do b_x = 1, nx/SZ
          !$omp parallel
          !$omp do collapse(2)
          do j = 1, SZ
            do i = 1, SZ
              tile(i, j) = u_x(i, j + (b_x - 1)*SZ, b_y + (ny/SZ)*(k - 1))
            end do
          end do
          !$omp end do
          !$omp do collapse(2)
          do j = 1, SZ
            do i = 1, SZ
              u_y(i, j + (b_y - 1)*SZ, b_x + (nx/SZ)*(k - 1)) = tile(j, i)
            end do
          end do
          !$omp end do
          !$omp end parallel
        end do
      end do
    end do
    !$omp end distribute
    !$omp end teams
    !$omp end target

  end subroutine reorder_omptgt_x2y

  subroutine reorder_omptgt_y2x(u_x, u_y, nx, ny, nz)
    real(dp), dimension(:, :, :), intent(inout) :: u_x
    real(dp), dimension(:, :, :), intent(in) :: u_y
    integer, intent(in) :: nx, ny, nz

    real(dp) :: tile(SZ, SZ)
    integer :: i, j, b_x, b_y, k

    !$omp target has_device_addr(u_x, u_y)
    !$omp teams private(tile)
    !$omp distribute collapse(3)
    do k = 1, nz
      do b_y = 1, ny/SZ
        do b_x = 1, nx/SZ
          !$omp parallel
          !$omp do collapse(2)
          do j = 1, SZ
            do i = 1, SZ
              tile(i, j) = u_y(i, j + (b_y - 1)*SZ, b_x + (nx/SZ)*(k - 1))
            end do
          end do
          !$omp end do
          !$omp do collapse(2)
          do j = 1, SZ
            do i = 1, SZ
              u_x(i, j + (b_x - 1)*SZ, b_y + (ny/SZ)*(k - 1)) = tile(j, i)
            end do
          end do
          !$omp end do
          !$omp end parallel
        end do
      end do
    end do
    !$omp end distribute
    !$omp end teams
    !$omp end target

  end subroutine reorder_omptgt_y2x

  subroutine reorder_omptgt_x2z(u_z, u_x, nx, ny, nz)
    real(dp), dimension(:, :, :), intent(inout) :: u_z
    real(dp), dimension(:, :, :), intent(in) :: u_x
    integer, intent(in) :: nx, ny, nz

    real(dp) :: tile(SZ, SZ)
    integer :: i, j, b_x, b_y, k

    !$omp target has_device_addr(u_z, u_x)
    !$omp teams private(tile)
    !$omp distribute collapse(3)
    do k = 1, nz
      do b_y = 1, ny/SZ
        do b_x = 1, nx/SZ
          !$omp parallel
          !$omp do collapse(2)
          do j = 1, SZ
            do i = 1, SZ
              tile(i, j) = u_x(i, j + (b_x - 1)*SZ, b_y + (ny/SZ)*(k - 1))
            end do
          end do
          !$omp end do
          !$omp do collapse(2)
          do j = 1, SZ
            do i = 1, SZ
              u_z(i, k, b_x + (nx/SZ)*(j + (b_y - 1)*SZ - 1)) = tile(j, i)
            end do
          end do
          !$omp end do
          !$omp end parallel
        end do
      end do
    end do
    !$omp end distribute
    !$omp end teams
    !$omp end target

  end subroutine reorder_omptgt_x2z

  subroutine reorder_omptgt_z2x(u_x, u_z, nx, ny, nz)
    real(dp), dimension(:, :, :), intent(inout) :: u_x
    real(dp), dimension(:, :, :), intent(in) :: u_z
    integer, intent(in) :: nx, ny, nz

    real(dp) :: tile(SZ, SZ)
    integer :: i, j, b_x, b_y, k

    !$omp target has_device_addr(u_x, u_z)
    !$omp teams private(tile)
    !$omp distribute collapse(3)
    do k = 1, nz
      do b_y = 1, ny/SZ
        do b_x = 1, nx/SZ
          !$omp parallel
          !$omp do collapse(2)
          do j = 1, SZ
            do i = 1, SZ
              tile(i, j) = u_z(i, k, b_x + (nx/SZ)*(j + (b_y - 1)*SZ - 1))
            end do
          end do
          !$omp end do
          !$omp do collapse(2)
          do j = 1, SZ
            do i = 1, SZ
              u_x(i, j + (b_x - 1)*SZ, b_y + (ny/SZ)*(k - 1)) = tile(j, i)
            end do
          end do
          !$omp end do
          !$omp end parallel
        end do
      end do
    end do
    !$omp end distribute
    !$omp end teams
    !$omp end target

  end subroutine reorder_omptgt_z2x

end module m_omptgt_kernels_reorder
