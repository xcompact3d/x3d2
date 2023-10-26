module m_derparams
   use m_common, only: dp

   implicit none

contains

   subroutine der_1_vv()
      implicit none
   end subroutine der_1_vv

   subroutine der_2_vv(coeffs, coeffs_s, coeffs_e, &
                       dist_fr, dist_bc, dist_af, dist_sa, dist_sc, &
                       n_halo, dx2, n, bcond)
      implicit none

      real(dp), allocatable, dimension(:), intent(out) :: coeffs, &
         dist_fr, dist_bc, dist_af, dist_sa, dist_sc
      real(dp), allocatable, dimension(:,:), intent(out) :: coeffs_s, coeffs_e
      integer, intent(out) :: n_halo
      real(dp), intent(in) :: dx2
      integer, intent(in) :: n
      character(len=*), intent(in) :: bcond

      real(dp), allocatable :: dist_b(:)
      real(dp) :: alfa, asi, bsi, csi, dsi
      integer :: i, n_stencil

      allocate(dist_sa(n), dist_sc(n), dist_b(n))

      alfa = 2._dp/11._dp
      asi = 12._dp/11._dp/dx2
      bsi = 3._dp/44._dp/dx2
      csi = 0._dp
      dsi = 0._dp

      n_halo = 4
      n_stencil = 2*n_halo+1

      coeffs = [dsi, csi, bsi, asi, &
               -2._dp*(asi+bsi+csi+dsi), &
               asi, bsi, csi, dsi]

      select case (bcond)
      case ('periodic')
         dist_sa(:) = alfa; dist_sc(:) = alfa; dist_b(:) = 1._dp
         allocate(coeffs_s(n_stencil, n_halo))
         allocate(coeffs_e(n_stencil, n_halo))
         do i = 1, n_halo
            coeffs_s(:, i) = coeffs(:)
            coeffs_e(:, i) = coeffs(:)
         end do
      case default
         print*, 'Boundary condition is not recognized :', bcond
      end select

      call process_dist(dist_fr, dist_bc, dist_af, dist_sa, dist_sc, dist_b, n)

   end subroutine der_2_vv

   subroutine process_dist(dist_fr, dist_bc, dist_af, &
                           dist_sa, dist_sc, dist_b, n)
      implicit none

      real(dp), allocatable, dimension(:), intent(out) :: dist_fr, dist_bc, &
                                                          dist_af
      real(dp), dimension(:), intent(inout) :: dist_sa, dist_sc, dist_b
      integer, intent(in) :: n

      integer :: i, nrank, nproc, m

      m = n
      nrank = 0; nproc = 1

      ! forward factors, backward factors, and auxiliary factor
      allocate(dist_fr(n), dist_bc(n), dist_af(n))

      do nrank = 0, nproc-1

         ! start the hybrid algorithm
         do i = 1+m*nrank, 2+m*nrank
            dist_sa(i) = dist_sa(i)/dist_b(i)
            dist_sc(i) = dist_sc(i)/dist_b(i)
            dist_bc(i) = dist_sc(i)
            dist_af(i) = 1._dp/dist_b(i)
         end do
         do i = 3+m*nrank, m+m*nrank
            dist_fr(i) = 1.d0/(dist_b(i)-dist_sa(i)*dist_sc(i-1))
            dist_af(i) = dist_sa(i)
            dist_sa(i) = -dist_fr(i)*dist_sa(i)*dist_sa(i-1)
            dist_sc(i) = dist_fr(i)*dist_sc(i)
            !dist_bc(i) = dist_sc(i)
         end do
         do i = m-2+m*nrank, 2+m*nrank, -1
            dist_sa(i) = dist_sa(i)-dist_sc(i)*dist_sa(i+1)
            dist_bc(i) = dist_sc(i)
            dist_sc(i) = -dist_sc(i)*dist_sc(i+1)
         end do
         ! dist_fr(1) is never used, so store this extra factor instead.
         dist_fr(1+m*nrank) = 1.d0/(1.d0-dist_sc(1+m*nrank)*dist_sa(2+m*nrank))

         dist_sa(1+m*nrank) = dist_fr(1+m*nrank)*dist_sa(1+m*nrank)
         dist_sc(1+m*nrank) = -dist_fr(1+m*nrank)*dist_sc(1+m*nrank) &
                               *dist_sc(2+m*nrank)
         !dist_bc(1+m*nrank) = dist_sc(1+m*nrank)
      end do

   end subroutine process_dist

end module m_derparams

