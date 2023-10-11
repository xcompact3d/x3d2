module m_omp_kernels_dist
   use omp_lib

   use m_common, only: dp
   use m_omp_common, only: SZ

   implicit none

contains

   subroutine der_univ_dist_omp( &
      du, send_u_b, send_u_e, u, u_b, u_e, coeffs_b, coeffs_e, coeffs, n, &
      alfa, ffr, fbc &
      )
      implicit none

      ! Arguments
      real(dp), intent(out), dimension(:, :) :: du, send_u_b, send_u_e
      real(dp), intent(in), dimension(:, :) :: u, u_b, u_e
      real(dp), intent(in), dimension(:) :: ffr, fbc
      real(dp), intent(in), dimension(:, :) :: coeffs_b, coeffs_e
      real(dp), intent(in), dimension(:) :: coeffs
      real(dp), intent(in) :: alfa
      integer, intent(in) :: n

      ! Local variables
      integer :: i, j!, b
      integer :: jm2, jm1, jp1, jp2
      integer :: n_s, n_m, n_b, n_e !stencil, middle, begin, end

      real(dp) :: temp_du, c_m4, c_m3, c_m2, c_m1, c_j, c_p1, c_p2, c_p3, c_p4

      !i = threadIdx%x
      !b = blockIdx%x
      !nblock = size(u, dim=3)

      n_s = (size(coeffs)-1)/2
      n_m = size(coeffs)
      n_b = size(coeffs_b, dim=2)
      n_e = size(coeffs_e, dim=2)

      ! store bulk coeffs in the registers
      c_m4 = coeffs(1); c_m3 = coeffs(2); c_m2 = coeffs(3); c_m1 = coeffs(4)
      c_j = coeffs(5)
      c_p1 = coeffs(6); c_p2 = coeffs(7); c_p3 = coeffs(8); c_p4 = coeffs(9)

      !$omp simd
      do i = 1, SZ
      du(i, 1) = coeffs(1)*u_b(i, 1) + coeffs(2)*u_b(i, 2) &
                 + coeffs(3)*u_b(i, 3) + coeffs(4)*u_b(i, 4) &
                 + coeffs(5)*u(i, 1) &
                 + coeffs(6)*u(i, 2) + coeffs(7)*u(i, 3) &
                 + coeffs(8)*u(i, 4) + coeffs(9)*u(i, 5)
      du(i, 2) = coeffs(1)*u_b(i, 2) + coeffs(2)*u_b(i, 3) &
                 + coeffs(3)*u_b(i, 4) + coeffs(4)*u(i, 1) &
                 + coeffs(5)*u(i, 2) &
                 + coeffs(6)*u(i, 3) + coeffs(7)*u(i, 4) &
                 + coeffs(8)*u(i, 5) + coeffs(9)*u(i, 6)
      du(i, 3) = coeffs(1)*u_b(i, 3) + coeffs(2)*u_b(i, 4) &
                 + coeffs(3)*u(i, 1) + coeffs(4)*u(i, 2) &
                 + coeffs(5)*u(i, 3) &
                 + coeffs(6)*u(i, 4) + coeffs(7)*u(i, 5) &
                 + coeffs(8)*u(i, 6) + coeffs(9)*u(i, 7)
      du(i, 3) = ffr(3)*(du(i, 3) - alfa*du(i, 2))
      du(i, 4) = coeffs(1)*u_b(i, 4) + coeffs(2)*u(i, 1) &
                 + coeffs(3)*u(i, 2) + coeffs(4)*u(i, 3) &
                 + coeffs(5)*u(i, 4) &
                 + coeffs(6)*u(i, 5) + coeffs(7)*u(i, 6) &
                 + coeffs(8)*u(i, 7) + coeffs(9)*u(i, 8)
      du(i, 4) = ffr(4)*(du(i, 4) - alfa*du(i, 3))
      end do
      !$omp end simd

      do j = n_s+1, n-n_s
         !$omp simd
         do i = 1, SZ
         temp_du = c_m4*u(i, j-4) + c_m3*u(i, j-3) &
                 + c_m2*u(i, j-2) + c_m1*u(i, j-1) &
                 + c_j*u(i, j) &
                 + c_p1*u(i, j+1) + c_p2*u(i, j+2) &
                 + c_p3*u(i, j+3) + c_p4*u(i, j+4)
         du(i, j) = ffr(j)*(temp_du - alfa*du(i, j-1))
         end do
         !$omp end simd
      end do

      !$omp simd
      do i = 1, SZ
      j = n-3
      du(i, j) = coeffs(1)*u(i, j-4) + coeffs(2)*u(i, j-3) &
               + coeffs(3)*u(i, j-2) + coeffs(4)*u(i, j-1) &
               + coeffs(5)*u(i, j) &
               + coeffs(6)*u(i, j+1) + coeffs(7)*u(i, j+2) &
               + coeffs(8)*u(i, j+3) + coeffs(9)*u_e(i, 1)
      du(i, j) = ffr(j)*(du(i, j) - alfa*du(i, j-1))
      j = n-2
      du(i, j) = coeffs(1)*u(i, j-4) + coeffs(2)*u(i, j-3) &
                 + coeffs(3)*u(i, j-2) + coeffs(4)*u(i, j-1) &
                 + coeffs(5)*u(i, j) &
                 + coeffs(6)*u(i, j+1) + coeffs(7)*u(i, j+2) &
                 + coeffs(8)*u_e(i, 1) + coeffs(9)*u_e(i, 2)
      du(i, j) = ffr(j)*(du(i, j) - alfa*du(i, j-1))
      j = n-1
      du(i, j) = coeffs(1)*u(i, j-4) + coeffs(2)*u(i, j-3) &
                 + coeffs(3)*u(i, j-2) + coeffs(4)*u(i, j-1) &
                 + coeffs(5)*u(i, j) &
                 + coeffs(6)*u(i, j+1) + coeffs(7)*u_e(i, 1) &
                 + coeffs(8)*u_e(i, 2) + coeffs(9)*u_e(i, 3)
      du(i, j) = ffr(j)*(du(i, j) - alfa*du(i, j-1))
      j = n
      du(i, j) = coeffs(1)*u(i, j-4) + coeffs(2)*u(i, j-3) &
                 + coeffs(3)*u(i, j-2) + coeffs(4)*u(i, j-1) &
                 + coeffs(5)*u(i, j) &
                 + coeffs(6)*u_e(i, 1) + coeffs(7)*u_e(i, 2) &
                 + coeffs(8)*u_e(i, 3) + coeffs(9)*u_e(i, 4)
      du(i, j) = ffr(j)*(du(i, j) - alfa*du(i, j-1))
      end do
      !$omp end simd

      !$omp simd
      do i = 1, SZ
      send_u_e(i, 1) = du(i, n)
      end do
      !$omp end simd

      ! Backward pass of the hybrid algorithm
      do j = n - 2, 2, -1
         !$omp simd
         do i = 1, SZ
         du(i, j) = du(i, j) - fbc(j)*du(i, j + 1)
         end do
         !$omp end simd
      end do
      !$omp simd
      do i = 1, SZ
      du(i, 1) = ffr(1)*(du(i, 1) - fbc(1)*du(i, 2))
      send_u_b(i, 1) = du(i, 1)
      end do
      !$omp end simd

   end subroutine der_univ_dist_omp

   subroutine der_univ_subs_omp(du, recv_u_b, recv_u_e, n, alfa, &
                                dist_sa, dist_sc)
      implicit none

      ! Arguments
      real(dp), intent(out), dimension(:, :) :: du
      real(dp), intent(in), dimension(:, :) :: recv_u_b, recv_u_e
      real(dp), intent(in), dimension(:) :: dist_sa, dist_sc
      real(dp), intent(in) :: alfa
      integer, intent(in) :: n

      ! Local variables
      integer :: i, j!, b
      real(dp) :: ur, bl, recp, du_1, du_n

      !i = threadIdx%x
      !b = blockIdx%x

      bl = dist_sa(1)
      ur = dist_sc(n)
      recp = 1._dp/(1._dp - ur*bl)

      !$omp simd
      do i = 1, SZ
      !du(i, 1) = recp*(du(i, 1) - bl*recv_u_b(i, 1))
      !du(i, n) = recp*(du(i, n) - ur*recv_u_e(i, 1))
      du_1 = recp*(du(i, 1) - bl*recv_u_b(i, 1))
      du_n = recp*(du(i, n) - ur*recv_u_e(i, 1))
      end do
      !$omp end simd

      !$omp simd
      do i = 1, SZ
      du(i, 1) = du_1
      end do
      !$omp end simd
      do j = 2, n-1
         !$omp simd
         do i = 1, SZ
         du(i, j) = (du(i, j) - dist_sa(j)*du_1 - dist_sc(j)*du_n)
         end do
         !$omp end simd
      end do
      !$omp simd
      do i = 1, SZ
      du(i, n) = du_n
      end do
      !$omp end simd

   end subroutine der_univ_subs_omp

end module m_omp_kernels_dist
