module m_cuda_kernels_dist
   use cudafor

   use m_common, only: dp

   implicit none

contains

   attributes(global) subroutine der_univ_dist( &
      du, send_u_s, send_u_e, u, u_s, u_e, coeffs_s, coeffs_e, coeffs, n, &
      ffr, fbc, faf &
      )
      implicit none

      ! Arguments
      real(dp), device, intent(out), dimension(:, :, :) :: du, send_u_s, &
                                                           send_u_e
      real(dp), device, intent(in), dimension(:, :, :) :: u, u_s, u_e
      real(dp), device, intent(in), dimension(:, :) :: coeffs_s, coeffs_e
      real(dp), device, intent(in), dimension(:) :: coeffs
      integer, value, intent(in) :: n
      real(dp), device, intent(in), dimension(:) :: ffr, fbc, faf

      ! Local variables
      integer :: i, j, b, k, lj
      integer :: jm2, jm1, jp1, jp2

      real(dp) :: c_m4, c_m3, c_m2, c_m1, c_j, c_p1, c_p2, c_p3, c_p4, &
                  temp_du, alpha, last_r

      i = threadIdx%x
      b = blockIdx%x

      ! store bulk coeffs in the registers
      c_m4 = coeffs(1); c_m3 = coeffs(2); c_m2 = coeffs(3); c_m1 = coeffs(4)
      c_j = coeffs(5)
      c_p1 = coeffs(6); c_p2 = coeffs(7); c_p3 = coeffs(8); c_p4 = coeffs(9)
      last_r = ffr(1)

      du(i, 1, b) = coeffs_s(1, 1)*u_s(i, 1, b) &
                    + coeffs_s(2, 1)*u_s(i, 2, b) &
                    + coeffs_s(3, 1)*u_s(i, 3, b) &
                    + coeffs_s(4, 1)*u_s(i, 4, b) &
                    + coeffs_s(5, 1)*u(i, 1, b) &
                    + coeffs_s(6, 1)*u(i, 2, b) &
                    + coeffs_s(7, 1)*u(i, 3, b) &
                    + coeffs_s(8, 1)*u(i, 4, b) &
                    + coeffs_s(9, 1)*u(i, 5, b)
      du(i, 1, b) = du(i, 1, b)*faf(1)
      du(i, 2, b) = coeffs_s(1, 2)*u_s(i, 2, b) &
                    + coeffs_s(2, 2)*u_s(i, 3, b) &
                    + coeffs_s(3, 2)*u_s(i, 4, b) &
                    + coeffs_s(4, 2)*u(i, 1, b) &
                    + coeffs_s(5, 2)*u(i, 2, b) &
                    + coeffs_s(6, 2)*u(i, 3, b) &
                    + coeffs_s(7, 2)*u(i, 4, b) &
                    + coeffs_s(8, 2)*u(i, 5, b) &
                    + coeffs_s(9, 2)*u(i, 6, b)
      du(i, 2, b) = du(i, 2, b)*faf(2)
      du(i, 3, b) = coeffs_s(1, 3)*u_s(i, 3, b) &
                    + coeffs_s(2, 3)*u_s(i, 4, b) &
                    + coeffs_s(3, 3)*u(i, 1, b) &
                    + coeffs_s(4, 3)*u(i, 2, b) &
                    + coeffs_s(5, 3)*u(i, 3, b) &
                    + coeffs_s(6, 3)*u(i, 4, b) &
                    + coeffs_s(7, 3)*u(i, 5, b) &
                    + coeffs_s(8, 3)*u(i, 6, b) &
                    + coeffs_s(9, 3)*u(i, 7, b)
      du(i, 3, b) = ffr(3)*(du(i, 3, b) - faf(3)*du(i, 2, b))
      du(i, 4, b) = coeffs_s(1, 4)*u_s(i, 4, b) &
                    + coeffs_s(2, 4)*u(i, 1, b) &
                    + coeffs_s(3, 4)*u(i, 2, b) &
                    + coeffs_s(4, 4)*u(i, 3, b) &
                    + coeffs_s(5, 4)*u(i, 4, b) &
                    + coeffs_s(6, 4)*u(i, 5, b) &
                    + coeffs_s(7, 4)*u(i, 6, b) &
                    + coeffs_s(8, 4)*u(i, 7, b) &
                    + coeffs_s(9, 4)*u(i, 8, b)
      du(i, 4, b) = ffr(4)*(du(i, 4, b) - faf(3)*du(i, 3, b))

      alpha = faf(5)

      do j = 5, n - 4
         temp_du = c_m4*u(i, j - 4, b) + c_m3*u(i, j - 3, b) &
                   + c_m2*u(i, j - 2, b) + c_m1*u(i, j - 1, b) &
                   + c_j*u(i, j, b) &
                   + c_p1*u(i, j + 1, b) + c_p2*u(i, j + 2, b) &
                   + c_p3*u(i, j + 3, b) + c_p4*u(i, j + 4, b)
         du(i, j, b) = ffr(j)*(temp_du - alpha*du(i, j - 1, b))
      end do

      j = n - 3
      du(i, j, b) = coeffs_e(1, 1)*u(i, j - 4, b) &
                    + coeffs_e(2, 1)*u(i, j - 3, b) &
                    + coeffs_e(3, 1)*u(i, j - 2, b) &
                    + coeffs_e(4, 1)*u(i, j - 1, b) &
                    + coeffs_e(5, 1)*u(i, j, b) &
                    + coeffs_e(6, 1)*u(i, j + 1, b) &
                    + coeffs_e(7, 1)*u(i, j + 2, b) &
                    + coeffs_e(8, 1)*u(i, j + 3, b) &
                    + coeffs_e(9, 1)*u_e(i, 1, b)
      du(i, j, b) = ffr(j)*(du(i, j, b) - faf(j)*du(i, j - 1, b))
      j = n - 2
      du(i, j, b) = coeffs_e(1, 2)*u(i, j - 4, b) &
                    + coeffs_e(2, 2)*u(i, j - 3, b) &
                    + coeffs_e(3, 2)*u(i, j - 2, b) &
                    + coeffs_e(4, 2)*u(i, j - 1, b) &
                    + coeffs_e(5, 2)*u(i, j, b) &
                    + coeffs_e(6, 2)*u(i, j + 1, b) &
                    + coeffs_e(7, 2)*u(i, j + 2, b) &
                    + coeffs_e(8, 2)*u_e(i, 1, b) &
                    + coeffs_e(9, 2)*u_e(i, 2, b)
      du(i, j, b) = ffr(j)*(du(i, j, b) - faf(j)*du(i, j - 1, b))
      j = n - 1
      du(i, j, b) = coeffs_e(1, 3)*u(i, j - 4, b) &
                    + coeffs_e(2, 3)*u(i, j - 3, b) &
                    + coeffs_e(3, 3)*u(i, j - 2, b) &
                    + coeffs_e(4, 3)*u(i, j - 1, b) &
                    + coeffs_e(5, 3)*u(i, j, b) &
                    + coeffs_e(6, 3)*u(i, j + 1, b) &
                    + coeffs_e(7, 3)*u_e(i, 1, b) &
                    + coeffs_e(8, 3)*u_e(i, 2, b) &
                    + coeffs_e(9, 3)*u_e(i, 3, b)
      du(i, j, b) = ffr(j)*(du(i, j, b) - faf(j)*du(i, j - 1, b))
      j = n
      du(i, j, b) = coeffs_e(1, 4)*u(i, j - 4, b) &
                    + coeffs_e(2, 4)*u(i, j - 3, b) &
                    + coeffs_e(3, 4)*u(i, j - 2, b) &
                    + coeffs_e(4, 4)*u(i, j - 1, b) &
                    + coeffs_e(5, 4)*u(i, j, b) &
                    + coeffs_e(6, 4)*u_e(i, 1, b) &
                    + coeffs_e(7, 4)*u_e(i, 2, b) &
                    + coeffs_e(8, 4)*u_e(i, 3, b) &
                    + coeffs_e(9, 4)*u_e(i, 4, b)
      du(i, j, b) = ffr(j)*(du(i, j, b) - faf(j)*du(i, j - 1, b))

      send_u_e(i, 1, b) = du(i, n, b)

      ! Backward pass of the hybrid algorithm
      do j = n - 2, 2, -1
         du(i, j, b) = du(i, j, b) - fbc(j)*du(i, j + 1, b)
      end do
      du(i, 1, b) = last_r*(du(i, 1, b) - fbc(1)*du(i, 2, b))
      send_u_s(i, 1, b) = du(i, 1, b)

   end subroutine der_univ_dist

   attributes(global) subroutine der_univ_subs(du, recv_u_s, recv_u_e, &
                                               n, dist_sa, dist_sc)
      implicit none

      ! Arguments
      real(dp), device, intent(out), dimension(:, :, :) :: du
      real(dp), device, intent(in), dimension(:, :, :) :: recv_u_s, recv_u_e
      real(dp), device, intent(in), dimension(:) :: dist_sa, dist_sc
      integer, value, intent(in) :: n

      ! Local variables
      integer :: i, j, b
      real(dp) :: ur, bl, recp, du_s, du_e

      i = threadIdx%x
      b = blockIdx%x

      ! A small trick we do here is valid for symmetric Toeplitz matrices.
      ! In our case our matrices satisfy this criteria in the (5:n-4) region
      ! and as long as a rank has around at least 20 entries the assumptions
      ! we make here are perfectly valid.

      ! bl is the bottom left entry in the 2x2 matrix
      ! ur is the upper right entry in the 2x2 matrix

      ! Start
      ! At the start we have the 'bl', and assume 'ur'
      bl = dist_sa(1)
      ur = dist_sa(1)
      recp = 1._dp/(1._dp - ur*bl)
      du_s = recp*(du(i, 1, b) - bl*recv_u_s(i, 1, b))

      ! End
      ! At the end we have the 'ur', and assume 'bl'
      bl = dist_sc(n)
      ur = dist_sc(n)
      recp = 1._dp/(1._dp - ur*bl)
      du_e = recp*(du(i, n, b) - ur*recv_u_e(i, 1, b))

      du(i, 1, b) = du_s
      do j = 2, n - 1
         du(i, j, b) = (du(i, j, b) - dist_sa(j)*du_s - dist_sc(j)*du_e)
      end do
      du(i, n, b) = du_e

   end subroutine der_univ_subs

end module m_cuda_kernels_dist
