module m_cuda_kernels_dist
   use cudafor

   use m_common, only: dp

   implicit none

contains

   attributes(global) subroutine der_univ_dist( &
      du, send_u_b, send_u_e, u, u_b, u_e, coeffs_b, coeffs_e, coeffs, n, &
      ffr, fbc, faf &
      )
      implicit none

      ! Arguments
      real(dp), device, intent(out), dimension(:, :, :) :: du, send_u_b, &
                                                           send_u_e
      real(dp), device, intent(in), dimension(:, :, :) :: u, u_b, u_e
      real(dp), device, intent(in), dimension(:, :) :: coeffs_b, coeffs_e
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

      du(i, 1, b) = coeffs(1)*u_b(i, 1, b) + coeffs(2)*u_b(i, 2, b) &
                    + coeffs(3)*u_b(i, 3, b) + coeffs(4)*u_b(i, 4, b) &
                    + coeffs(5)*u(i, 1, b) &
                    + coeffs(6)*u(i, 2, b) + coeffs(7)*u(i, 3, b) &
                    + coeffs(8)*u(i, 4, b) + coeffs(9)*u(i, 5, b)
      du(i, 1, b) = du(i, 1, b)*faf(1)
      du(i, 2, b) = coeffs(1)*u_b(i, 2, b) + coeffs(2)*u_b(i, 3, b) &
                    + coeffs(3)*u_b(i, 4, b) + coeffs(4)*u(i, 1, b) &
                    + coeffs(5)*u(i, 2, b) &
                    + coeffs(6)*u(i, 3, b) + coeffs(7)*u(i, 4, b) &
                    + coeffs(8)*u(i, 5, b) + coeffs(9)*u(i, 6, b)
      du(i, 2, b) = du(i, 2, b)*faf(2)
      du(i, 3, b) = coeffs(1)*u_b(i, 3, b) + coeffs(2)*u_b(i, 4, b) &
                    + coeffs(3)*u(i, 1, b) + coeffs(4)*u(i, 2, b) &
                    + coeffs(5)*u(i, 3, b) &
                    + coeffs(6)*u(i, 4, b) + coeffs(7)*u(i, 5, b) &
                    + coeffs(8)*u(i, 6, b) + coeffs(9)*u(i, 7, b)
      du(i, 3, b) = ffr(3)*(du(i, 3, b) - faf(3)*du(i, 2, b))
      du(i, 4, b) = coeffs(1)*u_b(i, 4, b) + coeffs(2)*u(i, 1, b) &
                    + coeffs(3)*u(i, 2, b) + coeffs(4)*u(i, 3, b) &
                    + coeffs(5)*u(i, 4, b) &
                    + coeffs(6)*u(i, 5, b) + coeffs(7)*u(i, 6, b) &
                    + coeffs(8)*u(i, 7, b) + coeffs(9)*u(i, 8, b)
      du(i, 4, b) = ffr(4)*(du(i, 4, b) - faf(3)*du(i, 3, b))

      alpha = faf(5)

      do j = 5, n-4
         temp_du = c_m4*u(i, j-4, b) + c_m3*u(i, j-3, b) &
                 + c_m2*u(i, j-2, b) + c_m1*u(i, j-1, b) &
                 + c_j*u(i, j, b) &
                 + c_p1*u(i, j+1, b) + c_p2*u(i, j+2, b) &
                 + c_p3*u(i, j+3, b) + c_p4*u(i, j+4, b)
         du(i, j, b) = ffr(j)*(temp_du - alpha*du(i, j-1, b))
      end do

      j = n-3
      du(i, j, b) = coeffs(1)*u(i, j-4, b) + coeffs(2)*u(i, j-3, b) &
                  + coeffs(3)*u(i, j-2, b) + coeffs(4)*u(i, j-1, b) &
                  + coeffs(5)*u(i, j, b) &
                  + coeffs(6)*u(i, j+1, b) + coeffs(7)*u(i, j+2, b) &
                  + coeffs(8)*u(i, j+3, b) + coeffs(9)*u_e(i, 1, b)
      du(i, j, b) = ffr(j)*(du(i, j, b) - faf(j)*du(i, j-1, b))
      j = n-2
      du(i, j, b) = coeffs(1)*u(i, j-4, b) + coeffs(2)*u(i, j-3, b) &
                    + coeffs(3)*u(i, j-2, b) + coeffs(4)*u(i, j-1, b) &
                    + coeffs(5)*u(i, j, b) &
                    + coeffs(6)*u(i, j+1, b) + coeffs(7)*u(i, j+2, b) &
                    + coeffs(8)*u_e(i, 1, b) + coeffs(9)*u_e(i, 2, b)
      du(i, j, b) = ffr(j)*(du(i, j, b) - faf(j)*du(i, j-1, b))
      j = n-1
      du(i, j, b) = coeffs(1)*u(i, j-4, b) + coeffs(2)*u(i, j-3, b) &
                    + coeffs(3)*u(i, j-2, b) + coeffs(4)*u(i, j-1, b) &
                    + coeffs(5)*u(i, j, b) &
                    + coeffs(6)*u(i, j+1, b) + coeffs(7)*u_e(i, 1, b) &
                    + coeffs(8)*u_e(i, 2, b) + coeffs(9)*u_e(i, 3, b)
      du(i, j, b) = ffr(j)*(du(i, j, b) - faf(j)*du(i, j-1, b))
      j = n
      du(i, j, b) = coeffs(1)*u(i, j-4, b) + coeffs(2)*u(i, j-3, b) &
                    + coeffs(3)*u(i, j-2, b) + coeffs(4)*u(i, j-1, b) &
                    + coeffs(5)*u(i, j, b) &
                    + coeffs(6)*u_e(i, 1, b) + coeffs(7)*u_e(i, 2, b) &
                    + coeffs(8)*u_e(i, 3, b) + coeffs(9)*u_e(i, 4, b)
      du(i, j, b) = ffr(j)*(du(i, j, b) - faf(j)*du(i, j-1, b))

      send_u_e(i, 1, b) = du(i, n, b)

      ! Backward pass of the hybrid algorithm
      do j = n - 2, 2, -1
         du(i, j, b) = du(i, j, b) - fbc(j)*du(i, j + 1, b)
      end do
      du(i, 1, b) = last_r*(du(i, 1, b) - fbc(1)*du(i, 2, b))
      send_u_b(i, 1, b) = du(i, 1, b)

   end subroutine der_univ_dist

   attributes(global) subroutine der_univ_subs(du, recv_u_b, recv_u_e, &
                                               n, dist_sa, dist_sc)
      implicit none

      ! Arguments
      real(dp), device, intent(out), dimension(:, :, :) :: du
      real(dp), device, intent(in), dimension(:, :, :) :: recv_u_b, recv_u_e
      real(dp), device, intent(in), dimension(:) :: dist_sa, dist_sc
      integer, value, intent(in) :: n

      ! Local variables
      integer :: i, j, b
      real(dp) :: ur, bl, recp, du_1, du_n

      i = threadIdx%x
      b = blockIdx%x

      bl = dist_sa(1)
      ur = dist_sc(n)
      recp = 1._dp/(1._dp - ur*bl)

      !du(i, 1, b) = recp*(du(i, 1, b) - bl*recv_u_b(i, 1, b))
      !du(i, n, b) = recp*(du(i, n, b) - ur*recv_u_e(i, 1, b))
      du_1 = recp*(du(i, 1, b) - bl*recv_u_b(i, 1, b))
      du_n = recp*(du(i, n, b) - ur*recv_u_e(i, 1, b))

      du(i, 1, b) = du_1
      do j = 2, n-1
         du(i, j, b) = (du(i, j, b) - dist_sa(j)*du_1 - dist_sc(j)*du_n)
      end do
      du(i, n, b) = du_n

   end subroutine der_univ_subs

end module m_cuda_kernels_dist
