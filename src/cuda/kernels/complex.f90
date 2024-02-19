module m_cuda_complex
   use cudafor

   use m_common, only: dp
   use m_cuda_common, only: SZ

   implicit none

contains

   attributes(global) subroutine processfftdiv( &
         div, waves, nx, ny, nz, ax, bx, ay, by, az, bz &
         )
      implicit none

      ! Arguments
      complex(dp), device, intent(inout), dimension(:, :, :) :: div
      complex(dp), device, intent(in), dimension(:, :, :) :: waves
      real(dp), device, intent(in), dimension(:) :: ax, bx, ay, by, az, bz
      integer, value, intent(in) :: nx, ny, nz

      ! Local variables
      integer :: i, j, b, ix, iy, iz
      real(dp) :: tmp_r, tmp_c, div_r, div_c

      i = threadIdx%x
      b = blockIdx%x

      do j = 1, nx
         ! normalisation
         div_r = real(div(i, j, b), kind=dp)/(nx*ny*nz)
         div_c = aimag(div(i, j, b))/(nx*ny*nz)

         ! get the indices for x, y, z directions
         ix = j; iy = i + (b-1)/(nz/2+1)*SZ; iz = mod(b-1, nz/2+1) + 1

         ! post-process forward
         ! post-process in z
         tmp_r = div_r
         tmp_c = div_c
         div_r = tmp_r*bz(iz) + tmp_c*az(iz)
         div_c = tmp_c*bz(iz) - tmp_r*az(iz)

         ! post-process in y
         tmp_r = div_r
         tmp_c = div_c
         div_r = tmp_r*by(iy) + tmp_c*ay(iy)
         div_c = tmp_c*by(iy) - tmp_r*ay(iy)
         if ( iy > ny/2 + 1 ) div_r = -div_r
         if ( iy > ny/2 + 1 ) div_c = -div_c

         ! post-process in x
         tmp_r = div_r
         tmp_c = div_c
         div_r = tmp_r*bx(ix) + tmp_c*ax(ix)
         div_c = tmp_c*bx(ix) - tmp_r*ax(ix)
         if ( ix > nx/2 + 1 ) div_r = -div_r
         if ( ix > nx/2 + 1 ) div_c = -div_c

         ! Solve Poisson
         tmp_r = real(waves(i, j, b), kind=dp)
         tmp_c = aimag(waves(i, j, b))
         if ((tmp_r < 1.e-16_dp) .or. (tmp_c < 1.e-16_dp)) then
            div_r = 0._dp; div_c = 0._dp
         else
            div_r =-div_r/tmp_r
            div_c =-div_c/tmp_c
         end if

         ! post-process backward
         ! post-process in z
         tmp_r = div_r
         tmp_c = div_c
         div_r = tmp_r*bz(iz) - tmp_c*az(iz)
         div_c =-tmp_c*bz(iz) - tmp_r*az(iz)

         ! post-process in y
         tmp_r = div_r
         tmp_c = div_c
         div_r = tmp_r*by(iy) + tmp_c*ay(iy)
         div_c = tmp_c*by(iy) - tmp_r*ay(iy)
         if ( iy > ny/2 + 1 ) div_r = -div_r
         if ( iy > ny/2 + 1 ) div_c = -div_c

         ! post-process in x
         tmp_r = div_r
         tmp_c = div_c
         div_r = tmp_r*bx(ix) + tmp_c*ax(ix)
         div_c =-tmp_c*bx(ix) + tmp_r*ax(ix)
         if ( ix > nx/2 + 1 ) div_r = -div_r
         if ( ix > nx/2 + 1 ) div_c = -div_c

         ! update the entry
         div(i, j, b) = cmplx(div_r, div_c, kind=dp)
      end do

   end subroutine processfftdiv

end module m_cuda_complex
