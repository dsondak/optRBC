module interpolation_pack

use global
use write_pack

contains

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
subroutine bilinear_interp(xc,yc,uyc,Tc)

real(dp),                  dimension(:),   intent(in) :: xc, yc
complex(C_DOUBLE_COMPLEX), dimension(:,:), intent(in) :: uyc, Tc
real(dp)                                              :: x1, x2, y1, y2
real(dp)                                              :: xi, eta
real(dp)                                              :: dx
integer                                               :: ii, jj
integer                                               :: ic, jc

ic = 1 ! Counter
! Initialize x1 and x2
x1 = xc(ic)
x2 = xc(ic+1)
dx = xp(2) - xp(1)
do ii = 1,Nx ! Loop over x
   if (ii == Nx) then
      ! Initialize x1 and x2
      x1 = xc(ic)
      x2 = abs(xc(1))
      xi = xp(ii) ! Locate x position on refined mesh
      jc = 1 ! Counter
      ! Initialize y1 and y2
      y1 = yc(jc)
      y2 = yc(jc+1)
      do jj = 1,Ny
         eta = yp(jj) ! Locate y position on refined mesh
         if (eta > y2) then! Check y interpolation points
            jc = jc + 1
            y1 = yc(jc)
            y2 = yc(jc+1)
         end if
         ! Set up coefficients for interpolation formula
         c0 = (x2 - x1)*(y2  - y1)
         c1 = (x2 - xi)*(y2  - eta)
         c2 = (x2 - xi)*(eta - y1)
         c3 = (xi - x1)*(y2  - eta)
         c4 = (xi - x1)*(eta - y1)
         ! Interpolate velocity
         uy(jj,ii) = (c1*uyc(jc,ic) + c2*uyc(jc+1,ic) + &
       &              c3*uyc(jc,1) + c4*uyc(jc+1,1)) / c0
         ! Interpolate temperature
         T(jj,ii) = (c1*Tc(jc,ic) + c2*Tc(jc+1,ic) + &
       &              c3*Tc(jc,1) + c4*Tc(jc+1,1)) / c0
      end do
   else
      xi = xp(ii) ! Locate x position on refined mesh
      !if (xi > x2) then ! Check interpolation points
      if ((xi - x2) >= dx/2.0_dp) then ! Check interpolation points
         ic = ic + 1
         x1 = xc(ic)
         x2 = xc(ic+1)
      end if
      jc = 1 ! Counter
      ! Initialize y1 and y2
      y1 = yc(jc)
      y2 = yc(jc+1)
      do jj = 1,Ny ! Loop over y
         eta = yp(jj) ! Locate y position on refined mesh
         if (eta > y2) then! Check y interpolation points
            jc = jc + 1
            y1 = yc(jc)
            y2 = yc(jc+1)
         end if
         ! Set up coefficients for interpolation formula
         c0 = (x2 - x1)*(y2  - y1)
         c1 = (x2 - xi)*(y2  - eta)
         c2 = (x2 - xi)*(eta - y1)
         c3 = (xi - x1)*(y2  - eta)
         c4 = (xi - x1)*(eta - y1)
         ! Interpolate velocity
         uy(jj,ii) = (c1*uyc(jc,ic) + c2*uyc(jc+1,ic) + &
       &              c3*uyc(jc,ic+1) + c4*uyc(jc+1,ic+1)) / c0
         ! Interpolate temperature
         T(jj,ii) = (c1*Tc(jc,ic) + c2*Tc(jc+1,ic) + &
       &              c3*Tc(jc,ic+1) + c4*Tc(jc+1,ic+1)) / c0
      end do
   end if 
end do

end subroutine bilinear_interp
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
subroutine linear_interp_x(xc,uyc,Tc)

real(dp),                  dimension(:),   intent(in) :: xc
complex(C_DOUBLE_COMPLEX), dimension(:,:), intent(in) :: uyc, Tc
real(dp)                                              :: x1, x2
real(dp)                                              :: xi
real(dp)                                              :: dx
integer                                               :: ii
integer                                               :: ic

ic = 1 ! Counter
! Initialize x1 and x2
x1 = xc(ic)
x2 = xc(ic+1)
dx = xp(2) - xp(1)
do ii = 1,Nx ! Loop over x
   xi = xp(ii) ! Locate x position on refined mesh
   if (ii == Nx) then
      ! Initialize x1 and x2
      x1 = xc(ic)
      x2 = abs(xc(1))
      ! Set up coefficients for interpolation formula
      c0 = (x2 - x1)
      c1 = (x2 - xi)
      c3 = (xi - x1)
      ! Interpolate velocity
      uy(:,ii) = (c1*uyc(:,ic) + c3*uyc(:,1)) / c0
      ! Interpolate temperature
      T(:,ii) = (c1*Tc(:,ic) + c3*Tc(:,1)) / c0
   else
      if ((xi - x2) >= dx/2.0_dp) then ! Check interpolation points
         ic = ic + 1
         x1 = xc(ic)
         x2 = xc(ic+1)
      end if
         ! Set up coefficients for interpolation formula
         c0 = (x2 - x1)
         c1 = (x2 - xi)
         c3 = (xi - x1)
         ! Interpolate velocity
         uy(:,ii) = (c1*uyc(:,ic) + c3*uyc(:,ic+1)) / c0
         ! Interpolate temperature
         T(:,ii)  = (c1*Tc(:,ic)  + c3*Tc(:,ic+1))  / c0
   end if 
end do

end subroutine linear_interp_x
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
subroutine linear_interp_y(yc,uyc,Tc)

real(dp),                  dimension(:),   intent(in) :: yc
complex(C_DOUBLE_COMPLEX), dimension(:,:), intent(in) :: uyc, Tc
real(dp)                                              :: y1, y2
real(dp)                                              :: eta
integer                                               :: jj
integer                                               :: jc

jc = 1 ! Counter
! Initialize x1 and x2
y1 = yc(jc)
y2 = yc(jc+1)
do jj = 1,Ny ! Loop over x
   eta = yp(jj) ! Locate x position on refined mesh
   if (eta > y2) then ! Check interpolation points
      jc = jc + 1
      y1 = yc(jc)
      y2 = yc(jc+1)
   end if
   ! Set up coefficients for interpolation formula
   c0 = (y2 - y1)
   c1 = (y2 - eta)
   c2 = (eta - y1)
   ! Interpolate velocity
   uy(jj,:) = (c1*uyc(jc,:) + c2*uyc(jc+1,:)) / c0
   ! Interpolate temperature
   T(jj,:)  = (c1*Tc(jc,:)  + c2*Tc(jc+1,:))  / c0
end do

end subroutine linear_interp_y
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

end module interpolation_pack
