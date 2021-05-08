module jacobians

use global
use write_pack

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function matrix_action(x) result(y)

implicit none

real(dp), intent(in)  :: x(:)
real(dp), allocatable :: y(:)
real(dp), allocatable :: x1(:), x2(:), x3(:)
real(dp) :: jj
integer :: n, ii

n = size(x)

allocate(x1(n), x2(n), x3(n), y(n), stat=alloc_err)
x1 = 0.0_dp
x2 = 0.0_dp
x3 = 0.0_dp
y  = 0.0_dp

x1(1) = 0.0_dp
x1(2:n) = x(1:n-1)

x2(1) = real((n-1) / 2, kind=dp)
jj = 1.0_dp
do ii = 2,n
   if (ii <= (n+1)/2) then
      x2(ii) = x2(ii-1) - 1.0_dp
   else
      x2(ii) = jj
      jj = jj + 1
   end if
end do

x2 = x2*x

x3 = x(2:n)

y = x1 + x2 + x3

end function matrix_action

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function jac_approx_flow_map(x) result(jact)

use global
use time_integrators

implicit none

real(dp), intent(in) :: x(:)
real(dp), allocatable :: jact(:)

real(dp), allocatable :: x0_delta(:), xT_delta(:), GT_delta(:)

real(dp) :: normx, dot_prod, eps
real(dp) :: dnrm2, ddot

integer :: n, ii
integer, parameter :: incx = 1, incy = 1

n = size(x)

allocate(jact(n), stat=alloc_err)
allocate(x0_delta(n), xT_Delta(n), GT_Delta(n), stat=alloc_err)
call check_alloc_err(alloc_err)

jact = 0.0_dp
x0_delta = 0.0_dp
xT_delta = 0.0_dp
GT_delta = 0.0_dp

normx = dnrm2(n,x,incx)
dot_prod = ddot(n,x,incx,x0,incy)

eps = sqrt(epsilon(1.0_dp)) / (normx**2.0_dp) * &
     &max(dot_prod, normx)

x0_delta = x0 + eps * x

! Upack fields
call unpackx(x0_delta)

! Calculate phi and ux from uy
do ii = 1,Nx
   if (kx(ii) /= 0.0_dp) then
      tmp_uy = uy(:,ii)
      ux(:,ii) = -CI*d1y(tmp_uy)/kx(ii)
   else if (kx(ii) == 0.0_dp) then
      ux(:,ii) = cmplx(0.0_dp, 0.0_dp, kind=C_DOUBLE_COMPLEX) ! Zero mean flow!
   end if
   phi(:,ii) = -kx(ii)**2.0_dp*uy(:,ii) + d2y(uy(:,ii))
end do

! Integrate out to some time
call imex_rk(.false., .false.)

! Form xT_delta from the fields
call packx(xT_delta)

! Compute flow map
GT_delta = (xT_delta - x0_delta) / t_final
jact = (GT_delta - GT) / eps

end function jac_approx_flow_map

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module jacobians
