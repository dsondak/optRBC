module nonlinear_solvers

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine flow_map(Nu, iter_nl, iter_gmres_ave)

use global
use write_pack
use allocate_vars
use imod
use statistics
use jacobians
use gmres_pack
use time_integrators

implicit none

real(dp), allocatable :: xT(:), Delta_x(:)
real(dp)              :: norm_GT, norm_x0, error
real(dp)              :: dnrm2                  ! BLAS L2-norm function
real(dp), intent(out) :: Nu                     ! Nu at top and bottom
real(dp), intent(out) :: iter_gmres_ave
integer , intent(out) :: iter_nl
real(dp), parameter   :: tol_nl = 1.0e-6_dp    ! Nonlinear iteration tolerance
real(dp), parameter   :: tol_gmres = 1.0e-6_dp  ! GMRES tolerance
integer , parameter   :: restart = 0            ! Restart GMRES = 1
integer               :: n, gmres_it
integer               :: iter, flag
integer, parameter    :: incx = 1
integer, parameter    :: nl_max = 25            ! Max nonlinear iterations
integer               :: nli, ii, jj, kk        ! Loop counters
integer               :: ri, rf

procedure(func), pointer :: infunc

n = size(x0)
gmres_it = 500

! Allocate flow map data
allocate(xT(n), Delta_x(n), stat=alloc_err)
if (alloc_err /= 0) then
   write(*,*) "ERROR:  Allocation in flow_map incomplete."
   stop
end if
xT = 0.0_dp
Delta_x = 0.0_dp

call init_random_seed()
call random_number(Delta_x)

! Make sure Delta_x is zero at the boundaries.
do ii = 1,size(b_pntr)
   ri = b_pntr(ii)*Nx - (Nx-1)
   rf = b_pntr(ii)*Nx
   Delta_x(ri:rf) = 0.0_dp
end do

! Point to function for approximating the jacobian based on the flow map
infunc => jac_approx_flow_map

iter_gmres_ave = 0.0_dp
norm_x0 = dnrm2(n,x0,incx)
do nli = 1,nl_max
   nli_global = nli
   ! Get solution at time T
   call imex_rk(.false., .false.)
   !  Form xT
   call packx(xT)
   ! Get flow map
   GT = (xT - x0) / t_final
   ! Compute norm of GT to see if we're at steady state
   norm_GT = dnrm2(n,GT,incx) / norm_x0
   if (norm_GT <= tol_nl) then

      iter_gmres_ave = real(iter_gmres_ave, kind=dp) / real((nli-1), kind=dp)
      iter_nl = nli-1
      if (iter_nl /= 0) then
        write(*,*) " "
        write(*,*) "Found a steady state at nonlinear iteration number ", nli-1
        write(*,*) "The nonlinear norm is ", norm_GT
        write(*,*) "Total number of GMRES iterations is", int(iter_gmres_ave *(nli-1))
        write(*,*) "The number of GMRES iterations per nonlinear step is ", iter_gmres_ave
        write(*,*) " "
      end if
      if (iter_nl == 0) then
              write(*,*) " "
              write(*,*) "Parabolic Fit has found local max of alpha"
              write(*,*) "The nonlinear norm is ", norm_GT
              write(*,*) " "
      end if
      flush(6)

      call nusselt(Nu, .true.)

      exit
   end if

   ! Use GMRES to find \Delta X.  Note that we want to force GT -> 0.
   ! Hence, we would like to take d GT / dX * \Delta X= -GT and find
   ! \Delta X.  The function jac_approx_flow_map does this by approximating
   ! the action of the Jacobian.
   call gmres(iter,error,flag, &
             &gmres_it, restart, tol_gmres, Delta_x, -GT, infunc)
   write(*,*) "nli, iter, gmres_error, nl_error= ", nli, iter, error, norm_GT
   iter_gmres_ave = iter_gmres_ave + iter
   flush(6)
   if (flag == 1) then
      write(*,*) "GMRES FAILED TO CONVERGE."
      write(*,*) "GMRES Iterations: ", iter
      write(*,*) "GMRES Error: ", error
      stop
   end if

   ! Update solution.
   x0 = x0 + Delta_x
   ! Unpack x0 into V, T, T0
   call unpackx(x0)

   ! Calculate phi and ux from uy
   do jj = 1,Nx
      if (kx(jj) /= 0.0_dp) then
         tuy = uy(:,jj)
         ux(:,jj) = -CI*d1y(tuy)/kx(jj)
      else if (kx(jj) == 0.0_dp) then
         ux(:,jj) = cmplx(0.0_dp, 0.0_dp, kind=C_DOUBLE_COMPLEX) ! Zero mean flow!
      end if
      phi(:,jj) = -kx(jj)**2.0_dp*uy(:,jj) + d2y(uy(:,jj))
   end do

!   dt = 0.05_dp

end do

deallocate(xT, Delta_x, stat=dealloc_err)
if (dealloc_err /= 0) then
   write(*,*) "ERROR:  Failed to deallocate xT and Delta_x in flow_map."
   stop
end if

end subroutine flow_map

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module nonlinear_solvers
