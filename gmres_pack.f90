module gmres_pack

use global
use imod
use jacobians

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine gmres(iter,error,flag, max_it,restart,toler, x, b, infunc, A)

implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!
!! This subroutine computes the solution to the linear system
!! Ax = b using the GMRES algorithm. Gram-Schmit is used for
!! orthonormalization.  Given's rotations are used for the
!! QR factorization.  The problem size is N.
!!
!! INPUT:  A (optional)      -- Square, possibly unsymmetric matrix; 
!!                              double precision
!!         infunc (optional) -- A function pointer.  This pointer points to
!!                              a function that returns the action of a
!!                              matrix on a vector.
!!         x                 -- Initial iterate; double precision
!!         b                 -- The RHS vector; double precision
!!         max_it(optional)  -- User-supplied maximum number of iterations.
!!                              If not present then defaults to 10; integer
!!         restart(optional) -- Number of GMRES iterations to perform 
!!                              before restarting the algorithm.  If the 
!!                              input value is 0 or N then the alorithm 
!!                              assumes no restarts.  If not present then
!!                              defaults to no restart; integer
!!         toler(optional)   -- User-supplied relative tolerance.  If not
!!                              present then defaults to 1e-6; 
!!                              double precision
!!
!! OUTPUT:  iter  -- Total number of GMRES iterations required.  Will be
!!                   the product of outer and inner iterations.  When
!!                   the un-restarted version is used, this number is
!!                   just the inner iterations as the outer iteration is
!!                   unity; integer
!!          error -- Relative error in the residual; double precision
!!          flag  -- Notifies the user if convergence was acheived or not.
!!                   If the algorithm converged, then flag = 0 otherwise
!!                   flag = 1; integer
!!          x     -- Solution vector; double precision
!!
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Declarations

real(dp), optional,  intent(in)       :: A(:,:)
real(dp), intent(in)                  :: b(:)
real(dp), intent(inout)               :: x(:)
real(dp), optional, intent(in)        :: toler
integer , optional, intent(in)        :: max_it
integer , optional, intent(in)        :: restart
procedure(func), optional, pointer    :: infunc
real(dp), intent(out)                 :: error
real(dp)                              :: dnrm2                  ! BLAS L2-norm function
integer , intent(out)                 :: iter, flag

real(dp), dimension(size(x))          :: r, Ax, w
real(dp), allocatable, dimension(:,:) :: P, H
real(dp), allocatable, dimension(:)   :: cs, sn
real(dp), allocatable, dimension(:)   :: s
real(dp), allocatable, dimension(:)   :: vt, e1
real(dp), allocatable, dimension(:)   :: y, yt

real(dp)                              :: bnrm2, ddot
real(dp)                              :: temp, k1, k2
real(dp)                              :: tol
real(dp), parameter                   :: zed = 1.0_dp
real(dp), parameter                   :: beta  = 0.0_dp

integer, parameter                    :: incx = 1, incy = 1
integer, dimension(2)                 :: sz
integer                               :: rows, cols
integer                               :: maxit
integer                               :: inner, outer
integer                               :: outiter, initer
integer                               :: ii

character(10) :: citer1
character(9)  :: fiter1
character(3)  :: citer2
character(2)  :: fiter2

iter = 0
flag = 0

rows = 0
cols = 0

w = 0.0_dp
Ax = 0.0_dp
r = 0.0_dp

if (present(A)) then
   if (present(infunc)) then
      write(*,*) "ERROR:  You must supply only one of A and infunc!"
      stop
   else
      sz = shape(A)
      rows = sz(1)
      cols = sz(2)
   end if
else if (present(infunc)) then
   if (present(A)) then
      write(*,*) "ERROR:  You must supply only one of A and infunc!"
      stop
   else
      rows = size(x)
   end if
else
   write(*,*) "ERROR:  You must supply either A or a function pointer!"
      stop
end if

if (present(max_it)) then
   if (max_it <= 0) then
      maxit = 10
   else
      maxit = max_it
   end if
else
   maxit = 10
end if

if (present(toler)) then
   tol = toler
else
   tol = 1.0e-6_dp
end if

if (present(restart)) then
   if ( (restart == 0).or.(restart == rows) ) then ! No restarts
      outer = 1
      inner = maxit
   else
      outer = maxit
      inner = restart
   end if
else ! Default to no restarts
   outer = 1
   inner = maxit
end if

allocate(P(rows,inner+1), stat = alloc_err)
allocate(H(inner+1,inner), stat = alloc_err)
allocate(cs(inner), sn(inner), stat = alloc_err)
allocate(e1(rows), vt(rows), stat = alloc_err)
allocate(s(rows+1), stat = alloc_err)

if (alloc_err /= 0) then
    write(*,*) "ERROR: Could not allocate space for arrays!"
    stop
end if

P = 0.0_dp
H = 0.0_dp
cs = 0.0_dp
sn = 0.0_dp
vt = 0.0_dp
s = 0.0_dp

bnrm2 = dnrm2(rows,b,1)
if (present(A)) then
   call dgemv('n', rows, cols, zed, A, rows, x, incx, beta, Ax, incy)
   r = b - Ax
else
   r = b - infunc(x)
end if
error  = norm(r) / bnrm2
e1 = unit_vec(rows,1)

gmres_outer: do outiter = 1,outer
                 if (present(A)) then
                    ! Compute A*x
                    call dgemv('n', rows, cols, zed, A, rows, x, incx, beta, Ax, incy)
                    ! Form residual 
                    r = b - Ax
                 else
                    r = b - infunc(x)
                 end if 
                 P(:,1) = r / norm(r)
                 s(1:rows) = norm(r) * e1
    gmres_inner: do initer = 1,inner
!                     write(*,*) " inner gmres iteration = ", initer
                     ! Use Arnoldi process to construct an orthonormal basis
                     vt = P(:,initer)
                     if (present(A)) then
                        call dgemv('n', rows, cols, zed, A, rows, vt, incx, beta, w, incy)
                     else
                        w = infunc(vt)
                     end if
                     do ii = 1,initer
                         H(ii,initer) = ddot(rows, w, incx, P(:,ii), incy)
                         w = w - H(ii,initer) * P(:,ii)
                     end do
                     H(initer+1,initer) = norm(w)
                     P(:,initer+1) = w / H(initer+1,initer)
                     ! Now apply the Givens rotation to get QR factorization.
                     ! This is the least-squares part of the problem.
                     do ii = 1,initer-1
                         temp           =  cs(ii)*H(ii,initer) + sn(ii)*H(ii+1,initer)
                         H(ii+1,initer) = -sn(ii)*H(ii,initer) + cs(ii)*H(ii+1,initer)
                         H(ii,initer)   =  temp
                     end do
                     ! Form the current rotation matrix
                     k1 = H(initer,initer)
                     k2 = H(initer+1,initer)
                     call drotg(k1, k2, cs(initer), sn(initer))
                     temp = cs(initer)*s(initer)
                     s(initer+1) = -sn(initer)*s(initer)
                     s(initer) = temp
                     H(initer,initer) = cs(initer)*H(initer,initer) + sn(initer)*H(initer+1,initer)
                     H(initer+1,initer) = 0.0_dp
                     ! Check residual norm
                     error = abs(s(initer+1)) / bnrm2
                     ! Think about updating the solution
                     if (error <= tol) then
                         ! Preparing to update
                         allocate(y(initer), stat = alloc_err)
                         allocate(yt(rows), stat = alloc_err)
                         if (alloc_err /= 0) then
                             write(*,*) "ERROR: Could not allocate space for arrays!"
                             stop
                         end if
                         ! Update solution
                         y = s(1:initer)
                         yt = 0.0_dp
                         ! Solve H*y = s for y
                         ! Triangular solve: Upper triangular, no transpose,
                         ! non-unitary diagonal, matrix is of order initer,
                         ! matrix is H(*,*), LDA = initer, RHS is y and result
                         ! will be returned as y, increment y by incx.
                         call dtrsv("u", "n", "n", initer, H(1:initer,1:initer), initer, y, incx)
                         ! Compute P*y
                         call dgemv('n', rows, initer, zed, P(:,1:initer), rows, y, incx, beta, yt, incy)
                         ! Update solution
                         x = x + yt
                         exit
                     end if
                 end do gmres_inner
                 ! Check error and exit outer loop if convergence is acheived
                 if (error <= tol) then
                    exit
                 end if
                 ! Prepare to update solution
                 allocate(y(inner), stat = alloc_err)
                 allocate(yt(rows), stat = alloc_err)
                 if (alloc_err /= 0) then
                     write(*,*) "ERROR: Could not allocate space for arrays!"
                     stop
                 end if
                 ! Update solution
                 y = s(1:inner)
                 yt = 0.0_dp
                 ! Solve H*y = s for y
                 call dtrsv("u", "n", "n", inner, H(1:inner,1:inner), inner, y, incx)
                 ! Compute P*y
                 call dgemv('n', rows, inner, zed, P(:,1:inner), rows, y, incx, beta, yt, incy)
                 ! Update solution
                 x = x + yt
                 ! Compute residual
                 if (present(A)) then
                    ! Compute A*x
                    call dgemv('n', rows, cols, zed, A, rows, x, incx, beta, Ax, incy)
                    ! Compute residual
                    r = b - Ax
                 else
                    r = b - infunc(x)
                 end if
                 s(inner+1) = norm(r)
                 error = s(inner+1) / bnrm2
                 ! Check error and exit if convergence is acheived
                 if (error <= tol) then
                     exit
                 end if
             end do gmres_outer

             ! Total number of GMRES iterations
             iter = initer * outiter

             ! If GMRES doesn't converge then throw a flag indicating this
             if (error > tol) then
                 flag = 1
             end if

             deallocate(P, stat = dealloc_err)
             deallocate(H, stat = dealloc_err)
             deallocate(cs, sn, stat = dealloc_err)
             deallocate(e1, vt, stat = dealloc_err)
             deallocate(y, stat = dealloc_err)
             deallocate(s, stat = dealloc_err)

             if (dealloc_err /= 0) then
                  write(*,*) "ERROR: Could not deallocate arrays!"
                  stop
             end if

2000 format(I3,E25.16E3                  )

end subroutine gmres

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real(dp) function norm(vec) result(n2)

implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This routine computes the L2 norm of a vector
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real(dp), intent(in) :: vec(:)
integer              :: i, M

M = size(vec)
n2 = 0.0_dp

do i = 1,M
 n2 = n2 + vec(i)**2.0_dp
end do
n2 = sqrt(n2)

end function norm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function unit_vec(n,pos) result(uvec)

implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! This routine returns a unit vector
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer, intent(in) :: n, pos
real(dp), allocatable :: uvec(:)
integer :: alloc_err

allocate(uvec(n), stat=alloc_err)
if (alloc_err /= 0) then
   write(*,*) "ERROR:  Allocation error in function unit_vec."
   stop
end if

uvec = 0.0_dp
uvec(pos) = 1.0_dp

end function unit_vec

end module gmres_pack
