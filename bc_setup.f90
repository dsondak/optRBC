module bc_setup

use global
use write_pack

implicit none

contains

subroutine init_bc(aii)

implicit none

real(dp),              intent(in)     :: aii
real(dp)                              :: pnu, qnu
real(dp)                              :: wavex

real(dp), allocatable, dimension(:)   :: d, dl, du
real(dp), allocatable, dimension(:)   :: phi1_b, phi1_t
real(dp), allocatable, dimension(:)   :: phi2_b, phi2_t
real(dp), allocatable, dimension(:,:) :: Fphi12, FV12

integer                               :: ii, jj
integer                               :: info

! Allocate local variables
allocate(Fphi12(Ny-2,2), FV12(Ny-2,2), stat=alloc_err)
call check_alloc_err(alloc_err)
allocate(d(Ny-2), stat=alloc_err)
call check_alloc_err(alloc_err)
allocate(dl(Ny-3), stat=alloc_err)
call check_alloc_err(alloc_err)
allocate(du(Ny-3), stat=alloc_err)
call check_alloc_err(alloc_err)
allocate(phi1_b(Nx), phi1_t(Nx), stat=alloc_err)
call check_alloc_err(alloc_err)
allocate(phi2_b(Nx), phi2_t(Nx), stat=alloc_err)
call check_alloc_err(alloc_err)

Fphi12 = 0.0_dp
FV12   = 0.0_dp
d      = 0.0_dp
dl     = 0.0_dp
du     = 0.0_dp
phi1_b = 1.0_dp
phi1_t = 0.0_dp
phi2_b = 0.0_dp
phi2_t = 1.0_dp

phi1_b(1) = 0.0_dp
phi2_t(1) = 0.0_dp

! Compute phi1 and phi2 from D2 phi_i = 0 where
! D2 = -(alpha*kx)^2 + dy^2 and phi1 satisfies the BCs
! phi1 = 1 at the bottom wall and phi1 = 0 at the top wall.
! phi2 satisfies phi2 = 0 at the bottom wall and phi2 = 1
! at the top wall.

! A comment on how these arrays are formed.  Since we are only
! dealing with Dirichlet conditions at the walls, we only need arrays of size
! Ny-2, i.e. we don't need the values at the walls.  Thus, in the calculations
! below, F(1) would correspond to point number 2 and F(Ny-2) would correspond
! to point Ny-1.

! Some parameters
pnu = nu0*dt*aii

do ii = 1,Nx

   Fphi12 = 0.0_dp

   if (abs(kx(ii)/alpha) > Nf/2) then
      phi1_b(ii) = 0.0_dp
      phi2_t(ii) = 0.0_dp
   end if

   wavex = kx(ii)

   qnu = 1.0_dp + pnu*wavex**2.0_dp

   do jj = 2,Ny-1
      d(jj-1) = qnu - pnu*g2(jj)
   end do

   do jj = 2,Ny-2
      du(jj-1) = -pnu*g3(jj)
   end do

   do jj = 3,Ny-1
      dl(jj-2) = -pnu*g1(jj)
   end do

   Fphi12(1,1)    = pnu*g1(2)*phi1_b(ii)
   Fphi12(Ny-2,2) = pnu*g3(Ny-1)*phi2_t(ii)

   !  Solve the system Aphi phi = Fphi
   call dgtsv(Ny-2, 2, dl, d, du, Fphi12, Ny-2, info)

   !  Put phi1 and phi2 together
   phi1(1,ii)      = phi1_b(ii)
   phi1(2:Ny-1,ii) = Fphi12(:,1)
   phi1(Ny,ii)     = 0.0_dp
   phi2(1,ii)      = 0.0_dp
   phi2(2:Ny-1,ii) = Fphi12(:,2)
   phi2(Ny,ii)     = phi2_t(ii)

   !  Calculate V1 and V2 from D2 V = phi
   !  Note that we used uy at top and bottom = 0 implicitly here.

   FV12(:,1) = phi1(2:Ny-1,ii)
   FV12(:,2) = phi2(2:Ny-1,ii)

   do jj = 2,Ny-1
      d(jj-1) = -wavex**2.0_dp + g2(jj)
   end do

   do jj = 2,Ny-2
      du(jj-1) = g3(jj)
   end do

   do jj = 3,Ny-1
      dl(jj-2) = g1(jj)
   end do

   call dgtsv(Ny-2, 2, dl, d, du, FV12, Ny-2, info)

   V1(2:Ny-1,ii) = FV12(:,1)
   V2(2:Ny-1,ii) = FV12(:,2)

   ! Calculate the wall derivatives.
   dyv1_B(ii) = h1(1)*V1(1,ii) + h2(1)*V1(2,ii) + h3(1)*V1(3,ii)
   dyv2_B(ii) = h1(1)*V2(1,ii) + h2(1)*V2(2,ii) + h3(1)*V2(3,ii)

   dyv1_T(ii) = h1(Ny)*V1(Ny-2,ii) + h2(Ny)*V1(Ny-1,ii) + h3(Ny)*V1(Ny,ii)
   dyv2_T(ii) = h1(Ny)*V2(Ny-2,ii) + h2(Ny)*V2(Ny-1,ii) + h3(Ny)*V2(Ny,ii)

end do ! kx loop

end subroutine init_bc

end module bc_setup
