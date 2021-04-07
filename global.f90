module global

use fftw

implicit none
save

! Get double-precision type
integer, parameter     :: dp = kind(1.0d0)
! Create fundamental constants
real(dp), parameter    :: pi = 3.141592653589793238462643383279502884197_dp
complex(dp), parameter :: CI = (0.0_dp, 1.0_dp)
complex(dp), parameter :: CR = (1.0_dp, 0.0_dp)
! Create fftw plans
type(C_PTR)            :: planT, planPhi, planux, planuy
type(C_PTR)            :: iplanT, iplanPhi, iplanux, iplanuy
type(C_PTR)            :: plannlT, plannlphi
type(C_PTR)            :: iplannlT, iplannlphi

! Set characters
!character(len=89) :: vtkloc="/fac/sondak/Documents/Research/ThermalConvection/Code/Working/Optimal_Box_Pr/Pr7/vtkdata/"
character(len=10) :: vtkloc="./vtkdata/"
character(len=2)  :: vtkname="Ra"
character(len=4)  :: gtype=".vtk"

! Set boolean parameters
logical :: no_slip = .true.

! Set physical parameters
real(dp) :: Pr
real(dp) :: Ra_global
real(dp) :: t_final

! Set numerical parameters
real(dp) :: ybot
real(dp) :: ytop
real(dp) :: dt, dt_init
integer  :: Nx
integer  :: Ny
integer  :: Nz

real(dp) :: dxmin, dymin
real(dp), parameter :: cfl = 2.0_dp

integer  :: Nf ! Number of modes to actually retain
integer  :: kmax ! Maximum wavenumber shell
real(dp) :: xL
real(dp) :: xR
real(dp) :: zF
real(dp) :: zB
real(dp) :: Lx ! Channel width
real(dp) :: Ly ! Channel height
real(dp) :: Lz ! Channel span
real(dp) :: nu0 ! Dimensionless viscosity
real(dp) :: kappa0 ! Dimensionless heat conductivity
real(dp) :: dx ! Grid spacing in the streamwise direction
real(dp) :: dy ! Grid spacing in the wall-normal direction
real(dp) :: dz ! Grid spacing in the spanwise direction
real(dp) :: alpha ! Fundamental streamwise wavelength
integer  :: nt ! Number of time-steps
integer  :: nli_global

! Allocatable variables
complex(C_DOUBLE_COMPLEX), allocatable :: T(:,:), uy(:,:), phi(:,:), ux(:,:)
complex(C_DOUBLE_COMPLEX), allocatable :: nlT(:,:), nlphi(:,:)
complex(C_DOUBLE_COMPLEX), allocatable :: tT(:), tuy(:), tux(:), tnlT(:), tnlphi(:)
complex(C_DOUBLE_COMPLEX), allocatable :: tphi(:)
complex(C_DOUBLE_COMPLEX), allocatable, dimension(:,:) :: Tptrb



real(dp), allocatable, dimension(:) :: x0, GT
real(dp), allocatable, dimension(:) :: xp, yp, zp
real(dp), allocatable, dimension(:) :: kx, kz
real(dp), allocatable, dimension(:) :: kx_modes

real(dp), allocatable, dimension(:) :: dynu
real(dp), allocatable, dimension(:) :: g1,g2,g3
real(dp), allocatable, dimension(:) :: h1,h2,h3

real(dp), allocatable :: phi1(:,:), phi2(:,:)
real(dp), allocatable :: V1(:,:), V2(:,:)

real(dp), allocatable :: dyv1_T(:), dyv2_T(:)
real(dp), allocatable :: dyv1_B(:), dyv2_B(:)

integer, dimension(8) :: b_pntr ! points to boundary terms
integer               :: alloc_err, dealloc_err

! Time integrator variables and arrays
complex(C_DOUBLE_COMPLEX), allocatable, dimension(:,:) :: phii, uyi, uxi, Ti, Tkappa
complex(C_DOUBLE_COMPLEX), allocatable, dimension(:,:) :: K1_phi, K1_T
complex(C_DOUBLE_COMPLEX), allocatable, dimension(:,:) :: K2_phi, K2_T
complex(C_DOUBLE_COMPLEX), allocatable, dimension(:,:) :: K3_phi, K3_T
complex(C_DOUBLE_COMPLEX), allocatable, dimension(:,:) :: K1hat_phi, K1hat_T
complex(C_DOUBLE_COMPLEX), allocatable, dimension(:,:) :: K2hat_phi, K2hat_T
complex(C_DOUBLE_COMPLEX), allocatable, dimension(:,:) :: K3hat_phi, K3hat_T
complex(C_DOUBLE_COMPLEX), allocatable, dimension(:,:) :: K4hat_phi, K4hat_T

real(dp), dimension(3)   :: b
real(dp), dimension(3,3) :: acoeffs
real(dp), dimension(4,4) :: ahatcoeffs
real(dp)                 :: gmma

complex(C_DOUBLE_COMPLEX), allocatable, dimension(:) :: tmp_phi, tmp_phi1, tmp_T, tmp_uy, tmp_uy1
complex(C_DOUBLE_COMPLEX), allocatable, dimension(:) :: tmp_K_phi, tmp_K_T



contains

subroutine global_params

  implicit none

  Nf = floor(2.0_dp*Nx/3.0_dp)
  kmax = ceiling(Nf/2.0_dp)

  Ly = ytop - ybot

  if (Ny == 1) then
     dy = 1.0_dp
  else
     dy = Ly / (real(Ny,kind=dp) - 1.0_dp)
  end if

  zF =  0.0_dp
  zB =  0.0_dp

  Lz = zF - zB

  if (Nz == 1) then
     dz = 1.0_dp
  else
     dz = Lz / (real(Nz,kind=dp) - 1.0_dp)
  end if

  nt = int(1.0_dp + t_final / dt)

  gmma = 0.4358665215_dp

  b(1) = -3.0_dp*gmma**2.0_dp/2.0_dp + 4.0_dp*gmma - 1.0_dp/4.0_dp
  b(2) =  3.0_dp*gmma**2.0_dp/2.0_dp - 5.0_dp*gmma + 5.0_dp/4.0_dp
  b(3) =  gmma

  ! a
  acoeffs(1,1) = gmma
  acoeffs(2,1) = (1.0_dp-gmma)/2.0_dp
  acoeffs(3,1) = b(1)

  acoeffs(1,2) = 0.0_dp
  acoeffs(2,2) = gmma
  acoeffs(3,2) = b(2)

  acoeffs(1,3) = 0.0_dp
  acoeffs(2,3) = 0.0_dp
  acoeffs(3,3) = gmma

  ! ahat
  ahatcoeffs(1,1) = 0.0_dp
  ahatcoeffs(2,1) = 0.4358665215_dp
  ahatcoeffs(3,1) = 0.3212788860_dp
  ahatcoeffs(4,1) = -0.105858296_dp

  ahatcoeffs(1,2) = 0.0_dp
  ahatcoeffs(2,2) = 0.0_dp
  ahatcoeffs(3,2) = 0.3966543747_dp
  ahatcoeffs(4,2) = 0.5529291479_dp

  ahatcoeffs(1,3) = 0.0_dp
  ahatcoeffs(2,3) = 0.0_dp
  ahatcoeffs(3,3) = 0.0_dp
  ahatcoeffs(4,3) = 0.5529291479_dp

  ahatcoeffs(1,4) = 0.0_dp
  ahatcoeffs(2,4) = 0.0_dp
  ahatcoeffs(3,4) = 0.0_dp
  ahatcoeffs(4,4) = 0.0_dp

end subroutine global_params

subroutine global_params_Ra(Ra)

  implicit none

  real(dp), optional, intent(in) :: Ra

! Derived parameters
  nu0 = sqrt(16.0_dp * Pr / Ra)
  kappa0 = sqrt(16.0_dp / (Pr * Ra))

end subroutine global_params_Ra

subroutine packx(xpack)

implicit none

real(dp), allocatable, intent(out) :: xpack(:)
integer                            :: n_flm, ii, jj, ri, rf, alloc_err

n_flm = 2*2*Ny*Nx

allocate(xpack(n_flm), stat=alloc_err)
if (alloc_err /= 0) then
   write(*,*) "ERROR:  Trouble allocating space in packx."
   stop
end if

do ii = 1,4
   select case (ii)
      case(1)
         do jj = 1,Ny
            rf = jj*Nx + (ii-1)*Nx*Ny
            ri = rf - Nx + 1
            xpack(ri:rf) = real(uy(jj,:))
         end do
      case(2)
         do jj = 1,Ny
            rf = jj*Nx + (ii-1)*Nx*Ny
            ri = rf - Nx + 1
            xpack(ri:rf) = aimag(uy(jj,:))
         end do
      case(3)
         do jj = 1,Ny
            rf = jj*Nx + (ii-1)*Nx*Ny
            ri = rf - Nx + 1
            xpack(ri:rf) = real(T(jj,:))
         end do
      case(4)
         do jj = 1,Ny
            rf = jj*Nx + (ii-1)*Nx*Ny
            ri = rf - Nx + 1
            xpack(ri:rf) = aimag(T(jj,:))
         end do
   end select
end do

end subroutine packx

subroutine unpackx(xsol)

implicit none

real(dp), allocatable, intent(in) :: xsol(:)
integer                            :: ii, jj, rre, rim

do ii = 1,Nx
   do jj = 1,Ny
      rre = jj*Nx + ii - Nx
      rim = rre + Nx*Ny
      uy(jj,ii) = cmplx(xsol(rre), xsol(rim), kind=C_DOUBLE_COMPLEX)
      rre = rre + 2*Nx*Ny
      rim = rre + Nx*Ny
      T(jj,ii) = cmplx(xsol(rre), xsol(rim), kind=C_DOUBLE_COMPLEX)
   end do
end do

end subroutine unpackx

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

function D2(vec) result(Lvec)

implicit none

real(dp), intent(in)  :: vec(:)
real(dp), allocatable :: Lvec(:)
integer               :: alloc_err, jj, n

n = size(vec)

allocate(Lvec(n), stat=alloc_err)
if (alloc_err /= 0) then
   write(*,*) "ERROR:  Allocation problem in D2."
   stop
end if

Lvec = 0.0_dp

do jj = 2,n-1
   Lvec(jj) = g1(jj)*vec(jj-1) + (-alpha**2.0_dp + g2(jj))*vec(jj) + g3(jj)*vec(jj+1)
end do

if (no_slip) then
   Lvec(1) = g2(2)*vec(2)   + g3(2)*vec(3)
   Lvec(n) = g1(n)*vec(n-2) + g2(n)*vec(n-1)
end if

end function D2

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

subroutine init_random_seed()
implicit none

integer :: i, n, clock
integer, dimension(:), allocatable :: seed

call random_seed(size = n)

allocate(seed(n))

call system_clock(count=clock)
seed = clock + 37 * (/ (i - 1, i = 1, n) /)
call random_seed(put = seed)

deallocate(seed)

end subroutine

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

function d1y(vec) result(d1yvec)

implicit none

complex(C_DOUBLE_COMPLEX), dimension(:), intent(in)  :: vec
complex(C_DOUBLE_COMPLEX), allocatable, dimension(:) :: d1yvec
integer               :: alloc_err, jj, n

n = size(vec)

allocate(d1yvec(n), stat=alloc_err)
if (alloc_err /= 0) then
   write(*,*) "ERROR:  Allocation problem in d1y."
   stop
end if

d1yvec = 0.0_dp

d1yvec(1) = h1(1)*vec(1) + h2(1)*vec(2) + h3(1)*vec(3)
do jj = 2,n-1
   d1yvec(jj) = h1(jj)*vec(jj-1) + h2(jj)*vec(jj) + h3(jj)*vec(jj+1)
end do
d1yvec(n) = h1(n)*vec(n-2) + h2(n)*vec(n-1) + h3(n)*vec(n)

end function d1y

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

function d2y(vec) result(d2yvec)

implicit none

complex(C_DOUBLE_COMPLEX), intent(in)  :: vec(:)
complex(C_DOUBLE_COMPLEX), allocatable :: d2yvec(:)
real(dp)                               :: d1, d2, d3, d4
real(dp)                               :: ell1, ell2, ell3, ell4
real(dp)                               :: dNm3, dNm2, dNm1, dN
real(dp)                               :: ellNm3, ellNm2, ellNm1, ellN
integer                                :: jj, n

n = size(vec)

allocate(d2yvec(n), stat=alloc_err)
if (alloc_err /= 0) then
   write(*,*) "ERROR:  Allocation problem in d2y."
   stop
end if

d2yvec = 0.0_dp

d1 = -dynu(1)*(dynu(1)+dynu(2))*(dynu(1)+dynu(2)+dynu(3))
d2 =  dynu(1)*dynu(2)*(dynu(2) + dynu(3))
d3 = -dynu(2)*dynu(3)*(dynu(1) + dynu(2))
d4 =  dynu(3)*(dynu(2)+dynu(3))*(dynu(1)+dynu(2)+dynu(3))

ell1 = -2.0_dp*(3.0_dp*dynu(1) + 2.0_dp*dynu(2) + dynu(3)) / d1
ell2 = -2.0_dp*(2.0_dp*(dynu(1) + dynu(2)) + dynu(3)) / d2
ell3 = -2.0_dp*(2.0_dp*dynu(1) + dynu(2) + dynu(3)) / d3
ell4 = -2.0_dp*(2.0_dp*dynu(1) + dynu(2)) / d4

dN   =  dynu(Ny-1)*(dynu(Ny-1)+dynu(Ny-2))*(dynu(Ny-1)+dynu(Ny-2)+dynu(Ny-3))
dNm1 = -dynu(Ny-1)*dynu(Ny-2)*(dynu(Ny-2)+dynu(Ny-3))
dNm2 =  dynu(Ny-2)*dynu(Ny-3)*(dynu(Ny-1)+dynu(Ny-2))
dNm3 = -dynu(Ny-3)*(dynu(Ny-2)+dynu(Ny-3))*(dynu(Ny-1)+dynu(Ny-2)+dynu(Ny-3))

ellN   = 2.0_dp*(3.0_dp*dynu(Ny-1) + 2.0_dp*dynu(Ny-2) + dynu(Ny-3)) / dN
ellNm1 = 2.0_dp*(2.0_dp*(dynu(Ny-1)+dynu(Ny-2)) + dynu(Ny-3)) / dNm1
ellNm2 = 2.0_dp*(2.0_dp*dynu(Ny-1) + dynu(Ny-2) + dynu(Ny-3)) / dNm2
ellNm3 = 2.0_dp*(2.0_dp*dynu(Ny-1) + dynu(Ny-2)) / dNm3

!d2yvec(1) = g1(1)*vec(1) + g2(1)*vec(2) + g3(1)*vec(3)
d2yvec(1) = ell1*vec(1) + ell2*vec(2) + ell3*vec(3) + ell4*vec(4)
do jj = 2,n-1
   d2yvec(jj) = g1(jj)*vec(jj-1) + g2(jj)*vec(jj) + g3(jj)*vec(jj+1)
end do
!d2yvec(n) = g1(n)*vec(n-2) + g2(n)*vec(n-1) + g3(n)*vec(n)
d2yvec(n) = ellNm3*vec(n-3) + ellNm2*vec(n-2) + &
           &ellNm1*vec(n-1) + ellN*vec(n)

end function d2y

subroutine check_alloc_err(err_mess)

integer, intent(in) :: err_mess

if (err_mess /= 0) then
  write(*,*) "Allocation error!"
  stop
end if

end subroutine check_alloc_err

end module global
