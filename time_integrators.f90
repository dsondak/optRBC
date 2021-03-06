module time_integrators

use fftw
use global
use write_pack
use allocate_vars
use bc_setup
use statistics

integer  :: it, jt, kkt
integer  :: info
real(dp) :: time, dtmax, dtmin, dt_old, dt_ramp, dt_final

contains

subroutine imex_rk(vtk_print, save_nusselt)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  This progam solves the equations of thermal convection using a Fourier
!!  spectral method in the x-direction and a 2nd order finite difference scheme in
!!  the y-direction.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


implicit none

integer, optional, intent(in)  :: vtk_print
logical, optional, intent(in)  :: save_nusselt

integer                        :: nti
integer                        :: nprint
logical                        :: wvtk
logical                        :: wnusselt
real(dp)                       :: nusselt_num

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if (present(vtk_print)) then
   wvtk = .true.
   nprint = vtk_print
else
   wvtk = .false.
   nprint = 100
end if

if (wvtk) then
   call write_to_vtk(0, .false.) ! false = Fourier space
end if

if (present(save_nusselt)) then
   wnusselt = .true.
else
   wnusselt = .false.
end if

dt = dt_init

call init_bc(acoeffs(1,1))

time = 0.0_dp

dtmax = 0.5_dp
dtmin = 1.0e-4_dp

dt_ramp = 1.1_dp

dt_old = dt

nti = 0

! Format for writing out single values.
1000 format(E25.16E3)

! Time integration
do ! while (time < t_final)

   dt_final = t_final - time

   if (dt_final <= dt) then
      time = t_final
   else
      time = time + dt
   end if

   !write(*,*) "time = ", time, "dt = ", dt

   nti = nti + 1

   !:::::::::::
   ! STAGE 1 ::
   !:::::::::::
   phii = phi
   Ti   = T
   uxi  = ux
   uyi  = uy
   call calc_explicit(1)
   do it = 1,Nx ! kx loop
      ! Compute phi1 and T1
      call calc_vari(tmp_phi, tmp_T, acoeffs(1,1), 1)
      ! Compute v1 from phi1 
      call calc_vi(tmp_uy, tmp_phi)
      ! BOUNDAY CONDITIONS!
      call update_bcs(tmp_phi1,tmp_uy1, tmp_phi,tmp_uy)
      tmp_phi = tmp_phi1
      tmp_uy  = tmp_uy1
      ! Compute K1_T and K1_phi
      call calc_implicit(tmp_K_phi,tmp_K_T, tmp_phi,tmp_T)
      K1_phi(:,it) = tmp_K_phi
      K1_T(:,it)   = tmp_K_T
      ! Compute u1 from v1
      if (kx(it) /= 0.0_dp) then
         !uxi(:,it) = -CI*d1y(tmp_uy)/kx(it)
         uxi(:,it) = CI*d1y(tmp_uy)/kx(it)
      else if (kx(it) == 0.0_dp) then
         uxi(:,it) = cmplx(0.0_dp, 0.0_dp, kind=C_DOUBLE_COMPLEX) ! Zero mean flow!
      end if
      phii(:,it) = tmp_phi
      Ti  (:,it) = tmp_T
      uyi (:,it) = tmp_uy
   end do
   ! Compute K2hat
   call calc_explicit(2)

   !:::::::::::
   ! STAGE 2 ::
   !:::::::::::
   do it = 1,Nx ! kx loop
      ! Compute phi2 and T2
      call calc_vari(tmp_phi, tmp_T, acoeffs(2,2), 2)
      ! Compute v1 from phi1 
      call calc_vi(tmp_uy, tmp_phi)
      ! BOUNDAY CONDITIONS!
      call update_bcs(tmp_phi1,tmp_uy1, tmp_phi,tmp_uy)
      tmp_phi = tmp_phi1
      tmp_uy  = tmp_uy1
      ! Compute K2_T and K2_phi
      call calc_implicit(tmp_K_phi,tmp_K_T, tmp_phi,tmp_T)
      K2_phi(:,it) = tmp_K_phi
      K2_T(:,it)   = tmp_K_T
      ! Compute u1 from v1
      if (kx(it) /= 0.0_dp) then
         !uxi(:,it) = -CI*d1y(tmp_uy)/kx(it)
         uxi(:,it) = CI*d1y(tmp_uy)/kx(it)
      else if (kx(it) == 0.0_dp) then
         uxi(:,it) = cmplx(0.0_dp, 0.0_dp, kind=C_DOUBLE_COMPLEX) ! Zero mean flow!
      end if
      phii(:,it) = tmp_phi
      Ti  (:,it) = tmp_T
      uyi (:,it) = tmp_uy
   end do
   ! Compute K3hat
   call calc_explicit(3)

   !:::::::::::
   ! STAGE 3 ::
   !:::::::::::
   do it = 1,Nx ! kx loop
      ! Compute phi3 and T3
      call calc_vari(tmp_phi, tmp_T, acoeffs(3,3), 3)
      ! Compute v1 from phi1
      call calc_vi(tmp_uy, tmp_phi)
      ! BOUNDAY CONDITIONS!
      call update_bcs(tmp_phi1,tmp_uy1, tmp_phi,tmp_uy)
      tmp_phi = tmp_phi1
      tmp_uy  = tmp_uy1
      ! Compute K3_T and K3_phi
      call calc_implicit(tmp_K_phi,tmp_K_T, tmp_phi,tmp_T)
      K3_phi(:,it) = tmp_K_phi
      K3_T(:,it)   = tmp_K_T
      ! Compute u1 from v1
      if (kx(it) /= 0.0_dp) then
         !uxi(:,it) = -CI*d1y(tmp_uy)/kx(it)
         uxi(:,it) = CI*d1y(tmp_uy)/kx(it)
      else if (kx(it) == 0.0_dp) then
         uxi(:,it) = cmplx(0.0_dp, 0.0_dp, kind=C_DOUBLE_COMPLEX) ! Zero mean flow!
      end if
      phii(:,it) = tmp_phi
      Ti  (:,it) = tmp_T
      uyi (:,it) = tmp_uy
   end do
   ! Compute K4hat
   call calc_explicit(4)

   ! UPDATE SOLUTIONS

   ! Get phi
   phi(2:Ny-1,:) = phi(2:Ny-1,:) + dt*(b(1)*(K1_phi(2:Ny-1,:) + K2hat_phi(2:Ny-1,:)) + &
                  &                    b(2)*(K2_phi(2:Ny-1,:) + K3hat_phi(2:Ny-1,:)) + &
                  &                    b(3)*(K3_phi(2:Ny-1,:) + K4hat_phi(2:Ny-1,:)))

   ! Get temperature
   T(2:Ny-1,:)   = T(2:Ny-1,:)  + dt*(b(1)*(K1_T(2:Ny-1,:) + K2hat_T(2:Ny-1,:)) + &
                  &                   b(2)*(K2_T(2:Ny-1,:) + K3hat_T(2:Ny-1,:)) + &
                  &                   b(3)*(K3_T(2:Ny-1,:) + K4hat_T(2:Ny-1,:)))

   ! Get ux and uy
   do it = 1,Nx
      ! Solve for v
      call calc_vi(tmp_uy, phi(:,it))
      uy(:,it) = tmp_uy
      ! Solve for u
      if (kx(it) /= 0.0_dp) then
         !ux(:,it) = -CI*d1y(tmp_uy)/kx(it)
         ux(:,it) = CI*d1y(tmp_uy)/kx(it)
      else if (kx(it) == 0.0_dp) then
         ux(:,it) = cmplx(0.0_dp, 0.0_dp, kind=C_DOUBLE_COMPLEX) ! Zero mean flow!
      end if
   end do

   if (time == t_final) then
      exit
   end if

   !call update_dt

   if (wvtk) then
      if (mod(nti,vtk_print) == 0) then
         call write_to_vtk(nti, .false.) ! false = Fourier space
      end if
   end if

   ! Calculate nusselt number.
   if (wnusselt) then
      call nusselt(nusselt_num, .true.) ! true = Fourier space
      write(8000, fmt=1000) nusselt_num
      flush(8000)
   end if

end do ! time loop

end subroutine imex_rk

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
subroutine calc_vari(phiout,Tout, aii, stage)

real(dp),                                             intent(in)  :: aii
integer,                                              intent(in)  :: stage
complex(C_DOUBLE_COMPLEX), allocatable, dimension(:), intent(out) :: phiout, Tout
complex(C_DOUBLE_COMPLEX), allocatable, dimension(:)              :: Fphi, FT
real(dp),                  allocatable, dimension(:,:)            :: phi_rhs, T_rhs
real(dp),                  allocatable, dimension(:)              :: dphi, duphi, dlphi
real(dp),                  allocatable, dimension(:)              :: ddT, duT, dlT

allocate(dphi(Ny-2), ddT(Ny-2), stat=alloc_err)
call check_alloc_err(alloc_err)
allocate(duphi(Ny-3), dlphi(Ny-3), duT(Ny-3), dlT(Ny-3), stat=alloc_err)
call check_alloc_err(alloc_err)
allocate(phiout(Ny), Tout(Ny), stat=alloc_err)
call check_alloc_err(alloc_err)
allocate(Fphi(Ny-2), FT(Ny-2), stat=alloc_err)
call check_alloc_err(alloc_err)
allocate(phi_rhs(Ny-2,2), T_rhs(Ny-2,2), stat=alloc_err)
call check_alloc_err(alloc_err)

dphi    = 0.0_dp
ddT     = 0.0_dp
duphi   = 0.0_dp
dlphi   = 0.0_dp
duT     = 0.0_dp
dlT     = 0.0_dp
phi_rhs = 0.0_dp
T_rhs   = 0.0_dp
Fphi    = (0.0_dp, 0.0_dp)
FT      = (0.0_dp, 0.0_dp)
phiout  = (0.0_dp, 0.0_dp)
Tout    = (0.0_dp, 0.0_dp)

! LHS Matrix (tridiagonal, not necessarily symmetric)
do jt = 2,Ny-1
   ddT (jt-1) = 1.0_dp - kappa0*dt*aii*(-kx(it)**2.0_dp + g2(jt))
   dphi(jt-1) = 1.0_dp - nu0   *dt*aii*(-kx(it)**2.0_dp + g2(jt))
end do

do jt = 2,Ny-2
   duT  (jt-1) = -kappa0*dt*g3(jt)*aii
   duphi(jt-1) = -nu0   *dt*g3(jt)*aii
end do

do jt = 3,Ny-1
   dlT  (jt-2) = -kappa0*dt*g1(jt)*aii
   dlphi(jt-2) = -nu0   *dt*g1(jt)*aii
end do

select case (stage)
  case(1)
    Fphi     = phi(2:Ny-1,it) + dt*ahatcoeffs(2,1)*K1hat_phi(2:Ny-1,it)
    FT       = T  (2:Ny-1,it) + dt*ahatcoeffs(2,1)*K1hat_T  (2:Ny-1,it)
    FT(1)    = FT(1) + kappa0*dt*aii*g1(2)*T(1,it) ! b/c Ti(y_1) = T(y_1)
    FT(Ny-2) = FT(Ny-2) + kappa0*dt*aii*g3(Ny-1)*T(Ny,it) ! b/c Ti(Ny) = T(Ny)

    phi_rhs(:,1) = real(Fphi)
    phi_rhs(:,2) = aimag(Fphi)

    T_rhs  (:,1) = real(FT)
    T_rhs  (:,2) = aimag(FT)

    call dgtsv(Ny-2, 2, dlT, ddT, duT, T_rhs, Ny-2, info)
    Tout(2:Ny-1) = cmplx(T_rhs(:,1), T_rhs(:,2), kind=C_DOUBLE_COMPLEX)
    ! Set temperature boundary conditions 
    Tout(1) = T(1,it)
    Tout(Ny) = T(Ny,it)

    call dgtsv(Ny-2, 2, dlphi, dphi, duphi, phi_rhs, Ny-2, info)
    phiout(2:Ny-1) = cmplx(phi_rhs(:,1), phi_rhs(:,2), kind=C_DOUBLE_COMPLEX)

  case(2)
    Fphi = phi(2:Ny-1,it) + dt*(acoeffs(2,1)*K1_phi(2:Ny-1,it)       + &
          &                     ahatcoeffs(3,1)*K1hat_phi(2:Ny-1,it) + &
          &                     ahatcoeffs(3,2)*K2hat_phi(2:Ny-1,it))
    FT   = T  (2:Ny-1,it) + dt*(acoeffs(2,1)*K1_T  (2:Ny-1,it)       + &
          &                     ahatcoeffs(3,1)*K1hat_T  (2:Ny-1,it) + &
          &                     ahatcoeffs(3,2)*K2hat_T  (2:Ny-1,it))
    FT(1)    = FT(1) + kappa0*dt*aii*g1(2)*T(1,it) ! b/c Ti(y_1) = T(y_1)
    FT(Ny-2) = FT(Ny-2) + kappa0*dt*aii*g3(Ny-1)*T(Ny,it) ! b/c Ti(Ny) = T(Ny)

    phi_rhs(:,1) = real(Fphi)
    phi_rhs(:,2) = aimag(Fphi)

    T_rhs  (:,1) = real(FT)
    T_rhs  (:,2) = aimag(FT)

    call dgtsv(Ny-2, 2, dlT, ddT, duT, T_rhs, Ny-2, info)
    Tout(2:Ny-1) = cmplx(T_rhs(:,1), T_rhs(:,2), kind=C_DOUBLE_COMPLEX)
    ! Set temperature boundary conditions 
    Tout(1) = T(1,it)
    Tout(Ny) = T(Ny,it)

    call dgtsv(Ny-2, 2, dlphi, dphi, duphi, phi_rhs, Ny-2, info)
    phiout(2:Ny-1) = cmplx(phi_rhs(:,1), phi_rhs(:,2), kind=C_DOUBLE_COMPLEX)

  case(3)
    Fphi = phi(2:Ny-1,it) + dt*(acoeffs(3,1)*K1_phi(2:Ny-1,it)       + &
                               &acoeffs(3,2)*K2_phi(2:Ny-1,it)       + &
                               &ahatcoeffs(4,1)*K1hat_phi(2:Ny-1,it) + &
                               &ahatcoeffs(4,2)*K2hat_phi(2:Ny-1,it) + &
                               &ahatcoeffs(4,3)*K3hat_phi(2:Ny-1,it))
    FT   = T  (2:Ny-1,it) + dt*(acoeffs(3,1)*K1_T  (2:Ny-1,it)       + &
                               &acoeffs(3,2)*K2_T  (2:Ny-1,it)       + &
                               &ahatcoeffs(4,1)*K1hat_T  (2:Ny-1,it) + &
                               &ahatcoeffs(4,2)*K2hat_T  (2:Ny-1,it) + &
                               &ahatcoeffs(4,3)*K3hat_T  (2:Ny-1,it))
    FT(1)    = FT(1) + kappa0*dt*aii*g1(2)*T(1,it) ! b/c Ti(y_1) = T(y_1)
    FT(Ny-2) = FT(Ny-2) + kappa0*dt*aii*g3(Ny-1)*T(Ny,it) ! b/c Ti(Ny) = T(Ny)

    phi_rhs(:,1) = real(Fphi)
    phi_rhs(:,2) = aimag(Fphi)
 
    T_rhs  (:,1) = real(FT)
    T_rhs  (:,2) = aimag(FT)

    call dgtsv(Ny-2, 2, dlT, ddT, duT, T_rhs, Ny-2, info)
    Tout(2:Ny-1) = cmplx(T_rhs(:,1), T_rhs(:,2), kind=C_DOUBLE_COMPLEX)
    ! Set temperature boundary conditions 
    Tout(1) = T(1,it)
    Tout(Ny) = T(Ny,it)

    call dgtsv(Ny-2, 2, dlphi, dphi, duphi, phi_rhs, Ny-2, info)
    phiout(2:Ny-1) = cmplx(phi_rhs(:,1), phi_rhs(:,2), kind=C_DOUBLE_COMPLEX)

end select

end subroutine calc_vari
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
subroutine calc_implicit(Kphi,KT, phiin,Tin)

complex(C_DOUBLE_COMPLEX),              dimension(:), intent(in)  :: phiin, Tin
complex(C_DOUBLE_COMPLEX), allocatable, dimension(:), intent(out) :: Kphi, KT

allocate(Kphi(Ny), KT(Ny), stat=alloc_err)
call check_alloc_err(alloc_err)

Kphi = cmplx(0.0_dp, 0.0_dp, kind=C_DOUBLE_COMPLEX)
KT   = cmplx(0.0_dp, 0.0_dp, kind=C_DOUBLE_COMPLEX)

Kphi = nu0   *(-kx(it)**2.0_dp*phiin + d2y(phiin))
KT   = kappa0*(-kx(it)**2.0_dp*Tin   + d2y(Tin))

end subroutine calc_implicit

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
subroutine calc_explicit(stage)

integer             :: i, j
integer, intent(in) :: stage

select case(stage)
   case (1)
      do i = 1,Nx
         K1hat_phi(:,i) = -kx(i)**2.0_dp*Ti(:,i)
      end do
   case (2)
      do i = 1,Nx
         K2hat_phi(:,i) = -kx(i)**2.0_dp*Ti(:,i)
      end do
   case (3)
      do i = 1,Nx
         K3hat_phi(:,i) = -kx(i)**2.0_dp*Ti(:,i)
      end do
   case (4)
      do i = 1,Nx
         K4hat_phi(:,i) = -kx(i)**2.0_dp*Ti(:,i)
      end do
end select

do i=1,Nx
   ! Compute dx(T) in Fourier space
   nlT  (:,i) =  kx(i)*Ti(:,i)
   ! Compute D2(ux)
   nlphi(:,i) = -kx(i)**2.0_dp*uxi(:,i) + d2y(uxi(:,i))
end do
!nlT = -CI*nlT
nlT = CI*nlT

do j = 1,Ny
   ! Bring everything to physical space
   tnlT   = nlT(j,:)
   tnlphi = nlphi(j,:)
   tT     = Ti(j,:)
   tux    = uxi(j,:)
   tuy    = uyi(j,:)
   tphi   = phii(j,:)
   call fftw_execute_dft(iplannlT, tnlT, tnlT)
   call fftw_execute_dft(iplannlphi, tnlphi, tnlphi)
   call fftw_execute_dft(iplanT, tT, tT)
   call fftw_execute_dft(iplanux, tux, tux)
   call fftw_execute_dft(iplanuy, tuy, tuy)
   call fftw_execute_dft(iplanphi, tphi, tphi)
   nlT(j,:)   = tnlT
   nlphi(j,:) = tnlphi
   Ti(j,:)   = tT
   uxi(j,:)  = tux
   uyi(j,:)  = tuy
   phii(j,:) = tphi
end do

! Calculate nonlinear term
do i = 1,Nx
   ! Temperature
   tmp_T = Ti(:,i)
   nlT(:,i) = uxi(:,i)*nlT(:,i) + uyi(:,i)*d1y(tmp_T)
   ! phi
   nlphi(:,i) = uxi(:,i)*phii(:,i) - uyi(:,i)*nlphi(:,i) 
end do

! Bring nonlinear terms back to Fourier space
do j = 1,Ny
   tnlT   = nlT(j,:)
   tnlphi = nlphi(j,:)
   call fftw_execute_dft(plannlT, tnlT, tnlT)
   call fftw_execute_dft(plannlphi, tnlphi, tnlphi)
   ! Dealias
   do i = 1,Nx
      if (abs(kx(i))/alpha >= Nf/2) then
         tnlT(i)   = cmplx(0.0_dp, 0.0_dp, kind=C_DOUBLE_COMPLEX)
         tnlphi(i) = cmplx(0.0_dp, 0.0_dp, kind=C_DOUBLE_COMPLEX)
      end if
   end do
   nlT(j,:)   = tnlT
   nlphi(j,:) = tnlphi
end do
nlT   = nlT   / real(Nx,kind=dp)
nlphi = nlphi / real(Nx,kind=dp)

select case (stage)
   case (1)
      do i = 1,Nx
        !K1hat_phi(:,i) = K1hat_phi(:,i) + CI*kx(i)*nlphi(:,i)
        K1hat_phi(:,i) = K1hat_phi(:,i) - CI*kx(i)*nlphi(:,i)
      end do
      K1hat_T = -nlT
   case (2)
      do i = 1,Nx
        !K2hat_phi(:,i) = K2hat_phi(:,i) + CI*kx(i)*nlphi(:,i)
        K2hat_phi(:,i) = K2hat_phi(:,i) - CI*kx(i)*nlphi(:,i)
      end do
      K2hat_T = -nlT
   case (3)
      do i = 1,Nx
        !K3hat_phi(:,i) = K3hat_phi(:,i) + CI*kx(i)*nlphi(:,i)
        K3hat_phi(:,i) = K3hat_phi(:,i) - CI*kx(i)*nlphi(:,i)
      end do
      K3hat_T = -nlT
   case (4)
      do i = 1,Nx
        !K4hat_phi(:,i) = K4hat_phi(:,i) + CI*kx(i)*nlphi(:,i)
        K4hat_phi(:,i) = K4hat_phi(:,i) - CI*kx(i)*nlphi(:,i)
      end do
      K4hat_T = -nlT
end select

end subroutine calc_explicit

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
subroutine calc_vi(vi, phiin)

complex(C_DOUBLE_COMPLEX),              dimension(:),    intent(in)  :: phiin
complex(C_DOUBLE_COMPLEX), allocatable, dimension(:),    intent(out) :: vi
real(dp),                  allocatable, dimension(:,:)               :: vi_rhs
real(dp),                  allocatable, dimension(:)                 :: dvi, dlvi, duvi
integer                                                              :: j

allocate(vi(Ny), stat=alloc_err)
call check_alloc_err(alloc_err)
allocate(vi_rhs(Ny-2,2), stat=alloc_err)
call check_alloc_err(alloc_err)
allocate(dvi(Ny-2), stat=alloc_err)
call check_alloc_err(alloc_err)
allocate(dlvi(Ny-3), duvi(Ny-3), stat=alloc_err)
call check_alloc_err(alloc_err)

vi     = cmplx(0.0_dp, 0.0_dp, kind=C_DOUBLE_COMPLEX)
vi_rhs = 0.0_dp
dvi    = 0.0_dp
dlvi   = 0.0_dp
duvi   = 0.0_dp

do j = 2,Ny-1
   dvi(j-1) = -kx(it)**2.0_dp + g2(j)
end do

do j = 2,Ny-2
   duvi(j-1) = g3(j)
end do

do j = 3,Ny-1
  dlvi(j-2) = g1(j)
end do

vi_rhs(:,1) = real (phiin(2:Ny-1))
vi_rhs(:,2) = aimag(phiin(2:Ny-1))

call dgtsv(Ny-2, 2, dlvi, dvi, duvi, vi_rhs, Ny-2, info)

vi(1)   = cmplx(0.0_dp, 0.0_dp, kind=C_DOUBLE_COMPLEX)
vi(2:Ny-1) = cmplx(vi_rhs(:,1), vi_rhs(:,2), kind=C_DOUBLE_COMPLEX)
vi(Ny)  = cmplx(0.0_dp, 0.0_dp, kind=C_DOUBLE_COMPLEX)

end subroutine calc_vi

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
subroutine update_bcs(phiout,vout, phiin,vin)

complex(C_DOUBLE_COMPLEX),              dimension(:),  intent(in)  :: phiin, vin
complex(C_DOUBLE_COMPLEX), allocatable, dimension(:),  intent(out) :: phiout, vout
real(dp),                               dimension(2,2)             :: C
real(dp)                                                           :: detC
complex(dp)                                                        :: c1, c2, c1t
complex(dp)                                                        :: dyV_T, dyV_B

allocate(phiout(Ny), vout(Ny), stat=alloc_err)
call check_alloc_err(alloc_err)
phiout = cmplx(0.0_dp, 0.0_dp, kind=C_DOUBLE_COMPLEX)
vout   = cmplx(0.0_dp, 0.0_dp, kind=C_DOUBLE_COMPLEX)
C      = 0.0_dp

C(1,1) = dyv1_T(it)
C(1,2) = dyv2_T(it)
C(2,1) = dyv1_B(it)
C(2,2) = dyv2_B(it)

detC = C(1,1)*C(2,2) - C(1,2)*C(2,1)

dyV_T = h1(Ny)*vin(Ny-2) + h2(Ny)*vin(Ny-1) + h3(Ny)*vin(Ny)
dyV_B = h1(1)*vin(1) + h2(1)*vin(2) + h3(1)*vin(3)

! Need to negate b/c want to solve Cx = -c12.
c1 = -dyV_T
c2 = -dyV_B

! Find c1 and c2.
if (detC == 0.0_dp) then
   c1 = (0.0_dp, 0.0_dp)
   c2 = (0.0_dp, 0.0_dp)
else
   c1t = (C(2,2)*c1 - C(1,2)*c2) / detC
   c2  = (C(1,1)*c2 - C(2,1)*c1) / detC
   c1  = c1t
end if

! Update uy and Phi.
vout   = vin   + c1*V1(:,it)   + c2*V2(:,it)
phiout = phiin + c1*phi1(:,it) + c2*phi2(:,it)

end subroutine update_bcs

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
subroutine update_dt

integer  :: ii, jj
real(dp) :: tmp, tmpmax

uxi = ux
uyi = uy

do jj = 1,Ny
   ! Bring everything to physical space
   tux    = uxi(jj,:)
   tuy    = uyi(jj,:)
   call fftw_execute_dft(iplanux, tux, tux)
   call fftw_execute_dft(iplanuy, tuy, tuy)
   uxi(jj,:)  = tux
   uyi(jj,:)  = tuy
end do

dt_old = dt

do jj = 1,Ny
   do ii = 1,Nx
      tmp = real(uxi(jj,ii)) / dxmin + real(uyi(jj,ii)) / dymin
      if (tmp > tmpmax) then
         tmpmax = tmp
      end if
   end do
end do

dt = cfl / tmpmax

if (dt > dt_ramp * dt_old) then
   dt = dt_ramp * dt_old
else if (dt < dt_old / dt_ramp) then
   dt = dt_old / dt_ramp
end if

if (dt > dtmax) then
   dt = dtmax
else if (dt < dtmin) then
   dt = dtmin
end if

call init_bc(acoeffs(1,1))

!write(*,*) "magu_max, dt:   ", magumax, dt
!write(*,*) "cfl, dxmin, magu_max, dt", cfl, dxmin, magumax, dt
!write(*,*) "cfl, tmpmax, dt", cfl, tmpmax, dt

end subroutine update_dt

end module time_integrators
