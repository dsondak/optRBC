program time_loop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use global
use strings
use write_pack
use interpolation_pack
use allocate_vars
use statistics
use jacobians
use gmres_pack
use mesh_pack
use time_integrators
use nonlinear_solvers

implicit none

integer                                                :: ntokens, ios
integer                                                :: Nxc, Nyc, Nzc
integer                                                :: nRa, ii, jj, n_flm, alphai
integer                                                :: ext
integer, parameter                                     :: alphai_max = 3
integer                                                :: refine_flag_x, refine_flag_y, refine_flag_z
integer                                                :: vtk_flag, rstrt_flag, opt_flag
integer                                                :: nl_iter_data
integer, parameter                                     :: iprint = 1
integer,                                dimension(1)   :: imin
logical                                                :: wvtk, save_restart, calc_opt, physical
logical                                                :: refine_x, refine_y, refine_z, refine_xy
logical                                                :: fTexist, fuyexist
logical                                                :: ex_Tptrb
real(dp)                                               :: temp, eta
real(dp)                                               :: gmres_iter_data
real(dp)                                               :: dxc
real(dp)                                               :: Ra, dRa, Nu, dalpha
real(dp)                                               :: Lmax, Numax, alphamax
real(dp),                               dimension(3)   :: Nui, Li
real(dp),                  allocatable, dimension(:)   :: Nus, Ras, Especu, EspecT
real(dp),                  allocatable, dimension(:)   :: Numaxs, Lmaxs
real(dp),                  allocatable, dimension(:)   :: xpc, ypc, zpc
complex(C_DOUBLE_COMPLEX), allocatable, dimension(:,:) :: uyc, Tc
character(10)                                          :: citer
character(9)                                           :: fiter
character(80)                                          :: line
character(80)                                          :: tokens(80)

open(2,file="input.data", status="old")
do ii = 1,7
   read(2,'(A)') line
   call parse(line,' ,', tokens, ntokens)
   if (ntokens > 0) then
      select case(ii)
        case (1)
          call value(tokens(1), Pr, ios)
          call value(tokens(2), alpha, ios)
          call value(tokens(3), dalpha, ios)
        case (2)
          call value(tokens(1), Ra,  ios)
          call value(tokens(2), nRa, ios)
          call value(tokens(3), dRa, ios)
        case (3)
          call value(tokens(1), t_final, ios)
          call value(tokens(2), dt, ios)
          dt_init = dt
        case (4)
          call value(tokens(1), ybot, ios)
          call value(tokens(2), ytop, ios)
        case (5)
          call value(tokens(1), Nx, ios)
          call value(tokens(2), Ny, ios)
          call value(tokens(3), Nz, ios)
        case (6)
          call value(tokens(1), refine_flag_x, ios)
          call value(tokens(2), refine_flag_y, ios)
          call value(tokens(3), refine_flag_z, ios)
          if (refine_flag_x == 1) then
             if (refine_flag_y == 1) then
                refine_xy = .true.
                refine_x = .false.
                refine_y = .false.
             else
               refine_x = .true.
               refine_y = .false.
               refine_xy = .false.
             end if
          else if (refine_flag_y == 1) then
             refine_x = .false.
             refine_y = .true.
             refine_xy = .false.
          else if (refine_flag_z == 1) then
             write(*,*) "Refinement in z-direction not available yet!"
             stop
          else
             refine_x = .false.
             refine_y = .false.
             refine_z = .false.
             refine_xy = .false.
          end if
        case (7)
          call value(tokens(1), vtk_flag, ios)
          call value(tokens(2), rstrt_flag, ios)
          call value(tokens(2), opt_flag, ios)
          if (vtk_flag == 1) then
             wvtk = .true.
          else
             wvtk = .false.
          end if
          if (rstrt_flag == 1) then
             save_restart = .true.
          else
             save_restart = .false.
          end if
          if (opt_flag == 1) then
             calc_opt = .true.
          else
             calc_opt = .false.
          end if
      end select
   end if
end do
close(unit=2)

! Create FFT plans
planuy = fftw_plan_dft_1d(Nx,tuy,tuy, FFTW_FORWARD,FFTW_ESTIMATE)
iplanuy = fftw_plan_dft_1d(Nx,tuy,tuy, FFTW_BACKWARD,FFTW_ESTIMATE)

planux = fftw_plan_dft_1d(Nx,tux,tux, FFTW_FORWARD,FFTW_ESTIMATE)
iplanux = fftw_plan_dft_1d(Nx,tux,tux, FFTW_BACKWARD,FFTW_ESTIMATE)

planphi = fftw_plan_dft_1d(Nx,tphi,tphi, FFTW_FORWARD,FFTW_ESTIMATE)
iplanphi = fftw_plan_dft_1d(Nx,tphi,tphi, FFTW_BACKWARD,FFTW_ESTIMATE)

planT = fftw_plan_dft_1d(Nx,tT,tT, FFTW_FORWARD,FFTW_ESTIMATE)
iplanT = fftw_plan_dft_1d(Nx,tT,tT, FFTW_BACKWARD,FFTW_ESTIMATE)

plannlT = fftw_plan_dft_1d(Nx,tnlT,tnlT, FFTW_FORWARD,FFTW_ESTIMATE)
iplannlT = fftw_plan_dft_1d(Nx,tnlT,tnlT, FFTW_BACKWARD,FFTW_ESTIMATE)

plannlphi = fftw_plan_dft_1d(Nx,tnlphi,tnlphi, FFTW_FORWARD,FFTW_ESTIMATE)
iplannlphi = fftw_plan_dft_1d(Nx,tnlphi,tnlphi, FFTW_BACKWARD,FFTW_ESTIMATE)

call global_params
call global_allocations

! Allocate some local arrays
allocate(Especu(kmax), EspecT(kmax), stat=alloc_err)
allocate(Nus(nRa), Ras(nRa), stat=alloc_err)
allocate(Numaxs(nRa), Lmaxs(nRa), stat=alloc_err)
call check_alloc_err(alloc_err)

Nus   = 0.0_dp
Ras   = 0.0_dp
Especu = 0.0_dp
EspecT = 0.0_dp

inquire(file="uy", exist=fuyexist)
if (fuyexist) then
   open(unit=2,file="uy",action="read", status="old", form="unformatted")
   if (refine_x) then
      Nxc =  Nx / 2
      allocate(uyc(Ny,Nxc), stat=alloc_err)
      call check_alloc_err(alloc_err)
      do ii = 1,Nxc
         do jj = 1,Ny
            read(2) temp
            uyc(jj,ii) = cmplx(temp, 0.0_dp, kind=C_DOUBLE_COMPLEX)
         end do
      end do
   else if (refine_y) then
      Nyc = (Ny+1) / 2
      allocate(uyc(Nyc,Nx), stat=alloc_err)
      call check_alloc_err(alloc_err)
      do ii = 1,Nx
         do jj = 1,Nyc
            read(2) temp
            uyc(jj,ii) = cmplx(temp, 0.0_dp, kind=C_DOUBLE_COMPLEX)
         end do
      end do
   else if (refine_xy) then
      Nxc =  Nx / 2
      Nyc = (Ny+1) / 2
      Nzc =  Nz
      allocate(uyc(Nyc,Nxc), stat=alloc_err)
      call check_alloc_err(alloc_err)
      do ii = 1,Nxc
         do jj = 1,Nyc
            read(2) temp
            uyc(jj,ii) = cmplx(temp, 0.0_dp, kind=C_DOUBLE_COMPLEX)
         end do
      end do
   else
      do ii = 1,Nx
         do jj = 1,Ny
            read(2) temp
            uy(jj,ii) = cmplx(temp, 0.0_dp, kind=C_DOUBLE_COMPLEX)
         end do
      end do
   end if
   close(unit=2)
end if

inquire(file="temperature", exist=fTexist)
if (fTexist) then
   open(unit=3,file="temperature",action="read", status="old", form="unformatted")
   if (refine_x) then
      Nxc =  Nx / 2
      allocate(Tc(Ny,Nxc), stat=alloc_err)
      call check_alloc_err(alloc_err)
      do ii = 1,Nxc
         do jj = 1,Ny
            read(3) temp
            Tc(jj,ii) = cmplx(temp, 0.0_dp, kind=C_DOUBLE_COMPLEX)
         end do
      end do
   else if (refine_y) then
      Nyc = (Ny+1) / 2
      allocate(Tc(Nyc,Nx), stat=alloc_err)
      call check_alloc_err(alloc_err)
      do ii = 1,Nx
         do jj = 1,Nyc
            read(3) temp
            Tc(jj,ii) = cmplx(temp, 0.0_dp, kind=C_DOUBLE_COMPLEX)
         end do
      end do
   else if (refine_xy) then
      Nxc =  Nx / 2
      Nyc = (Ny+1) / 2
      Nzc =  Nz
      allocate(Tc(Nyc,Nxc), stat=alloc_err)
      call check_alloc_err(alloc_err)
      do ii = 1,Nxc
         do jj = 1,Nyc
            read(3) temp
            Tc(jj,ii) = cmplx(temp, 0.0_dp, kind=C_DOUBLE_COMPLEX)
         end do
      end do
   else
      do ii = 1,Nx
         do jj = 1,Ny
            read(3) temp
            T(jj,ii) = cmplx(temp, 0.0_dp, kind=C_DOUBLE_COMPLEX)
         end do
      end do
   end if
   close(unit=3)
end if

xR =  pi / alpha
xL = -pi / alpha
Lx = xR - xL

if (Nx == 1) then
   dx = 1.0_dp
else
   dx = Lx / (real(Nx,kind=dp))
end if

! Create the grid
call cosine_mesh(xp,yp,zp, Nx,Ny,Nz) ! get coordinates
call dxdydz(dynu, xp,yp,zp) ! get mesh spacing for nonuniform grids
call y_mesh_params ! get metric coefficients for nonuniform grids
dymin = minval(dynu)
dxmin = sqrt(dx**2.0_dp + dymin**2.0_dp)

if (refine_xy) then
   if (fuyexist) then
      if (fTexist) then
         allocate(xpc(Nxc), stat=alloc_err)
         allocate(ypc(Nyc), stat=alloc_err)
         allocate(zpc(Nzc), stat=alloc_err)
         call check_alloc_err(alloc_err)
         dxc = Lx / real(Nxc,kind=dp)
         call cosine_mesh(xpc,ypc,zpc, Nxc,Nyc,Nzc, dxc)
         call bilinear_interp(xpc,ypc,uyc,Tc)
      else if (.not. fTexist) then
         write(*,*) "ERROR:  Need input data from uy and temperature.  Only have data from uy."
         stop
      end if
   else if (.not. fuyexist) then
      write(*,*) "ERROR:  Need input data from uy and temperature.  Not able to find uy data."
      stop
   end if
else if (refine_x) then
   if (fuyexist) then
      if (fTexist) then
         allocate(xpc(Nxc), stat=alloc_err)
         allocate(ypc(Ny), stat=alloc_err)
         allocate(zpc(Nz), stat=alloc_err)
         call check_alloc_err(alloc_err)
         dxc = Lx / real(Nxc,kind=dp)
         call cosine_mesh(xpc,ypc,zpc, Nxc,Ny,Nz, dxc)
         call linear_interp_x(xpc,uyc,Tc)
      else if (.not. fTexist) then
         write(*,*) "ERROR:  Need input data from uy and temperature.  Only have data from uy."
         stop
      end if
   else if (.not. fuyexist) then
      write(*,*) "ERROR:  Need input data from uy and temperature.  Not able to find uy data."
      stop
   end if
else if (refine_y) then
   if (fuyexist) then
      if (fTexist) then
         allocate(xpc(Nx), stat=alloc_err)
         allocate(ypc(Nyc), stat=alloc_err)
         allocate(zpc(Nz), stat=alloc_err)
         call check_alloc_err(alloc_err)
         call cosine_mesh(xpc,ypc,zpc, Nx,Nyc,Nz, dx)
         call linear_interp_y(ypc,uyc,Tc)
      else if (.not. fTexist) then
         write(*,*) "ERROR:  Need input data from uy and temperature.  Only have data from uy."
         stop
      end if
   else if (.not. fuyexist) then
      write(*,*) "ERROR:  Need input data from uy and temperature.  Not able to find uy data."
      stop
   end if
end if

do ii = 1,Nx/2
   kx_modes(ii) = real(ii,kind=dp) - 1.0_dp
end do
do ii = Nx/2+1, Nx
   kx_modes(ii) = real(ii-Nx,kind=dp) - 1.0_dp
end do
kx = alpha*kx_modes

! Write initial field to vtk
if (wvtk) then
   call write_to_vtk(int(Ra), .true.) ! true = already in physical space
end if

! Initialize fields.
call init_fields(ex_Tptrb, fTexist, Ra)
call init_to_fourier(ex_Tptrb)

! Allocate space for flow map variables
n_flm = 2*2*Ny*Nx ! # of fields X real/imag X Ny X Nx

! Form x0 and GT (flow map)
allocate(x0(n_flm), GT(n_flm), stat=alloc_err)
if (alloc_err /= 0) then
   write(*,*) "ERROR:  Allocation of x0 and GT incomplete."
   stop
end if
x0      = 0.0_dp
GT      = 0.0_dp

! p points to the boundaries.  Will be used to ensure that
! Delta_x is zero on the boundaries.
b_pntr(1)  = 1
b_pntr(2)  = Ny
b_pntr(3)  = Ny+1
b_pntr(4)  = 2*Ny
b_pntr(5)  = 2*Ny+1
b_pntr(6)  = 3*Ny
b_pntr(7)  = 3*Ny+1
b_pntr(8)  = 4*Ny

write(*,'(A70)')                  ' '
write(*,'(A70)')                  '*********************************************************************|'
write(*,'(A70)')                  '                TWO-DIMENSIONAL THERMAL CONVECTION                   |'
write(*,'(A70)')                  '*********************************************************************|'
write(*,'(5X,A,47X,A)')           'COMPUTATION SIZE:',                                                 '|'
write(*,'(A69,A)')                '                                                                  ','|'
write(*,'(20X,A23,I5,21X,A)')     'Nx                       = ', Nx,                                   '|'
write(*,'(20X,A23,I5,21X,A)')     'Ny                       = ', Ny,                                   '|'
write(*,'(20X,A23,I5,21X,A)')     'Nz                       = ', Nz,                                   '|'
write(*,'(20X,A23,I5,21X,A)')     'Number of time steps     = ', nt,                                   '|'
write(*,'(20X,A23,I5,21X,A)')     'Number of Fourier modes  = ', Nf,                                   '|'
write(*,'(A69,A)')                '                                                                  ','|'
write(*,'(A70)')                  '*********************************************************************|'
write(*,'(5X,A,45X,A)')           'PROBLEM PARAMETERS:',                                               '|'
write(*,'(A69,A)')                '                                                                  ','|'
write(*,'(10X,A32,ES16.8,11X,A)') 'Prandtl number  (Pr)                    = ', Pr,                    '|'
write(*,'(10X,A32,ES16.8,11X,A)') 'Initial Rayleigh number (Ra)            = ', Ra,                    '|'
write(*,'(10X,A32,ES16.8,11X,A)') 'Initial Reynolds number (Re)            = ', sqrt(Ra/(16.0_dp*Pr)), '|'
write(*,'(10X,A32,ES16.8,11X,A)') 'Ra Gradation (dRa)                      = ', dRa,                   '|'
write(*,'(10X,A32,I5,    21X,A)') 'Number of steps (nRa)                   = ', nRa,                   '|'
write(*,'(10X,A32,ES16.8,11X,A)') 'Final Rayleigh Number                   = ', (dRa**(nRa-1))*Ra,     '|'
write(*,'(A69,A)')                '                                                                  ','|'
write(*,'(A70)')                  '*********************************************************************|'
write(*,'(5X,A,42X,A)')           'PHYSICAL PROBLEM SIZE:',                                            '|'
write(*,'(A69,A)')                '                                                                  ','|'
write(*,'(10X,A32,ES16.8,11X,A)') 'alpha                     = ', alpha,                               '|'
write(*,'(10X,A32,ES16.8,11X,A)') 'Left coordinate of box    = ', xL,                                  '|'
write(*,'(10X,A32,ES16.8,11X,A)') 'Right coordinate of box   = ', xR,                                  '|'
write(*,'(10X,A32,ES16.8,11X,A)') 'Coordinate of bottom wall = ', ybot,                                '|'
write(*,'(10X,A32,ES16.8,11X,A)') 'Coordinate of top wall    = ', ytop,                                '|'
write(*,'(10X,A32,ES16.8,11X,A)') 'Box width                 = ', Lx,                                  '|'
write(*,'(10X,A32,ES16.8,11X,A)') 'Box height                = ', Ly,                                  '|'
write(*,'(10X,A32,ES16.8,11X,A)') 'Aspect Ratio              = ', Lx/Ly,                               '|'
write(*,'(A69,A)')                '                                                                  ','|'
write(*,'(10X,A32,ES16.8,11X,A)') 'Integration time (T)      = ', t_final,                              '|'
write(*,'(10X,A32,ES16.8,11X,A)') 'Time step size (Delta t)  = ', dt,                                   '|'
write(*,'(A69,A)')                '                                                                  ','|'
write(*,'(A70)')                  '*********************************************************************|'
write(*,'(A70)')                  '*********************************************************************|'
write(*,'(A70)')                  '                                                                      '

flush(6)

open(unit=8000, file="Nu_data.txt", action="write", status="unknown", position="append")
Raloop: do ii = 1, nRa

            Ras(ii) = Ra
            Ra_global = Ra

            if (mod((ii-1),iprint) == 0) then
               if (ii == 1) then
                  write(*,*) " "
                  write(*,*) "Ra, dt = ", Ra, dt
                  write(*,*) " "
               else
                  write(*,*) " "
                  write(*,*)"::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"
                  write(*,*) "Ra, dt = ", Ra, dt
                  write(*,*) " "
               end if
            end if

            ! Get nu0 and kappa0
            call global_params_Ra(Ra)

            ! Populate x0
            call packx(x0)

            ! Get solution with flow map
            !call flow_map(Nu, nl_iter_data, temp)
            call imex_rk(1)

            write(*,*) " "
            flush(6)

!            if (save_restart) then
!               call write_restart(.false.) ! false = Fourier space
!            end if
!            
!            if (wvtk) then
!               call write_to_vtk(int(Ra), .false.) ! false = Fourier space
!            end if
!
!            ! Increment Ra
!            Ra = dRa * Ra

        end do Raloop
        close(unit=8000)

!if (save_restart) then
!   call write_restart(.false.) ! false = Fourier space
!end if
!
!if (wvtk) then
!   call write_to_vtk(int(Ra), .false.) ! false = Fourier space
!end if

write(*,*) " "
write(*,*) "Done."

1000 format(E25.16E3                           )
2000 format(E25.16E3,E25.16E3                  )
3000 format(E25.16E3,E25.16E3,E25.16E3         )
4000 format(E25.16E3,E25.16E3,E25.16E3,E25.16E3)

end program time_loop

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
subroutine init_fields(ex_Tptrb, fTexist,Ra)

use global

implicit none

logical, intent(in)  :: fTexist
logical, intent(out) :: ex_Tptrb
real(dp), intent(in) :: Ra
integer              :: ii

! Initialize fields.
if (.not. fTexist) then
   if (Ra < 1710.0_dp) then
      ex_Tptrb = .true.
      do ii = 1,Nx
         Tptrb(:,ii) = 0.5_dp*2.0_dp*cos(alpha*xp(ii))*cos(pi*yp/2.0_dp)
         T(:,ii) = -yp + Tptrb(:,ii)
      end do
   else
      ex_Tptrb = .false.
      do ii = 1,Nx
         T(:,ii) = -yp + 0.5_dp*2.0_dp*cos(alpha*xp(ii))*cos(pi*yp/2.0_dp)
      end do
   end if
else if (fTexist) then   
   if (Ra < 1710.0_dp) then
      ex_Tptrb = .true.
      do ii = 1,Nx
         Tptrb(:,ii) = 0.5_dp*2.0_dp*cos(alpha*xp(ii))*cos(pi*yp/2.0_dp)
      end do
   else
      ex_Tptrb = .false.
   end if
end if

end subroutine init_fields
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
subroutine init_to_fourier(ex_Tptrb)

use global

implicit none

logical, intent(in) :: ex_Tptrb
integer             :: ii, jj

! Bring temperature and velocity to Fourier space.
do jj = 1,Ny
   tT = T(jj,:)
   tuy = uy(jj,:)
   call fftw_execute_dft(planT, tT, tT)
   call fftw_execute_dft(planuy, tuy, tuy)
   ! Truncate modes
   do ii = 1,Nx
      if (abs(kx(ii))/alpha >= Nf/2) then
         tT(ii)  = cmplx(0.0_dp, 0.0_dp, kind=C_DOUBLE_COMPLEX)
         tuy(ii) = cmplx(0.0_dp, 0.0_dp, kind=C_DOUBLE_COMPLEX)
      end if
   end do
   T(jj,:) = tT
   uy(jj,:) = tuy
   ! If temperature perturbation needed.
   if (ex_Tptrb) then
      tT = Tptrb(jj,:)
      call fftw_execute_dft(planT, tT, tT)
      ! Truncate modes
      do ii = 1,Nx
         if (abs(kx(ii))/alpha >= Nf/2) then
            tT(ii)  = cmplx(0.0_dp, 0.0_dp, kind=C_DOUBLE_COMPLEX)
         end if
      end do
      Tptrb(jj,:) = tT
   end if
end do
T = T / real(Nx, kind=dp)
uy = uy / real(Nx, kind=dp)

if (ex_Tptrb) then
   Tptrb = Tptrb / real(Nx, kind=dp)
end if

! Calculate phi and ux from uy
do ii = 1,Nx
   if (kx(ii) /= 0.0_dp) then
      tmp_uy = uy(:,ii)
      !ux(:,ii) = -CI*d1y(tmp_uy)/kx(ii)
      ux(:,ii) = CI*d1y(tmp_uy)/kx(ii)
   else if (kx(ii) == 0.0_dp) then
      ux(:,ii) = cmplx(0.0_dp, 0.0_dp, kind=C_DOUBLE_COMPLEX) ! Zero mean flow!
   end if
   phi(:,ii) = -kx(ii)**2.0_dp*uy(:,ii) + d2y(uy(:,ii))
end do

end subroutine init_to_fourier
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
subroutine fit_parab(xstar,ystar, x, y)

use global
use write_pack

implicit none

real(dp), dimension(3), intent(in)  :: x, y
real(dp), intent(out)               :: xstar, ystar
real(dp)                            :: den1, den2, den3
real(dp)                            :: c1, c2, c3
real(dp)                            :: x1,x2,x3,y1,y2,y3

x1 = x(1)
x2 = x(2)
x3 = x(3)

y1 = y(1)
y2 = y(2)
y3 = y(3)

den1 = (x1 - x2) * (x1 - x3)
den2 = (x2 - x1) * (x2 - x3)
den3 = (x3 - x1) * (x3 - x2)

c1 = y1 / den1 + y2 / den2 + y3 / den3

c2 = y1 * (x2 + x3) / den1 + &
&    y2 * (x1 + x3) / den2 + &
&    y3 * (x1 + x2) / den3

c3 = y1 * x2 * x3 / den1 + &
&    y2 * x1 * x3 / den2 + &
&    y3 * x1 * x2 / den3

xstar = c2 / (2.0_dp * c1)
ystar = c1 * xstar**2.0_dp - c2 * xstar + c3

end subroutine fit_parab

! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

subroutine calc_optimal(ave_num_gmres, alphain)

use global
use mesh_pack
use nonlinear_solvers

implicit none

integer               :: jj
integer               :: num_nonlinear
real(dp)              :: temp
real(dp), intent(in)  :: alphain
real(dp), intent(out) :: ave_num_gmres

alpha = alphain

xR =  pi / alpha
xL = -pi / alpha
Lx = xR - xL

if (Nx == 1) then
   dx = 1.0_dp
else
   dx = Lx / (real(Nx,kind=dp))
end if

kx = alpha * kx_modes

call cosine_mesh(xp,yp,zp, Nx,Ny,Nz) ! get coordinates

! Calculate phi and ux from uy
do jj = 1,Nx
   if (kx(jj) /= 0.0_dp) then
      tmp_uy = uy(:,jj)
      !ux(:,jj) = -CI*d1y(tmp_uy)/kx(jj)
      ux(:,jj) = CI*d1y(tmp_uy)/kx(jj)
   else if (kx(jj) == 0.0_dp) then
      ux(:,jj) = cmplx(0.0_dp, 0.0_dp, kind=C_DOUBLE_COMPLEX) ! Zero mean flow!
   end if
      phi(:,jj) = -kx(jj)**2.0_dp*uy(:,jj) + d2y(uy(:,jj))
end do

! Populate x0
call packx(x0)

! Get solution with flow map
call flow_map(temp, num_nonlinear, ave_num_gmres)

end subroutine calc_optimal
