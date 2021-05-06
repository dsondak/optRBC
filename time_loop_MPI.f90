program time_loop

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use global
use strings
use write_pack
use interpolation_pack
use allocate_vars
use statistics
use mesh_pack
use time_integrators_MPI

! include 'mpif.h'
implicit none

integer                                                :: ntokens, ios, i
integer                                                :: Nxc, Nyc, Nzc
integer                                                :: nRa, ii, jj
integer                                                :: refine_flag_x, refine_flag_y, refine_flag_z
integer                                                :: vtk_flag, rstrt_flag, opt_flag
logical                                                :: wvtk, save_restart, calc_opt
logical                                                :: refine_x, refine_y, refine_z, refine_xy
logical                                                :: fTexist, fuyexist
logical                                                :: ex_Tptrb
real(dp)                                               :: temp
real(dp)                                               :: dxc
real(dp)                                               :: Ra, dRa, Nu, dalpha
real(dp),                  allocatable, dimension(:)   :: xpc, ypc, zpc
complex(C_DOUBLE_COMPLEX), allocatable, dimension(:,:) :: uyc, Tc
character(10)                                          :: citer
character(9)                                           :: fiter
character(80)                                          :: line
character(80)                                          :: tokens(80)

! MPI specific variables
integer                                                :: mpierror, num_procs
integer                                                :: proc_id
character(4)                                           :: proc_id_str
real(dp)                                               :: mpi_spacing_y
integer                                                status(MPI_STATUS_SIZE)


! Initialize MPI.
call MPI_Init(mpierror)
call MPI_Comm_size(MPI_COMM_WORLD, num_procs, mpierror)
call MPI_Comm_rank(MPI_COMM_WORLD, proc_id, mpierror)

! Read file.
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

! Initialize MPI variables.
Ny = Ny / 4
mpi_spacing_y = (ytop - ybot) / num_procs
ybot = ybot + proc_id * mpi_spacing_y
ytop = ybot + mpi_spacing_y
write(proc_id_str, "(I3.3)") proc_id

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

xR =  pi / alpha
xL = -pi / alpha
Lx = xR - xL

if (Nx == 1) then
    dx = 1.0_dp
else
    dx = Lx / (real(Nx,kind=dp))
end if

call cosine_mesh_MPI(xp,yp,zp, Nx,Ny,Nz, proc_id, num_procs) ! get coordinates
call dxdydz_MPI(dynu, xp,yp,zp, proc_id, num_procs) ! get mesh spacing for nonuniform grids
call y_mesh_params_MPI(proc_id, num_procs) ! get metric coefficients for nonuniform grids

dymin = minval(dynu)
dxmin = sqrt(dx**2.0_dp + dymin**2.0_dp)

! Fourier modes
do ii = 1,Nx/2
    kx_modes(ii) = real(ii,kind=dp) - 1.0_dp
end do
do ii = Nx/2+1, Nx
    kx_modes(ii) = real(ii-Nx,kind=dp) - 1.0_dp
end do
kx = alpha*kx_modes

! Initialize fields.
call init_fields(ex_Tptrb, Ra)
call init_to_fourier(ex_Tptrb)

! Get nu0 and kappa0
call global_params_Ra(Ra)

if (proc_id == 0) then
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
    ! Print MPI division.
    write(*,*) "Ny = ", Ny*num_procs, " divided among ", num_procs, " processors -> ", &
               Ny, " rows per processor."
end if

call MPI_BARRIER(MPI_COMM_WORLD, mpierror)

write(*,*) "processor ", proc_id, "initialized with ", Ny, "rows."

! Get solution with time integration
call imex_rk_MPI(proc_id_str, wvtk, .true., proc_id, num_procs) ! true causes writing of nusselt number.

call MPI_Finalize(mpierror)

end program time_loop

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
subroutine init_fields(ex_Tptrb,Ra)

use global

implicit none

logical, intent(out) :: ex_Tptrb
real(dp), intent(in) :: Ra
integer              :: ii

! Initialize fields.
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