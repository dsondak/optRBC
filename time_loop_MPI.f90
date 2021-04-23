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
use time_integrators

implicit none
include 'mpif.h'

integer                                                :: ntokens, ios
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
real(dp)                                                :: mpi_spacing_y
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

! Print MPI division.
if (proc_id == 0) then
    write(*,*) "Ny = ", Ny, " divided among ", num_procs, " processors -> ", &
               Ny / num_procs, " rows per processor."
end if

! Initialize MPI variables.
Ny = Ny / 4
mpi_spacing_y = (ytop - ybot) / num_procs
ybot = ybot + proc_id * mpi_spacing_y
ytop = ybot + mpi_spacing_y

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


call cosine_mesh(xp,yp,zp, Nx,Ny,Nz) ! get coordinates
call dxdydz(dynu, xp,yp,zp) ! get mesh spacing for nonuniform grids
call y_mesh_params ! get metric coefficients for nonuniform grids
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

! Write initial field to vtk
if (wvtk) then
    write(proc_id_str, "(I3.3)") proc_id
    call write_to_vtk(int(Ra), .true., proc_id_str) ! true = already in physical space
 end if

write(*,*) "processor ", proc_id, "initialized with ", Ny, "rows."





call MPI_Finalize(mpierror)

end program time_loop

    