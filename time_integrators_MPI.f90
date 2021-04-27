module time_integrators_MPI

use fftw
use global
use write_pack
use allocate_vars
use bc_setup
use statistics
use omp_lib

integer  :: it, jt, kkt
integer  :: info
real(dp) :: time, dtmax, dtmin, dt_old, dt_ramp, dt_final

contains

subroutine imex_rk_MPI(proc_id_str, vtk_print, save_nusselt, proc_id, num_procs)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  This progam solves the equations of thermal convection using a Fourier
!!  spectral method in the x-direction and a 2nd order finite difference scheme in
!!  the y-direction.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


implicit none

character(3), intent(in)       :: proc_id_str
integer, optional, intent(in)  :: vtk_print
logical, optional, intent(in)  :: save_nusselt
integer,           intent(in)  :: proc_id, num_procs

integer                        :: nti, i, j
integer                        :: nprint
logical                        :: wvtk
real(dp)                       :: nusselt_num
real(dp)                       :: start, finish
real(dp)                       :: start_overall, finish_overall
integer                        :: nthreads, myid
integer, EXTERNAL              :: OMP_GET_THREAD_NUM, OMP_GET_NUM_THREADS
real(dp), EXTERNAL             :: OMP_GET_WTIME
EXTERNAL                       :: OMP_SET_NUM_THREADS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call OMP_SET_NUM_THREADS(1)

if (present(vtk_print)) then
    wvtk = .true.
    nprint = vtk_print
else
    wvtk = .false.
    nprint = 100
end if

write(*,*) "imex_rk_MPI from proc ", proc_id_str

if (wvtk) then
    call write_to_vtk(0, .false., proc_id_str) ! false = Fourier space
end if

dt = dt_init

call init_bc_MPI(acoeffs(1,1), proc_id, num_procs, proc_id_str)

if (proc_id == 0) then
    open(unit=9000, file="P000dyv1_b.txt", action="write", status="unknown")
    open(unit=9001, file="P000dyv2_b.txt", action="write", status="unknown")
    
    do i=1,Nx
        write (9000,*) dyv1_B(i)
        write (9001,*) dyv2_B(i)
    end do
    
    close(unit=9000)
    close(unit=9001)
end if

if (proc_id == num_procs-1) then
    open(unit=9000, file="P003dyv1_t.txt", action="write", status="unknown")
    open(unit=9001, file="P003dyv2_t.txt", action="write", status="unknown")
    
    do i=1,Nx
        write (9000,*) dyv1_T(i)
        write (9001,*) dyv2_T(i)
    end do
    
    close(unit=9000)
    close(unit=9001)
end if



end subroutine imex_rk_MPI


end module time_integrators_MPI
    