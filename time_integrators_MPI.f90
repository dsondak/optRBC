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
logical,           intent(in)  :: vtk_print
logical,           intent(in)  :: save_nusselt
integer,           intent(in)  :: proc_id, num_procs

integer                        :: nti, i, j, mpierror, total_ny
integer                        :: nprint, otherproc
logical                        :: wvtk
real(dp)                       :: nusselt_num
real(dp)                       :: start, finish
real(dp)                       :: start_overall, finish_overall
integer                        :: nthreads, myid, stind
complex(C_DOUBLE_COMPLEX), allocatable, dimension(:) :: tmp_phi_MPI, tmp_T_MPI, phi_MPI, phi_update_MPI
complex(C_DOUBLE_COMPLEX), allocatable, dimension(:) :: K1hat_phi_MPI, K2hat_phi_MPI, K3hat_phi_MPI
complex(C_DOUBLE_COMPLEX), allocatable, dimension(:) :: K1hat_T_MPI, K2hat_T_MPI, K3hat_T_MPI
complex(C_DOUBLE_COMPLEX), allocatable, dimension(:) :: K1_phi_MPI, K2_phi_MPI, K1_T_MPI, K2_T_MPI
complex(C_DOUBLE_COMPLEX), allocatable, dimension(:) :: T_MPI, tmp_uy_MPI, tmp_uy1_MPI, tmp_phi1_MPI
complex(C_DOUBLE_COMPLEX), allocatable, dimension(:) :: tmp_K_phi_MPI, tmp_K_T_MPI, uxi_MPI
real(dp), allocatable, dimension(:) :: phi1_MPI, phi2_MPI, V1_MPI, V2_MPI
complex(C_DOUBLE_COMPLEX), allocatable, dimension(:,:) :: single_phiout, single_Tout
real(dp), allocatable, dimension(:)   :: g1_total, g2_total, g3_total, dynu_MPI
real(dp), allocatable, dimension(:)   :: h1_total, h2_total, h3_total
real(dp)                       :: dyv1_T_it_MPI, dyv2_T_it_MPI, h1_end_MPI, h2_end_MPI, h3_end_MPI
integer, EXTERNAL              :: OMP_GET_THREAD_NUM, OMP_GET_NUM_THREADS
real(dp), EXTERNAL             :: OMP_GET_WTIME
EXTERNAL                       :: OMP_SET_NUM_THREADS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call OMP_SET_NUM_THREADS(1)

if (vtk_print) then
    wvtk = .true.
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

time = 0.0_dp

dtmax = 0.5_dp
dtmin = 1.0e-4_dp

dt_ramp = 1.1_dp

dt_old = dt

nti = 0

! Format for writing out single values.x
1000 format(E25.16E3)

if (proc_id == 0) then 
   total_ny = Ny * num_procs
   allocate(tmp_phi_MPI(total_ny), tmp_T_MPI(total_ny), phi_MPI(total_ny-2), phi_update_MPI(total_ny), stat=alloc_err)
   call check_alloc_err(alloc_err)
   allocate(K1hat_phi_MPI(total_ny-2), K2hat_phi_MPI(total_ny-2), K3hat_phi_MPI(total_ny-2), stat=alloc_err)
   call check_alloc_err(alloc_err)
   allocate(K1hat_T_MPI(total_ny-2), K2hat_T_MPI(total_ny-2), K3hat_T_MPI(total_ny-2), stat=alloc_err)
   call check_alloc_err(alloc_err)
   allocate(K1_phi_MPI(total_ny-2), K2_phi_MPI(total_ny-2), K1_T_MPI(total_ny-2), K2_T_MPI(total_ny-2), stat=alloc_err)
   call check_alloc_err(alloc_err)
   allocate(T_MPI(total_ny), stat=alloc_err)
   call check_alloc_err(alloc_err)
   allocate(g1_total(total_ny), g2_total(total_ny), g3_total(total_ny), dynu_MPI(total_ny-1), stat=alloc_err)
   call check_alloc_err(alloc_err)
   allocate(h1_total(total_ny), h2_total(total_ny), h3_total(total_ny), stat=alloc_err)
   call check_alloc_err(alloc_err)
   allocate(tmp_uy_MPI(total_ny), tmp_uy1_MPI(total_ny), tmp_phi1_MPI(total_ny), stat=alloc_err)
   call check_alloc_err(alloc_err)
   allocate(single_phiout(total_ny, Nx), single_Tout(total_ny,Nx), stat=alloc_err)
   call check_alloc_err(alloc_err)
   allocate(phi1_MPI(total_ny), phi2_MPI(total_ny), V1_MPI(total_ny), V2_MPI(total_ny), stat=alloc_err)
   call check_alloc_err(alloc_err)
   allocate(tmp_K_phi_MPI(total_ny), tmp_K_T_MPI(total_ny), uxi_MPI(total_ny), stat=alloc_err)
   call check_alloc_err(alloc_err)

    ! Receive g1,g2,g3 from all other nodes.
    do otherproc = 1,num_procs-1
        stind = otherproc * Ny + 1
        call MPI_RECV(g1_total(stind), Ny, MPI_DOUBLE, otherproc, 54, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
        call MPI_RECV(g2_total(stind), Ny, MPI_DOUBLE, otherproc, 55, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
        call MPI_RECV(g3_total(stind), Ny, MPI_DOUBLE, otherproc, 56, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
    end do
    g1_total(1:Ny) = g1
    g2_total(1:Ny) = g2
    g3_total(1:Ny) = g3

    ! Receive h1,h2,h3 from all other nodes.
    do otherproc = 1,num_procs-1
        stind = otherproc * Ny + 1
        call MPI_RECV(h1_total(stind), Ny, MPI_DOUBLE, otherproc, 67, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
        call MPI_RECV(h2_total(stind), Ny, MPI_DOUBLE, otherproc, 68, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
        call MPI_RECV(h3_total(stind), Ny, MPI_DOUBLE, otherproc, 69, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
    end do
    h1_total(1:Ny) = h1
    h2_total(1:Ny) = h2
    h3_total(1:Ny) = h3

    ! Receive dynu from all other nodes.
    do otherproc = 1,num_procs-1
        stind = otherproc * Ny
        call MPI_RECV(dynu_MPI(stind), Ny, MPI_DOUBLE, otherproc, 66, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
    end do
    dynu_MPI(1:Ny-1) = dynu(1:Ny-1)
else 
    ! Send g1, g2, g3 to first node.
    call MPI_SEND(g1(1), Ny, MPI_DOUBLE, 0, 54, MPI_COMM_WORLD, mpierror)
    call MPI_SEND(g2(1), Ny, MPI_DOUBLE, 0, 55, MPI_COMM_WORLD, mpierror)
    call MPI_SEND(g3(1), Ny, MPI_DOUBLE, 0, 56, MPI_COMM_WORLD, mpierror)

    ! Send h1, h2, h3 to first node.
    call MPI_SEND(h1(1), Ny, MPI_DOUBLE, 0, 67, MPI_COMM_WORLD, mpierror)
    call MPI_SEND(h2(1), Ny, MPI_DOUBLE, 0, 68, MPI_COMM_WORLD, mpierror)
    call MPI_SEND(h3(1), Ny, MPI_DOUBLE, 0, 69, MPI_COMM_WORLD, mpierror)

    ! Send dynu to main.
    call MPI_SEND(dynu(1), Ny, MPI_DOUBLE, 0, 66, MPI_COMM_WORLD, mpierror)

end if 

call MPI_BARRIER(MPI_COMM_WORLD, mpierror) 

! Time integration
do ! while (time < t_final)
    start_overall = OMP_GET_WTIME()
 
    dt_final = t_final - time
 
    if (dt_final <= dt) then
       time = t_final
    else
       time = time + dt
    end if
 
    write(*,*) "time = ", time, "dt = ", dt
 
    nti = nti + 1
 
    !:::::::::::
    ! STAGE 1 ::
    !:::::::::::
    phii = phi
    Ti   = T
    uxi  = ux
    uyi  = uy
    start = OMP_GET_WTIME()
    call calc_explicit_MPI(1, proc_id, num_procs, proc_id_str)
    finish = OMP_GET_WTIME()
    ! write(*,*) " - calc_explicit(1) timing: ", finish-start, "(s)"
    start = OMP_GET_WTIME()
    !$OMP PARALLEL DO private(tmp_phi, tmp_T, tmp_uy, tmp_phi1, tmp_uy1, tmp_K_phi, tmp_K_T) schedule(dynamic)
    do it = 1,Nx ! kx loop
        ! Compute phi1 and T1
        if (proc_id == 0) then
            ! Need to pull all of the variables into local main memory.
            phi_MPI(1:Ny-1) = phi(2:Ny,it)
            K1hat_phi_MPI(1:Ny-1) = K1hat_phi(2:Ny,it)
            K1hat_T_MPI(1:Ny-1) = K1hat_T(2:Ny,it)
            T_MPI(1:Ny) = T(1:Ny,it)

            ! Receive pieces from other nodes.
            do otherproc = 1,num_procs-2
                stind = otherproc * Ny
                call MPI_RECV(phi_MPI(stind), Ny, MPI_C_DOUBLE_COMPLEX, &
                              otherproc, 42, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
                call MPI_RECV(K1hat_phi_MPI(stind), Ny, MPI_C_DOUBLE_COMPLEX, &
                              otherproc, 43, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
                call MPI_RECV(K1hat_T_MPI(stind), Ny, MPI_C_DOUBLE_COMPLEX, & 
                              otherproc, 46, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
                call MPI_RECV(T_MPI(stind+1), Ny, MPI_C_DOUBLE_COMPLEX, & 
                              otherproc, 53, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
            end do
            
            ! Receive from last node.
            stind = (num_procs-1) * Ny
            call MPI_RECV(phi_MPI(stind), Ny-1, MPI_C_DOUBLE_COMPLEX, & 
                          num_procs-1, 42, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
            call MPI_RECV(K1hat_phi_MPI(stind), Ny-1, MPI_C_DOUBLE_COMPLEX, & 
                          num_procs-1, 43, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
            call MPI_RECV(K1hat_T_MPI(stind), Ny-1, MPI_C_DOUBLE_COMPLEX, & 
                          num_procs-1, 46, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
            call MPI_RECV(T_MPI(stind+1), Ny, MPI_C_DOUBLE_COMPLEX, & 
                          num_procs-1, 53, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)


            ! Call vari mod
            call calc_vari_mod_MPI(tmp_phi_MPI, tmp_T_MPI, acoeffs(1,1), 1, total_ny,&
                                kx(it), phi_MPI,&
                                K1hat_phi_MPI,K2hat_phi_MPI,K3hat_phi_MPI,&
                                K1hat_T_MPI,K2hat_T_MPI,K3hat_T_MPI,&
                                K1_phi_MPI, K2_phi_MPI, K1_T_MPI, K2_T_MPI,&
                                g1_total, g2_total, g3_total,&
                                T_MPI)
            
            ! Compute v1 from phi1
            call calc_vi_mod_MPI(tmp_uy_MPI, tmp_phi_MPI, kx(it), total_ny, g1_total, g2_total, g3_total)
            
            ! Receive dyv1 and dyv2
            call MPI_RECV(dyv1_T_it_MPI, 1, MPI_DOUBLE, num_procs-1, 57, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
            call MPI_RECV(dyv2_T_it_MPI, 1, MPI_DOUBLE, num_procs-1, 58, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)

            ! Receive V1, V2, phi1, phi2
            do otherproc = 1,num_procs-1
                stind = otherproc * Ny + 1
                call MPI_RECV(V1_MPI(stind), Ny, MPI_DOUBLE, otherproc, 59, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
                call MPI_RECV(V2_MPI(stind), Ny, MPI_DOUBLE, otherproc, 60, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
                call MPI_RECV(phi1_MPI(stind), Ny, MPI_DOUBLE, otherproc, 61, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
                call MPI_RECV(phi2_MPI(stind), Ny, MPI_DOUBLE, otherproc, 62, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
            end do
            ! Set local variables to MPI.
            V1_MPI(1:Ny) = V1(1:Ny, it)
            V2_MPI(1:Ny) = V2(1:Ny, it)
            phi1_MPI(1:Ny) = phi1(1:Ny, it)
            phi2_MPI(1:Ny) = phi2(1:Ny, it)

            ! BOUNDAY CONDITIONS!
            call update_bcs_mod_MPI(tmp_phi1_MPI,tmp_uy1_MPI, tmp_phi_MPI,tmp_uy_MPI,&
                                dyv1_T_it_MPI,dyv2_T_it_MPI,&
                                dyv1_B(it),dyv2_B(it),&
                                V1_MPI,V2_MPI,phi1_MPI,phi2_MPI,&
                                total_ny, h1_total(total_ny), h2_total(total_ny), h3_total(total_ny))
            
            tmp_phi_MPI = tmp_phi1_MPI
            tmp_uy_MPI  = tmp_uy1_MPI
            
            ! Implicit.
            call calc_implicit_mod_MPI(tmp_K_phi_MPI,tmp_K_T_MPI, tmp_phi_MPI,tmp_T_MPI,&
                                       kx(it), total_ny, dynu_MPI, g1_total, g2_total, g3_total)
            
            ! Compute u1 from v1
            if (kx(it) /= 0.0_dp) then
                uxi_MPI = CI*d1y_MPI(tmp_uy_MPI, h1_total, h2_total, h3_total)/kx(it)
            else if (kx(it) == 0.0_dp) then
                uxi_MPI = cmplx(0.0_dp, 0.0_dp, kind=C_DOUBLE_COMPLEX) ! Zero mean flow!
            end if

            ! Send K1_phi, K1_T, phii, Ti, uyi, uxi back to each node.
            do otherproc = 1,num_procs-1
                stind = otherproc * Ny + 1
                call MPI_SEND(tmp_K_phi_MPI(stind), Ny, MPI_C_DOUBLE_COMPLEX, otherproc, 70, MPI_COMM_WORLD, mpierror)
                call MPI_SEND(tmp_K_T_MPI(stind), Ny, MPI_C_DOUBLE_COMPLEX, otherproc, 71, MPI_COMM_WORLD, mpierror)
                call MPI_SEND(tmp_phi_MPI(stind), Ny, MPI_C_DOUBLE_COMPLEX, otherproc, 72, MPI_COMM_WORLD, mpierror)
                call MPI_SEND(tmp_T_MPI(stind), Ny, MPI_C_DOUBLE_COMPLEX, otherproc, 73, MPI_COMM_WORLD, mpierror)
                call MPI_SEND(tmp_uy_MPI(stind), Ny, MPI_C_DOUBLE_COMPLEX, otherproc, 74, MPI_COMM_WORLD, mpierror)
                call MPI_SEND(uxi_MPI(stind), Ny, MPI_C_DOUBLE_COMPLEX, otherproc, 75, MPI_COMM_WORLD, mpierror)
            end do

            ! Set local to main node variables.
            K1_phi(1:Ny,it) = tmp_K_phi_MPI(1:Ny)
            K1_T(1:Ny,it) = tmp_K_T_MPI(1:Ny)
            phii(1:Ny,it) = tmp_phi_MPI(1:Ny)
            Ti(1:Ny,it) = tmp_T_MPI(1:Ny)
            uyi(1:Ny,it) = tmp_uy_MPI(1:Ny)
            uxi(1:Ny,it) = uxi_MPI(1:Ny)

            
        else if (proc_id == num_procs - 1) then
            ! Send data to main node
            call MPI_SEND(phi(1:Ny-1,it), Ny-1, MPI_C_DOUBLE_COMPLEX, 0, 42, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(K1hat_phi(1:Ny-1,it), Ny-1, MPI_C_DOUBLE_COMPLEX, 0, 43, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(K1hat_T(1:Ny-1,it), Ny-1, MPI_C_DOUBLE_COMPLEX, 0, 46, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(T(1:Ny,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 53, MPI_COMM_WORLD, mpierror)

            ! Send dyv1 and dyv2.
            call MPI_SEND(dyv1_T(it), 1, MPI_DOUBLE, 0, 57, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(dyv2_T(it), 1, MPI_DOUBLE, 0, 58, MPI_COMM_WORLD, mpierror)

            ! Send V1, V2, phi1, phi2
            call MPI_SEND(V1(1:Ny,it), Ny, MPI_DOUBLE, 0, 59, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(V2(1:Ny,it), Ny, MPI_DOUBLE, 0, 60, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(phi1(1:Ny,it), Ny, MPI_DOUBLE, 0, 61, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(phi2(1:Ny,it), Ny, MPI_DOUBLE, 0, 62, MPI_COMM_WORLD, mpierror)

            
            ! Receive K1_phi, K1_T, phii, Ti, uyi, uxi from main node.
            call MPI_RECV(K1_phi(1:Ny,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 70,&
                          MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
            call MPI_RECV(K1_T(1:Ny,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 71,&
                          MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
            call MPI_RECV(phii(1:Ny,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 72,&
                          MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
            call MPI_RECV(Ti(1:Ny,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 73,&
                          MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
            call MPI_RECV(uyi(1:Ny,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 74,&
                          MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
            call MPI_RECV(uxi(1:Ny,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 75,&
                          MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
            
        else 
            ! Send data to main node
            call MPI_SEND(phi(1:Ny,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 42, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(K1hat_phi(1:Ny,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 43, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(K1hat_T(1:Ny,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 46, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(T(1:Ny,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 53, MPI_COMM_WORLD, mpierror)

            ! Send V1, V2, phi1, phi2
            call MPI_SEND(V1(1:Ny,it), Ny, MPI_DOUBLE, 0, 59, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(V2(1:Ny,it), Ny, MPI_DOUBLE, 0, 60, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(phi1(1:Ny,it), Ny, MPI_DOUBLE, 0, 61, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(phi2(1:Ny,it), Ny, MPI_DOUBLE, 0, 62, MPI_COMM_WORLD, mpierror)

            ! Receive K1_phi, K1_T, phii, Ti, uyi, uxi from main node.
            call MPI_RECV(K1_phi(1:Ny,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 70,&
                          MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
            call MPI_RECV(K1_T(1:Ny,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 71,&
                          MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
            call MPI_RECV(phii(1:Ny,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 72,&
                          MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
            call MPI_RECV(Ti(1:Ny,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 73,&
                          MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
            call MPI_RECV(uyi(1:Ny,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 74,&
                          MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
            call MPI_RECV(uxi(1:Ny,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 75,&
                          MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
        end if
    end do
    !$OMP END PARALLEL DO
    finish = OMP_GET_WTIME()
    ! write(*,*) " - stage 1 mid timing: ", finish-start, "(s)"

    ! Compute K2hat
    start = OMP_GET_WTIME()
    call calc_explicit_MPI(2, proc_id, num_procs, proc_id_str)
    finish = OMP_GET_WTIME()
    ! write(*,*) " - calc_explicit(2) timing: ", finish-start, "(s)"
    call MPI_BARRIER(MPI_COMM_WORLD, mpierror)

    !:::::::::::
    ! STAGE 2 ::
    !:::::::::::
    start = OMP_GET_WTIME()
    !$OMP PARALLEL DO private(tmp_phi, tmp_T, tmp_uy, tmp_phi1, tmp_uy1, tmp_K_phi, tmp_K_T) schedule(dynamic)
    do it = 1,Nx ! kx loop
        ! Compute phi1 and T1
        if (proc_id == 0) then
            ! Need to pull all of the variables into local main memory.
            phi_MPI(1:Ny-1) = phi(2:Ny,it)
            K1hat_phi_MPI(1:Ny-1) = K1hat_phi(2:Ny,it)
            K2hat_phi_MPI(1:Ny-1) = K2hat_phi(2:Ny,it)
            K1hat_T_MPI(1:Ny-1) = K1hat_T(2:Ny,it)
            K2hat_T_MPI(1:Ny-1) = K2hat_T(2:Ny,it)
            K1_phi_MPI(1:Ny-1) = K1_phi(2:Ny,it)
            K1_T_MPI(1:Ny-1) = K1_T(2:Ny,it)
            T_MPI(1:Ny) = T(1:Ny,it)

            ! Receive pieces from other nodes.
            do otherproc = 1,num_procs-2
                stind = otherproc * Ny
                call MPI_RECV(phi_MPI(stind), Ny, MPI_C_DOUBLE_COMPLEX, &
                              otherproc, 42, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
                call MPI_RECV(K1hat_phi_MPI(stind), Ny, MPI_C_DOUBLE_COMPLEX, &
                              otherproc, 43, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
                call MPI_RECV(K2hat_phi_MPI(stind), Ny, MPI_C_DOUBLE_COMPLEX, &
                              otherproc, 44, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
                call MPI_RECV(K1hat_T_MPI(stind), Ny, MPI_C_DOUBLE_COMPLEX, & 
                              otherproc, 46, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
                call MPI_RECV(K2hat_T_MPI(stind), Ny, MPI_C_DOUBLE_COMPLEX, & 
                              otherproc, 47, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
                call MPI_RECV(K1_phi_MPI(stind), Ny, MPI_C_DOUBLE_COMPLEX, & 
                              otherproc, 49, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
                call MPI_RECV(K1_T_MPI(stind), Ny, MPI_C_DOUBLE_COMPLEX, & 
                              otherproc, 51, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
                call MPI_RECV(T_MPI(stind+1), Ny, MPI_C_DOUBLE_COMPLEX, & 
                              otherproc, 53, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
            end do
            
            ! Receive from last node.
            stind = (num_procs-1) * Ny
            call MPI_RECV(phi_MPI(stind), Ny-1, MPI_C_DOUBLE_COMPLEX, & 
                          num_procs-1, 42, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
            call MPI_RECV(K1hat_phi_MPI(stind), Ny-1, MPI_C_DOUBLE_COMPLEX, & 
                          num_procs-1, 43, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
            call MPI_RECV(K2hat_phi_MPI(stind), Ny-1, MPI_C_DOUBLE_COMPLEX, & 
                          num_procs-1, 44, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
            call MPI_RECV(K1hat_T_MPI(stind), Ny-1, MPI_C_DOUBLE_COMPLEX, & 
                          num_procs-1, 46, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
            call MPI_RECV(K2hat_T_MPI(stind), Ny-1, MPI_C_DOUBLE_COMPLEX, & 
                          num_procs-1, 47, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
            call MPI_RECV(K1_phi_MPI(stind), Ny-1, MPI_C_DOUBLE_COMPLEX, & 
                          num_procs-1, 49, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
            call MPI_RECV(K1_T_MPI(stind), Ny-1, MPI_C_DOUBLE_COMPLEX, & 
                          num_procs-1, 51, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
            call MPI_RECV(T_MPI(stind+1), Ny, MPI_C_DOUBLE_COMPLEX, & 
                          num_procs-1, 53, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)

        
        

            ! Call vari mod
            call calc_vari_mod_MPI(tmp_phi_MPI, tmp_T_MPI, acoeffs(2,2), 2, total_ny,&
                                kx(it), phi_MPI,&
                                K1hat_phi_MPI,K2hat_phi_MPI,K3hat_phi_MPI,&
                                K1hat_T_MPI,K2hat_T_MPI,K3hat_T_MPI,&
                                K1_phi_MPI, K2_phi_MPI, K1_T_MPI, K2_T_MPI,&
                                g1_total, g2_total, g3_total,&
                                T_MPI)
            
            ! Compute v1 from phi1
            call calc_vi_mod_MPI(tmp_uy_MPI, tmp_phi_MPI, kx(it), total_ny, g1_total, g2_total, g3_total)
            
            ! Receive dyv1 and dyv2
            call MPI_RECV(dyv1_T_it_MPI, 1, MPI_DOUBLE, num_procs-1, 57, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
            call MPI_RECV(dyv2_T_it_MPI, 1, MPI_DOUBLE, num_procs-1, 58, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)

            ! Receive V1, V2, phi1, phi2
            do otherproc = 1,num_procs-1
                stind = otherproc * Ny + 1
                call MPI_RECV(V1_MPI(stind), Ny, MPI_DOUBLE, otherproc, 59, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
                call MPI_RECV(V2_MPI(stind), Ny, MPI_DOUBLE, otherproc, 60, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
                call MPI_RECV(phi1_MPI(stind), Ny, MPI_DOUBLE, otherproc, 61, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
                call MPI_RECV(phi2_MPI(stind), Ny, MPI_DOUBLE, otherproc, 62, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
            end do
            ! Set local variables to MPI.
            V1_MPI(1:Ny) = V1(1:Ny, it)
            V2_MPI(1:Ny) = V2(1:Ny, it)
            phi1_MPI(1:Ny) = phi1(1:Ny, it)
            phi2_MPI(1:Ny) = phi2(1:Ny, it)

            ! BOUNDAY CONDITIONS!
            call update_bcs_mod_MPI(tmp_phi1_MPI,tmp_uy1_MPI, tmp_phi_MPI,tmp_uy_MPI,&
                                dyv1_T_it_MPI,dyv2_T_it_MPI,&
                                dyv1_B(it),dyv2_B(it),&
                                V1_MPI,V2_MPI,phi1_MPI,phi2_MPI,&
                                total_ny, h1_total(total_ny), h2_total(total_ny), h3_total(total_ny))
            
            tmp_phi_MPI = tmp_phi1_MPI
            tmp_uy_MPI  = tmp_uy1_MPI

            
            
            ! Implicit.
            call calc_implicit_mod_MPI(tmp_K_phi_MPI,tmp_K_T_MPI, tmp_phi_MPI,tmp_T_MPI,&
                                       kx(it), total_ny, dynu_MPI, g1_total, g2_total, g3_total)
            
            ! Compute u1 from v1
            if (kx(it) /= 0.0_dp) then
                uxi_MPI = CI*d1y_MPI(tmp_uy_MPI, h1_total, h2_total, h3_total)/kx(it)
            else if (kx(it) == 0.0_dp) then
                uxi_MPI = cmplx(0.0_dp, 0.0_dp, kind=C_DOUBLE_COMPLEX) ! Zero mean flow!
            end if

            ! Send K2_phi, K2_T, phii, Ti, uyi, uxi back to each node.
            do otherproc = 1,num_procs-1
                stind = otherproc * Ny + 1
                call MPI_SEND(tmp_K_phi_MPI(stind), Ny, MPI_C_DOUBLE_COMPLEX, otherproc, 70, MPI_COMM_WORLD, mpierror)
                call MPI_SEND(tmp_K_T_MPI(stind), Ny, MPI_C_DOUBLE_COMPLEX, otherproc, 71, MPI_COMM_WORLD, mpierror)
                call MPI_SEND(tmp_phi_MPI(stind), Ny, MPI_C_DOUBLE_COMPLEX, otherproc, 72, MPI_COMM_WORLD, mpierror)
                call MPI_SEND(tmp_T_MPI(stind), Ny, MPI_C_DOUBLE_COMPLEX, otherproc, 73, MPI_COMM_WORLD, mpierror)
                call MPI_SEND(tmp_uy_MPI(stind), Ny, MPI_C_DOUBLE_COMPLEX, otherproc, 74, MPI_COMM_WORLD, mpierror)
                call MPI_SEND(uxi_MPI(stind), Ny, MPI_C_DOUBLE_COMPLEX, otherproc, 75, MPI_COMM_WORLD, mpierror)
            end do

            ! Set local to main node variables.
            K2_phi(1:Ny,it) = tmp_K_phi_MPI(1:Ny)
            K2_T(1:Ny,it) = tmp_K_T_MPI(1:Ny)
            phii(1:Ny,it) = tmp_phi_MPI(1:Ny)
            Ti(1:Ny,it) = tmp_T_MPI(1:Ny)
            uyi(1:Ny,it) = tmp_uy_MPI(1:Ny)
            uxi(1:Ny,it) = uxi_MPI(1:Ny)

            
        else if (proc_id == num_procs - 1) then
            ! Send data to main node
            call MPI_SEND(phi(1:Ny-1,it), Ny-1, MPI_C_DOUBLE_COMPLEX, 0, 42, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(K1hat_phi(1:Ny-1,it), Ny-1, MPI_C_DOUBLE_COMPLEX, 0, 43, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(K2hat_phi(1:Ny-1,it), Ny-1, MPI_C_DOUBLE_COMPLEX, 0, 44, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(K1hat_T(1:Ny-1,it), Ny-1, MPI_C_DOUBLE_COMPLEX, 0, 46, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(K2hat_T(1:Ny-1,it), Ny-1, MPI_C_DOUBLE_COMPLEX, 0, 47, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(K1_phi(1:Ny-1,it), Ny-1, MPI_C_DOUBLE_COMPLEX, 0, 49, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(K1_T(1:Ny-1,it), Ny-1, MPI_C_DOUBLE_COMPLEX, 0, 51, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(T(1:Ny,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 53, MPI_COMM_WORLD, mpierror)

            ! Send dyv1 and dyv2.
            call MPI_SEND(dyv1_T(it), 1, MPI_DOUBLE, 0, 57, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(dyv2_T(it), 1, MPI_DOUBLE, 0, 58, MPI_COMM_WORLD, mpierror)

            ! Send V1, V2, phi1, phi2
            call MPI_SEND(V1(1:Ny,it), Ny, MPI_DOUBLE, 0, 59, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(V2(1:Ny,it), Ny, MPI_DOUBLE, 0, 60, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(phi1(1:Ny,it), Ny, MPI_DOUBLE, 0, 61, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(phi2(1:Ny,it), Ny, MPI_DOUBLE, 0, 62, MPI_COMM_WORLD, mpierror)
            
           

            ! Receive K2_phi, K2_T, phii, Ti, uyi, uxi from main node.
            call MPI_RECV(K2_phi(1:Ny,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 70,&
                          MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
            call MPI_RECV(K2_T(1:Ny,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 71,&
                          MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
            call MPI_RECV(phii(1:Ny,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 72,&
                          MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
            call MPI_RECV(Ti(1:Ny,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 73,&
                          MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
            call MPI_RECV(uyi(1:Ny,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 74,&
                          MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
            call MPI_RECV(uxi(1:Ny,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 75,&
                          MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
            
        else 
            ! Send data to main node
            call MPI_SEND(phi(1:Ny,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 42, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(K1hat_phi(1:Ny,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 43, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(K2hat_phi(1:Ny,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 44, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(K1hat_T(1:Ny,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 46, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(K2hat_T(1:Ny,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 47, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(K1_phi(1:Ny,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 49, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(K1_T(1:Ny,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 51, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(T(1:Ny,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 53, MPI_COMM_WORLD, mpierror)

            ! Send V1, V2, phi1, phi2
            call MPI_SEND(V1(1:Ny,it), Ny, MPI_DOUBLE, 0, 59, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(V2(1:Ny,it), Ny, MPI_DOUBLE, 0, 60, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(phi1(1:Ny,it), Ny, MPI_DOUBLE, 0, 61, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(phi2(1:Ny,it), Ny, MPI_DOUBLE, 0, 62, MPI_COMM_WORLD, mpierror)
    

            ! Receive K2_phi, K2_T, phii, Ti, uyi, uxi from main node.
            call MPI_RECV(K2_phi(1:Ny,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 70,&
                          MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
            call MPI_RECV(K2_T(1:Ny,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 71,&
                          MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
            call MPI_RECV(phii(1:Ny,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 72,&
                          MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
            call MPI_RECV(Ti(1:Ny,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 73,&
                          MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
            call MPI_RECV(uyi(1:Ny,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 74,&
                          MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
            call MPI_RECV(uxi(1:Ny,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 75,&
                          MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
        end if
    end do
    !$OMP END PARALLEL DO
    finish = OMP_GET_WTIME()
    ! write(*,*) " - stage 2 mid timing: ", finish-start, "(s)"
    call MPI_BARRIER(MPI_COMM_WORLD, mpierror)
    ! Compute K3hat
    start = OMP_GET_WTIME()
    call calc_explicit_MPI(3, proc_id, num_procs, proc_id_str)
    finish = OMP_GET_WTIME()
    ! write(*,*) " - calc_explicit(3) timing: ", finish-start, "(s)"
    call MPI_BARRIER(MPI_COMM_WORLD, mpierror)

    !:::::::::::
    ! STAGE 3 ::
    !:::::::::::
    start = OMP_GET_WTIME()
    !$OMP PARALLEL DO private(tmp_phi, tmp_T, tmp_uy, tmp_phi1, tmp_uy1, tmp_K_phi, tmp_K_T) schedule(dynamic)
    do it = 1,Nx ! kx loop
        ! Compute phi1 and T1
        if (proc_id == 0) then
            ! Need to pull all of the variables into local main memory.
            phi_MPI(1:Ny-1) = phi(2:Ny,it)
            K1hat_phi_MPI(1:Ny-1) = K1hat_phi(2:Ny,it)
            K2hat_phi_MPI(1:Ny-1) = K2hat_phi(2:Ny,it)
            K3hat_phi_MPI(1:Ny-1) = K3hat_phi(2:Ny,it)
            K1hat_T_MPI(1:Ny-1) = K1hat_T(2:Ny,it)
            K2hat_T_MPI(1:Ny-1) = K2hat_T(2:Ny,it)
            K3hat_t_MPI(1:Ny-1) = K3hat_T(2:Ny,it)
            K1_phi_MPI(1:Ny-1) = K1_phi(2:Ny,it)
            K2_phi_MPI(1:Ny-1) = K2_phi(2:Ny,it)
            K1_T_MPI(1:Ny-1) = K1_T(2:Ny,it)
            K2_T_MPI(1:Ny-1) = K2_T(2:Ny,it)
            T_MPI(1:Ny) = T(1:Ny,it)

            ! Receive pieces from other nodes.
            do otherproc = 1,num_procs-2
                stind = otherproc * Ny
                call MPI_RECV(phi_MPI(stind), Ny, MPI_C_DOUBLE_COMPLEX, &
                              otherproc, 42, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
                call MPI_RECV(K1hat_phi_MPI(stind), Ny, MPI_C_DOUBLE_COMPLEX, &
                              otherproc, 43, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
                call MPI_RECV(K2hat_phi_MPI(stind), Ny, MPI_C_DOUBLE_COMPLEX, &
                              otherproc, 44, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
                call MPI_RECV(K3hat_phi_MPI(stind), Ny, MPI_C_DOUBLE_COMPLEX, & 
                              otherproc, 45, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
                call MPI_RECV(K1hat_T_MPI(stind), Ny, MPI_C_DOUBLE_COMPLEX, & 
                              otherproc, 46, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
                call MPI_RECV(K2hat_T_MPI(stind), Ny, MPI_C_DOUBLE_COMPLEX, & 
                              otherproc, 47, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
                call MPI_RECV(K3hat_T_MPI(stind), Ny, MPI_C_DOUBLE_COMPLEX, & 
                              otherproc, 48, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
                call MPI_RECV(K1_phi_MPI(stind), Ny, MPI_C_DOUBLE_COMPLEX, & 
                              otherproc, 49, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
                call MPI_RECV(K2_phi_MPI(stind), Ny, MPI_C_DOUBLE_COMPLEX, & 
                              otherproc, 50, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
                call MPI_RECV(K1_T_MPI(stind), Ny, MPI_C_DOUBLE_COMPLEX, & 
                              otherproc, 51, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
                call MPI_RECV(K2_T_MPI(stind), Ny, MPI_C_DOUBLE_COMPLEX, & 
                              otherproc, 52, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
                call MPI_RECV(T_MPI(stind+1), Ny, MPI_C_DOUBLE_COMPLEX, & 
                              otherproc, 53, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
            end do
            
            ! Receive from last node.
            stind = (num_procs-1) * Ny
            call MPI_RECV(phi_MPI(stind), Ny-1, MPI_C_DOUBLE_COMPLEX, & 
                          num_procs-1, 42, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
            call MPI_RECV(K1hat_phi_MPI(stind), Ny-1, MPI_C_DOUBLE_COMPLEX, & 
                          num_procs-1, 43, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
            call MPI_RECV(K2hat_phi_MPI(stind), Ny-1, MPI_C_DOUBLE_COMPLEX, & 
                          num_procs-1, 44, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
            call MPI_RECV(K3hat_phi_MPI(stind), Ny-1, MPI_C_DOUBLE_COMPLEX, & 
                          num_procs-1, 45, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
            call MPI_RECV(K1hat_T_MPI(stind), Ny-1, MPI_C_DOUBLE_COMPLEX, & 
                          num_procs-1, 46, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
            call MPI_RECV(K2hat_T_MPI(stind), Ny-1, MPI_C_DOUBLE_COMPLEX, & 
                          num_procs-1, 47, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
            call MPI_RECV(K3hat_T_MPI(stind), Ny-1, MPI_C_DOUBLE_COMPLEX, & 
                          num_procs-1, 48, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
            call MPI_RECV(K1_phi_MPI(stind), Ny-1, MPI_C_DOUBLE_COMPLEX, & 
                          num_procs-1, 49, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
            call MPI_RECV(K2_phi_MPI(stind), Ny-1, MPI_C_DOUBLE_COMPLEX, & 
                          num_procs-1, 50, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
            call MPI_RECV(K1_T_MPI(stind), Ny-1, MPI_C_DOUBLE_COMPLEX, & 
                          num_procs-1, 51, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
            call MPI_RECV(K2_T_MPI(stind), Ny-1, MPI_C_DOUBLE_COMPLEX, & 
                          num_procs-1, 52, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
            call MPI_RECV(T_MPI(stind+1), Ny, MPI_C_DOUBLE_COMPLEX, & 
                          num_procs-1, 53, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)

            ! Call vari mod
            call calc_vari_mod_MPI(tmp_phi_MPI, tmp_T_MPI, acoeffs(3,3), 3, total_ny,&
                                kx(it), phi_MPI,&
                                K1hat_phi_MPI,K2hat_phi_MPI,K3hat_phi_MPI,&
                                K1hat_T_MPI,K2hat_T_MPI,K3hat_T_MPI,&
                                K1_phi_MPI, K2_phi_MPI, K1_T_MPI, K2_T_MPI,&
                                g1_total, g2_total, g3_total,&
                                T_MPI)
            
            ! Compute v1 from phi1
            call calc_vi_mod_MPI(tmp_uy_MPI, tmp_phi_MPI, kx(it), total_ny, g1_total, g2_total, g3_total)
            
            ! Receive dyv1 and dyv2
            call MPI_RECV(dyv1_T_it_MPI, 1, MPI_DOUBLE, num_procs-1, 57, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
            call MPI_RECV(dyv2_T_it_MPI, 1, MPI_DOUBLE, num_procs-1, 58, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)

            ! Receive V1, V2, phi1, phi2
            do otherproc = 1,num_procs-1
                stind = otherproc * Ny + 1
                call MPI_RECV(V1_MPI(stind), Ny, MPI_DOUBLE, otherproc, 59, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
                call MPI_RECV(V2_MPI(stind), Ny, MPI_DOUBLE, otherproc, 60, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
                call MPI_RECV(phi1_MPI(stind), Ny, MPI_DOUBLE, otherproc, 61, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
                call MPI_RECV(phi2_MPI(stind), Ny, MPI_DOUBLE, otherproc, 62, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
            end do
            ! Set local variables to MPI.
            V1_MPI(1:Ny) = V1(1:Ny, it)
            V2_MPI(1:Ny) = V2(1:Ny, it)
            phi1_MPI(1:Ny) = phi1(1:Ny, it)
            phi2_MPI(1:Ny) = phi2(1:Ny, it)

            ! BOUNDAY CONDITIONS!
            call update_bcs_mod_MPI(tmp_phi1_MPI,tmp_uy1_MPI, tmp_phi_MPI,tmp_uy_MPI,&
                                dyv1_T_it_MPI,dyv2_T_it_MPI,&
                                dyv1_B(it),dyv2_B(it),&
                                V1_MPI,V2_MPI,phi1_MPI,phi2_MPI,&
                                total_ny, h1_total(total_ny), h2_total(total_ny), h3_total(total_ny))
            
            tmp_phi_MPI = tmp_phi1_MPI
            tmp_uy_MPI  = tmp_uy1_MPI

            
            ! Implicit.
            call calc_implicit_mod_MPI(tmp_K_phi_MPI,tmp_K_T_MPI, tmp_phi_MPI,tmp_T_MPI,&
                                       kx(it), total_ny, dynu_MPI, g1_total, g2_total, g3_total)

            ! Compute u1 from v1
            if (kx(it) /= 0.0_dp) then
                uxi_MPI = CI*d1y_MPI(tmp_uy_MPI, h1_total, h2_total, h3_total)/kx(it)
            else if (kx(it) == 0.0_dp) then
                uxi_MPI = cmplx(0.0_dp, 0.0_dp, kind=C_DOUBLE_COMPLEX) ! Zero mean flow!
            end if

            ! Send K3_phi, K3_T, phii, Ti, uyi, uxi back to each node.
            do otherproc = 1,num_procs-1
                stind = otherproc * Ny + 1
                call MPI_SEND(tmp_K_phi_MPI(stind), Ny, MPI_C_DOUBLE_COMPLEX, otherproc, 70, MPI_COMM_WORLD, mpierror)
                call MPI_SEND(tmp_K_T_MPI(stind), Ny, MPI_C_DOUBLE_COMPLEX, otherproc, 71, MPI_COMM_WORLD, mpierror)
                call MPI_SEND(tmp_phi_MPI(stind), Ny, MPI_C_DOUBLE_COMPLEX, otherproc, 72, MPI_COMM_WORLD, mpierror)
                call MPI_SEND(tmp_T_MPI(stind), Ny, MPI_C_DOUBLE_COMPLEX, otherproc, 73, MPI_COMM_WORLD, mpierror)
                call MPI_SEND(tmp_uy_MPI(stind), Ny, MPI_C_DOUBLE_COMPLEX, otherproc, 74, MPI_COMM_WORLD, mpierror)
                call MPI_SEND(uxi_MPI(stind), Ny, MPI_C_DOUBLE_COMPLEX, otherproc, 75, MPI_COMM_WORLD, mpierror)
            end do

            ! Set local to main node variables.
            K3_phi(1:Ny,it) = tmp_K_phi_MPI(1:Ny)
            K3_T(1:Ny,it) = tmp_K_T_MPI(1:Ny)
            phii(1:Ny,it) = tmp_phi_MPI(1:Ny)
            Ti(1:Ny,it) = tmp_T_MPI(1:Ny)
            uyi(1:Ny,it) = tmp_uy_MPI(1:Ny)
            uxi(1:Ny,it) = uxi_MPI(1:Ny)

            
        else if (proc_id == num_procs - 1) then
            ! Send data to main node
            call MPI_SEND(phi(1:Ny-1,it), Ny-1, MPI_C_DOUBLE_COMPLEX, 0, 42, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(K1hat_phi(1:Ny-1,it), Ny-1, MPI_C_DOUBLE_COMPLEX, 0, 43, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(K2hat_phi(1:Ny-1,it), Ny-1, MPI_C_DOUBLE_COMPLEX, 0, 44, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(K3hat_phi(1:Ny-1,it), Ny-1, MPI_C_DOUBLE_COMPLEX, 0, 45, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(K1hat_T(1:Ny-1,it), Ny-1, MPI_C_DOUBLE_COMPLEX, 0, 46, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(K2hat_T(1:Ny-1,it), Ny-1, MPI_C_DOUBLE_COMPLEX, 0, 47, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(K3hat_T(1:Ny-1,it), Ny-1, MPI_C_DOUBLE_COMPLEX, 0, 48, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(K1_phi(1:Ny-1,it), Ny-1, MPI_C_DOUBLE_COMPLEX, 0, 49, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(K2_phi(1:Ny-1,it), Ny-1, MPI_C_DOUBLE_COMPLEX, 0, 50, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(K1_T(1:Ny-1,it), Ny-1, MPI_C_DOUBLE_COMPLEX, 0, 51, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(K2_T(1:Ny-1,it), Ny-1, MPI_C_DOUBLE_COMPLEX, 0, 52, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(T(1:Ny,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 53, MPI_COMM_WORLD, mpierror)

            ! Send dyv1 and dyv2.
            call MPI_SEND(dyv1_T(it), 1, MPI_DOUBLE, 0, 57, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(dyv2_T(it), 1, MPI_DOUBLE, 0, 58, MPI_COMM_WORLD, mpierror)

            ! Send V1, V2, phi1, phi2
            call MPI_SEND(V1(1:Ny,it), Ny, MPI_DOUBLE, 0, 59, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(V2(1:Ny,it), Ny, MPI_DOUBLE, 0, 60, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(phi1(1:Ny,it), Ny, MPI_DOUBLE, 0, 61, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(phi2(1:Ny,it), Ny, MPI_DOUBLE, 0, 62, MPI_COMM_WORLD, mpierror)
            


            ! Receive K3_phi, K3_T, phii, Ti, uyi, uxi from main node.
            call MPI_RECV(K3_phi(1:Ny,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 70,&
                          MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
            call MPI_RECV(K3_T(1:Ny,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 71,&
                          MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
            call MPI_RECV(phii(1:Ny,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 72,&
                          MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
            call MPI_RECV(Ti(1:Ny,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 73,&
                          MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
            call MPI_RECV(uyi(1:Ny,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 74,&
                          MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
            call MPI_RECV(uxi(1:Ny,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 75,&
                          MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
            
        else 
            ! Send data to main node
            call MPI_SEND(phi(1:Ny,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 42, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(K1hat_phi(1:Ny,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 43, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(K2hat_phi(1:Ny,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 44, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(K3hat_phi(1:Ny,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 45, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(K1hat_T(1:Ny,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 46, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(K2hat_T(1:Ny,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 47, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(K3hat_T(1:Ny,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 48, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(K1_phi(1:Ny,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 49, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(K2_phi(1:Ny,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 50, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(K1_T(1:Ny,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 51, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(K2_T(1:Ny,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 52, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(T(1:Ny,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 53, MPI_COMM_WORLD, mpierror)

            ! Send V1, V2, phi1, phi2
            call MPI_SEND(V1(1:Ny,it), Ny, MPI_DOUBLE, 0, 59, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(V2(1:Ny,it), Ny, MPI_DOUBLE, 0, 60, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(phi1(1:Ny,it), Ny, MPI_DOUBLE, 0, 61, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(phi2(1:Ny,it), Ny, MPI_DOUBLE, 0, 62, MPI_COMM_WORLD, mpierror)

            ! Receive K3_phi, K3_T, phii, Ti, uyi, uxi from main node.
            call MPI_RECV(K3_phi(1:Ny,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 70,&
                          MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
            call MPI_RECV(K3_T(1:Ny,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 71,&
                          MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
            call MPI_RECV(phii(1:Ny,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 72,&
                          MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
            call MPI_RECV(Ti(1:Ny,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 73,&
                          MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
            call MPI_RECV(uyi(1:Ny,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 74,&
                          MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
            call MPI_RECV(uxi(1:Ny,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 75,&
                          MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
        end if
    end do
    !$OMP END PARALLEL DO
    finish = OMP_GET_WTIME()
    ! write(*,*) " - stage 3 mid timing: ", finish-start, "(s)"
    call MPI_BARRIER(MPI_COMM_WORLD, mpierror)
    ! Compute K4hat
    start = OMP_GET_WTIME()
    call calc_explicit_MPI(4, proc_id, num_procs, proc_id_str)
    finish = OMP_GET_WTIME()
    ! write(*,*) " - calc_explicit(4) timing: ", finish-start, "(s)"
    call MPI_BARRIER(MPI_COMM_WORLD, mpierror)

    start = OMP_GET_WTIME()
    ! UPDATE SOLUTIONS
    if (proc_id == 0) then
        ! Get phi
        phi(2:Ny,:) = phi(2:Ny,:) + dt*(b(1)*(K1_phi(2:Ny,:) + K2hat_phi(2:Ny,:)) + &
        &                    b(2)*(K2_phi(2:Ny,:) + K3hat_phi(2:Ny,:)) + &
        &                    b(3)*(K3_phi(2:Ny,:) + K4hat_phi(2:Ny,:)))

        ! Get temperature
        T(2:Ny,:)   = T(2:Ny,:)  + dt*(b(1)*(K1_T(2:Ny,:) + K2hat_T(2:Ny,:)) + &
        &                   b(2)*(K2_T(2:Ny,:) + K3hat_T(2:Ny,:)) + &
        &                   b(3)*(K3_T(2:Ny,:) + K4hat_T(2:Ny,:)))
    else if (proc_id == num_procs - 1) then
        ! Get phi
        phi(1:Ny-1,:) = phi(1:Ny-1,:) + dt*(b(1)*(K1_phi(1:Ny-1,:) + K2hat_phi(1:Ny-1,:)) + &
        &                    b(2)*(K2_phi(1:Ny-1,:) + K3hat_phi(1:Ny-1,:)) + &
        &                    b(3)*(K3_phi(1:Ny-1,:) + K4hat_phi(1:Ny-1,:)))

        ! Get temperature
        T(1:Ny-1,:)   = T(1:Ny-1,:)  + dt*(b(1)*(K1_T(1:Ny-1,:) + K2hat_T(1:Ny-1,:)) + &
        &                   b(2)*(K2_T(1:Ny-1,:) + K3hat_T(1:Ny-1,:)) + &
        &                   b(3)*(K3_T(1:Ny-1,:) + K4hat_T(1:Ny-1,:)))    
    else 
        ! Get phi
        phi(1:Ny,:) = phi(1:Ny,:) + dt*(b(1)*(K1_phi(1:Ny,:) + K2hat_phi(1:Ny,:)) + &
        &                    b(2)*(K2_phi(1:Ny,:) + K3hat_phi(1:Ny,:)) + &
        &                    b(3)*(K3_phi(1:Ny,:) + K4hat_phi(1:Ny,:)))
    
        ! Get temperature
        T(1:Ny,:)   = T(1:Ny,:)  + dt*(b(1)*(K1_T(1:Ny,:) + K2hat_T(1:Ny,:)) + &
        &                   b(2)*(K2_T(1:Ny,:) + K3hat_T(1:Ny,:)) + &
        &                   b(3)*(K3_T(1:Ny,:) + K4hat_T(1:Ny,:)))
    end if

    !$OMP PARALLEL DO private(tmp_uy, it) schedule(dynamic)
    do it = 1,Nx
        ! Only call from single node.
        if (proc_id == 0) then 

            do otherproc= 1, num_procs-1
                stind = otherproc * Ny + 1
                ! Receive phi
                call MPI_RECV(phi_update_MPI(stind), Ny, MPI_C_DOUBLE_COMPLEX, &
                              otherproc, 86, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
            end do
            ! Set local phi
            phi_update_MPI(1:Ny) = phi(1:Ny,it)

            ! Solve for v
            call calc_vi_mod_MPI(tmp_uy_MPI, phi_update_MPI, kx(it), total_ny, g1_total, g2_total, g3_total)
            ! Solve for u
            if (kx(it) /= 0.0_dp) then
                uxi_MPI = CI*d1y_MPI(tmp_uy_MPI, h1_total, h2_total, h3_total)/kx(it)
            else if (kx(it) == 0.0_dp) then
                uxi_MPI = cmplx(0.0_dp, 0.0_dp, kind=C_DOUBLE_COMPLEX) ! Zero mean flow!
            end if

            ! Send uy, ux back to each node.
            do otherproc = 1,num_procs-1
                stind = otherproc * Ny + 1
                call MPI_SEND(tmp_uy_MPI(stind), Ny, MPI_C_DOUBLE_COMPLEX, otherproc, 87, MPI_COMM_WORLD, mpierror)
                call MPI_SEND(uxi_MPI(stind), Ny, MPI_C_DOUBLE_COMPLEX, otherproc, 88, MPI_COMM_WORLD, mpierror)
            end do
            ! Set local versions.
            uy(1:Ny,it) = tmp_uy_MPI(1:Ny)
            ux(1:Ny,it) = uxi_MPI(1:Ny)

        else 
            ! Send phi
            call MPI_SEND(phi(1:Ny,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 86, MPI_COMM_WORLD, mpierror)

            ! Receive uy, ux from main node.
            call MPI_RECV(uy(1:Ny,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 87,&
                          MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
            call MPI_RECV(ux(1:Ny,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 88,&
                          MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
        end if 
    end do
    !$OMP END PARALLEL DO
    finish = OMP_GET_WTIME()
    ! write(*,*) " - update sols timing: ", finish-start, "(s)"
    call MPI_BARRIER(MPI_COMM_WORLD, mpierror)    

    if (time == t_final) then
        ! open(unit=9010, file="P"//proc_id_str//"T_real_update.txt", action="write", status="unknown")
        ! open(unit=9011, file="P"//proc_id_str//"T_im_update.txt", action="write", status="unknown")
        ! do i=1,Ny
        !     do j=1,Nx
        !         write (9010,*) REAL(T(i,j))
        !         write (9011,*) AIMAG(T(i,j))
        !     end do
        ! end do
        ! close(unit=9010)
        ! close(unit=9011)
        ! write(*,*) Ny
        ! write(*,*) "done writing T!"
        exit
    end if

    if (wvtk) then
        call write_to_vtk(nti, .false., proc_id_str) ! false = Fourier space
    end if

    call MPI_BARRIER(MPI_COMM_WORLD, mpierror)
    ! write(*,*) "proc ", proc_id_str, " write vtk complete for t", time 
 
 end do

end subroutine imex_rk_MPI

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine calc_explicit_MPI(stage, proc_id, num_procs, proc_id_str)

integer             :: i, j
integer, intent(in) :: stage, proc_id, num_procs
character(3), intent(in) :: proc_id_str
real(dp)            :: start, finish

start = OMP_GET_WTIME()
select case(stage)
    case (1)
        !$OMP PARALLEL DO schedule(dynamic)
        do i = 1,Nx
            K1hat_phi(:,i) = -kx(i)**2.0_dp*Ti(:,i)
        end do
        !$OMP END PARALLEL DO
    case (2)
        !$OMP PARALLEL DO schedule(dynamic)
        do i = 1,Nx
            K2hat_phi(:,i) = -kx(i)**2.0_dp*Ti(:,i)
        end do
        !$OMP END PARALLEL DO
    case (3)
        !$OMP PARALLEL DO schedule(dynamic)
        do i = 1,Nx
            K3hat_phi(:,i) = -kx(i)**2.0_dp*Ti(:,i)
        end do
        !$OMP END PARALLEL DO
    case (4)
        !$OMP PARALLEL DO schedule(dynamic)
        do i = 1,Nx
            K4hat_phi(:,i) = -kx(i)**2.0_dp*Ti(:,i)
        end do
        !$OMP END PARALLEL DO
end select
finish = OMP_GET_WTIME()
! write(*,*) " - - l1 timing: ", finish-start, "(s)"

start = OMP_GET_WTIME()
!$OMP PARALLEL DO schedule(dynamic)
do i=1,Nx
   ! Compute dx(T) in Fourier space
   nlT  (:,i) =  kx(i)*Ti(:,i)
   ! Compute D2(ux)
   nlphi(:,i) = -kx(i)**2.0_dp*uxi(:,i) + d2y_MPI2(uxi(:,i), proc_id, num_procs)
end do
!$OMP END PARALLEL DO
finish = OMP_GET_WTIME()

! write(*,*) " - - l2 timing: ", finish-start, "(s)"

!nlT = -CI*nlT
nlT = CI*nlT

start = OMP_GET_WTIME()
!$OMP PARALLEL DO private(tnlT, tnlphi, tT, tux, tuy, tphi) schedule(dynamic)
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
!$OMP END PARALLEL DO
finish = OMP_GET_WTIME()
! write(*,*) " - - l3 timing: ", finish-start, "(s)"

! Calculate nonlinear term
start = OMP_GET_WTIME()
!$OMP PARALLEL DO private(tmp_T) schedule(dynamic)
do i = 1,Nx
   ! Temperature
   tmp_T = Ti(:,i)
   nlT(:,i) = uxi(:,i)*nlT(:,i) + uyi(:,i)*d1y_MPI2(tmp_T, proc_id, num_procs)
   ! phi
   nlphi(:,i) = uxi(:,i)*phii(:,i) - uyi(:,i)*nlphi(:,i)
end do
!$OMP END PARALLEL DO
finish = OMP_GET_WTIME()
! write(*,*) " - - l4 timing: ", finish-start, "(s)"

! Bring nonlinear terms back to Fourier space
start = OMP_GET_WTIME()
!$OMP PARALLEL DO private(tnlT, tnlphi) schedule(dynamic)
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
!$OMP END PARALLEL DO
finish = OMP_GET_WTIME()
! write(*,*) " - - l5 timing: ", finish-start, "(s)"

nlT   = nlT   / real(Nx,kind=dp)
nlphi = nlphi / real(Nx,kind=dp)

start = OMP_GET_WTIME()
select case (stage)
   case (1)
      !$OMP PARALLEL DO schedule(dynamic)
      do i = 1,Nx
         !K1hat_phi(:,i) = K1hat_phi(:,i) + CI*kx(i)*nlphi(:,i)
         K1hat_phi(:,i) = K1hat_phi(:,i) - CI*kx(i)*nlphi(:,i)
      end do
      !$OMP END PARALLEL DO
      K1hat_T = -nlT
   case (2)
      !$OMP PARALLEL DO schedule(dynamic)
      do i = 1,Nx
         !K2hat_phi(:,i) = K2hat_phi(:,i) + CI*kx(i)*nlphi(:,i)
         K2hat_phi(:,i) = K2hat_phi(:,i) - CI*kx(i)*nlphi(:,i)
      end do
      !$OMP END PARALLEL DO
      K2hat_T = -nlT
   case (3)
      !$OMP PARALLEL DO schedule(dynamic)
      do i = 1,Nx
         !K3hat_phi(:,i) = K3hat_phi(:,i) + CI*kx(i)*nlphi(:,i)
         K3hat_phi(:,i) = K3hat_phi(:,i) - CI*kx(i)*nlphi(:,i)
      end do
      !$OMP END PARALLEL DO
      K3hat_T = -nlT
   case (4)
      !$OMP PARALLEL DO schedule(dynamic)
      do i = 1,Nx
         !K4hat_phi(:,i) = K4hat_phi(:,i) + CI*kx(i)*nlphi(:,i)
         K4hat_phi(:,i) = K4hat_phi(:,i) - CI*kx(i)*nlphi(:,i)
      end do
      !$OMP END PARALLEL DO
      K4hat_T = -nlT
end select
finish = OMP_GET_WTIME()
! write(*,*) " - - l6 timing: ", finish-start, "(s)"

end subroutine calc_explicit_MPI

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
subroutine calc_vari_mod_MPI(phiout,Tout, aii, stage, total_ny, kx_it, phi_in, &
    k1hat_phi_in, k2hat_phi_in, k3hat_phi_in,&
    k1hat_T_in, k2hat_T_in, k3hat_T_in,&
    k1_phi_in, k2_phi_in, k1_T_in, k2_T_in,&
    g1_total, g2_total, g3_total,&
    T_in)

real(dp),                                               intent(in)  :: aii
integer,                                                intent(in)  :: stage, total_ny
real(dp),                                               intent(in)  :: kx_it
complex(C_DOUBLE_COMPLEX),              dimension(total_ny-2), intent(in)  :: phi_in
complex(C_DOUBLE_COMPLEX),              dimension(total_ny-2), intent(in)  :: k1hat_phi_in
complex(C_DOUBLE_COMPLEX),              dimension(total_ny-2), intent(in)  :: k2hat_phi_in
complex(C_DOUBLE_COMPLEX),              dimension(total_ny-2), intent(in)  :: k3hat_phi_in
complex(C_DOUBLE_COMPLEX),              dimension(total_ny-2), intent(in)  :: k1hat_T_in
complex(C_DOUBLE_COMPLEX),              dimension(total_ny-2), intent(in)  :: k2hat_T_in
complex(C_DOUBLE_COMPLEX),              dimension(total_ny-2), intent(in)  :: k3hat_T_in
complex(C_DOUBLE_COMPLEX),              dimension(total_ny),   intent(in)  :: T_in
complex(C_DOUBLE_COMPLEX),              dimension(total_ny-2), intent(in)  :: k1_phi_in
complex(C_DOUBLE_COMPLEX),              dimension(total_ny-2), intent(in)  :: k2_phi_in
complex(C_DOUBLE_COMPLEX),              dimension(total_ny-2), intent(in)  :: k1_T_in
complex(C_DOUBLE_COMPLEX),              dimension(total_ny-2), intent(in)  :: k2_T_in
real(dp),              dimension(total_ny), intent(in)  :: g1_total, g2_total, g3_total
complex(C_DOUBLE_COMPLEX), allocatable, dimension(:),   intent(out) :: phiout, Tout
complex(C_DOUBLE_COMPLEX), allocatable, dimension(:)              :: Fphi, FT
real(dp),                  allocatable, dimension(:,:)            :: phi_rhs, T_rhs
real(dp),                  allocatable, dimension(:)              :: dphi, duphi, dlphi
real(dp),                  allocatable, dimension(:)              :: ddT, duT, dlT
integer :: i

allocate(dphi(total_ny-2), ddT(total_ny-2), stat=alloc_err)
call check_alloc_err(alloc_err)
allocate(duphi(total_ny-3), dlphi(total_ny-3), duT(total_ny-3), dlT(total_ny-3), stat=alloc_err)
call check_alloc_err(alloc_err)
allocate(phiout(total_ny), Tout(total_ny), stat=alloc_err)
call check_alloc_err(alloc_err)
allocate(Fphi(total_ny-2), FT(total_ny-2), stat=alloc_err)
call check_alloc_err(alloc_err)
allocate(phi_rhs(total_ny-2,2), T_rhs(total_ny-2,2), stat=alloc_err)
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
do jt = 2,total_ny-1
    ddT (jt-1) = 1.0_dp - kappa0*dt*aii*(-kx_it**2.0_dp + g2_total(jt))
    dphi(jt-1) = 1.0_dp - nu0   *dt*aii*(-kx_it**2.0_dp + g2_total(jt))
end do

do jt = 2,total_ny-2
    duT  (jt-1) = -kappa0*dt*g3_total(jt)*aii
    duphi(jt-1) = -nu0   *dt*g3_total(jt)*aii
end do

do jt = 3,total_ny-1
    dlT  (jt-2) = -kappa0*dt*g1_total(jt)*aii
    dlphi(jt-2) = -nu0   *dt*g1_total(jt)*aii
end do

select case (stage)
    case(1)
        Fphi     = phi_in + dt*ahatcoeffs(2,1)*k1hat_phi_in
        FT       = T_in(2:total_ny-1) + dt*ahatcoeffs(2,1)*k1hat_T_in
        FT(1)    = FT(1) + kappa0*dt*aii*g1_total(2)*T_in(1) ! b/c Ti(y_1) = T(y_1)
        FT(total_ny-2) = FT(total_ny-2) + kappa0*dt*aii*g3_total(total_ny-1)*T_in(total_ny) ! b/c Ti(total_ny) = T(total_ny)

        phi_rhs(:,1) = real(Fphi)
        phi_rhs(:,2) = aimag(Fphi)

        T_rhs  (:,1) = real(FT)
        T_rhs  (:,2) = aimag(FT)

        call dgtsv(total_ny-2, 2, dlT, ddT, duT, T_rhs, total_ny-2, info)
        Tout(2:total_ny-1) = cmplx(T_rhs(:,1), T_rhs(:,2), kind=C_DOUBLE_COMPLEX)
        ! Set temperature boundary conditions
        Tout(1) = T_in(1)
        Tout(total_ny) = T_in(total_ny)

        call dgtsv(total_ny-2, 2, dlphi, dphi, duphi, phi_rhs, total_ny-2, info)
        phiout(2:total_ny-1) = cmplx(phi_rhs(:,1), phi_rhs(:,2), kind=C_DOUBLE_COMPLEX)

    case(2)
        Fphi = phi_in + dt*(acoeffs(2,1)*k1_phi_in       + &
        &                     ahatcoeffs(3,1)*k1hat_phi_in + &
        &                     ahatcoeffs(3,2)*k2hat_phi_in)
        FT   = T_in(2:total_ny-1) + dt*(acoeffs(2,1)*k1_T_in      + &
        &                     ahatcoeffs(3,1)*k1hat_T_in + &
        &                     ahatcoeffs(3,2)*k2hat_T_in)
        FT(1)    = FT(1) + kappa0*dt*aii*g1_total(2)*T_in(1) ! b/c Ti(y_1) = T(y_1)
        FT(total_ny-2) = FT(total_ny-2) + kappa0*dt*aii*g3_total(total_ny-1)*T_in(total_ny) ! b/c Ti(total_ny) = T(total_ny)

        ! do i = 1,total_ny-2
        !     write(*,*) i, k2hat_phi_in(i)
        ! end do 

        phi_rhs(:,1) = real(Fphi)
        phi_rhs(:,2) = aimag(Fphi)

        T_rhs  (:,1) = real(FT)
        T_rhs  (:,2) = aimag(FT)

        call dgtsv(total_ny-2, 2, dlT, ddT, duT, T_rhs, total_ny-2, info)
        Tout(2:total_ny-1) = cmplx(T_rhs(:,1), T_rhs(:,2), kind=C_DOUBLE_COMPLEX)
        ! Set temperature boundary conditions
        Tout(1) = T_in(1)
        Tout(total_ny) = T_in(total_ny)

        call dgtsv(total_ny-2, 2, dlphi, dphi, duphi, phi_rhs, total_ny-2, info)
        phiout(2:total_ny-1) = cmplx(phi_rhs(:,1), phi_rhs(:,2), kind=C_DOUBLE_COMPLEX)

    case(3)
        Fphi = phi_in + dt*(acoeffs(3,1)*k1_phi_in       + &
                    &acoeffs(3,2)*k2_phi_in       + &
                    &ahatcoeffs(4,1)*k1hat_phi_in + &
                    &ahatcoeffs(4,2)*k2hat_phi_in + &
                    &ahatcoeffs(4,3)*k3hat_phi_in)
        FT   = T_in(2:total_ny-1) + dt*(acoeffs(3,1)*k1_T_in       + &
                    &acoeffs(3,2)*k2_T_in     + &
                    &ahatcoeffs(4,1)*k1hat_T_in + &
                    &ahatcoeffs(4,2)*k2hat_T_in + &
                    &ahatcoeffs(4,3)*k3hat_T_in)
        FT(1)    = FT(1) + kappa0*dt*aii*g1_total(2)*T_in(1) ! b/c Ti(y_1) = T(y_1)
        FT(total_ny-2) = FT(total_ny-2) + kappa0*dt*aii*g3_total(total_ny-1)*T_in(total_ny) ! b/c Ti(total_ny) = T(total_ny)

        phi_rhs(:,1) = real(Fphi)
        phi_rhs(:,2) = aimag(Fphi)

        T_rhs  (:,1) = real(FT)
        T_rhs  (:,2) = aimag(FT)

        call dgtsv(total_ny-2, 2, dlT, ddT, duT, T_rhs, total_ny-2, info)
        Tout(2:total_ny-1) = cmplx(T_rhs(:,1), T_rhs(:,2), kind=C_DOUBLE_COMPLEX)
        ! Set temperature boundary conditions
        Tout(1) = T_in(1)
        Tout(total_ny) = T_in(total_ny)

        call dgtsv(total_ny-2, 2, dlphi, dphi, duphi, phi_rhs, total_ny-2, info)
        phiout(2:total_ny-1) = cmplx(phi_rhs(:,1), phi_rhs(:,2), kind=C_DOUBLE_COMPLEX)

end select

end subroutine calc_vari_mod_MPI
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!

subroutine calc_vi_mod_MPI(vi, phiin, kx_it, total_ny, g1_total, g2_total, g3_total)
complex(C_DOUBLE_COMPLEX),              dimension(:),    intent(in)  :: phiin
complex(C_DOUBLE_COMPLEX),              dimension(:),    intent(out) :: vi
real(dp),                                                intent(in)  :: kx_it
real(dp),                  allocatable, dimension(:,:)               :: vi_rhs
real(dp),                  allocatable, dimension(:)                 :: dvi, dlvi, duvi
integer                                                              :: j
integer,                              intent(in)  :: total_ny
real(dp),              dimension(:), intent(in)  :: g1_total, g2_total, g3_total

allocate(vi_rhs(total_ny-2,2), stat=alloc_err)
call check_alloc_err(alloc_err)
allocate(dvi(total_ny-2), stat=alloc_err)
call check_alloc_err(alloc_err)
allocate(dlvi(total_ny-3), duvi(total_ny-3), stat=alloc_err)
call check_alloc_err(alloc_err)

vi     = cmplx(0.0_dp, 0.0_dp, kind=C_DOUBLE_COMPLEX)
vi_rhs = 0.0_dp
dvi    = 0.0_dp
dlvi   = 0.0_dp
duvi   = 0.0_dp

do j = 2,total_ny-1
    dvi(j-1) = -kx_it**2.0_dp + g2_total(j)
end do

do j = 2,total_ny-2
    duvi(j-1) = g3_total(j)
end do

do j = 3,total_ny-1
    dlvi(j-2) = g1_total(j)
end do

vi_rhs(:,1) = real (phiin(2:total_ny-1))
vi_rhs(:,2) = aimag(phiin(2:total_ny-1))

call dgtsv(total_ny-2, 2, dlvi, dvi, duvi, vi_rhs, total_ny-2, info)

vi(1)   = cmplx(0.0_dp, 0.0_dp, kind=C_DOUBLE_COMPLEX)
vi(2:total_ny-1) = cmplx(vi_rhs(:,1), vi_rhs(:,2), kind=C_DOUBLE_COMPLEX)
vi(total_ny)  = cmplx(0.0_dp, 0.0_dp, kind=C_DOUBLE_COMPLEX)

end subroutine calc_vi_mod_MPI

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
subroutine update_bcs_mod_MPI(phiout,vout, phiin,vin,dyv1_T_it,dyv2_T_it,dyv1_B_it,&
                              dyv2_B_it,V1_in, V2_in, phi1_in, phi2_in,&
                              total_ny, h1_end, h2_end, h3_end)

integer,                                               intent(in)  :: total_ny
complex(C_DOUBLE_COMPLEX),              dimension(:),  intent(in)  :: phiin, vin
complex(C_DOUBLE_COMPLEX), allocatable, dimension(:),  intent(out) :: phiout, vout
real(dp),                               dimension(2,2)             :: C
real(dp)                                                           :: detC
complex(dp)                                                        :: c1, c2, c1t
complex(dp)                                                        :: dyV_T, dyV_B
real(dp),                                                intent(in)  :: dyv1_T_it,dyv2_T_it,dyv1_B_it,dyv2_B_it
real(dp),              dimension(:),  intent(in)  :: V1_in, V2_in, phi1_in, phi2_in
real(dp),                                              intent(in)  :: h1_end, h2_end, h3_end
allocate(phiout(total_ny), vout(total_ny), stat=alloc_err)
call check_alloc_err(alloc_err)
phiout = cmplx(0.0_dp, 0.0_dp, kind=C_DOUBLE_COMPLEX)
vout   = cmplx(0.0_dp, 0.0_dp, kind=C_DOUBLE_COMPLEX)
C      = 0.0_dp

C(1,1) = dyv1_T_it
C(1,2) = dyv2_T_it
C(2,1) = dyv1_B_it
C(2,2) = dyv2_B_it

detC = C(1,1)*C(2,2) - C(1,2)*C(2,1)

dyV_T = h1_end*vin(total_ny-2) + h2_end*vin(total_ny-1) + h3_end*vin(total_ny)
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
vout   = vin   + c1*V1_in   + c2*V2_in
phiout = phiin + c1*phi1_in + c2*phi2_in

end subroutine update_bcs_mod_MPI

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

subroutine calc_implicit_mod_MPI(Kphi,KT, phiin,Tin, kx_it, total_ny, dynu_in, g1_total, g2_total, g3_total)

complex(C_DOUBLE_COMPLEX),              dimension(:), intent(in)  :: phiin, Tin
complex(C_DOUBLE_COMPLEX),              dimension(:), intent(out) :: Kphi, KT
real(dp),                                               intent(in)  :: kx_it
integer,                                              intent(in)  :: total_ny
real(dp),   dimension(:),  intent(in)  :: dynu_in, g1_total, g2_total, g3_total                

Kphi = cmplx(0.0_dp, 0.0_dp, kind=C_DOUBLE_COMPLEX)
KT   = cmplx(0.0_dp, 0.0_dp, kind=C_DOUBLE_COMPLEX)

Kphi = nu0   *(-kx_it**2.0_dp*phiin + d2y_MPI(phiin, dynu_in, g1_total, g2_total, g3_total, total_ny))
KT   = kappa0*(-kx_it**2.0_dp*Tin   + d2y_MPI(Tin, dynu_in, g1_total, g2_total, g3_total, total_ny))

end subroutine calc_implicit_mod_MPI

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!

end module time_integrators_MPI
    
