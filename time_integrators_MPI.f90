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

integer                        :: nti, i, j, mpierror, total_ny
integer                        :: nprint, otherproc
logical                        :: wvtk
real(dp)                       :: nusselt_num
real(dp)                       :: start, finish
real(dp)                       :: start_overall, finish_overall
integer                        :: nthreads, myid, stind
complex(C_DOUBLE_COMPLEX), allocatable, dimension(:) :: tmp_phi_MPI, tmp_T_MPI, phi_MPI
complex(C_DOUBLE_COMPLEX), allocatable, dimension(:) :: K1hat_phi_MPI, K2hat_phi_MPI, K3hat_phi_MPI
complex(C_DOUBLE_COMPLEX), allocatable, dimension(:) :: K1hat_T_MPI, K2hat_T_MPI, K3hat_T_MPI
complex(C_DOUBLE_COMPLEX), allocatable, dimension(:) :: K1_phi_MPI, K2_phi_MPI, K1_T_MPI, K2_T_MPI
complex(C_DOUBLE_COMPLEX), allocatable, dimension(:) :: T_MPI, demo_T_MPI
complex(C_DOUBLE_COMPLEX), allocatable, dimension(:,:) :: phiout, Tout
real(dp), allocatable, dimension(:)   :: g1_total, g2_total, g3_total
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
   allocate(tmp_phi_MPI(total_ny), tmp_T_MPI(total_ny), phi_MPI(total_ny-2), stat=alloc_err)
   call check_alloc_err(alloc_err)
   allocate(K1hat_phi_MPI(total_ny-2), K2hat_phi_MPI(total_ny-2), K3hat_phi_MPI(total_ny-2), stat=alloc_err)
   call check_alloc_err(alloc_err)
   allocate(K1hat_T_MPI(total_ny-2), K2hat_T_MPI(total_ny-2), K3hat_T_MPI(total_ny-2), stat=alloc_err)
   call check_alloc_err(alloc_err)
   allocate(K1_phi_MPI(total_ny-2), K2_phi_MPI(total_ny-2), K1_T_MPI(total_ny-2), K2_T_MPI(total_ny-2), stat=alloc_err)
   call check_alloc_err(alloc_err)
   allocate(T_MPI(total_ny), stat=alloc_err)
   call check_alloc_err(alloc_err)
   allocate(g1_total(total_ny), g2_total(total_ny), g3_total(total_ny), stat=alloc_err)
   call check_alloc_err(alloc_err)
end if

allocate(phiout(Ny,Nx), Tout(Ny,Nx), stat=alloc_err)
call check_alloc_err(alloc_err)

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
    write(*,*) " - calc_explicit(1) timing: ", finish-start, "(s)"
    start = OMP_GET_WTIME()
    !$OMP PARALLEL DO private(tmp_phi, tmp_T, tmp_uy, tmp_phi1, tmp_uy1, tmp_K_phi, tmp_K_T) schedule(dynamic)
    do it = 1,1 ! kx loop
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

            ! Call vari mod
            call calc_vari_mod_MPI(tmp_phi_MPI, tmp_T_MPI, acoeffs(1,1), 1, total_ny,&
                                kx(it), phi_MPI,&
                                K1hat_phi_MPI,K2hat_phi_MPI,K3hat_phi_MPI,&
                                K1hat_T_MPI,K2hat_T_MPI,K3hat_T_MPI,&
                                K1_phi_MPI, K2_phi_MPI, K1_T_MPI, K2_T_MPI,&
                                g1_total, g2_total, g3_total,&
                                T_MPI)
            do i = 1,total_ny
                write(*,*) i, tmp_T_MPI(i)
            end do 
            
        else if (proc_id == num_procs - 1) then
            ! Send data to main node
            call MPI_SEND(phi(1,it), Ny-1, MPI_C_DOUBLE_COMPLEX, 0, 42, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(K1hat_phi(1,it), Ny-1, MPI_C_DOUBLE_COMPLEX, 0, 43, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(K2hat_phi(1,it), Ny-1, MPI_C_DOUBLE_COMPLEX, 0, 44, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(K3hat_phi(1,it), Ny-1, MPI_C_DOUBLE_COMPLEX, 0, 45, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(K1hat_T(1,it), Ny-1, MPI_C_DOUBLE_COMPLEX, 0, 46, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(K2hat_T(1,it), Ny-1, MPI_C_DOUBLE_COMPLEX, 0, 47, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(K3hat_T(1,it), Ny-1, MPI_C_DOUBLE_COMPLEX, 0, 48, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(K1_phi(1,it), Ny-1, MPI_C_DOUBLE_COMPLEX, 0, 49, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(K2_phi(1,it), Ny-1, MPI_C_DOUBLE_COMPLEX, 0, 50, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(K1_T(1,it), Ny-1, MPI_C_DOUBLE_COMPLEX, 0, 51, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(K2_T(1,it), Ny-1, MPI_C_DOUBLE_COMPLEX, 0, 52, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(T(1:Ny,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 53, MPI_COMM_WORLD, mpierror)

            ! Send g1, g2, g3 to first node.
            call MPI_SEND(g1(1), Ny, MPI_DOUBLE, 0, 54, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(g2(1), Ny, MPI_DOUBLE, 0, 55, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(g3(1), Ny, MPI_DOUBLE, 0, 56, MPI_COMM_WORLD, mpierror)
        else 
            ! Send data to main node
            call MPI_SEND(phi(1,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 42, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(K1hat_phi(1,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 43, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(K2hat_phi(1,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 44, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(K3hat_phi(1,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 45, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(K1hat_T(1,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 46, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(K2hat_T(1,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 47, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(K3hat_T(1,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 48, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(K1_phi(1,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 49, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(K2_phi(1,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 50, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(K1_T(1,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 51, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(K2_T(1,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 52, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(T(1,it), Ny, MPI_C_DOUBLE_COMPLEX, 0, 53, MPI_COMM_WORLD, mpierror)

            ! Send g1, g2, g3 to first node.
            call MPI_SEND(g1(1), Ny, MPI_DOUBLE, 0, 54, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(g2(1), Ny, MPI_DOUBLE, 0, 55, MPI_COMM_WORLD, mpierror)
            call MPI_SEND(g3(1), Ny, MPI_DOUBLE, 0, 56, MPI_COMM_WORLD, mpierror)
        end if
    end do
    !$OMP END PARALLEL DO
    finish = OMP_GET_WTIME()
    write(*,*) " - stage 1 mid timing: ", finish-start, "(s)"
    call MPI_BARRIER(MPI_COMM_WORLD, mpierror)
    ! if (proc_id == 0) then 
    !     open(unit=9010, file="P"//proc_id_str//"phiout_real.txt", action="write", status="unknown")
    !     open(unit=9011, file="P"//proc_id_str//"phiout_im.txt", action="write", status="unknown")
    !     do i=1,total_ny
    !         do j=1,Nx
    !             write (9010,*) REAL(Tout(i,j))
    !             write (9011,*) AIMAG(Tout(i,j))
    !         end do
    !     end do
    !     close(unit=9010)
    !     close(unit=9011)
    ! end if 
   write(*,*) "done writing phi!"
    if (time == t_final) then
        exit
    end if
    
 
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
write(*,*) " - - l1 timing: ", finish-start, "(s)"

start = OMP_GET_WTIME()
!$OMP PARALLEL DO schedule(dynamic)
do i=1,Nx
   ! Compute dx(T) in Fourier space
   nlT  (:,i) =  kx(i)*Ti(:,i)
   ! Compute D2(ux)
   nlphi(:,i) = -kx(i)**2.0_dp*uxi(:,i) + d2y(uxi(:,i))
end do
!$OMP END PARALLEL DO
finish = OMP_GET_WTIME()

write(*,*) " - - l2 timing: ", finish-start, "(s)"

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
write(*,*) " - - l3 timing: ", finish-start, "(s)"

! Calculate nonlinear term
start = OMP_GET_WTIME()
!$OMP PARALLEL DO private(tmp_T) schedule(dynamic)
do i = 1,Nx
   ! Temperature
   tmp_T = Ti(:,i)
   nlT(:,i) = uxi(:,i)*nlT(:,i) + uyi(:,i)*d1y(tmp_T)
   ! phi
   nlphi(:,i) = uxi(:,i)*phii(:,i) - uyi(:,i)*nlphi(:,i)
end do
!$OMP END PARALLEL DO
finish = OMP_GET_WTIME()
write(*,*) " - - l4 timing: ", finish-start, "(s)"

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
write(*,*) " - - l5 timing: ", finish-start, "(s)"

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
write(*,*) " - - l6 timing: ", finish-start, "(s)"

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

! do i = 2,total_ny-1
!     write(*,*) i, ddT (i-1)
! end do 

! write(*,*) kappa0 ,dt , aii


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
        FT(1)    = FT(1) + kappa0*dt*aii*g1(2)*T_in(1) ! b/c Ti(y_1) = T(y_1)
        FT(total_ny-2) = FT(total_ny-2) + kappa0*dt*aii*g3(total_ny-1)*T_in(total_ny) ! b/c Ti(total_ny) = T(total_ny)

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
        FT(1)    = FT(1) + kappa0*dt*aii*g1(2)*T_in(1) ! b/c Ti(y_1) = T(y_1)
        FT(total_ny-2) = FT(total_ny-2) + kappa0*dt*aii*g3(total_ny-1)*T_in(total_ny) ! b/c Ti(total_ny) = T(total_ny)

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

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!

end module time_integrators_MPI
    