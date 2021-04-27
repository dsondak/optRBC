module bc_setup

use global
use write_pack

implicit none
include 'mpif.h'

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine init_bc_MPI(aii, proc_id, num_procs, proc_id_str)

implicit none

real(dp),              intent(in)     :: aii
integer,               intent(in)     :: proc_id, num_procs
character(3),  optional,        intent(in)     :: proc_id_str
real(dp)                              :: pnu, qnu
real(dp)                              :: wavex

real(dp), allocatable, dimension(:)   :: d, dl, du
real(dp), allocatable, dimension(:)   :: phi1_b, phi1_t
real(dp), allocatable, dimension(:)   :: g1_total, g2_total, g3_total
real(dp), allocatable, dimension(:)   :: recv_col, recv_col2, last_row_recv
real(dp), allocatable, dimension(:)   :: phi2_b, phi2_t
real(dp), allocatable, dimension(:,:) :: Fphi12, FV12


integer                               :: ii, jj, i, total_ny, otherproc, start
integer                               :: info, mpierror

! Only do allocation on main node.
if(proc_id == 0) then
   total_ny = Ny * num_procs

   ! Allocate local variables
   allocate(Fphi12(total_ny-2,2), FV12(total_ny-2,2), stat=alloc_err)
   call check_alloc_err(alloc_err)
   allocate(d(total_ny-2), stat=alloc_err)
   call check_alloc_err(alloc_err)
   allocate(dl(total_ny-3), stat=alloc_err)
   call check_alloc_err(alloc_err)
   allocate(du(total_ny-3), stat=alloc_err)
   call check_alloc_err(alloc_err)
   allocate(phi1_b(Nx), phi1_t(Nx), stat=alloc_err)
   call check_alloc_err(alloc_err)
   allocate(phi2_b(Nx), phi2_t(Nx), stat=alloc_err)
   call check_alloc_err(alloc_err)

   ! Need to pull all of g1, g2, g3 onto main node for calculations.
   allocate(g1_total(total_ny), g2_total(total_ny), g3_total(total_ny), stat=alloc_err)
   call check_alloc_err(alloc_err)

   ! Receive g1,g2,g3 from all other nodes.
   do otherproc = 1,num_procs-1
      start = otherproc * Ny + 1
      call MPI_RECV(g1_total(start), Ny, MPI_DOUBLE, otherproc, 42, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
      call MPI_RECV(g2_total(start), Ny, MPI_DOUBLE, otherproc, 43, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
      call MPI_RECV(g3_total(start), Ny, MPI_DOUBLE, otherproc, 44, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
   end do
   g1_total(1:Ny) = g1
   g2_total(1:Ny) = g2
   g3_total(1:Ny) = g3

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
   
      do jj = 2,total_ny-1
         d(jj-1) = qnu - pnu*g2_total(jj)
      end do
   
      do jj = 2,total_ny-2
         du(jj-1) = -pnu*g3_total(jj)
      end do
   
      do jj = 3,total_ny-1
         dl(jj-2) = -pnu*g1_total(jj)
      end do

      Fphi12(1,1)    = pnu*g1_total(2)*phi1_b(ii)
      Fphi12(total_ny-2,2) = pnu*g3_total(total_ny-1)*phi2_t(ii)

      !  Solve the system Aphi phi = Fphi
      call dgtsv(total_ny-2, 2, dl, d, du, Fphi12, total_ny-2, info)

      !  Put phi1 and phi2 together in main node.
      phi1(1,ii)      = phi1_b(ii)
      phi1(2:Ny,ii) = Fphi12(1:Ny-1,1)
      phi2(1,ii)      = 0.0_dp
      phi2(2:Ny,ii) = Fphi12(1:Ny-1,2)
      
      ! Send phi1 and phi2 columns to each node.
      do otherproc = 1,num_procs-1
         start = otherproc * Ny
         call MPI_SEND(Fphi12(start,1), Ny, MPI_DOUBLE, otherproc, 45, MPI_COMM_WORLD, mpierror)
         call MPI_SEND(Fphi12(start,2), Ny, MPI_DOUBLE, otherproc, 46, MPI_COMM_WORLD, mpierror)
      end do

      ! Getting V1 and V2 initialized.
      FV12(:,1) = Fphi12(:,1)
      FV12(:,2) = Fphi12(:,2)

      do jj = 2,total_ny-1
         d(jj-1) = -wavex**2.0_dp + g2_total(jj)
      end do

      do jj = 2,total_ny-2
         du(jj-1) = g3_total(jj)
      end do

      do jj = 3,total_ny-1
         dl(jj-2) = g1_total(jj)
      end do

      call dgtsv(total_ny-2, 2, dl, d, du, FV12, total_ny-2, info)

      ! Place V1, V2 in main node memory.
      V1(2:Ny,ii) = FV12(1:Ny-1,1)
      V2(2:Ny,ii) = FV12(1:Ny-1,2)

      ! Send V1 and V1 columns to each node.
      do otherproc = 1,num_procs-1
         start = otherproc * Ny
         call MPI_SEND(FV12(start,1), Ny, MPI_DOUBLE, otherproc, 48, MPI_COMM_WORLD, mpierror)
         call MPI_SEND(FV12(start,2), Ny, MPI_DOUBLE, otherproc, 49, MPI_COMM_WORLD, mpierror)
      end do

   end do

   ! Send last row to final proc.
   call MPI_SEND(phi2_t(1), Nx, MPI_DOUBLE, num_procs-1, 47, MPI_COMM_WORLD, mpierror)

else
   ! Send g1, g2, g3 to first node.
   call MPI_SEND(g1(1), 12, MPI_DOUBLE, 0, 42, MPI_COMM_WORLD, mpierror)
   call MPI_SEND(g2(1), 12, MPI_DOUBLE, 0, 43, MPI_COMM_WORLD, mpierror)
   call MPI_SEND(g3(1), 12, MPI_DOUBLE, 0, 44, MPI_COMM_WORLD, mpierror)

   allocate(recv_col(Ny), stat=alloc_err)
   call check_alloc_err(alloc_err)
   allocate(recv_col2(Ny), stat=alloc_err)
   call check_alloc_err(alloc_err)

   ! Receive cols of phi1 and phi2.
   do ii = 1,Nx
      call MPI_RECV(recv_col, Ny, MPI_DOUBLE, 0, 45, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
      call MPI_RECV(recv_col2, Ny, MPI_DOUBLE, 0, 46, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
      phi1(1:Ny, ii) = recv_col
      phi2(1:Ny, ii) = recv_col2
   end do

   ! Receive cols of V1 and V2.
   do ii = 1,Nx
      call MPI_RECV(recv_col, Ny, MPI_DOUBLE, 0, 48, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
      call MPI_RECV(recv_col2, Ny, MPI_DOUBLE, 0, 49, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
      V1(1:Ny, ii) = recv_col
      V2(1:Ny, ii) = recv_col2
   end do

   ! Handle boundary of last processor.
   if (proc_id == num_procs - 1) then
      allocate(last_row_recv(Nx), stat=alloc_err)
      call check_alloc_err(alloc_err)
      call MPI_RECV(last_row_recv, Nx, MPI_DOUBLE, 0, 47, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
      phi1(Ny,1:Nx)     = 0.0_dp
      phi2(Ny,1:Nx)     = last_row_recv
      V1(Ny, 1:Nx)      = 0.0_dp
      V2(Ny, 1:Nx)      = 0.0_dp
   end if 
end if 

! Calculate wall derivative on first and last processor.
if (proc_id == 0) then
   do ii = 1,Nx
      dyv1_B(ii) = h1(1)*V1(1,ii) + h2(1)*V1(2,ii) + h3(1)*V1(3,ii)
      dyv2_B(ii) = h1(1)*V2(1,ii) + h2(1)*V2(2,ii) + h3(1)*V2(3,ii)
   end do
else if (proc_id == num_procs - 1) then
   do ii = 1,Nx
      dyv1_T(ii) = h1(Ny)*V1(Ny-2,ii) + h2(Ny)*V1(Ny-1,ii) + h3(Ny)*V1(Ny,ii)
      dyv2_T(ii) = h1(Ny)*V2(Ny-2,ii) + h2(Ny)*V2(Ny-1,ii) + h3(Ny)*V2(Ny,ii)
   end do
end if

end subroutine init_bc_MPI

end module bc_setup
