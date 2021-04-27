module mesh_pack

use global
include 'mpif.h'

contains

subroutine cosine_mesh(xcoor,ycoor,zcoor, numx,numy,numz,dxin)

implicit none

integer,                             intent(in)  :: numx, numy, numz
integer                                          :: ii, jj, kk
real(dp)                                         :: c
real(dp), optional,                  intent(in)  :: dxin
real(dp)                                         :: dxj
real(dp), allocatable, dimension(:), intent(out) :: xcoor, ycoor, zcoor

if (present(dxin)) then
   dxj = dxin
else
   dxj = dx
end if

allocate(xcoor(numx), stat=alloc_err)
allocate(ycoor(numy), stat=alloc_err)
allocate(zcoor(numz), stat=alloc_err)
call check_alloc_err(alloc_err)

! Create the grid.
xcoor(1) = xL
zcoor(1) = 0.0_dp
do jj = 1,numy
    c = (real(jj,kind=dp) - 1.0_dp) / (real(numy,kind=dp) - 1.0_dp)
    ycoor(jj) = -cos(c*pi)
end do
do ii = 2,numx
   xcoor(ii) = xcoor(ii-1) + dxj
end do
do kk = 2,numz
   zcoor(kk) = 0.0_dp
end do

end subroutine cosine_mesh
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

subroutine cosine_mesh_MPI(xcoor,ycoor,zcoor, numx,numy,numz, proc_id,num_procs,dxin)

implicit none

integer,                             intent(in)  :: numx, numy, numz
integer,                             intent(in)  :: proc_id, num_procs
integer                                          :: mpierror
integer                                          :: ii, jj, kk
real(dp)                                         :: c
real(dp), optional,                  intent(in)  :: dxin
real(dp)                                         :: dxj
real(dp), allocatable, dimension(:), intent(out) :: xcoor, ycoor, zcoor
integer                                          :: total_ny, start_ny

if (present(dxin)) then
   dxj = dxin
else
   dxj = dx
end if

allocate(xcoor(numx), stat=alloc_err)
allocate(ycoor(numy), stat=alloc_err)
allocate(zcoor(numz), stat=alloc_err)
call check_alloc_err(alloc_err)

! Create the grid.
xcoor(1) = xL
zcoor(1) = 0.0_dp
total_ny = numy * num_procs
start_ny = proc_id * numy
do jj = 1,numy
      c = (real(jj+start_ny,kind=dp) - 1.0_dp) / (real(total_ny,kind=dp) - 1.0_dp)
      ycoor(jj) = -cos(c*pi)
end do
do ii = 2,numx
   xcoor(ii) = xcoor(ii-1) + dxj
end do
do kk = 2,numz
   zcoor(kk) = 0.0_dp
end do

end subroutine cosine_mesh_MPI
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

subroutine dxdydz(dyj, xcoor,ycoor,zcoor)

implicit none

real(dp),              dimension(:), intent(in)  :: xcoor, ycoor, zcoor
real(dp), allocatable, dimension(:), intent(out) :: dyj
integer                                          :: alloc_err
integer                                          :: numy, jj

numy = size(ycoor)

allocate(dyj(numy-1), stat=alloc_err)
call check_alloc_err(alloc_err)

! Calculate grid spacing
do jj = 2,numy
   dyj(jj-1) = ycoor(jj) - ycoor(jj-1)
end do

end subroutine dxdydz

! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

subroutine dxdydz_MPI(dyj, xcoor,ycoor,zcoor, proc_id, num_procs)

implicit none

real(dp),              dimension(:), intent(in)  :: xcoor, ycoor, zcoor
integer,                             intent(in)  :: proc_id, num_procs
real(dp), allocatable, dimension(:), intent(out) :: dyj
integer                                          :: alloc_err, mpierror
integer                                          :: numy, jj
real(dp)                                         :: prev_element

numy = size(ycoor)

if (proc_id == 0) then
   allocate(dyj(numy-1), stat=alloc_err)
   call check_alloc_err(alloc_err)
   ! Calculate grid spacing
   do jj = 2,numy
      dyj(jj-1) = ycoor(jj) - ycoor(jj-1)
   end do
   ! Send last element to next process to the right.
   call MPI_SEND(ycoor(numy), 1, MPI_DOUBLE, proc_id+1, 42, MPI_COMM_WORLD, mpierror)
else if (proc_id == num_procs - 1) then
   ! Recieve first element from next process to the left.
   call MPI_RECV(prev_element, 1, MPI_DOUBLE, proc_id-1, 42, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
   allocate(dyj(numy), stat=alloc_err)
   call check_alloc_err(alloc_err)
   ! Calculate grid spacing
   dyj(1) = ycoor(1) - prev_element
   do jj = 2,numy
      dyj(jj) = ycoor(jj) - ycoor(jj-1)
   end do
else
   ! Recieve first element from next process to the left.
   call MPI_RECV(prev_element, 1, MPI_DOUBLE, proc_id-1, 42, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
   ! Send last element to next process to the right.
   call MPI_SEND(ycoor(numy), 1, MPI_DOUBLE, proc_id+1, 42, MPI_COMM_WORLD, mpierror)
   allocate(dyj(numy), stat=alloc_err)
   call check_alloc_err(alloc_err)
   ! Calculate grid spacing
   dyj(1) = ycoor(1) - prev_element
   do jj = 2,numy
      dyj(jj) = ycoor(jj) - ycoor(jj-1)
   end do
end if


end subroutine dxdydz_MPI
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
subroutine y_mesh_params

! Calculate metric coefficients
! g is for second derivatives
! h is for first derivatives

g1(1) =  2.0_dp / (dynu(1)*(dynu(1)+dynu(2)))
g2(1) = -2.0_dp / (dynu(1)*dynu(2))
g3(1) =  2.0_dp / (dynu(2)*(dynu(1)+dynu(2)))

h1(1) = -(2.0_dp*dynu(1)+dynu(2)) / (dynu(1)*(dynu(1)+dynu(2)))
h2(1) = (dynu(1) + dynu(2)) / (dynu(1)*dynu(2))
h3(1) = -dynu(1) / (dynu(2)*(dynu(1)+dynu(2)))
do jj = 2,Ny-1
   g1(jj) =  2.0_dp / (dynu(jj-1)*(dynu(jj)+dynu(jj-1)))
   g2(jj) = -2.0_dp / (dynu(jj-1)*dynu(jj))
   g3(jj) =  2.0_dp / (dynu( jj )*(dynu(jj)+dynu(jj-1)))

   h1(jj) = -dynu(jj) / (dynu(jj-1)*(dynu(jj-1) + dynu(jj)))
   h2(jj) = (dynu(jj) - dynu(jj-1)) / (dynu(jj-1)*dynu(jj))
   h3(jj) = dynu(jj-1) / (dynu(jj)*(dynu(jj-1) + dynu(jj)))
end do
h1(Ny) = dynu(Ny-1) / (dynu(Ny-2)*(dynu(Ny-2)+dynu(Ny-1)))
h2(Ny) = -(dynu(Ny-2) + dynu(Ny-1)) / (dynu(Ny-2)*dynu(Ny-1))
h3(Ny) = (2.0_dp*dynu(Ny-1) + dynu(Ny-2)) / (dynu(Ny-1)*(dynu(Ny-2) + dynu(Ny-1)))

g1(Ny) =  2.0_dp / (dynu(Ny-2)*(dynu(Ny-2) + dynu(Ny-1)))
g2(Ny) = -2.0_dp / (dynu(Ny-1)*dynu(Ny-2))
g3(Ny) =  2.0_dp / (dynu(Ny-1)*(dynu(Ny-1) + dynu(Ny-2)))

end subroutine y_mesh_params
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
subroutine y_mesh_params_MPI (proc_id, num_procs)

implicit none

integer, intent(in)  :: proc_id, num_procs
real(dp)             :: recv_val
integer              :: mpierror, jj

! Calculate metric coefficients
! g is for second derivatives
! h is for first derivatives

! First element
if (proc_id == 0) then
   call MPI_RECV(recv_val, 1, MPI_DOUBLE, proc_id+1, 42, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)

   ! First elems
   g1(1) =  2.0_dp / (dynu(1)*(dynu(1)+dynu(2)))
   g2(1) = -2.0_dp / (dynu(1)*dynu(2))
   g3(1) =  2.0_dp / (dynu(2)*(dynu(1)+dynu(2)))

   h1(1) = -(2.0_dp*dynu(1)+dynu(2)) / (dynu(1)*(dynu(1)+dynu(2)))
   h2(1) = (dynu(1) + dynu(2)) / (dynu(1)*dynu(2))
   h3(1) = -dynu(1) / (dynu(2)*(dynu(1)+dynu(2)))

   ! Middle Ny-2 elements.
   do jj = 2,Ny-1
      g1(jj) =  2.0_dp / (dynu(jj-1)*(dynu(jj)+dynu(jj-1)))
      g2(jj) = -2.0_dp / (dynu(jj-1)*dynu(jj))
      g3(jj) =  2.0_dp / (dynu( jj )*(dynu(jj)+dynu(jj-1)))

      h1(jj) = -dynu(jj) / (dynu(jj-1)*(dynu(jj-1) + dynu(jj)))
      h2(jj) = (dynu(jj) - dynu(jj-1)) / (dynu(jj-1)*dynu(jj))
      h3(jj) = dynu(jj-1) / (dynu(jj)*(dynu(jj-1) + dynu(jj)))
   end do

   ! Last elems
   g1(Ny) =  2.0_dp / (dynu(Ny-1)*(recv_val+dynu(Ny-1)))
   g2(Ny) = -2.0_dp / (dynu(Ny-1)*recv_val)
   g3(Ny) =  2.0_dp / (recv_val*(recv_val+dynu(Ny-1)))

   h1(Ny) = -recv_val / (dynu(Ny-1)*(dynu(Ny-1) + recv_val))
   h2(Ny) = (recv_val - dynu(Ny-1)) / (dynu(Ny-1)*recv_val)
   h3(Ny) = dynu(Ny-1) / (recv_val*(dynu(Ny-1) + recv_val))
else if (proc_id == num_procs - 1) then
   call MPI_SEND(dynu(1), 1, MPI_DOUBLE, proc_id-1, 42, MPI_COMM_WORLD, mpierror)
   do jj = 1,Ny-1
      g1(jj) =  2.0_dp / (dynu(jj)*(dynu(jj+1)+dynu(jj)))
      g2(jj) = -2.0_dp / (dynu(jj)*dynu(jj+1))
      g3(jj) =  2.0_dp / (dynu(jj+1)*(dynu(jj+1)+dynu(jj)))

      h1(jj) = -dynu(jj+1) / (dynu(jj)*(dynu(jj) + dynu(jj+1)))
      h2(jj) = (dynu(jj+1) - dynu(jj)) / (dynu(jj)*dynu(jj+1))
      h3(jj) = dynu(jj) / (dynu(jj+1)*(dynu(jj) + dynu(jj+1)))
   end do
   h1(Ny) = dynu(Ny) / (dynu(Ny-1)*(dynu(Ny-1)+dynu(Ny)))
   h2(Ny) = -(dynu(Ny-1) + dynu(Ny)) / (dynu(Ny-1)*dynu(Ny))
   h3(Ny) = (2.0_dp*dynu(Ny) + dynu(Ny-1)) / (dynu(Ny)*(dynu(Ny-1) + dynu(Ny)))

   g1(Ny) =  2.0_dp / (dynu(Ny-1)*(dynu(Ny-1) + dynu(Ny)))
   g2(Ny) = -2.0_dp / (dynu(Ny)*dynu(Ny-1))
   g3(Ny) =  2.0_dp / (dynu(Ny)*(dynu(Ny) + dynu(Ny-1)))

else
   call MPI_SEND(dynu(1), 1, MPI_DOUBLE, proc_id-1, 42, MPI_COMM_WORLD, mpierror)
   call MPI_RECV(recv_val, 1, MPI_DOUBLE, proc_id+1, 42, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpierror)
   ! All but last elems.
   do jj = 1,Ny-1
      g1(jj) =  2.0_dp / (dynu(jj)*(dynu(jj+1)+dynu(jj)))
      g2(jj) = -2.0_dp / (dynu(jj)*dynu(jj+1))
      g3(jj) =  2.0_dp / (dynu(jj+1)*(dynu(jj+1)+dynu(jj)))

      h1(jj) = -dynu(jj+1) / (dynu(jj)*(dynu(jj) + dynu(jj+1)))
      h2(jj) = (dynu(jj+1) - dynu(jj)) / (dynu(jj)*dynu(jj+1))
      h3(jj) = dynu(jj) / (dynu(jj+1)*(dynu(jj) + dynu(jj+1)))
   end do
   ! Last elems
   g1(Ny) =  2.0_dp / (dynu(Ny)*(recv_val+dynu(Ny)))
   g2(Ny) = -2.0_dp / (dynu(Ny)*recv_val)
   g3(Ny) =  2.0_dp / (recv_val*(recv_val+dynu(Ny)))

   h1(Ny) = -recv_val / (dynu(Ny)*(dynu(Ny) + recv_val))
   h2(Ny) = (recv_val - dynu(Ny)) / (dynu(Ny)*recv_val)
   h3(Ny) = dynu(Ny) / (recv_val*(dynu(Ny) + recv_val))
end if 

end subroutine y_mesh_params_MPI
! :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

subroutine uniform_mesh(xcoor,ycoor,zcoor, numx,numy,numz)

implicit none

integer,                             intent(in)  :: numx, numy, numz
integer                                          :: ii, jj, kk
real(dp)                                         :: c
real(dp), allocatable, dimension(:), intent(out) :: xcoor, ycoor, zcoor

allocate(xcoor(numx), stat=alloc_err)
allocate(ycoor(numy), stat=alloc_err)
allocate(zcoor(numz), stat=alloc_err)
call check_alloc_err(alloc_err)

! Create the grid.
xcoor(1) = xL
ycoor(1) = ybot
zcoor(1) = 0.0_dp
do jj = 1,Ny
    ycoor(jj) = real((jj-1), kind=dp)*dy + ybot 
end do
do ii = 2,numx
   xcoor(ii) = xcoor(ii-1) + dx
end do
do kk = 2,numz
   zcoor(kk) = 0.0_dp
end do

end subroutine uniform_mesh

end module mesh_pack
