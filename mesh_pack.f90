module mesh_pack

use global

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
integer                                                :: mpierror
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
