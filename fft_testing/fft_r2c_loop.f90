program fft_r2c

use fftw

implicit none

integer :: alloc_err
integer, parameter :: dp = kind(1.0d0)
integer :: Nx
complex(dp), parameter :: CI = (0.0_dp,1.0_dp)
real(dp), allocatable :: xVals(:)
real(dp), allocatable :: kx(:)
real(dp) :: pi, xStep

integer :: i, j, n

type(C_PTR) :: tempPlan, itempPlan

real(dp), allocatable :: T(:)
real(dp), allocatable :: Tx(:)
complex(C_DOUBLE_COMPLEX), allocatable :: tempT(:)
complex(C_DOUBLE_COMPLEX), allocatable :: derivT(:)

pi = 4.0_dp*atan(1.0_dp)

do n=3,20

   Nx = 2**n
   write(*,*) "Nx = ", Nx

   allocate(T(Nx), Tx(Nx), derivT(Nx), kx(Nx), xVals(Nx), stat=alloc_err)

   !Compute wavenumbers
   do j=1,Nx/2
     kx(j) = real(j,kind=dp)-1.0_dp
   end do
   do j=Nx/2+1,Nx
     kx(j) = real(j-Nx,kind=dp) - 1.0_dp
   end do

   xStep = 2.0_dp * pi / real(Nx,kind=dp)
   do i=1,Nx
     xVals(i) = -pi + real(i-1, kind=dp) * xStep
   end do

   !T = 2.0_dp * sin(xVals) * cos(xVals)
   T = sin(xVals)

   !plans for taking T_x derivative
   tempPlan = fftw_plan_dft_r2c_1d(size(T),T,tempT,FFTW_ESTIMATE)
   itempPlan = fftw_plan_dft_c2r_1d(size(T),tempT,Tx,FFTW_ESTIMATE)

   allocate(tempT(Nx), stat=alloc_err)

   !go into Fourier space
   call fftw_execute_dft_r2c(tempPlan,T,tempT)
   tempT = tempT / real(Nx, kind=dp)
   tempT = CI * kx * tempT

   call fftw_execute_dft_c2r(itempPlan,tempT,Tx)

   open(unit=8000, file="errs_r2c.txt", action="write", status="unknown", position="append")
   write(8000, fmt=2000) Nx, maxval(abs(real(Tx - cos(xVals))))
   close(8000)

   call fftw_destroy_plan(tempPlan)
   call fftw_destroy_plan(itempPlan)

   deallocate(T, Tx, tempT, derivT, kx, xVals)

end do

2000 format(I7, E25.16E3)

end program fft_r2c
