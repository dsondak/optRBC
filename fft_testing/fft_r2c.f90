program fft_r2c

use fftw

implicit none

integer :: alloc_err
integer, parameter :: dp = kind(1.0d0)
integer, parameter :: Nx=1048576, Ny=10
complex(dp), parameter :: CI = (0.0_dp,1.0_dp)
real(dp), dimension(Nx) :: xVals
real(dp), dimension(Ny) :: yVals
!real(dp), dimension(Nx) :: derivT
real(dp) :: pi
real(dp), dimension(Nx) :: kx
real(dp) :: temp1
real(dp) :: maxError

real(dp), dimension(Nx) :: T
real(dp), dimension(Nx) :: Tx


integer :: i,j, maxIndex

type(C_PTR) :: tempPlan, itempPlan

complex(C_DOUBLE_COMPLEX), allocatable :: tempT(:)

allocate(tempT(Nx), stat=alloc_err)

write(*,*) "Nx: ", Nx

do j=1,Nx/2
  kx(j) = real(j,kind=dp)-1.0_dp
end do
do j=Nx/2+1,Nx
  kx(j) = real(j-Nx,kind=dp) - 1.0_dp
end do


call makeXArray(xVals,Nx)

tempT = (0.0_dp,0.0_dp)
T = 2.0_dp*sin(xVals)*cos(xVals)


!write(*,*) "---------------"

!write(*,*) "initial T array:", T

!write(*,*) "---------------"

!plans for taking Tx derivative
tempPlan = fftw_plan_dft_r2c_1d(size(T),T,tempT,FFTW_ESTIMATE)
itempPlan = fftw_plan_dft_c2r_1d(size(T),tempT,Tx,FFTW_ESTIMATE)


!write(*,*) "about to call r2c dft"
call fftw_execute_dft_r2c(tempPlan,T,tempT)
!write(*,*) "r2c dft complete"


tempT = tempT/Nx

!tempT = CI*kx*tempT

!write(*,*) "about to call c2r dft"

call fftw_execute_dft_c2r(itempPlan,tempT,Tx)

!write(*,*) "---------------"    

!write(*,*) "final T array:", Tx

!write(*,*) "---------------" 

!open(unit=1,file="deriv_r2c.txt",action="write",status="unknown")
!do i=1,Nx
!  write(1,*) xVals(i), Tx(i)
!end do
!close(1)

write(*,*) "error: ", maxval(abs(Tx-2.0_dp*sin(xVals)*cos(xVals)))

!write(*,*) "error: ", maxval(abs(Tx-2.0_dp*cos(2.0_dp*xVals)))

end program fft_r2c

subroutine makeXArray(xArray,N)
        implicit none
        integer, parameter :: dp=kind(1.0d0)
        integer,  intent(in) :: N
        real(dp) :: xStep
        integer :: i
        real(dp), dimension(N), intent(out) :: xArray
        real(dp), dimension(N) :: x1Array
        real(dp) :: pi
        pi = 4.0_dp*atan(1.0_dp)

        xStep = 2.0_dp*pi/real(N,kind=dp)
        do i=1,N
          xArray(i) = 0_dp
        end do
        do i=1,N
          xArray(i) = -pi + real(i-1,kind=dp)*xStep
        end do
end subroutine makeXArray
