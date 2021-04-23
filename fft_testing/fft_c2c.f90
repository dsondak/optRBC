program fft_c2c

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
complex(C_DOUBLE_COMPLEX), allocatable :: ckx(:)
real(dp) :: temp1
real(dp) :: maxError

integer :: i,j, maxIndex

type(C_PTR) :: tempPlan, itempPlan

complex(C_DOUBLE_COMPLEX), allocatable :: T(:)
complex(C_DOUBLE_COMPLEX), allocatable ::  Tx(:)
complex(C_DOUBLE_COMPLEX), allocatable :: tempT(:)
complex(C_DOUBLE_COMPLEX), allocatable :: derivT(:)



allocate(T(Nx),stat=alloc_err)
allocate(ckx(Nx),stat=alloc_err)


pi = 4.0_dp*atan(1.0_dp)

write(*,*) "----------------"

write(*,*) "Nx: ", Nx

write(*,*) "----------------"

!!Compute wavenumbers
do j=1,Nx/2
  kx(j) = real(j,kind=dp)-1.0_dp
  ckx(j) = cmplx(kx(j),kind=C_DOUBLE_COMPLEX)
end do
do j=Nx/2+1,Nx
  kx(j) = real(j-Nx,kind=dp) - 1.0_dp
  ckx(j) = cmplx(kx(j),kind=C_DOUBLE_COMPLEX)
end do

!write(*,*) "wavenumbes:"
!write(*,*) kx
!write(*,*) "----------"

!plans for taking T_x derivative
tempPlan = fftw_plan_dft_1d(Nx,tempT,tempT,FFTW_FORWARD,FFTW_ESTIMATE)
itempPlan = fftw_plan_dft_1d(Nx,tempT,tempT,FFTW_BACKWARD,FFTW_ESTIMATE)


call makeXArray(xVals,Nx)

do j=1,Nx
  T(j) = 2*sin(xVals(j))*cos(xVals(j))
  !T(j) = sin(2*xVals(j))
end do

!write(*,*) "initial T array"

!write(*,*) T

!write(*,*) "--------------"

tempT = T

!write(*,*) "tempT:"

!write(*,*) tempT

!write(*,*) "---------"

!!!!!!!!!!!!!!!!!!!
!!CAlculating T_x!!
!!!!!!!!!!!!!!!!!!!

!go into Fourier space
call fftw_execute_dft(tempPlan,tempT,tempT)

tempT = tempT/Nx

!write(*,*) "In Fourier space"

!write(*,*) tempT

!write(*,*) "------------"

!Compute derivative in Fourier space

!tempT = CI*ckx*tempT

!bring back to physical space

call fftw_execute_dft(itempPlan,tempT,tempT)

Tx = tempT

!write(*,*) "Tx:"

!write(*,*) Tx

!write(*,*) "-------------"

write(*,*) "error in Tx derivative:"

write(*,*) maxval(abs(real(Tx-2.0_dp*sin(xVals)*cos(xVals))))
!write(*,*) maxval(abs(real(Tx-2*cos(2*xVals))))
write(*,*) "-------------"

maxIndex = 0
maxError = 0

do i=1,Nx
  if (abs(Tx(i)-2*cos(2*xVals(i))) > maxError) then
    maxIndex = i
    maxError = abs(Tx(i)-2*cos(2*xVals(i)))
  end if
end do

write(*,*) "----------------"
write(*,*) "max error is at x=", xVals(maxIndex)
write(*,*) "----------------"

write(*,*) "----------------"
write(*,*) "Tx value at maxIndex is:" 
write(*,*) Tx(maxIndex)
write(*,*) "-------------------"

open(unit=1,file="deriv_c2c.txt",action="write",status="unknown")
do i=1,Nx
  write(1,*) xVals(i), Tx(i)
end do
close(1)

end program fft_c2c

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


