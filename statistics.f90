module statistics

use global

implicit none

contains

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
subroutine energy(k,u)

use global

implicit none

complex(dp), intent(in)  :: u(:,:)
real(dp),    intent(out) :: k
integer                  :: ii, jj

k = 0.0_dp
do ii = 1,Nx
   do jj = 2,Ny
      k = k  + conjg(u(jj,ii))*u(jj,ii) * dynu(jj-1)
   end do
end do
k = 0.5_dp * alpha * k
! Note:  Multiplied by alpha for the following reason:
!   k =0.5* \sum u*u \Delta kx \Delta y
!   but kx = \alpha l where l is the integer wavenumber.
!   So \Delta kx = \alpha l1 - \alpha l2
!                = \alpha (l1 - l2)
!                = \alpha b/c l1 - l2 = 1.
!   Thus \Delta kx is uniform and is equal to \alpha.

end subroutine energy
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
subroutine Espectra(E, u)

use global
use write_pack

implicit none
complex(dp),              intent(in)  :: u(:,:)
real(dp),    allocatable, intent(out) :: E(:)
real(dp)                              :: ubar
integer                               :: ii, jj
integer                               :: kmag

allocate(E(kmax), stat=alloc_err)

E    = 0.0_dp
ubar = 0.0_dp

do ii = 1,Nx
   kmag = abs(nint(kx(ii) / alpha))
   do jj = 2,Ny
      ubar = ubar + conjg(u(jj,ii))*u(jj,ii)*dynu(jj-1)
   end do
   ubar = ubar / Ly
   if (kmag < kmax) then
      E(int(kmag)+1) = E(int(kmag)+1) + 0.5_dp*ubar
   end if
end do

end subroutine Espectra

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::!
subroutine track_max(mags,u)

use global

implicit none

complex(dp), intent(in)  :: u(:,:)
real(dp), intent(out)    :: mags
real(dp)                 :: mag
integer                  :: ii, jj
integer                  :: imax, jmax

mag = 0.0_dp
mags = 0.0_dp

imax = size(u, 2)
jmax = size(u, 1)

do ii = 1,imax
   do jj = 1,jmax
      mag = sqrt(real(u(jj,ii))**2.0_dp + &
                &aimag(u(jj,ii))**2.0_dp)
      if (mag > mags) then
         mags = mag
      end if
   end do
end do

end subroutine track_max

subroutine nusselt(Nuave, in_fourier)

use global
use write_pack

implicit none

real(dp), intent(out) :: Nuave
logical, intent(in) :: in_fourier
real(dp) :: dyT_wall
integer :: ii

write(*,*) "Computing Nu..."


! Compute Nu at walls!
if (in_fourier) then
   dyT_wall = h1(1) * real(T(1,1)) + h2(1) * real(T(2,1)) + h3(1) * real(T(3,1))
   Nuave = -dyT_wall
else
   Nuave = 0.0_dp
   do ii = 1,Nx
      ! Compute temperature gradients at the walls
      dyT_wall = h1(1)*real(T(1,ii)) + h2(1)*real(T(2,ii)) + h3(1)*real(T(3,ii))
      Nuave = Nuave - dyT_wall ! Nu = -dy(T)
   end do
   Nuave = Nuave / Nx
end if
write(*,*) "Nu = ", Nuave

end subroutine nusselt

subroutine nusselt_Ti(Nuave, in_fourier)

   use global
   use write_pack
   
   implicit none
   
   real(dp), intent(out) :: Nuave
   logical, intent(in) :: in_fourier
   real(dp) :: dyT_wall
   integer :: ii
   
   write(*,*) "Computing Nu (Ti)..."
   
   
   ! Compute Nu at walls!
   if (in_fourier) then
      dyT_wall = h1(1) * real(T(1,1)) + h2(1) * real(T(2,1)) + h3(1) * real(T(3,1))
      Nuave = -dyT_wall
   else
      Nuave = 0.0_dp
      do ii = 1,Nx
         ! Compute temperature gradients at the walls
         dyT_wall = h1(1)*real(T(1,ii)) + h2(1)*real(T(2,ii)) + h3(1)*real(T(3,ii))
         Nuave = Nuave - dyT_wall ! Nu = -dy(T)
      end do
      Nuave = Nuave / Nx
   end if
   write(*,*) "Nu (Ti) = ", Nuave
   
end subroutine nusselt_Ti

end module statistics
