module imod

implicit none

abstract interface
   function func(x) result(jact)

      use global
      implicit none

      real(kind(1.0d0)), intent(in) :: x(:)
      real(kind(1.0d0)), allocatable :: jact(:)
      real(kind(1.0d0)), allocatable :: x0_delta(:), xT_delta(:), GT_delta(:)
      real(kind(1.0d0)) :: normx, dot_prod, eps
      real(kind(1.0d0)) :: dnrm2, ddot
      integer :: n
      integer, parameter :: incx = 1, incy = 1

   end function func
end interface

end module imod
