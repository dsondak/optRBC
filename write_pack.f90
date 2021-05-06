module write_pack

use global

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine write_vec(vname,place, vec)

character(len=*), intent(in) :: vname,place
real(dp), intent(in) :: vec(:)
integer :: ii, n

n = size(vec)

do ii = 1,n
   if (ii == 1) then
      write(*,*) vname//" = ", vec(ii)
   else
      write(*,*) place//"   ", vec(ii)
   end if
end do

end subroutine write_vec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine write_mat(vname,place, mat)

character(len=*), intent(in) :: vname,place
real(dp), intent(in)         :: mat(:,:)
integer                      :: jj, mx,my

mx = size(mat(1,:))
my = size(mat(:,1))

do jj = 1,my
   if (jj==1) then
      write(*,*) vname//" = ", mat(jj,:)
   else
      write(*,*) place//"   ", mat(jj,:)
   end if
end do

end subroutine write_mat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine write_vec_cmplx(vname,place,vec)

character(len=*), intent(in) :: vname, place
complex(dp), intent(in) :: vec(:)
integer :: ii, n

n = size(vec)

do ii = 1,n
   if (ii == 1) then
      write(*,*) vname//" = ", vec(ii)
   else
      write(*,*) place//"   ", vec(ii)
   end if
end do

end subroutine write_vec_cmplx

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine write_mat_cmplx(vname,place, mat)

character(len=*), intent(in) :: vname,place
complex(dp), intent(in) :: mat(:,:)
integer :: jj, mx,my

mx = size(mat(1,:))
my = size(mat(:,1))

do jj = 1,my
   if (jj==1) then
      write(*,*) vname//" = ", mat(jj,:)
   else
      write(*,*) place//"   ", mat(jj,:)
   end if
end do

end subroutine write_mat_cmplx

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine write_vtk_structured_grid(iter, process_id)

use global

implicit none

integer, intent(in) :: iter
character(10) :: citer
character(9)  :: fiter
integer :: n
integer :: ii, jj, kk
character(3), intent(in), optional  :: process_id

write(citer, "(I10)") 1000000000 + iter
fiter = citer(2:10)

n = (Nx+1)*Ny*Nz
if (PRESENT(process_id)) then
   open(unit=2, file=vtkloc//vtkname//"_P"//process_id//"_"//fiter//gtype, action="write", status="replace")
else
   open(unit=2, file=vtkloc//vtkname//fiter//gtype, action="write", status="replace")
endif
write(2, "(a)") "# vtk DataFile Version 3.1"
write(2, "(a)") "Thermal Convection Data"
write(2, "(a)") "ASCII"
write(2, "(a)") "DATASET STRUCTURED_GRID"
write(2, "(a, 3I5.1)") "DIMENSIONS", Nx+1, Ny, Nz
write(2, "(a, I15.1, a6)") "POINTS", n, "float"

do kk = 1,Nz
   do jj = 1,Ny
      do ii = 1,Nx+1
         if (ii < Nx+1) then
            write(2,"(3F18.14)") xp(ii), yp(jj), zp(kk)
         else
            write(2,"(3F18.14)") xR, yp(jj), zp(kk)
         end if
      end do
   end do
end do

write(2,"(a)") " "
write(2,"(a, a, I15.1)") "POINT_DATA", " ", n 
write(2,"(a, a, a, a, a)") "SCALARS", " ", "Temperature", " ", "float"
write(2,"(a, a, a)") "LOOKUP_TABLE", " ", "default"

do kk = 1,Nz
   do jj = 1,Ny
      do ii = 1,Nx+1
         if (ii < Nx+1) then
            write(2, "(F17.14)") real(T(jj,ii))
         else
            write(2, "(F17.14)") real(T(jj,1))
         end if
      end do
   end do
end do

write(2,"(a)") " "
write(2,"(a, a, a, a, a)") "SCALARS", " ", "ux", " ", "float"
write(2,"(a, a, a)") "LOOKUP_TABLE", " ", "default"

do kk = 1,Nz
   do jj = 1,Ny
      do ii = 1,Nx+1
         if (ii < Nx+1) then
            write(2, "(F17.14)") real(ux(jj,ii))
         else
            write(2, "(F17.14)") real(ux(jj,1))
         end if
      end do
   end do
end do

write(2,"(a)") " "
write(2,"(a, a, a, a, a)") "SCALARS", " ", "uy", " ", "float"
write(2,"(a, a, a)") "LOOKUP_TABLE", " ", "default"

do kk = 1,Nz
   do jj = 1,Ny
      do ii = 1,Nx+1
         if (ii < Nx+1) then
            write(2, "(F17.14)") real(uy(jj,ii))
         else
            write(2, "(F17.14)") real(uy(jj,1))
         end if
      end do
   end do
end do


close(unit=2)

end subroutine write_vtk_structured_grid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_to_vtk(step, physical, process_id)

integer, intent(in)            :: step
logical, intent(in)            :: physical
character(3), intent(in), optional  :: process_id
integer                        :: j


if (physical) then
   if (PRESENT(process_id)) then
      call write_vtk_structured_grid(step, process_id)
   else 
      call write_vtk_structured_grid(step)
   endif 
else
   do j = 1,Ny
      ! Bring everything to physical space
      tT  = T(j,:)
      tuy = uy(j,:)
      tux = ux(j,:)
      call fftw_execute_dft(iplanux, tux, tux)
      call fftw_execute_dft(iplanuy, tuy, tuy)
      call fftw_execute_dft(iplanT, tT, tT)
      T(j,:)   = tT
      ux(j,:)  = tux
      uy(j,:)  = tuy
   end do
   if (PRESENT(process_id)) then
      call write_vtk_structured_grid(step, process_id)
   else 
      call write_vtk_structured_grid(step)
   endif 
   do j = 1,Ny
      ! Bring everything to Fourier space
      tT  = T(j,:)
      tuy = uy(j,:)
      tux = ux(j,:)
      call fftw_execute_dft(planux, tux, tux)
      call fftw_execute_dft(planuy, tuy, tuy)
      call fftw_execute_dft(planT, tT, tT)
      T(j,:)   = tT
      ux(j,:)  = tux
      uy(j,:)  = tuy
   end do
   T  = T  / real(Nx, kind=dp)
   ux = ux / real(Nx, kind=dp)
   uy = uy / real(Nx, kind=dp)
end if

end subroutine write_to_vtk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine write_restart(physical)

integer :: j
logical, intent(in) :: physical

if (physical) then

   ! Write restart
   open(unit=2,file="uy_new",action="write", status="replace", form="unformatted")
   do ii = 1,Nx
      do jj = 1,Ny
         write(2) real(uy(jj,ii), kind=dp)
      end do
   end do
   close(unit=2)
   
   open(unit=3,file="temperature_new",action="write", status="replace", form="unformatted")
   do ii = 1,Nx
      do jj = 1,Ny
         write(3) real(T(jj,ii), kind=dp)
      end do
   end do
   close(unit=3)

else

   do j = 1,Ny
      ! Bring everything to physical space
      tT  = T(j,:)
      tuy = uy(j,:)
      call fftw_execute_dft(iplanuy, tuy, tuy)
      call fftw_execute_dft(iplanT, tT, tT)
      T(j,:)   = tT
      uy(j,:)  = tuy
   end do

   ! Write restart
   open(unit=2,file="uy_new",action="write", status="replace", form="unformatted")
   do ii = 1,Nx
      do jj = 1,Ny
         write(2) real(uy(jj,ii), kind=dp)
      end do
   end do
   close(unit=2)
   
   open(unit=3,file="temperature_new",action="write", status="replace", form="unformatted")
   do ii = 1,Nx
      do jj = 1,Ny
         write(3) real(T(jj,ii), kind=dp)
      end do
   end do
   close(unit=3)

   do j = 1,Ny
      ! Bring everything to Fourier space
      tT  = T(j,:)
      tuy = uy(j,:)
      call fftw_execute_dft(planuy, tuy, tuy)
      call fftw_execute_dft(planT, tT, tT)
      T(j,:)   = tT
      uy(j,:)  = tuy
   end do
   T  = T  / real(Nx, kind=dp)
   uy = uy / real(Nx, kind=dp)

end if

end subroutine write_restart

!::::::::::::::::::::::::::::::::::::::::::::
subroutine write_vtk_and_restart(step, physical)

integer, intent(in) :: step
logical, intent(in) :: physical
integer             :: j
character(10) :: citer
character(9)  :: fiter

write(citer, "(I10)") 1000000000 + step
fiter = citer(2:10)

if (physical) then

   ! Write to VTK
   call write_vtk_structured_grid(step)

   ! Write restart
   open(unit=2,file="uy"//fiter,action="write", status="replace", form="unformatted")
   do ii = 1,Nx
      do jj = 1,Ny
         write(2) real(uy(jj,ii), kind=dp)
      end do
   end do
   close(unit=2)
   
   open(unit=3,file="temperature"//fiter,action="write", status="replace", form="unformatted")
   do ii = 1,Nx
      do jj = 1,Ny
         write(3) real(T(jj,ii), kind=dp)
      end do
   end do
   close(unit=3)

else

   do j = 1,Ny
      ! Bring everything to physical space
      tT  = T(j,:)
      tuy = uy(j,:)
      tux = ux(j,:)
      call fftw_execute_dft(iplanux, tux, tux)
      call fftw_execute_dft(iplanuy, tuy, tuy)
      call fftw_execute_dft(iplanT, tT, tT)
      T(j,:)   = tT
      ux(j,:)  = tux
      uy(j,:)  = tuy
   end do

   ! Write to VTK
   call write_vtk_structured_grid(step)

   ! Write restart
   open(unit=2,file="uy"//fiter,action="write", status="replace", form="unformatted")
   do ii = 1,Nx
      do jj = 1,Ny
         write(2) real(uy(jj,ii), kind=dp)
      end do
   end do
   close(unit=2)
   
   open(unit=3,file="temperature"//fiter,action="write", status="replace", form="unformatted")
   do ii = 1,Nx
      do jj = 1,Ny
         write(3) real(T(jj,ii), kind=dp)
      end do
   end do
   close(unit=3)

   do j = 1,Ny
      ! Bring everything to Fourier space
      tT  = T(j,:)
      tuy = uy(j,:)
      tux = ux(j,:)
      call fftw_execute_dft(planux, tux, tux)
      call fftw_execute_dft(planuy, tuy, tuy)
      call fftw_execute_dft(planT, tT, tT)
      T(j,:)   = tT
      ux(j,:)  = tux
      uy(j,:)  = tuy
   end do
   T  = T  / real(Nx, kind=dp)
   ux = ux / real(Nx, kind=dp)
   uy = uy / real(Nx, kind=dp)

end if

end subroutine write_vtk_and_restart

end module write_pack
