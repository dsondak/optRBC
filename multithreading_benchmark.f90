program multithreading_benchmark

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use global
    use strings
    use write_pack
    use interpolation_pack
    use allocate_vars
    use statistics
    use mesh_pack
    use time_integrators

    implicit none

! http://www.fftw.org/fftw3_doc/Fortran-Examples.html

integer :: N_test
integer :: iret
integer :: N_threads, N_exp, curr_threads
complex(C_DOUBLE_COMPLEX), allocatable :: tu_in(:), tu_out(:)
real :: start, finish

! Create FFT plans
iret = fftw_init_threads()
! Uncomment below to check if init_threads working (should be nonzero)
PRINT *, "iret: ", iret 

do N_exp = 7,20
    N_test = 2**N_exp
    print '("Points: ",i12)', N_test
    curr_threads = 2
    do N_threads = 0,10
        curr_threads = 2**N_threads
        call fftw_plan_with_nthreads(curr_threads)
        plan = fftw_plan_dft_1d(N_test,tu_in,tu_out, FFTW_FORWARD,FFTW_ESTIMATE)

        allocate(tu_in(N_test), stat=alloc_err)
        allocate(tu_out(N_test), stat=alloc_err)
        tu_in    = (0.0_dp, 0.0_dp)
        tu_out    = (0.0_dp, 0.0_dp)

        call cpu_time(start)
        call fftw_execute_dft(plan, tu_in, tu_out)
        call cpu_time(finish)
        print '("Threads = ",i5,"   Time = ",f6.3," seconds.")',curr_threads,finish-start
        call fftw_destroy_plan(plan)
        deallocate(tu_in, tu_out)
    end do
end do
! Each output line contains results for an N
! results for each thread time for a single N comma separated

end program multithreading_benchmark