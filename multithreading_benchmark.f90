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

integer :: iret
! integer :: n_iter, thread_iter, num_threads
! integer :: count, count_rate, count_max
N_test = 32768
! Create FFT plans
iret = fftw_init_threads()
! Uncomment below to check if init_threads working (should be nonzero)
! PRINT *, "iret: ", iret 

call fftw_plan_with_nthreads(2048)
plan = fftw_plan_dft_1d(N_test,tu_in,tu_out, FFTW_FORWARD,FFTW_ESTIMATE)
call global_allocations

call fftw_execute_dft(plan, tu_in, tu_out)

call fftw_destroy_plan(plan)

! Each output line contains results for an N
! results for each thread time for a single N comma separated

end program multithreading_benchmark