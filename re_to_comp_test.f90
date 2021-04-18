program re_to_comp_test

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
    
    integer :: N_points
    ! go from real -> complex -> real
    real(C_DOUBLE), allocatable :: tu_in_real(:), tu_out_real(:)
    complex(C_DOUBLE_COMPLEX), allocatable :: tu_comp(:)
    real :: start, finish
    type(C_PTR) :: planre, planc
    N_points = 2048

    planre =  fftw_plan_dft_r2c_1d(N_points, tu_in_real, tu_comp, FFTW_FORWARD);
    planc = fftw_plan_dft_c2r_1d(N_points, tu_comp, tu_out_real, FFTW_BACKWARD);

    ! Allocate space for each array
    allocate(tu_in_real(N_points), stat=alloc_err)
    allocate(tu_out_real(N_points), stat=alloc_err)
    allocate(tu_comp(N_points/2 + 1), stat=alloc_err)
    tu_in_real = (0.0_dp, 0.0_dp)

    call cpu_time(start)
    call fftw_execute_dft_r2c(planre, tu_in_real, tu_comp)
    call fftw_execute_dft_c2r(planc, tu_comp, tu_out_real)
    call cpu_time(finish)

    call fftw_destroy_plan(planre)
    call fftw_destroy_plan(planc)

    deallocate(tu_in_real, tu_out_real)
    deallocate(tu_comp)
    
end program re_to_comp_test                