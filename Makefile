FC=mpif90
FFLAGS= -c -O3

LIBFLAGS2 = -L/usr/local/fftw/lib
LDFLAGS   = -I/usr/local/fftw/include


OBJECTS =            fftw.o global.o allocate_vars.o precmod.o stringmod.o write_pack.o interpolation_pack.o mesh_pack.o imod.o bc_setup.o statistics.o time_integrators.o time_integrators_MPI.o jacobians.o gmres_pack.o nonlinear_solvers.o time_loop.o
OBJECTS_MPI =        fftw.o global.o allocate_vars.o precmod.o stringmod.o write_pack.o interpolation_pack.o mesh_pack.o imod.o bc_setup.o statistics.o time_integrators.o time_integrators_MPI.o jacobians.o gmres_pack.o nonlinear_solvers.o time_loop_MPI.o
OBJECTS_FFT_THREAD = fftw.o global.o multithreading_benchmark.o
OBJECTS_FFT_RE =     fftw.o global.o re_to_comp_test.o
all : time_loop.exe time_loop_MPI.exe re_to_comp_test.exe multithreading_benchmark.exe
# all : re_to_comp_test.exe multithreading_benchmark.exe

time_loop.exe : $(OBJECTS)
	$(FC) -fopenmp -Wno-argument-mismatch $(LDFLAGS) -o time_loop.exe $(OBJECTS) $(LIBFLAGS1) -llapack -lblas $(LIBFLAGS2) -lfftw3 -lm

time_loop_MPI.exe : $(OBJECTS_MPI)
	$(FC) -fopenmp -Wno-argument-mismatch $(LDFLAGS) -o time_loop_MPI.exe $(OBJECTS_MPI) $(LIBFLAGS1) -llapack -lblas $(LIBFLAGS2) -lfftw3 -lm

re_to_comp_test.exe : $(OBJECTS_FFT_RE)
	$(FC) -fopenmp -Wno-argument-mismatch $(LDFLAGS) -o re_to_comp_test.exe $(OBJECTS_FFT_RE) $(LIBFLAGS1) -llapack -lblas $(LIBFLAGS2) -lfftw3 -lm

multithreading_benchmark.exe : $(OBJECTS_FFT_THREAD)
	$(FC) -Wno-argument-mismatch $(LDFLAGS) -o multithreading_benchmark.exe $(OBJECTS_FFT_THREAD) $(LIBFLAGS1) -llapack -lblas $(LIBFLAGS2) -lfftw3 -lm -lfftw3_threads --enable-threads

fftw.o : fftw.f90
	$(FC) $(FFLAGS) fftw.f90

global.o : global.f90
	$(FC) $(FFLAGS) global.f90

allocate_vars.o : allocate_vars.f90
	$(FC) $(FFLAGS) allocate_vars.f90

precmod.o : stringmod.f90
	$(FC) $(FFLAGS) precmod.f90

stringmod.o : stringmod.f90
	$(FC) $(FFLAGS) stringmod.f90

write_pack.o : write_pack.f90
	$(FC) $(FFLAGS) write_pack.f90

interpolation_pack.o : interpolation_pack.f90
	$(FC) $(FFLAGS) interpolation_pack.f90

mesh_pack.o : mesh_pack.f90
	$(FC) $(FFLAGS) mesh_pack.f90

imod.o : imod.f90
	$(FC) $(FFLAGS) imod.f90

bc_setup.o : bc_setup.f90
	$(FC) $(FFLAGS) bc_setup.f90

statistics.o : statistics.f90
	$(FC) $(FFLAGS) statistics.f90

time_integrators.o : time_integrators.f90
	$(FC) -fopenmp $(FFLAGS) time_integrators.f90

time_integrators_MPI.o : time_integrators_MPI.f90
	$(FC) -fopenmp -Wno-argument-mismatch  $(FFLAGS) time_integrators_MPI.f90

re_to_comp_test.o : re_to_comp_test.f90
	$(FC) -fopenmp $(FFLAGS) re_to_comp_test.f90

multithreading_benchmark.o : multithreading_benchmark.f90
	$(FC) -lfftw3_threads --enable-threads -lfftw3 -lm -Wno-argument-mismatch  $(FFLAGS) multithreading_benchmark.f90

jacobians.o : jacobians.f90
	$(FC) $(FFLAGS) jacobians.f90

gmres_pack.o : gmres_pack.f90
	$(FC) $(FFLAGS) gmres_pack.f90

nonlinear_solvers.o : nonlinear_solvers.f90
	$(FC) $(FFLAGS) nonlinear_solvers.f90

time_loop.o : time_loop.f90
	$(FC) $(FFLAGS) time_loop.f90

time_loop_MPI.o : time_loop_MPI.f90
	$(FC) $(FFLAGS) time_loop_MPI.f90

clean :
	rm -rf *.mod $(OBJECTS) $(PROGRAMS) 

cleanall :
	rm -rf *.mod *.txt $(OBJECTS) $(PROGRAMS) 
