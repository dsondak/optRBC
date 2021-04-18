FC=gfortran
FFLAGS= -c -O3
LIBFLAGS2 = -L/usr/local/Cellar/fftw/3.3.9/lib
LDFLAGS   = -I/usr/local/Cellar/fftw/3.3.9/include

#MAIN = Ra_loop
#MAIN = Ra_loop_no_opt
#MAIN = time_loop
MAIN = multithreading_benchmark

OBJECTS = fftw.o global.o allocate_vars.o precmod.o stringmod.o write_pack.o interpolation_pack.o mesh_pack.o imod.o bc_setup.o statistics.o time_integrators.o jacobians.o gmres_pack.o nonlinear_solvers.o $(MAIN).o
PROGRAMS = $(MAIN).exe

all: $(PROGRAMS)

$(PROGRAMS) : $(OBJECTS)
	$(FC) $(LDFLAGS) -o $(PROGRAMS) $(OBJECTS) $(LIBFLAGS1) -llapack -lblas $(LIBFLAGS2) -lfftw3 -lm -lfftw3_threads --enable-threads

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
	$(FC) $(FFLAGS) time_integrators.f90

jacobians.o : jacobians.f90
	$(FC) $(FFLAGS) jacobians.f90

gmres_pack.o : gmres_pack.f90
	$(FC) $(FFLAGS) gmres_pack.f90

nonlinear_solvers.o : nonlinear_solvers.f90
	$(FC) $(FFLAGS) nonlinear_solvers.f90

$(MAIN).o : $(MAIN).f90
	$(FC) $(FFLAGS) $(MAIN).f90

clean :
	rm -rf *.mod $(OBJECTS) $(PROGRAMS) 

cleanall :
	rm -rf *.mod *.txt $(OBJECTS) $(PROGRAMS) 
