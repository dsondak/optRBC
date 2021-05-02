FC=gfortran
FFLAGS= -c -O3

LIBFLAGS2 = -L/usr/local/fftw/lib
LDFLAGS   = -I/usr/local/fftw/include

MAIN = fft_c2c_loop
  
OBJECTS = fftw.o $(MAIN).o
PROGRAMS = $(MAIN).exe

all: $(PROGRAMS)

$(PROGRAMS) : $(OBJECTS)
	$(FC) $(LDFLAGS) -o $(PROGRAMS) $(OBJECTS) $(LIBFLAGS1) -llapack -lblas $(LIBFLAGS2) -lfftw3 -lm

fftw.o : fftw.f90
	$(FC) $(FFLAGS) fftw.f90


$(MAIN).o : $(MAIN).f90
	$(FC) $(FFLAGS) $(MAIN).f90

clean :
	rm -rf *.mod $(OBJECTS) $(PROGRAMS) 

cleanall :
	rm -rf *.mod *.txt $(OBJECTS) $(PROGRAMS) 
