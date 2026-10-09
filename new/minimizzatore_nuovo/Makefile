PROG =	min_multi

SRCS =	module_parameters.f90 module_min_mod.f90\
	drive_min.f90 read_min.f90 read_gupta.f90\
        xyz_writer.f90 read_initial.f90\
	minimization_gupta.f90 minimization_gupta_v.f90\
	l-bfgs-b.f90 bigvoi.f90 bigvoi_md.f90 voisin.f90 force_gupta.f90\
	minimization_in_loop.f90
	
	



OBJS =	module_parameters.o module_min_mod.o\
	drive_min.o read_min.o read_gupta.o\
        xyz_writer.o read_initial.o\
	minimization_gupta.o minimization_gupta_v.o\
	l-bfgs-b.o bigvoi.o bigvoi_md.o voisin.o force_gupta.o\
	minimization_in_loop.o

	


#LIBS = 

#FC=pgf90
#FC=gfortran
#FFLAGS=-static-libcxa -Xlinker -Bstatic -lsvml\
#-L/opt/intel/mkl/lib/32/ -lmkl_lapack -lmkl_p4 \
#-lguide -Vaxlib -Xlinker -Bdynamic -lpthread
#FC=ifort
#FC= pgf90
FC=gfortran -O3
FFLAGS= -std=legacy
F90=gfortran -O3
#F90=pgf90
#F90 =pgf90
#F90= ifort  -O1
F90FLAGS= -std=legacy
#LDFLAGS = -Wl,-R/usr/lib/gcc-lib/i486-linux/3.3.5 

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(LDFLAGS) $(F90FLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) *.mod

.SUFFIXES: $(SUFFIXES) .f90

.f90.o:
	$(F90) $(F90FLAGS) -c  $<

