ifeq ($(FCOMP),mpifort)
FCFLAGS = -O2 -freal-4-real-8 -march=native -flto -mavx
FCDEBUG = -g -fbacktrace -fcheck=all -fbounds-check -ffpe-trap=invalid,overflow,underflow,denormal
FCBUILD = -Wall -Wextra -pedantic
else ifeq ($(FCOMP),mpiifort)
FCFLAGS = -O2 -r8 -xHOST -ipo -assume byterecl
FCDEBUG = -g -debug -traceback -check all -CB -CU -CA -fpe0 -debug
FCBUILD = -warn all -stand f08
endif

PROGRAM =  mcgrid
SRCS =      ran2.f \
            constants.f90 \
            photon_vars.f90 \
            iarray.f90 \
            opt_prop.f90 \
            utils.f90 \
            subs.f90 \
            ch_opt.f90 \
            gridset.f90 \
            vector_class.f90 \
            inttau2.f90 \
            stokes.f90 \
            sourceph.f90 \
            writer.f90 \
            mcpolar.f90         
OBJECTS = $(SRCS:.f90=.o)

all:	$(PROGRAM)
debug:	FCFLAGS += $(FCDEBUG)
debug:	$(PROGRAM)
build:	FCFLAGS += $(FCBUILD)
build:	$(PROGRAM)

$(PROGRAM):	$(OBJECTS)
	$(FCOMP) $(FCFLAGS) -o $@ $^ 

%.o:  %.f90
	$(FCOMP)  $(FCFLAGS) -c $<

.PHONY:	clean

clean:
	rm -f *.o *.mod *.MOD mcgrid