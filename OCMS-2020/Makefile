#
# Makefile for the OCMS program package
#

# Include the file defining macro-variables for build.
# (Ex) make all params=$workspace/build-params.mk
include ${params}
include ${extras}

# Define SHELL variable.
# It avoids potential troubles in case of inherited
# SHELL variable from the environment.
SHELL = /bin/sh

# Define set of suffixes we use hereafter.
# The first command resets the current list and the second declares it.
.SUFFIX:
.SUFFIX: .F .F90 .o

# Define VPATH to allow automatic search of the prerequisites. If automatic
# variables are used, the matching will be added to the prerequisite name.
VPATH = src:src/def_cavities
FC = gfortran
CFLAGS = -DNBE=$(NBE) -DSYM=$(SYM) -DSIZE_PARAM=$(SIZE_PARAM)
ifdef DISPERSION
CFLAGS = -DNBE=$(NBE) -DSYM=$(SYM) -DSIZE_PARAM=$(SIZE_PARAM) -D$(DISPERSION)
endif
LFLAGS = -lslatec -llapack -lblas

DET = det.$(PLZ).NBE=$(NBE)
DET_DIRICHLET = det_dirichlet.$(PLZ).NBE=$(NBE)
WFUNC = wfunc.$(PLZ).NBE=$(NBE)
HUSIMI = husimi.$(PLZ).NBE=$(NBE)
CAVITYDATA = cavity_data
ESTIMATENBE = estimateNBE
ESTIMATERATIO = estimateRatio

# Compile the object pointed to by the CAVITY variable only once.
# It is to avoid that the object is compiled many times.
CAVITYOBJ = $(CAVITY:.F=.o)
OBJS = $(CAVITYOBJ) def_hankel_zbesh.o def_green_func.o def_matrix.o
OBJS4DET = $(OBJS) main.det.o
OBJS4DET_DIRICHLET = $(OBJS) main.det_dirichlet.o
OBJS4WF = $(OBJS) main.wfunc.o
OBJS4HUSIMI = $(OBJS) main.husimi.o
OBJS4CAVITY = $(CAVITYOBJ) main.cavity.o
OBJS4ESTIMATENBE = $(CAVITYOBJ) main.estimateNBE.o
OBJS4ESTIMATERATIO = $(CAVITYOBJ) main.estimateRatio.o

# A phony target is not related to any file; This avoids conflicts (e.g.
# there is a file named clean.F or all.F in the src directory) and it
# improves performances.
.PHONY: all clean cleanall distclean

default: det

all: det det_dirichlet wfunc husimi cavity estimateNBE estimateRatio clean

# Automatic variables
# $^ = return the list of all prerequisites
# $< = return the first prerequisite
# $@ = return the target
det: $(OBJS4DET)
	$(FC) $(CFLAGS) $^ -o $(DET) $(LFLAGS)

det_dirichlet: $(OBJS4DET_DIRICHLET)
	$(FC) $(CFLAGS) $^ -o $(DET_DIRICHLET) $(LFLAGS)

wfunc:  $(OBJS4WF)
	$(FC) $(CFLAGS) $^ -o $(WFUNC) $(LFLAGS)

husimi: $(OBJS4HUSIMI)
	$(FC) $(CFLAGS) $^ -o $(HUSIMI) $(LFLAGS)

cavity: $(OBJS4CAVITY)
	$(FC) $(CFLAGS) $^ -o $(CAVITYDATA) $(LFLAGS)

estimateNBE: $(OBJS4ESTIMATENBE)
	$(FC) $(CFLAGS) $^ -o $(ESTIMATENBE)

estimateRatio: $(OBJS4ESTIMATERATIO)
	$(FC) $(CFLAGS) $^ -o $(ESTIMATERATIO)

main.det.o: main.det.F
	$(FC) $(CFLAGS) -D$(PLZ) -c $^ -o $@

main.det_dirichlet.o: main.det_dirichlet.F
	$(FC) $(CFLAGS) -D$(PLZ) -c $^ -o $@

main.wfunc.o: main.wfunc.F
	$(FC) $(CFLAGS) -D$(PLZ) -c $^ -o $@

main.husimi.o: main.husimi.F
	$(FC) $(CFLAGS) -D$(PLZ) -DIXMAX=$(IXMAX) -DIYMAX=$(IYMAX) -c $^ -o $@

main.estimateNBE.o: main.estimateNBE.F90
	$(FC) $(CFLAGS) -c $^ -o $@

main.estimateRatio.o: main.estimateRatio.F90
	$(FC) $(CFLAGS) -c $^ -o $@

%.o: %.F
	$(FC) $(CFLAGS) -c $^ -o $@

clean:
	\rm -f *.o *~

cleanall: distclean

distclean:
	\rm -f det wfunc husimi estimateNBE estimateRatio cavity_data det*.T?.NBE=* wfunc.T?.NBE=* husimi.T?.NBE=* *.o *~
