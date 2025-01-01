# Makefile for Fortran programs

#FC = f90
#FC = g77
FC = f95
#FFLAGS = 
FFLAGS = -O3 -Wall -g
#FFLAGS = -Wall -g -pg
#FFLAGS = -C hopt -Wf"-L mrgmsg -pvctl fullmsg"
#FFLAGS = -g -p -C hopt -Wf"-L fmtlist"
#LDFLAGS =
#LIBS = 
#LIBS = -L. -lasl -lm
#LIBS = -lasl -lm

#OBJS = main.o cord.o cord_sph.o cord_time.o nlop.o mild_a_op.o b_op.o ch_fft4g.o fft4f2d.o setmix.o kainew.o spectrum.o melran.o

OBJS = main.o cord.o cord_time.o absorb.o refset.o setmix.o nlop.o prop.o fftsg.o fftsg2d.o spectrum.o plot2d_2.o sphere.o memory.o


default: a.out2

sb: $(OBJS)
	$(FC) $(FFLAGS) $(OBJS) $(LIBS) -o sb

a.out9: $(OBJS)
	$(FC) $(FFLAGS) $(OBJS) $(LIBS) -o a.out9

clean:
	\rm -f sb a.out *.o ftn* core *~ *.L

cleanall: distclean

distclean:
	\rm -f sb a.out ftn* core *.L fort.* mome.* wave.* ww.* smth.* memory-*.data *.o *~ \#*\#
