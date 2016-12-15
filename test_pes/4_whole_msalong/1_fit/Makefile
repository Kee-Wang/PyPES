FC = ifort -r8 -O
FFLAGS  = $(PES_MOD)
LIBS    = $(PES_LIB) -mkl=sequential


SRC1=1_fit.f90
OBJ1=$(SRC1:.f90=.o)


1_fit.x : $(OBJ1)
	$(FC) $(FFLAGS) -o fit.x $(OBJ1) $(LIBS)
	rm -rf *.o *.mod
%.o : %.f90
	$(FC) -c $(FFLAGS) $(LIBS) $<
clean:
	rm -rf *.o *.mod
