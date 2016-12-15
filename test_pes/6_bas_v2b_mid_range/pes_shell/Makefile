FC = ifort -r8  -O
FFLAGS  = $(PES_MOD)
LIBS    = $(PES_LIB)  -mkl=sequential


SRC=\
pes_shell.f90 \
xxralpha.f90 \
update.f90 

OBJ=$(SRC:.f90=.o)

update.x : $(OBJ)
	$(FC) $(FFLAGS) $(OBJ) $(LIBS) -o update.x
	rm update.x update.o
	echo 'update sucessfull'

%.o : %.f90
	$(FC) -c $(FFLAGS) $(LIBS) $<
clean:
	rm -rf *.o *.mod
