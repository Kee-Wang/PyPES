#compiles pes_fit and getpot program

iinput=getpot.f90
output=a.out

#----------------------------------
FC = ifort -r8  -O #-check
LIBS    = $(PES_LIB)  -mkl=sequential
PAATH=  ../pes_shell/
FFLAGS  = $(PES_MOD) -I $(PAATH)
SRC=  \
$(PAATH)pes_shell.o \
$(PAATH)xxralpha.o


%.o : %.f90
	$(FC) -c $(FFLAGS) $(LIBS) $<

OBJ=$(iinput:.f90=.o) $(SRC) 


$(output) : $(OBJ)
	$(FC) $(FFLAGS) $(OBJ) $(LIBS) -o $(output)
	rm -rf *.o *.mod

clean:
	rm -rf *.o *.mod
