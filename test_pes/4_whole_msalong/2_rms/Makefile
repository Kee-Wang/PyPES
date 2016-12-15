
#=====================User define area===================
input1=rms.f90
output1=rms.x
#========================================================

FC = ifort -r8 -O
LIBS = $(PES_LIB) -mkl=sequential
PAATH= ../pes_shell/
MMOD= $(PES_MOD)  -I $(PAATH)

SRC1 =\
$(PAATH)xxralpha.o\
$(PAATH)pes_shell.o


OBJECT1 := $(SRC1) $(input1:.f90=.o) 


$(output1) : $(OBJECT1)
	$(FC) -o $@ $^ $(LIBS)
	rm -f *.o *.mod 

%.o: %.f90
	$(FC) -c $< $(MMOD)
clean:
	rm -f *.o *.mod *.x	
