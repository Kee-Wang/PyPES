#=====================Info=================


#=====================User define area===================
input1=getpot.f90
output1=getpot.x
#========================================================

FC = ifort -r8
path0 = /home/kee/pes-lib/pes-co2h2o
path1 = /home/kee/pes-lib/pes-co2h2o-sr
path2 = /home/kee/pes-lib/pes-co2h2o-lr
path3 = /home/kee/pes-lib/pes-co2h2o-2b
path4 = /home/kee/pes-lib/pes-co2-huang
path5 = /home/kee/pes-lib/pes-h2o-ps
path6 = /home/kee/pes-lib/vnbpes
path7 = /home/kee/pes-lib/mods

FFLAGS =  $(PES_MOD) -I $(path0) -I $(path1) -I $(path2) -I $(path3) -I $(path4) -I $(path5) -I $(path6) -I $(path7)
LIB =  $(PES_LIB)

SRC1 := \
$(path0)/pes_co2h2o.o \
$(path1)/pes_co2h2o_sr.o\
$(path2)/basis.o \
$(path2)/pes_co2h2o_lr.o \
$(path3)/pes_co2h2o_2b.o \
$(path4)/pes_co2_huang.o \
$(path5)/pes_h2o_ps_mod.o \
$(path5)/pes_h2o_ps.o \
$(path6)/pes_tools.o\
$(path7)/constants.o


OBJECT1 := $(SRC1) $(input1:.f90=.o)


$(output1) : $(OBJECT1)
	$(FC) -o $@ $^ $(LIB) $(FFLAGS)
#	rm -f *.o *.mod 
$(output2) : $(OBJECT2)
	$(FC) -o $@ $^ $(LIB) $(FFLAGS)
#	rm -f *.o *.mod 

#get.x : $(OBJ) $(OBJGET)
#	$(FC) -o $@ $^ $(FFLAGS) $(LIB)

#xopt.x : $(OBJ) $(OBJOPT)
#	$(FC) -o $@ $^ $(FFLAGS) $(LIB)

%.o : %.f90
	$(FC) -c $(FFLAGS) $(LIB) $<

clean:
	rm -f *.x *.o *.mod

.PHONY: all clean
