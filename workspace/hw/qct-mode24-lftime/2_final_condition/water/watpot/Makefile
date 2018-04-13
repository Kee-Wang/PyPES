#custom water Makefile
MODDIR := /raid/PES-repository/water_WHBB/include

#Please make sure all libs also compliled using OPenMP

FC := /raid/opt/intel2013/bin/ifort

#LIBS := -static $(PES_LIB_MP)  -mkl=parallel
LIBS := -static $(PES_LIB)  -mkl=sequential


PATH1 :=  /home/kee/pes-lib/pes-hwww/src/src
FFLAGS := -u -O -r8 $(PES_MOD)  -I $(PATH1)

.PHONY : clean

OBJ := \
$(PATH1)/constants.o\
$(PATH1)/math.o\
$(PATH1)/smear.o\
$(PATH1)/ttm3f_mod2.o\
$(PATH1)/mnasa_mod.o\
$(PATH1)/mnasa.o\
$(PATH1)/ttm3f_mb.o\
$(PATH1)/pot_monomer_mod.o\
$(PATH1)/pot_monomer.o\
$(PATH1)/h2o_3bdy.o\
$(PATH1)/pes3b.o\
$(PATH1)/water_shell.o\
$(PATH1)/util.o\
$(PATH1)/hcl3_3bdy.o\
$(PATH1)/hcl_monomer.o\
$(PATH1)/hcl2_two.o\
$(PATH1)/hcl_shell.o\
$(PATH1)/w1h1_2bdy.o\
$(PATH1)/w2h1_3bdy.o\
$(PATH1)/w1h2_3bdy.o\
$(PATH1)/hw_shell.o\
$(PATH1)/pes_shell.o\
$(PATH1)/inv_share.o\
$(PATH1)/inv_mg42_j.o\
$(PATH1)/getr0.o\
$(PATH1)/getrvec.o\
$(PATH1)/getvec.o\
$(PATH1)/getx0.o\
$(PATH1)/h4o2.pes4.o\
#$(PATH1)/getpot.o

#main.o
SRC := eckart_math.o angular.f90 final_condition.o wat_pot.o
OBJ1 := $(SRC:.f90=.o)

%.o : %.f90
	$(FC) -c  $(FFLAGS) $(LIBS) $<
%.o : %.f
	$(FC) -c  $(FFLAGS) $(LIBS) $<

all : xget.x

xget.x : $(OBJ) $(OBJ1)
	$(FC) -o $@ $^  $(FFLAGS) $(LIBS)

water_shell.o: pot_monomer_mod.o nbody.f90
pot.o:      water_shell.o
pot_monomer.o: pot_monomer_mod.o
ttm3f_mb.o:    ttm3f_mod2.o
mnasa.o:       ttm3f_mod2.o mnasa_mod.o

clean: 
	rm -rf *.o *.x *.mod

