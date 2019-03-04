#################################################################### OVERVIEW
#  Makefile for HW #6b Code: Euler and Leapfrog Timestepping
#  PHYS:5905 Numerical Simulation of Plasmas
#  Prof. Gregory Howes
#
#  Version notes:
#
###############################################################################
# Intel Fortran compiler: ifort
FLAGS=  -O4 -r8 -double-size 128
# DEBUG
#FLAGS=  -g -r8 -double-size 128

# gfortran compiler: gfortran
#FLAGS= -O3 -DDOUBLE -ffree-form -ffree-line-length-512
#DEBUG 
#FLAGS= -g -DDOUBLE -ffree-form -ffree-line-length-512

COMP=   ifort
#COMP= gfortran

LIBS=	

HW6MOD= hw6_params.o hw6_funcs.o

HW6X=  	hw6.o

VPATH= 	src:include

###############################################################################
all: clean hw6 neat

hw6: $(HW6MOD) $(HW6X)
	$(COMP) -o hw6.e $(FLAGS) $(LIBS) $(HW6MOD) $(HW6X)

###############################################################################

###############################################################################

neat:	
	mv *.mod include/
	mv *.o include/

clean:
	rm -f include/*.o include/*.mod
	rm -f hw6*.e
	rm -f *.o

tar: 
	tar -cvf  hw6_`date +'%y%m%d'`.tar Makefile src/*.f90 include/* run_hw6.sh


#########Rules
%.o : %.f90
	$(COMP) -c $(FLAGS) $<

