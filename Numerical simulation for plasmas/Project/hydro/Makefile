#################################################################### OVERVIEW
#  Makefile for 
#
#  HYDRO: SAMPLE PARALLEL CODE FOR HYDRODYNAMIC SIMULATION
#  VERSION 1.0
#
#  Copyright, 2009
#  Gregory Howes
#
#  Version notes:
#  1. Code has been ported to BigBen, Cray XT3 at Pittsburgh Supercomputing Center
#	                      MacPro
#                             Moffett, SiCortex at Purdue
#			      Helium, at Univ of Iowa
#  LAST UPDATE:  25 MAY 2011
###############################################################################
# SET VARIABLES FOR SYSTEM, DEBUGGING, and PROFILING
# SYSTEM options: helium, SiCortex, XT3, MacPro, uiowa,
# SET export SYSTEM=argon
SYSTEM=ArchLinux
DEBUG= false
PROFILE= false

ifeq ($(SYSTEM), ArchLinux)
	COMP = mpif90
	FLAGS = -fdefault-real-8 -O3 -fopenmp
	ifeq ($(PROFILE), true)
		FLAGS += -g
	endif
endif

###############################################################################
###############################################################################
#SET FLAGS FOR U IOWA Shared Research Cluster neon.hpc.uiowa.edu
ifeq ($(SYSTEM),argon)
	COMP= mpif90
	FLAGS = -r8 -double-size 128 -O3 -fopenmp
	ifeq ($(PROFILE),true)
		FLAGS += -g
	endif
endif

#SET FLAGS FOR U IOWA Shared Research Cluster neon.hpc.uiowa.edu
ifeq ($(SYSTEM),neon)
	COMP= mpif90
	FLAGS = -r8 -double-size 128 -O3
	ifeq ($(PROFILE),true)
		FLAGS += -g
	endif
endif

#SET FLAGS FOR U IOWA Shared Research Cluster helium.hpc.uiowa.edu
ifeq ($(SYSTEM),helium)
	COMP= mpif90
	FLAGS = -r8 -double-size 128 -O3
	ifeq ($(PROFILE),true)
		FLAGS += -g
	endif
endif

#SET FLAGS FOR Moffett, SiCortex
ifeq ($(SYSTEM),SiCortex)
	COMP= mpif90
	FLAGS = -r8 -O3
#	ifeq ($(DEBUG),true)
	#For TotalView Debugger
#		FLAGS = -r8 -g -O0
#	endif		
	ifeq ($(PROFILE),true)
		FLAGS += -g
	endif
endif

#SET FLAGS FOR U IOWA 64-bit Research Cluster rs-003.its.uiowa.edu
ifeq ($(SYSTEM),uiowa)
	COMP= mpif90
	FLAGS = -fdefault-real-8 -O3
	ifeq ($(PROFILE),true)
		FLAGS += -g
	endif
endif

#SET FLAGS FOR Bigben, Cray XT3
ifeq ($(SYSTEM),XT3)
	COMP= ftn
	FLAGS = -fastsse -r8
	ifeq ($(DEBUG),true)
	#For TotalView Debugger
		FLAGS = -r8 -g -O0
	endif		
	ifeq ($(PROFILE),true)
	#For CrayPAT profiling
		FLAGS += -Mprof=func
	endif
endif

#SET FLAGS FOR MacPro
ifeq ($(SYSTEM),MacPro)
	COMP= mpif90
	FLAGS = -O3 -r8 -ipo -vec_report0
	ifeq ($(DEBUG),true)
		FLAGS = -g -implicitnone -warn all -check bounds -traceback -r8 	
	endif		
endif

LIBS=	

HYDX=  hydro_mpi_var.o hydro_var.o hydro_param.o hydro_hd.o hydro_grid.o \
	hydro_ic.o hydro_diag.o hydro_bcs.o hydro_dt.o hydro.o	

VPATH= src:include:/usr/include/

###############################################################################
all: hydro neat

hydro:  $(HYDX)
	$(COMP) -o hydro.e $(FLAGS) $(LIBS) $(HYDX)


###############################################################################

###############################################################################

neat:	
	mv *.o include/
	mv *.mod include/

clean:
	rm -f *.o include/*.o
	rm -f *.mod include/*.mod
	rm -f hydro.e

tar: 
	tar -cvf hyd`date +'%y%m%d'`.tar Makefile src/*.f90 include/ \
		data/ sample*in sample1_argon.sh hydro_documentation.pdf

#########Dependencies
%.o : %.f90
	$(COMP) -c $(FLAGS) $<


hydro_mpi_var.o:
hydro_var.o:
hydro_param.o: hydro_mpi_var.o
hydro_hd.o: hydro_var.o
hydro_grid.o: hydro_mpi_var.o hydro_param.o hydro_var.o
hydro_ic.o: hydro_mpi_var.o hydro_param.o hydro_var.o hydro_grid.o
hydro_diag.o: hydro_mpi_var.o hydro_param.o hydro_var.o hydro_grid.o
hydro_bcs.o: hydro_mpi_var.o hydro_var.o hydro_grid.o
hydro_dt.o: hydro_mpi_var.o hydro_var.o hydro_param.o hydro_grid.o \
	hydro_bcs.o hydro_hd.o hydro_ic.o 
hydro.o:   hydro_mpi_var.o hydro_param.o hydro_var.o  hydro_hd.o hydro_grid.o \
	hydro_ic.o hydro_diag.o hydro_bcs.o hydro_dt.o
