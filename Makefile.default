#	  /$$$$$$                                /$$ /$$                    
#	 /$$__  $$                              | $$|__/                    
#	| $$  \__/  /$$$$$$   /$$$$$$$  /$$$$$$ | $$ /$$ /$$$$$$$   /$$$$$$ 
#	| $$ /$$$$ |____  $$ /$$_____/ /$$__  $$| $$| $$| $$__  $$ /$$__  $$
#	| $$|_  $$  /$$$$$$$|  $$$$$$ | $$  \ $$| $$| $$| $$  \ $$| $$$$$$$$
#	| $$  \ $$ /$$__  $$ \____  $$| $$  | $$| $$| $$| $$  | $$| $$_____/
#	|  $$$$$$/|  $$$$$$$ /$$$$$$$/|  $$$$$$/| $$| $$| $$  | $$|  $$$$$$$
#	 \______/  \_______/|_______/  \______/ |__/|__/|__/  |__/ \_______/
                                                                                     
# Makefile for Gasoline (specify below)
# Then type make <platform>
# <platform> should be one of:
# null, pvm, pthread, ampi, lam_mpi, qmpi, mpi,
# charm
#

# Which C compiler should I use?
CC = gcc

#	 /$$      /$$                 /$$          
#	| $$$    /$$$                | $$          
#	| $$$$  /$$$$  /$$$$$$   /$$$$$$$  /$$$$$$ 
#	| $$ $$/$$ $$ /$$__  $$ /$$__  $$ /$$__  $$
#	| $$  $$$| $$| $$  \ $$| $$  | $$| $$$$$$$$
#	| $$\  $ | $$| $$  | $$| $$  | $$| $$_____/
#	| $$ \/  | $$|  $$$$$$/|  $$$$$$$|  $$$$$$$
#	|__/     |__/ \______/  \_______/ \_______/
	                                           
# What kind of problem are we running?  This will enable a default set of code components
# appropriate for different kinds of problems.  Select only one.
MODE_DEF = -DGALAXY            #Default mode: includes SPH, star formation, and feedback.
#MODE_DEF = -DCOLLISIONAL		#Collisional mode for simulating asteroids & planets
#MODE_DEF = -DNBODY 			#N-body only mode: gravity with nothing else.
#MODE_DEF = -DJUSTHYDRO			#SPH only mode: just do hydro/gravity and nothing else
#MODE_DEF = -DGLASS_IC 			#Generate a gravity glass (useful for making ICS)

#	 /$$   /$$                                         /$$
#	| $$  /$$/                                        | $$
#	| $$ /$$/   /$$$$$$   /$$$$$$  /$$$$$$$   /$$$$$$ | $$
#	| $$$$$/   /$$__  $$ /$$__  $$| $$__  $$ /$$__  $$| $$
#	| $$  $$  | $$$$$$$$| $$  \__/| $$  \ $$| $$$$$$$$| $$
#	| $$\  $$ | $$_____/| $$      | $$  | $$| $$_____/| $$
#	| $$ \  $$|  $$$$$$$| $$      | $$  | $$|  $$$$$$$| $$
#	|__/  \__/ \_______/|__/      |__/  |__/ \_______/|__/

# Which SPH Kernel should we use?
KERNEL_DEF = -DWENDLANDC4            #Use the Wendland C_4 Kernel (See Dehnen & Aly 2012)
#KERNEL_DEF = -DM4                    #The "classic" M4 Kernel
#KERNEL_DEF = -DWENDLAND              #Use the Wendland C_2 Kernel (See Dehnen & Aly 2012)
#KERNEL_DEF = -DHSHRINK               #M4 kernel times (pi/6)^(1/3) for a tighter kernel
#KERNEL_DEF = -DM43D                  #Use a fancier 3D derived M4 Kernel
#KERNEL_DEF = -DPEAKEDKERNEL          #Use a modified Peaked M4 kernel as per Thomas and Couchman 92 (Bibcode:1992MNRAS.257...11T)
#KERNEL_DEF = -DQUINTIC               #Use a quintic spline kernel from Dehnen & Aly 2012 (Bibcode:2012MNRAS.425.1068D)


#	  /$$$$$$                      /$$ /$$                    
#	 /$$__  $$                    | $$|__/                    
#	| $$  \__/  /$$$$$$   /$$$$$$ | $$ /$$ /$$$$$$$   /$$$$$$ 
#	| $$       /$$__  $$ /$$__  $$| $$| $$| $$__  $$ /$$__  $$
#	| $$      | $$  \ $$| $$  \ $$| $$| $$| $$  \ $$| $$  \ $$
#	| $$    $$| $$  | $$| $$  | $$| $$| $$| $$  | $$| $$  | $$
#	|  $$$$$$/|  $$$$$$/|  $$$$$$/| $$| $$| $$  | $$|  $$$$$$$
#	 \______/  \______/  \______/ |__/|__/|__/  |__/ \____  $$
#	                                                 /$$  \ $$
#	                                                |  $$$$$$/
#	                                                 \______/ 
#
# The following pairs of lines define which cooling method to use.  Uncomment
# only one pair

COOLING_OBJ = cooling_metal.o
COOLING_DEF = -DCOOLING_METAL

#COOLING_OBJ = cooling_cosmo.o
#COOLING_DEF = -DCOOLING_COSMO

#COOLING_OBJ = cooling_poly.o
#COOLING_DEF = -DCOOLING_POLY

#COOLING_OBJ = cooling_grackle.o
#COOLING_DEF = -DCOOLING_GRACKLE
#COOLING_LIB = libgrackle.so -lhdf5

#COOLING_OBJ = cooling_metal_H2.o
#COOLING_DEF = -DCOOLING_MOLECULARH

#COOLING_OBJ = cooling_planet.o
#COOLING_DEF = -DCOOLING_PLANET

#	 /$$$$$$$$             /$$                       
#	| $$_____/            | $$                       
#	| $$       /$$   /$$ /$$$$$$    /$$$$$$  /$$$$$$ 
#	| $$$$$   |  $$ /$$/|_  $$_/   /$$__  $$|____  $$
#	| $$__/    \  $$$$/   | $$    | $$  \__/ /$$$$$$$
#	| $$        >$$  $$   | $$ /$$| $$      /$$__  $$
#	| $$$$$$$$ /$$/\  $$  |  $$$$/| $$     |  $$$$$$$
#	|________/|__/  \__/   \___/  |__/      \_______/
#	                                                 
#	 /$$$$$$$$                    /$$                                            
#	| $$_____/                   | $$                                            
#	| $$     /$$$$$$   /$$$$$$  /$$$$$$   /$$   /$$  /$$$$$$   /$$$$$$   /$$$$$$$
#	| $$$$$ /$$__  $$ |____  $$|_  $$_/  | $$  | $$ /$$__  $$ /$$__  $$ /$$_____/
#	| $$__/| $$$$$$$$  /$$$$$$$  | $$    | $$  | $$| $$  \__/| $$$$$$$$|  $$$$$$ 
#	| $$   | $$_____/ /$$__  $$  | $$ /$$| $$  | $$| $$      | $$_____/ \____  $$
#	| $$   |  $$$$$$$|  $$$$$$$  |  $$$$/|  $$$$$$/| $$      |  $$$$$$$ /$$$$$$$/
#	|__/    \_______/ \_______/   \___/   \______/ |__/       \_______/|_______/ 
                                                                                                                              
EXTRA_DEF =                                                                  
## GRAVITY AND HYDRODYNAMICS (REQUIRE -DGASOLINE)
#EXTRA_DEF += -DCULLENDEHNEN          #Use the time-variable Cullen & Dehnen 2010 viscosity limit switch
#EXTRA_DEF += -DVARALPHA              #Use a Morris & Monaghan-style time variable viscosity
#EXTRA_DEF += -DINFLOWOUTFLOW         #Enable inflow/outflow boundary conditions
#EXTRA_DEF += -DOUTURBDRIVER          #Drive turbulence as in Price & Federrath 2010 with OU variables
#EXTRA_DEF += -DROT_FRAME             #Use a rotating reference frame with rotation rate given by dOmega
#EXTRA_DEF += -DSIMPLE_GAS_DRAG       #Apply a drag to the gas velocities using either the Epstein regime, or stopping time
#EXTRA_DEF += -DPEXT                  #Allow non-zero external pressure 
#EXTRA_DEF += -DJEANSSMOOTH 		  #Calculate presure floor to keep jeans length larger than SPH h only (default is larger of SPH h or eps)
#EXTRA_DEF += -DJEANSSOFTONLY         #Calculate presure floor to keep jeans length larger than eps only
#EXTRA_DEF += -DSINKINGAVERAGE        #Use a smoothed relative velocity rather than per-particle for SMBH accretion rates

## STAR FORMATION AND FEEDBACK (REQUIRE -DSTARFORM)
#EXTRA_DEF += -DBLASTWAVE			  #Use the Stinson+ 2006 Blastwave feedback 
#EXTRA_DEF += -DNONTHERMAL			  #Use an Agertz+ 2013/Teyssier+ 2013 non-thermal pressure feedback
#EXTRA_DEF += -DKROUPA01              #Use the newer Kroupa 01 IMF (DOI:10.1046/j.1365-8711.2001.04022.x)
#EXTRA_DEF += -DKROUPA                #Use the Kroupa IMF (from Kroupa et al 1993 Bibcode:1993MNRAS.262..545K)
#EXTRA_DEF += -DMILLERSCALO			  #Use the old-school Miller & Scalo 1979 IMF (DOI:10.1086/190629)
#EXTRA_DEF += -DFEEDBACKDIFFLIMIT     #Disable thermal diffusion for particles that have their cooling shut off (requires -DBLASTWAVE)
#EXTRA_DEF += -DTOPHATFEEDBACK        #Use a tophat kernel to smooth feedback over (default with superbubble feedback)
#EXTRA_DEF += -DVOLUMEFEEDBACK        #Smooth FB energy with a volume-weighted kernel
#EXTRA_DEF += -DMASSFEEDBACK          #Smooth FB energy over mass-weighted kernel (default with blastwave feedback)

## COLLISIONS AND PLANET FORMATION (REQUIRE -DCOLLISIONS)
#EXTRA_DEF += -DAGGS                  #Include support for aggregates (like asteroids and rubble piles) See aggs.c
#EXTRA_DEF += -DGROWMASS              #Allow collisional particles to grow
#EXTRA_DEF += -DN_DIM                 #Number of dimensions to use for solar system calculations (Defaults to 3)
#EXTRA_DEF += -DSURFACEAREA           #Calculate and output the surface area of all the particles.
#EXTRA_DEF += -DNORMAL                #Calculate and output surface normals for particles (Needs -DSURFACEAREA)
#EXTRA_DEF += -DRUBBLE_ZML            #Use "Rubble Pile" collision model 
#EXTRA_DEF += -DSAND_PILE             #Use the "Sand Pile" collision model 
#EXTRA_DEF += -DSLIDING_PATCH         #Alternative collision scheme (I don't really know what this does) 
#EXTRA_DEF += -DTUMBLER               #Generate a hard-walled cylinder boundary (for use with collisions and sand piles)

## OUTPUT VISUALIZATION
#EXTRA_DEF += -DGSS_DUMPFRAME         #Dumpframe that allows coloring by particle property (Greg's version)
#EXTRA_DEF += -DVOXEL                 #Dump voxels for volume imaging.

## INITIAL CONDITION BUILDING
#EXTRA_DEF += -DGLASSZ                #Build a glass only in the z-direction, useful for building disc initial conditions

## DEBUG
#EXTRA_DEF += -DSETTRAPFPE            #Enable Floating point exceptions

#  /$$       /$$ /$$                                    /$$                    
#  | $$      |__/| $$                                   |__/                    
#  | $$       /$$| $$$$$$$   /$$$$$$  /$$$$$$   /$$$$$$  /$$  /$$$$$$   /$$$$$$$
#  | $$      | $$| $$__  $$ /$$__  $$|____  $$ /$$__  $$| $$ /$$__  $$ /$$_____/
#  | $$      | $$| $$  \ $$| $$  \__/ /$$$$$$$| $$  \__/| $$| $$$$$$$$|  $$$$$$ 
#  | $$      | $$| $$  | $$| $$      /$$__  $$| $$      | $$| $$_____/ \____  $$
#  | $$$$$$$$| $$| $$$$$$$/| $$     |  $$$$$$$| $$      | $$|  $$$$$$$ /$$$$$$$/
#  |________/|__/|_______/ |__/      \_______/|__/      |__/ \_______/|_______/ 
#                                                                               

# If you do dump frame stuff with png you need these libs ( -DUSE_PNG )
#PNG_INCL = -I/usr/include
#PNG_LIB = -L/usr/lib -lpng -lz
#PNG_LIB = -lpng -lz
#PNG_OBJ = writepng.o
#PNG_DEF = $(PNG_INCL) -DUSE_PNG
PNG_INCL =
PNG_LIB =
PNG_OBJ =
PNG_DEF =

# if you use OUTURBDRIVER or TURBULENT you need libgsl
#TURB_OBJ = outurb.o
#GSL_LIB = -lgsl -lgslcblas
TURB_OBJ =
GSL_LIB =

BASE_LD_FLAGS = $(PNG_LIB) $(GSL_LIB)
BASE_LD_FLAGS = $(PNG_LIB) $(GSL_LIB) $(COOLING_LIB)

# Remove the object for NBODY runs
ifeq ($(strip $(MODE_DEF)),-DNBODY)
COOLING_OBJ =
endif

BASE_DEF = $(PNG_DEF) $(COOLING_DEF) $(MODE_DEF) $(KERNEL_DEF) $(EXTRA_DEF)

# The filename of the compiled executable
EXE = gasoline

#
#       NULL defines
#
NULL_MDL		= ../mdl/null
NULL_CFLAGS		= -O3 -I$(NULL_MDL) $(BASE_DEF)
NULL_LD_FLAGS	= $(BASE_LD_FLAGS) #-L/usr/lib -L/lib
NULL_XOBJ		= erf.o v_sqrt1.o
NULL_LIBMDL		= $(NULL_MDL)/mdl.o -lm

#
#       DBG defines
#
DBG_MDL		= ../mdl/null
DBG_CFLAGS		= -g -I$(DBG_MDL) $(BASE_DEF)
DBG_LD_FLAGS	= $(BASE_LD_FLAGS) -L/usr/lib -L/lib
DBG_XOBJ		= erf.o v_sqrt1.o
DBG_LIBMDL		= $(NULL_MDL)/mdl.o -lm

#
#       LINUX AMPI (Charm MPI) defines
#
CHARM=../charm/net-linux/bin/charmc
CHARMLINK=
AMPI_MDL			= ../mdl/ampi
AMPI_CFLAGS		= -O3 -malign-double -mpentiumpro -I$(AMPI_MDL) $(BASE_DEF)
AMPI_LD_FLAGS		=  $(BASE_LD_FLAGS)
AMPI_XOBJ                = erf.o v_sqrt1.o
AMPI_LIBMDL              = $(AMPI_MDL)/mdl.o -language ampi $(CHARMLINK) -lm

#
#       LINUX LAM MPI defines
#
LAM_MDL			= ../mdl/mpi
LAM_CFLAGS		= -fast -I$(LAM_MDL) $(BASE_DEF) -DMPI_LINUX
LAM_LD_FLAGS		=  $(BASE_LD_FLAGS)
LAM_XOBJ                = erf.o v_sqrt1.o
LAM_LIBMDL              = $(LAM_MDL)/mdl.o -lm
LAM_MDL_CFLAGS  = -fast -I$(LAM_MDL) $(BASE_DEF)

#
#       Quadrics MPI defines
#
QMPI_MDL                        = ../mdl/mpi
QMPI_CFLAGS             = -O3 -I$(QMPI_MDL) $(BASE_DEF) -DMPI_LINUX -DCCC
QMPI_LD_FLAGS           =  $(BASE_LD_FLAGS)
QMPI_XOBJ                = erf.o
QMPI_LIBMDL              = $(QMPI_MDL)/mdl.o -lmpi -lm
QMPI_MDL_CFLAGS  = -I$(QMPI_MDL) $(BASE_DEF) -DMPI_LINUX


#       Basic MPI defines
#
MPI_MDL			= ../mdl/mpi
MPI_CFLAGS		= -O3 -I$(MPI_MDL) $(BASE_DEF)
MPI_LD_FLAGS	= $(BASE_LD_FLAGS)
MPI_XOBJ		= v_sqrt1.o
MPI_LIBMDL		= $(MPI_MDL)/mdl.o -lm
MPI_MDL_CFLAGS	= -O3

#
#		PVM defines
#
PVMDIR	= $(PVM_ROOT)
PVMLIB	= $(PVMDIR)/lib/$(PVM_ARCH)/libpvm3.a
BDIR	= $(HOME)/pvm3/bin
XDIR	= $(BDIR)/$(PVM_ARCH)
PVM_MDL		= ../mdl/pvm
PVM_CFLAGS	= -O3 -I$(PVMDIR)/include -I$(PVM_MDL) $(BASE_DEF)
PVM_XOBJ	= v_sqrt1.o
PVM_LIBMDL	= $(PVM_MDL)/$(PVM_ARCH)/mdl.o $(PVMLIB) $(ARCHLIB) -lm
PVM_LD_FLAGS	= $(BASE_LD_FLAGS)

#
#       PTHREAD defines
#
PTHREAD_MDL			= ../mdl/pthread
PTHREAD_CFLAGS		= -O3  -D_REENTRANT -I$(PTHREAD_MDL) $(BASE_DEF)
PTHREAD_LD_FLAGS 	=  $(BASE_LD_FLAGS)
PTHREAD_XOBJ		= erf.o v_sqrt1.o
PTHREAD_LIBMDL 		= $(PTHREAD_MDL)/mdl.o -lm -lpthread

#
#       LINUX Charm  defines
#
CHARM=../charm/net-linux/bin/charmc
CHARMLINK=
CHARM_MDL			= ../mdl/charm
CHARM_CFLAGS		= -verbose -g -I$(CHARM_MDL) $(BASE_DEF)
CHARM_MDL_CFLAGS	= -verbose -g -Wall
CHARM_LD_FLAGS		=  $(BASE_LD_FLAGS) -language charm++ -memory os
CHARM_XOBJ                = erf.o v_sqrt1.o
CHARM_LIBMDL              = $(CHARM_MDL)/mdl.o $(CHARMLINK) -lm

OBJ	= main.o master.o param.o outtype.o pkd.o pst.o grav.o \
	  ewald.o walk.o eccanom.o fdl.o htable.o smooth.o \
	  smoothfcn.o collision.o qqsmooth.o $(COOLING_OBJ) cosmo.o romberg.o \
	  starform.o feedback.o millerscalo.o supernova.o supernovaia.o \
	  startime.o stiff.o runge.o dumpframe.o dffuncs.o dumpvoxel.o \
	  rotbar.o special.o ssio.o $(PNG_OBJ) $(TURB_OBJ) \
	  treezip.o log.o linalg.o aggs.o

EXTRA_OBJ = erf.o v_sqrt1.o 


default:
	@echo "Please tell me what architecture to make."
	@echo "Choices are null, pvm, pthread, ampi, lam_mpi, qmpi, mpi, charm"

install:
	@echo "No installation rules."

clean:
	-rm -f $(OBJ) $(EXTRA_OBJ)

spotless: clean
	@echo rm -rf $(EXE)
	@rm -f $(EXE)
	@rm -f ../mdl/*/mdl.o

depend:
	makedepend -Y -- $(BASE_DEF) -- *.c

$(XDIR):
	-mkdir $(BDIR)
	-mkdir $(XDIR)

null:
	cd $(NULL_MDL); make "CC=$(CC)" "CFLAGS=$(NULL_CFLAGS)"
	make $(EXE) "CFLAGS=$(NULL_CFLAGS)" "LD_FLAGS=$(NULL_LD_FLAGS)"\
		"MDL=$(NULL_MDL)" "XOBJ=$(NULL_XOBJ)" "LIBMDL=$(NULL_LIBMDL)"
dbg:
	cd $(DBG_MDL); make "CC=$(CC)" "CFLAGS=$(DBG_CFLAGS)"
	make $(EXE) "CFLAGS=$(DBG_CFLAGS)" "LD_FLAGS=$(DBG_LD_FLAGS)"\
		"MDL=$(DBG_MDL)" "XOBJ=$(DBG_XOBJ)" "LIBMDL=$(DBG_LIBMDL)"

pvm:
	cd $(PVM_MDL); aimk
	make $(EXE) "CFLAGS=$(PVM_CFLAGS)" "LD_FLAGS=$(PVM_LD_FLAGS)"\
		"MDL=$(PVM_MDL)" "XOBJ=$(PVM_XOBJ)" "LIBMDL=$(PVM_LIBMDL)"
	mv -f $(EXE) $(XDIR)

pthread:
	cd $(PTHREAD_MDL); make "CC=$(CC)" "CFLAGS=$(PTHREAD_CFLAGS)"
	make $(EXE) "EXE=$(EXE)" "CFLAGS=$(PTHREAD_CFLAGS)" "LD_FLAGS=$(PTHREAD_LD_FLAGS)"\
		"MDL=$(PTHREAD_MDL)" "XOBJ=$(PTHREAD_XOBJ)" "LIBMDL=$(PTHREAD_LIBMDL)"

ampi:
	cd $(AMPI_MDL); make CC=$(CHARM) "CFLAGS=$(AMPI_MDL_CFLAGS)"
	make $(EXE) CC=$(CHARM) "CFLAGS=$(AMPI_CFLAGS)" "LD_FLAGS=$(AMPI_LD_FLAGS)"\
		"MDL=$(AMPI_MDL)" "XOBJ=$(AMPI_XOBJ)" "LIBMDL=$(AMPI_LIBMDL)"

lam_mpi:
	cd $(LAM_MDL); make CC=mpicc "CC_FLAGS=$(LAM_MDL_CFLAGS)"
	make $(EXE) CC=mpicc "CFLAGS=$(LAM_CFLAGS)" "LD_FLAGS=$(LAM_LD_FLAGS)"\
		"MDL=$(LAM_MDL)" "XOBJ=$(LAM_XOBJ)" "LIBMDL=$(LAM_LIBMDL)"

qmpi:
	cd $(QMPI_MDL); make "CC_FLAGS=$(QMPI_MDL_CFLAGS)"
	make $(EXE) "CFLAGS=$(QMPI_CFLAGS)" "LD_FLAGS=$(QMPI_LD_FLAGS)"\
		"MDL=$(QMPI_MDL)" "XOBJ=$(QMPI_XOBJ)" "LIBMDL=$(QMPI_LIBMDL)"

mpi:
	cd $(MPI_MDL); make CC=mpicc "CC_FLAGS=$(MPI_MDL_CFLAGS)"
	make $(EXE) CC=mpicc "CFLAGS=$(MPI_CFLAGS)" "LD_FLAGS=$(MPI_LD_FLAGS)"\
		"MDL=$(MPI_MDL)" "XOBJ=$(MPI_XOBJ)" "LIBMDL=$(MPI_LIBMDL)"

charm:
	cd $(CHARM_MDL); make CC=$(CHARM) "CFLAGS=$(CHARM_MDL_CFLAGS)"
	make $(EXE) CC=$(CHARM) "CFLAGS=$(CHARM_CFLAGS)" "LD_FLAGS=$(CHARM_LD_FLAGS)"\
		"MDL=$(CHARM_MDL)" "XOBJ=$(CHARM_XOBJ)" "LIBMDL=$(CHARM_LIBMDL)"

$(EXE): $(OBJ) $(XOBJ)
	$(CC) $(CFLAGS) $(LD_FLAGS) -o $@ $(OBJ) $(XOBJ) $(LIBMDL)

$(OBJ) $(EXTRA_OBJ): Makefile

# DO NOT DELETE
cosmo.o: runge.h cosmo.h
ewald.o: ewald.h pkd.h floattype.h cooling.h meval.h qeval.h
fdl.o: htable.h fdl.h
grav.o: pkd.h floattype.h cooling.h grav.h meval.h qeval.h
htable.o: htable.h
integration.o: linalg.h floattype.h integration.h
aggs.o: aggs.c aggs.h
linalg.o: linalg.h floattype.h
main.o: master.h param.h pst.h pkd.h floattype.h cooling.h smoothfcn.h
main.o: parameters.h cosmo.h outtype.h
master.o: master.h param.h pst.h pkd.h floattype.h cooling.h smoothfcn.h
master.o: parameters.h cosmo.h tipsydefs.h opentype.h fdl.h htable.h
master.o: outtype.h
millerscalo.o: millerscalo.h
outtype.o: pkd.h floattype.h cooling.h outtype.h
outurb.o: pkd.h outurb.h
param.o: param.h
pkd.o: pkd.h floattype.h cooling.h ewald.h grav.h walk.h opentype.h
pkd.o: tipsydefs.h dumpframe.h
pst.o: pst.h pkd.h floattype.h cooling.h smoothfcn.h starform.h feedback.h
pst.o: outtype.h smooth.h dumpframe.h
romberg.o: floattype.h
smooth.o: smooth.h pkd.h floattype.h cooling.h smoothfcn.h
smoothfcn.o: smoothfcn.h pkd.h floattype.h cooling.h
supernova.o: pkd.h floattype.h cooling.h feedback.h supernova.h startime.h
supernova.o: millerscalo.h supernovaia.h
supernovaia.o: supernova.h startime.h millerscalo.h feedback.h pkd.h
supernovaia.o: floattype.h cooling.h supernovaia.h
dumpframe.o: dumpframe.c dumpframe.h
dumpvoxel.o: dumpvoxel.c dumpframe.h dumpvoxel.h
dffuncs.o: dffuncs.c dumpframe.h
treezip.o: treezip.c treezip.h treezipkey.h treeziptypes.h
writepng.o: writepng.c writepng.h
walk.o: walk.h pkd.h floattype.h cooling.h
rotbar.o: master.h rotbar.h
cooling_metal.o: cooling_metal.h cooling_metal.c
cooling_grackle.o: cooling_grackle.h cooling_grackle.c
log.o: log.h
