#
# Makefile for pkdgrav or gasoline (specify below)
#

#CC = gcc -Wall
EXE = pkdgrav
#EXE = gasoline
CODEDEF = -DCOLLISIONS -DSLIDING_PATCH
#CODEDEF = -DSUPERCOOL
#CODEDEF = -DGROWMASS
#CODEDEF = -DGASOLINE
#CODEDEF = -DDEBUG -HSHRINK

#ev6
#CC = ccc
#CODEDEF = -DGASOLINE -DCCC

CC = cc
# -DSUPERNOVA
CODEDEF = -DGASOLINE 

#
#       NULL defines
#
NULL_MDL		= ../mdl/null
NULL_CFLAGS		= -O3 -I$(NULL_MDL) $(CODEDEF)

#ev6 flags:
#NULL_CFLAGS		= -O3 -g3 -fast -arch ev6 -I$(NULL_MDL) $(CODEDEF)
#NULL_CFLAGS		= -O3 -fast -arch ev6 -I$(NULL_MDL) $(CODEDEF)

NULL_CFLAGS		= -O2 -I$(NULL_MDL) $(CODEDEF)
#NULL_CFLAGS		= -O2 -I$(NULL_MDL) $(CODEDEF)
#NULL_CFLAGS		= -g -I$(NULL_MDL) $(CODEDEF)
#NULL_LD_FLAGS	= -Wl,-s
NULL_XOBJ		= erf.o v_sqrt1.o
NULL_LIBMDL		= $(NULL_MDL)/mdl.o -lm

#
#       SGI defines
#
SGI_MDL			= ../mdl/mpi
SGI_CFLAGS		= -O2 -I$(SGI_MDL) $(CODEDEF) -mips4 -64 -r10000
SGI_LD_FLAGS	= -mips4 -64 -r10000
SGI_XOBJ		=
SGI_LIBMDL		= $(SGI_MDL)/mdl.o -lmpi -lm
SGI_MDL_CFLAGS	= -O2 -mips4 -64 -r10000

#
#       LINUX LAM MPI defines
#
LAM_MDL			= ../mdl/mpi
LAM_CFLAGS		= -O3 -I$(LAM_MDL) $(CODEDEF) -DMPI_LINUX -I/net/lam-6.3-b1/include
LAM_LD_FLAGS	= -L/net/lam-6.3-b1/lib -lmpi -ltstdio -lt -largs -ltrillium -ltstdio -lmpi++
LAM_XOBJ		=
LAM_LIBMDL		= $(LAM_MDL)/mdl.o -L/net/lam-6.3-b1/lib -lmpi -ltstdio -lt -largs -ltrillium -ltstdio -lmpi++
LAM_MDL_CFLAGS	= -O3 -I$(LAM_MDL) $(CODEDEF) -DMPI_LINUX  -I/net/lam-6.3-b1/include 

#
#       SP1/2 defines
#
SPX_MDL			= ../mdl/mpi
SPX_CFLAGS		= -O3 -I$(SPX_MDL) $(CODEDEF)
SPX_LD_FLAGS	=
SPX_XOBJ		= v_sqrt1.o
SPX_LIBMDL		= $(SPX_MDL)/mdl.o -lm
SPX_MDL_CFLAGS	= -g -O3

#
#		PVM defines
#
PVMDIR	= $(PVM_ROOT)
PVMLIB	= $(PVMDIR)/lib/$(PVM_ARCH)/libpvm3.a
BDIR	= $(HOME)/pvm3/bin
XDIR	= $(BDIR)/$(PVM_ARCH)

PVM_MDL		= ../mdl/pvm
PVM_CFLAGS	= -O3 -I$(PVMDIR)/include -I$(PVM_MDL) $(CODEDEF)
#PVM_CFLAGS	= -mips4 -g -I$(PVMDIR)/include -I$(PVM_MDL) $(CODEDEF)
PVM_XOBJ	= v_sqrt1.o
PVM_LIBMDL	= $(PVM_MDL)/$(PVM_ARCH)/mdl.o $(PVMLIB) $(ARCHLIB) -lm

#
#       PTHREAD defines
#
PTHREAD_MDL			= ../mdl/pthread
PTHREAD_CFLAGS		= -O3 -D_REENTRANT -I$(PTHREAD_MDL) $(CODEDEF)
PTHREAD_LD_FLAGS 	=
PTHREAD_XOBJ		= erf.o v_sqrt1.o
PTHREAD_LIBMDL 		= $(PTHREAD_MDL)/mdl.o -lm -lpthread

#
#       PTHREAD_SGI defines
#
PTHREAD_SGI_MDL			= ../mdl/pthread
PTHREAD_SGI_CFLAGS		= -O2 -D_REENTRANT -I$(PTHREAD_SGI_MDL) $(CODEDEF) -mips4 -64 -r10000
PTHREAD_SGI_LD_FLAGS 	= -mips4 -64 -r10000
PTHREAD_SGI_XOBJ		= 
PTHREAD_SGI_LIBMDL 		= $(PTHREAD_SGI_MDL)/mdl.o -lm -lpthread
PTHREAD_SGI_MDL_CFLAGS	= -O2 -mips4 -64 -r10000

#
#       T3D MPP defines
#
T3D_MDL	= ../mdl/mpp
V_SQRT	= ../v_sqrt/lib
V_SQRT1	= ../v_sqrt1/lib
RPC		= ../rpc

T3D_CFLAGS		= -O3 -g -DCRAY_T3D -I$(T3D_MDL) -I$(V_SQRT) -I$(RPC) $(CODEDEF)
T3D_XOBJ		= hyperlib.o
T3D_LIBMDL		= -O3 -g-L $(V_SQRT) -L $(V_SQRT1) -L $(RPC) \
					$(T3D_MDL)/mdl.o -lv_sqrtc -lv_sqrtc1 -lmpi -lrpc -lm
T3D_LD_FLAGS	=

#
#       T3DMPI MPP defines
#
T3DMPI_MDL	= ../mdl/t3dmpi
V_SQRT		= ../v_sqrt/lib
V_SQRT1		= ../v_sqrt1/lib
RPC			= ../rpc

T3DMPI_CFLAGS	= -O3 -DCRAY_T3D -I$(T3DMPI_MDL) -I$(V_SQRT) -I$(RPC) $(CODEDEF)
T3DMPI_XOBJ		= hyperlib.o
T3DMPI_LIBMDL	= -O3 -L $(V_SQRT) -L $(V_SQRT1) -L $(RPC) \
					$(T3DMPI_MDL)/mdl.o -lv_sqrtc -lv_sqrtc1 -lmpi -lrpc -lm
T3DMPI_LD_FLAGS	=

#
#       T3EMPI MPP defines
#
T3EMPI_MDL		= ../mdl/t3empi
T3EMPI_CFLAGS	= -O3 -DCRAY_T3D -I$(T3EMPI_MDL) -I$(RPC) $(CODEDEF)
T3EMPI_XOBJ		= hyperlib.o v_sqrt1.t3x.o
T3EMPI_LIBMDL	= $(T3EMPI_MDL)/mdl.o -lmpi -lm
T3DMPI_LD_FLAGS	=

#
#       KSR1 defines
#
KSR_MDL			= ../mdl/ksr
KSR_CFLAGS		= -O2 -w2 -I$(KSR_MDL) $(CODEDEF)
KSR_LD_FLAGS	= -para
KSR_XOBJ		= erf.o v_sqrt1.ksr.o
KSR_LIBMDL		= $(KSR_MDL)/mdl.o -lm -lrpc

OBJ	= main.o master.o param.o outtype.o pkd.o pst.o grav.o \
	  ewald.o walk.o eccanom.o hypanom.o fdl.o htable.o smooth.o \
	  smoothfcn.o collision.o qqsmooth.o cooling.o cosmo.o romberg.o

EXTRA_OBJ = erf.o hyperlib.o v_sqrt1.o v_sqrt1.ksr.o v_sqrt1.t3x.o

default:	
	@echo "Please tell me what architecture to make."
	@echo "Choices are null, sgi, pvm, pthread, spx, t3d, t3dmpi, t3empi, ksr."

install:
	@echo "No installation rules."

clean:
	-rm -f $(OBJ) $(EXTRA_OBJ)

spotless: clean
	-rm -f $(EXE)

depend:
	makedepend -Y -- $(CODEDEF) -- *.c

$(XDIR):
	-mkdir $(BDIR)
	-mkdir $(XDIR)

null:
	cd $(NULL_MDL); make "CC=$(CC)" "CFLAGS=$(NULL_CFLAGS)"
	make $(EXE) "CFLAGS=$(NULL_CFLAGS)" "LD_FLAGS=$(NULL_LD_FLAGS)"\
		"MDL=$(NULL_MDL)" "XOBJ=$(NULL_XOBJ)" "LIBMDL=$(NULL_LIBMDL)"

sgi:
	cd $(SGI_MDL); make CC=cc "CC_FLAGS=$(SGI_MDL_CFLAGS)"
	make $(EXE) "CFLAGS=$(SGI_CFLAGS)" "LD_FLAGS=$(SGI_LD_FLAGS)"\
		"MDL=$(SGI_MDL)" "XOBJ=$(SGI_XOBJ)" "LIBMDL=$(SGI_LIBMDL)"

pvm:
	cd $(PVM_MDL); aimk	
	make $(EXE) "CFLAGS=$(PVM_CFLAGS)" "LD_FLAGS=$(PVM_LD_FLAGS)"\
		"MDL=$(PVM_MDL)" "XOBJ=$(PVM_XOBJ)" "LIBMDL=$(PVM_LIBMDL)"
	mv -f $(EXE) $(XDIR)

pthread:
	cd $(PTHREAD_MDL); make
	make $(EXE) "CFLAGS=$(PTHREAD_CFLAGS)" "LD_FLAGS=$(PTHREAD_LD_FLAGS)"\
		"MDL=$(PTHREAD_MDL)" "XOBJ=$(PTHREAD_XOBJ)" "LIBMDL=$(PTHREAD_LIBMDL)"

pthread_sgi:
	cd $(PTHREAD_MDL); make CC=cc "CC_FLAGS=$(PTHREAD_SGI_MDL_CFLAGS)"
	make $(EXE) "CFLAGS=$(PTHREAD_SGI_CFLAGS)" "LD_FLAGS=$(PTHREAD_SGI_LD_FLAGS)"\
		"MDL=$(PTHREAD_SGI_MDL)" "XOBJ=$(PTHREAD_SGI_XOBJ)" "LIBMDL=$(PTHREAD_SGI_LIBMDL)"

lam_mpi:
	cd $(LAM_MDL); make CC=pgcc "CC_FLAGS=$(LAM_MDL_CFLAGS)"
	make $(EXE) CC=pgcc "CFLAGS=$(LAM_CFLAGS)" "LD_FLAGS=$(LAM_LD_FLAGS)"\
		"MDL=$(LAM_MDL)" "XOBJ=$(LAM_XOBJ)" "LIBMDL=$(LAM_LIBMDL)"


mpi: spx

spx:
	cd $(SPX_MDL); make CC=mpicc "CC_FLAGS=$(SPX_MDL_CFLAGS)"
	make $(EXE) CC=mpicc "CFLAGS=$(SPX_CFLAGS)" "LD_FLAGS=$(SPX_LD_FLAGS)"\
		"MDL=$(SPX_MDL)" "XOBJ=$(SPX_XOBJ)" "LIBMDL=$(SPX_LIBMDL)"

t3d:
	cd $(T3D_MDL); make
	make $(EXE) "CFLAGS=$(T3D_CFLAGS)" "LD_FLAGS=$(T3D_LD_FLAGS)"\
		"MDL=$(T3D_MDL)" "XOBJ=$(T3D_XOBJ)" "LIBMDL=$(T3D_LIBMDL)"

t3dmpi:
	cd $(T3DMPI_MDL); make
	make $(EXE) "CFLAGS=$(T3DMPI_CFLAGS)" "LD_FLAGS=$(T3DMPI_LD_FLAGS)"\
		"MDL=$(T3DMPI_MDL)" "XOBJ=$(T3DMPI_XOBJ)" "LIBMDL=$(T3DMPI_LIBMDL)"

t3empi:
	cd $(T3EMPI_MDL); make
	make $(EXE) "CFLAGS=$(T3EMPI_CFLAGS)" "LD_FLAGS=$(T3EMPI_LD_FLAGS)"\
		"MDL=$(T3EMPI_MDL)" "XOBJ=$(T3EMPI_XOBJ)" "LIBMDL=$(T3EMPI_LIBMDL)"

ksr:
	cd $(KSR_MDL); make
	make $(EXE) "CFLAGS=$(KSR_CFLAGS)" "LD_FLAGS=$(KSR_LD_FLAGS)"\
		"MDL=$(KSR_MDL)" "XOBJ=$(KSR_XOBJ)" "LIBMDL=$(KSR_LIBMDL)"

$(EXE): $(OBJ) $(XOBJ)
	$(CC) $(CFLAGS) $(LD_FLAGS) -o $@ $(OBJ) $(XOBJ) $(LIBMDL)

$(OBJ) $(EXTRA_OBJ): Makefile

# DO NOT DELETE

collision.o: pkd.h floattype.h collision.h
ewald.o: ewald.h pkd.h floattype.h meval.h qeval.h
fdl.o: htable.h fdl.h
grav.o: pkd.h floattype.h grav.h meval.h qeval.h
htable.o: htable.h
main.o: pst.h pkd.h floattype.h smoothfcn.h collision.h master.h param.h
main.o: parameters.h outtype.h
master.o: master.h param.h pst.h pkd.h floattype.h smoothfcn.h collision.h
master.o: parameters.h tipsydefs.h opentype.h fdl.h htable.h outtype.h
master.o: ssdefs.h ssio.h
outtype.o: pkd.h floattype.h outtype.h collision.h
param.o: param.h
pkd.o: pkd.h floattype.h ewald.h grav.h walk.h opentype.h tipsydefs.h
pkd.o: ssdefs.h ssio.h collision.h
pst.o: pst.h pkd.h floattype.h smoothfcn.h collision.h outtype.h smooth.h
qqsmooth.o: smooth.h pkd.h floattype.h smoothfcn.h
smooth.o: smooth.h pkd.h floattype.h smoothfcn.h
smoothfcn.o: smoothfcn.h pkd.h floattype.h ssdefs.h ssio.h collision.h
walk.o: walk.h pkd.h floattype.h
cooling.o: cooling.h





