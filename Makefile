#
# Makefile for pkdgrav
#
# PVM defines
#
PVMDIR	=	$(PVM_ROOT)
PVMLIB	=	$(PVMDIR)/lib/$(PVM_ARCH)/libpvm3.a
BDIR	=	$(HOME)/pvm3/bin
XDIR	=	$(BDIR)/$(PVM_ARCH)

PVM_MDL		=	../mdl/pvm
PVM_CFLAGS	=	-mips4 -O3 -I$(PVMDIR)/include -I$(PVM_MDL)
#PVM_CFLAGS	=	-mips4 -g -I$(PVMDIR)/include -I$(PVM_MDL)
PVM_LIBMDL	=	v_sqrt1.o $(PVM_MDL)/$(PVM_ARCH)/mdl.o $(PVMLIB) $(ARCHLIB) /local/lib/libmalloc.a -lm

#
#       KSR1 defines
#
KSR_MDL = ../mdl/ksr

KSR_CFLAGS = -O2 -w2 -I$(KSR_MDL)
KSR_LD_FLAGS = -para
KSR_LIBMDL = erf.o v_sqrt1.ksr.o $(KSR_MDL)/mdl.o -lm -lrpc

#
#       SP1/2 defines
#
SPX_MDL = ../mdl/mpi
SPX_CFLAGS = -O3 -I$(SPX_MDL)
SPX_LD_FLAGS = 
SPX_LIBMDL = v_sqrt1.o $(SPX_MDL)/mdl.o -lm

#
#       T3D MPP defines
#
T3D_MDL = ../mdl/mpp
V_SQRT = ../v_sqrt/lib
V_SQRT1 = ../v_sqrt1/lib
RPC = ../rpc

T3D_CFLAGS = -O3 -g -DCRAY_T3D -I$(T3D_MDL) -I$(V_SQRT) -I$(RPC)
T3D_LIBMDL = -O3 -g-L $(V_SQRT) -L $(V_SQRT1) -L $(RPC) hyperlib.o \
	$(T3D_MDL)/mdl.o -lv_sqrtc -lv_sqrtc1 -lmpi -lrpc -lm
T3D_LD_FLAGS = 

#
#       T3DMPI MPP defines
#
T3DMPI_MDL = ../mdl/t3dmpi
V_SQRT = ../v_sqrt/lib
V_SQRT1 = ../v_sqrt1/lib
RPC = ../rpc

T3DMPI_CFLAGS = -O3 -DCRAY_T3D -I$(T3DMPI_MDL) -I$(V_SQRT) -I$(RPC)
T3DMPI_LIBMDL = -O3 -L $(V_SQRT) -L $(V_SQRT1) -L $(RPC) hyperlib.o \
	$(T3DMPI_MDL)/mdl.o -lv_sqrtc -lv_sqrtc1 -lmpi -lrpc -lm
T3DMPI_LD_FLAGS = 

OBJS	= 	main.o master.o param.o outtype.o pkd.o pst.o grav.o \
		ewald.o walk.o eccanom.o hypanom.o fdl.o htable.o smooth.o
EXTRA_OBJ = 	erf.o v_sqrt1.o v_sqrt1.ksr.o hyperlib.o

default:	
	@echo "Please tell me what architecture to make."
	@echo "Choices are pvm, ksr, spx, t3d and t3dmpi."

$(XDIR):
	- mkdir $(BDIR)
	- mkdir $(XDIR)

clean:
	rm -f *.o

pvm:
	cd $(PVM_MDL); aimk	
	make pkdgrav "CFLAGS=$(PVM_CFLAGS)" "LD_FLAGS=$(PVM_LD_FLAGS)" "MDL=$(PVM_MDL)" "LIBMDL=$(PVM_LIBMDL)"
	mv -f pkdgrav $(XDIR)


ksr:
	cd $(KSR_MDL); make
	make pkdgrav "CFLAGS=$(KSR_CFLAGS)" "LD_FLAGS=$(KSR_LD_FLAGS)" "MDL=$(KSR_MDL)" "LIBMDL=$(KSR_LIBMDL)"

spx:
	cd $(SPX_MDL); make
	make pkdgrav CC=mpicc "CFLAGS=$(SPX_CFLAGS)" "LD_FLAGS=$(SPX_LD_FLAGS)" "MDL=$(SPX_MDL)" "LIBMDL=$(SPX_LIBMDL)"

t3d:
	cd $(T3D_MDL); make
	make pkdgrav "CFLAGS=$(T3D_CFLAGS)" "LD_FLAGS=$(T3D_LD_FLAGS)" "MDL=$(T3D_MDL)" "LIBMDL=$(T3D_LIBMDL)"

t3dmpi:
	cd $(T3DMPI_MDL); make
	make pkdgrav "CFLAGS=$(T3DMPI_CFLAGS)" "LD_FLAGS=$(T3DMPI_LD_FLAGS)" "MDL=$(T3DMPI_MDL)" "LIBMDL=$(T3DMPI_LIBMDL)"

pkdgrav: $(OBJS) $(EXTRA_OBJ)
	$(CC) $(CFLAGS) $(LD_FLAGS) -o pkdgrav $(OBJS) $(LIBMDL)

###
ewald.o: ewald.h pkd.h meval.h qeval.h
grav.o: pkd.h grav.h meval.h qeval.h
main.o: pst.h pkd.h master.h param.h parameters.h outtype.h
master.o: master.h param.h pst.h pkd.h parameters.h \
  tipsydefs.h opentype.h checkdefs.h
outtype.o: pkd.h outtype.h
param.o: param.c param.h
pkd.o: pkd.h ewald.h grav.h walk.h opentype.h tipsydefs.h \
  checkdefs.h parameters.h
pst.o: pst.h pkd.h outtype.h smooth.h
smooth.o: smooth.c smooth.h pkd.h
walk.o: walk.h pkd.h
