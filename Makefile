#
# Makefile for pkdgrav or gasoline (specify below)
#
#CC = gcc -Wall
EXE = pkdgrav
#CODEDEF = -DPLANETS
#CODEDEF = -DSUPERCOOL
#EXE = gasoline
#CODEDEF = -DGASOLINE
#
# PVM defines
#
PVMDIR	=	$(PVM_ROOT)
PVMLIB	=	$(PVMDIR)/lib/$(PVM_ARCH)/libpvm3.a
BDIR	=	$(HOME)/pvm3/bin
XDIR	=	$(BDIR)/$(PVM_ARCH)

PVM_MDL		=	../mdl/pvm
PVM_CFLAGS	=	-mips4 -O3 -I$(PVMDIR)/include -I$(PVM_MDL) $(CODEDEF)
#PVM_CFLAGS	=	-mips4 -g -I$(PVMDIR)/include -I$(PVM_MDL)
PVM_LIBMDL	=	v_sqrt1.o $(PVM_MDL)/$(PVM_ARCH)/mdl.o $(PVMLIB) $(ARCHLIB) /local/lib/libmalloc.a -lm

#
#       NULL MDL defines
#
NULL_MDL = ../mdl/null

NULL_CFLAGS = -O3 -I$(NULL_MDL) $(CODEDEF)
#NULL_CFLAGS = -O3 -g -malign-double -I$(NULL_MDL) $(CODEDEF)
NULL_LD_FLAGS = 
NULL_LIBMDL = erf.o v_sqrt1.o $(NULL_MDL)/mdl.o -lm

#
#       PTHREAD defines
#
PTHREAD_MDL = ../mdl/pthread

PTHREAD_CFLAGS = -O3 -malign-double -D_REENTRANT -I$(PTHREAD_MDL) $(CODEDEF)
PTHREAD_LD_FLAGS =
PTHREAD_LIBMDL = erf.o v_sqrt1.o $(PTHREAD_MDL)/mdl.o -lm -lpthread

#
#       KSR1 defines
#
KSR_MDL = ../mdl/ksr

KSR_CFLAGS = -O2 -w2 -I$(KSR_MDL) $(CODEDEF)
KSR_LD_FLAGS = -para
KSR_LIBMDL = erf.o v_sqrt1.ksr.o $(KSR_MDL)/mdl.o -lm -lrpc

#
#       SP1/2 defines
#
SPX_MDL = ../mdl/mpi
SPX_CFLAGS = -O3 -I$(SPX_MDL) $(CODEDEF)
SPX_LD_FLAGS = 
SPX_LIBMDL = v_sqrt1.o $(SPX_MDL)/mdl.o -lm

#
#       T3D MPP defines
#
T3D_MDL = ../mdl/mpp
V_SQRT = ../v_sqrt/lib
V_SQRT1 = ../v_sqrt1/lib
RPC = ../rpc

T3D_CFLAGS = -O3 -g -DCRAY_T3D -I$(T3D_MDL) -I$(V_SQRT) -I$(RPC) $(CODEDEF)
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

T3DMPI_CFLAGS = -O3 -DCRAY_T3D -I$(T3DMPI_MDL) -I$(V_SQRT) -I$(RPC) $(CODEDEF)
T3DMPI_LIBMDL = -O3 -L $(V_SQRT) -L $(V_SQRT1) -L $(RPC) hyperlib.o \
	$(T3DMPI_MDL)/mdl.o -lv_sqrtc -lv_sqrtc1 -lmpi -lrpc -lm
T3DMPI_LD_FLAGS = 

#
#       T3EMPI MPP defines
#
T3EMPI_MDL = ../mdl/t3empi

T3EMPI_CFLAGS = -O3 -DCRAY_T3D -I$(T3EMPI_MDL) -I$(RPC) $(CODEDEF)
T3EMPI_LIBMDL = hyperlib.o v_sqrt1.t3x.o \
	$(T3EMPI_MDL)/mdl.o -lmpi -lm
T3DMPI_LD_FLAGS = 

OBJS	= 	main.o master.o param.o outtype.o pkd.o pst.o grav.o \
		ewald.o walk.o eccanom.o hypanom.o fdl.o htable.o smooth.o \
		smoothfcn.o
EXTRA_OBJ = 	erf.o v_sqrt1.o v_sqrt1.ksr.o hyperlib.o

default:	
	@echo "Please tell me what architecture to make."
	@echo "Choices are pvm, null, pthread, ksr, spx, t3d and t3dmpi."

$(XDIR):
	- mkdir $(BDIR)
	- mkdir $(XDIR)

clean:
	rm -f *.o

depend:
	makedepend -- $(NULL_CFLAGS) -- *.c

cflow:
	-rm -f cflow.out
	cflow -I$(NULL_MDL) *.c > cflow.out

spotless:
	-rm -f $(EXE) $(OBJS) $(EXTRA_OBJ) cflow.out

pvm:
	cd $(PVM_MDL); aimk	
	make $(EXE) "CFLAGS=$(PVM_CFLAGS)" "LD_FLAGS=$(PVM_LD_FLAGS)" "MDL=$(PVM_MDL)" "LIBMDL=$(PVM_LIBMDL)"
	mv -f $(EXE) $(XDIR)

null:
	cd $(NULL_MDL); make
	make $(EXE) "CFLAGS=$(NULL_CFLAGS)" "LD_FLAGS=$(NULL_LD_FLAGS)" "MDL=$(NULL_MDL)" "LIBMDL=$(NULL_LIBMDL)"

pthread:
	cd $(PTHREAD_MDL); make
	make $(EXE) "CFLAGS=$(PTHREAD_CFLAGS)" "LD_FLAGS=$(PTHREAD_LD_FLAGS)" "MDL=$(PTHREAD_MDL)" "LIBMDL=$(PTHREAD_LIBMDL)"

ksr:
	cd $(KSR_MDL); make
	make $(EXE) "CFLAGS=$(KSR_CFLAGS)" "LD_FLAGS=$(KSR_LD_FLAGS)" "MDL=$(KSR_MDL)" "LIBMDL=$(KSR_LIBMDL)"

spx:
	cd $(SPX_MDL); make
	make $(EXE) CC=mpicc "CFLAGS=$(SPX_CFLAGS)" "LD_FLAGS=$(SPX_LD_FLAGS)" "MDL=$(SPX_MDL)" "LIBMDL=$(SPX_LIBMDL)"

t3d:
	cd $(T3D_MDL); make
	make $(EXE) "CFLAGS=$(T3D_CFLAGS)" "LD_FLAGS=$(T3D_LD_FLAGS)" "MDL=$(T3D_MDL)" "LIBMDL=$(T3D_LIBMDL)"

t3dmpi:
	cd $(T3DMPI_MDL); make
	make $(EXE) "CFLAGS=$(T3DMPI_CFLAGS)" "LD_FLAGS=$(T3DMPI_LD_FLAGS)" "MDL=$(T3DMPI_MDL)" "LIBMDL=$(T3DMPI_LIBMDL)"

t3empi:
	cd $(T3EMPI_MDL); make
	make $(EXE) "CFLAGS=$(T3EMPI_CFLAGS)" "LD_FLAGS=$(T3EMPI_LD_FLAGS)" "MDL=$(T3EMPI_MDL)" "LIBMDL=$(T3EMPI_LIBMDL)"

pkdgrav: $(OBJS) $(EXTRA_OBJ)
	$(CC) $(CFLAGS) $(LD_FLAGS) -o pkdgrav $(OBJS) $(LIBMDL)

gasoline: $(OBJS) $(EXTRA_OBJ)
	$(CC) $(CFLAGS) $(LD_FLAGS) -DGASOLINE -o gasoline $(OBJS) $(LIBMDL)

# DO NOT DELETE

eccanom.o: /usr/include/math.h /usr/include/assert.h
erf.o: /usr/include/math.h
ewald.o: /usr/include/stdio.h /usr/include/malloc.h /usr/include/math.h
ewald.o: /usr/include/assert.h ewald.h pkd.h /usr/include/sys/resource.h
ewald.o: /usr/include/sys/time.h /usr/include/sys/types.h
ewald.o: /usr/include/sgidefs.h /usr/include/sys/bsd_types.h
ewald.o: /usr/include/sys/select.h ../mdl/null/mdl.h meval.h qeval.h
fdl.o: /usr/include/stdio.h /usr/include/stdlib.h /usr/include/sgidefs.h
fdl.o: /usr/include/malloc.h /usr/include/assert.h /usr/include/ctype.h
fdl.o: /usr/include/string.h htable.h fdl.h
grav.o: /usr/include/stdio.h /usr/include/math.h /usr/include/assert.h
grav.o: ../mdl/null/mdl.h pkd.h /usr/include/sys/resource.h
grav.o: /usr/include/sys/time.h /usr/include/sys/types.h
grav.o: /usr/include/sgidefs.h /usr/include/sys/bsd_types.h
grav.o: /usr/include/sys/select.h grav.h meval.h qeval.h
htable.o: /usr/include/stdlib.h /usr/include/sgidefs.h /usr/include/malloc.h
htable.o: htable.h
hypanom.o: /usr/include/math.h /usr/include/assert.h
hyperlib.o: /usr/include/math.h
main.o: /usr/include/stdio.h /usr/include/math.h /usr/include/stdlib.h
main.o: /usr/include/sgidefs.h /usr/include/string.h /usr/include/assert.h
main.o: ../mdl/null/mdl.h pst.h pkd.h /usr/include/sys/resource.h
main.o: /usr/include/sys/time.h /usr/include/sys/types.h
main.o: /usr/include/sys/bsd_types.h /usr/include/sys/select.h smoothfcn.h
main.o: master.h param.h parameters.h outtype.h
master.o: /usr/include/stdio.h /usr/include/stdlib.h /usr/include/sgidefs.h
master.o: /usr/include/stddef.h /usr/include/string.h /usr/include/malloc.h
master.o: /usr/include/assert.h /usr/include/time.h /usr/include/math.h
master.o: /usr/include/rpc/types.h /usr/include/sys/types.h
master.o: /usr/include/sys/bsd_types.h /usr/include/sys/select.h
master.o: /usr/include/sys/time.h /usr/include/rpc/xdr.h master.h param.h
master.o: pst.h pkd.h /usr/include/sys/resource.h ../mdl/null/mdl.h
master.o: smoothfcn.h parameters.h tipsydefs.h opentype.h fdl.h htable.h
master.o: outtype.h
outtype.o: /usr/include/stdio.h /usr/include/math.h /usr/include/assert.h
outtype.o: pkd.h /usr/include/sys/resource.h /usr/include/sys/time.h
outtype.o: /usr/include/sys/types.h /usr/include/sgidefs.h
outtype.o: /usr/include/sys/bsd_types.h /usr/include/sys/select.h
outtype.o: ../mdl/null/mdl.h outtype.h
param.o: /usr/include/stdio.h /usr/include/stdlib.h /usr/include/sgidefs.h
param.o: /usr/include/stddef.h /usr/include/malloc.h /usr/include/string.h
param.o: /usr/include/assert.h /usr/include/ctype.h param.h
pkd.o: /usr/include/stdio.h /usr/include/stdlib.h /usr/include/sgidefs.h
pkd.o: /usr/include/stddef.h /usr/include/malloc.h /usr/include/math.h
pkd.o: /usr/include/assert.h /usr/include/sys/time.h /usr/include/sys/types.h
pkd.o: /usr/include/sys/bsd_types.h /usr/include/sys/select.h
pkd.o: /usr/include/rpc/types.h /usr/include/rpc/xdr.h pkd.h
pkd.o: /usr/include/sys/resource.h ../mdl/null/mdl.h ewald.h grav.h walk.h
pkd.o: opentype.h tipsydefs.h
pst.o: /usr/include/math.h /usr/include/stdlib.h /usr/include/sgidefs.h
pst.o: /usr/include/stddef.h /usr/include/string.h /usr/include/malloc.h
pst.o: /usr/include/assert.h ../mdl/null/mdl.h /usr/include/stdio.h pst.h
pst.o: pkd.h /usr/include/sys/resource.h /usr/include/sys/time.h
pst.o: /usr/include/sys/types.h /usr/include/sys/bsd_types.h
pst.o: /usr/include/sys/select.h smoothfcn.h outtype.h smooth.h
smooth.o: /usr/include/stdio.h /usr/include/stdlib.h /usr/include/sgidefs.h
smooth.o: /usr/include/stddef.h /usr/include/malloc.h /usr/include/math.h
smooth.o: /usr/include/assert.h smooth.h pkd.h /usr/include/sys/resource.h
smooth.o: /usr/include/sys/time.h /usr/include/sys/types.h
smooth.o: /usr/include/sys/bsd_types.h /usr/include/sys/select.h
smooth.o: ../mdl/null/mdl.h smoothfcn.h
smoothfcn.o: /usr/include/math.h /usr/include/assert.h smoothfcn.h pkd.h
smoothfcn.o: /usr/include/sys/resource.h /usr/include/sys/time.h
smoothfcn.o: /usr/include/sys/types.h /usr/include/sgidefs.h
smoothfcn.o: /usr/include/sys/bsd_types.h /usr/include/sys/select.h
smoothfcn.o: ../mdl/null/mdl.h /usr/include/stdio.h
v_sqrt1.o: /usr/include/math.h
v_sqrt1.t3x.o: /usr/include/math.h
walk.o: /usr/include/math.h /usr/include/stdlib.h /usr/include/sgidefs.h
walk.o: /usr/include/stddef.h /usr/include/malloc.h /usr/include/assert.h
walk.o: walk.h pkd.h /usr/include/sys/resource.h /usr/include/sys/time.h
walk.o: /usr/include/sys/types.h /usr/include/sys/bsd_types.h
walk.o: /usr/include/sys/select.h ../mdl/null/mdl.h /usr/include/stdio.h
