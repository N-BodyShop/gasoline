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
PVM_LIBMDL	=	v_sqrt1.o $(PVM_MDL)/$(PVM_ARCH)/mdl-lru.o $(PVMLIB) $(ARCHLIB) /local/lib/libmalloc.a -lm

#
#       KSR1 defines
#
KSR_MDL = ../mdl/ksr

KSR_CFLAGS = -O2 -w2 -I$(KSR_MDL)
KSR_LD_FLAGS = -para
KSR_LIBMDL = erf.o v_sqrt1.o $(KSR_MDL)/mdl.o -lm -lrpc

#
#       SP1/2 defines
#
MDL = ../mdl/spx
SPX_CFLAGS = -g -I$(MDL)
SPX_LD_FLAGS = 
SPX_LIBMDL = v_sqrt1.o $(MDL)/mdl.o -lm

#
#       T3D MPP defines
#
T3D_MDL = $(HOME)/ptools/mdl/mpp
V_SQRT = $(HOME)/ptools/v_sqrt/lib
V_SQRT1 = $(HOME)/ptools/v_sqrt1/lib
RPC = $(HOME)/rpc
RPCOBJS = $(RPC)/xdr.o $(RPC)/xdr_mem.o $(RPC)/xdr_rec.o \
          $(RPC)/xdr_reference.o $(RPC)/xdr_stdio.o $(RPC)/xdr_float.o

T3D_CFLAGS = -O3 -DCRAY_T3D -I$(T3D_MDL) -I$(V_SQRT)
T3D_LIBMDL = -O3 -L $(V_SQRT) -L $(V_SQRT1) hyperlib.o $(T3D_MDL)/mdl.o \
	-lv_sqrtc -lv_sqrtc1 -lm
T3D_LD_FLAGS = 

OBJS	= 	main.o master.o param.o outtype.o pkd.o pst.o grav.o \
		ewald.o walk.o eccanom.o hypanom.o
EXTRA_OBJ = 	erf.o v_sqrt1.o hyperlib.o

default:	
	@echo "Please tell me what architecture to make."
	@echo "Choices are pvm, ksr, spx, and t3d."

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
	make pkdgrav CC=mpcc "CFLAGS=$(SPX_CFLAGS)" "LD_FLAGS=$(SPX_LD_FLAGS)" "MDL=$(SPX_MDL)" "LIBMDL=$(SPX_LIBMDL)"

t3d:
	cd $(T3D_MDL); make
	make pkdgrav "CFLAGS=$(T3D_CFLAGS)" "LD_FLAGS=$(T3D_LD_FLAGS)" "MDL=$(T3D_MDL)" "LIBMDL=$(T3D_LIBMDL)"

pkdgrav: $(OBJS) $(EXTRA_OBJ)
	$(CC) $(CFLAGS) $(LD_FLAGS) -o pkdgrav $(OBJS) $(LIBMDL)






