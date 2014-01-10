#
# Custom section
# Set PVM_ARCH to your architecture type (SUN4, HP9K, RS6K, SGI, etc)
# if PVM_ARCH = BSD386 then set ARCHLIB = -lrpc
# if PVM_ARCH = SGI    then set ARCHLIB = -lsun
# if PVM_ARCH = I860   then set ARCHLIB = -lrpc -lsocket
# if PVM_ARCH = IPSC2  then set ARCHLIB = -lrpc -lsocket
# otherwise leave ARCHLIB blank
#
# PVM_ARCH and ARCHLIB are set for you if you use 'aimk'.
#
# PVM_ARCH	=	UNKNOWN
# ARCHLIB =
# END of custom section - leave this line here
#
PVMDIR	=	$(PVM_ROOT)
PVMLIB	=	$(PVMDIR)/lib/$(PVM_ARCH)/libpvm3.a
SDIR	=	..
BDIR	=	$(HOME)/pvm3/bin
XDIR	=	$(BDIR)/$(PVM_ARCH)

CFLAGS	=	-32 -O2 -I$(PVMDIR)/include
LIBS	=	$(PVMLIB) $(ARCHLIB) -lm

default:	mdl.o

mdl.o: $(SDIR)/mdl.c $(SDIR)/mdl.h
	$(CC) $(CFLAGS) -c $(SDIR)/mdl.c
