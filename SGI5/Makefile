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

MDL		=	$(SDIR)/../mdl/pvm
CC=gcc
CFLAGS	=	-g -O2 -I$(PVMDIR)/include -I$(MDL)
LIB	=	$(MDL)/$(PVM_ARCH)/mdl.o $(PVMLIB) $(ARCHLIB) -lm

default:	mdl pkdgrav

$(XDIR):
	- mkdir $(BDIR)
	- mkdir $(XDIR)

clean:
	rm -f *.o

mdl:
	cd $(MDL); aimk

pkdgrav: main.o master.o param.o outtype.o pkd.o pst.o grav.o \
	ewald.o walk.o v_sqrt1.o eccanom.o hypanom.o fdl.o htable.o smooth.o
	$(CC) $(CFLAGS) -o pkdgrav main.o master.o param.o outtype.o pkd.o \
	 pst.o grav.o ewald.o walk.o v_sqrt1.o eccanom.o hypanom.o fdl.o \
	htable.o smooth.o $(LIB)
	mv -f pkdgrav $(XDIR)

main.o: $(SDIR)/main.c $(SDIR)/master.h
	$(CC) $(CFLAGS) -c $(SDIR)/main.c

master.o: $(SDIR)/master.c
	$(CC) $(CFLAGS)	-c $(SDIR)/master.c

param.o: $(SDIR)/param.c
	$(CC) $(CFLAGS) -c $(SDIR)/param.c

outtype.o: $(SDIR)/outtype.c
	$(CC) $(CFLAGS) -c $(SDIR)/outtype.c

pkd.o: $(SDIR)/pkd.c
	$(CC) $(CFLAGS) -c $(SDIR)/pkd.c

pst.o: $(SDIR)/pst.c
	$(CC) $(CFLAGS) -c $(SDIR)/pst.c

grav.o: $(SDIR)/grav.c
	$(CC) $(CFLAGS) -c $(SDIR)/grav.c

ewald.o: $(SDIR)/ewald.c
	$(CC) $(CFLAGS) -c $(SDIR)/ewald.c

walk.o: $(SDIR)/walk.c
	$(CC) $(CFLAGS) -c $(SDIR)/walk.c

v_sqrt1.o: $(SDIR)/v_sqrt1.c
	$(CC) $(CFLAGS) -c $(SDIR)/v_sqrt1.c

eccanom.o: $(SDIR)/eccanom.c
	$(CC) $(CFLAGS) -c $(SDIR)/eccanom.c

hypanom.o: $(SDIR)/hypanom.c
	$(CC) $(CFLAGS) -c $(SDIR)/hypanom.c

fdl.o: $(SDIR)/fdl.c
	$(CC) $(CFLAGS) -c $(SDIR)/fdl.c

htable.o: $(SDIR)/htable.c
	$(CC) $(CFLAGS) -c $(SDIR)/htable.c

smooth.o: $(SDIR)/smooth.c
	$(CC) $(CFLAGS) -c $(SDIR)/smooth.c










