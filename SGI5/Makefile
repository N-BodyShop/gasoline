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

MDL		=	$(HOME)/bulk/ptools/mdl/mdl-0.2/pvm
CFLAGS	=	-32 -O2 -fullwarn -I$(PVMDIR)/include -I$(MDL)
LIBS	=	$(MDL)/$(PVM_ARCH)/mdl.o $(PVMLIB) $(ARCHLIB) -lm

default:	pkdgrav

$(XDIR):
	- mkdir $(BDIR)
	- mkdir $(XDIR)

clean:
	rm -f *.o


pkdgrav: main.o master.o param.o outtype.o pkd.o pst.o smooth.o grav.o \
	ewald.o walk.o
	$(CC) $(CFLAGS) -o pkdgrav main.o master.o param.o outtype.o pkd.o \
	 pst.o smooth.o grav.o ewald.o walk.o $(LIBS)
	mv -f pkdgrav $(XDIR)

main.o: $(SDIR)/main.c
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

smooth.o: $(SDIR)/smooth.c
	$(CC) $(CFLAGS) -c $(SDIR)/smooth.c

grav.o: $(SDIR)/grav.c
	$(CC) $(CFLAGS) -c $(SDIR)/grav.c

ewald.o: $(SDIR)/ewald.c
	$(CC) $(CFLAGS) -c $(SDIR)/ewald.c

walk.o: $(SDIR)/walk.c
	$(CC) $(CFLAGS) -c $(SDIR)/walk.c













