#
#       Charm makefile for MDL
#
CHARMC = $(CC)

default: mdl.o

mdl.o: mdl.ci mdl.C mdl.h mdlimpl.h mdl.def.h mdl.decl.h
	$(CHARMC) $(CFLAGS) -c mdl.C

mdl.def.h mdl.decl.h : mdl.ci
	$(CHARMC) mdl.ci
