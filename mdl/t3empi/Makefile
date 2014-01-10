#
#       SPX makefile for MDL
#
CC_FLAGS = -O3

default: mdl.o

mdl.o: mdl.c mdl.h
	cc -I/opt/ctl/craylibs/craylibs/include $(CC_FLAGS) -c mdl.c

#cc -I/usr/local/mpp/lib/include/mpp $(CC_FLAGS) -c mdl.c
