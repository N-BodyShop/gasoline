CC = gcc
testcool: testcool.c 
	${CC} -g -fopenmp -I../../../mdl/null -I../../ -DGASOLINE -DCOOLING_METAL ../../cooling_metal.c ../../stiff.c ../../log.c ../../param.c testcool.c -o testcool -lm 

