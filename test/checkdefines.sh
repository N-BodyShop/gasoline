#!/bin/bash

MAKEDEFS=`grep "CODE_DEF +=" ../Makefile | awk '{print $3;}' | sed -e 's/-D//'`
DEFSFOUND=`grep -h 'ifdef\|ifndef' ../*[ch] | awk '{print $2;}' |  sort | uniq | grep -v '^__i386__$' | grep -v '^__AGGS_H$' | grep -v 'HINCLUDED$' | grep -v '^COOLING_COSMO$' | grep -v '^COOLING_BATE$' | grep -v '^COOLING_METAL$'`
PASS=true
for DEF in $DEFSFOUND
do
	echo $MAKEDEFS | grep $DEF> /dev/null
	if [ $?  -ne 0 ]
	then
		PASS=false
		echo $DEF "Missing in Makefile!"
	fi
done
for DEF in $MAKEDEFS
do
	echo $DEFSFOUND | grep $DEF> /dev/null
	if [ $?  -ne 0 ]
	then
		PASS=false
		echo $DEF "Found in Makefile, but isn't in the code!"
	fi
done
if $PASS
then
	echo "Makefile defines are good!"
fi
