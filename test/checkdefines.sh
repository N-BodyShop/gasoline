#!/bin/bash

MAKEDEFS=`grep "CODE_DEF +=" ../Makefile | sed -e 's/=/ /g' | awk '{print $3;}' | sed -e 's/-D//' | sort | uniq`
DEFSFOUND=`grep 'ifdef\|ifndef' ../*.[ch] | awk '{print $2;}' |  sort | uniq | grep -v '^__i386__$' | grep -v '^__AGGS_H$' | grep -v 'HINCLUDED$' | grep -v '^COOLING_COSMO$' | grep -v '^COOLING_BATE$' | grep -v '^COOLING_METAL$'  | grep -v "^[/][*]" | grep -v "^TRUE$" | grep -v "^FALSE$" | grep -v "^\#ifdef$" | grep -v "^Y_EMIN$" | grep -v "^AMPI$" | grep -v "^NOCOOLING$" | grep -v "^COOLING_MOLECULARH$" | grep -v "^_CRAYMPP$" | grep -v "^CRAY_T3D$" | grep -v "^CRAY_XT3$" | grep -v "^BOOLEAN$" | grep -v "^MAX$" | grep -v "^COOLING_BOLEY$" | grep -v "^COOLING_DISK$" | grep -v "^COOLING_PLANET$" | grep -v "^COOLING_POLY$" | grep -v "^__crayx1$" | grep -v "^__DATE__$" | grep -v "^__TIME__$" | grep -v "^DBL_MAX$" | grep -v "^FLT_MAX$" | grep -v "^lint$" | grep -v "^__FAST_MATH__$" | grep -v "^__LINALG_H$" | grep -v "^_REENTRANT$" | grep -v "^TWO_PI$" | grep -v "^SINGLE$" | grep -v "^MAXPATHLEN$" | grep -v "^NEWTIME$" | grep -v "^SGN$" | grep -v "^_LARGE_FILES$" | grep -v "^MAXLISTLEN$" | grep -v "^N_DIM$" | grep -v "^SUPPRESSMASSCHECKREPORT$" | grep -v "^DEBUGFORCE$" | grep -v "^NEWGLASS$" | grep -v "^LARGEFBALL$" | grep -v "^NEWINTEG$" | grep -v "^TESTRATE$"`
PASS=true
for DEF in $DEFSFOUND
do
    echo $MAKEDEFS | grep $DEF> /dev/null
    if [ $?  -ne 0 ]
    then
        PASS=false
        echo $DEF "Missing in Makefile!"
        for FNAME in `grep -l $DEF ../*.[ch]`
        do
            echo $DEF "is in" $FNAME
        done
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
