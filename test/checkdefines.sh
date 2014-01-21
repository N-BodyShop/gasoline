#!/bin/bash

#These are defines used for internal purposes that don't need to be user-enabled.
EXCLUDES='^__i386__$|^__AGGS_H$|HINCLUDED$|^COOLING_COSMO$|^COOLING_BATE$|^COOLING_METAL$|^[/][*]|^TRUE$|^FALSE$|^\#ifdef$|^Y_EMIN$|^AMPI$|^NOCOOLING$|^COOLING_MOLECULARH$|^_CRAYMPP$|^CRAY_T3D$|^CRAY_XT3$|^BOOLEAN$|^MAX$|^COOLING_BOLEY$|^COOLING_DISK$|^COOLING_PLANET$|^COOLING_POLY$|^__crayx1$|^__DATE__$|^__TIME__$|^DBL_MAX$|^FLT_MAX$|^lint$|^__FAST_MATH__$|^__LINALG_H$|^_REENTRANT$|^TWO_PI$|^SINGLE$|^MAXPATHLEN$|^NEWTIME$|^SGN$|^_LARGE_FILES$|^MAXLISTLEN$|^N_DIM$|^SUPPRESSMASSCHECKREPORT$|^DEBUGFORCE$|^NEWGLASS$|^LARGEFBALL$|^NEWINTEG$|^TESTRATE$'
#Get the defines that are in the Makefile
MAKEDEFS=`grep "CODE_DEF +=" ../Makefile | sed -e 's/=/ /g' | awk '{print $3;}' | sed -e 's/-D//' | sort | uniq`
#Get the defines that are in the code
DEFSFOUND=`grep 'ifdef\|ifndef' ../*.[ch] | awk '{print $2;}' |  sort | uniq | grep -v -E $EXCLUDES`
PASS=true
#Check for defines missing in the makefile but found in the C
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

#Vice Versa
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
