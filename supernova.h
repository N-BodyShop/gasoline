#include "pkd.h"

typedef struct snContext 
{
    PKD pkd;                    /* pkdContext needed for pkdDeleteParticle 
                                   in DeleteGas (smoothfcn.c) */
    double dMinMassFrac;        /* minimum fraction of average mass of neighbour 
                                   particles required for particle to avoid deletion */
    double dDeltaT;		/* timestep in system units */
    double dGmUnit;		/* system mass in grams */
    double dGmPerCcUnit;	/* system density in gm/cc */
    double dSecUnit;		/* system time in seconds */
    double dErgUnit;		/* system energy in ergs */

    double dESN;		/* SN energy released per SN */
    double dMSNrem;		/* mass of supernovae remnant */

    double dMSNIImin;           /* minimum mass of star that goes SNII */
    double dMSNIImax;           /* maximum mass of star that goes SNII */

    double dMBmin;              /* minimum mass of binary that goes SNIa */
    double dMBmax;              /* maximum mass of binary that goes SNIa */
    
    double dFracBinSNIa;        /* fraction of binaries that go SNIa */

    double dMFeconst;           /* normalization constant for formula
                                   for mass of ejected Fe as a
                                   function of stellar mass */
    double dMFeexp;             /* power of stellar mass in formula
                                   for mass of ejected Fe as a
                                   function of stellar mass */
    double dMOxconst;           /* same for oxygen */
    double dMOxexp;
    
    int bInitialized;
} * SN;

void snInitialize(SN *psn);

void snInitConstants(SN sn);

void pkdCalcSNEnergy(PKD pkd, SN sn, double dTime);
void snCalcSNEnergy(SN sn, PKD pkd, PARTICLE *p, double dTime);

        /* seconds per year = 3600*24*365.25 */
#define SEC_YR 3.15576e07;
	/* mass of hydrogen atom in grams */
#define MHYDR 1.67e-24
	/* solar mass in grams */
#define MSOLG 1.99e33

#ifndef STARFORM_HINCLUDED
#define STARFORM_HINCLUDED
#endif
