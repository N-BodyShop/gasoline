#ifdef STARFORM
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "pkd.h"
#include "millerscalo.h"
#include "supernova.h"
#include "startime.h"
#include "supernovaia.h"

#define max(A,B) ((A) > (B) ? (A) : (B))
#define min(A,B) ((A) < (B) ? (A) : (B))

/*
 * Star forming module for GASOLINE
 */

void snInitialize(SN *psn)
{
    SN sn;
    
    sn = (SN) malloc(sizeof(struct snContext));
    assert(sn != NULL);
    
    sn->bInitialized = 0;
    sn->dSecUnit = 0;
    sn->dGmPerCcUnit = 0;
    sn->dGmUnit = 0;
    sn->dErgUnit = 0;
    sn->dESN = 0.0;
    sn->dMSNrem = 0.0;
    sn->dMSNIImin = 0.0;
    sn->dMSNIImax = 0.0;
    sn->dMBmin = 0.0;
    sn->dMBmax = 0.0;
    sn->dFracBinSNIa = 0.0;
    sn->dMFeconst = 0.0;
    sn->dMFeexp = 0.0;
    sn->dMOxconst = 0.0;
    sn->dMOxexp = 0.0;
    *psn = sn;
    }

/*
     Subroutine to initialize parameters dealing with supernovae feedback
*/
void snInitConstants(SN sn) 
{
    assert(sn->dSecUnit != 0.0);
    assert(sn->dGmPerCcUnit != 0.0);
    assert(sn->dGmUnit != 0.0);
    
    sn->dMSNrem = 1.4*MSOLG/sn->dGmUnit; /* mass of supernova remnant in solar masses */
    
    sn->dESN = 1e51/sn->dErgUnit; /* Energy of Supernovae in ergs. */
    sn->dMSNIImin = 8.*MSOLG/sn->dGmUnit; /* Mass above which stars
                                            supernova in solar
                                            masses */
    sn->dMSNIImax = 40.*MSOLG/sn->dGmUnit; /* Mass below which stars
                                             supernova in solar masses */
    sn->dMBmin = 3.*MSOLG/sn->dGmUnit; /* Minimum mass of binary that
                                          can go SNIa */
    
    sn->dMBmax = 16.*MSOLG/sn->dGmUnit; /* Maximum mass of binary that
                                           can go SNIa */
    sn->dFracBinSNIa = 0.16; /* fraction of binary systems in
                                appropriate mass range that go SNIa =
                                0.16 (van den Bergh & McClure, ApJ
                                425, 205, 1994) */
    /* normalization constant and exponent in formulae for masses of
       ejected Fe and O16 as a function of stellar mass taken from
       Raiteri, Villata and Navarro, A&A 315, 105, 1996 */
    sn->dMFeexp = 1.864;    
    sn->dMFeconst = 2.802e-4; 
    sn->dMOxexp = 2.721;    
    sn->dMOxconst = 4.586e-4; 
/* normalize to system mass units */
    sn->dMFeconst *= pow (sn->dGmUnit/MSOLG, sn->dMFeexp-1);    
    sn->dMOxconst *= pow (sn->dGmUnit/MSOLG, sn->dMOxexp-1);    

    sn->bInitialized = 1;
    }

void pkdCalcSNEnergy(PKD pkd, SN sn, double dTime)
{
    int i;
    PARTICLE *p;
    int n = pkdLocal(pkd);
    
    for(i = 0; i < n; ++i) {
        p = &pkd->pStore[i];
	if(pkdIsStar(pkd, p))
	    snCalcSNEnergy(sn, pkd, p, dTime);
	else if(pkdIsGas(pkd, p))
	    p->fESNrate = 0; /* reset SN heating rate of gas to zero */
    }
}

void snCalcSNEnergy(SN sn, PKD pkd, PARTICLE *p, double dTime)
{
    MSPARAM MSparam;
    PDVAPARAM ppdva;
    struct inMSIMFSec mssn;    
    
    double dStarLtimeMin, dStarLtimeMax;
    double dMStarMin, dMStarMax;
    double dMStarMinII, dMStarMaxII;
    double dMStarMinIa, dMStarMaxIa;
    double dCumMMin, dCumMMax, dMtot;
    double dCumMMinII, dCumMMaxII;
    
    double dMSNTypeII, dESNTypeII, dNSNTypeII;
    double dMSNTypeIa, dESNTypeIa, dNSNTypeIa;
    double dCumNMin, dCumNMax;
    double dCumNMinII, dCumNMaxII;
    double dMSNtot, dNSNtot, dESNtot;

    double  dDeltaMSNrem;

    MSInitialize(&MSparam);
    PadovaInitialize(&ppdva);
    mssn.ms = *MSparam; /* needed by SN Ia functions */
    mssn.sn = *sn;      /* needed by SN Ia functions */

/* stellar lifetimes corresponding to beginning and end of 
   current timestep with respect to starbirth time in yrs */
    dStarLtimeMin = dTime - p->fTimeForm; 
    dStarLtimeMax = dStarLtimeMin + sn->dDeltaT;
    dStarLtimeMin *= sn->dSecUnit/SEC_YR;
    dStarLtimeMax *= sn->dSecUnit/SEC_YR;

    dMStarMin = dSTMStarLtime(ppdva, sn, dStarLtimeMax, p); 
    dMStarMax = dSTMStarLtime(ppdva, sn, dStarLtimeMin, p); 

    dMtot = dMSCumMass(MSparam, 0.0); /* total mass in stars integrated over IMF */

/* masses corresponding to these stellar lifetimes in solar masses */
    if (dMStarMin > sn->dMSNIImax || dMStarMax < sn->dMSNIImin) {
        dMSNTypeII = 0;
        dNSNTypeII = 0;
        dESNTypeII = 0;
    } else {
        dMStarMinII = max (sn->dMSNIImin, dMStarMin); 
        dMStarMaxII = min (sn->dMSNIImax, dMStarMax); 

        assert (dMStarMinII < dMStarMaxII && dMStarMinII >0 && dMStarMaxII > 0);

/* cumulative mass of stars with mass greater than dMStarMinII and dMStarMaxII
   in solar masses */
        dCumMMin = dMSCumMass (MSparam, dMStarMinII);
        dCumMMax = dMSCumMass (MSparam, dMStarMaxII);

        dMSNTypeII = dCumMMin - dCumMMax; /* total mass of stars that go SN II
                                             in this timestep in solar masses */
/* cumulative number of stars with mass greater than dMStarMinII and
   dMstarMaxII in solar masses */
        dCumNMinII = dMSCumNumber (MSparam, dMStarMinII); 
        dCumNMaxII = dMSCumNumber (MSparam, dMStarMaxII);

        dMSNTypeII /= dMtot; /* convert to mass fraction of stars */
        dMSNTypeII *= p->fMass; /* convert to mass of stars that go SNIa */

        dNSNTypeII = dCumNMinII - dCumNMaxII;
        dNSNTypeII /= dMtot; /* convert to number per solar mass of stars */
        dNSNTypeII *= p->fMass; /* convert to absolute number of SNII */

        dESNTypeII = dNSNTypeII * sn->dESN;
    }
    
    if (dMStarMin > sn->dMBmin && dMStarMax < sn->dMBmax/2.) {

        dMStarMinIa = max (sn->dMBmin, dMStarMin); 
        dMStarMaxIa = min (sn->dMBmax/2., dMStarMax); 

        assert (dMStarMinIa < dMStarMaxIa && dMStarMinIa >0 && dMStarMaxIa > 0);
        
        dMSNTypeIa = dMSNIa (&mssn, dMStarMinIa, dMStarMaxIa); /* mass of
                                                               stars that
                                                               go SNIa */
        dNSNTypeIa = dNSNIa (&mssn, dMStarMinIa, dMStarMaxIa); /* number of
                                                                stars that go
                                                                SNIa */
        dMSNTypeIa /= dMtot; /* convert to mass fraction of stars */
        dMSNTypeIa *= p->fMass; /* convert to mass of stars that go SNIa */

        dNSNTypeIa /= dMtot; /* convert to number per solar mass of stars */
        dNSNTypeIa *= p->fMass; /* convert to absolute number of SNIa */

        dESNTypeIa = dNSNTypeIa * sn->dESN;
    } else {
        dNSNTypeIa = 0;
        dMSNTypeIa = 0;
        dESNTypeIa = 0;
    }
    
    dESNtot = dESNTypeII + dESNTypeIa;
    dNSNtot = dNSNTypeII + dNSNTypeIa;
    dMSNtot = dMSNTypeII + dMSNTypeIa;


    dDeltaMSNrem = dNSNtot*sn->dMSNrem; /* mass in SN remnants */
    p->fMass += dDeltaMSNrem - dMSNtot; /* decrement mass of star
                                           particle by mass of stars
                                           that go SN plus mass of SN
                                           remnants */
    p->fMSN = dMSNtot - dDeltaMSNrem;   /* store mass lost by star
                                           particle for re-distribution
                                           as gas to neighbouring gas
                                           particles */
    p->fESNrate = dESNtot/p->fMass/sn->dDeltaT;    
/* SN specific energy rate to be re-distributed among neighbouring gas
   particles */

    p->fMetals += sn->dMFeconst * pow(dMSNtot, sn->dMFeexp); 
/* mass in Fe to be re-distributed among neighbouring gas particles */

    free(MSparam);
    free(ppdva);
}
#endif
