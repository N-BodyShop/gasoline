#ifdef STARFORM
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "pkd.h"
#include "feedback.h"
#include "supernova.h"
#include "supernovaia.h"

#define max(A,B) ((A) > (B) ? (A) : (B))
#define min(A,B) ((A) < (B) ? (A) : (B))

/*
 * Supernova module for GASOLINE
 */

void snInitialize(SN *psn)
{
    SN sn;
    
    sn = (SN) malloc(sizeof(struct snContext));
    assert(sn != NULL);
    
    sn->bInitialized = 0;
    sn->dESN = 0.0;
    sn->dMSNrem = 0.0;
    sn->dMSNIImin = 0.0;
    sn->dMSNIImax = 0.0;
    sn->dMBmin = 0.0;
    sn->dMBmax = 0.0;
    sn->dFracBinSNIa = 0.0;
    sn->dMEjexp = 0;    
    sn->dMEjconst = 0; 
    sn->dMFeconst = 0.0;
    sn->dMFeexp = 0.0;
    sn->dMOxconst = 0.0;
    sn->dMOxexp = 0.0;
    sn->dSNIaMetals = 0.0;
    *psn = sn;
    }

void snFree(SN sn)
{
/*    free(sn->MSparam);
    free(sn->ppdva);*/
    
    free(sn);
    }

/*
   Subroutine to initialize parameters dealing with supernovae feedback
   */
void snInitConstants(SN sn) 
{
    MSPARAM tempMS;
    PDVAPARAM tempPDVA;
    
    sn->dMSNrem = 1.4;			/* mass of supernova remnant in solar masses */
    
    /*sn->dESN = 1e51;		*/	/* Energy of Supernovae in ergs. */
    sn->dMSNIImin = 8.0;		/* Mass above which stars
								   supernova in solar
								   masses */
    sn->dMSNIImax = 40.;		/* Mass below which stars
								   supernova in solar masses */
    sn->dMBmin = 3.0;			/* Minimum mass of binary that
								   can go SNIa */
    
    sn->dMBmax = 16.0;			/* Maximum mass of binary that
								   can go SNIa */
    sn->dFracBinSNIa = 0.16;	/* fraction of binary systems in
								   appropriate mass range that go SNIa =
								   0.16 (van den Bergh & McClure, ApJ
								   425, 205, 1994) */
    /* normalization constant and exponent in formulae for masses of
       ejected Fe and O16 as a function of stellar mass taken from
       Raiteri, Villata and Navarro, A&A 315, 105, 1996 */
    sn->dMEjexp = 1.056;    
    sn->dMEjconst = 0.7682; 
    sn->dMFeexp = 2.864;    
    sn->dMFeconst = 9.783e-5; 
    sn->dMOxexp = 3.721;    
    sn->dMOxconst = 1.2325e-4; 
/*    sn->dMFeexp = 1.864;    
    sn->dMFeconst = 2.802e-4;       The values from the article, but
    sn->dMOxexp = 2.721;            we want to integrate the functions
    sn->dMOxconst = 4.586e-4; */
    sn->dSNIaMetals = 0.76;  /* 0.63 Msol Fe + 0.13 Msol Ox (Theilemann 1986)*/
    
    MSInitialize(&tempMS);
    sn->MSparam = *tempMS;
    free(tempMS);
    
    PadovaInitialize(&tempPDVA);
    sn->ppdva = *tempPDVA;
    free(tempPDVA);

    sn->bInitialized = 1;
    }

void snCalcSNIIFeedback(SN sn, SFEvent sfEvent,
						double dTime, /* current time in years */
						double dDelta, /* length of timestep (years) */
						FBEffects *fbEffects)
{
    double dStarLtimeMin, dStarLtimeMax;
    double dMStarMin, dMStarMax;
    double dMStarMinII, dMStarMaxII;
    double dCumMMin, dCumMMax, dMtot;
    
    double dMSNTypeII, dESNTypeII, dNSNTypeII;
    double dCumNMinII, dCumNMaxII;
    double dMeanMStar;			/* Average mass of star going
								   supernovea for metalicity calculation */
    double  dDeltaMSNrem;


	/* stellar lifetimes corresponding to beginning and end of 
	   current timestep with respect to starbirth time in yrs */
    dStarLtimeMin = dTime - sfEvent.dTimeForm; 
    dStarLtimeMax = dStarLtimeMin + dDelta;

    dMStarMin = dSTMStarLtime(&sn->ppdva, dStarLtimeMax, sfEvent.dMetals); 
    dMStarMax = dSTMStarLtime(&sn->ppdva, dStarLtimeMin, sfEvent.dMetals); 

    dMtot = dMSCumMass(&sn->MSparam, 0.0); /* total mass in stars integrated over IMF */

	/* masses corresponding to these stellar lifetimes in solar masses */
    if (dMStarMin > sn->dMSNIImax || dMStarMax < sn->dMSNIImin) {
		fbEffects->dMassLoss = 0.0;
		fbEffects->dEnergy = 0.0;
		fbEffects->dMetals = 0.0;
                fbEffects->dMIron = 0.0;
                fbEffects->dMOxygen = 0.0;
		return;
		} else {
			dMStarMinII = max (sn->dMSNIImin, dMStarMin); 
			dMStarMaxII = min (sn->dMSNIImax, dMStarMax); 

			assert (dMStarMinII < dMStarMaxII && dMStarMinII >0 && dMStarMaxII > 0);

			/* cumulative mass of stars with mass greater than dMStarMinII and dMStarMaxII
			   in solar masses */
			dCumMMin = dMSCumMass (&sn->MSparam, dMStarMinII);
			dCumMMax = dMSCumMass (&sn->MSparam, dMStarMaxII);

			dMSNTypeII = dCumMMin - dCumMMax; /* total mass of stars that go SN II
												 in this timestep in solar masses */
			/* cumulative number of stars with mass greater than dMStarMinII and
			   dMstarMaxII in solar masses */
			dCumNMinII = dMSCumNumber (&sn->MSparam, dMStarMinII); 
			dCumNMaxII = dMSCumNumber (&sn->MSparam, dMStarMaxII);

			dMSNTypeII /= dMtot; /* convert to mass fraction of stars */
			dMSNTypeII *= sfEvent.dMass; /* convert to mass of stars that go SNIa */

			dNSNTypeII = dCumNMinII - dCumNMaxII;
			dNSNTypeII /= dMtot; /* convert to number per solar mass of stars */
			dNSNTypeII *= sfEvent.dMass; /* convert to absolute number of SNII */

			/* Average Star mass for metalicity calculation. */
			dMeanMStar = dMSNTypeII/dNSNTypeII;

			dESNTypeII = dNSNTypeII * sn->dESN;
			}
    
    /* decrement mass of star particle by mass of stars that go SN
       plus mass of SN remnants */
    dDeltaMSNrem = dNSNTypeII*sn->dMSNrem; /* mass in SN remnants */
    fbEffects->dMassLoss = dMSNTypeII - dDeltaMSNrem;

    /* SN specific energy rate to be re-distributed among neighbouring gas
       particles */
    fbEffects->dEnergy = dESNTypeII/(MSOLG*fbEffects->dMassLoss);   

	/* fraction of mass in metals to be re-distributed among neighbouring
	 * gas particles.  Formula 6-8 of  Raiteri, Villata and Navarro, A&A
	 * 315, 105, 1996  are used to calculate SNII yields
         * Integrate over power law from lowest mass progenitor to highest mass.
	 */
    fbEffects->dMIron = 
        sn->dMFeconst * (pow(dMStarMax, sn->dMFeexp) - pow(dMStarMin, sn->dMFeexp));
    fbEffects->dMOxygen = 
        sn->dMOxconst * (pow(dMStarMax, sn->dMOxexp) - pow(dMStarMin, sn->dMOxexp));
    fbEffects->dMetals = ( fbEffects->dMIron + fbEffects->dMOxygen )/ fbEffects->dMassLoss;
/*		/ (sn->dMEjconst * pow(dMeanMStar, sn->dMEjexp));*/
    
	}

void snCalcSNIaFeedback(SN sn, SFEvent sfEvent,
						double dTime, /* current time in years */
						double dDelta, /* length of timestep (years) */
						FBEffects *fbEffects)
{
    double dMStarMin, dMStarMax;
    double dMStarMinIa, dMStarMaxIa;
    double dMSNTypeIa, dESNTypeIa, dNSNTypeIa;
    struct inMSIMFSec mssn;    
    double dStarLtimeMin, dStarLtimeMax;
    double dMtot;
    double  dDeltaMSNrem;

    mssn.ms = sn->MSparam;		/* needed by SN Ia functions */
    mssn.sn = *sn;				/* needed by SN Ia functions */

	/* stellar lifetimes corresponding to beginning and end of 
	 * current timestep with respect to starbirth time in yrs */
    dStarLtimeMin = dTime - sfEvent.dTimeForm; 
    dStarLtimeMax = dStarLtimeMin + dDelta;

    dMStarMin = dSTMStarLtime(&sn->ppdva, dStarLtimeMax, sfEvent.dMetals); 
    dMStarMax = dSTMStarLtime(&sn->ppdva, dStarLtimeMin, sfEvent.dMetals); 

    dMtot = dMSCumMass(&sn->MSparam, 0.0); /* total mass in stars
											 integrated over IMF */

    if (dMStarMin > sn->dMBmin && dMStarMax < sn->dMBmax/2.) {

        dMStarMinIa = max (sn->dMBmin, dMStarMin); 
        dMStarMaxIa = min (sn->dMBmax/2., dMStarMax); 

        assert (dMStarMinIa < dMStarMaxIa && dMStarMinIa >0 && dMStarMaxIa > 0);
        
        /* mass of stars that go SNIa */
        dMSNTypeIa = dMSNIa (&mssn, dMStarMinIa, dMStarMaxIa); 
        
        /* number of stars that go SNIa */        
        dNSNTypeIa = dNSNIa (&mssn, dMStarMinIa, dMStarMaxIa); 
        dMSNTypeIa /= dMtot;	/* convert to mass fraction of stars */
        
        /* convert to mass of stars that go SNIa */ 
        dMSNTypeIa *= sfEvent.dMass; 

        dNSNTypeIa /= dMtot;	/* convert to number per solar mass of stars */
        dNSNTypeIa *= sfEvent.dMass; /* convert to absolute number of SNIa */

        dESNTypeIa = dNSNTypeIa * sn->dESN;
		} else {
			fbEffects->dMassLoss = 0.0;
			fbEffects->dEnergy = 0.0;    
			fbEffects->dMetals = 0.0; 
                        fbEffects->dMIron = 0.0;
                        fbEffects->dMOxygen = 0.0;
			return;
			}
    
    /* decrement mass of star particle by mass of stars that go SN
       plus mass of SN remnants */
    dDeltaMSNrem = dNSNTypeIa*sn->dMSNrem; /* mass in SN remnants */
    fbEffects->dMassLoss = dMSNTypeIa - dDeltaMSNrem;

    /* SN specific energy rate to be re-distributed among neighbouring gas
       particles */
    fbEffects->dEnergy = dESNTypeIa/(MSOLG*fbEffects->dMassLoss);   

    /* Following Raiteri 1996 who follows Thielemann's 1986 W7 model for
     * SNIa explosions, the same mass of iron and oxygen is released in
     * every single explosion.  A comparable amount of Silicon is ejected
     * to the amount of Oxygen that is ejected.
     */
    fbEffects->dMIron = dNSNTypeIa*0.63;
    fbEffects->dMOxygen = dNSNTypeIa*0.13;
	/* Fraction of mass in metals to be re-distributed among neighbouring
	 * gas particles: assumes fixed amount of metals per supernovea independent of
	 * mass. See Raiteri, Villata and Navarro, page 108.
	 */
    fbEffects->dMetals = dNSNTypeIa*sn->dSNIaMetals/fbEffects->dMassLoss; 
	}
#endif
