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
  /*  sn->dMFeconst = 0.0;
      sn->dMFeexp = 0.0;*/
  sn->dMOxconst = 0.0;
  sn->dMOxexp = 0.0;
  sn->dMIronslope = 0.0;
  sn->dMIronintercept = 0.0;
  sn->dMOxygenslope = 0.0;
  sn->dMOxygenintercept = 0.0;
#ifdef MORE_METALS
  sn->dMCslope = 0.0;
  sn->dMCintercept = 0.0;
  sn->dMNslope = 0.0;
  sn->dMNintercept = 0.0;
  sn->dMNeslope = 0.0;
  sn->dMNeintercept = 0.0;
  sn->dMMgslope = 0.0;
  sn->dMMgintercept = 0.0;
  sn->dMSislope = 0.0;
  sn->dMSiintercept = 0.0;
#endif /*MORE_METALS*/
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
  AGBPARAM tempAGB;
    
  sn->dMSNrem = 1.4;			/* mass of supernova remnant in solar masses 
                                         * Also used for SNIa ejected mass */
    
  /*sn->dESN = 1e51;		*/	/* Energy of Supernovae in ergs. */
  sn->dMSNIImin = 8.0;		/* Mass above which stars
				   supernova in solar
				   masses */
  sn->dMSNIImax = 40.;		/* Mass below which stars
				   supernova in solar masses */
  sn->dMBmin = 0.8;			/* Minimum mass of secondary that
                                           can cause SNIa */
    
  sn->dMBmax = 16.0;			/* Maximum mass of binary that
                                           can go SNIa */
  sn->dFracBinSNIa = 0.05;/* fraction of binary systems in
			    appropriate mass range that go SNIa.
			    .05 is inline with chemical evolution
			    models of the Milky Way (Francois et al 2004) */

  /* normalization constant and exponent in formulae for masses of
     ejected Fe and O16 as a function of stellar mass taken from
     Raiteri, Villata and Navarro, A&A 315, 105, 1996 */
  sn->dMEjexp = 1.056;    
  sn->dMEjconst = 0.7682; 
  /*  sn->dMFeexp = 1.864;    
      sn->dMFeconst = 2.802e-4;*/
  sn->dMOxexp = 2.721;
  sn->dMOxconst = 4.586e-4; 
  sn->dMIronslope = 0.003;
  sn->dMIronintercept = 0.0;
  sn->dMOxygenslope = 0.210714;
  sn->dMOxygenintercept = -2.0;
#ifdef MORE_METALS
  sn->dMCslope = 0.01;
  sn->dMCintercept = 0.0;
  sn->dMNslope = 3.67e-3;
  sn->dMNintercept = -6.7e-3;
  sn->dMNeslope = 0.08;
  sn->dMNeintercept = -2.0;
  sn->dMMgslope = 0.01;
  sn->dMMgintercept = -0.1;
  sn->dMSislope = 0.01;
  sn->dMSiintercept = -0.1;
#endif /*MORE_METALS*/
  sn->dSNIaMetals = 1.4;  /* Entire Chandrashekar mass goes to metals
			     (lots of Ni) */
    
  MSInitialize(&tempMS);
  sn->MSparam = *tempMS;
  free(tempMS);
    
  PadovaInitialize(&tempPDVA);
  sn->ppdva = *tempPDVA;
  free(tempPDVA);

  AGBInitialize(&tempAGB);
  sn->AGBtabs = *tempAGB;
  free(tempAGB);

  sn->bInitialized = 1;
}

#define YIELD_LIN(min, max, m, b) (0.5*m*(max*max - min*min) + b*(max-min) )
#define MET_YIELD_LIN(min, max, m, b,met) met*(0.5*m*(max*max - min*min) + b*(max-min) )

void snCalcSNIIFeedback(SN sn, SFEvent sfEvent,
			double dTime, /* current time in years */
			double dDelta, /* length of timestep (years) */
			FBEffects *fbEffects)
{
  double dStarLtimeMin, dStarLtimeMax;
  double dMStarMin, dMStarMax;
  double dMStarMinII, dMStarMaxII;
  double dCumMMin, dCumMMax, dMtot, dMej, dMmet;
    
  double dMSNTypeII, dESNTypeII, dNSNTypeII;
  double dCumNMinII, dCumNMaxII;
  double dMeanMStar;/* Average mass of star going
		       supernovea for metalicity calculation */
  double  dDeltaMSNrem;
  double dLuminosity, dOneMtot;

  /* stellar lifetimes corresponding to beginning and end of 
     current timestep with respect to starbirth time in yrs */
  dStarLtimeMin = dTime - sfEvent.dTimeForm; 
  dStarLtimeMax = dStarLtimeMin + dDelta;

  /* masses corresponding to these stellar lifetimes in solar masses */
  dMStarMin = dSTMStarLtime(&sn->ppdva, dStarLtimeMax, sfEvent.dMetals); 
  dMStarMax = dSTMStarLtime(&sn->ppdva, dStarLtimeMin, sfEvent.dMetals); 

  dMtot = dMSCumMass(&sn->MSparam, 0.0); /* total mass in stars integrated over IMF */

  if (sn->iNSNIIQuantum && dMStarMin > sn->dMSNIImax && dStarLtimeMin){
    /* blow wind before SN */

    fbEffects->dMetals = 0.0; /* zero enrichment from winds */
    fbEffects->dMIron = 0.0;
    fbEffects->dMOxygen = 0.0;
#ifdef MORE_METALS
    fbEffects->dMC = 0.0;
    fbEffects->dMN = 0.0;
    fbEffects->dMNe = 0.0;
    fbEffects->dMMg = 0.0;
    fbEffects->dMSi = 0.0;
#endif

    dNSNTypeII = dMSCumNumber (&sn->MSparam, sn->dMSNIImin);
    dNSNTypeII *= sfEvent.dMass/dMtot; /* normalize to star particle mass */
    fbEffects->dMassLoss = dNSNTypeII * dDelta / 1e6; /* 1 M_sol / Myr / SN */
    fbEffects->dEnergy = dNSNTypeII * 3e42 * dDelta/ /* 1e35 erg/s /SN */
      (MSOLG*fbEffects->dMassLoss); 
    
  } 
#ifndef JAMESFB
  else if (sn->dRadPresFrac && dMStarMin > sn->dMSNIImax && dStarLtimeMin){
    /* use radiation pressure before SN */

    fbEffects->dMetals = 0.0; /* zero enrichment from photons */
    fbEffects->dMIron = 0.0;
    fbEffects->dMOxygen = 0.0;
#ifdef MORE_METALS
    fbEffects->dMC = 0.0;
    fbEffects->dMN = 0.0;
    fbEffects->dMNe = 0.0;
    fbEffects->dMMg = 0.0;
    fbEffects->dMSi = 0.0;
#endif

    dLuminosity = imfCumLuminosity(&sn->MSparam, 1.0);
    dOneMtot = dMSCumMass(&sn->MSparam, 1.0); /* total mass in stars
						 integrated over IMF
						 from 1 M_sun up */
    dLuminosity *= sfEvent.dMass/dOneMtot; /* normalize to star particle
					   mass */
    fbEffects->dMassLoss = 1.0; /* in M_sol: minimal mass loss */
#define SECONDSPERYEAR   31557600.
    fbEffects->dEnergy = sn->dRadPresFrac * dLuminosity * 2e33 *
	SECONDSPERYEAR* dDelta/(MSOLG*fbEffects->dMassLoss); 
    
      } 
#endif /* JAMESFB */
  else if (dMStarMin > sn->dMSNIImax  || dMStarMax < sn->dMSNIImin) {
    /* do nothing */
    fbEffects->dMassLoss = 0.0;
    fbEffects->dEnergy = 0.0;
    fbEffects->dMetals = 0.0;
    fbEffects->dMIron = 0.0;
    fbEffects->dMOxygen = 0.0;
#ifdef MORE_METALS
    fbEffects->dMC = 0.0;
    fbEffects->dMN = 0.0;
    fbEffects->dMNe = 0.0;
    fbEffects->dMMg = 0.0;
    fbEffects->dMSi = 0.0;
#endif
    return;
  } else {
    /* supernova */
    dMStarMinII = max (sn->dMSNIImin, dMStarMin); 
    dMStarMaxII = min (sn->dMSNIImax, dMStarMax);

    assert (dMStarMinII < dMStarMaxII && 
	    dMStarMinII >0.0 && dMStarMaxII > 0.0);

    /* cumulative mass of stars with mass greater than dMStarMinII 
       and dMStarMaxII in solar masses */
    dCumMMin = dMSCumMass (&sn->MSparam, dMStarMinII);
    dCumMMax = dMSCumMass (&sn->MSparam, dMStarMaxII);

    if(dCumMMax > dCumMMin || dCumMMax < 0) dMSNTypeII = dCumMMin;
    else dMSNTypeII = dCumMMin - dCumMMax; /* mass of stars that go SN II
					 in this timestep FAKE MASS 
					 (normalized to IMF stuff, 
					 REAL MASS is below).*/

    dMSNTypeII *= sfEvent.dMass/dMtot; /* REAL MASS of stars that go SNII */

    /* cumulative number of stars with mass greater than dMStarMinII and
       less than dMstarMaxII in solar masses */
    dCumNMinII = dMSCumNumber (&sn->MSparam, dMStarMinII); 
    dCumNMaxII = dMSCumNumber (&sn->MSparam, dMStarMaxII);

    if(dCumNMaxII > dCumNMinII || dCumNMaxII < 0) dNSNTypeII = dCumNMinII;
    else dNSNTypeII = dCumNMinII - dCumNMaxII;
    dNSNTypeII *= sfEvent.dMass/dMtot; /* number of SNII in star particle */

    /* Average Star mass for metalicity calculation. */
    if (dNSNTypeII > 0) dMeanMStar = dMSNTypeII/dNSNTypeII;
    else dMeanMStar =0;

    dESNTypeII = dNSNTypeII * sn->dESN;
    
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
     * Integrate over power law from lowest mass progenitor to highest
     * mass.
     * See Gibson (2002) IAUS conf. proc. Fig 1 for Oxygen + Fe
     * behavior in progenitors up to 100 M_sun 
     * -O not M^3 past 30 M_sun!
     * -Fe doesn't go much above 0.3 M_sun
     */
    /*    fbEffects->dMIron = sn->dMFeconst * pow(dMeanMStar, sn->dMFeexp)
      / (sn->dMEjconst * pow(dMeanMStar, sn->dMEjexp));
    fbEffects->dMOxygen = sn->dMOxconst*pow(dMeanMStar, sn->dMOxexp)
      / (sn->dMEjconst * pow(dMeanMStar, sn->dMEjexp));
        printf("old Oxygen: %g\n",fbEffects->dMOxygen);
	  printf("dMeanMStar: %g\n",dMeanMStar);*/
    dMej = (sn->dMEjconst) / (sn->dMEjexp + 1 ) * 
	(pow(dMStarMaxII,sn->dMEjexp + 1) - pow(dMStarMinII,sn->dMEjexp + 1));
    fbEffects->dMIron = YIELD_LIN(dMStarMinII,dMStarMaxII,sn->dMIronslope,
				  sn->dMIronintercept) / dMej;
    dMmet = max(0,YIELD_LIN(dMStarMinII,dMStarMaxII,sn->dMOxygenslope,
			    sn->dMOxygenintercept) );
    fbEffects->dMOxygen = dMmet / dMej;
    /*	/ (sn->dMEjconst * pow(dMeanMStar, sn->dMEjexp));*/
    /*    printf("new Oxygen: %g\n",fbEffects->dMOxygen);
	  printf("dMStarMinII: %g; dMStarMaxII: %g\n",dMStarMinII,dMStarMaxII);*/

#ifdef MORE_METALS
    /* Approximate C, N, Ne, Mg, and Si yields using linear fits to 
     * a combination of Woosley & Weaver (1995) and Kobayashi, C et al (2006)
     */
    dMmet = max(0, YIELD_LIN(dMStarMinII,dMStarMaxII,sn->dMCslope,
			     sn->dMCintercept));
    fbEffects->dMC = dMmet / dMej;

    dMmet = max(0, MET_YIELD_LIN(dMStarMinII,dMStarMaxII,sn->dMNslope,
				 sn->dMNintercept,sfEvent.dMetals/0.02) );
    fbEffects->dMN = dMmet / dMej;

    dMmet = max(0, YIELD_LIN(dMStarMinII,dMStarMaxII,sn->dMNeslope,
			     sn->dMNeintercept) );
    fbEffects->dMNe = dMmet / dMej;

    dMmet = max(0, YIELD_LIN(dMStarMinII,dMStarMaxII,sn->dMMgslope,
			     sn->dMMgintercept) );
    fbEffects->dMMg = dMmet / dMej;

    dMmet = max(0, YIELD_LIN(dMStarMinII,dMStarMaxII,sn->dMSislope,
			     sn->dMSiintercept) );
    fbEffects->dMSi = dMmet / dMej;

    fbEffects->dMetals = ( fbEffects->dMIron + fbEffects->dMOxygen +
			   fbEffects->dMC + fbEffects->dMN + fbEffects->dMNe +
			   fbEffects->dMMg + fbEffects->dMSi);
#else
    /* Use ratio of Fe to total iron group and O to total non-iron
       group derived from Asplund et al 2009 */
    fbEffects->dMetals = ( 1.06*fbEffects->dMIron + 2.09*fbEffects->dMOxygen );
#endif /* MORE_METALS */
  }
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

  dMtot = dMSCumMass(&sn->MSparam, 0.0); /* total mass in stars integrated over IMF */

  if (dMStarMin < sn->dMBmax/2.) {

    dMStarMinIa = dMStarMin; 
    dMStarMaxIa = min (sn->dMBmax/2., dMStarMax); 

    assert (dMStarMinIa < dMStarMaxIa && dMStarMinIa >0.0 && dMStarMaxIa > 0.0);
        
    /* mass of stars that go SNIa based on normalized IMF */
    /*        dMSNTypeIa = dMSNIa (&mssn, dMStarMinIa, dMStarMaxIa); 
	      dMSNTypeIa /= dMtot;	/* convert to mass fraction of stars /
	      /* convert to mass of stars that go SNIa / 
	      dMSNTypeIa *= sfEvent.dMass; */

    /* number of stars that go SNIa */        
    dNSNTypeIa = dNSNIa (&mssn, dMStarMinIa, dMStarMaxIa); 
    dNSNTypeIa /= dMtot;	/* convert to number per solar mass of stars */
    dNSNTypeIa *= sfEvent.dMass; /* convert to absolute number of SNIa */

    dESNTypeIa = dNSNTypeIa * sn->dESN;
    /* decrement mass of star particle by mass of stars that go SN
       and SN Ia have no remnants.   
       The MSNrem is the same Chandrasekhar mass that explodes as a SNIa.*/
    fbEffects->dMassLoss = dNSNTypeIa*sn->dMSNrem;
    
    /* SN specific energy rate to be re-distributed among neighbouring gas
       particles */
    fbEffects->dEnergy = dESNTypeIa/(MSOLG*fbEffects->dMassLoss);   
    
    /* Following Raiteri 1996 who follows Thielemann's 1986 W7 model for
     * SNIa explosions, the same mass of iron and oxygen is released in
     * every single explosion.  A comparable amount of Silicon is ejected
     * to the amount of Oxygen that is ejected.
     */
    fbEffects->dMIron = dNSNTypeIa*0.74/fbEffects->dMassLoss;
    fbEffects->dMOxygen = dNSNTypeIa*0.143/fbEffects->dMassLoss;
    /* Fraction of mass in metals to be re-distributed among neighbouring
     * gas particles: assumes fixed amount of metals per supernovea independent of
     * mass. See Raiteri, Villata and Navarro, page 108.
     */
#ifdef MORE_METALS
    /* should use W7 column in Table 1 from Nomoto et al (1997) 
     * http://adsabs.harvard.edu/abs/1997NuPhA.621..467N
     */
    fbEffects->dMC = dNSNTypeIa*0.0483/fbEffects->dMassLoss;
    fbEffects->dMN = dNSNTypeIa*1.16e-6/fbEffects->dMassLoss;
    fbEffects->dMNe = dNSNTypeIa*0.0045/fbEffects->dMassLoss;
    fbEffects->dMMg = dNSNTypeIa*0.009/fbEffects->dMassLoss;
    fbEffects->dMSi = dNSNTypeIa*0.15/fbEffects->dMassLoss;
#endif
  fbEffects->dMetals = dNSNTypeIa*sn->dSNIaMetals/fbEffects->dMassLoss; 
  } else {
    fbEffects->dMassLoss = 0.0;
    fbEffects->dEnergy = 0.0;    
    fbEffects->dMetals = 0.0; 
    fbEffects->dMIron = 0.0;
    fbEffects->dMOxygen = 0.0;
#ifdef MORE_METALS
    fbEffects->dMC = 0.0;
    fbEffects->dMN = 0.0;
    fbEffects->dMNe = 0.0;
    fbEffects->dMMg = 0.0;
    fbEffects->dMSi = 0.0;
#endif
    return;
  }
    
}
#endif
