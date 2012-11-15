#ifdef STARFORM
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "pkd.h"
#include "supernova.h"

void AGBInitialize(AGBPARAM *pagb) 
{
    AGBPARAM agb;
    struct AGBContext initagb = {
	{
	    { -3.61000e-05, -4.09000e-05, 0.0027925, 0.0027922, 0.0031905, 
	      0.0034096, 0.004608, 0.0049273, 0.00469660, 0.001313, 0.0006325,  
	      0.000825,  0.0006888},
	    { -7.96100e-05,  0.00068539, 0.00171888, 0.00171888, 0.00219961,   
	      0.00350988, 0.0049408, 0.0059123, 0.0073926, 0.001108, 0.0003818, 
	      0.000509, 0.000461},
	    { -0.0001346, -0.0001546, 0.001297, 0.0012971, 0.0025483, 0.002929, 
	      0.0044115, 0.0055536, 0.0068345, 0.000958, 9.81400e-05, 0.0001758,  
	      0.000106},
	    { -0.0003841, -7.50000e-05, 3.36000e-05, -0.0004449, 0.0008553, 
	      0.0013303, 0.0024097, 0.0040841, 0.0051768, 0.0046898, -0.0006,
	      -0.000616,-0.000747},
	    { -0.0009011, -0.000234, -0.000677, -0.0013618, -0.001561, -0.001665, 
	      -0.000938, 0.000761, 0.002048, 0.004701, 0.004848, -0.001604, 
	      -0.00167300}
	},

	{
	{2.64000e-05, 3.04000e-05, 3.64000e-05, 3.71000e-05, 4.24000e-05,  
	 4.59000e-05, 8.11000e-05, 8.82000e-05, 8.02000e-05, 0.00472, 0.00657,
	 0.00869000,   0.00794000},
	{0.000105, 0.000124, 0.00016, 0.00016, 0.000198, 0.000215, 0.000337, 
	 0.000418, 0.000444, 0.00447, 0.00639, 0.00858, 0.00939},
	{0.000198, 0.000241, 0.000306, 0.000308, 0.000356000,  0.000394000,  
	 0.000576, 0.000733, 0.000784, 0.00443, 0.00643, 0.00881, 0.0111},
	{0.000485, 0.000773, 0.00122, 0.000867, 0.000986, 0.00102000, 
	 0.00127, 0.00154, 0.00173, 0.00182, 0.00729, 0.00976, 0.0108},
	{0.00114, 0.00171, 0.00214, 0.00171, 0.00205, 0.0023, 0.00244, 0.0028,
	 0.00317, 0.0033, 0.00359, 0.0116, 0.0124}
       
	},

	{
	{-1.97000e-06,-2.23000e-06, 0.000249,0.000249, 0.000284, 0.000304, 
	 0.0004, 0.000429, 0.000417, 0.000431, 0.000412, 6.29000e-05,  
	 9.60000e-05},
	{-7.22000e-07, 6.36000e-05, 0.000144, 0.000143, 0.000207, 0.000309,
	 0.000361, 0.000396, 0.000512, 0.000276, 0.000257, -0.000254,-0.000321},
	{2.85000e-05, 5.36000e-05, 0.000139, 0.000148, 0.000216, 0.000249,  
	 0.000254, 0.000235, 0.000305, 5.72000e-05,  6.60000e-05, -0.000632, 
	 -0.000992},
	{ -3.87000e-05, 0.000981, 0.00187,0.000416, 0.000371, 0.000178, 
	  7.77000e-05, -4.19000e-05, -1.47000e-06, -0.000216, -0.000291, 
	  -0.00128,  -0.0016},
	{ -8.05000e-05, 0.00191, 0.00174,-0.000107, 2.28000e-05, 0.000201,  
	  8.78000e-05, -7.98000e-05, 9.64000e-05, -0.000417, -0.00059,  
	  -0.00125, -0.0013}
	}
    };

    agb = (AGBPARAM) malloc(sizeof(struct AGBContext));
    assert(agb != NULL);
    *agb = initagb;
    *pagb = agb;
    }

/* returns interpolated metal yield based on initial star mass and
 * it's metallicty
 */
double InterpolateMet(float fMass, float fMetals, double table[5][13])
{
    int i,j;
    double massFac, metFac;
    double t00,t01,t10,t11,t0,t1;
    float mass[13] = {0.9,1.0,1.3,1.3,1.5,1.7,2.0,2.5,3.0,4.0,5.0,7.0,8.0};
    float zs[5] = {0.001,0.004,0.008,0.02,0.04};
    for (i=0; mass[i] < fMass && i < 12; i++);
    for (j=0; zs[j] < fMetals && j < 4; j++);
    if(i > 0 && j > 0) {
	massFac = (fMass - mass[i-1]) / (mass[i] - mass[i-1]);
	metFac = (fMetals - zs[j-1]) / (zs[j] - zs[j-1]);
	t11 = table[j][i];
	t01 = table[j-1][i];
	t10 = table[j][i-1];
	t00 = table[j-1][i-1];
	t1 = metFac*(t11-t10) + t10;
	t0 = metFac*(t01 - t00) + t00;
	/*	printf("AGB fMass:  %g; i: %d; j: %d; t1: %g t0: %g; \nAGB massFac: %g; metFac: %g\n",
		fMass,i,j,t1,t0,massFac,metFac);*/
	return massFac*(t1-t0) + t0;
	} else if(i > 0) {
	massFac = (fMass - mass[i-1]) / (mass[i] - mass[i-1]);
	t1 = table[0][i];
	t0 = table[0][i-1];
	/*	printf("AGB fMass:  %g; i: %d; j: %d; t1: %g t0: %g; \nAGB massFac: %g; metFac: %g\n",
		fMass,i,j,t1,t0,massFac,metFac);*/
	return massFac*(t1-t0) + t0;
	} else if(j > 0) {
	metFac = (fMetals - zs[j-1]) / (zs[j] - zs[j-1]);
	t1 = table[j][0];
	t0 = table[j-1][0];
	/*	printf("AGB fMass:  %g; i: %d; j: %d; t1: %g t0: %g; \nAGB massFac: %g; metFac: %g\n",
		fMass,i,j,t1,t0,massFac,metFac);*/
	return metFac*(t1-t0) + t0;
	} else return table[0][0];
    }

void snCalcWindFeedback(SN sn, SFEvent sfEvent,
                        double dTime, /* current time in years */
                        double dDelta, /* length of timestep (years) */
                        FBEffects *fbEffects)
{
    double dMStarMin, dMStarMax;
    double dStarLtimeMin, dStarLtimeMax;
    double dMCumMin, dMCumMax,dMTot;
    double dMmin, dMmax;
    double dMassFracReturned;
    double dMDying;
    double dMFracC, dMFracN, dMFracO;
    AGBPARAM pagb = &(sn->AGBtabs);

    /* First determine if dying stars are between 0.9-8 Msolar

    * stellar lifetimes corresponding to beginning and end of 
    * current timestep with respect to starbirth time in yrs */
    dMmin=0.9;
    dMmax=8.0;
    dStarLtimeMin = dTime - sfEvent.dTimeForm; 
    dStarLtimeMax = dStarLtimeMin + dDelta;

    dMStarMin = dSTMStarLtime(&sn->ppdva, dStarLtimeMax, sfEvent.dMetals); 
    dMStarMax = dSTMStarLtime(&sn->ppdva, dStarLtimeMin, sfEvent.dMetals); 
    assert(dMStarMax >= dMStarMin);

    if (((dMStarMin < dMmax) && (dMStarMax > dMmin)) && dMStarMax > dMStarMin) {

	/* Mass Fraction returned to ISM taken from Weidemann, 1987, A&A 188 74 
	   then fit to function: MFreturned = 0.86 - exp(-Mass/1.1) 
	   More recent observations/theory Weidemann (2000), Kalirai
	   et al (2005 + 2008) don't seem to affect this much.
	*/
	dMassFracReturned=0.86-exp(-((dMStarMax+dMStarMin)/2.)/1.1);
	/* an attempt to integrate this failed badly
	  dMassFracReturned=0.86*(dMStarMax-dMStarMin) + 
	  1.1*(exp(-dMStarMax/1.1) - exp(-dMStarMin/1.1));*/
#ifdef MORE_METALS
	/* fits to van den Hoek + Groenewegen (1997) or Karakas (2010) */
	/*dMC = 5e-2*exp(-4.0 * pow(log((dMStarMax+dMStarMin)/2.0)-1.3,2.0) )/
	    (dMStarMax+dMStarMin)/2.0; 
	/* integrating this Carbon expression is not fun, so using
	   mean star mass.*/
	/*dMN = 1.7e-3 * (exp(dMStarMax/1.7) - exp(dMStarMin/1.7));
	dMO = 1e-3 * (dMStarMax-dMStarMin);
	printf("AGB dMC: %g; dMN: %g; dMO: %g\n",dMC,dMN,dMO);*/
	/* The above analytical expressions will be slow to
	   calculate, so interpolate tables instead. */

	/* Yield tables are mass fractions, */
	dMFracC = (InterpolateMet(dMStarMin, sfEvent.dMetals, pagb->Cej) + 
		   InterpolateMet(dMStarMax, sfEvent.dMetals, pagb->Cej))/2.0
	    /dMassFracReturned;
	dMFracN = (InterpolateMet(dMStarMin, sfEvent.dMetals, pagb->Nej) +
		   InterpolateMet(dMStarMax, sfEvent.dMetals, pagb->Nej))/2.0
	    /dMassFracReturned;
#endif
	dMFracO = (InterpolateMet(dMStarMin, sfEvent.dMetals, pagb->Oej) + 
		   InterpolateMet(dMStarMax, sfEvent.dMetals, pagb->Oej))/2.0
	    /dMassFracReturned;
	/*printf("AGB dMC: %g; dMN: %g; dMO: %g\n",dMC,dMN,dMO);*/

	dMCumMin = dMSCumMass(&sn->MSparam, dMStarMin);
	dMCumMax = dMSCumMass(&sn->MSparam, dMStarMax);
	dMTot = dMSCumMass(&sn->MSparam,0.0);
	/* Find out mass fraction of dying stars, then multiply by the original
	   mass of the star particle */
	if (dMTot > 0.0) dMDying = (dMCumMin - dMCumMax)/dMTot;
	else dMDying = 0.0;
	dMDying *= sfEvent.dMass;

	/* Figure out feedback effects */
	fbEffects->dMassLoss = dMDying * dMassFracReturned;
	fbEffects->dEnergy = 0.0;    
	/* Use star's metallicity for gas returned */
	fbEffects->dMIron = sfEvent.dMFracIron; 
	/* These will be multiplied by dMassLoss, so fraction is
	   appropriate quantity.*/
	fbEffects->dMOxygen = sfEvent.dMFracOxygen + dMFracO;
#ifdef MORE_METALS
	fbEffects->dMC = sfEvent.dMFracC + dMFracC;
	fbEffects->dMN = sfEvent.dMFracN + dMFracN;
	/* printf("AGB dMassLoss:  %g, dMDying: %g, fbC: %g; fbN: %g; fbO: %g\n",
	       dMDying,fbEffects->dMassLoss,
	       fbEffects->dMC, fbEffects->dMN,fbEffects->dMOxygen);
	*/
	fbEffects->dMNe = sfEvent.dMFracNe;
	fbEffects->dMMg = sfEvent.dMFracMg; 
	fbEffects->dMSi = sfEvent.dMFracSi; 
	fbEffects->dMetals = ( fbEffects->dMIron + fbEffects->dMOxygen +
			       fbEffects->dMC + fbEffects->dMN + 
			       fbEffects->dMNe + fbEffects->dMMg + 
			       fbEffects->dMSi);
#else
	fbEffects->dMetals = ( fbEffects->dMIron + fbEffects->dMOxygen );
#endif /* MORE_METALS */
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
	    }
    }
#endif
