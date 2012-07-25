#ifdef STARFORM
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "pkd.h"
#include "feedback.h"
#include "supernova.h"

void snCalcUVFeedback(SN sn, SFEvent sfEvent,
		      double dTime, /* current time in years */
		      double dDelta, /* length of timestep (years) */
		      FBEffects *fbEffects);

/*
 * Feedback module for GASOLINE
 */

double mod(double a, int b) {
    /* printf("MOD: a=%g, b=%d, floor(a/b)=%g\n",a,b,floor(a/b));*/
    return (a-b*floor(a/b));
    }

void fbInitialize(FB *pfb)
{
    FB fb;
    
    fb = (FB) malloc(sizeof(struct fbContext));
    assert(fb != NULL);
    
    fb->dGmUnit = 0.0;
    fb->dSecUnit = 0.0;
    fb->dErgPerGmUnit = 0.0;
    *pfb = fb;
    }

void pkdFeedback(PKD pkd, FB fb, SN sn, double dTime, double dDelta,
		 FBEffects *fbTotals)
{
    int i;
    PARTICLE *p;
    int n = pkdLocal(pkd);
    SFEvent sfEvent;
    FBEffects fbEffects;
    double dTotMassLoss;
    double dTotMetals;
    double dTotMOxygen;
    double dTotMIron;
#ifdef MORE_METALS
    double dTotMC, dTotMN, dTotMNe, dTotMMg, dTotMSi;
#endif
    double dDeltaYr;
    double dSNIaMassStore;
    double dNSNII, dProb, dStarAge, dMinAge;
    double dRandomNum;
    int j;
        
    for(i = 0; i < FB_NFEEDBACKS; i++) {
	fbTotals[i].dMassLoss = 0.0;
	fbTotals[i].dEnergy = 0.0;
	fbTotals[i].dMetals = 0.0;
	fbTotals[i].dMIron = 0.0;
	fbTotals[i].dMOxygen = 0.0;
#ifdef MORE_METALS
	fbTotals[i].dMC = 0.0;
	fbTotals[i].dMN = 0.0;
	fbTotals[i].dMNe = 0.0;
	fbTotals[i].dMMg = 0.0;
	fbTotals[i].dMSi = 0.0;
#endif
	}
    dTime *= fb->dSecUnit/SEC_YR;
    dDeltaYr = dDelta*fb->dSecUnit/SEC_YR;
    
    for(i = 0; i < n; ++i) {
	p = &pkd->pStore[i];
	if(pkdIsStar(pkd, p) && p->fTimeForm >= 0.0) {
	    dTotMassLoss = 0.0;
	    dTotMetals = 0.0;
	    dTotMOxygen = 0.0;
	    dTotMIron = 0.0;
#ifdef MORE_METALS
	    dTotMC = 0.0;
	    dTotMN = 0.0;
	    dTotMNe = 0.0;
	    dTotMMg = 0.0;
	    dTotMSi = 0.0;
#endif
	    p->fESNrate = 0.0;
	    p->fNSN = 0.0;

	    sfEvent.dMass = p->fMassForm*fb->dGmUnit/MSOLG;
	    sfEvent.dTimeForm = p->fTimeForm*fb->dSecUnit/SEC_YR;
	    dStarAge = dTime - sfEvent.dTimeForm;
	    sfEvent.dMetals = p->fMetals;
	    sfEvent.dMFracOxygen = p->fMFracOxygen;
	    sfEvent.dMFracIron = p->fMFracIron;
#ifdef MORE_METALS
	    sfEvent.dMFracC = p->fMFracC;
	    sfEvent.dMFracN = p->fMFracN;
	    sfEvent.dMFracNe = p->fMFracNe;
	    sfEvent.dMFracMg = p->fMFracMg;
	    sfEvent.dMFracSi = p->fMFracSi;
#endif
	    /*
	     * Call all the effects in order and accumulate them.
	     */
	    dSNIaMassStore=0.0;  /* Stores mass loss of Ia so as
				    not to double count it in 
				    wind feedback */
                                    
	    for(j = 0; j < FB_NFEEDBACKS; j++) {
		dNSNII = 0;
		switch (j) {
		case FB_SNII:
		    snCalcSNIIFeedback(sn, sfEvent, dTime,
				       dDeltaYr, &fbEffects);
		    if( sn->dESN > 0.0)
			dNSNII = fbEffects.dEnergy * MSOLG*fbEffects.dMassLoss/
			    sn->dESN;
		    /* Blow winds before SN (power of winds ~ power of SN) */
		    dMinAge = dSTLtimeMStar(&sn->ppdva, sn->dMSNIImax, 
					    sfEvent.dMetals); 
		    if (dNSNII > 0 && sn->iNSNIIQuantum > 0 && 
			dStarAge > dMinAge) {
			/* Make sure only a iNSNIIQuantum number of
			 * SNII go off at a time */
			dProb = mod(dNSNII, sn->iNSNIIQuantum)/
			    sn->iNSNIIQuantum;
			dRandomNum = (rand()/((double) RAND_MAX));
			/*	  printf("Random Number = %g\n",dRandomNum);*/
			if(dRandomNum < dProb) /* SN occurred */ {
			    /* Adds missing part to make up quantum */
			    p->fNSN = dNSNII + (1.-dProb)*sn->iNSNIIQuantum;
			    /*printf("NSN: +factor=%g   dNSNII=%g  result=%g fNSN=%g\n",(1.-dProb)*sn->iNSNIIQuantum,dNSNII,dNSNII + (1.-dProb)*sn->iNSNIIQuantum,p->fNSN);*/
			    } else {
				p->fNSN = dNSNII - dProb*sn->iNSNIIQuantum;
				/*printf("NSN: -factor=%g   dNSNII=%g  result=%g fNSN=%g\n",dProb*sn->iNSNIIQuantum,dNSNII,dNSNII - dProb*sn->iNSNIIQuantum,p->fNSN);*/
				}
			if(p->fNSN < sn->iNSNIIQuantum) p->fNSN = 0;
			fbEffects.dEnergy = p->fNSN*sn->dESN/(MSOLG*fbEffects.dMassLoss);   
			} 
		    else if(dStarAge < dMinAge && sn->iNSNIIQuantum) 
			p->fNSN = dNSNII;
		    else p->fNSN += dNSNII;
		    break;
		case FB_SNIA:
		    snCalcSNIaFeedback(sn, sfEvent, dTime,
				       dDeltaYr, &fbEffects);
		    dSNIaMassStore=fbEffects.dMassLoss;
		    break;
		case FB_WIND:
		    snCalcWindFeedback(sn, sfEvent, dTime,
				       dDeltaYr, &fbEffects);
		    if(dSNIaMassStore < fbEffects.dMassLoss)
			fbEffects.dMassLoss -= dSNIaMassStore;
		    break;
		case FB_UV:
		    snCalcUVFeedback(sn, sfEvent, dTime, dDeltaYr,
				     &fbEffects);
		    break;
		default:
		    assert(0);
		    }

		fbEffects.dMassLoss *= MSOLG/fb->dGmUnit;
		fbEffects.dEnergy /= fb->dErgPerGmUnit;

		dTotMassLoss += fbEffects.dMassLoss;
		p->fESNrate += fbEffects.dEnergy*fbEffects.dMassLoss;
		dTotMetals += fbEffects.dMetals*fbEffects.dMassLoss;
		dTotMOxygen += fbEffects.dMOxygen*fbEffects.dMassLoss;
		dTotMIron += fbEffects.dMIron*fbEffects.dMassLoss;
#ifdef MORE_METALS
		dTotMC += fbEffects.dMC*fbEffects.dMassLoss;
		dTotMN += fbEffects.dMN*fbEffects.dMassLoss;
		dTotMNe += fbEffects.dMNe*fbEffects.dMassLoss;
		dTotMMg += fbEffects.dMMg*fbEffects.dMassLoss;
		dTotMSi += fbEffects.dMSi*fbEffects.dMassLoss;
#endif

		fbTotals[j].dMassLoss += fbEffects.dMassLoss;
		fbTotals[j].dEnergy += fbEffects.dEnergy*fbEffects.dMassLoss;
		fbTotals[j].dMetals += fbEffects.dMetals*fbEffects.dMassLoss;
		fbTotals[j].dMIron += fbEffects.dMIron*fbEffects.dMassLoss;
		fbTotals[j].dMOxygen += fbEffects.dMOxygen*fbEffects.dMassLoss;
#ifdef MORE_METALS
		fbTotals[j].dMC += fbEffects.dMC*fbEffects.dMassLoss;
		fbTotals[j].dMN += fbEffects.dMN*fbEffects.dMassLoss;
		fbTotals[j].dMNe += fbEffects.dMNe*fbEffects.dMassLoss;
		fbTotals[j].dMMg += fbEffects.dMMg*fbEffects.dMassLoss;
		fbTotals[j].dMSi += fbEffects.dMSi*fbEffects.dMassLoss;
#endif
		}

	    /*
	     * Modify star particle
	     */
	    assert(p->fMass > dTotMassLoss);

	    p->fMass -= dTotMassLoss;
	    p->fMSN = dTotMassLoss;
	    /* The SNMetals and ESNrate used to be specific
	       quantities, but we are making them totals as
	       they leave the stars so that they are easier
	       to divvy up among the gas particles in 
	       distSNEnergy in smoothfcn.c.  These quantities
	       will be converted back to specific quantities when
	       they are parts of gas particles. */
	    p->fSNMetals = dTotMetals;
	    p->fMIronOut = dTotMIron;
	    p->fMOxygenOut = dTotMOxygen;
#ifdef MORE_METALS
	    p->fMCOut = dTotMC;
	    p->fMNOut = dTotMN;
	    p->fMNeOut = dTotMNe;
	    p->fMMgOut = dTotMMg;
	    p->fMSiOut = dTotMSi;
#endif
	    p->fESNrate /= dDelta; /* convert to rate */
#ifdef  RADIATIVEBOX /* Calculates LW radiation from a stellar particle of a given age and mass (assumes Kroupa IMF), CC */
	    double  a0 = -84550.812,
	      a1 =  54346.066,
	      a2 = -13934.144,
	      a3 =  1782.1741,
	      a4 = -113.68717,
	      a5 =  2.8930795;

	    double a0old =  70.908586,
	      a1old = -4.0643123;
	    
	    double Alog101, Alog102, dLW1, dLW2;
	    Alog101 = dStarAge;
	    if (Alog101 < 1e7) Alog101 = 7; /*used Alog10 = 1e6, at one point*/
	    else Alog101 = log10(Alog101);
	    if (Alog101 < 9.0) dLW1 = pow(10,a0
		       + a1*Alog101
		       + a2*Alog101*Alog101
		       + a3*Alog101*Alog101*Alog101
		       + a4*Alog101*Alog101*Alog101*Alog101
		       + a5*Alog101*Alog101*Alog101*Alog101*Alog101 - 30.0);
	    else dLW1 = pow(10,a0old + a1old*Alog101 - 30.0);

	    Alog102 = dStarAge + dDeltaYr;
	    if (Alog102 < 1e7) Alog102 = 7; /*used Alog10 = 1e6, at one point*/
	    else Alog102 = log10(Alog102);
	    if (Alog102 < 9.0) dLW2 = pow(10,a0
		       + a1*Alog102
		       + a2*Alog102*Alog102
		       + a3*Alog102*Alog102*Alog102
		       + a4*Alog102*Alog102*Alog102*Alog102
		       + a5*Alog102*Alog102*Alog102*Alog102*Alog102 - 30.0); /*Divide by 1e30 to keep things in bounds*/
	    else dLW2 = pow(10,a0old + a1old*Alog102 - 30.0);

	    p->CoolParticle.dLymanWerner = pkd->Cool->dMsolUnit*pkd->Cool->dInitStarMass*(dLW1 + dLW2)/2; 
#endif
	    }
	else if(pkdIsGas(pkd, p)){
	    assert(p->u >= 0.0);
	    assert(p->uPred >= 0.0);
	    p->fESNrate = 0.0;	/* reset SN heating rate of gas to zero */
	    }
	}

    }

void snCalcUVFeedback(SN sn, SFEvent sfEvent,
		      double dTime, /* current time in years */
		      double dDelta, /* length of timestep (years) */
		      FBEffects *fbEffects)
{
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

#endif
