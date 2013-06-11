#ifdef STARFORM
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "pkd.h"
#include "feedback.h"
#include "supernova.h"

void snCalcWindFeedback(SN sn, SFEvent sfEvent,
                        double dTime, /* current time in years */
                        double dDelta, /* length of timestep (years) */
                        FBEffects *fbEffects);

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
	    p->uDotFB = 0.0;
	    p->fNSN = 0.0;
#ifdef FBPARTICLE
            {
            double tstart = 4e6, tend = 30e6; /* yr */
            double mdotonmstar = 10*10/100./(tend-tstart); /* gm / gm / yr */
            double edotonmstar = 1e51/(100*2e33)/(tend-tstart); /* erg / gm / yr */
            double mFB = 2.5e-2*p->fMassForm; /* mass of fb particles (code units) */
            double tFB0,tFB1,nFac;
            int nFB0,nFB1;

            tFB1 =  dTime-p->fTimeForm*fb->dSecUnit/SEC_YR;
            tFB0 =  tFB1 - dDeltaYr;
            if (tFB1 > tstart && tFB0 < tend) {
                nFac = mdotonmstar*p->fMassForm/mFB;
                nFB0 = floor(nFac*(tFB0 > tstart ? tFB0-tstart : 0)+0.5);
                nFB1 = floor(nFac*(tFB1 < tend ? tFB1-tstart : tend-tstart)+0.5);
                printf("FBP: %d %g %g %g %d %d\n",p->iOrder,dTime,tFB0,tFB1,nFB0,nFB1);
                while (nFB1 > nFB0) {
                    PARTICLE pNew;
                    pNew = *p;
                    TYPEReset(&pNew, TYPE_STAR|TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE|TYPE_ACTIVE);
                    TYPESet(&pNew, TYPE_GAS);
                    
                    /* For mass conservation would need to do this ... */
                    /*p->fMass -= mFB;
                      assert(p->fMass > 0);*/

                    pNew.fMass = mFB;
                    pNew.u = (edotonmstar/mdotonmstar)/fb->dErgPerGmUnit;
                    pNew.uPred = pNew.u;
                    pNew.c = sqrt(pNew.u);
                    pNew.PoverRho2 = 0;
                    pNew.uDot = 0;
                    pNew.uDotPdV = 0;
                    pNew.uDotAV = 0;
                    pNew.uDotDiff = 0;
                    pNew.uDotFB = 0;
                    pkdNewParticle(pkd, pNew);
                    printf("FBP: Particle made %g %g\n",mFB,pNew.u);
                    nFB1--;
                    }
                }
            }
        continue; /* skip normal feedback */
#endif                

	    sfEvent.dMass = p->fMassForm*fb->dGmUnit/MSOLG;
	    sfEvent.dTimeForm = p->fTimeForm*fb->dSecUnit/SEC_YR;
	    dStarAge = dTime - sfEvent.dTimeForm;
	    sfEvent.dMetals = p->fMetals;
	    sfEvent.dMFracOxygen = p->fMFracOxygen;
	    sfEvent.dMFracIron = p->fMFracIron;

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
/*			    printf("NSN: +factor=%g   dNSNII=%g  result=%g fNSN=%g\n",(1.-dProb)*sn->iNSNIIQuantum,dNSNII,dNSNII + (1.-dProb)*sn->iNSNIIQuantum,p->fNSN);*/
			    } 
			else {
			    p->fNSN = dNSNII - dProb*sn->iNSNIIQuantum;
/*			    printf("NSN: -factor=%g   dNSNII=%g  result=%g fNSN=%g\n",dProb*sn->iNSNIIQuantum,dNSNII,dNSNII - dProb*sn->iNSNIIQuantum,p->fNSN);*/
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
		p->uDotFB += fbEffects.dEnergy*fbEffects.dMassLoss;
		dTotMetals += fbEffects.dMetals*fbEffects.dMassLoss;
		dTotMOxygen += fbEffects.dMOxygen*fbEffects.dMassLoss;
		dTotMIron += fbEffects.dMIron*fbEffects.dMassLoss;

		fbTotals[j].dMassLoss += fbEffects.dMassLoss;
		fbTotals[j].dEnergy += fbEffects.dEnergy*fbEffects.dMassLoss;
		fbTotals[j].dMetals += fbEffects.dMetals*fbEffects.dMassLoss;
		fbTotals[j].dMIron += fbEffects.dMIron*fbEffects.dMassLoss;
		fbTotals[j].dMOxygen += fbEffects.dMOxygen*fbEffects.dMassLoss;
		}

	    /*
	     * Modify star particle
	     */
//        fprintf(stderr,"Mass dTotMassLoss %d %g %g  %g %g\n",p->iOrder,p->fMass,dTotMassLoss,dTime/(fb->dSecUnit/SEC_YR),p->fTimeForm);
	    assert(p->fMass > dTotMassLoss);

	    p->fMass -= dTotMassLoss;
	    p->fMSN = dTotMassLoss;
	    /* The SNMetals and uDotFB (SN rate) used to be specific
	       quantities, but we are making them totals as
	       they leave the stars so that they are easier
	       to divvy up among the gas particles in 
	       distSNEnergy in smoothfcn.c.  These quantities
	       will be converted back to specific quantities when
	       they are parts of gas particles. */
	    p->fSNMetals = dTotMetals;
	    p->fMIronOut = dTotMIron;
	    p->fMOxygenOut = dTotMOxygen;
	    p->uDotFB /= dDelta; /* convert to rate */
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
	    if (dStarAge < 0) p->CoolParticle.dLymanWerner = 0; /*If BH, no LW radiation*/
	    else {
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
	    }
#endif
	    }
	else if(pkdIsGas(pkd, p)){
	    assert(p->u >= 0.0);
	    assert(p->uPred >= 0.0);
	    p->uDotFB = 0.0;	/* reset SN heating rate of gas to zero */
	    }
	}

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

    /* First determine if dying stars are between 1-8 Msolar

    * stellar lifetimes corresponding to beginning and end of 
    * current timestep with respect to starbirth time in yrs */
    dMmin=1.0;
    dMmax=8.0;
    dStarLtimeMin = dTime - sfEvent.dTimeForm; 
    dStarLtimeMax = dStarLtimeMin + dDelta;

    dMStarMin = dSTMStarLtime(&sn->ppdva, dStarLtimeMax, sfEvent.dMetals); 
    dMStarMax = dSTMStarLtime(&sn->ppdva, dStarLtimeMin, sfEvent.dMetals); 
    assert(dMStarMax >= dMStarMin);

    if (((dMStarMin < dMmax) && (dMStarMax > dMmin)) && dMStarMax > dMStarMin) {

	/* Mass Fraction returned to ISM taken from Weidemann, 1987, A&A 188 74 
	   then fit to function: MFreturned = 0.86 - exp(-Mass/1.1) */
	dMassFracReturned=0.86-exp(-((dMStarMax+dMStarMin)/2.)/1.1);

	dMCumMin = dMSCumMass(&sn->MSparam, dMStarMin);
	dMCumMax = dMSCumMass(&sn->MSparam, dMStarMax);
	dMTot = dMSCumMass(&sn->MSparam,0.0);
	/* Find out mass fraction of dying stars, then multiply by the original
	   mass of the star particle */
	if (dMTot == 0.0){
	    dMDying = 0.0;
	    } else { 
		dMDying = (dMCumMin - dMCumMax)/dMTot;
		}
	dMDying *= sfEvent.dMass;

	/* Figure out feedback effects */
	fbEffects->dMassLoss = dMDying * dMassFracReturned;
	fbEffects->dEnergy = 0.0;    
	/* Use star's metallicity for gas returned */
	fbEffects->dMetals = sfEvent.dMetals; 
	fbEffects->dMIron = sfEvent.dMFracIron; 
	fbEffects->dMOxygen = sfEvent.dMFracOxygen; 

	} else {
	    fbEffects->dMassLoss = 0.0;
	    fbEffects->dEnergy = 0.0;    
	    fbEffects->dMetals = 0.0; 
	    fbEffects->dMIron = 0.0;
	    fbEffects->dMOxygen = 0.0;
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
    }

#endif
