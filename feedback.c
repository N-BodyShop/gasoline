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
            p->fESNrate = 0.0;
            p->fNSN = 0.0;

            if(fb->dInitStarMass > 0.0)
                sfEvent.dMass = fb->dInitStarMass*fb->dGmUnit/MSOLG;
            else
                sfEvent.dMass = p->fMassForm*fb->dGmUnit/MSOLG;
            sfEvent.dTimeForm = p->fTimeForm*fb->dSecUnit/SEC_YR;;
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
                switch (j) {
                case FB_SNII:
                    snCalcSNIIFeedback(sn, sfEvent, dTime,
                                        dDeltaYr, &fbEffects);
                    break;
                case FB_SNIA:
                    snCalcSNIaFeedback(sn, sfEvent, dTime,
                                        dDeltaYr, &fbEffects);
                    dSNIaMassStore=fbEffects.dMassLoss;

                    break;
                case FB_WIND:
                    snCalcWindFeedback(sn, sfEvent, dTime,
                                        dDeltaYr, &fbEffects);
                    fbEffects.dMassLoss -= dSNIaMassStore;
                    /*printf("Wind, SNaI Mass Loss: %d   %d\n",fbEffects.dMassLoss,dSNIaMassStore); */
                    
                    break;
                case FB_UV:
                    snCalcUVFeedback(sn, sfEvent, dTime, dDeltaYr,
                                                     &fbEffects);
                    break;
                default:
                    assert(0);
                    }

                if( sn->dESN > 0.0) p->fNSN += fbEffects.dEnergy * MSOLG*fbEffects.dMassLoss/sn->dESN;
                else p->fNSN += 0.0;
                fbEffects.dMassLoss *= MSOLG/fb->dGmUnit;
                fbEffects.dEnergy /= fb->dErgPerGmUnit;

                dTotMassLoss += fbEffects.dMassLoss;
                p->fESNrate += fbEffects.dEnergy*fbEffects.dMassLoss;
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
            p->fESNrate /= dDelta; /* convert to rate */
            }
	    
	
		else if(pkdIsGas(pkd, p)){
                        assert(p->u >= 0.0);
                        assert(p->uPred >= 0.0);
			p->fESNrate = 0.0;	/* reset SN heating rate of gas to zero */
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
  /*printf("");  Intel optimizer needs this to avoid floating point exception on sharks. */

    /* Mass Fraction returned to ISM taken from Weidermann, 1987, A&A 188 74 
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
    /* Use stars metallicity for gas returned */
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
