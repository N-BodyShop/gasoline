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

void pkdFeedback(PKD pkd, FB fb, double dTime, double dDelta,
				 FBEffects *fbTotals)
{
    int i;
    PARTICLE *p;
    int n = pkdLocal(pkd);
    SFEvent sfEvent;
    FBEffects fbEffects;
    double dTotMassLoss;
    double dTotMetals;
    double dDeltaYr;
    SN sn;
    int j;
    
    snInitialize(&sn);
    snInitConstants(sn);
    
    for(i = 0; i < FB_NFEEDBACKS; i++) {
		fbTotals[i].dMassLoss = 0.0;
		fbTotals[i].dEnergy = 0.0;
		fbTotals[i].dMetals = 0.0;
		}
    dTime *= fb->dSecUnit/SEC_YR;
    dDeltaYr = dDelta*fb->dSecUnit/SEC_YR;
    
    for(i = 0; i < n; ++i) {
        p = &pkd->pStore[i];
#ifdef COOLDEBUG
		if (p->iOrder == 842079) {
			fprintf(stderr,"Particle %i in pStore[%i]\n",p->iOrder,(int) (p-pkd->pStore));
			}
		assert(p->u >= 0);
		assert(p->uPred >= 0);
#endif
		if(pkdIsStar(pkd, p)) {
			dTotMassLoss = 0.0;
			dTotMetals = 0.0;
			p->fESNrate = 0.0;
	    
			if(fb->dInitStarMass > 0.0)
			    sfEvent.dMass = fb->dInitStarMass*fb->dGmUnit/MSOLG;
			else
			    sfEvent.dMass = p->fMass*fb->dGmUnit/MSOLG;
			sfEvent.dTimeForm = p->fTimeForm*fb->dSecUnit/SEC_YR;;
			sfEvent.dMetals = p->fMetals;
	    
			/*
			 * Call all the effects in order and accumulate them.
			 */
			for(j = 0; j < FB_NFEEDBACKS; j++) {
				switch (j) {
				case FB_SNII:
					snCalcSNIIFeedback(sn, sfEvent, dTime,
									   dDeltaYr, &fbEffects);
					break;
				case FB_SNIA:
					snCalcSNIaFeedback(sn, sfEvent, dTime,
									   dDeltaYr, &fbEffects);
					break;
				case FB_WIND:
					snCalcWindFeedback(sn, sfEvent, dTime,
									   dDeltaYr, &fbEffects);
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
	    
				fbTotals[j].dMassLoss += fbEffects.dMassLoss;
				fbTotals[j].dEnergy += fbEffects.dEnergy*fbEffects.dMassLoss;
				fbTotals[j].dMetals += fbEffects.dMetals*fbEffects.dMassLoss;
				}
	    
			/*
			 * Modify star particle
			 */
			assert(p->fMass > dTotMassLoss);
	    
			p->fMass -= dTotMassLoss;
			p->fMSN = dTotMassLoss;
			if(dTotMassLoss != 0.0) {
				p->fSNMetals = dTotMetals/dTotMassLoss;
				p->fESNrate /= dTotMassLoss*dDelta; /* convert to rate */
				}
			else {
				p->fSNMetals = 0.0;
				p->fESNrate = 0.0;
				}
			}
	    
	
		else if(pkdIsGas(pkd, p))
			assert(p->u >= 0);
		    assert(p->uPred >= 0);
			p->fESNrate = 0;	/* reset SN heating rate of gas to zero */
		}
    snFree(sn);

#ifdef COOLDEBUG
    for(i=0;i<pkdLocal(pkd);++i) {
		p = &pkd->pStore[i];
		if (p->iOrder == 842079) fprintf(stderr,"Particle %i in pStore[%i] after\n",p->iOrder,(int) (p-pkd->pStore));
		assert(p->u >= 0);
		assert(p->uPred >= 0);
		}
#endif
	}

void snCalcWindFeedback(SN sn, SFEvent sfEvent,
						double dTime, /* current time in years */
						double dDelta, /* length of timestep (years) */
						FBEffects *fbEffects)
{
    fbEffects->dMassLoss = 0.0;
    fbEffects->dEnergy = 0.0;
    fbEffects->dMetals = 0.0;
    }

void snCalcUVFeedback(SN sn, SFEvent sfEvent,
					  double dTime, /* current time in years */
					  double dDelta, /* length of timestep (years) */
					  FBEffects *fbEffects)
{
    fbEffects->dMassLoss = 0.0;
    fbEffects->dEnergy = 0.0;
    fbEffects->dMetals = 0.0;
    }

#endif
