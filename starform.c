#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "pkd.h"
#include "starform.h"
#include "millerscalo.h"

#define max(A,B) ((A) > (B) ? (A) : (B))
#define min(A,B) ((A) < (B) ? (A) : (B))

#ifdef STARFORM

#ifdef NOCOOLING
#error "STARFORM requires Cooling on"
#endif


/*
 * Star forming module for GASOLINE
 */

void stfmInitialize(STFM *pstfm)
{
    STFM stfm;
    
    stfm = (STFM) malloc(sizeof(struct stfmContext));
    assert(stfm != NULL);
    
    stfm->dPhysDenMin = 0;
    stfm->dOverDenMin = 0;
    stfm->dTempMax = 0;
    stfm->dCStar = 0;
    stfm->dSecUnit = 0;
    stfm->dGmPerCcUnit = 0;
    stfm->dGmUnit = 0;
    stfm->dStarEff = 0.0;
	stfm->dInitStarMass = 0.0;
    stfm->dMinGasMass = 0.0;
    stfm->dMinMassFrac = 0.0;
    stfm->dMaxStarMass = 0.0;
    *pstfm = stfm;
    }

/*
     taken from TREESPH and modified greatly.
     Uses the following formula for the star formation rate:

              d(ln(rhostar))/dt=cstar/tdyn

*/


void stfmFormStars(STFM stfm, PKD pkd, PARTICLE *p,
		   double dTime, /* current time */
		   int *nFormed, /* number of stars formed */
		   double *dMassFormed,	/* mass of stars formed */
		   int *nDeleted) /* gas particles deleted */
{
    double tdyn;
    double tform;
    double tsound;
    double tcool;
    double E,T;
    COOL *cl = pkd->Cool;
    double dExp = 1.0/(1.0 + cl->z);
    double dCosmoFac = dExp*dExp*dExp;
    PARTICLE starp;
    double dMprob;
    double dDeltaM;
    double l_jeans2;
    int small_jeans = 0;
    int j;
    
    /*  This version of the code has only three conditions unless 
	-D SFCONDITIONS is set:
	  converging flow (p->divv < 0)
	  T < dTempMax
	  density > dOverdenmin && density > dPhysDenMin
	Anil - Nov. 2003
    */

    /*
     * Is particle in convergent part of flow?
     */
    
#ifndef DIVVOFF
    if(p->divv >= 0.0)
	return;
#endif /*DIVVOFF*/
    /*
     * Determine dynamical time.
     */
    tdyn = 1.0/sqrt(4.0*M_PI*p->fDensity/dCosmoFac);

    /*
     * Determine cooling time.
     */


    T = CoolCodeEnergyToTemperature( cl, &p->CoolParticle, p->u );
#if (0)
    E = CoolCodeEnergyToErgPerGm( cl, p->u );
    tcool = E/(-CoolHeatingRate( cl, &p->CoolParticle, T, 
		 CodeDensityToComovingGmPerCc(cl,p->fDensity ), p->fMetals )
                    -CoolCodeWorkToErgPerGmPerSec( cl, p->PdV ));
    tcool = CoolSecondsToCodeTime( cl, tcool ); 
    printf("tcool %i: %g %g %g\n",p->iOrder,T,p->fDensity,tcool);
#endif
    tcool = p->u/(
#ifdef DENSITYU
	-CoolEdotInstantCode( cl, &p->CoolParticle, p->u, p->fDensityU, p->fMetals, p->r )
#else
	-CoolEdotInstantCode( cl, &p->CoolParticle, p->u, p->fDensity, p->fMetals, p->r )
#endif
	-p->PdV );
#ifdef CHECKSF
    p->tOff = CoolCodeTimeToSeconds( cl, p->fTimeCoolIsOffUntil - dTime)/3.1557e7;  /* years */
    p->tcool = CoolCodeTimeToSeconds( cl, tcool)/3.1557e7;
    p->tdyn = CoolCodeTimeToSeconds( cl, tdyn)/3.1557e7;
    p->ratiosounddyn = sqrt(0.25*p->fBall2)/p->c/tdyn;
    p->l_jeans = sqrt(M_PI*p->c*p->c/p->fDensity*dCosmoFac); /* Why not comoving? */
    p->small_jeans = small_jeans;
#endif
#ifdef SFCONDITIONS
    if(tcool < 0.0 && T > stfm->dTempMax) return;
    /*
     * Determine sound crossing time.
     */
    tsound = sqrt(0.25*p->fBall2)/p->c;

    /* 
     * criteria that stars form if the Jean's length is less than the
     * softening
     */
    l_jeans2 = M_PI*p->c*p->c/p->fDensity*dCosmoFac;
#if (0)
/* Old code: problem -- compares physics L_J to comoving softening */
    if (l_jeans2 < p->fSoft*p->fSoft*stfm->dSoftMin*stfm->dSoftMin) 
        small_jeans = 1;
    /*printf("tsound:  %g  c_sound:  %g  Temp:  %g  tdyn:  %g\n",tsound,p->c,T,tdyn);*/
#ifdef CHECKSF
    p->small_jeans = small_jeans;
#endif
    if (!small_jeans && tsound <= tdyn)
        return;

#else
/* New code: physical L_J vs. physics smoothing length (with multiplier) */
    if (l_jeans2 >= 0.25*p->fBall2*dExp*dExp*stfm->dSoftMin*stfm->dSoftMin) return;

#ifdef CHECKSF
    p->small_jeans = 1;
#endif
#endif

#else
    if(T > stfm->dTempMax) return;
#endif /*SFCONDITIONS*/

    /*
     * Determine if this particle satisfies all conditions.
     */
    
    if(p->fDensity < stfm->dOverDenMin ||
       p->fDensity/dCosmoFac < stfm->dPhysDenMin)
	return;

    if(tcool < 0.0 || tdyn > tcool || T < stfm->dTempMax)
	tform = tdyn;
    else
	tform = tcool;
    
    dMprob = 1.0 - exp(-stfm->dCStar*stfm->dDeltaT/tform);
    
    /*
     * Decrement mass of particle.
     */

    if (stfm->dInitStarMass > 0) 
        dDeltaM = stfm->dInitStarMass;
    else 
        dDeltaM = p->fMass*stfm->dStarEff;

    /* No negative or very tiny masses please! */
    if ( (dDeltaM > p->fMass) ) dDeltaM = p->fMass;

    if(dMprob*p->fMass < dDeltaM*(rand()/((double) RAND_MAX)))
	return;

    p->fMass -= dDeltaM;
	assert(p->fMass >= 0.0);

    /* 
     * Note on number of stars formed:
     * n = log(dMinGasMass/dInitMass)/log(1-dStarEff) = max no. stars 
     * formed per gas particle, e.g. if min gas mass = 10% initial mass,
     * dStarEff = 1/3, max no. stars formed = 6 (round up so gas mass 
     * goes below min gas mass)
     */

    starp = *p; 		/* grab copy before possible deletion */

    if(p->fMass < stfm->dMinGasMass) {
		(*nDeleted)++;
		pkdDeleteParticle(pkd, p);
		}

    /*
     * form star
     */

    starp.fMass = dDeltaM;
    starp.fTimeForm = dTime;
    starp.fMassForm = dDeltaM;
    starp.fBallMax = 0.0;
    /*
     * Save quantities
     */
    for(j = 0; j < 3; j++) {
	starp.rForm[j] = starp.r[j];
	starp.vForm[j] = starp.v[j];
	}
    starp.u = T;
    starp.fNSNtot = 0.0;
    starp.iGasOrder = starp.iOrder; /* iOrder gets reassigned in
				       NewParticle() */

	/* NB: It is important that the star inherit special properties of the gas
	   particle such as being a target for movies or other tracing
	   Thus: Do not remove all the TYPE properties -- just the gas specific ones */
    TYPEReset(&starp, TYPE_GAS|TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE|TYPE_ACTIVE);
    TYPESet(&starp, TYPE_STAR);
    TYPEReset(&starp, TYPE_NbrOfACTIVE); /* just a precaution */
    
    (*nFormed)++;
    *dMassFormed += dDeltaM;
    
    pkdNewParticle(pkd, starp);    

}

void pkdFormStars(PKD pkd, STFM stfm, double dTime, int *nFormed,
		  double *dMassFormed, int *nDeleted)
{
    int i;
    PARTICLE *p;
    int n = pkdLocal(pkd);
    
    *nFormed = 0;
    *nDeleted = 0;
    *dMassFormed = 0.0;
    
    for(i = 0; i < n; ++i) {
        p = &pkd->pStore[i];
        if(pkdIsGas(pkd, p))
            stfmFormStars(stfm, pkd, p, dTime, nFormed, dMassFormed, nDeleted);
        assert(p->u >= 0);
        assert(p->uPred >= 0);
        assert(p->fMass >= 0);
        }
    }

#endif

#ifdef SIMPLESF
void pkdSimpleStarForm(PKD pkd, double dRateCoeff, double dTMax, double dDenMin, double dDelta, double dTime,
					   double dInitStarMass, double dESNPerStarMass, double dtCoolingShutoff, int bdivv,
					   int *nFormed, /* number of stars formed */
					   double *dMassFormed,	/* mass of stars formed */
					   int *nDeleted) /* gas particles deleted */
{
    int i,j;
    PARTICLE *p;
    int n = pkdLocal(pkd);

    double T;
    COOL *cl = pkd->Cool;

	double mstardot;
    PARTICLE starp;
    
    *nFormed = 0;
    *nDeleted = 0;
    *dMassFormed = 0.0;

    for(i = 0; i < n; ++i) {
        p = &pkd->pStore[i];
/*		if(TYPEFilter(p,TYPE_GAS|TYPE_ACTIVE,TYPE_GAS|TYPE_ACTIVE)) { */
		if(TYPEFilter(p,TYPE_GAS,TYPE_GAS)) {
			/* Make sure cool again up to date */
			if (p->fTimeForm < dTime) p->fTimeForm = dTime;

			/* Ref: stfmFormStars(stfm, pkd, p, dTime, nFormed, dMassFormed, nDeleted); */
			/* Is particle in convergent part of flow?  */
			if (p->fDensity < dDenMin || (bdivv && p->divv >= 0.0)) continue;
			
			if ((T = CoolCodeEnergyToTemperature(pkd->Cool,&p->CoolParticle, p->u)) > dTMax) continue; 
			
			mstardot = dRateCoeff*sqrt(p->fDensity)*(p->fMass-p->fMassStar); /* Predictor corrector for second order? */

			p->fMassStar += mstardot*dDelta; /* sanity checks occur later */

			/* Star formation event? */
			if (p->fMassStar > dInitStarMass) { 
				starp = *p; 		/* grab copy before possible deletion */
				starp.fESN = p->u + dESNPerStarMass; /* ESN per unit mass -- includes gas internal energy */

				if (p->fMassStar > p->fMass-0.5*dInitStarMass) {
					starp.fMass = p->fMass;
					p->fMassStar = 0;
					(*nDeleted)++;
					pkdDeleteParticle(pkd, p);
					}
				else {
					starp.fMass = dInitStarMass;
					p->fMass -= dInitStarMass;
					p->fMassStar -= dInitStarMass;
					assert(p->fMass > 0);
					}

				starp.PdV = dtCoolingShutoff; /* Max local Cooling shutoff period */

				/*
				 * Save quantities -- as per old STARFORM
				 */
				for(j = 0; j < 3; j++) {
					starp.rForm[j] = starp.r[j];
					starp.vForm[j] = starp.v[j];
					}
				starp.u = T;
				starp.iGasOrder = starp.iOrder; /* iOrder gets reassigned in NewParticle() */

				starp.fTimeForm = dTime;
				starp.fBallMax = 0.0;
    
				/* NB: It is important that the star inherit special properties of the gas
				   particle such as being a target for movies or other tracing
				   Thus: Do not remove all the TYPE properties -- just gas specific ones */
				TYPEReset(&starp, TYPE_GAS);
				TYPESet(&starp, TYPE_STAR);

				/* Energy distribution */
				TYPESet(&starp, TYPE_SMOOTHACTIVE);

				(*nFormed)++;
				*dMassFormed += starp.fMass;
				
				pkdNewParticle(pkd, starp);    
				}
			}
		}
	}

#endif
