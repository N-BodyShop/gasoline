#ifdef STARFORM
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "pkd.h"
#include "starform.h"
#include "millerscalo.h"

#define max(A,B) ((A) > (B) ? (A) : (B))
#define min(A,B) ((A) < (B) ? (A) : (B))

/*
 * Star forming module for GASOLINE
 */

void stfmInitialize(STFM *pstfm)
{
    STFM stfm;
    
    stfm = (STFM) malloc(sizeof(struct stfmContext));
    assert(stfm != NULL);
    
    stfm->bInitialized = 0;
    stfm->dPhysDenMin = 0;
    stfm->dOverDenMin = 0;
    stfm->dTempMax = 0;
    stfm->dCStar = 0;
    stfm->dSecUnit = 0;
    stfm->dGmPerCcUnit = 0;
    stfm->dGmUnit = 0;
    stfm->dStarEff = 0.0;
    stfm->dMinGasMass = 0.0;
    stfm->dMinMassFrac = 0.0;
    stfm->dMaxStarMass = 0.0;
    *pstfm = stfm;
    }

/*
     Subroutine to initialize parameters dealing with star formation
     taken from TREESPH and modified greatly.
     Uses the following formula for the star formation rate:

              d(ln(rhostar))/dt=cstar/tdyn

*/

void stfmInitConstants(STFM stfm) 
{
    assert(stfm->dSecUnit != 0.0);
    assert(stfm->dGmPerCcUnit != 0.0);
    assert(stfm->dGmUnit != 0.0);
    
    stfm->dCStar = 0.1;
    stfm->dTempMax = 3.0e4;
    stfm->dSoftMin = 1.0;
    
    stfm->bInitialized = 1;
    }

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
    double E;
    PERBARYON Y;
    RATE Rate;
    double T;
    CL *cl = &(pkd->cl);
    double dExp = 1.0/(1.0 + cl->z);
    double dCosmoFac = dExp*dExp*dExp;
    PARTICLE starp;
    double dMprob;
    double dDeltaM;
    double l_jeans2;
    int small_jeans = 0;

    /*
     * Is particle in convergent part of flow?
     */
    
    if(p->divv >= 0.0)
	return;

    /*
     * Determine dynamical time.
     */
    tdyn = 1.0/sqrt(4.0*M_PI*p->fDensity/dCosmoFac);

    /*
     * Determine cooling time.
     */
    E = p->u * cl->dErgPerGmUnit;
    pkdPARTICLE2PERBARYON(&Y, p, cl->Y_H, cl->Y_He);
    T = clTemperature(Y.Total, E); 
    clRates(cl, &Rate, T);
    tcool = E/(clCoolTotal(cl, &Y, &Rate, p->fDensity*cl->dComovingGmPerCcUnit)
	       - clHeatTotal(cl, &Y)
	       - p->PdV*cl->dErgPerGmPerSecUnit);
    tcool /= cl->dSecUnit;

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
    if (l_jeans2 < p->fSoft*p->fSoft*stfm->dSoftMin*stfm->dSoftMin) 
        small_jeans = 1;

    if (!small_jeans && tsound <= tdyn)
        return;

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
    if(dMprob < (rand()/((double) RAND_MAX)))
	return;
    
    /*
     * Decrement mass of particle.
     */

    dDeltaM = p->fMass*stfm->dStarEff;
    p->fMass -= dDeltaM;
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
    starp.fBallMax = 0.0;
    TYPEReset(&starp, TYPE_GAS);
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
	}
    }

#endif
