#ifdef STARFORM
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "pkd.h"
#include "starform.h"
#include "millerscalo.h"

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
    stfm->dMaxStarMass = 0.0;
    
    *pstfm = stfm;
    }

/*
     Subroutine to initialize parameters dealing with star formation
     taken from TREESPH and modified greatly.
     Uses the following formula for the star formation rate:

              d(ln(rhostar))/dt=cstar/tdyn

*/

	/* mass of hydrogen atom in grams */
#define MHYDR 1.67e-24
	/* solar mass in grams */
#define MSOLG 1.99e33

void stfmInitConstants(STFM stfm) 
{
    const double dMRem=1.4;		/* mass of the supernovae
					   remnant in solar masses. */
    const double mSN=8.;	/* Mass above which stars supernova */
    double eff;			/* fraction of gas mass going into
				   FORMING stars */
    double dEsn=1.e51;		/* Energy of Supernovae in ergs. */
    const double dHDenMin=0.1; /* minimum density for star
				     formation in atoms/CC */
    const double dPSNTime=2.e7;	/* exponential decay time for
				   supernova in years */
    double dSNTime;		/* decay time in system units */
    const double dYield=0.02;	/* Fraction of metals returned to ISM
				   per unit mass of star formation */
    double dMgtMSN;		/* mass of stars greater than SN mass */
    double dMtot;		/* total mass of stars */
    double dNgtMSN;		/* number of stars greater than SN mass */
    MSPARAM MSparam;
    
    assert(stfm->dSecUnit != 0.0);
    assert(stfm->dGmPerCcUnit != 0.0);
    assert(stfm->dGmUnit != 0.0);

    stfm->dCStar = 0.1;
    stfm->dTempMax = 3.0e4;
    stfm->dSoftMin = 1.0e30;
    
    stfm->dPhysDenMin = dHDenMin/stfm->dGmPerCcUnit*MHYDR;
    dEsn /= stfm->dErgUnit;
    
    dSNTime = dPSNTime*3600*24*365.25/stfm->dSecUnit;
    stfm->dSNFrac = exp(-stfm->dDeltaT/dSNTime);
    
    MSInitialize(&MSparam);
    dNgtMSN = dMSCumNumber(MSparam, mSN);
    dMtot = dMSCumMass(MSparam, 0.0);
    dMgtMSN = dMSCumMass(MSparam, mSN);
    
    /* adjust for mass lost to SN, and mass gained from SN remnants */
    eff = stfm->dStarEff;
    stfm->dStarEff = eff*(1.0 - dMgtMSN/dMtot + dMRem/dMtot*dNgtMSN);

    stfm->dESNdM = eff*stfm->dGmUnit/MSOLG/dMtot*dNgtMSN*dEsn/(1.0 -
							     stfm->dStarEff);
    stfm->dESNdM *= (1.0/stfm->dSNFrac - 1.0)*stfm->dSNFrac;
    stfm->dYlddE = dYield/(stfm->dGmUnit/MSOLG/dMtot*dNgtMSN*dEsn);
    stfm->bInitialized = 1;
    free(MSparam);
    }

void stfmFormStars(STFM stfm, PKD pkd, PARTICLE *p, double dTime)
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
    PARTICLE starp;
    double dMprob;
    double dDeltaM;
    
    /*
     * Adjust supernovea energy.
     */
    p->fSNRes *= stfm->dSNFrac;
    if(p->fSNRes < 1.0e-4*stfm->dESNdM)
	p->fSNRes = 0.0;

    /*
     * Is particle in convergent part of flow?
     */
    if(p->divv >= 0.0)
	return;

    /*
     * Determine dynamical time.
     */
    tdyn = 1.0/sqrt(4.0*M_PI*dExp*dExp*dExp/p->fDensity);

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

    if(tcool < 0.0) return;
    /*
     * Determine sound crossing time.
     */
    tsound = sqrt(0.25*p->fBall2)/p->c;
    
#if 0
    if(p->fBall2 <= stfm->dSoftMin) /* XXX There should probably be a
				       criteria such that stars form
				       if the Jean's length is less
				       than the softening.  C.F.
				       the code in
				       tipsy:starform_func() */
	too_small = 1;
#endif
    /*
     * Determine if this particle satisfies all conditions.
     */
    if(!too_small && tsound <= tdyn)
	return;

    
    if(p->fDensity < stfm->dOverDenMin ||
       p->fDensity*cl->dComovingGmPerCcUnit/MHYDR < stfm->dPhysDenMin)
	return;

    if(tcool < 0.0 || tdyn > tcool || T < stfm->dTempMax)
	tform = tdyn;
    else
	tform = tcool;
    
    dMprob = 1.0 - exp(-stfm->dCStar*stfm->dDeltaT/tform);
    if(dMprob > (rand()/((double) RAND_MAX)))
	return;
    
    /*
     * Decrement mass of particle.
     */
    dDeltaM = p->fGasMass*stfm->dStarEff;
    p->fGasMass -= dDeltaM;
    p->fSNRes /= 1.0 - stfm->dStarEff;
    p->fSNRes += stfm->dESNdM*dDeltaM;
    /*
     * XXX also adjust metals and age.
     */
    /* 
     * XXX do this deletion later.
     */
    if(p->fGasMass < stfm->dMinGasMass) {
	pkdDeleteParticle(pkd, p - pkd->pStore);
	}
    /*
     * Calculate metalicity, and add into star part of metalicity.
     * Include Fe?
     */
    /*
     * form star
     */
    if(p->fGasMass < stfm->dMinGasMass
       || p->fMass - p->fGasMass >  stfm->dMaxStarMass) {
	starp = *p;
	if(p->fGasMass < stfm->dMinGasMass)
	    starp.fMass = p->fMass;
	else
	    starp.fMass = p->fMass - p->fGasMass;
	p->fMass = p->fGasMass;	/* reduce mass of gas particle */
	starp.fTimeForm = dTime;
	pkdNewParticle(pkd, starp);
	}
    }

void pkdFormStars(PKD pkd, STFM stfm, double dTime)
{
    int i;
    PARTICLE *p;
    int n = pkdLocal(pkd);
    
    for(i = 0; i < n; ++i) {
        p = &pkd->pStore[i];
	if(pkdIsGas(pkd, p))
	    stfmFormStars(stfm, pkd, p, dTime);
	}
    }

#endif
