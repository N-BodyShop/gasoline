#ifndef STARFORM_HINCLUDED
#define STARFORM_HINCLUDED

typedef struct stfmContext 
{
    double dDeltaT;		/* timestep in system units */
    double dGmUnit;		/* system mass in grams */
    double dGmPerCcUnit;	/* system density in gm/cc */
    double dSecUnit;		/* system time in seconds */
    double dErgUnit;		/* system energy in ergs */
    double dPhysDenMin;		/* Physical density minimum for star
				   formation (in system units) */
    double dOverDenMin;		/* Overdensity minimum for star formation */
    double dTempMax;		/* Use dynamical time when T is below
				   this for forming stars. */
    double dSoftMin;		/* */
    double dCStar;		/* Star formation constant */
    double dStarEff;		/* Fraction of gas mass converted into
				 star mass per timestep. */
    double dMinGasMass;		/* minimum mass gas before we delete
				   the particle. */
    double dMaxStarMass;	/* maximum mass star particle to form */
    double dESNdM;		/* SN energy released per mass of
				   stars formed. */
    double dSNFrac;		/* fraction of star mass that goes
				   supernovae per timestep. */
    double dYlddE;		/* Metal yield per energy */
    
    int bInitialized;
} * STFM;
	
void stfmInitialize(STFM *pstfm);

void stfmInitConstants(STFM stfm);

void stfmFormStars(STFM stfm, PKD pkd, PARTICLE *p, double dTime);
void pkdFormStars(PKD pkd, STFM stfm, double dTime);

#endif
