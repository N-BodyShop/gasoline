#ifdef STARFORM
#ifndef STARFORM_HINCLUDED
#define STARFORM_HINCLUDED

typedef struct stfmContext 
{
    double dMinMassFrac;        /* minimum fraction of average mass of neighbour 
                                   particles required for particle to avoid deletion */
    double dDeltaT;		/* timestep in system units */
    double dGmUnit;		/* system mass in grams */
    double dGmPerCcUnit;	/* system density in gm/cc */
    double dSecUnit;		/* system time in seconds */
    double dErgUnit;		/* system energy in ergs */
    double dPhysDenMin;		/* Physical density minimum for star
				   formation (in system units) */
    double dOverDenMin;		/* Overdensity minimum for star formation */
    double dTempMax;		/* Form stars below this temperature
				   EVEN IF the gas is not cooling. */
    double dSoftMin;		/* Jean's length as a fraction of
				   softening at which to form stars*/
    double dCStar;		/* Star formation constant */
    double dStarEff;		/* Fraction of gas mass converted into
				 star mass per timestep. */
	double dInitStarMass;    /* Fixed Initial Star Mass */
    double dMinGasMass;		/* minimum mass gas before we delete
				   the particle. */
    double dMaxStarMass;	/* maximum mass star particle to form */

} * STFM;
	
void stfmInitialize(STFM *pstfm);


void stfmFormStars(STFM stfm, PKD pkd, PARTICLE *p,
		   double dTime, /* current time */
		   int *nFormed, /* number of stars formed */
		   double *dMassFormed,	/* mass of stars formed */
		   int *nDeleted); /* gas particles deleted */
void pkdFormStars(PKD pkd, STFM stfm, double dTime,
		   int *nFormed, /* number of stars formed */
		   double *dMassFormed,	/* mass of stars formed */
		   int *nDeleted); /* gas particles deleted */

#endif
#endif

#ifdef SIMPLESF
void pkdSimpleStarForm(PKD pkd, double dRateCoeff, double dTMax, double dDenMin, double dDelta, double dTime,
					   double dInitStarMass, double dESNPerStarMass, double dtCoolingShutoff, int bdivv,
                                           int *nFormed, /* number of stars formed */
                                           double *dMassFormed, /* mass of stars formed */
                                           int *nDeleted); /* gas particles deleted */
#endif


