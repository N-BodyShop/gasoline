
#ifdef STARFORM
#ifndef STARFORM_HINCLUDED
#define STARFORM_HINCLUDED

typedef struct stfmContext 
    {
    double dMinMassFrac;        /* minimum fraction of average mass of neighbour 
                                   particles required for particle to avoid deletion */
    double dDeltaT;		/* timestep in system units */
    double dExp;        /* a */
    double dCosmoFac;   /* a^3 */
    double dGmUnit;		/* system mass in grams */
    double dGmPerCcUnit;	/* system density in gm/cc */
    double dSecUnit;		/* system time in seconds */
    double dErgUnit;		/* system energy in ergs */
    double dPhysDenMin;		/* Physical density minimum for star
				   formation (in system units) */
    double dOverDenMin;		/* Overdensity minimum for star formation */
    int bTempInclHot;   /* Include uHot in temp estimate for TempMax */
    double dTempMax;		/* Form stars below this temperature
				   EVEN IF the gas is not cooling. */
    double dSoftMin;		/* Jean's length as a fraction of
				   softening at which to form stars*/
    double dCStar;		/* Star formation constant */
    double dStarEff;		/* Fraction of gas mass converted into
				 star mass per timestep. */
    double dInitStarMass;       /* Fixed Initial Star Mass */
    double dMinSpawnStarMass;   /* Minimum Initial Star Mass */
    double dMinGasMass;		/* minimum mass gas before we delete
				   the particle. */
    double dMaxGasMass;		/* maxmimum mass gas so that it gets
				   no more feedback. */
    double dMaxStarMass;	/* maximum mass star particle to form */
    double dZAMSDelayTime;      /* Time Delay between star particle formation and ZAMS time of stars in it */
    double dBHFormProb;         /* Probability star will become a BH */
    int bBHForm;                /* are BH seeds allowed to form */ 
    double dInitBHMass;         /* Initial BH mass */
    double dStarClusterMass;   /* Target Star Cluster Mass */
    double dStarClusterRatio;   /* Target Star Cluster Mass: default 0.5*3/5 */
#ifdef COOLING_MOLECULARH
  double dStarFormEfficiencyH2; /* Star formation efficiency, CStar, is multiplied by dStarFormEfficiencyH2 times the fraction of hydrogen in molecular form  */
#endif

} * STFM;
	
void stfmInitialize(STFM *pstfm);

void pkdStarLogInit(PKD pkd);
void pkdStarLogFlush(PKD pkd, char *pszFileName);

double stfmFormStarProb(STFM stfm, PKD pkd, PARTICLE *p,
    double dTime /* current time */);
void stfmFormStarParticle(STFM stfm, PKD pkd, PARTICLE *p,
    double dDeltaM, /* star mass */
    double dTime, /* current time */
    int *nFormed, /* number of stars formed */
    double *dMassFormed,	/* mass of stars formed */
    int *nDeleted); /* gas particles deleted */
void stfmFormStars(STFM stfm, PKD pkd, PARTICLE *p,
    double dTime, /* current time */
    int *nFormed, /* number of stars formed */
    double *dMassFormed,	/* mass of stars formed */
    int *nDeleted); /* gas particles deleted */
void pkdFormStars(PKD pkd, STFM stfm, double dTime,
    int *nFormed, /* number of stars formed */
    double *dMassFormed,	/* mass of stars formed */
    int *nDeleted); /* gas particles deleted */

struct inStarClusterFormPrecondition
{
    double dTime;
    struct stfmContext stfm;
    };
struct outStarClusterFormPrecondition 
{
    int n;
    };
void pkdStarClusterFormPrecondition(PKD pkd, struct inStarClusterFormPrecondition in);

#endif
#endif

#ifdef SIMPLESF
void pkdSimpleStarForm(PKD pkd, double dRateCoeff, double dTMax, double dDenMin, double dDelta, double dTime,
					   double dInitStarMass, double dESNPerStarMass, double dtCoolingShutoff, int bdivv,
                                           int *nFormed, /* number of stars formed */
                                           double *dMassFormed, /* mass of stars formed */
                                           int *nDeleted); /* gas particles deleted */
#endif


