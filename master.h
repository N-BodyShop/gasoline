#ifndef MASTER_HINCLUDED
#define MASTER_HINCLUDED

#include "param.h"
#include "pst.h"
#include "mdl.h"
#include "parameters.h"
#include "floattype.h"

#define MSR_INIT_ECOSMO		1
#define MSR_STEP_ECOSMO		0

/*
 ** An integer marking the type of tree currently in use.
 ** MSR_TREE_NONE: undefined tree type
 ** MSR_TREE_SPATIAL: spatial binary tree
 ** MSR_TREE_DENSITY: density binary tree (the old style KD-tree!)
 ** MSR_TREE_QQ: perihelion-aphelion tree for planets.
 */
#define MSR_TREE_NONE		0
#define MSR_TREE_SPATIAL	1
#define MSR_TREE_DENSITY	2
#define MSR_TREE_QQ		3

typedef struct msrContext {
	PRM prm;
	PST pst;
	MDL mdl;
	LCL lcl;
	FLOAT fCenter[3];
	/*
	 ** Parameters.
	 */
	struct parameters param;	   
	/*
	 ** Other stuff...
	 */
	int nThreads;
	int N;
	int nDark;
	int nGas;
	int nStar;
        int nMaxOrder;		/* Order number of last particle */
	int nMaxOrderGas;
	int nMaxOrderDark;
	int iCurrMaxRung;
	int bOpenSpec;	/* was an opening parameter specified (used by +restart) */
	int iOpenType;
	double dCrit;
	/*
	 ** Comoving coordinate variables.
	 */
	double dEcosmo;
	double dUOld;
	double dTimeOld;
	/*
	 ** Redshift output points.
	 */
	int nMaxOuts;
	int nOuts;
	double *pdOutTime;
	int iOut;
	/*
	 ** Processor mapping for one-node-output functions.
	 */
	int *pMap;
	/*
	 ** An integer marking the type of tree currently in use.
	 ** MSR_TREE_NONE: undefined tree type
	 ** MSR_TREE_SPATIAL: spatial binary tree
	 ** MSR_TREE_DENSITY: density binary tree (the old style KD-tree!)
	 */
	int iTreeType;
	int bGravityTree;
        /*
         * Domain Decomposition Done
         */
        int bDoneDomainDecomp;
        int nActive;
        int nTreeActive;
        int nSmoothActive;
	} * MSR;



void msrInitialize(MSR *,MDL,int,char **);
void msrLogParams(MSR msr, FILE *fp);
int msrGetLock(MSR msr);
int msrCheckForStop(MSR msr);
void msrFinish(MSR);
int msrReadASCII(MSR, char *, int, double *);
double msrReadTipsy(MSR);
void msrWriteTipsy(MSR,char *,double);
void msrSetSoft(MSR msr,double);
void msrDomainDecomp(MSR);
void msrBuildTree(MSR,int,double,int);
void msrDomainColor(MSR);
void msrReorder(MSR);
void msrOutArray(MSR,char *,int);
void msrOutVector(MSR,char *,int);
void msrGetGasPressure(MSR);
void msrSmooth(MSR,double,int,int);
void msrReSmooth(MSR,double,int,int);
void msrGravity(MSR,double,int,int *,double *,double *,double *,int *);
void msrCalcE(MSR,int,double,double *,double *,double *,double *);
void msrDrift(MSR,double,double);
void msrKick(MSR,double,double);
double msrReadCheck(MSR,int *);
void msrWriteCheck(MSR,double,int);
int msrOutTime(MSR,double);
void msrReadOuts(MSR,double);
double msrMassCheck(MSR,double,char *);
void msrTopStepDKD(MSR msr, double dStep, double dTime, double dDelta, 
				double *pdMultiEff);
void msrTopStepKDK(MSR msr,
		   double dStep, /* Current step */
		   double dTime, /* Current time */
		   double dDelta, /* Time step */
		   int iRung,	/* Rung level */
		   int iKickRung, /* Gravity on all rungs from iRung
				     to iKickRung */
		   int iAdjust,	/* Do an adjust? */
		   double *pdActiveSum,
		   double *pdWMax,
		   double *pdIMax,
		   double *pdEMax,
		   int *piSec);

void msrRungStats(MSR);

/*------------------*/
/* Active Functions */
/*------------------*/
void msrActiveRung(MSR msr, int iRung, int bGreater);
void msrActiveOrder(MSR msr);

/* Deprecated functions */
/*
void msrSmoothActiveRung(MSR msr, int iRung, int bGreater);
void msrActiveGas(MSR msr);
void msrActiveDark(MSR msr);
void msrActiveStar(MSR msr);
void msrTreeActiveGas(MSR msr);
void msrTreeActiveDark(MSR msr);
void msrTreeActiveStar(MSR msr);
void msrSmoothActiveGas(MSR msr);
void msrSmoothActiveDark(MSR msr);
void msrSmoothActiveStar(MSR msr);
void msrActiveOrder(MSR msr);
void msrTreeActiveOrder(MSR msr);
#ifdef SUPERCOOL
void msrActiveCool(MSR msr);
void msrTreeActiveCool(MSR msr);
void msrSmoothActiveCool(MSR msr);
#endif
*/

/* Replacement functions */
void msrResetTouchRung(MSR msr, unsigned int iTestMask, unsigned int iSetMask);
void msrActiveExactType(MSR msr, unsigned int iFilterMask, unsigned int iTestMask, unsigned int iSetMask);
void msrActiveType(MSR msr, unsigned int iTestMask, unsigned int iSetMask);
void msrSetType(MSR msr, unsigned int iTestMask, unsigned int iSetMask);
void msrResetType(MSR msr, unsigned int iTestMask, unsigned int iSetMask);
int  msrCountType(MSR msr, unsigned int iFilterMask, unsigned int iTestMask);
void msrActiveMaskRung(MSR msr, unsigned int iSetMask, int iRung, int bGreater);
void msrActiveTypeRung(MSR msr, unsigned int iTestMask, unsigned int iSetMask, int iRung, int bGreater);
void msrActiveTypeOrder(MSR msr, unsigned int iTestMask );
/*------------------*/
/* Active Functions */
/*------------------*/

void msrVelocityRung(MSR msr, int iRung, double dDelta, double dTime,
		     int bAll);
void msrCoolVelocity(MSR,double,double);
void msrGrowMass(MSR msr, double dTime, double dDelta);
void msrCalcWriteStart(MSR);
void msrAddDelParticles(MSR msr);
void msrAccelStep(MSR msr, double dTime);
void msrInitDt(MSR msr);
void msrDtToRung(MSR msr, int iRung, double dDelta, int bAll);

/*
 ** Interface functions.
 */
int msrSteps(MSR);
char *msrOutName(MSR);
double msrDelta(MSR);
int msrLogInterval(MSR);
int msrCheckInterval(MSR);
int msrOutInterval(MSR);
int msrRestart(MSR);
int msrComove(MSR);
int msrKDK(MSR);
int msrDoSun(MSR);
double msrSoft(MSR);
int msrDoDensity(MSR);
int msrDoGravity(MSR msr);
int msrDoGas(MSR msr);
int msrFastGas(MSR msr);
void msrInitStep(MSR msr);
void msrSetRung(MSR msr, int iRung);
void msrInitAccel(MSR msr);
void msrSwitchTheta(MSR msr,double);
int msrMaxOrder(MSR msr);

void msrInitTimeSteps(MSR,double,double);

#ifdef GASOLINE
void msrInitSph(MSR,double);
int msrSphCurrRung(MSR msr, int iRung, int bGreater);
void msrSphStep(MSR msr, double dTime);
void msrSphViscosityLimiter(MSR msr, int bOn, double dTime);
void msrInitCooling(MSR msr);
#endif
#ifdef GLASS
void msrInitGlass(MSR);
#endif
#ifdef COLLISIONS
void msrFindRejects(MSR msr);
double msrReadSS(MSR msr);
void msrWriteSS(MSR msr, char *pszFileName, double dTime);
void msrCalcHill(MSR msr);
void msrHillStep(MSR msr);
void msrPlanetsKDK(MSR msr, double dStep, double dTime, double dDelta,
				   double *pdWMax, double *pdIMax, double *pdEMax, int *piSec);
void msrPlanetsDrift(MSR msr, double dStep, double dTime, double dDelta);
void msrFindEncounter(MSR msr, double dStart, double dEnd, double *dNext);
void msrMarkEncounters(MSR msr, double dTmax);
void msrLinearKDK(MSR msr, double dStep, double dTime, double dDelta);
void msrDoCollisions(MSR msr, double dTime, double dDelta);
#endif /* COLLISIONS */

#endif




