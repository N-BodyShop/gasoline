#ifndef PKD_HINCLUDED
#define PKD_HINCLUDED

#include <sys/resource.h>
#include "mdl.h"
#include "floattype.h"

#ifdef GASOLINE
#include "cooling.h"
#endif

/*
 ** The following sort of definition should really be in a global
 ** configuration header file -- someday...
 */

#if defined(GASOLINE) || defined(ROT_FRAME) || defined(SLIDING_PATCH)
#define NEED_VPRED
#endif

#define CID_TOP			0
#define CID_PARTICLE	0
#define CID_CELL		1

#define ROOT		1
#define LOWER(i)	(i<<1)
#define UPPER(i)	((i<<1)+1)
#define PARENT(i)	(i>>1)
#define SIBLING(i) 	((i&1)?i-1:i+1)
#define SETNEXT(i)\
{\
	while (i&1) i=i>>1;\
	++i;\
	}

#define MAX_TIMERS		10

typedef struct particle {
	int iOrder;
	unsigned int iActive;  
	int iRung;
	int cpStart;
	FLOAT fWeight;
	FLOAT fMass;
	FLOAT fSoft;
#ifdef CHANGESOFT
	FLOAT fSoft0;
#endif
	FLOAT r[3];
	FLOAT v[3];
	FLOAT a[3];
	FLOAT fPot;
	FLOAT fBall2;
	FLOAT fDensity;
	FLOAT dt;			/* a time step suggestion */
	FLOAT dtGrav;		/* suggested 1/dt^2 from gravity */
#ifdef SUPERCOOL
	FLOAT vMean[3];
#endif
#ifdef COLORCODE
	FLOAT fColor;
#endif
	FLOAT fBallMax;		/* SPH 2h Max value */
#ifdef GASOLINE
	FLOAT uPred;		/* predicted thermal energy */
	FLOAT PoverRho2;	/* P/rho^2 */
	FLOAT u;	        /* thermal energy */ 
	FLOAT c;			/* sound speed */
	FLOAT mumax;		/* sound speed like viscosity term */
	FLOAT PdV;	        /* P dV heating (includes shocking) */
#ifdef PDVDEBUG
	FLOAT PdVvisc;		/* P dV from shock (testing) */
	FLOAT PdVpres;		/* P dV from adiabatic compression (testing) */
#endif
	FLOAT divv;		/* Balsara viscosity reduction */
	FLOAT curlv[3];         
	FLOAT BalsaraSwitch;
#ifdef SHOCKTRACK
        FLOAT aPres[3];
        FLOAT ShockTracker;     /* Shock tracker */
        FLOAT divrhov;          /* debug */
        FLOAT gradrho[3];       /* debug */
#endif
/*	FLOAT fDensSave;*/	/* Used by diagnostic DensCheck funcs */
#ifndef NOCOOLING
	FLOAT uDot;			/* Rate of change of u -- for predicting u */
	COOLPARTICLE CoolParticle;  /* Abundances and any other cooling internal variables */
#endif
	FLOAT fMetals;
	FLOAT fTimeForm;
#ifdef SIMPLESF
    FLOAT fMassStar;
	FLOAT fESN;
	FLOAT rForm[3];		/* record pos and vel of star formation */
	FLOAT vForm[3];
	int iGasOrder;		/* gas from which star formed */
#endif
#ifdef STARFORM
	FLOAT fESNrate;
	FLOAT fMSN;
	FLOAT fMSNII;
	FLOAT fSNMetals;
        FLOAT fTimeCoolIsOffUntil;
	FLOAT rForm[3];		/* record pos and vel of star formation */
	FLOAT vForm[3];
	FLOAT fMassForm;	/* record original mass of star */
	int iGasOrder;		/* gas from which star formed */
#endif
#endif  /* GASOLINE */
#ifdef COLLISIONS
	int iOrgIdx;		/* for tracking of mergers, aggregates etc. */
	FLOAT w[3];			/* spin vector */
	int iColor;			/* handy color tag */
	int iDriftType;		/* either NORMAL or KEPLER */
	double dtCol;		/* time to next encounter or collision */
	int iOrderCol;		/* neighbour or collider iOrder */
	double dtPrevCol;	/* time of previous collision */
	int iPrevCol;		/* iOrder of previous collider */
	int bTinyStep;		/* flag for imminent collapse */
#endif
#ifdef SAND_PILE
	int bStuck;
#endif
#ifdef NEED_VPRED
	FLOAT vPred[3];		/* predicted velocity (time centered) */
#endif
#ifdef AGGS
	/* Position in principal frame of the aggregate (normally).  We
	 *  temporarily store the COM frame position here during the
	 *  process of computing the aggregate parameters. */
	FLOAT r_agg[3];
#endif
	} PARTICLE;

/* Active Type Masks */

/* Active: -- eg. Calculate new acceleration, PdV, etc... for this particle */
#define TYPE_ACTIVE            (1<<0)
/* In the Tree: */
#define TYPE_TREEACTIVE        (1<<1)
/* Gather to/Scatter from this particle with in smooths: */
#define TYPE_SMOOTHACTIVE      (1<<2)
/* Smooth has processed this particle */
#define TYPE_SMOOTHDONE        (1<<3)

/* Types used for Fast Density only (so far) */
/* Sum Fast Density on this particle */
#define TYPE_DensACTIVE        (1<<4)
/* Neighbour of ACTIVE (incl. ACTIVE): */
#define TYPE_NbrOfACTIVE       (1<<5)
/* Potential Scatter Neighbour */
#define TYPE_Scatter           (1<<6)
/* Density set to zero already */
#define TYPE_DensZeroed        (1<<7)

/* Particle Type Masks */

#define TYPE_GAS               (1<<8)
#define TYPE_DARK              (1<<9)
#define TYPE_STAR              (1<<10)
#define TYPE_SUPERCOOL         (1<<11)

/* Particle marked for deletion.  Will be deleted in next
   msrAddDelParticles(); */
#define TYPE_DELETED           (1<<12)

#define TYPE_PHOTOGENIC        (1<<13)

/* Combination Masks */
#define TYPE_ALLACTIVE			(TYPE_ACTIVE|TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE)
#define TYPE_ALL				(TYPE_GAS|TYPE_DARK|TYPE_STAR)

/* Type Macros */
int TYPEQueryACTIVE      ( PARTICLE *a );
int TYPEQueryTREEACTIVE  ( PARTICLE *a );
int TYPEQuerySMOOTHACTIVE( PARTICLE *a );
int TYPETest  ( PARTICLE *a, unsigned int mask );
int TYPEFilter( PARTICLE *a, unsigned int filter, unsigned int mask );
int TYPESet   ( PARTICLE *a, unsigned int mask );
int TYPEReset ( PARTICLE *a, unsigned int mask );
/* This retains Particle Type and clears all flags: */
int TYPEClearACTIVE( PARTICLE *a ); 

/* Warning: This erases Particle Type */
int TYPEClear( PARTICLE *a ); 

#define TYPEQueryACTIVE(a)       ((a)->iActive & TYPE_ACTIVE)
#define TYPEQueryTREEACTIVE(a)   ((a)->iActive & TYPE_TREEACTIVE)
#define TYPEQuerySMOOTHACTIVE(a) ((a)->iActive & TYPE_SMOOTHACTIVE)
#define TYPETest(a,b)            ((a)->iActive & (b))
#define TYPEFilter(a,b,c)        (((a)->iActive & (b))==(c))
#define TYPESet(a,b)             ((a)->iActive |= (b))
#define TYPEReset(a,b)           ((a)->iActive &= (~(b)))
#define TYPEClearACTIVE(a)       ((a)->iActive &= (TYPE_ALL|TYPE_SUPERCOOL))
#define TYPEClear(a)             ((a)->iActive = 0)


#define CHECKPOINT_VERSION 7

typedef struct chkParticle {
	int iOrder;
	FLOAT fMass;
	FLOAT fSoft;
	FLOAT r[3];
	FLOAT v[3];
#ifdef GASOLINE
	FLOAT u;
	FLOAT fMetals;
#ifndef NOCOOLING
	COOLPARTICLE CoolParticle;
#endif
#ifdef STARFORM
        FLOAT fTimeCoolIsOffUntil;
	FLOAT fTimeForm;
	FLOAT rForm[3];		/* record pos and vel of star formation */
	FLOAT vForm[3];
	FLOAT fMassForm;	/* record original mass of star */
	FLOAT fDenForm;
	int iGasOrder;
#endif
#ifdef SIMPLESF
	FLOAT fMassStar;
	FLOAT fTimeForm;
	FLOAT rForm[3];		/* record pos and vel of star formation */
	FLOAT vForm[3];
	FLOAT fDenForm;
	int iGasOrder;
#endif
#endif
#ifdef COLLISIONS
	int iOrgIdx; /* added for version 7 */
	FLOAT w[3];
	int iColor;
#endif /* COLLISIONS */
	} CHKPART;

typedef struct bndBound {
	FLOAT fMin[3];
	FLOAT fMax[3];
	} BND;

struct pkdCalcCellStruct {
	double Qxx,Qyy,Qzz,Qxy,Qxz,Qyz;
	/*
	 ** Reduced multipole moments for l>2 !!!
	 */
	double Oxxx,Oxyy,Oxxy,Oyyy,Oxxz,Oyyz,Oxyz;
	double Oxzz, Oyzz, Ozzz;
	double Hxxxx,Hxyyy,Hxxxy,Hyyyy,Hxxxz,Hyyyz,Hxxyy,Hxxyz,Hxyyz;
	double Hxxzz, Hxyzz, Hxzzz, Hyyzz, Hyzzz, Hzzzz;
	double Bmax,B2,B3,B4,B5,B6;
	};

typedef struct kdNode {
	int iDim;
	double fSplit;
	BND bnd;
	BND bndBall;	/* Bound including fBall*(1+changemax) */
	int pLower;		/* also doubles as thread id for the LTT */
	int pUpper;		/* pUpper < 0 indicates no particles in tree! */
	int iLower;
	int iUpper;
	double fMass;
	double fSoft;
	FLOAT r[3];
	struct pkdCalcCellStruct mom;
	double fOpen2;
	} KDN;

typedef struct ilPart {
	double m,h;
	double x,y,z;
	} ILP;

typedef struct ilCellSoft {
	double m,h;
	double x,y,z;
	double xx,yy,zz,xy,xz,yz;
	} ILCS;

/*
 ** moment tensor components.
 */
typedef struct ilCellNewt {
	double m;
	double x,y,z;
	double xx,yy,xy,xz,yz;
	double zz;
	double xxx,xyy,xxy,yyy,xxz,yyz,xyz;
	double xzz,yzz,zzz;
	double xxxx,xyyy,xxxy,yyyy,xxxz,yyyz,xxyy,xxyz,xyyz;
	double xxzz,xyzz,xzzz,yyzz,yzzz,zzzz;
	} ILCN;

/* IBM brain damage */
#undef hz

typedef struct ewaldTable {
	double hx,hy,hz;
	double hCfac,hSfac;
	} EWT;

typedef struct pkdContext {
	MDL mdl;
	int idSelf;
	int nThreads;
	int nStore;
	int nRejects;
	int nLocal;
	int nActive;
	int nTreeActive;
	int nSmoothActive;
	int nDark;
	int nGas;
	int nStar;
	int nMaxOrderDark;
	int nMaxOrderGas;
	int nBucket;
	int nLevels;
	int nSplit;
	int nNodes;
	int iExtraBucket;
	int iOrder;
	int iFreeCell;
	int iRoot;
	FLOAT fPeriod[3];
	int *piLeaf;
	KDN *kdTop;
	KDN *kdNodes;
	PARTICLE *pStore;
	/*
	 ** gravitational interaction lists
	 */
	int nMaxPart;
	int nMaxCellSoft;
	int nMaxCellNewt;
	int nPart;
	int nCellSoft;
	int nCellNewt;
	int nSqrtTmp;
	ILP *ilp;
	ILCS *ilcs;
	ILCN *ilcn;
	double *sqrttmp;
	double *d2a;
	/*
	 ** Ewald summation setup.
	 */
	ILCN ilcnRoot;
	int nMaxEwhLoop;
	int nEwhLoop;
	EWT *ewt;
	/*
	 ** Timers stuff.
	 */
	struct timer {
		double sec;
		double stamp;
		double system_sec;
		double system_stamp;
		double wallclock_sec;
		double wallclock_stamp;
		int iActive;
		} ti[MAX_TIMERS];
#ifdef GASOLINE
	/* 
	 ** Cooling 
	 */
#ifndef NOCOOLING
	COOL *Cool;
#endif
#endif
#ifdef SLIDING_PATCH
	/*
	 ** Rotating frame info...
	 */
	double dOrbFreq;
	double dTime;
#endif
	} * PKD;

int pkdIsGas(PKD,PARTICLE *);
#define pkdIsGas( pkd, pTMP) TYPETest( (pTMP), TYPE_GAS )
int pkdIsDark(PKD,PARTICLE *);
#define pkdIsDark( pkd, pTMP) TYPETest( (pTMP), TYPE_DARK )
int pkdIsStar(PKD,PARTICLE *);
#define pkdIsStar( pkd, pTMP) TYPETest( (pTMP), TYPE_STAR )

int pkdIsGasByOrder(PKD pkd,PARTICLE *p);
int pkdIsDarkByOrder(PKD pkd,PARTICLE *p);
int pkdIsStarByOrder(PKD pkd,PARTICLE *p);


typedef struct CacheStatistics {
	double dpNumAccess;
	double dpMissRatio;
	double dpCollRatio;
	double dpMinRatio;
	double dcNumAccess;
	double dcMissRatio;
	double dcCollRatio;
	double dcMinRatio;
	} CASTAT;

/* JW: */

#define GASMODEL_UNSET -1
enum GasModel {
	GASMODEL_ADIABATIC, 
	GASMODEL_ISOTHERMAL, 
	GASMODEL_COOLING, 
	GASMODEL_GLASS
	}; 

#define PKD_ORDERTEMP	256

#define pkdRoot(iCell,id)\
{\
	id = -1;\
	iCell = ROOT;\
	}

#define pkdIsRoot(iCell,id)		((id==-1)?((iCell==ROOT)?1:0):0)

/*
 * There is now a slight inefficency here when going from the top tree to
 * a node tree in that we visit the root cell twice (once as the leaf of the
 * top tree and once as the root of the node tree).  This is necessary to
 * check if the root cell is a bucket.
 */

#define pkdLower(pkd,iCell,id)\
{\
	if (id == -1) {\
		id = pkd->kdTop[iCell].pLower;\
		if (id != -1) iCell = ROOT;\
		else iCell = LOWER(iCell);\
		}\
	else iCell = LOWER(iCell);\
	}

#define pkdUpper(pkd,iCell,id)\
{\
	if (id == -1) {\
		id = pkd->kdTop[iCell].pLower;\
		if (id != -1) iCell = ROOT;\
		else iCell = UPPER(iCell);\
		}\
	else iCell = UPPER(iCell);\
	}

#define pkdParent(pkd,iCell,id)\
{\
	iCell = PARENT(iCell);\
	if (iCell == ROOT) {\
		if (id != -1) {\
			iCell = pkd->piLeaf[id];\
			id = -1;\
			}\
		}\
	}

#define pkdNext(pkd,iCell,id)\
{\
	SETNEXT(iCell);\
	if (iCell == ROOT) {\
		if (id != -1) {\
			iCell = pkd->piLeaf[id];\
			id = -1;\
			SETNEXT(iCell);\
			}\
		}\
	}

double pkdGetTimer(PKD,int);
double pkdGetSystemTimer(PKD,int);
double pkdGetWallClockTimer(PKD,int);
void pkdClearTimer(PKD,int);
void pkdStartTimer(PKD,int);
void pkdStopTimer(PKD,int);
void pkdInitialize(PKD *,MDL,int,int,int,FLOAT *,int,int,int);
void pkdFinish(PKD);
void pkdReadTipsy(PKD,char *,int,int,int,int,double,double);
void pkdSetSoft(PKD pkd,double dSoft);
#ifdef CHANGESOFT
void pkdPhysicalSoft(PKD pkd,double, double, int);
void pkdPreVariableSoft(PKD pkd);
void pkdPostVariableSoft(PKD pkd,double dSoftMax,int bSoftMaxMul);
#endif
void pkdCalcBound(PKD,BND *,BND *,BND *,BND *);
void pkdGasWeight(PKD);
void pkdRungDDWeight(PKD, int, double);
int pkdWeight(PKD,int,FLOAT,int,int,int,int *,int *,FLOAT *,FLOAT *);
int pkdLowerPart(PKD,int,FLOAT,int,int);
int pkdUpperPart(PKD,int,FLOAT,int,int);
int pkdWeightWrap(PKD,int,FLOAT,FLOAT,int,int,int,int *,int *,FLOAT *,FLOAT *);
int pkdLowerPartWrap(PKD,int,FLOAT,FLOAT,int,int);
int pkdUpperPartWrap(PKD,int,FLOAT,FLOAT,int,int);
int pkdLowerOrdPart(PKD,int,int,int);
int pkdUpperOrdPart(PKD,int,int,int);
int pkdActiveTypeOrder(PKD, unsigned int);
int pkdActiveOrder(PKD);
int pkdColRejects(PKD,int,FLOAT,FLOAT,int);
int pkdSwapRejects(PKD,int);
int pkdSwapSpace(PKD);
int pkdFreeStore(PKD);
int pkdLocal(PKD);
int pkdActive(PKD);
int pkdTreeActive(PKD);
int pkdInactive(PKD);
int pkdTreeInactive(PKD);
int pkdNodes(PKD);
void pkdDomainColor(PKD);
int pkdColOrdRejects(PKD,int,int);
void pkdLocalOrder(PKD);
void pkdWriteTipsy(PKD,char *,int,int,double,double,int);
void pkdCombine(KDN *,KDN *,KDN *);
void pkdCalcCell(PKD,KDN *,FLOAT *,int,struct pkdCalcCellStruct *);
double pkdCalcOpen(KDN *,int,double,int);
void pkdBuildLocal(PKD,int,int,double,int,int,int,KDN *);
void pkdBuildBinary(PKD,int,int,double,int,int,int,KDN *);
void pkdThreadTree(PKD pkd,int iCell,int iNext);
void pkdGravAll(PKD,int,int,int,int,int,double,double,int,double *,int *,
				double *,double *,double *,CASTAT *,double *);
void pkdCalcEandL(PKD,double *,double *,double *,double []);
void pkdCalcEandLExt(PKD,double *,double[],double [],double *);
void pkdDrift(PKD,double,FLOAT *,int,int,FLOAT);
void pkdUpdateUdot(PKD pkd,double,double,double,int,int);
void pkdKick(PKD pkd,double,double, double, double, double, double, int, double, double);
void pkdReadCheck(PKD,char *,int,int,int,int);
void pkdWriteCheck(PKD,char *,int,int);
void pkdDistribCells(PKD,int,KDN *);
void pkdCalcRoot(PKD,struct ilCellNewt *);
void pkdDistribRoot(PKD,struct ilCellNewt *);
void pkdSwapAll(PKD pkd, int idSwap);
double pkdMassCheck(PKD pkd);
void pkdMassMetalsEnergyCheck(PKD pkd, double *dTotMass, 
                    double *dTotMetals, double *dTotEnergy);
void pkdSetRung(PKD pkd, int iRung);
void pkdBallMax(PKD pkd, int iRung, int bGreater, double ddHonHLimit);
int pkdActiveRung(PKD pkd, int iRung, int bGreater);
int pkdCurrRung(PKD pkd, int iRung);
void pkdGravStep(PKD pkd, double dEta);
void pkdAccelStep(PKD pkd, double dEta, double dVelFac, double
				  dAccFac, int bDoGravity, int bEpsAcc, int bSqrtPhi, double dhMinOverSoft);
void pkdDensityStep(PKD pkd, double dEta, double dRhoFac);
int pkdDtToRung(PKD pkd, int iRung, double dDelta, int iMaxRung, int bAll,
		int *pnMaxRung, int *piMaxRungIdeal);
void pkdInitDt(PKD pkd, double dDelta);
int pkdRungParticles(PKD,int);
void pkdCoolVelocity(PKD,int,double,double,double);
void pkdGrowMass(PKD pkd,int nGrowMass, double dDeltaM);
void pkdInitAccel(PKD);
int pkdOrdWeight(PKD,int,int,int,int,int *,int *);
void pkdUnDeleteParticle(PKD pkd, PARTICLE *p);
void pkdDeleteParticle(PKD pkd, PARTICLE *p);
void pkdNewParticle(PKD pkd, PARTICLE p);
int pkdResetTouchRung(PKD pkd, unsigned int iTestMask, unsigned int iSetMask);
int pkdActiveExactType(PKD pkd, unsigned int iFilterMask, unsigned int iTestMask, unsigned int iSetMask);
int pkdActiveType(PKD pkd, unsigned int iTestMask, unsigned int iSetMask);
int pkdSetType(PKD pkd, unsigned int iTestMask, unsigned int iSetMask);
int pkdResetType(PKD pkd, unsigned int iTestMask, unsigned int iSetMask);
int pkdCountType(PKD pkd, unsigned int iFilterMask, unsigned int iTestMask);
int pkdActiveMaskRung(PKD pkd, unsigned int iSetMask, int iRung, int bGreater );
int pkdActiveTypeRung(PKD pkd, unsigned int iTestMask, unsigned int iSetMask, int iRung, int bGreater);
int pkdSetTypeFromFile(PKD pkd, int iSetMask, char *file, int *niOrder, int *nSet);

void pkdSetParticleTypes(PKD pkd, int nSuperCool);
void pkdColNParts(PKD pkd, int *pnNew, int *nDeltaGas, int *nDeltaDark,
				  int *nDeltaStar);
void pkdNewOrder(PKD pkd, int nStart);

struct outGetNParts { 
	int n;
    int nGas;
    int nDark;
    int nStar;
    int iMaxOrderGas;
    int iMaxOrderDark;
    int iMaxOrderStar;
    };

void pkdGetNParts(PKD pkd, struct outGetNParts *out );
void pkdSetNParts(PKD pkd, int nGas, int nDark, int nStar, int, int nMaxOrderGas,
				  int nMaxOrderDark);
void pkdSunIndirect(PKD,double *,int,double);
void pkdLogHalo(PKD);
void pkdHernquistSpheroid(PKD pkd);
void pkdNFWSpheroid(PKD pkd);
void pkdHomogSpheroid(PKD pkd);
void pkdBodyForce(PKD pkd);
void pkdMiyamotoDisk(PKD pkd);
void pkdTimeVarying(PKD pkd,double dTime);
#ifdef ROT_FRAME
void pkdRotFrame(PKD pkd, double dOmega, double dOmegaDot);
#endif

#ifdef GASOLINE

void pkdUpdateuDot(PKD,double,double,double,int,int);
void pkdUpdateShockTracker(PKD,double, double, double);
void pkdAdiabaticGasPressure(PKD, double gammam1, double gamma);
void pkdCoolingGasPressure(PKD, double gammam1, double gamma);
void pkdLowerSoundSpeed(PKD, double);
void pkdInitEnergy(PKD pkd, double dTuFac, double z, double dTime );
void pkdKickRhopred(PKD pkd, double dHubbFac, double dDelta);
int pkdSphCurrRung(PKD pkd, int iRung, int bGreater);
void pkdSphStep(PKD pkd, double dCosmoFac, double dEtaCourant, double dEtauDot, int bViscosityLimitdt);
void pkdSphViscosityLimiter(PKD pkd, int bOn, int bShockTracker);

void pkdDensCheck(PKD pkd, int iRung, int bGreater, int iMeasure, void *data);

#endif /* GASOLINE */
#ifdef GLASS
void pkdGlassGasPressure(PKD, void *in);
void pkdRandomVelocities(PKD pkd, double dMaxVL, double dMaxVR);
#endif

#ifdef COLLISIONS
int pkdNumRejects(PKD pkd);
void pkdReadSS(PKD pkd, char *pszFileName, int nStart, int nLocal);
void pkdWriteSS(PKD pkd, char *pszFileName, int nStart);
void pkdKickUnifGrav(PKD pkd, double dvx, double dvy, double dvz);
void pkdNextEncounter(PKD pkd, double *dt);
void pkdMarkEncounters(PKD pkd, double dt);
#ifdef SIMPLE_GAS_DRAG
void pkdSimpleGasDrag(PKD pkd,int iFlowOpt,int bEpstein,double dGamma,
					  double dOmegaZ);
#endif
#ifdef OLD_KEPLER/*DEBUG*/
int pkdLowerQQPart(PKD pkd, int d, FLOAT fSplit, int i, int j);
int pkdUpperQQPart(PKD pkd, int d, FLOAT fSplit, int i, int j);
void pkdQQCalcBound(PKD pkd, BND *pbnd, BND *pbndActive);
void pkdQQBuild(PKD pkd, int nBucket, int bActiveOnly, KDN *pRoot);
#endif
#endif /* COLLISIONS */

#ifdef SLIDING_PATCH
#define SHEAR(ix,lx,ly,w,t)\
	((ix) < 0 ? fmod(0.5*(ly) - 1.5*(ix)*(w)*(lx)*(t),(ly)) - 0.5*(ly):\
	 (ix) > 0 ? 0.5*(ly) - fmod(0.5*(ly) + 1.5*(ix)*(w)*(lx)*(t),(ly)): 0.0)
void pkdPatch(PKD pkd, double dOrbFreqZ2);
#endif

#ifdef NEED_VPRED
#ifdef GASOLINE
void pkdKickVpred(PKD pkd, double dvFacOne, double dvFacTwo, double duDelta,
				  int iGasModel, double z, double duDotLimit);
#else
void pkdKickVpred(PKD pkd, double dvFacOne, double dvFacTwo);
#endif
#endif

#endif
void pkdCOM(PKD pkd, double *com);
void pkdCOMByType(PKD pkd, int type, double *com);
void pkdOldestStar(PKD pkd, double *com);
