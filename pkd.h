#ifndef PKD_HINCLUDED
#define PKD_HINCLUDED

#include <sys/resource.h>
#include "mdl.h"
#include "floattype.h"
#include "cooling.h"

#define CID_TOP		0
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
        unsigned int iTouchRung; /* To be removed */
	int cpStart;
	FLOAT fWeight;
	FLOAT fMass;
	FLOAT fSoft;
	FLOAT r[3];
	FLOAT v[3];
	FLOAT a[3];
	FLOAT fPot;
	FLOAT fBall2;
        FLOAT fDensity;
	FLOAT dt;			/* a time step suggestion */
#ifdef SUPERCOOL
	FLOAT vMean[3];
#endif
#ifdef COLORCODE
	FLOAT fColor;
#endif
#ifdef GASOLINE
        FLOAT h;                /* SPH h */
	FLOAT vPred[3];		/* predicted velocity (time centered) */
	FLOAT uPred;		/* predicted thermal energy */
	FLOAT PoverRho2;	/* P/rho^2 */
	FLOAT u;	        /* thermal energy */ 
        FLOAT c;                /* sound speed */
        FLOAT mumax;            /* sound speed like viscosity term */
        FLOAT PdV;	        /* P dV heating (includes shocking) */
#ifdef DEBUG
        FLOAT PdVvisc;          /* P dV from shock (testing) */
        FLOAT PdVpres;          /* P dV from adiabatic compression (testing) */
#endif
        FLOAT divv;             /* Balsara viscosity reduction */
        FLOAT curlv[3];         
        FLOAT BalsaraSwitch;
  
        FLOAT udot;             /* Rate of change of u -- for predicting */
        FLOAT Y_HI, Y_HeI, Y_HeII;  /* Abundance of ions */

	FLOAT fMetals;
	FLOAT fTimeForm;
#endif
#ifdef COLLISIONS
	FLOAT fRedHill;		/* radius of reduced Hill sphere, times HILL_SCALE */
	FLOAT w[3];			/* spin vector */
	int iColor;			/* handy color tag */
	int iDriftType;		/* either NORMAL or KEPLER */
	double dTEnc;		/* time to next encounter */
/*DEBUG	PARTICLE_ID idNbr;*/	/* encounter neighbor id */
#ifdef SAND_PILE
	int bStuck;
#endif /* SAND_PILE */
#endif /* COLLISIONS */

	} PARTICLE;

/* Active Type Masks */

/* Active: -- eg. Calculate new acceleration, PdV, etc... for this particle */
#define TYPE_ACTIVE            1
/* In the Tree: */
#define TYPE_TREEACTIVE        2
/* Gather to/Scatter from this particle with in smooths: */
#define TYPE_SMOOTHACTIVE      4
/* Smooth has processed this particle */
#define TYPE_SMOOTHDONE        8

/* Types used for Fast Density only (so far) */
/* Sum Fast Density on this particle */
#define TYPE_DensACTIVE        16
/* Neighbour of active particle: */
/* Neighbour of DENSACTIVE (incl. DENSACTIVE): */
#define TYPE_NbrOfDensACTIVE           32
/* Neighbour of Neighbour of DENSACTIVE (incl. DENSACTIVE & NACTIVE): */
#define TYPE_NbrOfNbrOfDensACTIVE          64
/* Neighbour of NOT-DENSACTIVE: */
#define TYPE_NbrOfNonDensACTIVE        128

/* Particle Type Masks */

#define TYPE_GAS               256
#define TYPE_DARK              512
#define TYPE_STAR              1024
#define TYPE_SUPERCOOL         2048

/* Combination Masks */
#define TYPE_ALL               (TYPE_GAS|TYPE_DARK|TYPE_STAR)

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
#define TYPEClearACTIVE(a)       ((a)->iActive &= (TYPE_GAS|TYPE_DARK|TYPE_STAR|TYPE_SUPERCOOL))
#define TYPEClear(a)             ((a)->iActive = 0)

#define CHECKPOINT_VERSION 6

typedef struct chkParticle {
	int iOrder;
	FLOAT fMass;
	FLOAT fSoft;
	FLOAT r[3];
	FLOAT v[3];
#ifdef GASOLINE
	FLOAT u;
#endif
#ifdef COLLISIONS
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

typedef struct ewaldTable {
	double hx,hy,hz;
	double hCfac,hSfac;
	} EWT;

#ifdef COLLISIONS

typedef struct {
	int iPid;
	int iIndex;
	int iOrder;
	} PARTICLE_ID;

typedef struct {
	PARTICLE_ID id;
	FLOAT fMass;
	FLOAT fRadius; /*DEBUG = 2*fSoft so spline not used in force calcs*/
	FLOAT r[3];
	FLOAT v[3];
	FLOAT w[3];
	FLOAT dt;
	} COLLIDER;

#endif /* COLLISIONS */

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
		} ti[MAX_TIMERS];
        /* 
	 ** Cooling 
	 */
        CL  cl;
#ifdef COLLISIONS
	double dImpactTime;
	COLLIDER Collider1,Collider2;
#endif /* COLLISIONS */
	} * PKD;


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

enum GasModel { GASMODEL_ADIABATIC, 
		GASMODEL_ISOTHERMAL, 
                GASMODEL_COOLING, 
		GASMODEL_COOLING_NONEQM,
   		GASMODEL_GLASS }; 


#ifdef GASOLINE

#endif

#define PKD_ORDERTEMP	256

#define pkdRoot(iCell,id)\
{\
	id = -1;\
	iCell = ROOT;\
	}

#define pkdIsRoot(iCell,id)		((id==-1)?((iCell==ROOT)?1:0):0)

#define pkdLower(pkd,iCell,id)\
{\
	if (id == -1) {\
		id = pkd->kdTop[iCell].pLower;\
		if (id != -1) iCell = ROOT;\
		}\
	iCell = LOWER(iCell);\
	}

#define pkdUpper(pkd,iCell,id)\
{\
	if (id == -1) {\
		id = pkd->kdTop[iCell].pLower;\
		if (id != -1) iCell = ROOT;\
		}\
	iCell = UPPER(iCell);\
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
void pkdClearTimer(PKD,int);
void pkdStartTimer(PKD,int);
void pkdStopTimer(PKD,int);
void pkdInitialize(PKD *,MDL,int,int,int,FLOAT *,int,int,int);
void pkdFinish(PKD);
void pkdReadTipsy(PKD,char *,int,int,int,double,double);
void pkdSetSoft(PKD pkd,double dSoft);
void pkdCalcBound(PKD,BND *,BND *,BND *);
int pkdWeight(PKD,int,FLOAT,int,int,int,int *,int *,FLOAT *,FLOAT *);
int pkdLowerPart(PKD,int,FLOAT,int,int);
int pkdUpperPart(PKD,int,FLOAT,int,int);
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
void pkdGravAll(PKD,int,int,int,int,double,double,int,double *,int *,
				double *,double *,double *,CASTAT *,double *); 
void pkdCalcE(PKD,double *,double *,double *);
void pkdDrift(PKD,double,FLOAT *,int,int,FLOAT);
void pkdKick(PKD pkd,double,double, double, double, double, double, int, double, double);
void pkdReadCheck(PKD,char *,int,int,int,int);
void pkdWriteCheck(PKD,char *,int,int);
void pkdDistribCells(PKD,int,KDN *);
void pkdCalcRoot(PKD,struct ilCellNewt *);
void pkdDistribRoot(PKD,struct ilCellNewt *);
void pkdSwapAll(PKD pkd, int idSwap);
double pkdMassCheck(PKD pkd);
void pkdSetRung(PKD pkd, int iRung);
int pkdActiveRung(PKD pkd, int iRung, int bGreater);
int pkdCurrRung(PKD pkd, int iRung);
void pkdDensityStep(PKD pkd, double dEta, double
		    dRhoFac);
void pkdAccelStep(PKD pkd, double dEta, double dVelFac, double
		     dAccFac, int bDoGravity, int bEpsVel, int bSqrtPhi);
int pkdDtToRung(PKD pkd, int iRung, double dDelta, int iMaxRung, int
		bAll);
void pkdInitDt(PKD pkd, double dDelta);

int pkdRungParticles(PKD,int);
void pkdCoolVelocity(PKD,int,double,double,double);
void pkdGrowMass(PKD pkd,int nGrowMass, double dDeltaM);
void pkdInitAccel(PKD);
int pkdOrdWeight(PKD,int,int,int,int,int *,int *);
void pkdDeleteParticle(PKD pkd, int i);
void pkdNewParticle(PKD pkd, PARTICLE p);

int pkdResetTouchRung(PKD pkd, unsigned int iTestMask, unsigned int iSetMask);
int pkdActiveExactType(PKD pkd, unsigned int iFilterMask, unsigned int iTestMask, unsigned int iSetMask);
int pkdActiveType(PKD pkd, unsigned int iTestMask, unsigned int iSetMask);
int pkdSetType(PKD pkd, unsigned int iTestMask, unsigned int iSetMask);
int pkdResetType(PKD pkd, unsigned int iTestMask, unsigned int iSetMask);
int pkdCountType(PKD pkd, unsigned int iFilterMask, unsigned int iTestMask);
int pkdActiveMaskRung(PKD pkd, unsigned int iSetMask, int iRung, int bGreater );
int pkdActiveTypeRung(PKD pkd, unsigned int iTestMask, unsigned int iSetMask, int iRung, int bGreater);

int pkdIsGas(PKD,PARTICLE *);
int pkdIsDark(PKD,PARTICLE *);
int pkdIsStar(PKD,PARTICLE *);
void pkdSetParticleTypes(PKD pkd, int nSuperCool);

void pkdColNParts(PKD pkd, int *pnNew, int *nDeltaGas, int *nDeltaDark,
		  int *nDeltaStar);
void pkdNewOrder(PKD pkd, int nStart);
void pkdSetNParts(PKD pkd, int nGas, int nDark, int nStar, int nMaxOrderGas,
		  int nMaxOrderDark);
void pkdSunIndirect(PKD,double *,int,double);
void pkdLogHalo(PKD);
void pkdHernquistSpheroid(PKD pkd);
void pkdMiyamotoDisk(PKD pkd);

#ifdef GASOLINE
void pkdAdiabaticGasPressure(PKD, double gammam1, double gamma);
void pkdInitEnergy(PKD pkd, double dTuFac, double z);
#ifdef GLASS
void pkdGlassGasPressure(PKD, void *in);
#endif
void pkdKickVpred(PKD pkd, double dvFacOne, double dvFacTwo, double duDelta,int iGasModel, double z, double duDotLimit);
void pkdKickRhopred(PKD pkd, double dHubbFac, double dDelta);
int pkdSphCurrRung(PKD pkd, int iRung, int bGreater);
void pkdSphStep(PKD pkd, double dCosmoFac, double dEtaCourant, double dEtauDot);
void pkdSphViscosityLimiter(PKD pkd, int bOn);

void pkdPARTICLE2PERBARYON(PERBARYON *Y, PARTICLE *p, double HTotal, double HeTotal);

#define pkdPARTICLE2PERBARYON(Y, p, HTotal, HeTotal) { \
    (Y)->HI = (p)->Y_HI; \
    (Y)->HII = HTotal - (Y)->HI; \
    (Y)->HeI = (p)->Y_HeI; \
    (Y)->HeII = (p)->Y_HeII; \
    (Y)->HeIII = HeTotal - (Y)->HeI - (Y)->HeII; \
    (Y)->e = (Y)->HII + (Y)->HeII + 2*(Y)->HeIII; \
    (Y)->Total = (Y)->e + HTotal + HeTotal; }

void pkdPERBARYON2PARTICLE(PERBARYON *Y, PARTICLE *p);

#define pkdPERBARYON2PARTICLE(Y, p) { \
    (p)->Y_HI = (Y)->HI; \
    (p)->Y_HeI = (Y)->HeI; \
    (p)->Y_HeII = (Y)->HeII; }

#endif
#ifdef GLASS
void pkdRandomVelocities(PKD pkd, double dMaxVL, double dMaxVR);
#endif

#ifdef COLLISIONS
int pkdNumRejects(PKD pkd);
void pkdReadSS(PKD pkd, char *pszFileName, int nStart, int nLocal);
void pkdWriteSS(PKD pkd, char *pszFileName, int nStart);
void pkdCalcHill(PKD pkd, double dCentMass);
void pkdHillStep(PKD pkd, double dEta);
void pkdFindEncounter(PKD pkd, double *dNext);
void pkdMarkEncounters(PKD pkd, double dt);
int pkdLowerQQPart(PKD pkd, int d, FLOAT fSplit, int i, int j);
int pkdUpperQQPart(PKD pkd, int d, FLOAT fSplit, int i, int j);
void pkdQQCalcBound(PKD pkd, BND *pbnd, BND *pbndActive);
void pkdQQBuild(PKD pkd, int nBucket, int bActiveOnly, KDN *pRoot);
#endif /* COLLISIONS */

#endif




