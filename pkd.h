
#ifndef PKD_HINCLUDED
#define PKD_HINCLUDED

#include <sys/time.h>
#include <sys/resource.h>
#include "mdl.h"
#include "floattype.h"
#include "treezip.h"

#ifdef GASOLINE
#include "cooling.h"
#endif
#include "rotbar.h"

#ifdef SLIDING_PATCH
#include "patch.h"
#endif

/* Allow for creation of this many new gas particles */
#if defined(INFLOWOUTFLOW) || defined(FBPARTICLE) || defined(PARTICLESPLIT)
#define NIORDERGASBUFFER 1000000000
#else
#define NIORDERGASBUFFER 0
#endif

/*
** The following sort of definition should really be in a global
** configuration header file -- someday...
*/

#if defined(GASOLINE) || defined(ROT_FRAME) || defined(SIMPLE_GAS_DRAG) || defined(GR_DRAG)
#define NEED_VPRED
#endif

//There are a handful of TYPE bits we never want leaving a local thread or a smooth call.
#define TYPE_MASK (~TYPE_RESMOOTHINNER & ~TYPE_MARK)

/* SPH variable ALPHA */
#define ALPHAMIN 0.01
#define ALPHAMAX 1

#ifdef RTF
#define RTDENSITY
#define RTFORCE
#endif

#ifdef SUPERBUBBLE
#define PROMOTE
#define THERMALCOND
#define TWOPHASE
#define TOPHATFEEDBACK
#endif

#if defined(TWOPHASE)  && !defined(UNONCOOL)
#define UNONCOOL
#endif

#if defined(DENSITYUNOTP) && !defined(DENSITYU)
#define DENSITYU
#endif

#ifdef DIFFUSION

#if defined(FEEDBACKDIFFLIMIT) && !defined(DIFFUSIONHARMONIC)
#define DIFFUSIONHARMONIC
#endif

#define DIFFRATE(p_) ((p_)->diff)
#else
#define DIFFRATE(p_) (1e-30)

#endif

/* Note: UDOT_HYDRO is only correct if there is only thermal pressure (no UNONCOOL or Jeans Floor) */
#define UDOT_HYDRO(p_)   ((p_)->uDotPdV+(p_)->uDotAV+(p_)->uDotDiff)
#ifndef PONRHOFLOOR
#define PONRHOFLOOR 0
#endif

#define StarClusterFormfBall2Save(p) (p->curlv[0])
#define StarClusterFormiOrder(p) (p->curlv[1])

/* (note bVWarnings still applies) */

#define CID_TOP                 0
#define CID_PARTICLE    0
#define CID_CELL                1

#define ROOT            1
#define LOWER(i)        (i<<1)
#define UPPER(i)        ((i<<1)+1)
#define PARENT(i)       (i>>1)
#define SIBLING(i)      ((i&1)?i-1:i+1)
#define SETNEXT(i)                              \
    {                                           \
    while (i&1) i=i>>1;                         \
    ++i;                                        \
    }

#define MAX_TIMERS              10

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
    FLOAT dt;               /* a time step suggestion */
    FLOAT dtNew;            /* SPH new dt estimate */
    FLOAT dtOld;            /* SPH Old dt */
    FLOAT dtGrav;           /* suggested 1/dt^2 from gravity */
#ifdef SLIDING_PATCH
    FLOAT dPy;              /* Canonical momentum for Hill eqn. */
#endif
#ifdef SUPERCOOL
    FLOAT vMean[3];
#endif
    FLOAT fBallMax;         /* SPH 2h Max value */
#ifdef GASOLINE
    FLOAT c;                /* sound speed */
    FLOAT PoverRho2;        /* P/rho^2 */
    FLOAT mumax;            /* sound speed like viscosity term (OBSOLETE) */
    FLOAT u;                /* thermal energy  */ 
    FLOAT uPred;            /* predicted thermal energy, */
    FLOAT uDotPdV;          /* PdV heating                          [Sink Lx] */
    FLOAT uDotAV;           /* Shock Heating (Artificial Viscosity) [Sink Ly] */
    FLOAT uDotDiff;         /* Thermal Energy diffusion             [Sink Lz] */
#ifndef NOCOOLING
    FLOAT uDot;                 /* Rate of change of u -- for predicting u */
    COOLPARTICLE CoolParticle;  /* Abundances and any other cooling internal variables */
#endif
#ifdef TWOPHASE
    FLOAT fMassHot;
    COOLPARTICLE CoolParticleHot;  /* Abundances and any other cooling internal variables */
#endif
#ifdef UNONCOOL
    FLOAT uHot;
    FLOAT uHotPred;
    FLOAT uHotDot;
    FLOAT uHotDotDiff;  /* Hot Energy diffusion */
#endif
    FLOAT divv;             
    FLOAT dvds;
#ifdef VARALPHA
    FLOAT alpha;
    FLOAT alphaPred;
#endif
#ifdef CULLENDEHNEN
    FLOAT alpha;
    FLOAT dTime_divv;
    FLOAT divv_old; // stored old value for checking that nbrs also all compressing
#endif
    FLOAT curlv[3];         /* Note this is used as workspace and value is not preserved */
    FLOAT BalsaraSwitch;    /* Balsara viscosity reduction */
#ifdef THERMALCOND
    FLOAT fThermalCond;
    FLOAT fThermalLength;
#endif
#ifdef DIFFUSION
    FLOAT diff;
    FLOAT fMetalsDot;
    FLOAT fMetalsPred;
#endif
#ifdef DENSITYU
    FLOAT fDensityU;
#endif
    FLOAT fDivv_t;
    FLOAT fDivv_Corrector;
#ifdef SINKING
    FLOAT rSinking0Unit[3];
    FLOAT rSinking0Mag;
    FLOAT vSinkingTang0Unit[3];
    FLOAT vSinkingTang0Mag;
    FLOAT vSinkingr0;
    FLOAT fSinkingTime;  
    FLOAT fTrueMass;
    int iSinkingOnto;
#endif
#ifdef SURFACEAREA
    FLOAT fArea; 
#ifdef NORMAL
    FLOAT normal[3];
#endif
#endif
/*      FLOAT fDensSave;*/      /* Used by diagnostic DensCheck funcs */
    FLOAT fMetals;  /* mass fraction in metals, a.k.a, Z */
    FLOAT fTimeForm;
#ifdef SFBOUND
    FLOAT fSigma2;
#endif
#ifdef STARFORM
    FLOAT uDotFB;
    FLOAT uDotESF;
    FLOAT fMSN;
    FLOAT fNSN;           
    FLOAT fMOxygenOut;
    FLOAT fMIronOut;
    FLOAT fMFracOxygen;
    FLOAT fMFracIron;
#ifdef DIFFUSION
    FLOAT fMFracOxygenDot;
    FLOAT fMFracIronDot;
    FLOAT fMFracOxygenPred;
    FLOAT fMFracIronPred;
#endif
    FLOAT fSNMetals;
    FLOAT fNSNtot;
    FLOAT fTimeCoolIsOffUntil;
    FLOAT fMassForm;        /* record original mass of star */
    int iGasOrder;          /* gas from which star formed */
#endif
#endif  /* GASOLINE */
#ifdef COLLISIONS
    int iOrgIdx;            /* for tracking of mergers, aggregates etc. */
    FLOAT w[3];                     /* spin vector */
    int iColor;                     /* handy color tag */
    int iDriftType;         /* either NORMAL or KEPLER */
    double dtCol;           /* time to next encounter or collision */
    int iOrderCol;          /* neighbour or collider iOrder */
    double dtPrevCol;       /* time of previous collision */
    int iPrevCol;           /* iOrder of previous collider */
    int bTinyStep;          /* flag for imminent collapse */
    FLOAT mindist2;         /* record min dist for all encounters */
    int bGhostExclude;      /* particle not included in ghost cells */
#endif /* COLLISIONS */
#ifdef SLIDING_PATCH
    int bAzWrap;        /* flag set on azimuthal boundary wrap */
#endif
#ifdef SAND_PILE
    int bStuck;
#endif
#ifdef NEED_VPRED
    FLOAT vPred[3];         /* predicted velocity (time centered) */
#endif
#ifdef AGGS
    /*
    ** Position of particle in principal frame of the aggregate
    ** (normally).  We temporarily store the COM frame position
    ** here during the process of computing the aggregate
    ** parameters.
    */
    FLOAT r_agg[3];
#endif
#ifdef RUBBLE_ZML
    double dDustMass;       /* predicted mass increase from dust */
    int iBin;                               /* dust bin that planetesimal is in */
    int bMayCollide;        /* true if planetesimal is predicted to
                               collide with another planetesimal during
                               the top step interval */
#endif
    } PARTICLE;

#define GAMMA_JEANS    (2.0)
#define GAMMA_NONCOOL  (5./3.)

#ifdef GLASS
struct GlassData {
    /* Glass */
	double dGlassPoverRhoL;
	double dGlassPoverRhoR;
	double dGlassxL;
	double dGlassxR;
	double dxBoundL;
	double dxBoundR;
    double dGamma;
    };
#endif

struct GasPressureContext {
    int iGasModel;
  /* Adiabatic */
	double gamma;
	double gammam1;
	double dResolveJeans;
	double dCosmoFac;
    double dtFacCourant;
  /* Isothermal */
  /* Ion evolving */
#ifdef GLASS
    struct GlassData g;
#endif
#if defined(THERMALCOND) || defined(TWOPHASE)
    double dThermalCondCoeffCode;
    double dThermalCondSatCoeff;
    double dThermalCond2CoeffCode;
    double dThermalCond2SatCoeff;
#endif
#if defined(PROMOTE) || defined(TWOPHASE)
	double dEvapCoeffCode;
	double dEvapMinTemp;
#endif
    };

typedef struct uHotContext {
    double dHotConvRate;
    double dHotConvRateMul;
    double dHotConvRateMax;
    double dHotConvUMin;
    struct GasPressureContext gpc;
#ifdef TWOPHASE
    double dMultiPhaseMinTemp;
    double dMultiPhaseMaxFrac;
    double dMultiPhaseMaxTime;
#endif
    } UHC;

#define SINK_Lx(_a) (((PARTICLE *) (_a))->uDotPdV)
#define SINK_Ly(_a) (((PARTICLE *) (_a))->uDotAV)
#define SINK_Lz(_a) (((PARTICLE *) (_a))->uDotDiff)

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
#define TYPE_SINK              (1<<14)
#define TYPE_SINKING           (1<<15)
#define TYPE_NEWSINKING        (1<<16)
#define TYPE_INFLOW            (1<<17)
#define TYPE_OUTFLOW           (1<<18)
#define TYPE_FEEDBACK          (1<<19)
#define TYPE_PROMOTED          (1<<20)
#define TYPE_DENMAX            (1<<21)
#define TYPE_STARFORM          (1<<22)
#define TYPE_MARK              (1<<23)
#define TYPE_RESMOOTHINNER     (1<<24)
#define TYPE_TWOPHASE          (1<<25)

/* Combination Masks */
#define TYPE_ALLACTIVE                  (TYPE_ACTIVE|TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE)
#define TYPE_ALL                                (TYPE_GAS|TYPE_DARK|TYPE_STAR)

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

enum CheckSanityProblem {
    PROBLEM_IORDER,
    PROBLEM_MASS,
    PROBLEM_SOFT,
    PROBLEM_POSITION,
    PROBLEM_VELOCITY,
    PROBLEM_U,
    PROBLEM_METALS,

    PROBLEM_ENDLIST
    }; 

/* Alternate if no detailed sanity check is done */
#define CHECKSANEALT(xxx_) xxx_
#define CHECKSANE(w_,x_,y_,z_)
#define CHECKSANEFIX(nProb_,iVar_,iTest_)

#ifdef RUBBLE_ZML
/* RUBBLE_ZML puts extra stuff in the checkpoint file, changes version to indicate this  ZML 01.08.04 */
#define CHECKPOINT_VERSION 81
#else
#define CHECKPOINT_VERSION 8
#endif

typedef struct chkParticle {
    int iOrder;
    int iActive;
    FLOAT fMass;
    FLOAT fSoft;
    FLOAT r[3];
    FLOAT v[3];
#ifdef GASOLINE
    FLOAT u;
#ifdef TWOPHASE
    FLOAT fMassHot;
    COOLPARTICLE CoolParticleHot;  /* Abundances and any other cooling internal variables */
#endif
#ifdef UNONCOOL
    FLOAT uHot;
#endif
#ifdef STARSINK
    FLOAT Lx,Ly,Lz;
#endif
#ifdef VARALPHA
    FLOAT alpha;
#endif
    FLOAT fMetals;
#ifndef NOCOOLING
    COOLPARTICLE CoolParticle;
#endif
    FLOAT fTimeForm;
#ifdef STARFORM
    FLOAT fTimeCoolIsOffUntil;
    FLOAT fMassForm;        /* record original mass of star */
    FLOAT fNSN;
    FLOAT fMFracOxygen;
    FLOAT fMFracIron;
    int iGasOrder;
#endif
#ifdef SINKING
    FLOAT rSinking0Unit[3];
    FLOAT rSinking0Mag;
    FLOAT vSinkingTang0Unit[3];
    FLOAT vSinkingTang0Mag;
    FLOAT vSinkingr0;
    FLOAT fSinkingTime;  
    FLOAT fTrueMass;
    int iSinkingOnto; /* used for nSinkingOnto for sink itself */
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

/* Used by bLongRangeStep, and -D ONGRANGESTEP */
typedef struct bndDt {
    FLOAT vMin[3],vMax[3];
    FLOAT cMax,drMax2;
    } BNDDT;

#define DIAGDIST2(fDist2,rMin,rMax) {                   \
        FLOAT DD_dx,DD_dy,DD_dz;                        \
        DD_dx = (rMax)[0] - (rMin)[0];                  \
        DD_dy = (rMax)[1] - (rMin)[1];                  \
        DD_dz = (rMax)[2] - (rMin)[2];                  \
        fDist2 = DD_dx*DD_dx+DD_dy*DD_dy+DD_dz*DD_dz; }

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
#ifdef  RADIATIVEBOX
    double fLW;
    double gmass;
    double gmom;
    FLOAT cLumLW[3];
#endif
    };


typedef struct kdNode {
    int iDim;
    double fSplit;
	BND bnd;
	BND bndBall;	/* Bound including fBall*(1+changemax) */
    BNDDT bndDt;
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

typedef struct sfEvent 		/* Holds statistics of the star
                               formation event */
    {
    int iOrdStar;
    int iOrdGas;
    double timeForm;
    double rForm[3];
    double vForm[3];
    double massForm;
    double rhoForm;
    double TForm;
#ifdef SFEVENTCRIT
    double tcool;
    double tdyn;
#endif
#ifdef COOLING_MOLECULARH
    double H2fracForm;
#endif
    } SFEVENT;

typedef struct starLog
    {
    int nLog;			/* number of events in buffer */
    int nMaxLog;		/* max size of buffer; increase when needed */
    int nOrdered;		/* The number of events that have been
                           globally ordered, incremented by
                           pkdNewOrder() */
    SFEVENT *seTab;		/* The actual table */
    } STARLOG;

enum SinkEventType {
    SINK_EVENT_NULL,
    SINK_EVENT_ACCRETE_AT_FORMATION,
    SINK_EVENT_FORM,
    SINK_EVENT_ACCRETE,
    SINK_EVENT_ACCRETE_UPDATE,
    SINK_EVENT_MERGER
    }; 

typedef struct sinkEvent 	/* Holds statistics of the sink/accretion/merger
                               formation event */
    {
    int iOrdSink;
    int iOrdVictim;
    double time;
    double mass;
    double r[3];
    double v[3];
    double L[3];
    } SINKEVENT;

typedef struct sinkLog
    {
    int nLog;			/* number of events in buffer */
    int nMaxLog;		/* max size of buffer; increase when needed */
    int nLogOrdered;		/* The number of events that have been
                               globally ordered, incremented by
                               pkdNewOrder() prior to flush */
    int nFormOrdered;           /* Number of formation events ordered ... */
    int nForm;                  /* Number of new sink formation events in table */
    int nAccrete;               /* Gas Accrete events */
    int nMerge;                 /* Sink Merger events */
    SINKEVENT *SinkEventTab;		/* The actual table */
    } SINKLOG;

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
    int nMaxOrder;
	int nBucket;
	int nLevels;
	int nSplit;
	int nNodes;
	int iExtraBucket;
	int iOrder;
	int iFreeCell;
	int iRoot;
	FLOAT fPeriod[3];
	FLOAT dxInflow, dxOutflow;
	int *piLeaf;
	KDN *kdTop;
	KDN *kdNodes;
	PARTICLE *pStore;
    double duTFac;
    double dvFac;
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
#ifdef OUTURBDRIVER
    void *outurb;
#endif
    STARLOG starLog;
#endif
    SINKLOG sinkLog;
#ifdef SLIDING_PATCH
    /*
    ** Info needed for sliding patch model...
    */
    double dTime;
    PATCH_PARAMS *PP;
#endif
    ROTBAR  rotbar;
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

#define pkdRoot(iCell,id)                       \
    {                                           \
	id = -1;                                    \
	iCell = ROOT;                               \
	}

#define pkdIsRoot(iCell,id)		((id==-1)?((iCell==ROOT)?1:0):0)

/*
 * There is now a slight inefficency here when going from the top tree to
 * a node tree in that we visit the root cell twice (once as the leaf of the
 * top tree and once as the root of the node tree).  This is necessary to
 * check if the root cell is a bucket.
 */

#define pkdLower(pkd,iCell,id)                  \
    {                                           \
	if (id == -1) {                             \
		id = pkd->kdTop[iCell].pLower;          \
		if (id != -1) iCell = ROOT;             \
		else iCell = LOWER(iCell);              \
		}                                       \
	else iCell = LOWER(iCell);                  \
	}

#define pkdUpper(pkd,iCell,id)                  \
    {                                           \
	if (id == -1) {                             \
		id = pkd->kdTop[iCell].pLower;          \
		if (id != -1) iCell = ROOT;             \
		else iCell = UPPER(iCell);              \
		}                                       \
	else iCell = UPPER(iCell);                  \
	}

#define pkdParent(pkd,iCell,id)                 \
    {                                           \
	iCell = PARENT(iCell);                      \
	if (iCell == ROOT) {                        \
		if (id != -1) {                         \
			iCell = pkd->piLeaf[id];            \
			id = -1;                            \
			}                                   \
		}                                       \
	}

#define pkdNext(pkd,iCell,id)                   \
    {                                           \
	SETNEXT(iCell);                             \
	if (iCell == ROOT) {                        \
		if (id != -1) {                         \
			iCell = pkd->piLeaf[id];            \
			id = -1;                            \
			SETNEXT(iCell);                     \
			}                                   \
		}                                       \
	}

double pkdGetTimer(PKD,int);
double pkdGetSystemTimer(PKD,int);
double pkdGetWallClockTimer(PKD,int);
void pkdClearTimer(PKD,int);
void pkdStartTimer(PKD,int);
void pkdStopTimer(PKD,int);
void pkdInitialize(PKD *,MDL,int,int,int,FLOAT *,FLOAT,FLOAT,int,int,int);
void pkdFinish(PKD);
void pkdReadTipsy(PKD,char *,int,int,int,int,double,double);
void pkdOutputBlackHoles(PKD pkd,char *pszFileName, double dvFac);
void pkdSetSoft(PKD pkd,double dSoft);
#ifdef CHANGESOFT
void pkdPhysicalSoft(PKD pkd,double, double, int);
void pkdPreVariableSoft(PKD pkd,int iVariableSoftType);
void pkdPostVariableSoft(PKD pkd,double dSoftMax,int bSoftMaxMul,int iVariableSoftType);
#endif
void pkdCalcBound(PKD,BND *,BND *,BND *,BND *, BNDDT *);
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
void pkdTotals(PKD pkd, int *nDark, int *nGas, int *nStar);
int pkdActive(PKD);
int pkdTreeActive(PKD);
int pkdInactive(PKD);
int pkdTreeInactive(PKD);
int pkdNodes(PKD);
void pkdDomainColor(PKD);
int pkdColOrdRejects(PKD,int,int);
void pkdLocalOrder(PKD);
void pkdWriteTipsy(PKD,char *,int,int,double,double,int);
void pkdTreeZip(PKD pkd,char *pszFileName, double *dmin, double *dmax);
void pkdCombine(KDN *,KDN *,KDN *);
void pkdCalcCell(PKD,KDN *,FLOAT *,int,struct pkdCalcCellStruct *);
double pkdCalcOpen(KDN *,int,double,int);
void pkdBuildLocal(PKD,int,int,double,int,int,int,KDN *);
void pkdBuildBinary(PKD,int,int,double,int,int,int,KDN *);
void pkdThreadTree(PKD pkd,int iCell,int iNext);
void pkdGravAll(PKD pkd,int nReps,int bPeriodic,int iOrder, int bEwald,int iEwOrder,
    double fEwCut,double fEwhCut, int bComove, double dRhoFac,
    int bDoSun,double dSunSoft, double *aSun,int *nActive,
    double *pdPartSum,double *pdCellSum,double *pdSoftSum,CASTAT *pcs,
    double *pdFlop);
void pkdCalcEandL(PKD,double *,double *,double *,double []);
void pkdCalcEandLExt(PKD,double *,double[],double [],double *);
void pkdDrift(PKD,double,FLOAT *,int,int,int,FLOAT,double);

double pkduHotConvRate(PKD pkd, UHC uhc, FLOAT fBall2, double uHotPred, double uPred);
void pkdUpdateuDot(PKD pkd, double duDelta, double dTime, double z, UHC uhc, int iGasModel, int bUpdateState );
void pkdKick(PKD pkd, double dvFacOne, double dvFacTwo, double dvPredFacOne,
    double dvPredFacTwo, double duDelta, double duPredDelta, int iGasModel,
    double z, double duDotLimit, double dTimeEnd, UHC uhc );
void pkdEmergencyAdjust(PKD pkd, int iRung, int iMaxRung, double dDelta, double dDeltaThresh, int *pnUn, int *piMaxRungIdeal, int *pnMaxRung, int *piMaxRungOut);
void pkdKickPatch(PKD pkd, double dvFacOne, double dvFacTwo,
    double dOrbFreq, int bOpen);
void pkdGravInflow(PKD pkd, double r);
void pkdCreateInflow(PKD pkd, int Ny, int iGasModel, double dTuFac, double pmass, double x, double vx, double density, double temp, double metals, double eps, double dt, int iRung);
void pkdReadCheck(PKD,char *,int,int,int,int);
void pkdWriteCheck(PKD,char *,int,int);
void pkdDistribCells(PKD,int,KDN *);
void pkdCalcRoot(PKD,struct ilCellNewt *);
void pkdDistribRoot(PKD,struct ilCellNewt *);
void pkdSwapAll(PKD pkd, int idSwap);
double pkdMassCheck(PKD pkd);
void pkdMassMetalsEnergyCheck(PKD pkd, double *dTotMass, double *dTotMetals, 
    double *dTotOx, double *dTotFe, double *dTotEnergy);
void pkdSetRung(PKD pkd, int iRung);
void pkdBallMax(PKD pkd, int iRung, int bGreater, double ddHonHLimit);
int pkdActiveRung(PKD pkd, int iRung, int bGreater);
int pkdCurrRung(PKD pkd, int iRung);
void pkdGravStep(PKD pkd, double dEta);
void pkdAccelStep(PKD pkd, double dEta, double dVelFac, double
    dAccFac, int bDoGravity, int bEpsAcc, int bSqrtPhi, double dhMinOverSoft);
void pkdDensityStep(PKD pkd, double dEta, double dRhoFac);

int pkdOneParticleDtToRung( int iRung,double dDelta,double dt);
int pkdDtToRung(PKD pkd, int iRung, double dDelta, int iMaxRung, int bAll, int bDiagExceed,
    int *pnMaxRung, int *piMaxRungIdeal );
void pkdInitDt(PKD pkd, double dDelta);
int pkdRungParticles(PKD,int);
void pkdCoolVelocity(PKD,int,double,double,double);
void pkdGrowMass(PKD pkd,int nGrowMass, int iGrowType, double dDeltaM, double dMinM, double dMaxM);
void pkdInitAccel(PKD);
void pkdModifyAccel(PKD pkd, double);
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
int pkdSetTypeFromFile(PKD pkd, int iSetMask, int biGasOrder, char *file, int *pniOrder, int *pnSet, int *pnSetiGasOrder);

void pkdSetParticleTypes(PKD pkd, int nSuperCool);

struct SoughtParticle {
    int iOrder;
    int iActive; 
    double x,y,z,m;
    };

int pkdSoughtParticleList(PKD pkd, int iTypeSought, int nMax, int *n, struct SoughtParticle *sp);

void pkdCoolUsingParticleList(PKD pkd, int nList, struct SoughtParticle *sp);

void pkdColNParts(PKD pkd, int *pnNew, int *nAddGas, int *nAddDark,
    int *nAddStar, int *nDelGas, int *nDelDark, int *nDelStar);

void pkdNewOrder(PKD pkd, int nStartGas, int nStarDark, int nStartStar);
void pkdMoveParticle(PKD pkd, double *xcenter,double *xoffset,int iOrder);


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
void pkdSunIndirect(PKD,double *,int,double,double);
void pkdLogHalo(PKD, double, double, double, double);
void pkdHernquistSpheroid(PKD pkd);
void pkdNFWSpheroid(PKD pkd, double M_200, double r_200, double c, double dSoft);
void pkdElliptical(PKD pkd, int bEllipticalDarkNFW);
void pkdHomogSpheroid(PKD pkd, double M_s, double R_s);
void pkdBodyForce(PKD pkd, double dConst);
void pkdGalaxyDiskVerticalPotentialForce(PKD pkd, double Vc, double R, double StarSigma, double StarH, double GasSigma, double GasH, double Gasa, double Gasb, double Gasc);
void pkdMiyamotoDisk(PKD pkd);
void pkdTimeVarying(PKD pkd,double dTime);
double pkdDtFacCourant( double dEtaCourant, double dCosmoFac );
#ifdef ROT_FRAME
void pkdRotFrame(PKD pkd, double dOmega, double dOmegaDot);
#endif

#ifdef GASOLINE

void pkdUpdateuDot(PKD pkd, double duDelta, double dTime, double z, UHC uhc, int iGasModel, int bUpdateState );
void pkdUpdateShockTracker(PKD,double, double, double);
double pkdPoverRhoFloorJeansParticle(PKD pkd, double dResolveJeans, PARTICLE *p);
void pkdSetThermalCond(PKD pkd, struct GasPressureContext *pgpc, PARTICLE *p);
void pkdGasPressureParticle(PKD pkd, struct GasPressureContext *pgpc, PARTICLE *p, 
    double *pPoverRhoFloorJeans, double *pPoverRhoHot, double *pPoverRhoGas, double *pcGas );
void pkdGasPressure(PKD, struct GasPressureContext *pgpc);
void pkdGetDensityU(PKD, double);
void pkdLowerSoundSpeed(PKD, double);
void pkdInitEnergy(PKD pkd, double dTuFac, double z, double dTime );
void pkdKickRhopred(PKD pkd, double dHubbFac, double dDelta);
int pkdSphCurrRung(PKD pkd, int iRung, int bGreater);
void pkdSphStep(PKD pkd, double dCosmoFac, double dEtaCourant, double dEtauDot, double duMindt, double dDiffCoeff, double dEtaDiffusion, double dResolveJeans, int bViscosityLimitdt, double *pdtMinGas);
void pkdSinkStep(PKD pkd, double dtMax );
void pkdSetSphStep(PKD pkd, double dt );
void pkdSphViscosityLimiter(PKD pkd, int bOn, int bShockTracker);

void pkdDensCheck(PKD pkd, int iRung, int bGreater, int iMeasure, void *data);

#endif /* GASOLINE */
#ifdef GLASS
void pkdGlassGasPressure(PKD, void *in);
void pkdRandomVelocities(PKD pkd, double dMaxVL, double dMaxVR);
#endif

#ifdef SLIDING_PATCH
double SHEAR(int,double,PATCH_PARAMS *);
#define SHEAR(ix,t,pp)                                                  \
	((ix) < 0 ? fmod(0.5*(pp)->dLength - 1.5*(ix)*(pp)->dOrbFreq*(pp)->dWidth*(t),(pp)->dLength) - 0.5*(pp)->dLength: \
    (ix) > 0 ? 0.5*(pp)->dLength - fmod(0.5*(pp)->dLength + 1.5*(ix)*(pp)->dOrbFreq*(pp)->dWidth*(t),(pp)->dLength): 0.0)
void pkdPatch(PKD pkd);
int pkdRandAzWrap(PKD pkd);
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
    double dTime);
#endif
#endif /* COLLISIONS */

void pkdMassInR(PKD pkd, double R, double *pdMass, FLOAT *com);

#ifdef NEED_VPRED
#ifdef GASOLINE
void pkdKickVpred(PKD pkd,double dvFacOne,double dvFacTwo,double duDelta,
    int iGasModel,double z,double duDotLimit, double dTimeEnd,UHC uhc);
#else
void pkdKickVpred(PKD pkd, double dvFacOne, double dvFacTwo);
#endif
#endif

#ifdef RADIATIVEBOX
void pkdLocalFinishLWTree(PKD pkd,int iCell, double fPrevLW, FLOAT *PrevcLumLW);
void pkdFinishLWTree(PKD pkd);
#endif

void pkdInitRotBar(PKD pkd, ROTBAR rotbar);
void pkdRotatingBar(PKD pkd, double amp, /* relative amplitude of bar */
    double posang, /* position angle of bar */
    double b5,	/* radial scale length (^5) */
    FLOAT *aCom, /* Center of mass */
    double *accCom, /* acceleration (returned) */
    double *dTorque); /* acceleration (returned) */


void pkdCOM(PKD pkd, double *com);
void pkdCOMByType(PKD pkd, int type, double *com);
void pkdOldestStar(PKD pkd, double *com);
int pkdSetSink(PKD pkd, double dSinkMassMin);
void pkdFormSinks(PKD pkd, int bJeans, double dJConst2, int bDensity, double dDensityCut, double dTime, int iKickRung, int bSimple, int *nCandidates, double *Jvalmin);

void pkdSinkLogInit(PKD pkd);
void pkdSinkLogFlush(PKD pkd, char *pszFileName);
#ifdef PARTICLESPLIT
void pkdSplitGas(PKD pkd, double dInitGasMass);
#endif
#endif /* PKD_HINCLUDED */
