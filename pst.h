#ifndef PST_HINCLUDED
#define PST_HINCLUDED

#include "pkd.h"
#include "mdl.h"
#include "smoothfcn.h"
#include "floattype.h"
#include "dumpframe.h"

#ifdef COLLISIONS
#include "collision.h"
#endif
 
#ifdef AGGS
#include "aggs.h"
#endif

#ifdef RUBBLE_ZML
#include "rubble.h"
#endif

#ifdef STARFORM
#include "starform.h"
#include "feedback.h"
#include "supernova.h"
#endif

typedef struct lclBlock {
	char *pszDataPath;
	PKD	pkd;
	int nPstLvl;
	int iWtFrom;
	int iWtTo;
	int iPart;
	int iOrdSplit;
	FLOAT fSplit;
	FLOAT fWtLow;
	FLOAT fWtHigh;
	FLOAT fLow;
	FLOAT fHigh;
	int nWriteStart;
	int nDarkWriteStart;
	int nGasWriteStart;
	int nStarWriteStart;
	} LCL;

typedef struct pstContext {
	struct pstContext *pstLower;
	MDL mdl;
	LCL *plcl;
	int idSelf;
	int idUpper;
	int nLeaves;
	int nLower;
	int nUpper;
	int iLvl;
	BND bnd;
	BND bndActive;
	BND bndTreeActive;
	int iSplitDim;
	int iOrdSplit;
	FLOAT fSplit;
	FLOAT fSplitInactive;
	int nTotal;
	int nDark;
	int nGas;
	int nStar;
	/*
	 ** The PST node is also a valid cell for the tree.
	 */
	KDN kdn;
	} * PST;


#define PST_FILENAME_SIZE	512

enum pst_service {
      PST_SRV_STOP,
      PST_SETADD,
      PST_LEVELIZE,
      PST_READTIPSY,
      PST_MOVEPART,
      PST_DOMAINDECOMP,
      PST_CALCBOUND,
      PST_GASWEIGHT,
      PST_RUNGDDWEIGHT,
      PST_WEIGHT,
      PST_WEIGHTWRAP,
      PST_FREESTORE,
      PST_COLREJECTS,
      PST_SWAPREJECTS,
      PST_DOMAINCOLOR,
      PST_COLORDREJECTS,
      PST_DOMAINORDER,
      PST_LOCALORDER,
      PST_OUTARRAY,
      PST_OUTVECTOR,
      PST_OUTNCVECTOR,
      PST_WRITETIPSY,
      PST_TREEZIP,
      PST_BUILDTREE,
      PST_SMOOTH,
      PST_GRAVITY,
      PST_GRAVEXTERNAL,
      PST_CALCEANDL,
	  PST_CALCEANDLEXT,
      PST_DRIFT,
      PST_KICK,
      PST_KICKPATCH,
      PST_READCHECK,
      PST_WRITECHECK,
      PST_OUTPUTBLACKHOLES,
      PST_SETSOFT,
      PST_PHYSICALSOFT,
      PST_PREVARIABLESOFT,
      PST_POSTVARIABLESOFT,
      PST_SETTOTAL,
      PST_SETTOTALS,
      PST_CALCCELL,
      PST_COLCELLS,
      PST_DISTRIBCELLS,
      PST_CALCROOT,
      PST_DISTRIBROOT,
      PST_ONENODEREADINIT,
      PST_SWAPALL,
      PST_MASSCHECK,
      PST_MASSMETALSENERGYCHECK,
      PST_ACTIVETYPEORDER,
      PST_ACTIVEORDER,
      PST_SETRUNG,
      PST_ACTIVERUNG,
      PST_CURRRUNG,
      PST_GRAVSTEP,
      PST_ACCELSTEP,
      PST_DENSITYSTEP,
      PST_RUNGSTATS,
      PST_GETMAP,
      PST_COOLVELOCITY,
      PST_RESETTOUCHRUNG,
      PST_ACTIVEEXACTTYPE,
      PST_ACTIVETYPE,
      PST_SETTYPE,
      PST_RESETTYPE,
      PST_COUNTTYPE,
      PST_ACTIVEMASKRUNG,
      PST_ACTIVETYPERUNG,
      PST_SETPARTICLETYPES,
	  PST_SOUGHTPARTICLELIST,
	  PST_COOLUSINGPARTICLELIST,
	  PST_SETTYPEFROMFILE,
      PST_MARKSMOOTH,
      PST_RESMOOTH,
      PST_INITACCEL,
      PST_MODIFYACCEL,
      PST_DTTORUNG,
      PST_INITDT,
      PST_ORDWEIGHT,
      PST_SETWRITESTART,
      PST_SETNCWRITESTART,
      PST_COLNPARTS,
      PST_NEWORDER,
      PST_GETNPARTS,
      PST_SETNPARTS,
      PST_UPDATEUDOT,
      PST_UPDATESHOCKTRACKER,
      PST_GETGASPRESSURE,
      PST_GETDENSITYU,
      PST_LOWERSOUNDSPEED,
      PST_INITENERGY,
      PST_KICKVPRED,
      PST_KICKRHOPRED,
      PST_SPHCURRRUNG,
      PST_RANDOMVELOCITIES,
      PST_DENSCHECK,
	  PST_GETSPECIALPARTICLES,
	  PST_DOSPECIALPARTICLES,
      PST_GHOSTEXCLUDE,
	  /* following for COLLISIONS */
      PST_NUMREJECTS,
      PST_READSS,
      PST_WRITESS,
	  PST_KICKUNIFGRAV,
      PST_NEXTENCOUNTER,
      PST_MARKENCOUNTERS,
      PST_NEXTCOLLISION,
	  PST_GETCOLLIDERINFO,
      PST_DOCOLLISION,
	  PST_RESETCOLLIDERS,
      PST_FINDBINARY,
      PST_MERGEBINARY,
	  /* following for OLD_KEPLER */
#ifdef OLD_KEPLER
      PST_QQCALCBOUND,
      PST_QQDOMAINDECOMP,
      PST_QQBUILDTREE,
      PST_QQSMOOTH,
#endif
#ifdef SLIDING_PATCH
      PST_FINDLM,
      PST_GETNEIGHBORS,
#endif
	  /* following for AGGS */
	  PST_AGGSFIND,
	  PST_AGGSCONFIRM,
	  PST_AGGSMERGE,
	  PST_AGGSBACKDRIFT,
	  PST_AGGSGETCOM,
	  PST_AGGSGETAXESANDSPIN,
	  PST_AGGSSETBODYPOS,
	  PST_AGGSSETSPACEPOS,
	  PST_AGGSSETSPACEVEL,
	  PST_AGGSSETSPACESPINS,
	  PST_AGGSDELETE,
	  PST_AGGSGETACCEL,
	  PST_AGGSCHECKSTRESS,
	  PST_AGGSGETTORQUE,
	  /* following for SLIDING_PATCH */
      PST_RANDAZWRAP,
	  /* following for RUBBLE_ZML */
	  PST_DUSTBINSGETMASS,
	  PST_DUSTBINSAPPLY,
	  PST_RUBBLERESETCOLFLAG,
	  PST_RUBBLECHECKFORKDKRESTART,
	  PST_RUBBLESTEP,
	  PST_RUBCLEANUP,
	  PST_RUBINTERPCLEANUP,
	  /* following for SPH, etc. */
      PST_SPHSTEP,
      PST_SINKSTEP,
      PST_SETSPHSTEP,
      PST_SETBALL,
      PST_SPHVISCOSITYLIMITER,
      PST_INITCOOLING,
      PST_COOLTABLEREAD,
      PST_GROWMASS,
      PST_CLEARTIMER,
      PST_MASSINR,
      PST_ROTBARINIT,
      PST_BALLMAX,
      PST_FORMSINKS,
      PST_FORMSTARS,
      PST_FEEDBACK,
      PST_GRAVINFLOW,
      PST_CREATEINFLOW,
      PST_SIMPLESTARFORM,
	  PST_SSFCREATESTARS,
      PST_DUMPFRAME,
      PST_DUMPVOXEL,
      PST_COM,
	  PST_COMBYTYPE,
	  PST_OLDESTSTAR,
	  PST_SETSINK,
      PST_INITSTARLOG,
      PST_FLUSHSTARLOG
      };

void pstAddServices(PST,MDL);
void pstInitialize(PST *,MDL,LCL *);
void pstFinish(PST);

/* PST_SETADD */
struct inSetAdd {
	int id;
	};
void pstSetAdd(PST,void *,int,void *,int *);

/* PST_LEVELIZE */
struct inLevelize {
	int iLvl;
	};
void pstLevelize(PST,void *,int,void *,int *);

/* PST_READTIPSY */
struct inReadTipsy {
	int nFileStart;
	int nFileEnd;
	int nDark;	
	int nGas;
	int nStar;
	int iOrder;
	int nIOProcessor;
	float fExtraStore;
	FLOAT fPeriod[3];
        FLOAT dxInflow, dxOutflow;
	int bStandard;
	int iReadIOrder;
	double dvFac;
	double dTuFac;
	char achInFile[PST_FILENAME_SIZE];
	};
void pstReadTipsy(PST,void *,int,void *,int *);

/* PST_DOMAINDECOMP */
struct inDomainDecomp {
    int bDoRootFind;
    int bDoSplitDimFind;
    int bSplitWork;
    };

void pstDomainDecomp(PST,void *,int,void *,int *);

/* PST_CALCBOUND */
struct outCalcBound {
	BND bnd;
	BND bndActive;
	BND bndTreeActive;
	BND bndBall;
	};
void pstCalcBound(PST,void *,int,void *,int *);

void pstGasWeight(PST,void *,int,void *,int *);

struct inRungDDWeight {
	int iMaxRung;
        double dWeight; 
        };

void pstRungDDWeight(PST,void *,int,void *,int *);

/* PST_WEIGHT */
struct inWeight {
	int iSplitDim;
	FLOAT fSplit;
	int iSplitSide;
	int ittr;
	int pFlag;
	};
struct outWeight {
	int nLow;
	int nHigh;
	FLOAT fLow;
	FLOAT fHigh;
	};
void pstWeight(PST,void *,int,void *,int *);
struct inWeightWrap {
	int iSplitDim;
	FLOAT fSplit;
	FLOAT fSplit2;
	int iSplitSide;
	int ittr;
	int pFlag;
	};
void pstWeightWrap(PST,void *,int,void *,int *);

/* PST_FREESTORE */
struct outFreeStore {
	int nFreeStore;
	};
void pstFreeStore(PST,void *,int,void *,int *);

/*
 ** This structure is used by reject collectors and SwapRejects
 */
typedef struct outReject {
	int id;
	int nRejects;
	int nSpace;
	int nLocal;
	} OREJ;

/* PST_COLREJECTS */
struct inColRejects {
	int iSplitDim;
	FLOAT fSplit;
	FLOAT fSplitInactive;
	int iSplitSide;
	};
void pstColRejects(PST,void *,int,void *,int *);

/* PST_SWAPREJECTS */
void pstSwapRejects(PST,void *,int,void *,int *);

/* PST_DOMAINCOLOR */
void pstDomainColor(PST,void *,int,void *,int *);

/* PST_COLORDREJECTS */
struct inColOrdRejects {
	int iOrdSplit;
	int iSplitSide;
	};
void pstColOrdRejects(PST,void *,int,void *,int *);

/* PST_DOMAINORDER */
struct inDomainOrder {
	int iMaxOrder;
	};
void pstDomainOrder(PST,void *,int,void *,int *);

/* PST_LOCALORDER */
void pstLocalOrder(PST,void *,int,void *,int *);

/* PST_OUTARRAY */
struct inOutArray {
	char achOutFile[PST_FILENAME_SIZE];
        int nStart;
	int iType;
	int iBinaryOutput;
	int bStandard;
	};
void pstOutArray(PST,void *,int,void *,int *);

/*PST_OUTNCVECTOR*/
struct outNC {
    float min[3][3];
    float max[3][3];
    };
void pstOutNCVector(PST,void *,int,void *,int *);

/* PST_OUTVECTOR */
struct inOutput {
	char achOutFile[PST_FILENAME_SIZE];
	int iDim;
	int iType;
	int iBinaryOutput;
	int N;
	int bStandard;
	int nIOProcessor;
	double duTFac;
	double dvFac;
	};
void pstOutVector(PST,void *,int,void *,int *);

/* PST_WRITETIPSY */
struct inWriteTipsy {
	int bStandard;
	int nIOProcessor;
	double dvFac;
	double duTFac;
	int iGasModel;
	char achOutFile[PST_FILENAME_SIZE];
	};
void pstWriteTipsy(PST,void *,int,void *,int *);

struct inTreeZip {
        char achFile[PST_FILENAME_SIZE];
        int bLocal;
        BND bnd;
        };
void pstTreeZip(PST,void *,int,void *,int *);

/* PST_BUILDTREE */
struct inBuildTree {
	int nBucket;
	int iOpenType;
	int iOrder;
	double dCrit;
	int bBinary;
        int bActiveOnly;
	int bTreeActiveOnly;
	int bGravity;
	};
struct outBuildTree {
	KDN kdn;
	};
void pstBuildTree(PST,void *,int,void *,int *);

/* PST_SMOOTH */
struct inSmooth {
	int nSmooth;
        int bPeriodic;
	int bSymmetric;
	int iSmoothType;
        double dfBall2OverSoft2;
	SMF smf;
	};
struct outSmooth {
    int iSmoothFlags;  /* Warning Flags for need to smooth again, etc... */
	/*
	 ** Cache Statistics.
	 */
	double dpASum;
	double dpMSum;
	double dpCSum;
	double dpTSum;
	double dcASum;
	double dcMSum;
	double dcCSum;
	double dcTSum;
	};
void pstSmooth(PST,void *,int,void *,int *);

/* PST_MOVEPART */
struct inMoveParticle {
    PARTICLE p;
    FLOAT dOrigin[3]; /* Origin of new frame */
    FLOAT dRelx[3]; /* Displacement from newx */
#ifdef SLIDING_PATCH
  PATCH_PARAMS PP;
#endif
    };
void pstMoveParticle(PST,void *,int,void *,int *);

/* PST_GRAVITY */
struct inGravity {
	int nReps;
	int bPeriodic;
	int iOrder;
	int bEwald;
	int iEwOrder;
	int bComove;
	int bDoSun;
	double dSunSoft;
	double dEwCut;
	double dEwhCut;
	double dRhoFac;
#ifdef SLIDING_PATCH
  double dTime;
  PATCH_PARAMS PP;
#endif
	};
struct outGravity {
	int nActive;
	int nTreeActive;
    double aSun[3];
	double dPartSum;
	double dCellSum;
	double dSoftSum;
	double dFlop;
	/*	
	 ** Collected CPU time stats.
	 */
	double dWSum;
	double dWMax;
	double dWMin;
	double dISum;
	double dIMax;
	double dIMin;
	double dESum;
	double dEMax;
	double dEMin;
	/*
	 ** Cache Statistics.
	 */
	double dpASum;
	double dpMSum;
	double dpCSum;
	double dpTSum;
	double dcASum;
	double dcMSum;
	double dcCSum;
	double dcTSum;
	};
void pstGravity(PST,void *,int,void *,int *);

/* PST_GRAVEXTERNAL */
struct inGravExternal {
    int bIndirect;
    int bDoSun;
	double dTime;
    double dSunMass;
	double dSunSoft;
    double aSun[3];
	/*
	 ** For external galaxy potential stuff
	 */
	int bLogHalo;
	int bHernquistSpheroid;
    int bNFWSpheroid;
        double dNFWm200;
        double dNFWr200;
        double dNFWconc;
        double dNFWsoft;
    int bElliptical;
    int bEllipticalDarkNFW;
	int bHomogSpheroid;
	int bBodyForce;
	double dBodyForceConst;
	int bMiyamotoDisk;
	int bTimeVarying;
	int bRotatingBar;
	double dRotBarAmp;
	double dRotBarPosAng;
	double dRotBarB5;
    FLOAT aCom[3];
    
#ifdef ROT_FRAME
	int bRotFrame;
	double dOmega;
	double dOmegaDot;
#endif
#ifdef SLIDING_PATCH
  int bPatch;
  PATCH_PARAMS PP;
#endif
#ifdef SIMPLE_GAS_DRAG
	int bSimpleGasDrag;
	int iFlowOpt;
	int bEpstein;
	double dGamma;
	double dTime;
#endif
	};
struct outGravExternal {
    double dAcc[3];
    double dTorque[3];
    };

void pstGravExternal(PST,void *,int,void *,int *);

/* PST_CALCEANDL */
struct outCalcEandL {
	double T;
	double U;
	double Eth;
	double L[3];
	};
void pstCalcEandL(PST,void *,int,void *,int *);

/* PST_CALCEANDLEXT */
struct inCalcEandLExt {
    int bHeliocentric;
	};
struct outCalcEandLExt {
	double dMass;
    double dSumMR[3];
    double dSumMV[3];
	double dPot;
	};
void pstCalcEandLExt(PST,void *,int,void *,int *);

/* PST_DRIFT */
struct inDrift {
	double dDelta;
	FLOAT fCenter[3];
	int bPeriodic;
        int bInflowOutflow;
	int bFandG;
	FLOAT fCentMass;
        double dTime;
#ifdef SLIDING_PATCH
  PATCH_PARAMS PP;
#endif
	};
void pstDrift(PST,void *,int,void *,int *);

/* PST_UPDATEUDOT */
struct inUpdateuDot {
	double duDelta;
	double dTime;	
	double z;
	int iGasModel;
	int bUpdateState;
	};
struct outUpdateuDot {
	double Time;
	double MaxTime;
	double SumTime;
	int nSum;
	};

void pstUpdateuDot(PST,void *,int,void *,int *);

/* PST_KICK */
struct inKick {
	double dvFacOne;
	double dvFacTwo;
	double dvPredFacOne;
	double dvPredFacTwo;
	double duDelta;
	double duPredDelta;
	double duDotLimit;
	int iGasModel;
        double z;
        double dTimeEnd;
	};
struct outKick {
	double Time;
	double MaxTime;
	double SumTime;
	int nSum;
	};

void pstKick(PST,void *,int,void *,int *);

/* PST_KICKPATCH */
struct inKickPatch {
    int bOpen;
    double dOrbFreq;
    double dvFacOne;
    double dvFacTwo;
    };
void pstKickPatch(PST,void *,int,void *,int *);

/* PST_READCHECK */
struct inReadCheck {
	int iVersion;
	int iOffset;
	int nFileStart;
	int nFileEnd;
	int nDark;
	int nGas;
	int nStar;
	int iOrder;
	int nIOProcessor;
	float fExtraStore;
	FLOAT fPeriod[3];
        FLOAT dxInflow, dxOutflow;
	char achInFile[PST_FILENAME_SIZE];
	};
void pstReadCheck(PST,void *,int,void *,int *);

/* PST_WRITECHECK */
struct inWriteCheck {
	int iOffset;
	int nIOProcessor;
	char achOutFile[PST_FILENAME_SIZE];
	};
void pstWriteCheck(PST,void *,int,void *,int *);

/* PST_OUTPUTBLACKHOLES */
void pstOutputBlackHoles(PST pst,void *vin,int nIn,void *vout,int *pnOut);

/* PST_SETSOFT */
struct inSetSoft {
	double dSoft;
	};
void pstSetSoft(PST,void *,int,void *,int *);

#ifdef CHANGESOFT
/* PST_PHYSICALSOFT */
struct inPhysicalSoft {
        double dSoftMax;
        double dFac;
        int bSoftMaxMul;
        };
void pstPhysicalSoft(PST,void *,int,void *,int *);

/* PST_PREVARIABLESOFT */
struct inPreVariableSoft {
        int iVariableSoftType;
        };
void pstPreVariableSoft(PST,void *,int,void *,int *);

/* PST_POSTVARIABLESOFT */
struct inPostVariableSoft {
        double dSoftMax;
        int bSoftMaxMul;
        int iVariableSoftType;
        };
void pstPostVariableSoft(PST,void *,int,void *,int *);
#endif

/* PST_SETTOTAL */
struct outSetTotal {
	int nTotal;
	};
void pstSetTotal(PST,void *,int,void *,int *);

/* PST_SETTOTALS */
struct outSetTotals {
	int nGas;
	int nStar;
	int nDark;
	};
void pstSetTotals(PST,void *,int,void *,int *);

/* PST_CALCCELL */
struct inCalcCell {
	int iOrder;
	FLOAT rcm[3];
	};
struct outCalcCell {
	struct pkdCalcCellStruct mom;
	};
void pstCalcCell(PST,void *,int,void *,int *);

/* PST_COLCELLS */
struct inColCells {
	int iCell;
	int nCell;
	};
void pstColCells(PST,void *,int,void *,int *);

/* PST_DISTRIBCELLS */
void pstDistribCells(PST,void *,int,void *,int *);

/* PST_CALCROOT */
struct ioCalcRoot {
	struct ilCellNewt ilcn;
	};
void pstCalcRoot(PST,void *,int,void *,int *);

/* PST_DISTRIBROOT */
void pstDistribRoot(PST,void *,int,void *,int *);

/* PST_ONENODEREADINIT */
void pstOneNodeReadInit(PST pst,void *vin,int nIn,void *vout,int *pnOut);

/* PST_SWAPALL */
void pstSwapAll(PST pst,void *vin,int nIn,void *vout,int *pnOut);

/* PST_MASSCHECK */
struct outMassCheck {
	double dMass;
	};
void pstMassCheck(PST,void *,int,void *,int *);

struct outMassMetalsEnergyCheck {
	double dTotMass;
        double dTotMetals;
        double dTotOx;
        double dTotFe;
        double dTotEnergy;
	};
void pstMassMetalsEnergyCheck(PST,void *,int,void *,int *);

struct inActiveTypeOrder {
        unsigned int iTestMask;
        };

void pstActiveTypeOrder(PST,void *,int,void *,int *);

/* PST_ACTIVEORDER */
void pstActiveOrder(PST,void *,int,void *,int *);

/* PST_SETRUNG */
struct inSetRung {
    int iRung;
    };
void pstSetRung(PST,void *,int,void *,int *);

struct inDensCheck {
    int iRung;
    int bGreater;
    int iMeasure;
    };
struct outDensCheck {
    double dMaxDensError;
    double dAvgDensError;
    int nError;
    int nTotal;
    };

void pstDensCheck(PST,void *,int,void *,int *);

struct inBallMax {
    int iRung;
    int bGreater;
    double dhFac;
    };

void pstBallMax(PST,void *,int,void *,int *);

/* PST_ACTIVERUNG */
struct inActiveRung {
    int iRung;
    int bGreater;
    };

void pstActiveRung(PST,void *,int,void *,int *);

/* PST_CURRRUNG */
struct inCurrRung {
    int iRung;
    };
struct outCurrRung {
    int iCurrent;
    };
void pstCurrRung(PST,void *,int,void *,int *);

/* PST_GRAVSTEP */
struct inGravStep {
    double dEta;
    };
void pstGravStep(PST,void *,int,void *,int *);

/* PST_ACCELSTEP */
struct inAccelStep {
    double dEta;
    double dVelFac;
    double dAccFac;
    int    bDoGravity;
    int    bEpsAcc;
    int    bSqrtPhi;
    double dhMinOverSoft;
    };
void pstAccelStep(PST,void *,int,void *,int *);

/* PST_DENSITYSTEP */
struct inDensityStep {
    double dEta;
    double dRhoFac;
    };
void pstDensityStep(PST,void *,int,void *,int *);

/* PST_RUNGSTATS */
struct inRungStats {
	int iRung;
	};
struct outRungStats {
	int nParticles;
	};
void pstRungStats(PST,void *,int,void *,int *);

/* PST_GETMAP */
struct inGetMap {
	int nStart;
	};
void pstGetMap(PST,void *,int,void *,int *);

/* PST_COOLVELOCITY */
struct inCoolVelocity {
	int nSuperCool;
	double dCoolFac;
	double dCoolDens;
	double dCoolMaxDens;
	};
void pstCoolVelocity(PST,void *,int,void *,int *);

struct inActiveType {
	unsigned int iFilterMask;
	unsigned int iTestMask;
	unsigned int iSetMask;
	int iRung;
	int bGreater;
	};

void pstActiveExactType(PST,void *,int,void *,int *);
void pstActiveType(PST,void *,int,void *,int *);
void pstSetType(PST,void *,int,void *,int *);
void pstResetType(PST,void *,int,void *,int *);
void pstCountType(PST,void *,int,void *,int *);
void pstActiveMaskRung(PST,void *,int,void *,int *);
void pstActiveTypeRung(PST,void *,int,void *,int *);

#define PST_SETTYPEFROMFILEMAXLEN 160
struct inSetTypeFromFile {
  int iSetMask;
  int biGasOrder;
  char file[PST_SETTYPEFROMFILEMAXLEN];
};

struct outSetTypeFromFile {
  int niOrder;
  int nSet;
  int nSetiGasOrder;
};

void pstSetTypeFromFile(PST,void *,int,void *,int *);

struct inSetParticleTypes {
	int nSuperCool;
	};

void pstSetParticleTypes(PST,void *,int,void *,int *);

#define MAXSOUGHTPARTICLELIST 10

struct inSoughtParticleList {
	int iTypeSought;
    int nMax;
	};

struct inoutParticleList {
    int n;
    struct SoughtParticle p[MAXSOUGHTPARTICLELIST];
	};

void pstSoughtParticleList(PST,void *,int,void *,int *);
void pstCoolUsingParticleList(PST,void *,int,void *,int *);

/* PST_RESMOOTH */
struct inMarkSmooth {
	int nSmooth;
	int bPeriodic;
	int bSymmetric;
	int iSmoothType;
	int iMarkType;
	SMF smf;
	};
void pstMarkSmooth(PST,void *,int,void *,int *);

struct inReSmooth {
	int nSmooth;
	int bPeriodic;
	int bSymmetric;
	int iSmoothType;
        double dfBall2OverSoft2;
	SMF smf;
	};
struct outReSmooth {
    int iSmoothFlags;  /* Warning Flags for need to smooth again, etc... */
	/*
	 ** Cache Statistics.
	 */
	double dpASum;
	double dpMSum;
	double dpCSum;
	double dpTSum;
	double dcASum;
	double dcMSum;
	double dcCSum;
	double dcTSum;
	};
void pstReSmooth(PST,void *,int,void *,int *);

/* PST_INITACCEL */
void pstInitAccel(PST,void *,int,void *,int *);

struct inModifyAccel {
    double a;
    };

void pstModifyAccel(PST,void *,int,void *,int *);

/* PST_DTTORUNG */
struct inDtToRung {
    int iRung;
    double dDelta;
    int iMaxRung;
    int bAll;
    };
struct outDtToRung {
    int iMaxRung;
    int nMaxRung;
    int iMaxRungIdeal;
    };
void pstDtToRung(PST,void *,int,void *,int *);

/* PST_INITDT */
struct inInitDt {
    double dDelta;
    };
void pstInitDt(PST,void *,int,void *,int *);

/* PST_ORDWEIGHT */
struct inOrdWeight {
	int iOrdSplit;
	int iSplitSide;
	int ittr;
	};
struct outOrdWeight {
	int nLow;
	int nHigh;
	};
void pstOrdWeight(PST,void *,int,void *,int *);

/* PST_SETWRITESTART */
struct inSetWriteStart {
	int nWriteStart;
	};
void pstSetWriteStart(PST,void *,int,void *,int *);

/* PST_SETNCWRITESTART */
struct inSetNCWriteStart {
	int nDarkWriteStart;
	int nGasWriteStart;
	int nStarWriteStart;
	};
void pstSetNCWriteStart(PST,void *,int,void *,int *);

/* PST_COLNPARTS */
struct outColNParts {
    int nNew;
    int nAddGas;
    int nAddDark;
    int nAddStar;
    int nDelGas;
    int nDelDark;
    int nDelStar;
    };
void pstColNParts(PST, void *, int, void *, int *);

/* PST_NEWORDER */
struct inNewOrder {
    int Gas, Dark, Star;
    };
void pstNewOrder(PST, void *, int, void *, int *);

/* PST_GETNPARTS */
/* see pkd.h
 struct outGetNParts { 
	int n;
    int nGas;
    int nDark;
    int nStar;
    int iMaxOrderGas;
    int iMaxOrderDark;
    int iMaxOrderStar;
    };
*/
void pstGetNParts(PST, void *, int, void *, int *);

/* PST_SETNPARTS */
struct inSetNParts {
    int nGas;
    int nDark;
    int nStar;
	int nMaxOrder;
    int nMaxOrderGas;
    int nMaxOrderDark;
    };
void pstSetNParts(PST, void *, int, void *, int *);

#ifdef GASOLINE

struct inGetGasPressure {
	enum GasModel iGasModel; 
  /* Adiabatic */
	double gamma;
	double gammam1;
	double dResolveJeans;
	double dCosmoFac;
    
  /* Isothermal */

  /* Ion evolving */


#ifdef GLASS
    struct GlassData g;
#endif
	};



/* PST_GETGASPRESSURE */
void pstGetGasPressure(PST, void *,int,void *,int *);

void pstGetDensityU(PST, void *,int,void *,int *);

struct inLowerSoundSpeed {
	double dhMinOverSoft;
	};

void pstLowerSoundSpeed(PST, void *,int,void *,int *);

/* See cooling_*.h for details on input structs */
void pstInitEnergy(PST, void *,int,void *,int *);
void pstInitCooling(PST,void *,int,void *,int *);
void pstCoolTableRead(PST,void *,int,void *,int *);

/* PST_KICKVPRED */
struct inKickRhopred {
	double dHubbFac;
	double dDelta;
	};
void pstKickRhopred(PST,void *,int,void *,int *);

/* PST_SPHCURRRUNG */
struct inSphCurrRung {
    int iRung;
    int bGreater;
    };
struct outSphCurrRung {
    int iCurrent;
    };
void pstSphCurrRung(PST,void *,int,void *,int *);

#endif

#ifdef GLASS
/* PST_RANDOMVELOCITIES */
struct inRandomVelocities {
    double dMaxVelocityL;
    double dMaxVelocityR;
    };
void pstRandomVelocities(PST,void *,int,void *,int *);
#endif

#ifdef SPECIAL_PARTICLES

#include "special.h"

/* PST_GETSPECIALPARTICLES */
struct inGetSpecial {
	int nSpecial;
	int iId[MAX_NUM_SPECIAL_PARTICLES];
	SPECIAL_MASTER_INFO mInfo;
	};
struct outGetSpecial {
	SPECIAL_PARTICLE_INFO sInfo[MAX_NUM_SPECIAL_PARTICLES];
	};
void pstGetSpecialParticles(PST,void *,int,void *,int *);

/* PST_DOSPECIALPARTICLES */
struct inDoSpecial {
	int nSpecial;
	SPECIAL_MASTER_INFO mInfo;
	SPECIAL_PARTICLE_DATA sData[MAX_NUM_SPECIAL_PARTICLES];
	SPECIAL_PARTICLE_INFO sInfo[MAX_NUM_SPECIAL_PARTICLES];
	};
struct outDoSpecial {
	FLOAT aFrame[3];
	};
void pstDoSpecialParticles(PST,void *,int,void *,int *);

#endif

#ifdef COLLISIONS

struct inSetBall {
  double dDelta;
  double dBallVelFact;
};
void pstSetBall(PST,void *,int,void *,int *);

/* PST_NUMREJECTS */
struct outNumRejects {
	int nRej;
	};
void pstNumRejects(PST,void *,int,void *,int *);

/* PST_READSS */
struct inReadSS {
	int nFileStart;
	int nFileEnd;
	int nDark;
	int nGas;			/* always zero */
	int nStar;			/* always zero */
	int iOrder;
	float fExtraStore;
	FLOAT fPeriod[3];	/* for compatability */
        FLOAT dxInflow, dxOutflow;
	char achInFile[PST_FILENAME_SIZE];
	};
void pstReadSS(PST,void *,int,void *,int *);

/* PST_WRITESS */
struct inWriteSS {
	char achOutFile[PST_FILENAME_SIZE];
	};
void pstWriteSS(PST,void *,int,void *,int *);

/* #define PST_FINDENCOUNTER	67 */
struct outFindEncounter {
	double dNext;
	};
void pstFindEncounter(PST,void *,int,void *,int *);

/* PST_KICKUNIFGRAV */
struct inKickUnifGrav {
	double dvx,dvy,dvz;
	};
void pstKickUnifGrav(PST,void *,int,void *,int *);

/* #define PST_MARKENCOUNTERS	68 */
/* PST_NEXTENCOUNTER */
struct outNextEncounter {
	double dt;
	};
void pstNextEncounter(PST,void *,int,void *,int *);

/* PST_MARKENCOUNTERS */

struct inMarkEncounters {
	double dt;
	};
void pstMarkEncounters(PST,void *,int,void *,int *);

/* PST_NEXTCOLLISION */
struct outNextCollision {
	double dt;
	int iOrder1,iOrder2;
	};
void pstNextCollision(PST,void *,int,void *,int *);

/* PST_GETCOLLIDERINFO */
struct inGetColliderInfo {
	int iOrder;
	};
struct outGetColliderInfo {
	COLLIDER Collider;
	};
void pstGetColliderInfo(PST,void *,int,void *,int *);

/* PST_DOCOLLISION */
struct inDoCollision {
	double dt;
	COLLIDER Collider1,Collider2;
	int bPeriodic;
	COLLISION_PARAMS CP;
	double dBaseTime;
	int iTime0;
#ifdef SLIDING_PATCH
  PATCH_PARAMS PP;
#endif
	double dTime;
#ifdef AGGS
	int iAggNewIdx;
#endif
	};
struct outDoCollision {
	COLLIDER Out[MAX_NUM_FRAG];
	double dT;
	int iOutcome,nOut;
	};
void pstDoCollision(PST,void *,int,void *,int *);

/* PST_RESETCOLLIDERS */
struct inResetColliders {
	int iOrder1,iOrder2;
	};
void pstResetColliders(PST,void *,int,void *,int *);

/* PST_FINDTIGHTESTBINARY */
struct inBinary {
        int n;
        };
struct outBinary {
        int n;
        double dBindEn;
        int iOrder1;
        int iOrder2;
        };
void pstFindTightestBinary(PST,void *,int,void *,int *);

/* PST_MERGEBINARY */
struct inMrgBnry {
       COLLIDER c1;
       COLLIDER c2;
       int bPeriodic;
       double dBaseStep;
       double dTimeNow;
       double dDensity;
#ifdef SLIDING_PATCH
  PATCH_PARAMS PP;
       double dTime;
#endif
       int iTime0;
       int iPad; /* Cheat to try to make the code owrk in Linux-ia64 */    
       };
struct outMrgBnry {
       COLLIDER cOut;
       int n;
       };
void pstMergeBinary(PST,void *,int,void *,int *);

#ifdef SLIDING_PATCH
struct inLargeMass {
    double dMass;
    FLOAT fNumHillSphere;
    double dCentMass;
    double dOrbRad;    
    };
struct outLargeMass
{
       PARTICLE p[MAXLARGEMASS];           /* large mass particles */
       double dRadius[MAXLARGEMASS];        /* radius of sphere to be exchanged */    
       int n;                 /* # of large masses */
    };
void pstFindLargeMasses(PST pst,void *vin,int nIn,void *vout,int *pnOut);

struct inGetNeighbors
{
    FLOAT x[3];
    double dDist;
#ifdef SLIDING_PATCH
  PATCH_PARAMS PP;
    double dTime;
#endif
    int id;
    };
struct outGetNeighbors 
{
    PARTICLE p[MAXNEIGHBORS];
    double dSep2[MAXNEIGHBORS];    
    int n;
    };
void pstGetNeighborParticles(PST pst,void *vin,int nIn,void *vout,int *pnOut);

#endif /* SLIDING_PATCH */

#ifdef OLD_KEPLER
/* PST_QQCALCBOUND */
void pstQQCalcBound(PST,void *,int,void *,int *);

/* PST_QQDOMAINDECOMP */
void pstQQDomainDecomp(PST,void *,int,void *,int *);

/* PST_QQBUILDTREE */
void pstQQBuildTree(PST,void *,int,void *,int *);

/* PST_QQSMOOTH */
void pstQQSmooth(PST,void *,int,void *,int *);
#endif

#endif /* COLLISIONS */

#ifdef AGGS

/* PST_AGGSFIND */
struct outAggsFind {
	int iMaxIdx;
	};
void pstAggsFind(PST pst,void *vIn,int nIn,void *vout,int *pnOut);

/* PST_AGGSCONFIRM */
struct inAggsConfirm {
	int iAggIdx;
	};
struct outAggsConfirm {
	int bAssigned;
	};
void pstAggsConfirm(PST pst,void *vIn,int nIn,void *vout,int *pnOut);

/* PST_AGGSMERGE */
struct inAggsMerge {
	int iOldIdx,iNewIdx;
	};
void pstAggsMerge(PST pst,void *vIn,int nIn,void *vout,int *pnOut);

/* PST_AGGSBACKDRIFT */
struct inAggsBackDrift {
	int iAggIdx;
	double dt;
	};
void pstAggsBackDrift(PST pst,void *vIn,int nIn,void *vout,int *pnOut);

/* PST_AGGSGETCOM */
struct inAggsGetCOM {
	int iAggIdx;
	};
struct outAggsGetCOM {
	Scalar m;
	Vector mr,mv; /* position & velocity moments */
	};
void pstAggsGetCOM(PST pst,void *vIn,int nIn,void *vout,int *pnOut);
 
/* PST_AGGSGETAXESANDSPIN */
struct inAggsGetAxesAndSpin {
	int iAggIdx;
	Vector r_com,v_com;
	};
struct outAggsGetAxesAndSpin {
	Matrix I;
	Vector L;
	};
void pstAggsGetAxesAndSpin(PST pst,void *vIn,int nIn,void *vout,int *pnOut);
 
/* PST_AGGSSETBODYPOS */
struct inAggsSetBodyPos {
	int iAggIdx;
	Matrix spaceToBody; /* transpose of lambda */
	};
void pstAggsSetBodyPos(PST pst,void *vIn,int nIn,void *vout,int *pnOut);

/* PST_AGGSSETSPACEPOS */
struct inAggsSetSpacePos {
	int iAggIdx;
	Vector r_com;
	Matrix lambda;
	};
void pstAggsSetSpacePos(PST pst,void *vIn,int nIn,void *vout,int *pnOut);

/* PST_AGGSSETSPACEVEL */
struct inAggsSetSpaceVel {
	int iAggIdx;
	Vector v_com,omega;
	Matrix lambda;
	};
void pstAggsSetSpaceVel(PST pst,void *vIn,int nIn,void *vout,int *pnOut);

/* PST_AGGSSETSPACESPINS */
struct inAggsSetSpaceSpins {
	int iAggIdx;
	Vector omega;
	};
void pstAggsSetSpaceSpins(PST pst,void *vIn,int nIn,void *vout,int *pnOut);

/* PST_AGGSDELETE */
struct inAggsDelete {
	int iAggIdx;
	};
struct outAggsDelete {
	int bFound;
	};
void pstAggsDelete(PST pst,void *vIn,int nIn,void *vout,int *pnOut);

/* PST_AGGSGETACCEL */
struct inAggsGetAccel {
	int iAggIdx;
	};
struct outAggsGetAccel {
	Scalar m;
	Vector ma; /* acceleration moment */
	};
void pstAggsGetAccel(PST pst,void *vIn,int nIn,void *vout,int *pnOut);

/* PST_AGGSCHECKSTRESS */
struct inAggsCheckStress {
	int iAggIdx;
	Vector r_com,a_com,omega;
	FLOAT fTensileStrength,fShearStrength;
	};
struct outAggsCheckStress {
	int nLost,nLeft;
	};
void pstAggsCheckStress(PST pst,void *vIn,int nIn,void *vout,int *pnOut);

/* PST_AGGSGETTORQUE */
struct inAggsGetTorque {
	int iAggIdx;
	Vector r_com,a_com;
	};
struct outAggsGetTorque {
	Vector torque;
	};
void pstAggsGetTorque(PST pst,void *vIn,int nIn,void *vout,int *pnOut);
 
#endif /* AGGS */

#ifdef SLIDING_PATCH

/* PST_RANDAZWRAP */
struct inRandAzWrap {
  PATCH_PARAMS PP;
};
struct outRandAzWrap {
	int nRandomized;
};
void pstRandAzWrap(PST,void *,int,void *,int *);

#endif /* SLIDING_PATCH */

#ifdef RUBBLE_ZML

struct inDustBinsGetMass {
	DUST_BINS_PARAMS DB;
	DustBins aDustBins[DUST_BINS_MAX_NUM];
	double dTimeInt;
	double dCentMass;
	};
struct outDustBinsGetMass {
	double aDustBinsMassLoss[DUST_BINS_MAX_NUM];
	};
void pstDustBinsGetMass(PST pst,void *vIn,int nIn,void *vout,int *pnOut);

struct inDustBinsApply {
	double dCentMass;
	double aMassIncrFrac[DUST_BINS_MAX_NUM];
	int nBins;
	};
void pstDustBinsApply(PST pst,void *vIn,int nIn,void *vout,int *pnOut);

void pstRubbleResetColFlag(PST pst,void *vIn,int nIn,void *vout,int *pnOut);

struct outRubbleCheckForKDKRestart {
	int bRestart;
	};
void pstRubbleCheckForKDKRestart(PST pst,void *vIn,int nIn,void *vout,int *pnOut);

struct inRubbleStep {
	double dMaxStep,dMinStep;
	};
void pstRubbleStep(PST pst,void *vIn,int nIn,void *vout,int *pnOut);

struct inRubCleanup {
	int iColor;
	DUST_BINS_PARAMS DB;
	double dCentMass;
	};
struct outRubCleanup {
	double aDustBinsMassGain[DUST_BINS_MAX_NUM];
	double dDustBinsRubTrash;
	};
void pstRubCleanup(PST pst,void *vIn,int nIn,void *vout,int *pnOut);

struct inRubInterpCleanup {
	DUST_BINS_PARAMS DB;
	double dCentMass;
	int iOrder;
	};
struct outRubInterpCleanup {
	double dDustBinsInterpMass;
	int iBin;
	};
void pstRubInterpCleanup(PST pst,void *vIn,int nIn,void *vout,int *pnOut);

#endif /* RUBBLE_ZML */

struct inGravInflow
    {
    double r;
    };

void pstGravInflow(PST,void *,int,void *,int *);

struct inCreateInflow
    {
    int Ny,iGasModel,iRung;
    double dTuFac;
    double pmass, x, vx, density, temp, metals, eps;
    double dt;
    };

void pstCreateInflow(PST,void *,int,void *,int *);

#ifdef GASOLINE

/* PST_SPHSTEP */
struct inSphStep {
    double dCosmoFac;
    double dEtaCourant;
    double dEtauDot;
    int bViscosityLimitdt;
    };
void pstSphStep(PST,void *,int,void *,int *);

struct inSinkStep {
    double dtMax;
    };
void pstSinkStep(PST,void *,int,void *,int *);

struct inSetSphStep {
    double dt;
    };
void pstSetSphStep(PST,void *,int,void *,int *);

struct inSphViscosityLimiter {
    int bOn;
    int bShockTracker;
    };
void pstSphViscosityLimiter(PST,void *,int,void *,int *);

struct inUpdateShockTracker {
    double dDelta;
    double dShockTrackerA;
    double dShockTrackerB;
};

void pstUpdateShockTracker(PST,void *,int,void *,int *);

#ifdef STARFORM
/* PST_FORMSTARS */
struct inFormStars
{
    struct stfmContext stfm;
    double dTime;
    };

struct outFormStars 
{
    int nFormed;
    int nDeleted;
    double dMassFormed;
    };

void pstFormStars(PST,void *,int,void *,int *);

struct inFeedback
{
    struct fbContext fb;
    struct snContext sn;
    double dTime;
    double dDelta;
    };

struct outFeedback
{
    FBEffects fbTotals[FB_NFEEDBACKS];
    };

void pstFeedback(PST,void *,int,void *,int *);

void pstInitStarLog(PST,void *,int,void *,int *);

struct inFlushStarLog 
{
    char achStarLogFile[PST_FILENAME_SIZE];
    };

void pstFlushStarLog(PST,void *,int,void *,int *);

#endif

#ifdef SIMPLESF
struct inSimpleStarForm
{
    double dRateCoeff;
    double dTMax;
    double dDenMin;
    double dDelta;

	double dTime;
	double dInitStarMass;
	double dESNPerStarMass;
	double dtCoolingShutoff;
    int bdivv;
    };

struct outSimpleStarForm 
{
    int nFormed;
    int nDeleted;
    double dMassFormed;
    };

void pstSimpleStarForm(PST,void *,int,void *,int *);

#endif

#endif

/* Return is pixmap */
void pstDumpFrame(PST,void *,int,void *,int *);

void pstCOM(PST,void *,int,void *,int *);

void pstCOMByType(PST,void *,int,void *,int *);

void pstOldestStar(PST,void *,int,void *,int *);

struct inSetSink 
{
  double dSinkMassMin;
  };

struct outSetSink 
{
  int nSink;
  };

void pstSetSink(PST,void *,int,void *,int *);

struct inFormSinks
{
    int bJeans;
    int bDensity;
    int iKickRung;
    int bSimple;
    double dDensityCut;
    double dTime;
    double dJConst2;
    };

struct outFormSinks 
{
    int nCandidates;
    double Jvalmin;
    };

void pstFormSinks(PST,void *,int,void *,int *);

/* Return is pixmap */
void pstDumpVoxel(PST,void *,int,void *,int *);

/* PST_GROWMASS */
struct inGrowMass 
{
    int nGrowMass;
    int iGrowType;
    double dDeltaM;
    double dMinM;
    double dMaxM;
    };

void pstGrowMass(PST,void *,int,void *,int *);

/* PST_CLEARTIMER */
struct inClearTimer 
{
    int iTimer;
    };

void pstClearTimer(PST,void *,int,void *,int *);

/* PST_MASSINR */
struct inMassInR 
{
    double R;
    };

struct outMassInR
{
    double dMass;
    FLOAT com[3];
    };
void pstMassInR(PST,void *,int,void *,int *);

struct inRotBar
{
    struct rotbarContext rotbar;
    };
void pstInitRotBar(PST,void *,int,void *,int *);

/* PST_KICKVPRED */
#ifdef NEED_VPRED
struct inKickVpred {
	double dvFacOne;
	double dvFacTwo;
	double duDelta;
	double duDotLimit;
	int iGasModel;
	double z;
        double dTimeEnd;
	};
void pstKickVpred(PST,void *,int,void *,int *);
#endif
#endif
