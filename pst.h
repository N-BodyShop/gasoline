#ifndef PST_HINCLUDED
#define PST_HINCLUDED

#include "pkd.h"
#include "mdl.h"
#include "smoothfcn.h"

typedef struct lclBlock {
	char *pszDataPath;
	PKD	pkd;
	int nPstLvl;
	int iWtFrom;
	int iWtTo;
	int iPart;
	float fSplit;
	float fWtLow;
	float fWtHigh;
	float fLow;
	float fHigh;
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
	BND	bnd;
	BND bndActive;
	int iSplitDim;
	float fSplit;
	float fSplitInactive;
	int nStart;
	int nEnd;
	int nOrdSplit;
	int nTotal;
	/*
	 ** The PST node is also a valid cell for the tree.
	 */
	KDN kdn;
	} * PST;


#define PST_SERVICES		51
#define PST_FILENAME_SIZE	512

void pstAddServices(PST,MDL);
void pstInitialize(PST *,MDL,LCL *);
void pstFinish(PST);

#define PST_SETADD			2
struct inSetAdd {
	int id;
	};
void pstSetAdd(PST,void *,int,void *,int *);

#define PST_LEVELIZE		3
struct inLevelize {
	int iLvl;
	};
void pstLevelize(PST,void *,int,void *,int *);

#define PST_READTIPSY		5
struct inReadTipsy {
	int nStart;
	int nEnd;
	int nDark;	
	int nGas;
	int nStar;
	int iOrder;
	float fExtraStore;
	float fPeriod[3];
	int bStandard;
	double dvFac;
	double dTuFac;
	char achInFile[PST_FILENAME_SIZE];
	};
void pstReadTipsy(PST,void *,int,void *,int *);

#define PST_DOMAINDECOMP	6
void pstDomainDecomp(PST,void *,int,void *,int *);

#define PST_CALCBOUND		7
struct outCalcBound {
	BND bnd;
	BND bndActive;
	};
void pstCalcBound(PST,void *,int,void *,int *);

#define PST_WEIGHT			8
struct inWeight {
	int iSplitDim;
	float fSplit;
	int iSplitSide;
	int ittr;
	int pFlag;
	};
struct outWeight {
	int nLow;
	int nHigh;
	float fLow;
	float fHigh;
	};
void pstWeight(PST,void *,int,void *,int *);

#define PST_FREESTORE		9
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

#define PST_COLREJECTS		10
struct inColRejects {
	int iSplitDim;
	float fSplit;
	float fSplitInactive;
	int iSplitSide;
	};
void pstColRejects(PST,void *,int,void *,int *);

#define PST_SWAPREJECTS		11
void pstSwapRejects(PST,void *,int,void *,int *);

#define PST_DOMAINCOLOR		12
void pstDomainColor(PST,void *,int,void *,int *);

#define PST_COLORDREJECTS	13
struct inColOrdRejects {
	int nOrdSplit;
	int iSplitSide;
	};
void pstColOrdRejects(PST,void *,int,void *,int *);

#define PST_DOMAINORDER		14
void pstDomainOrder(PST,void *,int,void *,int *);

#define PST_LOCALORDER		15
void pstLocalOrder(PST,void *,int,void *,int *);

#define PST_OUTARRAY		16
struct inOutArray {
	char achOutFile[PST_FILENAME_SIZE];
	int iType;
	};
void pstOutArray(PST,void *,int,void *,int *);

#define PST_OUTVECTOR		17
struct inOutVector {
	char achOutFile[PST_FILENAME_SIZE];
	int iDim;
	int iType;
	};
void pstOutVector(PST,void *,int,void *,int *);

#define PST_WRITETIPSY		18
struct inWriteTipsy {
	int bStandard;
	double dvFac;
	double duTFac;
	char achOutFile[PST_FILENAME_SIZE];
	};
void pstWriteTipsy(PST,void *,int,void *,int *);

#define PST_BUILDTREE		19
struct inBuildTree {
	int nBucket;
	int iOpenType;
	int iOrder;
	double dCrit;
	int bBinary;
	int bActiveOnly;
	};
struct outBuildTree {
	KDN kdn;
	};
void pstBuildTree(PST,void *,int,void *,int *);

#define PST_SMOOTH			20
struct inSmooth {
	int nSmooth;
    int bPeriodic;
	int bSymmetric;
	int iSmoothType;
	SMF smf;
	};
void pstSmooth(PST,void *,int,void *,int *);

#define PST_GRAVITY			21
struct inGravity {
	int nReps;
	int bPeriodic;
	int iOrder;
	int iEwOrder;
	double dEwCut;
	double dEwhCut;
	};
struct outGravity {
	int nActive;
	double dPartSum;
	double dCellSum;
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

#define PST_CALCE			22
struct outCalcE {
	double T;
	double U;
	};
void pstCalcE(PST,void *,int,void *,int *);

#define PST_DRIFT			23
struct inDrift {
	double dDelta;
	float fCenter[3];
	int bPeriodic;
	};
void pstDrift(PST,void *,int,void *,int *);

#define PST_KICK			24
struct inKick {
	double dvFacOne;
	double dvFacTwo;
	double dvPredFacOne;
	double dvPredFacTwo;
	};
void pstKick(PST,void *,int,void *,int *);

#define PST_READCHECK		25
struct inReadCheck {
	int iVersion;
	int iOffset;
	int nStart;
	int nEnd;
	int nDark;
	int nGas;
	int nStar;
	int iOrder;
	float fExtraStore;
	float fPeriod[3];
	char achInFile[PST_FILENAME_SIZE];
	};
void pstReadCheck(PST,void *,int,void *,int *);

#define PST_WRITECHECK		26
struct inWriteCheck {
	int iOffset;
	int nStart;
	char achOutFile[PST_FILENAME_SIZE];
	};
void pstWriteCheck(PST,void *,int,void *,int *);

#define PST_SETSOFT		27
struct inSetSoft {
	double dSoft;
	};
void pstSetSoft(PST,void *,int,void *,int *);

#define PST_SETTOTAL		28
struct outSetTotal {
	int nTotal;
	};
void pstSetTotal(PST,void *,int,void *,int *);

#define PST_CALCCELL		30
struct inCalcCell {
	int iOrder;
	float rcm[3];
	};
struct outCalcCell {
	struct pkdCalcCellStruct mom;
	};
void pstCalcCell(PST,void *,int,void *,int *);

#define PST_COLCELLS		31
struct inColCells {
	int iCell;
	int nCell;
	};
void pstColCells(PST,void *,int,void *,int *);

#define PST_DISTRIBCELLS	32
void pstDistribCells(PST,void *,int,void *,int *);

#define PST_CALCROOT		33
struct ioCalcRoot {
	struct ilCellNewt ilcn;
	};
void pstCalcRoot(PST,void *,int,void *,int *);

#define PST_DISTRIBROOT		34
void pstDistribRoot(PST,void *,int,void *,int *);

#define PST_ONENODEREADINIT	35
void pstOneNodeReadInit(PST pst,void *vin,int nIn,void *vout,int *pnOut);

#define PST_SWAPALL		36
void pstSwapAll(PST pst,void *vin,int nIn,void *vout,int *pnOut);

#define PST_MASSCHECK		37
struct outMassCheck {
	double dMass;
	};
void pstMassCheck(PST,void *,int,void *,int *);

#define PST_ACTIVEORDER		38
void pstActiveOrder(PST,void *,int,void *,int *);

#define PST_SETRUNG		39
struct inSetRung {
    int iRung;
    };
void pstSetRung(PST,void *,int,void *,int *);

#define PST_ACTIVERUNG		40
struct inActiveRung {
    int iRung;
    int bGreater;
    };
void pstActiveRung(PST,void *,int,void *,int *);

#define PST_CURRRUNG		41
struct inCurrRung {
    int iRung;
    };
struct outCurrRung {
    int iCurrent;
    };
void pstCurrRung(PST,void *,int,void *,int *);

#define PST_DENSITYRUNG		42
struct inDensityRung {
    int iRung;
    double dDelta;
    double dEta;
    double dRhoFac;
    int bAll;
    };
struct outDensityRung {
    int iMaxRung;
    };
void pstDensityRung(PST,void *,int,void *,int *);

#define PST_RUNGSTATS		43
struct inRungStats {
	int iRung;
	};
struct outRungStats {
	int nParticles;
	};
void pstRungStats(PST,void *,int,void *,int *);

#define PST_GETMAP			44
struct inGetMap {
	int nStart;
	};
void pstGetMap(PST,void *,int,void *,int *);

#define PST_VELOCITYRUNG		45
struct inVelocityRung {
    int iRung;
    double dDelta;
    double dEta;
    int iMaxRung;
    double dVelFac;
    double dAccFac;
    int bAll;
    };
struct outVelocityRung {
    int iMaxRung;
    };
void pstVelocityRung(PST,void *,int,void *,int *);

#define PST_COOLVELOCITY		46
struct inCoolVelocity {
	int nSuperCool;
	double dCoolFac;
	double dCoolDens;
	double dCoolMaxDens;
	};
void pstCoolVelocity(PST,void *,int,void *,int *);

#define PST_ACTIVECOOL			47
struct inActiveCool {
	int nSuperCool;
	};
void pstActiveCool(PST,void *,int,void *,int *);

#define PST_RESMOOTH			48
struct inReSmooth {
	int nSmooth;
    int bPeriodic;
	int bSymmetric;
	int iSmoothType;
	SMF smf;
	};
void pstReSmooth(PST,void *,int,void *,int *);

#ifdef GASOLINE

#define PST_ACTIVEGAS			49
void pstActiveGas(PST,void *,int,void *,int *);

#define PST_CALCETHDOT			50
void pstCalcEthdot(PST,void *,int,void *,int *);

#define PST_KICKVPRED			51
struct inKickVpred {
	double dvFacOne;
	double dvFacTwo;
	};
void pstKickVpred(PST,void *,int,void *,int *);

#define PST_SPHCURRRUNG		52
struct inSphCurrRung {
    int iRung;
    };
struct outSphCurrRung {
    int iCurrent;
    };
void pstSphCurrRung(PST,void *,int,void *,int *);

#endif

#endif

