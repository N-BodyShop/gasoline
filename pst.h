#ifndef PST_HINCLUDED
#define PST_HINCLUDED

#include "pkd.h"
#include "mdl.h"

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
	int iSplitDim;
	float fSplit;
	int nStart;
	int nEnd;
	int nOrdSplit;
	int nTotal;
	/*
	 ** The PST node is also a valid cell for the tree.
	 */
	KDN kdn;
	} * PST;


#define PST_SERVICES		50
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
	int iOrder;
	float fExtraStore;
	float fPeriod[3];
	int bParaRead;
	double dvFac;
	char achInFile[PST_FILENAME_SIZE];
	};
void pstReadTipsy(PST,void *,int,void *,int *);

#define PST_DOMAINDECOMP	6
void pstDomainDecomp(PST,void *,int,void *,int *);

#define PST_CALCBOUND		7
struct outCalcBound {
	BND bnd;
	};
void pstCalcBound(PST,void *,int,void *,int *);

#define PST_WEIGHT			8
struct inWeight {
	int iSplitDim;
	float fSplit;
	int iSplitSide;
	int ittr;
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
	int bParaWrite;
	double dvFac;
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
	};
struct outBuildTree {
	KDN kdn;
	};
void pstBuildTree(PST,void *,int,void *,int *);

#define PST_DENSITY			20
struct inDensity {
	int nSmooth;
	int bGatherScatter;
	};
void pstDensity(PST,void *,int,void *,int *);

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
	double dPartSum;
	double dCellSum;
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
	};
void pstKick(PST,void *,int,void *,int *);

#define PST_READCHECK		25
struct inReadCheck {
	int iVersion;
	int iOffset;
	int nStart;
	int nEnd;
	int iOrder;
	float fExtraStore;
	float fPeriod[3];
	int bParaRead;
	char achInFile[PST_FILENAME_SIZE];
	};
void pstReadCheck(PST,void *,int,void *,int *);

#define PST_WRITECHECK		26
struct inWriteCheck {
	int iOffset;
	int nStart;
	int bParaWrite;
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

#define PST_MASSCHECK		35
struct outMassCheck {
	double dMass;
	};
void pstMassCheck(PST,void *,int,void *,int *);

#endif











