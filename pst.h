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
	int nTopNodes;
	KDN *kdTop;
	int *piLeaf;
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
	} * PST;


#define PST_SERVICES		50
#define PST_FILENAME_SIZE	512

void pstInitialize(PST *,MDL,LCL *);
void pstFinish(PST);

#define PST_SETADD			2
struct inSetAdd {
	int id;
	};
void pstSetAdd(PST,char *,int,char *,int *);

#define PST_LEVELIZE		3
struct inLevelize {
	int iLvl;
	};
void pstLevelize(PST,char *,int,char *,int *);

#define PST_SHOWPST			4
void pstShowPst(PST,char *,int,char *,int *);

#define PST_READTIPSY		5
struct inReadTipsy {
	int nStart;
	int nEnd;
	int iOrder;
	float fExtraStore;
	float fPeriod[3];
	char achInFile[PST_FILENAME_SIZE];
	};
void pstReadTipsy(PST,char *,int,char *,int *);

#define PST_DOMAINDECOMP	6
void pstDomainDecomp(PST,char *,int,char *,int *);

#define PST_CALCBOUND		7
struct outCalcBound {
	BND bnd;
	};
void pstCalcBound(PST,char *,int,char *,int *);

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
void pstWeight(PST,char *,int,char *,int *);

#define PST_FREESTORE		9
struct outFreeStore {
	int nFreeStore;
	};
void pstFreeStore(PST,char *,int,char *,int *);

/*
 ** This structure is used by reject collectors and SwapRejects
 */
typedef struct outReject {
	int id;
	int nRejects;
	int nSpace;
	} OREJ;

#define PST_COLREJECTS		10
struct inColRejects {
	int iSplitDim;
	float fSplit;
	int iSplitSide;
	};
void pstColRejects(PST,char *,int,char *,int *);

#define PST_SWAPREJECTS		11
void pstSwapRejects(PST,char *,int,char *,int *);

#define PST_DOMAINCOLOR		12
void pstDomainColor(PST,char *,int,char *,int *);

#define PST_COLORDREJECTS	13
struct inColOrdRejects {
	int nOrdSplit;
	int iSplitSide;
	};
void pstColOrdRejects(PST,char *,int,char *,int *);

#define PST_DOMAINORDER		14
void pstDomainOrder(PST,char *,int,char *,int *);

#define PST_LOCALORDER		15
void pstLocalOrder(PST,char *,int,char *,int *);

#define PST_OUTARRAY		16
struct inOutArray {
	char achOutFile[PST_FILENAME_SIZE];
	int iType;
	};
void pstOutArray(PST,char *,int,char *,int *);

#define PST_OUTVECTOR		17
struct inOutVector {
	char achOutFile[PST_FILENAME_SIZE];
	int iDim;
	int iType;
	};
void pstOutVector(PST,char *,int,char *,int *);

#define PST_WRITETIPSY		18
struct inWriteTipsy {
	char achOutFile[PST_FILENAME_SIZE];
	};
void pstWriteTipsy(PST,char *,int,char *,int *);

#define PST_BUILDTREE		19
struct inBuildTree {
	int iCell;
	int nBucket;
	int iOpenType;
	double dCrit;
	};
void pstBuildTree(PST,char *,int,char *,int *);

#define PST_DENSITY			20
struct inDensity {
	int nSmooth;
	int bGatherScatter;
	};
void pstDensity(PST,char *,int,char *,int *);

#define PST_GRAVITY			21
struct inGravity {
	int nReps;
	int bPeriodic;
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
	};
void pstGravity(PST,char *,int,char *,int *);

#define PST_CALCE			22
struct outCalcE {
	double T;
	double U;
	};
void pstCalcE(PST,char *,int,char *,int *);

#define PST_DRIFT			23
struct inDrift {
	double dDelta;
	float fCenter[3];
	};
void pstDrift(PST,char *,int,char *,int *);

#define PST_KICK			24
struct inKick {
	double dvFacOne;
	double dvFacTwo;
	};
void pstKick(PST,char *,int,char *,int *);

#define PST_READCHECK		25
struct inReadCheck {
	int bNewCheck;
	int nStart;
	int nEnd;
	int iOrder;
	float fExtraStore;
	float fPeriod[3];
	char achInFile[PST_FILENAME_SIZE];
	};
void pstReadCheck(PST,char *,int,char *,int *);

#define PST_WRITECHECK		26
struct inWriteCheck {
	int bNewCheck;
	int nStart;
	char achOutFile[PST_FILENAME_SIZE];
	};
void pstWriteCheck(PST,char *,int,char *,int *);

#define PST_SETSOFT		27
struct inSetSoft {
	double dSoft;
	};
void pstSetSoft(PST,char *,int,char *,int *);

#define PST_SETTOTAL		28
struct outSetTotal {
	int nTotal;
	};
void pstSetTotal(PST,char *,int,char *,int *);

#endif











