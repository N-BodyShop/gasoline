#ifndef PKD_HINCLUDED
#define PKD_HINCLUDED

#include <sys/resource.h>
#include "mdl.h"

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
	int iActive;
	int iRung;
	int cpStart;
	float fWeight;
	float fMass;
	float fSoft;
	float r[3];
	float v[3];
	float a[3];
	float fPot;
	float fBall2;
	float fDensity;
	float dt;			/* a time step suggestion */
#ifdef SUPERCOOL
	float vMean[3];
#endif
#ifdef COLORCODE
	float fColor;
#endif
#ifdef GASOLINE
	float vPred[3];		/* predicted velocity (time centered) */
	float fHsmDivv;		/* 0.5*sqrt(fBall2)*div(vPred) */
	float fRhoDivv;		/* <fDensity*div(vPred)_j> */
	float fCutVisc;
	float u;			/* New thermal energy */ 
	float du;			/* time derivative of thermal energy */
	float uOld;			/* Old thermal energy */ 
	float A;			/* sqrt(uNew_i) prefactor in duNew/dt */
	float B;			/* constant term in duNew/dt */
	float fMetals;
	float fTimeForm;
#endif
#ifdef PLANETS
	float dCollTime;	/* time to next collision */
	int	iCollider;		/* iOrder of collider */
#endif /* PLANETS */
	} PARTICLE;


#define CHECKPOINT_VERSION		4

typedef struct chkParticle {
	int iOrder;
	float fMass;
	float fSoft;
	float r[3];
	float v[3];
	} CHKPART;


typedef struct bndBound {
	float fMin[3];
	float fMax[3];
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
	int pUpper;
	int iLower;
	int iUpper;
	double fMass;
	double fSoft;
	float r[3];
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

typedef struct pkdContext {
	MDL mdl;
	int idSelf;
	int nThreads;
	int nStore;
	int nRejects;
	int nLocal;
	int nActive;
	int nDark;
	int nGas;
	int nStar;
	int nBucket;
	int nLevels;
	int nSplit;
	int nNodes;
	int iOrder;
	int iFreeCell;
	int iRoot;
	float fPeriod[3];
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
void pkdInitialize(PKD *,MDL,int,int,int,float *,int,int,int);
void pkdFinish(PKD);
void pkdReadTipsy(PKD,char *,int,int,int,double,double);
void pkdSetSoft(PKD pkd,double dSoft);
void pkdCalcBound(PKD,BND *,BND *);
int pkdWeight(PKD,int,float,int,int,int,int *,int *,float *,float *);
int pkdLowerPart(PKD,int,float,int,int);
int pkdUpperPart(PKD,int,float,int,int);
int pkdLowerOrdPart(PKD,int,int,int);
int pkdUpperOrdPart(PKD,int,int,int);
void pkdActiveOrder(PKD);
int pkdColRejects(PKD,int,float,float,int);
int pkdSwapRejects(PKD,int);
int pkdSwapSpace(PKD);
int pkdFreeStore(PKD);
int pkdLocal(PKD);
int pkdActive(PKD);
int pkdInactive(PKD);
int pkdNodes(PKD);
void pkdDomainColor(PKD);
int pkdColOrdRejects(PKD,int,int);
void pkdLocalOrder(PKD);
void pkdWriteTipsy(PKD,char *,int,int,double,double);
void pkdCombine(KDN *,KDN *,KDN *);
void pkdCalcCell(PKD,KDN *,float *,int,struct pkdCalcCellStruct *);
double pkdCalcOpen(KDN *,int,double,int);
void pkdBuildLocal(PKD,int,int,double,int,int,KDN *);
void pkdBuildBinary(PKD,int,int,double,int,int,KDN *);
void pkdGravAll(PKD,int,int,int,int,double,double,int *,
				double *,double *,CASTAT *,double *); 
void pkdCalcE(PKD,double *,double *);
void pkdDrift(PKD,double,float *,int);
void pkdKick(PKD pkd,double,double, double, double);
void pkdReadCheck(PKD,char *,int,int,int,int);
void pkdWriteCheck(PKD,char *,int,int);
void pkdDistribCells(PKD,int,KDN *);
void pkdCalcRoot(PKD,struct ilCellNewt *);
void pkdDistribRoot(PKD,struct ilCellNewt *);
void pkdSwapAll(PKD pkd, int idSwap);
double pkdMassCheck(PKD pkd);
void pkdSetRung(PKD pkd, int iRung);
void pkdActiveRung(PKD pkd, int iRung, int bGreater);
int pkdCurrRung(PKD pkd, int iRung);
void pkdDensityStep(PKD pkd, double dEta, double
		    dRhoFac);
void pkdVelocityStep(PKD pkd, double dEta, double dVelFac, double
		     dAccFac);
int pkdDtToRung(PKD pkd, int iRung, double dDelta, int iMaxRung, int
		bAll);
void pkdInitDt(PKD pkd, double dDelta);

int pkdRungParticles(PKD,int);
void pkdCoolVelocity(PKD,int,double,double,double);
void pkdActiveCool(PKD,int);
void pkdInitAccel(PKD);
int pkdOrdWeight(PKD,int,int,int,int,int *,int *);
int pkdIsGas(PKD,PARTICLE *);
int pkdIsDark(PKD,PARTICLE *);
int pkdIsStar(PKD,PARTICLE *);

#ifdef GASOLINE

void pkdActiveGas(PKD);
void pkdCalcEthdot(PKD);
void pkdKickVpred(PKD pdk, double dvFacOne, double dvFacTwo);
int pkdSphCurrRung(PKD pkd, int iRung);

#endif

#endif

