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
	float fWeight;
	float fMass;
	float fSoft;
	float r[3];
	float v[3];
	float a[3];
	float fPot;
#ifdef SMOOTH_CODE
	float fBall2;
	float fDensity;
#endif
	} PARTICLE;

typedef struct chkParticle {
	int iOrder;
	float fMass;
	float fSoft;
	float r[3];
	float v[3];
	} CHKPART;


typedef struct ochkParticle {
	int iOrder;
	float fWeight;
	float fMass;
	float fSoft;
	float r[3];
	float v[3];
	float a[3];
	float fPot;
	float fBall2;
	float fDensity;
	float fColor;
	} OCHKPART;


typedef struct bndBound {
	float fMin[3];
	float fMax[3];
	} BND;

typedef struct kdNode {
	int iDim;
	float fSplit;
	BND bnd;
	int pLower;		/* also doubles as thread id for the LTT */
	int pUpper;
	float fMass;
	float fSoft;
	float r[3];
	float fQxx,fQyy,fQzz,fQxy,fQxz,fQyz;
	float fBmax,fB2,fB3,fB4;
	float fOpen2;
	} KDN;

typedef struct ilPart {
	float m,h;
	float x,y,z;
	} ILP;

typedef struct ilCellSoft {
	float m,h;
	float x,y,z;
	float xx,yy,zz,xy,xz,yz;
	} ILCS;

/*
 ** Reduced Quadrupole moment tensor, Q_ij = 3.0*M_ij - delta_ij*tr(M_ij).
 */
typedef struct ilCellNewt {
	float m;
	float x,y,z;
	float xx,yy,xy,xz,yz;
	} ILCN;

typedef struct ewaldTable {
	float hx,hy,hz;
	float k2,k3x,k3y,k3z;
	} EWT;

typedef struct pkdContext {
	MDL mdl;
	int idSelf;
	int nThreads;
	int nStore;
	int nRejects;
	int nLocal;
	int nBucket;
	int nLevels;
	int nSplit;
	int nNodes;
	int iOrder;
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
	ILP *ilp;
	ILCS *ilcs;
	ILCN *ilcn;
	/*
	 ** Ewald summation setup.
	 */
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
void pkdInitialize(PKD *,MDL,int,int,int,float *);
void pkdFinish(PKD);
void pkdReadTipsy(PKD,char *,int,int);
void pkdSetSoft(PKD pkd,double dSoft);
void pkdCalcBound(PKD,BND *);
int pkdWeight(PKD,int,float,int,int,int,int *,int *,float *,float *);
int pkdLowerPart(PKD,int,float,int,int);
int pkdUpperPart(PKD,int,float,int,int);
int pkdLowerOrdPart(PKD,int,int,int);
int pkdUpperOrdPart(PKD,int,int,int);
int pkdColRejects(PKD,int,float,int);
int pkdSwapRejects(PKD,int);
int pkdSwapSpace(PKD);
int pkdFreeStore(PKD);
int pkdLocal(PKD);
int pkdNodes(PKD);
void pkdDomainColor(PKD);
int pkdColOrdRejects(PKD,int,int);
void pkdLocalOrder(PKD,int);
void pkdWriteTipsy(PKD,char *,int,int,int);
void pkdBuildLocal(PKD,int,int,double,KDN *);
void pkdBuildTop(PKD,int,double,KDN *,int,int *);
void pkdGravAll(PKD,int,int,float,float,double *,double *,CASTAT *,CASTAT *); 
void pkdCalcE(PKD,double *,double *);
void pkdDrift(PKD,double,float *);
void pkdKick(PKD pkd,double,double);
void pkdReadCheckOld(PKD,char *,int,int);
void pkdReadCheckNew(PKD,char *,int,int);
void pkdWriteCheckNew(PKD,char *,int);

#endif

