#ifndef SMOOTH_HINCLUDED
#define SMOOTH_HINCLUDED

#include "pkd.h"
#include "smoothfcn.h"
#include "floattype.h"

#define PQ_LOAD_FACTOR 0.50


typedef struct pqNode {
	FLOAT fKey;
	struct pqNode *pqLoser;
	struct pqNode *pqFromInt;
	struct pqNode *pqFromExt;
	struct pqNode *pqWinner;	/* Only used when building initial tree */
	struct pqNode *pNext;
	struct pqNode *pPrev;
	int p;
	int id;
	PARTICLE *pPart;
	FLOAT dx;
	FLOAT dy;
	FLOAT dz;
	FLOAT ax;
	FLOAT ay;
	FLOAT az;
	} PQ;


typedef struct smContext {
	PKD pkd;
	int nSmooth;
	int bPeriodic;
	void (*fcnSmooth)(PARTICLE *,int,NN *,SMF *);
	void (*fcnPost)(PARTICLE *,SMF *);
	int *piMark;
	int nListSize;
	NN *nnList;
	int *pbRelease;
	PQ *pq;
	PQ *pqHead;
	int nHash;
	PQ **pqHash;
	} * SMX;


#define PQ_INIT(pq,n)\
{\
	int j;\
	if ((n) == 1) {\
		(pq)[0].pqFromInt = NULL;\
		(pq)[0].pqFromExt = NULL;\
		}\
	else {\
	    for (j=0;j<(n);++j) {\
		    if (j < 2) (pq)[j].pqFromInt = NULL;\
		    else (pq)[j].pqFromInt = &(pq)[j>>1];\
		    (pq)[j].pqFromExt = &(pq)[(j+(n))>>1];\
		    }\
	    }\
    }


#define PQ_BUILD(pq,n,q)\
{\
    int i,j;\
	PQ *t,*lt;\
	for (j=(n)-1;j>0;--j) {\
		i = (j<<1);\
		if (i < (n)) t = (pq)[i].pqWinner;\
		else t = &(pq)[i-(n)];\
		++i;\
		if (i < (n)) lt = (pq)[i].pqWinner;\
		else lt = &(pq)[i-(n)];\
		if (t->fKey < lt->fKey) {\
			(pq)[j].pqLoser = t;\
			(pq)[j].pqWinner = lt;\
			}\
		else {\
			(pq)[j].pqLoser = lt;\
			(pq)[j].pqWinner = t;\
			}\
		}\
    if ((n) == 1) (q) = (pq);\
	else (q) = (pq)[1].pqWinner;\
	}


#define PQ_REPLACE(q)\
{\
	PQ *t,*lt;\
	t = (q)->pqFromExt;\
	while (t) {\
		if (t->pqLoser->fKey > (q)->fKey) {\
			lt = t->pqLoser;\
			t->pqLoser = (q);\
			(q) = lt;\
			}\
		t = t->pqFromInt;\
		}\
	}


#define PQ_HASHADD(pqHash,nHash,pq)\
{\
	int i;\
	PQ *t;\
	i = (pq)->p%(nHash);\
	t = (pqHash)[i];\
/*\
	assert((pq) != t);\
*/\
	if (t) t->pPrev = (pq);\
	(pq)->pPrev = NULL;\
	(pq)->pNext = t;\
	(pqHash)[i] = (pq);\
	}

#define PQ_HASHDEL(pqHash,nHash,pq)\
{\
	int i;\
	PQ *t,*lt;\
	t = (pq)->pNext;\
	lt = (pq)->pPrev;\
	if (t) t->pPrev = lt;\
	if (lt) lt->pNext = t;\
	else {\
		i = (pq)->p%(nHash);\
/*\
		assert((pqHash)[i] == (pq));\
*/\
		(pqHash)[i] = t;\
		}\
	}

#define PQ_INQUEUE(pqHash,nHash,pRef,idRef,goto_label)\
{\
	PQ *t;\
	t = (pqHash)[(pRef)%(nHash)];\
	while (t) {\
		if (t->p == (pRef) && t->id == (idRef)) goto goto_label;\
		t = t->pNext;\
		}\
	}

#define INTERSECTNP(pkdn,fBall2,x,y,z,label)\
{\
	FLOAT INTRSCT_dx,INTRSCT_dy,INTRSCT_dz;\
	FLOAT INTRSCT_dx1,INTRSCT_dy1,INTRSCT_dz1;\
    FLOAT INTRSCT_fDist2;\
	INTRSCT_dx = (pkdn)->bnd.fMin[0]-x;\
	INTRSCT_dx1 = x-(pkdn)->bnd.fMax[0];\
	INTRSCT_dy = (pkdn)->bnd.fMin[1]-y;\
	INTRSCT_dy1 = y-(pkdn)->bnd.fMax[1];\
	INTRSCT_dz = (pkdn)->bnd.fMin[2]-z;\
	INTRSCT_dz1 = z-(pkdn)->bnd.fMax[2];\
	if (INTRSCT_dx > 0.0) INTRSCT_fDist2 = INTRSCT_dx*INTRSCT_dx;\
	else if (INTRSCT_dx1 > 0.0) INTRSCT_fDist2 = INTRSCT_dx1*INTRSCT_dx1;\
	else INTRSCT_fDist2 = 0.0;\
	if (INTRSCT_dy > 0.0) INTRSCT_fDist2 += INTRSCT_dy*INTRSCT_dy;\
	else if (INTRSCT_dy1 > 0.0) INTRSCT_fDist2 += INTRSCT_dy1*INTRSCT_dy1;\
	if (INTRSCT_dz > 0.0) INTRSCT_fDist2 += INTRSCT_dz*INTRSCT_dz;\
	else if (INTRSCT_dz1 > 0.0) INTRSCT_fDist2 += INTRSCT_dz1*INTRSCT_dz1;\
    if (INTRSCT_fDist2 > fBall2) goto label;\
	}


#define INTERSECT(pkdn,fBall2,lx,ly,lz,x,y,z,sx,sy,sz,bPeriodic,label)\
{\
	FLOAT INTRSCT_dx,INTRSCT_dy,INTRSCT_dz;\
	FLOAT INTRSCT_dx1,INTRSCT_dy1,INTRSCT_dz1;\
    FLOAT INTRSCT_fDist2;\
	INTRSCT_dx = (pkdn)->bnd.fMin[0]-x;\
	INTRSCT_dx1 = x-(pkdn)->bnd.fMax[0];\
	INTRSCT_dy = (pkdn)->bnd.fMin[1]-y;\
	INTRSCT_dy1 = y-(pkdn)->bnd.fMax[1];\
	INTRSCT_dz = (pkdn)->bnd.fMin[2]-z;\
	INTRSCT_dz1 = z-(pkdn)->bnd.fMax[2];\
	if (INTRSCT_dx > 0.0) {\
		INTRSCT_dx1 += lx;\
		if (INTRSCT_dx1 < INTRSCT_dx) {\
			INTRSCT_fDist2 = INTRSCT_dx1*INTRSCT_dx1;\
			sx = x+lx;\
			bPeriodic = 1;\
			}\
		else {\
			INTRSCT_fDist2 = INTRSCT_dx*INTRSCT_dx;\
			sx = x;\
			}\
		if (INTRSCT_fDist2 > fBall2) goto label;\
		}\
	else if (INTRSCT_dx1 > 0.0) {\
		INTRSCT_dx += lx;\
		if (INTRSCT_dx < INTRSCT_dx1) {\
			INTRSCT_fDist2 = INTRSCT_dx*INTRSCT_dx;\
			sx = x-lx;\
			bPeriodic = 1;\
			}\
		else {\
			INTRSCT_fDist2 = INTRSCT_dx1*INTRSCT_dx1;\
			sx = x;\
			}\
		if (INTRSCT_fDist2 > fBall2) goto label;\
		}\
	else {\
		INTRSCT_fDist2 = 0.0;\
		sx = x;\
		}\
	if (INTRSCT_dy > 0.0) {\
		INTRSCT_dy1 += ly;\
		if (INTRSCT_dy1 < INTRSCT_dy) {\
			INTRSCT_fDist2 += INTRSCT_dy1*INTRSCT_dy1;\
			sy = y+ly;\
			bPeriodic = 1;\
			}\
		else {\
			INTRSCT_fDist2 += INTRSCT_dy*INTRSCT_dy;\
			sy = y;\
			}\
		if (INTRSCT_fDist2 > fBall2) goto label;\
		}\
	else if (INTRSCT_dy1 > 0.0) {\
		INTRSCT_dy += ly;\
		if (INTRSCT_dy < INTRSCT_dy1) {\
			INTRSCT_fDist2 += INTRSCT_dy*INTRSCT_dy;\
			sy = y-ly;\
			bPeriodic = 1;\
			}\
		else {\
			INTRSCT_fDist2 += INTRSCT_dy1*INTRSCT_dy1;\
			sy = y;\
			}\
		if (INTRSCT_fDist2 > fBall2) goto label;\
		}\
	else {\
		sy = y;\
		}\
	if (INTRSCT_dz > 0.0) {\
		INTRSCT_dz1 += lz;\
		if (INTRSCT_dz1 < INTRSCT_dz) {\
			INTRSCT_fDist2 += INTRSCT_dz1*INTRSCT_dz1;\
			sz = z+lz;\
			bPeriodic = 1;\
			}\
		else {\
			INTRSCT_fDist2 += INTRSCT_dz*INTRSCT_dz;\
			sz = z;\
			}\
		if (INTRSCT_fDist2 > fBall2) goto label;\
		}\
	else if (INTRSCT_dz1 > 0.0) {\
		INTRSCT_dz += lz;\
		if (INTRSCT_dz < INTRSCT_dz1) {\
			INTRSCT_fDist2 += INTRSCT_dz*INTRSCT_dz;\
			sz = z-lz;\
			bPeriodic = 1;\
			}\
		else {\
			INTRSCT_fDist2 += INTRSCT_dz1*INTRSCT_dz1;\
			sz = z;\
			}\
		if (INTRSCT_fDist2 > fBall2) goto label;\
		}\
	else {\
		sz = z;\
		}\
	}


#define INTERSECTSCATTERNP(pkdn,x,y,z,label)\
{\
        if ( x<(pkdn)->bndBall.fMin[0] || x>(pkdn)->bndBall.fMax[0] ||\
             y<(pkdn)->bndBall.fMin[1] || y>(pkdn)->bndBall.fMax[1] ||\
             z<(pkdn)->bndBall.fMin[2] || z>(pkdn)->bndBall.fMax[2] ) goto label; \
	}


#define INTERSECTSCATTER(pkdn,lx,ly,lz,x,y,z,sx,sy,sz,bPeriodic,label)\
{\
        if (x < (pkdn)->bndBall.fMin[0]) { \
                sx=x+lx; \
                if (sx > (pkdn)->bndBall.fMax[0]) goto label; \
                bPeriodic = 1; \
                } \
        else if (x > (pkdn)->bndBall.fMax[0]) { \
                sx=x-lx; \
                if (sx < (pkdn)->bndBall.fMin[0]) goto label; \
                bPeriodic = 1; \
                } \
        else { \
                sx=x; \
	        } \
        if (y < (pkdn)->bndBall.fMin[1]) { \
                sy=y+ly; \
                if (sy > (pkdn)->bndBall.fMax[1]) goto label; \
                bPeriodic = 1; \
                } \
        else if (y > (pkdn)->bndBall.fMax[1]) { \
                sy=y-ly; \
                if (sy < (pkdn)->bndBall.fMin[1]) goto label; \
                bPeriodic = 1; \
                } \
        else { \
                sy=y; \
	        } \
        if (z < (pkdn)->bndBall.fMin[2]) { \
                sz=z+lz; \
                if (sz > (pkdn)->bndBall.fMax[2]) goto label; \
                bPeriodic = 1; \
                } \
        else if (z > (pkdn)->bndBall.fMax[2]) { \
                sz=z-lz; \
                if (sz < (pkdn)->bndBall.fMin[2]) goto label; \
                bPeriodic = 1; \
                } \
        else { \
                sz=z; \
	        } \
	}


int smInitialize(SMX *,PKD,SMF *,int,int,int,int,int);
void smFinish(SMX,SMF *);
void smSmooth(SMX,SMF *);
void smMarkSmooth(SMX,SMF *, int);
void smReSmooth(SMX,SMF *);
void smGrowList(SMX smx);

#ifdef COLLISIONS
void smQQSmooth(SMX smx, SMF *smf);
#endif /* COLLISIONS */

#endif











