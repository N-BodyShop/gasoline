#include "define.h"
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

#define LOWHFIX_NONE         0
#define LOWHFIX_HOVERSOFT    1 
#define LOWHFIX_SINKRADIUS   2
#define LOWHFIX_SINKRADIUS_BUFF   3

typedef struct smContext {
	PKD pkd;
	int nSmooth;
    int nSmoothMax;
    int nSmoothInner;
    double fBall2InnerFrac;
	int bPeriodic;
	int iLowhFix;
    int bUseBallMax;
	double dfBall2OverSoft2;
    FLOAT lx,ly,lz;
	void (*fcnSmooth)(PARTICLE *,int,NN *,SMF *);
	void (*fcnPost)(PARTICLE *,SMF *);
	int *piMark; /* deprecated MARK code */
	int nListSize;
	NN *nnList;
	int *pbRelease;
	PQ *pq;
	PQ *pqHead;
	int nHash;
	PQ **pqHash;
#ifdef SLIDING_PATCH
    double dTime;
    PATCH_PARAMS *PP;
#endif
	} * SMX;


#define PQ_INIT(pq,n)                               \
    {                                               \
	int j;                                          \
	if ((n) == 1) {                                 \
		(pq)[0].pqFromInt = NULL;                   \
		(pq)[0].pqFromExt = NULL;                   \
		}                                           \
	else {                                          \
	    for (j=0;j<(n);++j) {                       \
		    if (j < 2) (pq)[j].pqFromInt = NULL;    \
		    else (pq)[j].pqFromInt = &(pq)[j>>1];   \
		    (pq)[j].pqFromExt = &(pq)[(j+(n))>>1];  \
		    }                                       \
	    }                                           \
    }


#define PQ_BUILD(pq,n,q)                        \
    {                                           \
    int i,j;                                    \
	PQ *t,*lt;                                  \
	for (j=(n)-1;j>0;--j) {                     \
		i = (j<<1);                             \
		if (i < (n)) t = (pq)[i].pqWinner;      \
		else t = &(pq)[i-(n)];                  \
		++i;                                    \
		if (i < (n)) lt = (pq)[i].pqWinner;     \
		else lt = &(pq)[i-(n)];                 \
		if (t->fKey < lt->fKey) {               \
			(pq)[j].pqLoser = t;                \
			(pq)[j].pqWinner = lt;              \
			}                                   \
		else {                                  \
			(pq)[j].pqLoser = lt;               \
			(pq)[j].pqWinner = t;               \
			}                                   \
		}                                       \
    if ((n) == 1) (q) = (pq);                   \
	else (q) = (pq)[1].pqWinner;                \
	}


#define PQ_REPLACE(q)                           \
    {                                           \
	PQ *t,*lt;                                  \
	t = (q)->pqFromExt;                         \
	while (t) {                                 \
		if (t->pqLoser->fKey > (q)->fKey) {     \
			lt = t->pqLoser;                    \
			t->pqLoser = (q);                   \
			(q) = lt;                           \
			}                                   \
		t = t->pqFromInt;                       \
		}                                       \
	}


#define PQ_HASHADD(pqHash,nHash,pq)             \
    {                                           \
	int i;                                      \
	PQ *t;                                      \
	i = (pq)->p%(nHash);                        \
	t = (pqHash)[i];                            \
/*                                              \
  assert((pq) != t);                            \
*/                                              \
	if (t) t->pPrev = (pq);                     \
	(pq)->pPrev = NULL;                         \
	(pq)->pNext = t;                            \
	(pqHash)[i] = (pq);                         \
	}

#define PQ_HASHDEL(pqHash,nHash,pq)             \
    {                                           \
	int i;                                      \
	PQ *t,*lt;                                  \
	t = (pq)->pNext;                            \
	lt = (pq)->pPrev;                           \
	if (t) t->pPrev = lt;                       \
	if (lt) lt->pNext = t;                      \
	else {                                      \
		i = (pq)->p%(nHash);                    \
/*                                              \
  assert((pqHash)[i] == (pq));                  \
*/                                              \
		(pqHash)[i] = t;                        \
		}                                       \
	}

#define PQ_INQUEUE(pqHash,nHash,pRef,idRef,goto_label)              \
    {                                                               \
	PQ *t;                                                          \
	t = (pqHash)[(pRef)%(nHash)];                                   \
	while (t) {                                                     \
		if (t->p == (pRef) && t->id == (idRef)) goto goto_label;    \
		t = t->pNext;                                               \
		}                                                           \
	}

#define INTERSECTNP(pkdn,fBall2,x,y,z,label)                            \
    {                                                                   \
	FLOAT INTRSCT_dx,INTRSCT_dy,INTRSCT_dz;                             \
	FLOAT INTRSCT_dx1,INTRSCT_dy1,INTRSCT_dz1;                          \
    FLOAT INTRSCT_fDist2;                                               \
	INTRSCT_dx = (pkdn)->bnd.fMin[0]-x;                                 \
	INTRSCT_dx1 = x-(pkdn)->bnd.fMax[0];                                \
	INTRSCT_dy = (pkdn)->bnd.fMin[1]-y;                                 \
	INTRSCT_dy1 = y-(pkdn)->bnd.fMax[1];                                \
	INTRSCT_dz = (pkdn)->bnd.fMin[2]-z;                                 \
	INTRSCT_dz1 = z-(pkdn)->bnd.fMax[2];                                \
	if (INTRSCT_dx > 0.0) INTRSCT_fDist2 = INTRSCT_dx*INTRSCT_dx;       \
	else if (INTRSCT_dx1 > 0.0) INTRSCT_fDist2 = INTRSCT_dx1*INTRSCT_dx1; \
	else INTRSCT_fDist2 = 0.0;                                          \
	if (INTRSCT_dy > 0.0) INTRSCT_fDist2 += INTRSCT_dy*INTRSCT_dy;      \
	else if (INTRSCT_dy1 > 0.0) INTRSCT_fDist2 += INTRSCT_dy1*INTRSCT_dy1; \
	if (INTRSCT_dz > 0.0) INTRSCT_fDist2 += INTRSCT_dz*INTRSCT_dz;      \
	else if (INTRSCT_dz1 > 0.0) INTRSCT_fDist2 += INTRSCT_dz1*INTRSCT_dz1; \
    if (INTRSCT_fDist2 > fBall2) goto label;                            \
	}

#ifdef SLIDING_PATCH
#define INTERSECT(pkdn,fBall2,lx,ly,lz,x,y,z,sx,sy,sz,bPeriodic,label)  \
    {                                                                   \
	FLOAT INTRSCT_dx,INTRSCT_dy,INTRSCT_dz;                             \
	FLOAT INTRSCT_dx1,INTRSCT_dy1,INTRSCT_dz1;                          \
	FLOAT INTRSCT_fDist2;                                               \
	sy = y;                                                             \
	INTRSCT_dx = (pkdn)->bnd.fMin[0]-x;                                 \
	INTRSCT_dx1 = x-(pkdn)->bnd.fMax[0];                                \
	if (INTRSCT_dx > 0.0) {                                             \
		INTRSCT_dx1 += lx;                                              \
		if (INTRSCT_dx1 < INTRSCT_dx) {                                 \
			INTRSCT_fDist2 = INTRSCT_dx1*INTRSCT_dx1;                   \
			sx = x+lx;                                                  \
			sy += SHEAR(1,smx->dTime,smx->PP);                          \
			if (sy >= 0.5*ly) sy -= ly;                                 \
			else if (sy < - 0.5*ly) sy += ly;                           \
			bPeriodic = 1;                                              \
			}                                                           \
		else {                                                          \
			INTRSCT_fDist2 = INTRSCT_dx*INTRSCT_dx;                     \
			sx = x;                                                     \
			}                                                           \
		if (INTRSCT_fDist2 > fBall2) goto label;                        \
		}                                                               \
	else if (INTRSCT_dx1 > 0.0) {                                       \
		INTRSCT_dx += lx;                                               \
		if (INTRSCT_dx < INTRSCT_dx1) {                                 \
			INTRSCT_fDist2 = INTRSCT_dx*INTRSCT_dx;                     \
			sx = x-lx;                                                  \
			sy += SHEAR(-1,smx->dTime,smx->PP);                         \
			if (sy >= 0.5*ly) sy -= ly;                                 \
			else if (sy < - 0.5*ly) sy += ly;                           \
			bPeriodic = 1;                                              \
			}                                                           \
		else {                                                          \
			INTRSCT_fDist2 = INTRSCT_dx1*INTRSCT_dx1;                   \
			sx = x;                                                     \
			}                                                           \
		if (INTRSCT_fDist2 > fBall2) goto label;                        \
		}                                                               \
	else {                                                              \
		INTRSCT_fDist2 = 0.0;                                           \
		sx = x;                                                         \
		}                                                               \
	INTRSCT_dy = (pkdn)->bnd.fMin[1]-sy;                                \
	INTRSCT_dy1 = sy-(pkdn)->bnd.fMax[1];                               \
	if (INTRSCT_dy > 0.0) {                                             \
		INTRSCT_dy1 += ly;                                              \
		if (INTRSCT_dy1 < INTRSCT_dy) {                                 \
			INTRSCT_fDist2 += INTRSCT_dy1*INTRSCT_dy1;                  \
			sy += ly;                                                   \
			bPeriodic = 1;                                              \
			}                                                           \
		else {                                                          \
			INTRSCT_fDist2 += INTRSCT_dy*INTRSCT_dy;                    \
			}                                                           \
		if ((pkdn)->bnd.fMax[0] - (pkdn)->bnd.fMin[0] < 0.5*lx &&       \
			INTRSCT_fDist2 > fBall2) goto label;                        \
		}                                                               \
	else if (INTRSCT_dy1 > 0.0) {                                       \
		INTRSCT_dy += ly;                                               \
		if (INTRSCT_dy < INTRSCT_dy1) {                                 \
			INTRSCT_fDist2 += INTRSCT_dy*INTRSCT_dy;                    \
			sy -= ly;                                                   \
			bPeriodic = 1;                                              \
			}                                                           \
		else {                                                          \
			INTRSCT_fDist2 += INTRSCT_dy1*INTRSCT_dy1;                  \
			}                                                           \
		if ((pkdn)->bnd.fMax[0] - (pkdn)->bnd.fMin[0] < 0.5*lx &&       \
			INTRSCT_fDist2 > fBall2) goto label;                        \
		}                                                               \
	else {                                                              \
		}                                                               \
	INTRSCT_dz = (pkdn)->bnd.fMin[2]-z;                                 \
	INTRSCT_dz1 = z-(pkdn)->bnd.fMax[2];                                \
	if (INTRSCT_dz > 0.0) {                                             \
		INTRSCT_dz1 += lz;                                              \
		if (INTRSCT_dz1 < INTRSCT_dz) {                                 \
			INTRSCT_fDist2 += INTRSCT_dz1*INTRSCT_dz1;                  \
			sz = z+lz;                                                  \
			bPeriodic = 1;                                              \
			}                                                           \
		else {                                                          \
			INTRSCT_fDist2 += INTRSCT_dz*INTRSCT_dz;                    \
			sz = z;                                                     \
			}                                                           \
		if (INTRSCT_fDist2 > fBall2) goto label;                        \
		}                                                               \
	else if (INTRSCT_dz1 > 0.0) {                                       \
		INTRSCT_dz += lz;                                               \
		if (INTRSCT_dz < INTRSCT_dz1) {                                 \
			INTRSCT_fDist2 += INTRSCT_dz*INTRSCT_dz;                    \
			sz = z-lz;                                                  \
			bPeriodic = 1;                                              \
			}                                                           \
		else {                                                          \
			INTRSCT_fDist2 += INTRSCT_dz1*INTRSCT_dz1;                  \
			sz = z;                                                     \
			}                                                           \
		if (INTRSCT_fDist2 > fBall2) goto label;                        \
		}                                                               \
	else {                                                              \
		sz = z;                                                         \
		}                                                               \
	}
#else
#define INTERSECT(pkdn,fBall2,lx,ly,lz,x,y,z,sx,sy,sz,bPeriodic,label)  \
    {                                                                   \
	FLOAT INTRSCT_dx,INTRSCT_dy,INTRSCT_dz;                             \
	FLOAT INTRSCT_dx1,INTRSCT_dy1,INTRSCT_dz1;                          \
    FLOAT INTRSCT_fDist2;                                               \
	INTRSCT_dx = (pkdn)->bnd.fMin[0]-x;                                 \
	INTRSCT_dx1 = x-(pkdn)->bnd.fMax[0];                                \
	INTRSCT_dy = (pkdn)->bnd.fMin[1]-y;                                 \
	INTRSCT_dy1 = y-(pkdn)->bnd.fMax[1];                                \
	INTRSCT_dz = (pkdn)->bnd.fMin[2]-z;                                 \
	INTRSCT_dz1 = z-(pkdn)->bnd.fMax[2];                                \
	if (INTRSCT_dx > 0.0) {                                             \
		INTRSCT_dx1 += lx;                                              \
		if (INTRSCT_dx1 < INTRSCT_dx) {                                 \
			INTRSCT_fDist2 = INTRSCT_dx1*INTRSCT_dx1;                   \
			sx = x+lx;                                                  \
			bPeriodic = 1;                                              \
			}                                                           \
		else {                                                          \
			INTRSCT_fDist2 = INTRSCT_dx*INTRSCT_dx;                     \
			sx = x;                                                     \
			}                                                           \
		if (INTRSCT_fDist2 > fBall2) goto label;                        \
		}                                                               \
	else if (INTRSCT_dx1 > 0.0) {                                       \
		INTRSCT_dx += lx;                                               \
		if (INTRSCT_dx < INTRSCT_dx1) {                                 \
			INTRSCT_fDist2 = INTRSCT_dx*INTRSCT_dx;                     \
			sx = x-lx;                                                  \
			bPeriodic = 1;                                              \
			}                                                           \
		else {                                                          \
			INTRSCT_fDist2 = INTRSCT_dx1*INTRSCT_dx1;                   \
			sx = x;                                                     \
			}                                                           \
		if (INTRSCT_fDist2 > fBall2) goto label;                        \
		}                                                               \
	else {                                                              \
		INTRSCT_fDist2 = 0.0;                                           \
		sx = x;                                                         \
		}                                                               \
	if (INTRSCT_dy > 0.0) {                                             \
		INTRSCT_dy1 += ly;                                              \
		if (INTRSCT_dy1 < INTRSCT_dy) {                                 \
			INTRSCT_fDist2 += INTRSCT_dy1*INTRSCT_dy1;                  \
			sy = y+ly;                                                  \
			bPeriodic = 1;                                              \
			}                                                           \
		else {                                                          \
			INTRSCT_fDist2 += INTRSCT_dy*INTRSCT_dy;                    \
			sy = y;                                                     \
			}                                                           \
		if (INTRSCT_fDist2 > fBall2) goto label;                        \
		}                                                               \
	else if (INTRSCT_dy1 > 0.0) {                                       \
		INTRSCT_dy += ly;                                               \
		if (INTRSCT_dy < INTRSCT_dy1) {                                 \
			INTRSCT_fDist2 += INTRSCT_dy*INTRSCT_dy;                    \
			sy = y-ly;                                                  \
			bPeriodic = 1;                                              \
			}                                                           \
		else {                                                          \
			INTRSCT_fDist2 += INTRSCT_dy1*INTRSCT_dy1;                  \
			sy = y;                                                     \
			}                                                           \
		if (INTRSCT_fDist2 > fBall2) goto label;                        \
		}                                                               \
	else {                                                              \
		sy = y;                                                         \
		}                                                               \
	if (INTRSCT_dz > 0.0) {                                             \
		INTRSCT_dz1 += lz;                                              \
		if (INTRSCT_dz1 < INTRSCT_dz) {                                 \
			INTRSCT_fDist2 += INTRSCT_dz1*INTRSCT_dz1;                  \
			sz = z+lz;                                                  \
			bPeriodic = 1;                                              \
			}                                                           \
		else {                                                          \
			INTRSCT_fDist2 += INTRSCT_dz*INTRSCT_dz;                    \
			sz = z;                                                     \
			}                                                           \
		if (INTRSCT_fDist2 > fBall2) goto label;                        \
		}                                                               \
	else if (INTRSCT_dz1 > 0.0) {                                       \
		INTRSCT_dz += lz;                                               \
		if (INTRSCT_dz < INTRSCT_dz1) {                                 \
			INTRSCT_fDist2 += INTRSCT_dz*INTRSCT_dz;                    \
			sz = z-lz;                                                  \
			bPeriodic = 1;                                              \
			}                                                           \
		else {                                                          \
			INTRSCT_fDist2 += INTRSCT_dz1*INTRSCT_dz1;                  \
			sz = z;                                                     \
			}                                                           \
		if (INTRSCT_fDist2 > fBall2) goto label;                        \
		}                                                               \
	else {                                                              \
		sz = z;                                                         \
		}                                                               \
	}
#endif


#define INTERSECTSCATTERNP(pkdn,x,y,z,label)                            \
    {                                                                   \
    if ( x<(pkdn)->bndBall.fMin[0] || x>(pkdn)->bndBall.fMax[0] ||      \
        y<(pkdn)->bndBall.fMin[1] || y>(pkdn)->bndBall.fMax[1] ||       \
        z<(pkdn)->bndBall.fMin[2] || z>(pkdn)->bndBall.fMax[2] ) goto label; \
	}


#define INTERSECTSCATTER(pkdn,lx,ly,lz,x,y,z,sx,sy,sz,bPeriodic,label)  \
    {                                                                   \
    if (x < (pkdn)->bndBall.fMin[0]) {                                  \
        sx=x+lx;                                                        \
        if (sx > (pkdn)->bndBall.fMax[0]) goto label;                   \
        bPeriodic = 1;                                                  \
        }                                                               \
    else if (x > (pkdn)->bndBall.fMax[0]) {                             \
        sx=x-lx;                                                        \
        if (sx < (pkdn)->bndBall.fMin[0]) goto label;                   \
        bPeriodic = 1;                                                  \
        }                                                               \
    else {                                                              \
        sx=x;                                                           \
        }                                                               \
    if (y < (pkdn)->bndBall.fMin[1]) {                                  \
        sy=y+ly;                                                        \
        if (sy > (pkdn)->bndBall.fMax[1]) goto label;                   \
        bPeriodic = 1;                                                  \
        }                                                               \
    else if (y > (pkdn)->bndBall.fMax[1]) {                             \
        sy=y-ly;                                                        \
        if (sy < (pkdn)->bndBall.fMin[1]) goto label;                   \
        bPeriodic = 1;                                                  \
        }                                                               \
    else {                                                              \
        sy=y;                                                           \
        }                                                               \
    if (z < (pkdn)->bndBall.fMin[2]) {                                  \
        sz=z+lz;                                                        \
        if (sz > (pkdn)->bndBall.fMax[2]) goto label;                   \
        bPeriodic = 1;                                                  \
        }                                                               \
    else if (z > (pkdn)->bndBall.fMax[2]) {                             \
        sz=z-lz;                                                        \
        if (sz < (pkdn)->bndBall.fMin[2]) goto label;                   \
        bPeriodic = 1;                                                  \
        }                                                               \
    else {                                                              \
        sz=z;                                                           \
        }                                                               \
	}

#define DTINTERSECTNP(pkdn,dt2,rMin2,x,y,z,vx,vy,vz,label)              \
    {                                                                   \
	FLOAT DTIS_dx,DTIS_dy,DTIS_dz;                                      \
	FLOAT DTIS_dx1,DTIS_dy1,DTIS_dz1;                                   \
	FLOAT DTIS_dvx,DTIS_dvy,DTIS_dvz;                                   \
	FLOAT DTIS_dvx1,DTIS_dvy1,DTIS_dvz1;                                \
	FLOAT DTIS_dv2, DTIS_MaxDist2;                                      \
    DTIS_dvx = fabs((pkdn)->bndDt.vMin[0]-vx);                          \
    DTIS_dvx1 = fabs((pkdn)->bndDt.vMax[0]-vx);                         \
    DTIS_dv2 = (DTIS_dvx > DTIS_dvx1 ? DTIS_dvx*DTIS_dvx : DTIS_dvx1*DTIS_dvx1); \
    DTIS_dvy = fabs((pkdn)->bndDt.vMin[1]-vy);                          \
    DTIS_dvy1 = fabs((pkdn)->bndDt.vMax[1]-vy);                         \
    DTIS_dv2 += (DTIS_dvy > DTIS_dvy1 ? DTIS_dvy*DTIS_dvy : DTIS_dvy1*DTIS_dvy1); \
    DTIS_dvz = fabs((pkdn)->bndDt.vMin[2]-vz);                          \
    DTIS_dvz1 = fabs((pkdn)->bndDt.vMax[2]-vz);                         \
    DTIS_dv2 += (DTIS_dvz > DTIS_dvz1 ? DTIS_dvz*DTIS_dvz : DTIS_dvz1*DTIS_dvz1); \
	DTIS_MaxDist2 = (DTIS_dv2 + (pkdn)->bndDt.cMax*(pkdn)->bndDt.cMax)*dt2; \
	DTIS_dx = (pkdn)->bnd.fMin[0]-x;                                    \
	DTIS_dx1 = x-(pkdn)->bnd.fMax[0];                                   \
	DTIS_dy = (pkdn)->bnd.fMin[1]-y;                                    \
	DTIS_dy1 = y-(pkdn)->bnd.fMax[1];                                   \
	DTIS_dz = (pkdn)->bnd.fMin[2]-z;                                    \
	DTIS_dz1 = z-(pkdn)->bnd.fMax[2];                                   \
	if (DTIS_dx > 0.0) rMin2 = DTIS_dx*DTIS_dx;                         \
	else if (DTIS_dx1 > 0.0) rMin2 = DTIS_dx1*DTIS_dx1;                 \
	else rMin2 = 0.0;                                                   \
	if (DTIS_dy > 0.0) rMin2 += DTIS_dy*DTIS_dy;                        \
	else if (DTIS_dy1 > 0.0) rMin2 += DTIS_dy1*DTIS_dy1;                \
	if (DTIS_dz > 0.0) rMin2 += DTIS_dz*DTIS_dz;                        \
	else if (DTIS_dz1 > 0.0) rMin2 += DTIS_dz1*DTIS_dz1;                \
    if (rMin2 > DTIS_MaxDist2) goto label;                              \
	}

//	if ((x+.5)*(x+.5)+y*y+z*z < 0.033*0.033) printf("MAXDIST2 %g %g %g %f %f  %f  ",(pkdn)->bndDt.cMax,(pkdn)->bndDt.vMin[0],(pkdn)->bndDt.vMax[0],DTIS_dv2,dt2,DTIS_MaxDist2); \

#define DTINTERSECT(pkdn,dt2,rMin2,lx,ly,lz,x,y,z,sx,sy,sz,bPeriodic,vx,vy,vz,label) \
    {                                                                   \
	FLOAT DTIS_dx,DTIS_dy,DTIS_dz;                                      \
	FLOAT DTIS_dx1,DTIS_dy1,DTIS_dz1;                                   \
	FLOAT DTIS_dvx,DTIS_dvy,DTIS_dvz;                                   \
	FLOAT DTIS_dvx1,DTIS_dvy1,DTIS_dvz1;                                \
	FLOAT DTIS_dv2, DTIS_MaxDist2;                                      \
    DTIS_dvx = fabs((pkdn)->bndDt.vMin[0]-vx);                          \
    DTIS_dvx1 = fabs((pkdn)->bndDt.vMax[0]-vx);                         \
    DTIS_dv2 = (DTIS_dvx > DTIS_dvx1 ? DTIS_dvx*DTIS_dvx : DTIS_dvx1*DTIS_dvx1); \
    DTIS_dvy = fabs((pkdn)->bndDt.vMin[1]-vy);                          \
    DTIS_dvy1 = fabs((pkdn)->bndDt.vMax[1]-vy);                         \
    DTIS_dv2 += (DTIS_dvy > DTIS_dvy1 ? DTIS_dvy*DTIS_dvy : DTIS_dvy1*DTIS_dvy1); \
    DTIS_dvz = fabs((pkdn)->bndDt.vMin[2]-vz);                          \
    DTIS_dvz1 = fabs((pkdn)->bndDt.vMax[2]-vz);                         \
    DTIS_dv2 += (DTIS_dvz > DTIS_dvz1 ? DTIS_dvz*DTIS_dvz : DTIS_dvz1*DTIS_dvz1); \
	DTIS_MaxDist2 = (DTIS_dv2 + (pkdn)->bndDt.cMax*(pkdn)->bndDt.cMax)*dt2; \
	DTIS_dx = (pkdn)->bnd.fMin[0]-x;                                    \
	DTIS_dx1 = x-(pkdn)->bnd.fMax[0];                                   \
	DTIS_dy = (pkdn)->bnd.fMin[1]-y;                                    \
	DTIS_dy1 = y-(pkdn)->bnd.fMax[1];                                   \
	DTIS_dz = (pkdn)->bnd.fMin[2]-z;                                    \
	DTIS_dz1 = z-(pkdn)->bnd.fMax[2];                                   \
	if (DTIS_dx > 0.0) {                                                \
		DTIS_dx1 += lx;                                                 \
		if (DTIS_dx1 < DTIS_dx) {                                       \
			rMin2 = DTIS_dx1*DTIS_dx1;                                  \
			sx = x+lx;                                                  \
			bPeriodic = 1;                                              \
			}                                                           \
		else {                                                          \
			rMin2 = DTIS_dx*DTIS_dx;                                    \
			sx = x;                                                     \
			}                                                           \
		if (rMin2 > DTIS_MaxDist2) goto label;                          \
		}                                                               \
	else if (DTIS_dx1 > 0.0) {                                          \
		DTIS_dx += lx;                                                  \
		if (DTIS_dx < DTIS_dx1) {                                       \
			rMin2 = DTIS_dx*DTIS_dx;                                    \
			sx = x-lx;                                                  \
			bPeriodic = 1;                                              \
			}                                                           \
		else {                                                          \
			rMin2 = DTIS_dx1*DTIS_dx1;                                  \
			sx = x;                                                     \
			}                                                           \
		if (rMin2 > DTIS_MaxDist2) goto label;                          \
		}                                                               \
	else {                                                              \
		rMin2 = 0.0;                                                    \
		sx = x;                                                         \
		}                                                               \
	if (DTIS_dy > 0.0) {                                                \
		DTIS_dy1 += ly;                                                 \
		if (DTIS_dy1 < DTIS_dy) {                                       \
			rMin2 += DTIS_dy1*DTIS_dy1;                                 \
			sy = y+ly;                                                  \
			bPeriodic = 1;                                              \
			}                                                           \
		else {                                                          \
			rMin2 += DTIS_dy*DTIS_dy;                                   \
			sy = y;                                                     \
			}                                                           \
		if (rMin2 > DTIS_MaxDist2) goto label;                          \
		}                                                               \
	else if (DTIS_dy1 > 0.0) {                                          \
		DTIS_dy += ly;                                                  \
		if (DTIS_dy < DTIS_dy1) {                                       \
			rMin2 += DTIS_dy*DTIS_dy;                                   \
			sy = y-ly;                                                  \
			bPeriodic = 1;                                              \
			}                                                           \
		else {                                                          \
			rMin2 += DTIS_dy1*DTIS_dy1;                                 \
			sy = y;                                                     \
			}                                                           \
		if (rMin2 > DTIS_MaxDist2) goto label;                          \
		}                                                               \
	else {                                                              \
		sy = y;                                                         \
		}                                                               \
	if (DTIS_dz > 0.0) {                                                \
		DTIS_dz1 += lz;                                                 \
		if (DTIS_dz1 < DTIS_dz) {                                       \
			rMin2 += DTIS_dz1*DTIS_dz1;                                 \
			sz = z+lz;                                                  \
			bPeriodic = 1;                                              \
			}                                                           \
		else {                                                          \
			rMin2 += DTIS_dz*DTIS_dz;                                   \
			sz = z;                                                     \
			}                                                           \
		if (rMin2 > DTIS_MaxDist2) goto label;                          \
		}                                                               \
	else if (DTIS_dz1 > 0.0) {                                          \
		DTIS_dz += lz;                                                  \
		if (DTIS_dz < DTIS_dz1) {                                       \
			rMin2 += DTIS_dz*DTIS_dz;                                   \
			sz = z-lz;                                                  \
			bPeriodic = 1;                                              \
			}                                                           \
		else {                                                          \
			rMin2 += DTIS_dz1*DTIS_dz1;                                 \
			sz = z;                                                     \
			}                                                           \
		if (rMin2 > DTIS_MaxDist2) goto label;                          \
		}                                                               \
	else {                                                              \
		sz = z;                                                         \
		}                                                               \
	}

#define DTESTIMATOR(pkdn,dt_est,sx,sy,sz,vx,vy,vz)                      \
    {                                                                   \
	FLOAT DTEST_dx,DTEST_dy,DTEST_dz;                                   \
	FLOAT DTEST_dvx,DTEST_dvy,DTEST_dvz;                                \
	FLOAT DTEST_dvx1,DTEST_dvy1,DTEST_dvz1;                             \
	FLOAT DTEST_fDist, DTEST_dvdr;                                      \
	DTEST_dx = sx - (pkdn)->r[0];                                       \
	DTEST_dvx = (vx-(pkdn)->bndDt.vMin[0])*DTEST_dx;                    \
    DTEST_dvx1 = (vx-(pkdn)->bndDt.vMax[0])*DTEST_dx;                   \
    DTEST_dvdr = (DTEST_dvx < DTEST_dvx1 ? DTEST_dvx : DTEST_dvx1);     \
    DTEST_dy = sy - (pkdn)->r[1];                                       \
	DTEST_dvy = (vy-(pkdn)->bndDt.vMin[1])*DTEST_dy;                    \
    DTEST_dvy1 = (vy-(pkdn)->bndDt.vMax[1])*DTEST_dy;                   \
    DTEST_dvdr += (DTEST_dvy < DTEST_dvy1 ? DTEST_dvy : DTEST_dvy1);    \
	DTEST_dz = sz - (pkdn)->r[2];                                       \
	DTEST_dvz = (vz-(pkdn)->bndDt.vMin[2])*DTEST_dz;                    \
    DTEST_dvz1 = (vz-(pkdn)->bndDt.vMax[2])*DTEST_dz;                   \
    DTEST_dvdr += (DTEST_dvz < DTEST_dvz1 ? DTEST_dvz : DTEST_dvz1);    \
    DTEST_fDist = sqrt(DTEST_dx*DTEST_dx + DTEST_dy*DTEST_dy + DTEST_dz*DTEST_dz); \
    dt_est = DTEST_fDist/((pkdn)->bndDt.cMax - (DTEST_dvdr < 0 ? DTEST_dvdr/DTEST_fDist : 0)); }

void smInitList(SMX smx, int nSmooth);
void smFinishList(SMX smx);
void smGrowList(SMX smx);
int smInitialize(SMX *,PKD,SMF *,int,int,int,int,int,double);
void smFinish(SMX smx,SMF *smf, CASTAT *pcs);
void smLargefBallCheck(SMX smx,PARTICLE *p,FLOAT lx, FLOAT ly, FLOAT lz);

PQ *smLoadPQ(SMX smx, int pi, int nSmooth, int nTree, FLOAT x, FLOAT y, FLOAT z,  FLOAT lx, FLOAT ly, FLOAT lz);
PQ *smRecentrePQ(SMX smx, PQ *pqn, int nSmooth, FLOAT x, FLOAT y, FLOAT z);
PQ *smBallSearch(SMX smx,PQ *pq,FLOAT *ri,int *cpStart);
PQ *smBallSearchNP(SMX smx,PQ *pq,FLOAT *ri,int *cpStart);
PQ *smBallSearchAll(SMX smx, PQ *pq, int pi, FLOAT x, FLOAT y, FLOAT z,  FLOAT lx, FLOAT ly, FLOAT lz);
void smLoadNNFromPQ(NN *nnList,PQ *pqi);

int smBallGather(SMX smx,FLOAT fBall2,FLOAT *ri);
int smBallGatherNP(SMX smx,FLOAT fBall2,FLOAT *ri,int cp);
void smBallScatter(SMX smx,FLOAT *ri,int iMarkType);
void smBallScatterNP(SMX smx,FLOAT *ri,int iMarkType,int cp);
void smDtBall(SMX smx,FLOAT *ri,FLOAT *vi,FLOAT *pdt2,FLOAT fBall2);
void smDtBallNP(SMX smx,FLOAT *ri,FLOAT *vi,FLOAT *pdt2,FLOAT fBall2,int cp);

int smResmoothParticle(SMX smx,SMF *smf,int pi,FLOAT fBall2,int cp,FLOAT lx,FLOAT ly,FLOAT lz);

void smSmooth(SMX,SMF *);
void smMarkSmooth(SMX,SMF *, int);
void smReSmooth(SMX,SMF *);
void smDtSmooth(SMX,SMF *);

#ifdef OLD_KEPLER
void smQQSmooth(SMX smx, SMF *smf);
#endif

#endif
