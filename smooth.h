#ifndef SMOOTH_HINCLUDED
#define SMOOTH_HINCLUDED

#include "pkd.h"

#define RESMOOTH_SAFE  10
#define PQ_LOAD_FACTOR 0.50


typedef struct pqNode {
	float fKey;
	struct pqNode *pqLoser;
	struct pqNode *pqFromInt;
	struct pqNode *pqFromExt;
	struct pqNode *pqWinner;	/* Only used when building initial tree */
	struct pqNode *pNext;
	struct pqNode *pPrev;
	int p;
	int id;
	PARTICLE *pPart;
	float ax;
	float ay;
	float az;
	} PQ;


typedef struct nNeighbor {
	PARTICLE *pPart;
	float fDist2;
	} NN;

typedef struct smContext {
	PKD pkd;
	int nSmooth;
        int bPeriodic;
	int *piMark;
	int nListSize;
	NN *nnList;
	PQ *pq;
	PQ *pqHead;
	int nHash;
	PQ **pqHash;
	} * SMX;


#define PQ_INIT(pq,n)\
{\
	int j;\
	for (j=0;j<(n);++j) {\
		if (j < 2) (pq)[j].pqFromInt = NULL;\
		else (pq)[j].pqFromInt = &(pq)[j>>1];\
		(pq)[j].pqFromExt = &(pq)[(j+(n))>>1];\
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
	(q) = (pq)[1].pqWinner;\
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

int smInitialize(SMX *,PKD,int,int, int);
void smFinish(SMX);
void smSmooth(SMX,void (*)(SMX,int,int,NN *));
/*
 ** Smoothing functions.
 */
void smDensity(SMX,int,int,NN *);
void smDensitySym(SMX,int,int,NN *);


#endif











