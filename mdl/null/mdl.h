#ifndef MDL_HINCLUDED
#define MDL_HINCLUDED
#include <stdio.h>
#include <assert.h>

#define MDL_VERSION_NUMBER 2.21

#if defined(__osf__) || defined(__sgi)
#define vsnprintf(a,b,c,d) vsprintf((a),(c),(d))
#endif

#define SRV_STOP		0

typedef struct workSpace {
    char *cWorkList;         /* Pointer to this PEs work list */
    int iWorkEltSize;        /* Number of bytes per work list element */
    int nWorkElts;           /* Number of elemets in cWorkList */
    int iNextLocalElt;       /* The next element for the local PE to work on */
    int iNextRemoteElt;      /* The next element for a remote PE to work on, if requested */
    int nLocalWorkRemaining; /* Total number of unassigned work elements remaining */
} WORK;

typedef struct cacheSpace {
	int iType;
	char *pData;
	int iDataSize;
	int nData;
	void (*init)(void *);
	void (*combine)(void *,void *);
	/*	
	 ** Statistics stuff.
	 */
	int nAccess;
	int nAccHigh;
	long nMiss;
	long nColl;
	long nMin;
	} CACHE;


typedef struct serviceRec {
	int nInBytes;
	int nOutBytes;
	void *p1;
	void (*fcnService)(void *,void *,int,void *,int *);	
	} SERVICE;


typedef struct mdlContext {
	int nThreads;
	int idSelf;
	int bDiag;
	int *atid;
	FILE *fpDiag;
	/*
	 ** Services stuff!
	 */
	int nMaxServices;
	int nMaxInBytes;
	int nMaxOutBytes;
	SERVICE *psrv;
	char *pszIn;
	char *pszOut;
	/*
	 ** Caching stuff!
	 */
	unsigned long uRand;
	int iMaxDataSize;
	int nMaxCacheIds;
	CACHE *cache;
    /*
    ** Work management stuff!
    */
    WORK work;

	} * MDL;


/*
** These are for reduction operations, which need MPI datatypes
*/
/* MPI REDUCTION OPERATIONS */
enum ntropy_reduction {
    MDL_REDUCE_MAX,
    MDL_REDUCE_MIN,
    MDL_REDUCE_SUM,
    MDL_REDUCE_PROD,
    MDL_REDUCE_LAND,
    MDL_REDUCE_BAND,
    MDL_REDUCE_LOR,
    MDL_REDUCE_BOR,
    MDL_REDUCE_LXOR,
    MDL_REDUCE_BXOR,
    MDL_REDUCE_MAXLOC,
    MDL_REDUCE_MINLOC
};

/* MPI DATATYPES */
enum ntropy_datatypes {
    MDL_TYPE_INT,
    MDL_TYPE_LONG,
    MDL_TYPE_SHORT,
    MDL_TYPE_UNSIGNED_SHORT,
    MDL_TYPE_UNSIGNED,
    MDL_TYPE_UNSIGNED_LONG,
    MDL_TYPE_FLOAT,
    MDL_TYPE_DOUBLE,
    MDL_TYPE_LONG_DOUBLE,
    MDL_TYPE_BYTE
};


/*
 * MDL debug and Timer macros and prototypes 
 */
/* 
 * Compile time mdl debugging options
 *
 * mdl asserts: define MDLASSERT
 * Probably should always be on unless you want no mdlDiag output at all
 *
 * NB: defining NDEBUG turns off all asserts so MDLASSERT will not assert
 * however it will output uding mdlDiag and the code continues.
 */
#define MDLASSERT
/* 
 * Debug functions active: define MDLDEBUG
 * Adds debugging mdldebug prints and mdldebugassert asserts
 */
#define MDLDEBUG
/* 
 * Timer functions active: define MDLTIMER
 * Makes mdl timer functions active
 */
#ifndef _CRAYMPP
#define MDLTIMER
#endif


void mdlprintf( MDL mdl, const char *format, ... );

#ifdef MDLASSERT
#ifndef __STRING
#define __STRING( arg )   (("arg"))
#endif
#define mdlassert(mdl,expr) \
    { \
      if (!(expr)) { \
             mdlprintf( mdl, "%s:%d Assertion `%s' failed.\n", __FILE__, __LINE__, __STRING(expr) ); \
             assert( expr ); \
             } \
    }
#else
#define mdlassert(mdl,expr)  assert(expr)
#endif

#ifdef MDLDEBUG
#define mdldebugassert(mdl,expr)   mdlassert(mdl,expr)
void mdldebug( MDL mdl, const char *format, ... );
#else
#define mdldebug
#define mdldebugassert
#endif

typedef struct {
  double wallclock;
  double cpu;
  double system;
} mdlTimer;

#ifdef MDLTIMER
void mdlZeroTimer(MDL mdl,mdlTimer *);
void mdlGetTimer(MDL mdl,mdlTimer *,mdlTimer *);
void mdlPrintTimer(MDL mdl,char *message,mdlTimer *);
#else
#define mdlZeroTimer
#define mdlGetTimer
#define mdlPrintTimer
#endif

/*
 ** General Functions
 */

double mdlVersion(MDL);
double mdlCpuTimer(MDL);
int mdlInitialize(MDL *,char **,void (*)(MDL));
void mdlFinish(MDL);
int mdlThreads(MDL);
int mdlSelf(MDL);
int mdlSwap(MDL,int,size_t,void *,size_t,size_t *,size_t *);
void mdlDiag(MDL,char *);
void mdlAddService(MDL,int,void *,void (*)(void *,void *,int,void *,int *),
				   int,int);
void mdlReqService(MDL,int,int,void *,int);
void mdlGetReply(MDL,int,void *,int *);
void mdlHandler(MDL);
/*
 ** Caching functions.
 */
void *mdlMalloc(MDL,int);
void *mdlMallocMax(MDL mdl, size_t iSize, int *iMaxSize);
void *mdlMallocShared(MDL mdl, size_t iSize, int *iEltSize);
void mdlFree(MDL,void *);
void mdlROcache(MDL,int,void *,int,int);
void mdlCOcache(MDL,int,void *,int,int,
				void (*)(void *),void (*)(void *,void *));
void mdlDUMcache(MDL,int);
void mdlFinishCache(MDL,int);
void mdlCacheCheck(MDL);
void *mdlAquire(MDL,int,int,int);
void mdlRelease(MDL,int,void *);
/*
 ** Cache statistics functions.
 */
double mdlNumAccess(MDL,int);
double mdlMissRatio(MDL,int);
double mdlCollRatio(MDL,int);
double mdlMinRatio(MDL,int);
double mdlWaitReplace(MDL,int);
double mdlWaitFlush(MDL,int);
/* 
** Work management functions.
*/
void mdlInitWork(MDL mdl, void *pWorkList, int iWorkEltSize, int nWorkElts, int bDynamic);
void mdlFinishWork(MDL mdl, void *pWorkList);
void *mdlRequestWork(MDL mdl, void *pWorkList, int *pidHome);
/*
** Collectives
*/
void mdlCollectShared(MDL mdl, void *, int);
void mdlBroadcast(MDL mdl, int iRoot, void *vBuf, int iEltSize);
void mdlAllReduce(MDL mdl, int iType, int iReduce, void *pSendArray, 
                  void *pReceiveArray, int iEltSize, int nElements);
void mdlGather(MDL mdl, void *pSentElt, void *pReceiveArray, int iEltSize);
void *mdlGatherv(MDL mdl, void *vOutElt, int *iEltInSizes, int iEltOutSize,
                 int *piInArrayStride);
void mdlScatter(MDL mdl, void *vOutArray, void *vInElt, int iEltSize);
void mdlAllToAll(MDL mdl, void *pSendArray, void *pReceiveArray, int iEltSize);
int mdlComputeEltSizeFromType(MDL mdl, int iType);

#endif

