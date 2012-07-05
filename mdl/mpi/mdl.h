#ifndef MDL_HINCLUDED
#define MDL_HINCLUDED
#include <stdio.h>
#include <assert.h>
#include "mpi.h"

#define MDL_VERSION_NUMBER 2.21

#define SRV_STOP		0

#define MDL_CACHE_SIZE		64000000
#define MDL_CACHELINE_BITS	3
#define MDL_CACHELINE_ELTS	(1<<MDL_CACHELINE_BITS)
#define MDL_CACHE_MASK		(MDL_CACHELINE_ELTS-1)
#define MDL_INDEX_MASK		(~MDL_CACHE_MASK)


/*
 * Work management structs
 */

/* Priority queue structs */

/* Size of the PQ array (should be somewhat larger than number of threads)*/
#define MDL_PQHEAP_SIZE  2048

typedef struct
{
    int priority;
    int data;
} Packet;

typedef struct {
	Packet packets[MDL_PQHEAP_SIZE];
	unsigned int size;
} PacketHeap;


/* This stuct is used by every thread to manage its local work elements
 * and to send/receive work related communications. */
typedef struct workSpace {
    char *cWorkList;         /* Pointer to this PEs work list */
    int iWorkEltSize;        /* Number of bytes per work list element */
    int nWorkElts;           /* Number of elemets in cWorkList */
    int iNextLocalElt;       /* The next element for the local PE to work on */
    int iNextRemoteElt;      /* The next element for a remote PE to work on, if requested */
    int nLocalWorkRemaining; /* Total number of unassigned work elements remaining */
    int idScheduler;         /* MPI id of the master scheduler node */
    char *pbufRemoteWork;    /* Buffer containing previously fetched remote work assignments */
    int nbufRemoteWork;      /* Number of elements in remote work buffer */
    int ibufRemoteWork;      /* "Index" for the current element in remote work buffer */
    int idRemoteWork;        /* The thread that the remote work came from */
    char *pbufIncomingWork;  /* Temporary buffer to store incoming work assignments */
    int nbufIncomingWork;    /* Number of elements in incoming work buffer */
    int retIncomingWork;     /* The mdlWorkReceive return code for the incoming work buffer */
    int idIncomingWork;      /* The thread that the incoming work came from */
    /* MPI Buffers, handles, etc */
    int iBufSize;            /* Size (in bytes) of the buffer to send/receive work comms */
    char *pbufSnd;           /* Buffer for regular sends */
    char *pbufRcv;           /* Buffer in which to receive all MDL_TAG_WORK messages */
    MPI_Request handleRcv;   /* MPI handle for the request dumping into pbufRcv */
    char *pbufISnd;          /* Buffer for Isends */
    MPI_Request handleISnd;  /* MPI handle for Isends */
    int bRequestedMoreWork;  /* Have I send out a request for more work already? */
    int iRidCurrent;         /* The current outstanding work request ID */
    int bDynamic;            /* Do we use dynamic load balancing or not? */
} WORK;

/* This struct is only used by the master scheduler thread.  All other
 * threads simply set bIAmScheduler to 0. */
typedef struct scheduleSpace {
    int bIAmScheduler;
    int nCheckIn;
    int nCheckOut;
    int *iWorkRemainingList; /* An array 0..nthreads-1 of local work elements that
                              * each thread has still unassigned */
    PacketHeap pWorkRemainingHeap;
} SCHEDULE;

typedef struct workHeader {
    int rid; /* Request ID */
    int oid; /* Request originator ID (not necessarily the thread that sent this message) */
    int id;  /* The actual thread that sent this specific message */
    int nWorkUnits; /* The number of work units contained in this message */
    int nLocalWorkRemaining; /* The work remaining for thread "id" */
    int iMessage; /* The "message" */
} WORKHEAD;
    
typedef struct cacheTag {
	int iKey;
	int nLock;
	int nLast;
	int iLink;
	} CTAG;

/*
 ** This structure should be "maximally" aligned, with 4 ints it
 ** should align up to at least QUAD word, which should be enough.
 */
typedef struct cacheHeader {
	int cid;   /* Cache ID */
	int mid;   /* Message ID */
	int id;    /* ID of sender */
	int iLine; /* Line requested */
	} CAHEAD;


typedef struct cacheSpace {
	int iType;
	char *pData;
	int iDataSize;
	int nData;
	int iLineSize;
	int nLines;
    int iLine;
	int nTrans;
	int iTransMask;
        int iKeyShift;
        int iInvKeyShift;
        int iIdMask;
	int *pTrans;
	CTAG *pTag;
	char *pLine;
	int nCheckIn;
	int nCheckOut;
	CAHEAD caReq;
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
	int nKeyMax;
	char *pbKey;
    /* Timers */
    double dTimerWaitReplace;
    double dTimerWaitFlush;
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
	FILE *fpDiag;
	int dontcare;
	int allgrp;
	/*
	 ** Services stuff!
	 */
	int nMaxServices;
	int nMaxSrvBytes;
	SERVICE *psrv;
	char *pszIn;
	char *pszOut;
	char *pszBuf;
	/*
	 ** Swapping buffer.
	 */
	char *pszTrans;
	/*
	 ** Caching stuff!
	 */
	unsigned long uRand;
	int iMaxDataSize;
	int iCaBufSize;
	char *pszRcv;
	int *pmidRpl;
	MPI_Request *pReqRpl;
	MPI_Request ReqRcv;
	char **ppszRpl;
	char *pszFlsh;
	int nMaxCacheIds;
	CACHE *cache;
    /*
    ** Work management stuff!
    */
    WORK work;
    SCHEDULE sch;
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
 * however it will output using mdlDiag and the code continues.
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
#define MDLTIMER


void mdlprintf( MDL mdl, const char *format, ... );

#ifdef MDLASSERT
#ifdef __ANSI_CPP__
#define mdlassert(mdl,expr) \
    { \
      if (!(expr)) { \
             mdlprintf( mdl, "%s:%d Assertion `%s' failed.\n", __FILE__, __LINE__, # expr ); \
             assert( expr ); \
             } \
    }
#else
#define mdlassert(mdl,expr) \
    { \
      if (!(expr)) { \
             mdlprintf( mdl, "%s:%d Assertion `%s' failed.\n", __FILE__, __LINE__, "expr" ); \
             assert( expr ); \
             } \
    }
#endif
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

/* MDL Timer struct */
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
void *mdlMalloc(MDL,size_t);
void *mdlMallocMax(MDL mdl, size_t, int *);
void *mdlMallocShared(MDL mdl, size_t, int *);
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
int  mdlWorkReceive(MDL mdl, char *pWork);
void mdlInitWork(MDL mdl, void *pWorkList, int iWorkEltSize, int nWorkElts, int bDynamic);
void mdlFinishWork(MDL mdl, void *pWorkList);
void *mdlRequestWork(MDL mdl, void *pWorkList, int *pidHome);
/*
** Collectives
*/
int mdlComputeEltSizeFromType(MDL mdl, int iType);
void mdlCollectShared(MDL mdl, void *, int);
void mdlBroadcast(MDL mdl, int iRoot, void *vBuf, int iEltSize);
void mdlGather(MDL mdl, void *vOutElt, void *vInArray, int iEltSize);
void *mdlGatherv(MDL mdl, void *vOutElt, int *iEltInSizes, int iEltOutSize,
                 int *piInArrayStride);
void mdlScatter(MDL mdl, void *vOutArray, void *vInElt, int iEltSize);
void mdlAllReduce(MDL mdl, int iType, int iReduce, void *pSendArray, 
                  void *pReceiveArray, int iEltSize, int nElements);
void mdlAllToAll(MDL mdl, void *pScatterArray, void *pGatherArray, int iEltSize);

#endif
