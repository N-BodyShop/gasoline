#ifndef MDLIMPL_HINCLUDED
#define MDLIMPL_HINCLUDED

#include "pup_stl.h"

#include "mdl.decl.h"

#define MDL_CACHE_SIZE		16000000
#define MDL_CACHELINE_BITS	4
#define MDL_CACHELINE_ELTS	(1<<MDL_CACHELINE_BITS)
#define MDL_CACHE_MASK		(MDL_CACHELINE_ELTS-1)
#define MDL_INDEX_MASK		(~MDL_CACHE_MASK)

/*
 ** This structure should be "maximally" aligned, with 4 ints it
 ** should align up to at least QUAD word, which should be enough.
 */
typedef struct srvHeader {
	int idFrom;
	int sid;
	int nInBytes;
	int nOutBytes;
	} SRVHEAD;

/*
 ** This structure should be "maximally" aligned, with 4 ints it
 ** should align up to at least QUAD word, which should be enough.
 */
typedef struct cacheHeader {
    int cid;	// Which cache
    int rid;	// proc. id to which request is being made
    int id; 	// id of requesting processor
    int iNode;  // node of requesting processor
    int iLine;
    } CAHEAD;

class MdlMsg : public CMessage_MdlMsg 
{
 public:
    SRVHEAD ph;
    char *pszBuf;

    static void *alloc(int mnum, size_t size, int *sizes, int priobits);
    static void *pack(MdlMsg *msg);
    static MdlMsg *unpack(void *buf);  
    };

class MdlSwapMsg : public CMessage_MdlSwapMsg 
{
 public:
    int nBytes;
    char *pszBuf;

    static void *alloc(int mnum, size_t size, int *sizes, int priobits);
    static void *pack(MdlSwapMsg *msg);
    static MdlSwapMsg *unpack(void *buf);  
    };

class MdlCacheMsg : public CMessage_MdlCacheMsg 
{
 public:
    CAHEAD ch;
    char *pszBuf;

    static void *alloc(int mnum, size_t size, int *sizes, int priobits);
    static void *pack(MdlCacheMsg *msg);
    static MdlCacheMsg *unpack(void *buf);  
    };

class MdlCacheFlshMsg : public CMessage_MdlCacheFlshMsg 
{
 public:
    CAHEAD ch;
    int nLines;
    int *pLine;
    char *pszBuf;

    static void *alloc(int mnum, size_t size, int *sizes, int priobits);
    static void *pack(MdlCacheFlshMsg *msg);
    static MdlCacheFlshMsg *unpack(void *buf);  
    };

extern "C"
void AMPI_Main(int argc, char **);

class Main : public Chare
{
    int nfinished;
    int nThreads;
    
public:
    Main(CkArgMsg* m);
    void startMain();
    void done(void);
};

typedef struct cacheTag {
	int iKey;
	int nLock;
	int nLast;
	int iLink;
	int iIdLock;
	int bFetching;
	} CTAG;

typedef	struct procDATA {
	    char *pData;
	    int nData;
	    } PDATA;

typedef struct cacheSpace {
	int iType;
        PDATA *procData;
	int iDataSize;
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
	void (*init)(void *);
	void (*combine)(void *,void *);
	int nOut;
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
    } CACHE;

#define MAXELEM 32

class grpCache : public NodeGroup 
{
 public:
    int nMaxCacheIds;
    CACHE *cache;
    CmiNodeLock lock;
    CthThreadStruct * threadBarrier;
    int nFlush;
    int idFlushing;
    int nElem;			/* number of array elements in this
				   group */
    int vecIndex[MAXELEM];
    
    grpCache();
    void CacheInitialize(int cid,void *pData,int iDataSize,int nData,
			 void (*init)(void *),void (*combine)(void *,void *),
			 int iRank, int nThread);
    void CacheRequest(MdlCacheMsg *mesg);
    void CacheReply(MdlCacheMsg *mesg);
    MdlCacheMsg *waitCache(int iRank) ;
    void flushreply();
    void waitflush();
    void FinishCache(int cid, int idSelf, int nThreads);
    int elementRegister(int index);
    inline int indexRank(int index) 
	{
	    for(int i = 0; i < nElem; i++) {
		if(index == vecIndex[i])
		    return i;
		}
	    return -1;
	    }
};

void mdlSetup(MDL *pmdl, int bDiag, const char *);

PUPbytes(void *);

// class AMdl : public ArrayElement1D
class AMdl : public CBase_AMdl
{
public:
    struct {			/* state data for mdlSwap() */
	size_t nInBytes;
	size_t nOutBytes;
	size_t nBufBytes;
	size_t nOutBufBytes;
	size_t nRcvBytes;
	size_t nSndBytes;
	int id;
	char *pszOut;
	char *pszIn;
	int done;
	} swapData;
    CthThreadStruct * threadSrvWait;
    CACHE *cache;		/* pointer to nodegroup cache */
    CmiNodeLock *lock;		/* pointer to nodegroup lock */
    CkCallback * cbSwap;
    CkCallback ** cbService;
    CkCallback * cbCache;
    CkCallback * cbBarrier;
    int idReplyWait;
    int nInBar;
    int nFlush;
    int iMyRank;
    int *vecIndex;
    int nElem;
    
    MDL mdl;
    AMdl(int bDiag, const std::string& progname);

    void AMdlInit(void *fcnPtr);
    AMdl(CkMigrateMessage*) {}
    void swapInit(size_t, size_t);
    void swapSendMore();
    void swapGetMore(MdlSwapMsg *);
    void swapDone();
    void waitSrvStop();
    void stopSrv();
    void reqReply(MdlMsg * mesg);
    void reqHandle(MdlMsg * mesg);
    void unblockCache(MdlCacheMsg *);
    void CacheFlush(MdlCacheMsg *mesg);
    void CacheFlushAll(MdlCacheFlshMsg *mesg);
    void barrierRel();
    void barrierEnter(CkReductionMsg *);
    void barrier();
    void waitflush();
    void waitflushAwaken();
    inline int indexRank(int index) 
	{
	    for(int i = 0; i < nElem; i++) {
		if(index == vecIndex[i])
		    return i;
		}
	    return -1;
	    }
#ifdef CHARM_MIGRATE
    void *pst;  // pointer to pst data to be pupped
    void pup(PUP::er& p);
#endif
};

#endif
