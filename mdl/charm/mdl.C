/*
 ** Charm++ version of MDL.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <malloc.h>
#include <math.h>
#include <limits.h>
#include <assert.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <values.h>
#include "mdl.h"
#include "mdlimpl.h"

CProxy_AMdl aId;
CProxy_Main MainId;
CProxy_grpCache CacheId;

#define MDL_NOCACHE			0
#define MDL_ROCACHE			1
#define MDL_COCACHE			2

#define MDL_DEFAULT_BYTES		80000
#define MDL_DEFAULT_SERVICES	50
#define MDL_DEFAULT_CACHEIDS	5

#define MDL_TRANS_SIZE		50000

void _srvNull(void *p1,void *vin,int nIn,void *vout,int *pnOut)
{
	return;
	}

double mdlCpuTimer(MDL mdl)
{
#ifndef _CRAYMPP
	struct rusage ru;

	getrusage(0,&ru);
	return((double)ru.ru_utime.tv_sec + 1e-6*(double)ru.ru_utime.tv_usec);
#else
	return( ((double) clock())/CLOCKS_PER_SEC);
#endif
	}

/* 
 * MDL debug and Timer functions 
 */
#define MDLPRINTF_STRING_MAXLEN 256
void mdlprintf( MDL mdl, const char *format, ... )
{
     static char ach[MDLPRINTF_STRING_MAXLEN];
     va_list args;

     if (mdl->bDiag) {	
         va_start( args, format);
         vsnprintf( ach, MDLPRINTF_STRING_MAXLEN, format, args);
         mdlDiag( mdl, ach);
         va_end( args);
         }
}

#ifdef MDLDEBUG
void mdldebug( MDL mdl, const char *format, ... )
{
     static char ach[MDLPRINTF_STRING_MAXLEN];
     va_list args;

     if (mdl->bDiag) {	
         va_start( args, format);
	 vsnprintf( ach, MDLPRINTF_STRING_MAXLEN, format, args);
	 mdlDiag( mdl, ach);
	 va_end( args);
         }
}
#endif

#ifdef MDLTIMER
void mdlZeroTimer(MDL mdl, mdlTimer *t)
{
  struct timezone tz;
  struct timeval tv;
  struct rusage ru;
  tz.tz_minuteswest = 0;
  tz.tz_dsttime = 0;
  gettimeofday(&tv,&tz);
  t->wallclock = tv.tv_sec + 1e-6*(double) tv.tv_usec;
  getrusage(0,&ru);
  t->cpu = (double)ru.ru_utime.tv_sec + 1e-6*(double)ru.ru_utime.tv_usec;
  t->system = (double)ru.ru_stime.tv_sec + 1e-6*(double)ru.ru_stime.tv_usec;
}

void mdlGetTimer(MDL mdl, mdlTimer *t0, mdlTimer *t)
{
  struct timezone tz;
  struct timeval tv;
  struct rusage ru;

  getrusage(0,&ru);
  t->cpu = (double)ru.ru_utime.tv_sec + 1e-6*(double)ru.ru_utime.tv_usec - t0->cpu;
  t->system = (double)ru.ru_stime.tv_sec + 1e-6*(double)ru.ru_stime.tv_usec - t0->system;
  tz.tz_minuteswest = 0;
  tz.tz_dsttime = 0;
  gettimeofday(&tv,&tz);
  t->wallclock = tv.tv_sec + 1e-6*(double) tv.tv_usec - t0->wallclock;
}

void mdlPrintTimer(MDL mdl,char *message, mdlTimer *t0) 
{
  mdlTimer lt;

  if (mdl->bDiag) {	
      mdlGetTimer(mdl,t0,&lt);
      mdlprintf(mdl,"%s %f %f %f\n",message,lt.wallclock,lt.cpu,lt.system);
      }
}
#endif

/*
 * Charm start-up strategy:
 * The mainchare Main() is called first.
 * It will then call the main() for PKDGRAV which will call 
 * mdlInitialize()
 * Slight change: a threaded routine has to start main to keep things
 * from blocking startMain() serves this purpose.
 * mdlInitialize will then call proxies on other processors which will
 * invoke fcnChild()
 * mdlInitialize will then return so that the main() in PKDGRAV will
 * continue.
 */


// Plan here: have a helper function do all the MDL structure setup
// For main thread return it via pmdl here
// For child threads pass it to fcnChild in the AMdl() constructor

// here is the helper functoin
void mdlSetup(MDL *pmdl, int bDiag, const char* progname)
{
	MDL mdl;
	int i;
	char *p,ach[256],achDiag[256];
    
	*pmdl = NULL;
	mdl = (mdlContext *) malloc(sizeof(struct mdlContext));
	assert(mdl != NULL);
	/*
	 ** Set default "maximums" for structures. These are NOT hard
	 ** maximums, as the structures will be realloc'd when these
	 ** values are exceeded.
	 */
	mdl->nMaxServices = MDL_DEFAULT_SERVICES;
	mdl->nMaxSrvBytes = MDL_DEFAULT_BYTES;
	/*
	 ** Now allocate the initial service slots.
	 */
	mdl->psrv = (SERVICE *)malloc(mdl->nMaxServices*sizeof(SERVICE));
	assert(mdl->psrv != NULL);
	/*
	 ** Initialize the new service slots.
	 */
	for (i=0;i<mdl->nMaxServices;++i) {
		mdl->psrv[i].p1 = NULL;
		mdl->psrv[i].nInBytes = 0;
		mdl->psrv[i].nOutBytes = 0;
		mdl->psrv[i].fcnService = NULL;
		}
	/*
	 ** Provide a 'null' service for sid = 0, so that stopping the 
	 ** service handler is well defined!
	 */
	mdl->psrv[0].p1 = NULL;
	mdl->psrv[0].nInBytes = 0;
	mdl->psrv[0].nOutBytes = 0;
	mdl->psrv[0].fcnService = _srvNull;
	/*
	 ** Allocate service buffers.
	 */
	mdl->pszOut = (char *)malloc(mdl->nMaxSrvBytes+sizeof(SRVHEAD));
	assert(mdl->pszOut != NULL);

	/*
	 ** Do some low level argument parsing for number of threads, and
	 ** diagnostic flag!
	 */
	if (bDiag) {
	    p = getenv("MDL_DIAGNOSTIC");
	    if (!p) p = getenv("HOME");
	    if (!p) sprintf(ach,"/tmp");
	    else sprintf(ach,"%s",p);
	    }

	mdl->bDiag = bDiag;
	*pmdl = mdl;
	if (mdl->bDiag) {
		sprintf(achDiag,"%s/%s.%d",ach,progname,mdl->idSelf);
		mdl->fpDiag = fopen(achDiag,"w");
		assert(mdl->fpDiag != NULL);
	    }
    }

/*
 * Count in bit-reversed order
 */
static
int revadd(int i, int nbits) {
        int result;
        int ibit;

        result = i;
        for(ibit = nbits - 1; ibit >= 0; ibit--) {
             result = result ^ (1 << ibit);
             if(result & (1 << ibit))
                break;
        }
        return result;
}

class treeMap : public CkArrayMap
{
private:
public:
    int *pmap;
    treeMap(int nThreads) 
    {
	int ilog2; int rev;
	int i;
	int nPes = CkNumPes();
	int nThdPe = nThreads/nPes;
	
	// find number of bits needed
	for(i = nThreads-1, ilog2 = 0; i > 0; i >>= 1)
		ilog2++;
           
	pmap = new int[nThreads];
	
	rev = 0;
	for (i=0; i < nThreads; ++i) {
	    pmap[i] = rev/nThdPe;
	    while((rev = revadd(rev, ilog2)) >= nThreads);
	    // pmap[i] = i;
	    }
    }
    ~treeMap() { delete[] pmap; }
    
    int procNum(int , const CkArrayIndex &element) 
    {
	int myPe = *(element.data());
	    
	return pmap[myPe];
    }
};

Main::Main(CkArgMsg* m)
{
      char **argv = m->argv;
	int i,bDiag,bThreads;

	/*
	 ** Do some low level argument parsing for number of threads, and
	 ** diagnostic flag!
	 */
	bDiag = 0;
	bThreads = 0;
	i = 1;
	nThreads = 1;
	while (argv[i]) {
		if (!strcmp(argv[i],"-sz") && !bThreads) {
			++i;
			if (argv[i]) {
			    bThreads = 1;
			    nThreads = atoi(argv[i]);
			    }
			}
		if (!strcmp(argv[i],"-d") && !bDiag) {
			bDiag = 1;
			}
		++i;
		}
	if (bThreads) {
		fprintf(stderr,"Warning: -sz refers to threads per processor\n");
		fflush(stderr);
		}

	char *tmp = strrchr(argv[0],'/');
	if (!tmp) tmp = argv[0];
	else ++tmp;

#if 0
	// LiveVis stuff -polling version
	liveVisConfig cfg(true, false);
	CkArrayOptions opts(CkNumPes());
	liveVisPollInit(cfg, opts);
#endif
	
      nfinished = 0;
      MainId = thishandle;

      nThreads *= CkNumPes();
      
      CProxy_treeMap procMap = CProxy_treeMap::ckNew(nThreads);
      CkArrayOptions opts(nThreads);
      opts.setMap(procMap);
      CacheId = CProxy_grpCache::ckNew();
      aId = CProxy_AMdl::ckNew(bDiag, std::string(tmp), opts);
      
#if 0
	// LiveVis stuff  non- polling version
	liveVisConfig cfg(true, false);
	liveVisInit(cfg, aId, drawing_callback);
#endif

      delete m;
      MainId.startMain();
    };

void
Main::startMain()
{
      // get command-line arguments
      int argc = CkGetArgc();
      char **argv = CkGetArgv();


      AMPI_Main(argc, argv);
};

void
Main::done(void)
  {
      nfinished++;
      
      if(nfinished >= nThreads)
	  CkExit();
  };

int
grpCache::elementRegister(int index)
{
    int thisRank;
    
    CmiLock(lock);
    assert(nElem < MAXELEM);
    thisRank = nElem;
    vecIndex[nElem] = index;
    nElem++;
    CmiUnlock(lock);
    return thisRank;
    }

AMdl::AMdl(int bDiag, const std::string& progname)
{
    int i;
    
	mdlSetup(&mdl, bDiag, progname.c_str());
	mdl->nThreads = numElements;
	mdl->idSelf = thisIndex;
	mdl->iNodeSelf = CkMyNode();
	mdl->pSelf = this;
	
	cbService = (CkCallback **) malloc(mdl->nThreads*sizeof(CkCallback *));
	for(i = 0; i < mdl->nThreads; i++)
	    cbService[i] = NULL;
	nInBar = 0;
	swapData.id = -1;

	CProxy_grpCache proxyCache(CacheId);
	iMyRank = proxyCache.ckLocalBranch()->elementRegister(thisIndex);
	vecIndex = proxyCache.ckLocalBranch()->vecIndex;
    }

void
AMdl::AMdlInit(void *fcnPtr)
    {
	
	void (*fcnChild)(MDL) = (void (*)(MDL)) fcnPtr;
	
	if(thisIndex == 0) return;
	(*fcnChild)(mdl);
	mdlFinish(mdl);
    }

extern "C"
int mdlInitialize(MDL *pmdl,char **argv,void (*fcnChild)(MDL))
{

	void *fcnPtr = (void *)fcnChild;
	
	CProxy_AMdl amdlProxy(aId);

	*pmdl = amdlProxy[0].ckLocal()->mdl;

	for(int i = 0; i < (*pmdl)->nThreads; i++) {
	    amdlProxy[i].AMdlInit(fcnPtr);
	    }
	
	
	return (*pmdl)->nThreads;
	
    }

// Creation routines for variable size messages.
void* MdlMsg::alloc(int mnum, size_t size, int *sizes, int priobits){
    int total_size = size + sizes[0] * sizeof(char);
    MdlMsg * mesg = (MdlMsg *)CkAllocMsg(mnum, total_size, priobits);
    mesg->pszBuf = (char *)((char*)mesg + sizeof(MdlMsg));
    return (void *)mesg;
}

void* MdlMsg::pack(MdlMsg *mesg){
    return (void *) mesg;
}

MdlMsg* MdlMsg::unpack(void *buf){
    MdlMsg *mesg = (MdlMsg*)buf;
    mesg->pszBuf = (char *)((char*)mesg + sizeof(MdlMsg));
    return mesg;
}

// Creation routines for swap messages.
void* MdlSwapMsg::alloc(int mnum, size_t size, int *sizes, int priobits){
    int total_size = size + sizes[0] * sizeof(char);
    MdlSwapMsg * mesg = (MdlSwapMsg *)CkAllocMsg(mnum, total_size, priobits);
    mesg->pszBuf = (char *)((char*)mesg + sizeof(MdlSwapMsg));
    return (void *)mesg;
}

void* MdlSwapMsg::pack(MdlSwapMsg *mesg){
    return (void *) mesg;
}

MdlSwapMsg* MdlSwapMsg::unpack(void *buf){
    MdlSwapMsg *mesg = (MdlSwapMsg*)buf;
    mesg->pszBuf = (char *)((char*)mesg + sizeof(MdlSwapMsg));
    return mesg;
}

// Creation routines for cache messages.
void* MdlCacheMsg::alloc(int mnum, size_t size, int *sizes, int priobits){
    int total_size = size + sizes[0] * sizeof(char);
    MdlCacheMsg * mesg = (MdlCacheMsg *)CkAllocMsg(mnum, total_size, priobits);
    mesg->pszBuf = (char *)((char*)mesg + sizeof(MdlCacheMsg));
    return (void *)mesg;
}

void* MdlCacheMsg::pack(MdlCacheMsg *mesg){
    return (void *) mesg;
}

MdlCacheMsg* MdlCacheMsg::unpack(void *buf){
    MdlCacheMsg *mesg = (MdlCacheMsg*)buf;
    mesg->pszBuf = (char *)((char*)mesg + sizeof(MdlCacheMsg));
    return mesg;
}

// Creation routines for cache flush all messages.
void* MdlCacheFlshMsg::alloc(int mnum, size_t size, int *sizes, int priobits){
    int total_size = size + sizes[0] * sizeof(int)
	+ sizes[0] * sizes[1] * sizeof(char);
    MdlCacheFlshMsg * mesg = (MdlCacheFlshMsg *)CkAllocMsg(mnum, total_size,
						       priobits);
    mesg->nLines = sizes[0];
    mesg->pLine = (int *)((char*)mesg + sizeof(MdlCacheFlshMsg));
    mesg->pszBuf = (char *)((char*)mesg + sizeof(MdlCacheFlshMsg)
			    + sizes[0]*sizeof(int));
    return (void *)mesg;
}

void* MdlCacheFlshMsg::pack(MdlCacheFlshMsg *mesg){
    return (void *) mesg;
}

MdlCacheFlshMsg* MdlCacheFlshMsg::unpack(void *buf){
    MdlCacheFlshMsg *mesg = (MdlCacheFlshMsg*)buf;
    mesg->pLine = (int *)((char*)mesg + sizeof(MdlCacheFlshMsg));
    mesg->pszBuf = (char *)((char*)mesg + sizeof(MdlCacheFlshMsg)
			    + mesg->nLines*sizeof(int));
    return mesg;
}

extern "C"
void mdlFinish(MDL mdl)
{
	/*
	 ** Close Diagnostic file.
	 */
	if (mdl->bDiag) {
		fclose(mdl->fpDiag);
		}
	/*
	 ** Deallocate storage.
	 */
	free(mdl->psrv);
	free(mdl->pszOut);
	free(mdl);
	
	CProxy_Main proxyMain(MainId);
	proxyMain.done();
	CthSuspend(); 
	}


/*
 ** This function returns the number of threads in the set of 
 ** threads.
 */
extern "C"
int mdlThreads(MDL mdl)
{
	return(mdl->nThreads);
	}


/*
 ** This function returns this threads 'id' number within the specified
 ** MDL Context. Parent thread always has 'id' of 0, where as children
 ** have 'id's ranging from 1..(nThreads - 1).
 */
extern "C"
int mdlSelf(MDL mdl)
{
	return(mdl->idSelf);
	}


/*
 ** This is a tricky function. It initiates a bilateral transfer between
 ** two threads. Both threads MUST be expecting this transfer. The transfer
 ** occurs between idSelf <---> 'id' or 'id' <---> idSelf as seen from the
 ** opposing thread. It is designed as a high performance non-local memory
 ** swapping primitive and implementation will vary in non-trivial ways 
 ** between differing architectures and parallel paradigms (eg. message
 ** passing and shared address space). A buffer is specified by 'pszBuf'
 ** which is 'nBufBytes' in size. Of this buffer the LAST 'nOutBytes' are
 ** transfered to the opponent, in turn, the opponent thread transfers his
 ** nBufBytes to this thread's buffer starting at 'pszBuf'.
 ** If the transfer completes with no problems the function returns 1.
 ** If the function returns 0 then one of the players has not received all
 ** of the others memory, however he will have successfully transfered all
 ** of his memory.
 */
// Plan here:
// invoke method to give "OutBytes" and "BufBytes"
// That method invokes sendMore()
// Which invokes GetMore()
// repeat until done

extern "C"
int mdlSwap(MDL mdl,int id,
	    size_t nBufBytes,	// available free space on my processor
	    void *vBuf,		// beginning of buffer space
	    size_t nOutBytes,	// Bytes I want to transfer
	    size_t *pnSndBytes, // bytes actually sent
	    size_t *pnRcvBytes  // bytes received
	    )
{
	AMdl *mmdl = mdl->pSelf;
	char *pszBuf = (char *) vBuf;
        CProxy_AMdl proxyAmdl(aId);
	
	mmdl->cbSwap = new CkCallback(CkCallback::resumeThread);
	mmdl->swapData.nOutBytes = nOutBytes;
	mmdl->swapData.nBufBytes = nBufBytes;
	mmdl->swapData.nRcvBytes = 0;
	mmdl->swapData.nSndBytes = 0;
	mmdl->swapData.id = id;
	mmdl->swapData.pszIn = pszBuf;
	mmdl->swapData.pszOut = &pszBuf[nBufBytes - nOutBytes];
	mmdl->swapData.done = 0;
	
	assert(nBufBytes >= nOutBytes);
	proxyAmdl[id].swapInit(nOutBytes, nBufBytes);
	(mmdl->cbSwap)->thread_delay();
	mmdl->swapData.id = -1;
	delete mmdl->cbSwap;

	*pnRcvBytes = mmdl->swapData.nRcvBytes;
	*pnSndBytes = mmdl->swapData.nSndBytes;

	if (mmdl->swapData.nOutBytes) return(0);
	else if (mmdl->swapData.nInBytes) return(0);
	else return(1);
    }

void
AMdl::swapInit(size_t nInBytes, size_t nBufBytes)
{
    swapData.nInBytes = nInBytes;
    swapData.nOutBufBytes = nBufBytes;
    
    while(swapData.id == -1)
	CthYield();
    
    swapSendMore();
    }

void
AMdl::swapSendMore() 
{
    CProxy_AMdl proxyAmdl(aId);

    if(swapData.nOutBytes && swapData.nOutBufBytes) {
	int nOutMax = (swapData.nOutBytes < MDL_TRANS_SIZE)
	    ? swapData.nOutBytes : MDL_TRANS_SIZE;
	nOutMax = (nOutMax < swapData.nOutBufBytes)
	    ? nOutMax : swapData.nOutBufBytes;

	assert(nOutMax > 0);
	
	MdlSwapMsg *mesg = new(&nOutMax, 0) MdlSwapMsg;
	mesg->nBytes = nOutMax;
	
	memcpy(mesg->pszBuf, swapData.pszOut, nOutMax);
	
	// Adjust counts for next iteration
	swapData.pszOut += nOutMax;
	swapData.nOutBytes -= nOutMax;
	swapData.nOutBufBytes -= nOutMax;
	swapData.nSndBytes += nOutMax;
	
	proxyAmdl[swapData.id].swapGetMore(mesg);
	}
    else {
	proxyAmdl[thisIndex].swapDone();
	proxyAmdl[swapData.id].swapDone();
	}
    }

void
AMdl::swapGetMore(MdlSwapMsg *mesg) 
{
    CProxy_AMdl proxyAmdl(aId);
    int nBytes = mesg->nBytes;	// temporary for bytes transferred
    
    while(swapData.pszIn + nBytes > swapData.pszOut) {
	CthYield();		// pause while buffer is transferred out
	}

    memcpy(swapData.pszIn, mesg->pszBuf, nBytes);

    swapData.pszIn += nBytes;
    swapData.nInBytes -= nBytes;
    swapData.nBufBytes -= nBytes;
    swapData.nRcvBytes += nBytes;

    delete mesg;
    proxyAmdl[swapData.id].swapSendMore();
    }

void
AMdl::swapDone()
{
    swapData.done++;
    if(swapData.done == 2) {
	cbSwap->send();
	}
    
    }
	
extern "C"
void mdlDiag(MDL mdl,char *psz)
{
	if (mdl->bDiag) {	
		fputs(psz,mdl->fpDiag);
		fflush(mdl->fpDiag);
		}
	}

extern "C"
void mdlAddService(MDL mdl,int sid,void *p1,
				   void (*fcnService)(void *,void *,int,void *,int *),
				   int nInBytes,int nOutBytes)
{
	int i,nMaxServices,nMaxBytes;

	assert(sid > 0);
	if (sid >= mdl->nMaxServices) {
		/*
		 ** reallocate service buffer, adding space for 8 new services
		 ** including the one just defined.
		 */
		nMaxServices = sid + 9;
		mdl->psrv = (SERVICE *) realloc(mdl->psrv,
						nMaxServices*sizeof(SERVICE));
		assert(mdl->psrv != NULL);
		/*
		 ** Initialize the new service slots.
		 */
		for (i=mdl->nMaxServices;i<nMaxServices;++i) {
			mdl->psrv[i].p1 = NULL;
			mdl->psrv[i].nInBytes = 0;
			mdl->psrv[i].nOutBytes = 0;
			mdl->psrv[i].fcnService = NULL;
			}
		mdl->nMaxServices = nMaxServices;
		}
	/*
	 ** Make sure the service buffers are big enough!
	 */
	nMaxBytes = (nInBytes > nOutBytes)?nInBytes:nOutBytes;
	if (nMaxBytes > mdl->nMaxSrvBytes) {
		mdl->pszOut = (char *) realloc(mdl->pszOut,nMaxBytes+sizeof(SRVHEAD));
		assert(mdl->pszOut != NULL);
		mdl->nMaxSrvBytes = nMaxBytes;
		}
	mdl->psrv[sid].p1 = p1;
	mdl->psrv[sid].nInBytes = nInBytes;
	mdl->psrv[sid].nOutBytes = nOutBytes;
	mdl->psrv[sid].fcnService = fcnService;
	}


extern "C"
void mdlReqService(MDL mdl,int id,int sid,void *vin,int nInBytes)
{
	char *pszIn = (char *) vin;
	MdlMsg *mesg;
	int i;
	
	mesg = new(&nInBytes, sizeof(int)) MdlMsg;
	// *((int *)CkPriorityPtr(mesg)) = MAXINT;
	// CkSetQueueing(mesg, CK_QUEUEING_IFIFO);
	
	mesg->ph.idFrom = mdl->idSelf;
	mesg->ph.sid = sid;
	if (!pszIn) mesg->ph.nInBytes = 0;
	else mesg->ph.nInBytes = nInBytes;
	if (nInBytes > 0 && pszIn != NULL) {
		for (i=0;i<nInBytes;++i) mesg->pszBuf[i] = pszIn[i];
		}
	CProxy_AMdl aProxy(aId);
	
	assert(mdl->pSelf->cbService[id] == NULL);

	mdl->pSelf->cbService[id] = new CkCallback(CkCallback::resumeThread);
	
	aProxy[id].reqHandle(mesg);
	}

extern "C"
void mdlGetReply(MDL mdl,int id,void *vout,int *pnOutBytes)
{
	char *pszOut = (char *) vout;
	int i, nOutBytes;
	MdlMsg *mesg;

	mesg = (MdlMsg *)(mdl->pSelf->cbService[id])->thread_delay();

	nOutBytes = mesg->ph.nOutBytes;
	
	if (nOutBytes > 0 && pszOut != NULL) {
		for (i=0; i< nOutBytes; ++i) pszOut[i] = mesg->pszBuf[i];
		}
	if (pnOutBytes) *pnOutBytes = nOutBytes;
	delete mesg;
	delete mdl->pSelf->cbService[id];
	mdl->pSelf->cbService[id] = NULL;
	}

void
AMdl::reqReply(MdlMsg *mesg) 
{
    assert(cbService[mesg->ph.idFrom] != NULL);
    
    cbService[mesg->ph.idFrom]->send(mesg);
}


// In the CHARM case, this doesn't do much:
// It needs to just wait for the SRV_STOP function to happen, then return.

extern "C"
void mdlHandler(MDL mdl)
{
    mdl->pSelf->waitSrvStop();
    }

void
AMdl::waitSrvStop()
{
    threadSrvWait = CthSelf();
    CthSuspend();
}

void
AMdl::stopSrv()
{
    assert(threadSrvWait != NULL);
    CthAwaken(threadSrvWait);
    }

void
AMdl::reqHandle(MdlMsg * mesg) 
{
	char *pszIn = mesg->pszBuf;
	int sid,id,nInBytes, nOutBytes;
	char *pszOut = &mdl->pszOut[sizeof(SRVHEAD)];
	int i;

	id = mesg->ph.idFrom;
	sid = mesg->ph.sid;
	nInBytes = mesg->ph.nInBytes;
	mdlassert(mdl, sid < mdl->nMaxServices);
	while(mdl->psrv[sid].fcnService == NULL)
	    CthYield();		// Wait for service to be registered
	
	mdlassert(mdl, nInBytes <= mdl->psrv[sid].nInBytes);
	nOutBytes = 0;
	assert(mdl->psrv[sid].fcnService != NULL);
	
	(*mdl->psrv[sid].fcnService)(mdl->psrv[sid].p1, pszIn,
				     nInBytes, pszOut, &nOutBytes);
	delete mesg;
	assert(nOutBytes <= mdl->psrv[sid].nOutBytes);

	mesg = new(&nOutBytes, 0) MdlMsg;
	mesg->ph.idFrom = mdl->idSelf;
	mesg->ph.sid = sid;
	mesg->ph.nInBytes = nInBytes;
	mesg->ph.nOutBytes = nOutBytes;
	for(i = 0; i < nOutBytes; i++)
	    mesg->pszBuf[i] = pszOut[i];

	CProxy_AMdl proxyAMdl(aId);
	proxyAMdl[id].reqReply(mesg);

	if(sid == SRV_STOP) {
	    // Stop handler
	    CProxy_AMdl proxyAMdl(aId);
	    proxyAMdl[mdl->idSelf].stopSrv();
	    }
	}

#define MDL_CHECK_MASK  	0x7f
#define BILLION				1000000000

grpCache::grpCache()
{
    int i;
    
    lock = CmiCreateLock();
    nElem = 0;
    for(i = 0; i < MAXELEM; i++) vecIndex[i] = -1;
    
    threadBarrier = 0;
    nMaxCacheIds = MDL_DEFAULT_CACHEIDS;

    /*
     ** Allocate initial cache spaces.
     */
    cache = (CACHE *) malloc(nMaxCacheIds*sizeof(CACHE));
    assert(cache != NULL);
    /*
     ** Initialize caching spaces.
     */
    for (i=0;i<nMaxCacheIds;++i) {
	    cache[i].iType = MDL_NOCACHE;
	    }
    }

void
grpCache::CacheRequest(MdlCacheMsg *mesg)
{
	CACHE *c;
	char *pszRpl;
	char *t;
	int i;
	int iLineSize;
        // CProxy_AMdl proxyAMdl(aId);
	CProxy_grpCache proxyCache(CacheId);
	char *pData;
	int nData;

	c = &(cache[mesg->ch.cid]);
	assert(c->iType != MDL_NOCACHE);
	assert(indexRank(mesg->ch.rid) != -1);
	
	MdlCacheMsg *mesgRpl;

	pData = c->procData[indexRank(mesg->ch.rid)].pData;
	nData = c->procData[indexRank(mesg->ch.rid)].nData;
	
	t = &pData[mesg->ch.iLine*c->iLineSize];
	if(t+c->iLineSize > pData + nData*c->iDataSize)
		iLineSize = pData + nData*c->iDataSize - t;
	else
		iLineSize = c->iLineSize;
	mesgRpl = new(&c->iLineSize, 0) MdlCacheMsg;

	pszRpl = mesgRpl->pszBuf;
	mesgRpl->ch.cid = mesg->ch.cid;
	mesgRpl->ch.id = mesg->ch.rid;
	mesgRpl->ch.rid = mesg->ch.id;
	for (i=0;i<iLineSize;++i) pszRpl[i] = t[i];

	proxyCache[mesg->ch.iNode].CacheReply(mesgRpl);
	delete mesg;
    }

void
grpCache::CacheReply(MdlCacheMsg *mesg)
{
	/*
	 ** For now assume no prefetching!
	 ** This means that this WILL be the reply to this Aquire
	 ** request.
	 */

    CProxy_AMdl proxyAMdl(aId);
    proxyAMdl[mesg->ch.rid].unblockCache(mesg);
    }

void
AMdl::unblockCache(MdlCacheMsg *mesg)
{
    cbCache->send(mesg);
    }

void
AMdl::CacheFlush(MdlCacheMsg *mesg)
{
	CACHE *c;
	char *pszRcv = mesg->pszBuf;
	char *t;
	int n,i;
	int iDataSize;
	char *pData;
	int nData;
	CProxy_grpCache proxyCache(CacheId);

	c = &(cache[mesg->ch.cid]);
	
	assert(c->iType == MDL_COCACHE);
	i = mesg->ch.iLine*MDL_CACHELINE_ELTS;
	pData = c->procData[iMyRank].pData;
	nData = c->procData[iMyRank].nData;
	t = &pData[i*c->iDataSize];
	/*
	 ** Make sure we don't combine beyond the number of data elements!
	 */
	n = i + MDL_CACHELINE_ELTS;
	if (n > nData) n = nData;
	n -= i;
	n *= c->iDataSize;
	iDataSize = c->iDataSize;
	for (i=0;i<n;i+=iDataSize) {
		(*c->combine)(&t[i],&pszRcv[i]);
		}
	proxyCache[mesg->ch.id].flushreply();
	delete mesg;
    }

void
AMdl::CacheFlushAll(MdlCacheFlshMsg *mesg)
{
	CACHE *c;
	int *pLine = mesg->pLine;
	char *pszRcv = mesg->pszBuf;
	int nLines = mesg->nLines;
	char *t;
	int n,i;
	int iDataSize;
	int iLine;
	CProxy_grpCache proxyCache(CacheId);
	char *pData;
	int nData;

	c = &(cache[mesg->ch.cid]);
	
	assert(c->iType == MDL_COCACHE);
	pData = c->procData[iMyRank].pData;
	nData = c->procData[iMyRank].nData;

	for(iLine = 0; iLine < nLines; iLine++) {
	    i = pLine[iLine]*MDL_CACHELINE_ELTS;
	    t = &pData[i*c->iDataSize];
	    /*
	     ** Make sure we don't combine beyond the number of data elements!
	     */
	    n = i + MDL_CACHELINE_ELTS;
	    if (n > nData) n = nData;
	    n -= i;
	    n *= c->iDataSize;
	    iDataSize = c->iDataSize;
	    for (i=0;i<n;i+=iDataSize) {
		    (*c->combine)(&t[i],&pszRcv[i]);
		}
	    pszRcv += c->iLineSize;
	    }
	proxyCache[mesg->ch.id].flushreply();
	delete mesg;
    }

/*
 ** Special MDL memory allocation functions for allocating memory 
 ** which must be visible to other processors thru the MDL cache 
 ** functions.
 ** mdlMalloc() is defined to return a pointer to AT LEAST iSize bytes 
 ** of memory. This pointer will be passed to either mdlROcache or 
 ** mdlCOcache as the pData parameter.
 ** For PVM and most machines these functions are trivial, but on the 
 ** T3D and perhaps some future machines these functions are required.
 */
extern "C"
void *mdlMalloc(MDL mdl,size_t iSize)
{	
	return(malloc(iSize));
	}

extern "C"
void mdlFree(MDL mdl,void *p)
{
	free(p);
	}

extern "C"
void *mdlRealloc(MDL mdl,void *p, size_t iSize)
{	
	return(realloc(p, iSize));
	}

/*
 ** Common initialization for all types of caches.
 */
void
grpCache::CacheInitialize(int cid,void *pData,int iDataSize,int nData,
			  void (*init)(void *),void (*combine)(void *,void *),
			  int iRank, int nThreads)
{
    CACHE *c;
    int i,newMaxCacheIds;

    assert(cid >= 0);

    CmiLock(lock);		// single thread at a time

    /*
     ** Allocate more cache spaces if required!
     */
    if (cid >= nMaxCacheIds) {
	    /*
	     ** reallocate cache spaces, adding space for 2 new cache spaces
	     ** including the one just defined.
	     */
	    newMaxCacheIds = cid + 3;
	    cache = (CACHE *) realloc(cache,newMaxCacheIds*sizeof(CACHE));
	    assert(cache != NULL);
	    /*
	     ** Initialize the new cache slots.
	     */
	    for (i=nMaxCacheIds;i<newMaxCacheIds;++i) {
		    cache[i].iType = MDL_NOCACHE;
		    }
	    nMaxCacheIds = newMaxCacheIds;
	    }
    c = &cache[cid];
    if(c->iType == MDL_NOCACHE) { // first time for this group
	c->procData = (PDATA *)malloc(nElem*sizeof(*(c->procData)));
	assert(c->procData != NULL);
	c->iDataSize = iDataSize;
	}
    else {
	assert(c->iDataSize == iDataSize);
	}

    c->procData[iRank].pData = (char *)pData;
    c->procData[iRank].nData = nData;
    if(c->iType == MDL_NOCACHE) { // first time for this group
	c->iLineSize = MDL_CACHELINE_ELTS*c->iDataSize;
	c->iKeyShift = 0;
	while((1 << c->iKeyShift) < nThreads) ++c->iKeyShift;
	c->iIdMask = (1 << c->iKeyShift) - 1;

	if(c->iKeyShift < MDL_CACHELINE_BITS) {
	  /*
	   * Key will be (index & MDL_INDEX_MASK) | id.
	   */
	    c->iInvKeyShift = MDL_CACHELINE_BITS;
	    c->iKeyShift = 0;
	    }
	else {
	  /*
	   * Key will be (index & MDL_INDEX_MASK) << KeyShift | id.
	   */
	    c->iInvKeyShift = c->iKeyShift;
	    c->iKeyShift -= MDL_CACHELINE_BITS;
	    }

	/*
	 ** Determine the number of cache lines to be allocated.
	 */
	c->nLines = (MDL_CACHE_SIZE/c->iDataSize)
	    >> MDL_CACHELINE_BITS;
	c->iLine = 1;
	c->nTrans = 1;
	while(c->nTrans < c->nLines) c->nTrans *= 2;
	c->nTrans *= 4;
	c->iTransMask = c->nTrans-1;
	/*
	 **	Set up the translation table.
	 */
	c->pTrans = (int *) malloc(c->nTrans*sizeof(int));	
	assert(c->pTrans != NULL);
	for (i=0;i<c->nTrans;++i) c->pTrans[i] = 0;
	/*
	 ** Set up the tags. Note pTag[0] is a Sentinel!
	 */
	c->pTag = (CTAG *) malloc(c->nLines*sizeof(CTAG));
	assert(c->pTag != NULL);
	for (i=0;i<c->nLines;++i) {
		c->pTag[i].iKey = -1;	/* invalid */	
		c->pTag[i].nLock = 0;
		c->pTag[i].nLast = -1;	/* invalid */
		c->pTag[i].iLink = 0;
		c->pTag[i].iIdLock = -1; /* invalid */	
		c->pTag[i].bFetching = 0;
		}
	c->pTag[0].nLock = 1;		/* always locked */
	c->pTag[0].nLast = INT_MAX;  	/* always Most Recently Used */
	c->nAccess = 0;
	c->nAccHigh = 0;
	c->nMiss = 0;				/* !!!, not NB */
	c->nColl = 0;				/* !!!, not NB */
	c->nMin = 0;				/* !!!, not NB */	
	c->nKeyMax = 500;				/* !!!, not NB */
	c->pbKey = (char *) malloc(c->nKeyMax);	/* !!!, not NB */
	assert(c->pbKey != NULL);			/* !!!, not NB */
	for (i=0;i<c->nKeyMax;++i) c->pbKey[i] = 0;	/* !!!, not NB */
	/*
	 ** Allocate cache data lines.
	 */
	c->pLine = (char *) malloc(c->nLines*c->iLineSize);
	assert(c->pLine != NULL);
	if(init == NULL)
	    c->iType = MDL_ROCACHE;
	else {
	    c->iType = MDL_COCACHE;
	    c->init = init;
	    c->combine = combine;
	    }
	c->nOut = 0;
	nFlush = 0;
	}
    else {
	if(init == NULL)
	    assert(c->iType == MDL_ROCACHE);
	else
	    assert(c->iType == MDL_COCACHE);
	}
    CmiUnlock(lock);
    }

void
AMdl::barrierEnter(CkReductionMsg *msg)
{

    delete msg;
    
    cbBarrier->send();
    }

void
AMdl::barrier()
{
    // Use a null reduction as a barrier.

    cbBarrier = new CkCallback(CkCallback::resumeThread);
    CkCallback cb(CkIndex_AMdl::barrierEnter(NULL), thisProxy);
    contribute(0, 0,  CkReduction::concat, cb);
    cbBarrier->thread_delay();
    delete cbBarrier;
}

void
grpCache::flushreply()
{
    CProxy_AMdl proxyAMdl(aId);

    CmiLock(lock);		// single thread at a time
    nFlush--;
    if(nFlush == 0) {
	if(threadBarrier)
	    proxyAMdl[idFlushing].waitflushAwaken();
	}
    CmiUnlock(lock);		// single thread at a time
    }

void
AMdl::waitflushAwaken()
{
    CProxy_grpCache proxyCache(CacheId);

    CthAwaken(proxyCache.ckLocalBranch()->threadBarrier);
}

void
grpCache::waitflush()
{
    CmiLock(lock);
    if(nFlush) {
	threadBarrier = CthSelf();
	CmiUnlock(lock);
	CthSuspend();
	threadBarrier = 0;
	}
    else {
	CmiUnlock(lock);
	}
    }

/*
 ** Initialize a Read-Only caching space.
 */
extern "C"
void mdlROcache(MDL mdl,int cid,void *pData,int iDataSize,int nData)
{
	CProxy_grpCache proxyCache(CacheId);

	proxyCache.ckLocalBranch()->CacheInitialize(cid, pData, iDataSize,
						    nData, NULL, NULL,
						    mdl->pSelf->iMyRank,
						    mdl->nThreads);

	mdl->pSelf->cache = 
	    proxyCache.ckLocalBranch()->cache;
	mdl->pSelf->lock = 
	    &(proxyCache.ckLocalBranch()->lock);
	mdl->pSelf->nElem = 
	    proxyCache.ckLocalBranch()->nElem;

	mdl->pSelf->barrier();
	}

/*
 ** Initialize a Combiner caching space.
 */
extern "C"
void mdlCOcache(MDL mdl,int cid,void *pData,int iDataSize,int nData,
				void (*init)(void *),void (*combine)(void *,void *))
{
	CProxy_grpCache proxyCache(CacheId);

	assert(init != NULL);
	assert(combine != NULL);
	proxyCache.ckLocalBranch()->CacheInitialize(cid, pData, iDataSize,
						    nData, init, combine,
						    mdl->pSelf->iMyRank,
						    mdl->nThreads);

	mdl->pSelf->cache = 
	    proxyCache.ckLocalBranch()->cache;
	mdl->pSelf->lock = 
	    &(proxyCache.ckLocalBranch()->lock);
	mdl->pSelf->nElem = 
	    proxyCache.ckLocalBranch()->nElem;

	mdl->pSelf->barrier();
	}

void
grpCache::FinishCache(int cid, int idSelf, int nThreads)
{
    int i,id;
    char *t;
    int j, iKey;
    CACHE *c = &cache[cid];
    CProxy_AMdl proxyAMdl(aId);

    CmiLock(lock);		// single thread at a time
    c->nOut++;
    if(c->nOut < nElem) {
	CmiUnlock(lock);
	return;
	}
    idFlushing = idSelf;
    CmiUnlock(lock);

    if (c->iType == MDL_COCACHE) {
	/*
	 ** Must flush all valid data elements.
	 */
	/*
	 * Do it a processor at a time
	 */
	int idFlush;

	for(idFlush = 0; idFlush < nThreads; idFlush++) {
	    int iLine = 0;

	    /* first count lines */
	    for (i=1;i<c->nLines;++i) {
		iKey = c->pTag[i].iKey;
		if (iKey >= 0) {
		    id = iKey & c->iIdMask;
		    if(id == idFlush)
			iLine++;
		    }
		}
	    if(iLine == 0)
		continue;

	    int iFlushSize[2];
	    iFlushSize[0] = iLine;
	    iFlushSize[1] = c->iLineSize;

	    MdlCacheFlshMsg *mesgFlsh = new(iFlushSize,0) MdlCacheFlshMsg;
	    char *pszFlsh = mesgFlsh->pszBuf;
	    int *pLine = mesgFlsh->pLine;
	    mesgFlsh->ch.cid = cid;
	    mesgFlsh->ch.id = CmiMyNode();

	    iLine = 0;
	    for (i=1;i<c->nLines;++i) {
		    iKey = c->pTag[i].iKey;
		    if (iKey >= 0) {
			    /*
			     ** Flush element since it is valid!
			     */
			id = iKey & c->iIdMask;
			if(id == idFlush) {
			    pLine[iLine] = iKey >> c->iInvKeyShift;
			    t = &c->pLine[i*c->iLineSize];
			    for(j = 0; j < c->iLineSize; ++j)
				pszFlsh[j] = t[j];
			    iLine++;
			    pszFlsh += c->iLineSize;
			    }
			}
		    }

	    assert(mesgFlsh->nLines == iLine);

	    CmiLock(lock);		// single thread at a time
	    nFlush++;
	    CmiUnlock(lock);		// single thread at a time
	    proxyAMdl[idFlush].CacheFlushAll(mesgFlsh);
	    }
	waitflush();
	}
    /*
     ** Free up storage and finish.
     */
    free(c->pTrans);
    free(c->pTag);
    free(c->pbKey);
    free(c->pLine);
    }

extern "C"
void mdlFinishCache(MDL mdl,int cid)
{
	CProxy_grpCache proxyCache(CacheId);

	proxyCache.ckLocalBranch()->FinishCache(cid, mdl->idSelf, mdl->nThreads);
	mdl->pSelf->barrier();
	mdl->pSelf->cache[cid].iType = MDL_NOCACHE;
	mdl->pSelf->barrier();
	}

#ifndef mdlCacheCheck

extern "C"
void mdlCacheCheck(MDL mdl)
{
    //    int dummy = CmiDeliverMsgs(0);
    CthYield();
    
    }
#endif

extern "C"
void *mdlAquire(MDL mdl,int cid,int iIndex,int id)
{
        CProxy_AMdl proxyAMdl(aId);
	CProxy_grpCache proxyCache(CacheId);
	CACHE *c = &(mdl->pSelf->cache[cid]);
	char *pLine;
	int iElt,iLine,i,iKey,iKeyVic,nKeyNew;
	int idVic;
	int iVictim,*pi;
	char ach[80];
	char *pszFlsh;
	CmiNodeLock *lock = mdl->pSelf->lock;
	int bLocalFetch;	// Can we get the line locally?
	char* pLocalAddress;	// Address of line if so
	int peId;
	int iRank = proxyAMdl[mdl->idSelf].ckLocal()->indexRank(id);
	
	if((c->iType == MDL_ROCACHE && iRank != -1) || id == mdl->idSelf) {
	    // It's on node or purely local: use shared memory
	    return(&(c->procData[iRank].pData[iIndex*c->iDataSize]));
	    }
	
	// if (!(c->nAccess & MDL_CHECK_MASK))
	//        mdlCacheCheck(mdl);
	/*
	 ** Determine memory block key value and cache line.
	iLine = iIndex >> MDL_CACHELINE_BITS;
	iKey = iLine*mdl->nThreads + id;
	 */
	iKey = ((iIndex&MDL_INDEX_MASK) << c->iKeyShift)| id;

	/*
	 ** Check for a match!
	 */
	while(CmiTryLock(*lock))
	    CthYield();
	
	++c->nAccess;

	bLocalFetch = 0;
	CTAG *copyTag = NULL;
	
	/*
	 ** Collision chain search.
	 */
	for(i = c->pTrans[iKey & c->iTransMask]; i;
	    i = c->pTag[i].iLink, ++c->nColl) {
	    if (c->pTag[i].iKey == iKey) {
		// Also match on processor for CO cache
		if(c->iType == MDL_ROCACHE
		   || c->pTag[i].nLock == 0
		   || c->pTag[i].iIdLock == mdl->idSelf) {
		    ++c->pTag[i].nLock;
		    c->pTag[i].iIdLock = mdl->idSelf;
		    c->pTag[i].nLast = c->nAccess;
		    CmiUnlock(*lock);
		    while(c->pTag[i].bFetching == 1)
			CthYield();
		    pLine = &c->pLine[i*c->iLineSize];
		    iElt = iIndex & MDL_CACHE_MASK;
		    return(&pLine[iElt*c->iDataSize]);
		    }
		// matched but another processor has it locked
		bLocalFetch = 1;
		copyTag = &c->pTag[i];
		copyTag->nLock++;
		pLocalAddress = &c->pLine[i*c->iLineSize];
		break;
		}
	    }

	/*
	 ** Cache Miss.
	 */
	iLine = iIndex >> MDL_CACHELINE_BITS;
	//
	// Check if it is locally available
	//
	if(c->iType == MDL_COCACHE && iRank != -1) {
	    bLocalFetch = 1;
	    pLocalAddress = &(c->procData[iRank].pData[iLine*c->iLineSize]);
	    }
	    
        if(!bLocalFetch) {
	    
	    int nRequestBytes = 0;	// Just a place holder

	    peId = proxyAMdl.ckLocalBranch()->lastKnown(CkArrayIndex1D(id));
	    MdlCacheMsg *mesg = new(&nRequestBytes, sizeof(int)) MdlCacheMsg;
	    *((int *)CkPriorityPtr(mesg)) = - INT_MAX;
	    CkSetQueueing(mesg, CK_QUEUEING_IFIFO);
	    mesg->ch.cid = cid;
	    mesg->ch.rid = id;
	    mesg->ch.id = mdl->idSelf;
	    mesg->ch.iNode = mdl->iNodeSelf;
	    mesg->ch.iLine = iLine;

	    // CkArrayIndex1D aidxId(id);

	    mdl->pSelf->cbCache = new CkCallback(CkCallback::resumeThread);
	    
	    proxyCache[CkNodeOf(peId)].CacheRequest(mesg);
	    ++c->nMiss;
	    }
	
	/*
	 **	LRU Victim Search!
	 ** If nAccess > BILLION then we reset all LRU counters.
	 ** This *should* be sufficient to prevent overflow of the 
	 ** Access counter, but it *is* cutting corners a bit. 
	 */
	iElt = iIndex & MDL_CACHE_MASK;
	if (c->nAccess > BILLION) {
		for (i=1;i<c->nLines;++i) c->pTag[i].nLast = 0;
		c->nAccess -= BILLION;
		c->nAccHigh += 1;
		}
	iVictim = 0;
	for (i=c->iLine;i<c->nLines;++i) {
	  if (!c->pTag[i].nLock) { 
		c->iLine = i+1;
		iVictim = i;
		goto GotVictim;
	  }
	}
	for (i=1;i<c->iLine;++i) {
	  if (!c->pTag[i].nLock) { 
		c->iLine = i+1;
		iVictim = i;
		goto GotVictim;
	  }
	}
	if (!iVictim) {
		/*
		 ** Cache Failure!
		 */
		sprintf(ach,"MDL CACHE FAILURE: cid == %d, no unlocked lines!\n",cid);
		mdlDiag(mdl,ach);
		assert(0);
		}
GotVictim:
	iKeyVic = c->pTag[iVictim].iKey;
	/*
	 ** 'pLine' will point to the actual data line in the cache.
	 */
	pLine = &c->pLine[iVictim*c->iLineSize];

	if (iKeyVic >= 0) {
		if (c->iType == MDL_COCACHE) {
			/*
			 ** Flush element since it is valid!
			 */
		        idVic = iKeyVic&c->iIdMask;
			MdlCacheMsg *mesgFlsh = new(&c->iLineSize,0) MdlCacheMsg;
			
			pszFlsh = mesgFlsh->pszBuf;
		        mesgFlsh->ch.cid = cid;
			mesgFlsh->ch.id = CkMyNode();
			mesgFlsh->ch.iLine = iKeyVic >> c->iInvKeyShift;
			for(i = 0; i < c->iLineSize; ++i)
			    pszFlsh[i] = pLine[i];
			proxyCache.ckLocalBranch()->nFlush++;
			proxyAMdl[idVic].CacheFlush(mesgFlsh);
			}
		/*
		 ** If valid iLine then "unlink" it from the cache.
		 */
		pi = &c->pTrans[iKeyVic & c->iTransMask];
		while (*pi != iVictim) pi = &c->pTag[*pi].iLink;
		*pi = c->pTag[iVictim].iLink;
		}
	c->pTag[iVictim].iKey = iKey;
	c->pTag[iVictim].nLock = 1;
	c->pTag[iVictim].iIdLock = mdl->idSelf;
	c->pTag[iVictim].nLast = c->nAccess;
	c->pTag[iVictim].bFetching = 1;
	/*
	 **	Add the modified victim tag back into the cache.
	 ** Note: the new element is placed at the head of the chain.
	 */
	pi = &c->pTrans[iKey & c->iTransMask];
	c->pTag[iVictim].iLink = *pi;
	*pi = iVictim;
	/*
	 ** Figure out whether this is a "new" miss.
	 ** This is for statistics only!
	 */
	if (iKey >= c->nKeyMax) {			/* !!! */
		nKeyNew = iKey+500;
		c->pbKey = (char *) realloc(c->pbKey,nKeyNew);
		assert(c->pbKey != NULL);
		for (i=c->nKeyMax;i<nKeyNew;++i) c->pbKey[i] = 0;
		c->nKeyMax = nKeyNew;
		}
	if (!c->pbKey[iKey]) {
		c->pbKey[iKey] = 1;
		++c->nMin;
		}					/* !!! */
	/*
	 ** At this point 'pLine' is the recipient cache line for the 
	 ** data requested from processor 'id'.
	 */

	// I can unlock here since each line in the CO Cache is
	// processor specific.
	// The "bFetching" bit in the RO Cache also lets me unlock here.
	CmiUnlock(*lock);

	char *pszLine;
	MdlCacheMsg *mesg;
	if(!bLocalFetch) {
	    mesg = (MdlCacheMsg *)(mdl->pSelf->cbCache)->thread_delay();
	    assert(mesg->ch.id == id);
	    assert(mesg->ch.cid == cid);
	    delete mdl->pSelf->cbCache;

	    pszLine = mesg->pszBuf;
	    }
	else {
	    pszLine = pLocalAddress;
	    if(copyTag) {
		while(copyTag->bFetching)
		    CthYield();
		}
	    }
	   
	   
	for(i = 0; i < c->iLineSize; i++)
	    pLine[i] = pszLine[i];

	if(!bLocalFetch)
	    delete mesg;
	if(copyTag) {
	    CmiLock(*lock);
	    copyTag->nLock--;
	    CmiUnlock(*lock);
	    }
	
	if (c->iType == MDL_COCACHE && c->init) {
	    /*
	    ** Call the initializer function for all elements in 
	    ** the cache line.
	    */
	    for (i=0;i<c->iLineSize;i+=c->iDataSize) {
		(*c->init)(&pLine[i]);
		}
	    }
	c->pTag[iVictim].bFetching = 0;
	return(&pLine[iElt*c->iDataSize]);
	}

extern "C"
void mdlRelease(MDL mdl,int cid,void *p)
{
	CACHE *c = &(mdl->pSelf->cache[cid]);
	int iLine,iData;
	
	assert(c != NULL);
	
	iLine = ((char *)p - c->pLine) / c->iLineSize;
	/*
	 ** Check if the pointer fell in a cache line, otherwise it
	 ** must have been a local pointer.
	 */
	if (iLine > 0 && iLine < c->nLines) {
	    while(CmiTryLock(*mdl->pSelf->lock))
		  CthYield();
		--c->pTag[iLine].nLock;
		assert(c->pTag[iLine].nLock >= 0);
		CmiUnlock(*mdl->pSelf->lock);
		return;
		}
	if(c->iType == MDL_ROCACHE) { // assume we're correct with
				      // shared memory
		return;
		}
	else {
		char *pData = c->procData[mdl->pSelf->iMyRank].pData;
		int nData = c->procData[mdl->pSelf->iMyRank].nData;

		iData = ((char *)p - pData) / c->iDataSize;
		assert(iData >= 0 && iData < nData);
		}
	}

extern "C"
double mdlNumAccess(MDL mdl,int cid)
{
	CACHE *c = &(mdl->pSelf->cache[cid]);

	return((c->nAccHigh*1e9 + c->nAccess)/CkNodeSize(CkMyNode()));
	}


extern "C"
double mdlMissRatio(MDL mdl,int cid)
{
	CACHE *c = &(mdl->pSelf->cache[cid]);
	double dAccess = c->nAccHigh*1e9 + c->nAccess;
	
	if (dAccess > 0.0) return(c->nMiss/dAccess);
	else return(0.0);
	}


extern "C"
double mdlCollRatio(MDL mdl,int cid)
{
	CACHE *c = &(mdl->pSelf->cache[cid]);
	double dAccess = c->nAccHigh*1e9 + c->nAccess;

	if (dAccess > 0.0) return(c->nColl/dAccess);
	else return(0.0);
	}


extern "C"
double mdlMinRatio(MDL mdl,int cid)
{
	CACHE *c = &(mdl->pSelf->cache[cid]);
	double dAccess = c->nAccHigh*1e9 + c->nAccess;

	if (dAccess > 0.0) return(c->nMin/dAccess);
	else return(0.0);
	}

#include "mdl.def.h"
