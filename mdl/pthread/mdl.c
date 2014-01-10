/*
 ** A very low level Machine model. For homogeneous systems only!
 ** This is the pthreads mdl module.
 **
 */
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <limits.h>
#include <stdarg.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <pthread.h>
#include "mdl.h"

#if (__i486__)
#include <linux/smp.h>		/* FOR LINUX ONLY!!! */
#else
int smp_num_cpus = 1;
#endif

#define MDL_NOCACHE			0
#define MDL_ROCACHE			1
#define MDL_COCACHE			2


#define MDL_DEFAULT_BYTES		4096
#define MDL_DEFAULT_SERVICES	50
#define MDL_DEFAULT_CACHEIDS	5

#define MDL_TRANS_SIZE			50000 

/*
 ** GLOBAL BARRIER VARIABLES! All threads must see these.
 */
int MDLnInBar = 0;
int MDLnEpisode = 0;
pthread_mutex_t MDLmuxBar = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t MDLsigBar = PTHREAD_COND_INITIALIZER;

void mdlBarrier(MDL mdl)
{
    int iEpisode;
    
    pthread_mutex_lock(&MDLmuxBar);
    iEpisode = MDLnEpisode;
    ++MDLnInBar;
    if (MDLnInBar == mdl->nThreads) {
	++MDLnEpisode;
	MDLnInBar = 0;
	pthread_cond_broadcast(&MDLsigBar);
	}
    else {
	while(MDLnEpisode == iEpisode) {
	    pthread_cond_wait(&MDLsigBar,&MDLmuxBar);
	    }
	}
    pthread_mutex_unlock(&MDLmuxBar);
    }

void srvNull(void *p1,void *vin,int nIn,void *vout,int *pnOut)
{
    return;
    }


double mdlCpuTimer(MDL mdl)
{
	struct rusage ru;

	getrusage(0,&ru);
	return((double)ru.ru_utime.tv_sec + 1e-6*(double)ru.ru_utime.tv_usec);
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
 ** This function performs basic initialization, common to all
 ** MDL contexts.
 */
void BasicInit(MDL mdl)
{
    int i;

    /*
     ** Set default "maximums" for structures. These are NOT hard
     ** maximums, as the structures will be realloc'd when these
     ** values are exceeded.
     */
    mdl->nMaxServices = MDL_DEFAULT_SERVICES;
    mdl->nMaxInBytes = MDL_DEFAULT_BYTES;
    mdl->nMaxOutBytes = MDL_DEFAULT_BYTES;
    mdl->nMaxCacheIds = MDL_DEFAULT_CACHEIDS;
    /*
     ** Now allocate the initial service slots.
     */
    mdl->psrv = malloc(mdl->nMaxServices*sizeof(SERVICE));
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
    mdl->psrv[0].fcnService = srvNull;
    /*
     ** Allocate initial cache spaces.
     */
    mdl->cache = malloc(mdl->nMaxCacheIds*sizeof(CACHE));
    assert(mdl->cache != NULL);
    /*
     ** Initialize caching spaces.
     */
    for (i=0;i<mdl->nMaxCacheIds;++i) {
		mdl->cache[i].iType = MDL_NOCACHE;
		} 
    /*
     ** Initialize the mailboxes.
     */
    mdl->mbxOwn.pszIn = malloc(mdl->nMaxInBytes);
    assert(mdl->mbxOwn.pszIn != NULL);
    mdl->mbxOwn.pszOut = malloc(mdl->nMaxOutBytes);
    assert(mdl->mbxOwn.pszOut != NULL);
    mdl->mbxOwn.bReq = 0;
    mdl->mbxOwn.bRpl = 0;
    /*
     ** Will have to release once the handler is entered, this
     ** way we make sure that all services are set up before
     ** any thread can start transferring. The main reason for
     ** this change is that the buffers are dynamically 
     ** sized in mdlAddService().
     */
    mdl->mbxOwn.bRel = 0;
    pthread_mutex_init(&mdl->mbxOwn.mux,NULL);
	pthread_cond_init(&mdl->mbxOwn.sigReq,NULL);
	pthread_cond_init(&mdl->mbxOwn.sigRpl,NULL);
	pthread_cond_init(&mdl->mbxOwn.sigRel,NULL);
    /*
     ** Initialize the swapbox, this buffer has fixed size but must
     ** still be malloc'd to ensure that no allignment problems can
     ** occur.
     */
    mdl->swxOwn.pszBuf = malloc(MDL_TRANS_SIZE);
    assert(mdl->swxOwn.pszBuf != NULL);
    mdl->swxOwn.bRec = 0;
    pthread_mutex_init(&mdl->swxOwn.mux,NULL);
    pthread_cond_init(&mdl->swxOwn.sigRec,NULL);
    pthread_cond_init(&mdl->swxOwn.sigSnd,NULL);
    /*
    ** Initialize Cache mailboxes.
    */
    mdl->iMaxDataSize = 0;
    mdl->iCaBufSize = sizeof(CAHEAD);
    pthread_mutex_init(&mdl->muxRing,NULL);
    mdl->iRingHd = 0;
    mdl->iRingTl = 0;
    for(i = 0; i < MDL_MBX_RING_SZ; i++) {
	mdl->mbxCache[i].pszIn = malloc(mdl->iCaBufSize);
	assert(mdl->mbxCache[i].pszIn != NULL);
	mdl->mbxCache[i].bRel = 1;
	mdl->mbxCache[i].bReq = 0;
	pthread_mutex_init(&mdl->mbxCache[i].mux,NULL);
	pthread_cond_init(&mdl->mbxCache[i].sigReq,NULL);
	pthread_cond_init(&mdl->mbxCache[i].sigRel,NULL);
	}
    }


void BasicDestroy(MDL mdl)
{
    int i;
    
    pthread_mutex_destroy(&mdl->mbxOwn.mux);
	pthread_cond_destroy(&mdl->mbxOwn.sigReq);
	pthread_cond_destroy(&mdl->mbxOwn.sigRpl);
	pthread_cond_destroy(&mdl->mbxOwn.sigRel);
    pthread_mutex_destroy(&mdl->swxOwn.mux);
	pthread_cond_destroy(&mdl->swxOwn.sigRel);
	pthread_cond_destroy(&mdl->swxOwn.sigRec);
	pthread_cond_destroy(&mdl->swxOwn.sigSnd);
    pthread_mutex_destroy(&mdl->muxRing);
    for(i = 0; i < MDL_MBX_RING_SZ; i++) {
	pthread_mutex_destroy(&mdl->mbxCache[i].mux);
	pthread_cond_destroy(&mdl->mbxCache[i].sigReq);
	pthread_cond_destroy(&mdl->mbxCache[i].sigRel);
	free(mdl->mbxCache[i].pszIn);
	}
    
    free(mdl->swxOwn.pszBuf);
    free(mdl->psrv);
    free(mdl->mbxOwn.pszIn);
    free(mdl->mbxOwn.pszOut);
    free(mdl->cache);
    }


int mdlInitialize(MDL *pmdl,char **argv,void (*fcnChild)(MDL))
{
    MDL mdl,tmdl;
    int i,nThreads=1,bThreads,bDiag;
    char *p,ach[256],achDiag[256];
    pthread_attr_t attr;
#ifdef TINY_PTHREAD_STACK
    static int first = 1;
    static MDL *save_mdl;
    extern void *(main)(void *);
    
    if(!first) {
        *pmdl = *save_mdl;
        return (*pmdl)->nThreads;
	}
    first = 0;
#endif
    
    *pmdl = NULL;
    /*
     ** Do some low level argument parsing for number of threads, and
     ** diagnostic flag!
     */
    bDiag = 0;
    bThreads = 0;
    i = 1;
    while (argv[i]) {
		if (!strcmp(argv[i],"-sz") && !bThreads) {
			++i;
			if (argv[i]) {
				nThreads = atoi(argv[i]);
				bThreads = 1;
				}
			}
		if (!strcmp(argv[i],"+d") && !bDiag) {
			p = getenv("MDL_DIAGNOSTIC");
			if (!p) p = getenv("HOME");
			if (!p) sprintf(ach,"/tmp");
			else sprintf(ach,"%s",p);
			bDiag = 1;
			}
		++i;
		}
    if (!bThreads) {
		/*
		 ** Default the number of threads to the number of CPUs in the 
		 ** SMP system. See #include <linux/smp.h> above, that
		 ** will vary with system type. Currently defaults to 1 on 
		 ** systems not running LINUX.
		 */
		nThreads = smp_num_cpus;
		}
    mdl = malloc(sizeof(struct mdlContext));
    assert(mdl != NULL);
    mdl->pmdl = malloc(nThreads*sizeof(MDL));
    assert(mdl->pmdl != NULL);
    mdl->pmdl[0] = mdl;			/* that's me! */
    mdl->pt = (pthread_t *)malloc(nThreads*sizeof(pthread_t));
    assert(mdl->pt != NULL);
    *pmdl = mdl;
    pthread_attr_init(&attr);
#ifdef TINY_PTHREAD_STACK
    /*
     * Create 1 Meg stack.
     */
    pthread_attr_setstacksize(&attr, 10*1024*1024);
#endif
    if (nThreads > 1) {
		for (i=1;i<nThreads;++i) {
			/*
			 ** Allocate the children mdl data structures.
			 */
			tmdl = malloc(sizeof(struct mdlContext));
			assert(tmdl != NULL);
			mdl->pmdl[i] = tmdl;
			tmdl->pt = NULL;
			}
		for (i=0;i<nThreads;++i) {
			/*
			 ** Set up all the mdl data structures.
			 */
			tmdl = mdl->pmdl[i];
            BasicInit(tmdl);
			tmdl->pmdl = mdl->pmdl;
			tmdl->idSelf = i;
			tmdl->bDiag = bDiag;
			tmdl->nThreads = nThreads;
			if (tmdl->bDiag) {
				char *tmp = strrchr(argv[0],'/');
				if (!tmp) tmp = argv[0];
				else ++tmp;
				sprintf(achDiag,"%s/%s.%d",ach,tmp,tmdl->idSelf);
				tmdl->fpDiag = fopen(achDiag,"w");
				assert(tmdl->fpDiag != NULL);
				}
			}
		for (i=1;i<nThreads;++i) {
			/*
			 ** Start all the children.
			 */
			pthread_create(&mdl->pt[i],&attr,
				       (void *(*)(void *))fcnChild,
				       mdl->pmdl[i]);
			}
		}
    else {
		/*
		 ** A unik!
		 */
        BasicInit(mdl);
		mdl->bDiag = bDiag;
		mdl->nThreads = 1;
		mdl->idSelf = 0;
		if (mdl->bDiag) {
			char *tmp = strrchr(argv[0],'/');
			if (!tmp) tmp = argv[0];
			else ++tmp;
			sprintf(achDiag,"%s/%s.%d",ach,tmp,mdl->idSelf);
			mdl->fpDiag = fopen(achDiag,"w");
			assert(mdl->fpDiag != NULL);
			}
		}
#ifdef TINY_PTHREAD_STACK
    save_mdl = &mdl;
    /*
     * Restart myself with a bigger stack.
     */
    pthread_create(mdl->pt,&attr,main,mdl->pmdl[0]);
    pthread_exit(0);
#endif
    return(nThreads);
    }


void mdlFinish(MDL mdl)
{
    MDL tmdl;
    int i;
    
    for (i=1;i<mdl->nThreads;++i) {
        pthread_join(mdl->pt[i],0);
        }
    for (i=0;i<mdl->nThreads;++i) {
        tmdl = mdl->pmdl[i];
        BasicDestroy(tmdl);
        /*
         ** Close Diagnostic file.
         */
        if (tmdl->bDiag) {
            fclose(tmdl->fpDiag);
            }
        }
    for (i=1;i<mdl->nThreads;++i) {
        free(mdl->pmdl[i]);
        }
    free(mdl->pmdl);
    free(mdl->pt);
    free(mdl);
    }


/*
 ** This function returns the number of threads in the set of 
 ** threads.
 */
int mdlThreads(MDL mdl)
{
    return(mdl->nThreads);
    }


/*
 ** This function returns this threads 'id' number within the specified
 ** MDL Context. Parent thread always has 'id' of 0, where as children
 ** have 'id's ranging from 1..(nThreads - 1).
 */
int mdlSelf(MDL mdl)
{
    return(mdl->idSelf);
    }


void mdlDiag(MDL mdl,char *psz)
{
    if (mdl->bDiag) {	
		fputs(psz,mdl->fpDiag);
		fflush(mdl->fpDiag);
		}
    }


void mdlAddService(MDL mdl,int sid,void *p1,
				   void (*fcnService)(void *,void *,int,void *,int *),
				   int nInBytes,int nOutBytes)
{
    MBX *pmbx;
    int i,nMaxServices;
    
    assert(sid > 0);
    if (sid >= mdl->nMaxServices) {
		/*
		 ** reallocate service buffer, adding space for 8 new services
		 ** including the one just defined.
		 */
		nMaxServices = sid + 9;
		mdl->psrv = realloc(mdl->psrv,nMaxServices*sizeof(SERVICE));
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
    pmbx = &mdl->mbxOwn;
    if (nInBytes > mdl->nMaxInBytes) {
		/*
		 ** Don't need to aquire lock here, because only 2 cases can
		 ** occur. 1) We are adding service outside of the handler,
		 ** and this okay (usual). 2) We are trying to allow a service 
		 ** to add a new service, this is also okay since while the 
		 ** service is running it already has the lock on the mailbox.
		 */
		pmbx->pszIn = realloc(pmbx->pszIn,nInBytes);
		mdl->nMaxInBytes = nInBytes;
		}
    if (nOutBytes > mdl->nMaxOutBytes) {
		pmbx->pszOut = realloc(pmbx->pszOut,nOutBytes);
		mdl->nMaxOutBytes = nOutBytes;
		}
    mdl->psrv[sid].p1 = p1;
    mdl->psrv[sid].nInBytes = nInBytes;
    mdl->psrv[sid].nOutBytes = nOutBytes;
    mdl->psrv[sid].fcnService = fcnService;
    }


void mdlReqService(MDL mdl,int id,int sid,void *vIn,int nInBytes)
{
    MBX *pmbx;
    char *pszIn = vIn;
    int i;
    
    pmbx = &mdl->pmdl[id]->mbxOwn;
	pthread_mutex_lock(&pmbx->mux);
	while (!pmbx->bRel) {
		pthread_cond_wait(&pmbx->sigRel,&pmbx->mux);
		}
    pmbx->bRel = 0;
    pmbx->sid = sid;
    pmbx->nBytes = nInBytes;
    for (i=0;i<nInBytes;++i) pmbx->pszIn[i] = pszIn[i];
	pmbx->bReq = 1;
	pthread_cond_signal(&pmbx->sigReq);
    pthread_mutex_unlock(&pmbx->mux);
    }


void mdlGetReply(MDL mdl,int id,void *vOut,int *pnOutBytes)
{
    MBX *pmbx;
    char *pszOut = vOut;
    int i;
    
    pmbx = &mdl->pmdl[id]->mbxOwn;
	pthread_mutex_lock(&pmbx->mux);
	while (!pmbx->bRpl) {
		pthread_cond_wait(&pmbx->sigRpl,&pmbx->mux);
		}
    pmbx->bRpl = 0;
    if (pnOutBytes) *pnOutBytes = pmbx->nBytes;
    if (pmbx->nBytes > 0 && pszOut != NULL) {
		for (i=0;i<pmbx->nBytes;++i) pszOut[i] = pmbx->pszOut[i];
		}
    /*
     ** If a STOP service was requested lock out any future requests 
     ** until the handler is restarted! In other words don't release 
     ** requests to this mailbox.
     ** Have to allow the restarted handler to aquire the lock, hence set
     ** bReq if this was a SRV_STOP service.
     */
    if (pmbx->sid != SRV_STOP) {
		pmbx->bRel = 1;
		pthread_cond_signal(&pmbx->sigRel);
		}
    pthread_mutex_unlock(&pmbx->mux);
    }


void mdlHandler(MDL mdl)
{
    MBX *pmbx;
    int nInBytes,sid;
    
    /*
     ** First open the floodgates!
     */
    pmbx = &mdl->mbxOwn;
    pthread_mutex_lock(&pmbx->mux);
    pmbx->bRel = 1;
	pthread_cond_signal(&pmbx->sigRel);
    pthread_mutex_unlock(&pmbx->mux);
    sid = 1;
    while (sid != SRV_STOP) {
		pthread_mutex_lock(&pmbx->mux);
		while (!pmbx->bReq) {
			pthread_cond_wait(&pmbx->sigReq,&pmbx->mux);
			}
		pmbx->bReq = 0;
        sid = pmbx->sid;
        assert(sid < mdl->nMaxServices);
        nInBytes = pmbx->nBytes;
        assert(nInBytes <= mdl->psrv[sid].nInBytes);
		pmbx->nBytes = 0;
        assert(mdl->psrv[sid].fcnService != NULL);
	(*mdl->psrv[sid].fcnService)(mdl->psrv[sid].p1,
				     pmbx->pszIn,nInBytes,
				     pmbx->pszOut,&pmbx->nBytes);
        assert(pmbx->nBytes <= mdl->psrv[sid].nOutBytes);
        pmbx->bRpl = 1;
		pthread_cond_signal(&pmbx->sigRpl);
        pthread_mutex_unlock(&pmbx->mux);
		}
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
int mdlSwap(MDL mdl,int id,size_t nBufBytes,void *vBuf,size_t nOutBytes,
	    size_t *pnSndBytes,size_t *pnRcvBytes)
{
    SWX *pout,*pin;
    size_t i,nInBytes,nOutBufBytes,nInMax,nOutMax;
    char *pszBuf = vBuf;
    char *pszIn,*pszOut;
    
    *pnRcvBytes = 0;
    *pnSndBytes = 0;
    /*
     **	Send number of rejects to target thread amount of free space
     */ 
    pout = &mdl->pmdl[id]->swxOwn;
    pin = &mdl->swxOwn;
	
	pthread_mutex_lock(&pout->mux);
    pout->nInBytes = nOutBytes;
    pout->nOutBufBytes = nBufBytes;
	pout->bRec = 1;
	pthread_cond_signal(&pout->sigRec);
    pthread_mutex_unlock(&pout->mux);
    /*
     ** Receive the number of target thread rejects and target free space
     */
	pthread_mutex_lock(&pin->mux);
    while (!pin->bRec) {
		pthread_cond_wait(&pin->sigRec,&pin->mux);
		}
    nInBytes = pin->nInBytes;
    nOutBufBytes = pin->nOutBufBytes;
    pin->bRec = 0;
	pthread_cond_signal(&pin->sigSnd);
    pthread_mutex_unlock(&pin->mux);
    pszIn = pszBuf;
    pszOut = &pszBuf[nBufBytes-nOutBytes];
    /*
     ** Start bilateral transfers. Note: One processor is GUARANTEED to 
     ** complete all its transfers.
     */
    while (nOutBytes && nInBytes) {
		/*
		 ** nOutMax is the maximum number of bytes allowed to be sent
		 ** nInMax is the number of bytes which will be received.
		 */
        nOutMax = (nOutBytes < MDL_TRANS_SIZE)?nOutBytes:MDL_TRANS_SIZE;
		nOutMax = (nOutMax < nOutBufBytes)?nOutMax:nOutBufBytes;
        nInMax = (nInBytes < MDL_TRANS_SIZE)?nInBytes:MDL_TRANS_SIZE;
		nInMax = (nInMax < nBufBytes)?nInMax:nBufBytes;
		/*
		 ** Transfer...
		 */
		pthread_mutex_lock(&pout->mux);
        while (pout->bRec) {
			pthread_cond_wait(&pout->sigSnd,&pout->mux);
            }
        for (i=0;i<nOutMax;++i) pout->pszBuf[i] = pszOut[i];
        pout->bRec = 1;
		pthread_cond_signal(&pout->sigRec);
        pthread_mutex_unlock(&pout->mux);
		pthread_mutex_lock(&pin->mux);
        while (!pin->bRec) {
			pthread_cond_wait(&pin->sigRec,&pin->mux);
            }
		for (i=0;i<nInMax;++i) pszIn[i] = pin->pszBuf[i];
        pin->bRec = 0;
		pthread_cond_signal(&pin->sigSnd);
        pthread_mutex_unlock(&pin->mux);
		/*
		 ** Adjust pointers and counts for next itteration.
		 */
		pszOut = &pszOut[nOutMax];
		nOutBytes -= nOutMax;
		nOutBufBytes -= nOutMax;
		*pnSndBytes += nOutMax;
		pszIn = &pszIn[nInMax];
		nInBytes -= nInMax;
		nBufBytes -= nInMax;
		*pnRcvBytes += nInMax;
		}
    /*
     ** At this stage we perform only unilateral transfers, and here we
     ** could exceed the opponent's storage capacity.
     */
    while (nOutBytes && nOutBufBytes) {
        nOutMax = (nOutBytes < MDL_TRANS_SIZE)?nOutBytes:MDL_TRANS_SIZE;
		nOutMax = (nOutMax < nOutBufBytes)?nOutMax:nOutBufBytes;
		pthread_mutex_lock(&pout->mux);
        while (pout->bRec) {
			pthread_cond_wait(&pout->sigSnd,&pout->mux);
            }
        for (i=0;i<nOutMax;++i) pout->pszBuf[i] = pszOut[i];
        pout->bRec = 1;
		pthread_cond_signal(&pout->sigRec);
        pthread_mutex_unlock(&pout->mux);
		/*
		 ** Adjust pointers and counts.
		 */
		pszOut = &pszOut[nOutMax];
		nOutBytes -= nOutMax;
		nOutBufBytes -= nOutMax;
		*pnSndBytes += nOutMax;
		}
    while (nInBytes && nBufBytes) {
        nInMax = (nInBytes < MDL_TRANS_SIZE)?nInBytes:MDL_TRANS_SIZE;
		nInMax = (nInMax < nBufBytes)?nInMax:nBufBytes;
		pthread_mutex_lock(&pin->mux);
        while (!pin->bRec) {
			pthread_cond_wait(&pin->sigRec,&pin->mux);
            }
		for (i=0;i<nInMax;++i) pszIn[i] = pin->pszBuf[i];
        pin->bRec = 0;
		pthread_cond_signal(&pin->sigSnd);
        pthread_mutex_unlock(&pin->mux);
        pszIn = &pszIn[nInMax];
		nInBytes -= nInMax;
		nBufBytes -= nInMax;
		*pnRcvBytes += nInMax;
		}
    if (nOutBytes) return(0);
    else if (nInBytes) return(0);
    else return(1);
    }

/*
 ** START OF CACHE CODE
 */
#define MDL_MID_CACHEOUT	4
#define MDL_MID_CACHEFLSH	5

#define MDL_RANDMOD		1771875
#define MDL_RAND(mdl) (mdl->uRand = (mdl->uRand*2416+374441)%MDL_RANDMOD)
#define MDL_CHECK_MASK  	0x7f
#define BILLION				1000000000

#define MDL_ADVANCE_RING(X) (X)++; if((X) >= MDL_MBX_RING_SZ) (X) = 0

int mdlCacheReceive(MDL mdl)
{
	CACHE *c;
	MBX *pmbx; /* Cache mailbox */
	CAHEAD *ph;
	char *pszRcv;
	char *t;
	int n,i,nBytes;

	pthread_mutex_lock(&mdl->muxRing);
	assert(mdl->iRingHd != mdl->iRingTl);
	pmbx = &mdl->mbxCache[mdl->iRingHd]; /* Cache mailbox */
	MDL_ADVANCE_RING(mdl->iRingHd);
	pthread_mutex_unlock(&mdl->muxRing);

	ph = (CAHEAD *)pmbx->pszIn;
	pszRcv = &pmbx->pszIn[sizeof(CAHEAD)];
	pthread_mutex_lock(&pmbx->mux);

	assert(ph->iSeq == mdl->iRecSeq[ph->id]);
	
	mdl->iRecSeq[ph->id]++;

	c = &mdl->cache[ph->cid];
	switch (ph->mid) {
	case MDL_MID_CACHEOUT:
		++c->nCheckOut;
		pthread_mutex_unlock(&pmbx->mux);
		return(0);
	case MDL_MID_CACHEFLSH:
		assert(c->iType == MDL_COCACHE);
		/*
		 ** Unpack the data into the 'sentinel-line' cache data.
		 */
		nBytes = c->iLineSize;
		for (i=0;i<nBytes;++i) c->pLine[i] = pszRcv[i];
		i = ph->iLine*MDL_CACHELINE_ELTS;
		/*
		 * Data is out; unlock it
		 */
		pthread_mutex_unlock(&pmbx->mux);
		t = &c->pData[i*((size_t)c->iDataSize)];
		/*
		 ** Make sure we don't combine beyond the number of data elements!
		 */
		n = i + MDL_CACHELINE_ELTS;
		if (n > c->nData) n = c->nData;
		n -= i;
		n *= c->iDataSize;
		for (i=0;i<n;i+=c->iDataSize) {
			(*c->combine)(&t[i],&c->pLine[i]);
			}
		return(0);
	default:
	    assert(0);
	    }
	}

void AdjustDataSize(MDL mdl)
{
	int i,iMaxDataSize;

	/*
	 ** Change buffer size?
	 */
	iMaxDataSize = 0;
	for (i=0;i<mdl->nMaxCacheIds;++i) {
		if (mdl->cache[i].iType == MDL_NOCACHE) continue;
		if (mdl->cache[i].iDataSize > iMaxDataSize) {
			iMaxDataSize = mdl->cache[i].iDataSize;
			}
		}
	if (iMaxDataSize != mdl->iMaxDataSize) {
		/*
		 ** Create new buffer with realloc?
		 ** Be very careful when reallocing buffers in other libraries
		 ** (not PVM) to be sure that the buffers are not in use!
		 ** A pending non-blocking receive on a buffer which is realloced
		 ** here will cause problems, make sure to take this into account!
		 ** This is certainly true in using the MPL library.
		 */
		mdl->iMaxDataSize = iMaxDataSize;
		mdl->iCaBufSize = sizeof(CAHEAD) + 
			iMaxDataSize*(1 << MDL_CACHELINE_BITS);
		for(i = 0; i < MDL_MBX_RING_SZ; i++) {
		    mdl->mbxCache[i].pszIn = realloc(mdl->mbxCache[i].pszIn,
						  mdl->iCaBufSize);
		    assert(mdl->mbxCache[i].pszIn != NULL);
		    }
		}
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
void *mdlMalloc(MDL mdl,size_t iSize)
{
    return(malloc(iSize));
	}

void mdlFree(MDL mdl,void *p)
{
	free(p);
	}

/*
 ** Common initialization for all types of caches.
 */
CACHE *CacheInitialize(MDL mdl,int cid,void *pData,int iDataSize,int nData)
{
	CACHE *c;
	int i,nMaxCacheIds;

	/*
	 ** Allocate more cache spaces if required!
	 */
	assert(cid >= 0);
	if (cid >= mdl->nMaxCacheIds) {
		/*
		 ** reallocate cache spaces, adding space for 2 new cache spaces
		 ** including the one just defined.
		 */
		nMaxCacheIds = cid + 3;
		mdl->cache = realloc(mdl->cache,nMaxCacheIds*sizeof(CACHE));
		assert(mdl->cache != NULL);
		/*
		 ** Initialize the new cache slots.
		 */
		for (i=mdl->nMaxCacheIds;i<nMaxCacheIds;++i) {
			mdl->cache[i].iType = MDL_NOCACHE;
			}
		mdl->nMaxCacheIds = nMaxCacheIds;
		}
	c = &mdl->cache[cid];
	assert(c->iType == MDL_NOCACHE);
	c->pData = pData;
	c->iDataSize = iDataSize;
	c->nData = nData;
	c->pDataMax = nData*((size_t)iDataSize);

	c->pTrans = NULL;
	c->pTag = NULL;
	c->pbKey = NULL;
	c->pLine = NULL;
	
	return(c);
	}

/*
 ** Initialize a caching space.
 */
void mdlROcache(MDL mdl,int cid,void *pData,int iDataSize,int nData)
{
	CACHE *c;

	c = CacheInitialize(mdl,cid,pData,iDataSize,nData);
	c->iType = MDL_ROCACHE;
	/*
	 ** For an ROcache these two functions are not needed.
	 */
	c->init = NULL;
	c->combine = NULL;
	/*
	 ** THIS IS A SYNCHRONIZE!!!
	 */
	mdlBarrier(mdl);
	}

/*
 ** Initialize a combiner caching space.
 */
void mdlCOcache(MDL mdl,int cid,void *pData,int iDataSize,int nData,
				void (*init)(void *),void (*combine)(void *,void *))
{
	CACHE *c;
	int i;

	c = CacheInitialize(mdl,cid,pData,iDataSize,nData);

	c->iType = MDL_COCACHE;
	c->init = init;
	c->combine = combine;
	c->iLineSize = MDL_CACHELINE_ELTS*c->iDataSize;
	c->iKeyShift = 0;
	while((1 << c->iKeyShift) < mdl->nThreads) ++c->iKeyShift;
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
	c->nLines = (MDL_CACHE_SIZE/c->iDataSize) >> MDL_CACHELINE_BITS;
	assert(c->nLines < MDL_RANDMOD);
	c->nTrans = 1;
	while(c->nTrans < c->nLines) c->nTrans *= 2;
	c->nTrans *= 2;
	c->iTransMask = c->nTrans-1;
	/*
	 **	Set up the translation table.
	 */
	c->pTrans = malloc(c->nTrans*sizeof(int));	
	assert(c->pTrans != NULL);
	for (i=0;i<c->nTrans;++i) c->pTrans[i] = 0;
	/*
	 ** Set up the tags. Note pTag[0] is a Sentinel!
	 */
	c->pTag = malloc(c->nLines*sizeof(CTAG));
	assert(c->pTag != NULL);
	for (i=0;i<c->nLines;++i) {
		c->pTag[i].iKey = -1;	/* invalid */	
		c->pTag[i].nLock = 0;
		c->pTag[i].nLast = 0;	/* !!! */
		c->pTag[i].iLink = 0;
		}
	c->pTag[0].nLock = 1;		/* always locked */
	c->pTag[0].nLast = INT_MAX;  	/* always Most Recently Used */
	c->nAccess = 0;
	c->nAccHigh = 0;
	c->nMiss = 0;				/* !!!, not NB */
	c->nColl = 0;				/* !!!, not NB */
	c->nMin = 0;				/* !!!, not NB */	
	c->nKeyMax = 500;				/* !!!, not NB */
	c->pbKey = malloc(c->nKeyMax);			/* !!!, not NB */
	assert(c->pbKey != NULL);			/* !!!, not NB */
	for (i=0;i<c->nKeyMax;++i) c->pbKey[i] = 0;	/* !!!, not NB */
	/*
	 ** Allocate cache data lines.
	 */
	c->pLine = malloc(c->nLines*c->iLineSize);
	assert(c->pLine != NULL);
	c->nCheckOut = 0;
	for(i = 0; i < mdl->nThreads; i++) {
	    mdl->iRecSeq[i] = 0;
	    mdl->iSndSeq[i] = 0;
	    }
	/*
	 ** THIS IS A SYNCHRONIZE!!!
	 */
	AdjustDataSize(mdl);
	mdlBarrier(mdl);
	}

void mdlCacheRequest(MDL mdl, int id, int cid, int mid, char *pszData,
		     int iLine, int iLineSize)
{
    pthread_mutex_t *pmux_ring = &mdl->pmdl[id]->muxRing;
    MBX *pmbx;
    CAHEAD *caFlsh;
    char *pszFlsh;
    int iRingTl;
    int iOldRingTl;

    assert(id != mdl->idSelf);
    /*
     * Grab my place in the ring.
     */
    pthread_mutex_lock(pmux_ring);
    iRingTl = mdl->pmdl[id]->iRingTl;
    MDL_ADVANCE_RING(iRingTl);
    while(iRingTl == mdl->pmdl[id]->iRingHd) {
	pthread_mutex_unlock(pmux_ring);
	mdlCacheCheck(mdl);
	pthread_mutex_lock(pmux_ring);
	iRingTl = mdl->pmdl[id]->iRingTl;
	MDL_ADVANCE_RING(iRingTl);
	}
    
    iOldRingTl = mdl->pmdl[id]->iRingTl;
    mdl->pmdl[id]->iRingTl = iRingTl;
    
    pmbx = &mdl->pmdl[id]->mbxCache[iOldRingTl];

    caFlsh = (CAHEAD *)pmbx->pszIn;
    pszFlsh = &pmbx->pszIn[sizeof(CAHEAD)];
    pthread_mutex_lock(&pmbx->mux);

    pthread_mutex_unlock(pmux_ring);

    caFlsh->cid = cid;
    caFlsh->mid = mid;
    caFlsh->id = mdl->idSelf;
    caFlsh->iSeq = mdl->iSndSeq[id];
    mdl->iSndSeq[id]++;
    if(mid == MDL_MID_CACHEFLSH) {
	int j;
	
	caFlsh->iLine = iLine;
	for(j = 0; j < iLineSize; ++j)
	    pszFlsh[j] = pszData[j];
	}
    pthread_mutex_unlock(&pmbx->mux);
    }

void mdlFinishCache(MDL mdl,int cid)
{
	CACHE *c = &mdl->cache[cid];
	int i,id;
	int iKey;

	/*
	 ** THIS IS A SYNCHRONIZE!!!
	 */
	if(c->iType == MDL_COCACHE) {
		/*
		 ** Must flush all valid data elements.
		 */
		for (i=1;i<c->nLines;++i) {
			iKey = c->pTag[i].iKey;
			if (iKey >= 0) {
				/*
				 ** Flush element since it is valid!
				 */
				int iLine = iKey >> c->iInvKeyShift;
			    
				id = iKey & c->iIdMask;
				mdlCacheRequest(mdl, id, cid,
						MDL_MID_CACHEFLSH, 
						&c->pLine[i*c->iLineSize],
						iLine, c->iLineSize);
						
				mdlCacheCheck(mdl); /* service incoming */
				}
			}
		if(mdl->idSelf == 0) {
		    ++c->nCheckOut;
		    while(c->nCheckOut < mdl->nThreads) {
			mdlCacheCheck(mdl);
			}
		    }
		else {
		    mdlCacheRequest(mdl, 0, cid, MDL_MID_CACHEOUT,
				    NULL, 0, 0);
		    }
		if(mdl->idSelf == 0) {
		    for(id = 1; id < mdl->nThreads; id++) {
			mdlCacheRequest(mdl, id, cid, MDL_MID_CACHEOUT,
				    NULL, 0, 0);
			}
		    }
		else {
		    c->nCheckOut = 0;
		    while (c->nCheckOut == 0) {
			mdlCacheCheck(mdl);
			}	
		    }	
	    }
	else {
	    mdlBarrier(mdl);
	    }
	
	/*
	 ** Free up storage and finish.
	 */
	if(c->iType == MDL_COCACHE) {
	    free(c->pTrans);
	    free(c->pTag);
	    free(c->pbKey);
	    free(c->pLine);
	    }
	c->iType = MDL_NOCACHE;
	AdjustDataSize(mdl);
	}

void *doMiss(MDL mdl, int cid, int iIndex, int id, int iKey, int lock);

void *mdlAquire(MDL mdl,int cid,int iIndex,int id)
{
	CACHE *c = &mdl->cache[cid];
	int iKey, i;
	char *pLine;		/* matched cache line */
	int iElt;		/* Element in line */

	if(c->iType == MDL_ROCACHE || id == mdl->idSelf) {
	    CACHE *cc = &mdl->pmdl[id]->cache[cid];
	    return(&cc->pData[iIndex*((size_t)c->iDataSize)]);
	    }

	++c->nAccess;
	if (!(c->nAccess & MDL_CHECK_MASK)) {
	    mdlCacheCheck(mdl);
	    }
	
	/*
	 ** Determine memory block key value and cache line.
	 */
	iKey = ((iIndex&MDL_INDEX_MASK) << c->iKeyShift)| id;

	i = c->pTrans[iKey & c->iTransMask];
	/*
	 ** Check for a match!
	 */
	if (c->pTag[i].iKey == iKey) {
		++c->pTag[i].nLock;
		pLine = &c->pLine[i*c->iLineSize];
		iElt = iIndex & MDL_CACHE_MASK;
		return(&pLine[iElt*c->iDataSize]);
		}
	i = c->pTag[i].iLink;
	/*
	 ** Collision chain search.
	 */
	while (i) {
		++c->nColl;
		if (c->pTag[i].iKey == iKey) {
			++c->pTag[i].nLock;
			pLine = &c->pLine[i*c->iLineSize];
			iElt = iIndex & MDL_CACHE_MASK;
			return(&pLine[iElt*c->iDataSize]);
			}
		i = c->pTag[i].iLink;
		}
	return(doMiss(mdl, cid, iIndex, id, iKey, 1));
    }

void *doMiss(MDL mdl, int cid, int iIndex, int id, int iKey, int lock)
{
	CACHE *c = &mdl->cache[cid];
	CACHE *cc;
	char *pLine;
	int iElt,iLine,i,iKeyVic,nKeyNew;
	int idVic;
	int iVictim,*pi;
	int iLineSize;

	/*
	 ** Cache Miss.
	 */
	++c->nMiss;
	iLine = iIndex >> MDL_CACHELINE_BITS;
	/*
	 **	Victim Search!
	 ** Note: if more than 1771875 cache lines are present this random
	 ** number generation may have to be changed, although none of the
	 ** code will break in this case. The only problem may be non-optimal
	 ** cache line replacement. Maybe give a warning at initialization?
	 */
	iVictim = MDL_RAND(mdl)%c->nLines;
	iElt = iIndex & MDL_CACHE_MASK;
	for (i=0;i<c->nLines;++i) {
	    if (!c->pTag[iVictim].nLock) {
		/*
		 ** Found victim.
		 */
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
			    int iLine = iKeyVic >> c->iInvKeyShift;

			    idVic = iKeyVic&c->iIdMask;

#if 0
				fprintf(stderr, "%d %d %d\n",
					mdl->idSelf, idVic, iLine);
#endif
			    
			    mdlCacheRequest(mdl, idVic, cid,
					    MDL_MID_CACHEFLSH, pLine,
					    iLine, c->iLineSize);
			    }
			/*
			 ** If valid iLine then "unlink" it from the cache.
			 */
			pi = &c->pTrans[iKeyVic & c->iTransMask];
			while (*pi != iVictim) pi = &c->pTag[*pi].iLink;
			*pi = c->pTag[iVictim].iLink;
			}
		c->pTag[iVictim].iKey=iKey;
		if(lock)
		    c->pTag[iVictim].nLock = 1;
		/*
		 **	Add the modified victim tag back into the cache.
		 ** Note: the new element is placed at the head of the chain.
		 */
		pi = &c->pTrans[iKey & c->iTransMask];
		c->pTag[iVictim].iLink = *pi;
		*pi = iVictim;
		goto Await;
		}
	    if (++iVictim == c->nLines) iVictim = 0;
	    }
	/*
	 ** Cache Failure!
	 */
	mdlprintf(mdl, "MDL CACHE FAILURE: cid == %d, no unlocked lines!\n",cid);
	assert(0);
 Await:
	/*
	 ** Figure out whether this is a "new" miss.
	 ** This is for statistics only!
	 */
	if (iKey >= c->nKeyMax) {			/* !!! */
		nKeyNew = iKey+500;
		c->pbKey = realloc(c->pbKey,nKeyNew);
		assert(c->pbKey != NULL);
		for (i=c->nKeyMax;i<nKeyNew;++i) c->pbKey[i] = 0;
		c->nKeyMax = nKeyNew;
		}
	if (!c->pbKey[iKey]) {
		c->pbKey[iKey] = 1;
		++c->nMin;
		}
	/*
	 ** At this point 'pLine' is the recipient cache line for the 
	 ** data requested from processor 'id'.
	 */	
	cc = &mdl->pmdl[id]->cache[cid];

	if ((iLine+1)*c->iLineSize > cc->pDataMax) {
	    iLineSize = cc->pDataMax - iLine*c->iLineSize;
	}
	else {
	    iLineSize = c->iLineSize;
	} 

	for(i = 0; i < iLineSize; i++)
	    pLine[i] = cc->pData[iLine*c->iLineSize + i];
	
	if (c->iType == MDL_COCACHE) {
		/*
		 ** Call the initializer function for all elements in 
		 ** the cache line.
		 */
		for (i=0;i<iLineSize;i+=c->iDataSize) {
			(*c->init)(&pLine[i]);
			}
		}

	return(&pLine[iElt*c->iDataSize]);
	}

void mdlRelease(MDL mdl,int cid,void *p)
{
	CACHE *c = &mdl->cache[cid];
	int iLine,iData;
	
	if(c->iType == MDL_ROCACHE)
	    return;
	
	iLine = ((char *)p - c->pLine) / c->iLineSize;
	/*
	 ** Check if the pointer fell in a cache line, otherwise it
	 ** must have been a local pointer.
	 */
	if (iLine > 0 && iLine < c->nLines) {
		--c->pTag[iLine].nLock;
		assert(c->pTag[iLine].nLock >= 0);
		}
	else {
		iData = ((char *)p - c->pData) / c->iDataSize;
		assert(iData >= 0 && iData < c->nData);
		}
    }

void mdlCacheCheck(MDL mdl)
{
    while (1) {
	pthread_mutex_lock(&mdl->muxRing);
	if(mdl->iRingHd == mdl->iRingTl) {
	    pthread_mutex_unlock(&mdl->muxRing);
	    break;
	    }
	pthread_mutex_unlock(&mdl->muxRing);
	mdlCacheReceive(mdl);
	}
    }

double mdlNumAccess(MDL mdl,int cid)
{
	CACHE *c = &mdl->cache[cid];

	return(c->nAccHigh*1e9 + c->nAccess);
	}


double mdlMissRatio(MDL mdl,int cid)
{
	CACHE *c = &mdl->cache[cid];
	double dAccess = c->nAccHigh*1e9 + c->nAccess;
	
	if (dAccess > 0.0) return(c->nMiss/dAccess);
	else return(0.0);
	}


double mdlCollRatio(MDL mdl,int cid)
{
	CACHE *c = &mdl->cache[cid];
	double dAccess = c->nAccHigh*1e9 + c->nAccess;

	if (dAccess > 0.0) return(c->nColl/dAccess);
	else return(0.0);
	}


double mdlMinRatio(MDL mdl,int cid)
{
	CACHE *c = &mdl->cache[cid];
	double dAccess = c->nAccHigh*1e9 + c->nAccess;

	if (dAccess > 0.0) return(c->nMin/dAccess);
	else return(0.0);
	}
