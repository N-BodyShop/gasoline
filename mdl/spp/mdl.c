/*
 ** A very low level Machine model. For homogeneous systems only!
 ** This is the KSR-pthreads mdl module.
 **
 ** WARNING: Do not use pthread_cond_wait() anywhere!
 */
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>
#include <cps.h>
#include <sys/times.h>
#include <unistd.h>
#include <assert.h>
#include <stdarg.h>
#include "mdl.h"


#define MDL_NOCACHE			0
#define MDL_ROCACHE			1
#define MDL_COCACHE			2


#define MDL_DEFAULT_BYTES		4096
#define MDL_DEFAULT_SERVICES	50
#define MDL_DEFAULT_CACHEIDS	5

#define MDL_TRANS_SIZE			50000 


void srvNull(void *p1,void *vin,int nIn,void *vout,int *pnOut)
{
    return;
    }


double mdlCpuTimer(MDL mdl)
{
/*    return(user_seconds()); */
        double t;
	struct tms timeinfo;

/*	t=(double) sysconf(_SC_CLK_TCK);
	times(&timeinfo);
	t=timeinfo.tms_utime/t; */
	t=.001;
	return t;
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
    cps_mutex_alloc(&mdl->mbxOwn.mux);
    /*
     ** Initialize the swapbox, this buffer has fixed size but must
     ** still be malloc'd to ensure that no allignment problems can
     ** occur.
     */
    mdl->swxOwn.pszBuf = malloc(MDL_TRANS_SIZE);
    assert(mdl->swxOwn.pszBuf != NULL);
    mdl->swxOwn.bSet = 0;
    cps_mutex_alloc(&mdl->swxOwn.mux);
    }


void BasicDestroy(MDL mdl)
{
    cps_mutex_free(&mdl->mbxOwn.mux);
    cps_mutex_free(&mdl->swxOwn.mux);
    free(mdl->swxOwn.pszBuf);
    free(mdl->psrv);
    free(mdl->mbxOwn.pszIn);
    free(mdl->mbxOwn.pszOut);
    free(mdl->cache);
    }


int mdlInitialize(MDL *pmdl,char **argv,void (*fcnChild)(MDL))
{
    MDL mdl,tmdl;
    int i,nThreads,bThreads,bDiag;
    char *p,ach[256],achDiag[256];
    const int node=CPS_ANY_NODE;
    
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
		nThreads = cps_complex_cpus();
		}
    mdl = malloc(sizeof(struct mdlContext));
    assert(mdl != NULL);
    mdl->pmdl = malloc(nThreads*sizeof(MDL));
    assert(mdl->pmdl != NULL);
    mdl->pmdl[0] = mdl;			/* that's me! */
    mdl->pt = (ktid_t *)malloc(nThreads*sizeof(ktid_t));
    assert(mdl->pt != NULL);
    /*
     ** Initialize caching barrier.
     */
    cps_barrier_alloc(&mdl->bar);
    *pmdl = mdl;
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
				sprintf(achDiag,"%s.%d",ach,tmdl->idSelf);
				tmdl->fpDiag = fopen(achDiag,"w");
				assert(tmdl->fpDiag != NULL);
				}
			}
		for (i=1;i<nThreads;++i) {
			/*
			 ** Start all the children.
			 */
		    mdl->pt[i]=cps_thread_create(&node,fcnChild,mdl->pmdl[i]);
		    /*pthread_create(&mdl->pt[i],pthread_attr_default,fcnChild,
						   mdl->pmdl[i]);*/
			}
		return(nThreads);
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
			sprintf(achDiag,"%s.%d",ach,mdl->idSelf);
			mdl->fpDiag = fopen(achDiag,"w");
			assert(mdl->fpDiag != NULL);
			}
		return(nThreads);
		}
    }


void mdlFinish(MDL mdl)
{
    MDL tmdl;
    int i;
    int flag=0,nthreads;
    
    nthreads=cps_thread_wait(&flag);
    fprintf(stderr,"%d threads still active\n",nthreads);
    flag=1;
    cps_thread_wait(&flag);
	
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
    cps_barrier_free(&mdl->bar);
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
int mdlSwap(MDL mdl,int id,int nBufBytes,void *vBuf,int nOutBytes,
	    int *pnSndBytes,int *pnRcvBytes)
{
    SWX *pout,*pin;
    int i,nInBytes,nOutBufBytes,nInMax,nOutMax;
    char *pszBuf = vBuf;
    char *pszIn,*pszOut;
    
    *pnRcvBytes = 0;
    *pnSndBytes = 0;
    /*
     **	Send number of rejects to target thread amount of free space
     */ 
    pout = &mdl->pmdl[id]->swxOwn;
    pin = &mdl->swxOwn;
    while (1) {
	cps_mutex_lock(&pout->mux);
	if (!pout->bSet) break;
        cps_mutex_unlock(&pout->mux);
	}
    pout->bSet = 1;
    pout->nInBytes = nOutBytes;
    pout->nOutBufBytes = nBufBytes;
    cps_mutex_unlock(&pout->mux);
    /*
     ** Receive the number of target thread rejects and target free space
     */
    while (1) {
	cps_mutex_lock(&pin->mux);
	if (pin->bSet) break;
        cps_mutex_unlock(&pin->mux);
	}
    pin->bSet = 0;
    nInBytes = pin->nInBytes;
    nOutBufBytes = pin->nOutBufBytes;
    cps_mutex_unlock(&pin->mux);
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
        while (1) {
            cps_mutex_lock(&pout->mux);
            if (!pout->bSet) break;
            cps_mutex_unlock(&pout->mux);
            }
        pout->bSet = 1;
        for (i=0;i<nOutMax;++i) pout->pszBuf[i] = pszOut[i];
        cps_mutex_unlock(&pout->mux);
        while (1) {
            cps_mutex_lock(&pin->mux);
            if (pin->bSet) break;
            cps_mutex_unlock(&pin->mux);
            }
        pin->bSet = 0;
	for (i=0;i<nInMax;++i) pszIn[i] = pin->pszBuf[i];
        cps_mutex_unlock(&pin->mux);
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
        while (1) {
            cps_mutex_lock(&pout->mux);
            if (!pout->bSet) break;
            cps_mutex_unlock(&pout->mux);
            }
        pout->bSet = 1;
        for (i=0;i<nOutMax;++i) pout->pszBuf[i] = pszOut[i];
        cps_mutex_unlock(&pout->mux);
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
        while (1) {
            cps_mutex_lock(&pin->mux);
            if (pin->bSet) break;
            cps_mutex_unlock(&pin->mux);
            }
        pin->bSet = 0;
		for (i=0;i<nInMax;++i) pszIn[i] = pin->pszBuf[i];
        cps_mutex_unlock(&pin->mux);
        pszIn = &pszIn[nInMax];
		nInBytes -= nInMax;
		nBufBytes -= nInMax;
		*pnRcvBytes += nInMax;
		}
    if (nOutBytes) return(0);
    else if (nInBytes) return(0);
    else return(1);
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
    while (1) {
		cps_mutex_lock(&pmbx->mux);
		if (pmbx->bRel) break;
        cps_mutex_unlock(&pmbx->mux);
		}
    pmbx->bRel = 0;
    pmbx->bReq = 1;
    pmbx->sid = sid;
    pmbx->nBytes = nInBytes;
    for (i=0;i<nInBytes;++i) pmbx->pszIn[i] = pszIn[i];
    cps_mutex_unlock(&pmbx->mux);
    }


void mdlGetReply(MDL mdl,int id,void *vOut,int *pnOutBytes)
{
    MBX *pmbx;
    char *pszOut = vOut;
    int i;
    
    pmbx = &mdl->pmdl[id]->mbxOwn;
    while (1) {
		cps_mutex_lock(&pmbx->mux);
		if (pmbx->bRpl) break;
        cps_mutex_unlock(&pmbx->mux);
		}
    pmbx->bRpl = 0;
    /*
     ** If a STOP service was requested lock out any future requests 
     ** until the handler is restarted! In other words don't release 
     ** requests to this mailbox.
     ** Have to allow the restarted handler to aquire the lock, hence set
     ** bReq if this was a SRV_STOP service.
     */
    if (pmbx->sid != SRV_STOP) pmbx->bRel = 1;
    if (pnOutBytes) *pnOutBytes = pmbx->nBytes;
    if (pmbx->nBytes > 0 && pszOut != NULL) {
		for (i=0;i<pmbx->nBytes;++i) pszOut[i] = pmbx->pszOut[i];
		}
    cps_mutex_unlock(&pmbx->mux);
    }


void mdlHandler(MDL mdl)
{
    MBX *pmbx;
    int i,nInBytes,sid;
    
    /*
     ** First open the floodgates!
     */
    pmbx = &mdl->mbxOwn;
    cps_mutex_lock(&pmbx->mux);
    pmbx->bRel = 1;
    cps_mutex_unlock(&pmbx->mux);
    sid = 1;
    while (sid != SRV_STOP) {
		while (1) {
			cps_mutex_lock(&pmbx->mux);
			if (pmbx->bReq) break;
            cps_mutex_unlock(&pmbx->mux);
			}
        pmbx->bReq = 0;
        pmbx->bRpl = 1;
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
        cps_mutex_unlock(&pmbx->mux);
		}
    }


#define MDL_ROCACHE			1
#define MDL_COCACHE			2

#define MDL_RANDMOD		1771875
#define MDL_RAND(mdl) (mdl->uRand = (mdl->uRand*2416+374441)%MDL_RANDMOD)


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
		 */
		mdl->iMaxDataSize = iMaxDataSize;
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
void *mdlMalloc(MDL mdl,int iSize)
{	
	return(malloc(iSize));
	}


void mdlFree(MDL mdl,void *p)
{
	free(p);
	}


/*
 ** Initialize a caching space.
 */
void mdlROcache(MDL mdl,int cid,void *pData,int iDataSize,int nData)
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
    c->iType = MDL_ROCACHE;
    c->pData = pData;
    c->iDataSize = iDataSize;
    c->nData = nData;
    /*
     ** Must synchronize.
     */
/*    pthread_barrier_checkout(&mdl->pmdl[0]->bar,mdl->idSelf);
    pthread_barrier_checkin(&mdl->pmdl[0]->bar,mdl->idSelf); */
    cps_barrier(&mdl->pmdl[0]->bar, &mdl->nThreads);
    AdjustDataSize(mdl);
    }


#define MDL_CACHE_SIZE		500000
#define MDL_CACHELINE_BITS	4
#define MDL_CACHE_MASK		((1<<MDL_CACHELINE_BITS)-1)


void mdlCOcache(MDL mdl,int cid,void *pData,int iDataSize,int nData,
				void (*init)(void *),void (*combine)(void *,void *))
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
    c->iType = MDL_COCACHE;
    c->pData = pData;
    c->iDataSize = iDataSize;
    c->nData = nData;
    c->nLineElts = (1 << MDL_CACHELINE_BITS); 
    c->iLineSize = c->nLineElts*c->iDataSize;
    c->nLines = (MDL_CACHE_SIZE/c->iDataSize) >> MDL_CACHELINE_BITS;
    assert(c->nLines < MDL_RANDMOD);
    c->nTrans = 1;
    while(c->nTrans < c->nLines) c->nTrans *= 2;
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
    c->pTag[0].nLock = 1;		/* always locked */
    for (i=1;i<c->nLines;++i) {
		c->pTag[i].id = -1;		/* invalid */	
		c->pTag[i].iLine = -1;	/* invalid */	
		c->pTag[i].nLock = 0;
		c->pTag[i].iLink = 0;
		}
    /*
     ** Allocate cache data lines.
     */
    c->pLine = malloc(c->nLines*c->iLineSize);
    assert(c->pLine != NULL);
    c->init = init;
    c->combine = combine;
    /*
     ** Must allocate an array of mutex for the actual data.
     */
    c->nMux = (nData >> MDL_CACHELINE_BITS) + 1;
    c->pMux = malloc(c->nMux*sizeof(cps_mutex_t));
    assert(c->pMux != NULL);
    for (i=0;i<c->nMux;++i) {
		cps_mutex_alloc(&c->pMux[i]);
		}
    /*
     ** Must synchronize.
     */
/*    pthread_barrier_checkout(&mdl->pmdl[0]->bar,mdl->idSelf);
    pthread_barrier_checkin(&mdl->pmdl[0]->bar,mdl->idSelf);*/
    cps_barrier(&mdl->pmdl[0]->bar,&mdl->nThreads);
    AdjustDataSize(mdl);
    }


void mdlFinishCache(MDL mdl,int cid)
{
    CACHE *c = &mdl->cache[cid];
    CACHE *cc;
    char *srce,*dest;
    int i,j,n,iLine,id;
    
    if (c->iType == MDL_ROCACHE) {
		/*
		 ** Must synchronize.
		 */
/*		pthread_barrier_checkout(&mdl->pmdl[0]->bar,mdl->idSelf);
		pthread_barrier_checkin(&mdl->pmdl[0]->bar,mdl->idSelf); */
	cps_barrier(&mdl->pmdl[0]->bar,&mdl->nThreads);
        c->iType = MDL_NOCACHE;
        AdjustDataSize(mdl);
		return;
		}
    /*
     ** Must flush all valid data elements.
     */
    for (i=1;i<c->nLines;++i) {
		iLine = c->pTag[i].iLine;
		id = c->pTag[i].id;
		if (iLine >= 0) {
			/*
			 ** Flush element since it is valid!
			 */
			cc = &mdl->pmdl[id]->cache[cid];	
			srce = &c->pLine[i*c->iLineSize];
			dest = &cc->pData[iLine*c->iLineSize];
			/*
			 ** Make sure we don't combine beyond the number of data elements!
			 */
			j = iLine*c->nLineElts;
			n = j + c->nLineElts;
			if (n > cc->nData) n = cc->nData;
			n -= j;
			n *= c->iDataSize;
			/*
			 ** Lock data line and combine!
			 */
			cps_mutex_lock(&cc->pMux[iLine]);
			for (j=0;j<n;j+=c->iDataSize) {
				(*c->combine)(&dest[j],&srce[j]);
				}
			cps_mutex_unlock(&cc->pMux[iLine]);
			}
		}
    /*
     ** Must synchronize.
     */
/*    pthread_barrier_checkout(&mdl->pmdl[0]->bar,mdl->idSelf);
    pthread_barrier_checkin(&mdl->pmdl[0]->bar,mdl->idSelf); */
    cps_barrier(&mdl->pmdl[0]->bar,&mdl->nThreads);
    /*
     ** Destroy array of mutex.
     */
    for (i=0;i<c->nMux;++i) {
		cps_mutex_free(&c->pMux[i]);
		}
    /*
     ** Free up storage and finish.
     */
    free(c->pMux);
    free(c->pTrans);
    free(c->pTag);
    free(c->pLine);
    c->iType = MDL_NOCACHE;
    AdjustDataSize(mdl);
    }


void *mdlAquire(MDL mdl,int cid,int iIndex,int id)
{
    CACHE *cc = &mdl->pmdl[id]->cache[cid];
    CACHE *c,*cflsh;
    char *pLine,*srce,*dest;
    int iElt,iLine,i,j,n;
    int iVictim,iLineVic,idVic,*pi;
    char ach[80];
    
    if (cc->iType == MDL_ROCACHE) {
		return(&cc->pData[iIndex*cc->iDataSize]);
		}
    c = &mdl->cache[cid];
    /*
     ** Note that for KSR even local requests are cached, because
     ** even they have to aquire a lock to access the memory.
     ** This excludes local memory "cheats" where the COcache is
     ** being used!
     ** Determine memory block key value and cache line.
     */
    iElt = iIndex & MDL_CACHE_MASK;
    iLine = iIndex >> MDL_CACHELINE_BITS;
    i = c->pTrans[iLine % c->iTransMask];
    /*
     ** Check for a match!
     */
    while (i) {
		if (c->pTag[i].id == id) {
			if (c->pTag[i].iLine == iLine) {
				++c->pTag[i].nLock;
				pLine = &c->pLine[i*c->iLineSize];
				return(&pLine[iElt*c->iDataSize]);
				}
			}
		i = c->pTag[i].iLink;
		}
    /*
     ** Cache Miss.
     **	Victim Search!
     ** Note: if more than 1771875 cache lines are present this random
     ** number generation may have to be changed, although none of the
     ** code will break in this case. The only problem may be non-optimal
     ** cache line replacement. Maybe give a warning at initialization?
     */
    iVictim = MDL_RAND(mdl)%c->nLines;
    for (i=0;i<c->nLines;++i) {
		if (!c->pTag[iVictim].nLock) {
			/*
			 ** Found victim.
			 */
			iLineVic = c->pTag[iVictim].iLine;
			idVic = c->pTag[iVictim].id;
			/*
			 ** 'pLine' will point to the actual data line in the cache.
			 */
			pLine = &c->pLine[iVictim*c->iLineSize];
			if (iLineVic >= 0) {
				/*
				 ** Flush element since it is valid!
				 */
				cflsh = &mdl->pmdl[idVic]->cache[cid];	
				dest = &cflsh->pData[iLineVic*c->iLineSize];
				/*
				 ** Make sure we don't combine beyond the number
				 ** of data elements!
				 */
				j = iLineVic*c->nLineElts;
				n = j + c->nLineElts;
				if (n > cflsh->nData) n = cflsh->nData;
				n -= j;
				n *= c->iDataSize;
				/*
				 ** Lock data line and combine!
				 */
				cps_mutex_lock(&cflsh->pMux[iLineVic]);
				for (j=0;j<n;j+=c->iDataSize) {
					(*c->combine)(&dest[j],&pLine[j]);
					}
				cps_mutex_unlock(&cflsh->pMux[iLineVic]);
				/*
				 ** If valid iLine then "unlink" it from the cache.
				 */
				pi = &c->pTrans[iLineVic % c->iTransMask];
				while (*pi != iVictim) pi = &c->pTag[*pi].iLink;
				*pi = c->pTag[iVictim].iLink;
				}
			c->pTag[iVictim].id = id;	
			c->pTag[iVictim].iLine = iLine;
			c->pTag[iVictim].nLock = 1;
			/*
			 **	Add the modified victim tag back into the cache.
			 ** Note: the new element is placed at the head of the chain.
			 */
			pi = &c->pTrans[iLine % c->iTransMask];
			c->pTag[iVictim].iLink = *pi;
			*pi = iVictim;
			/*
			 ** Grab new cache line, don't really need a lock here.
			 */
			srce = &cc->pData[iLine*c->iLineSize];
			for (j=0;j<c->iLineSize;++j) {
				pLine[j] = srce[j];
				}
			/*
			 ** Call the initializer function for all elements in 
			 ** the cache line.
			 */
			for (j=0;j<c->iLineSize;j+=c->iDataSize) {
				(*c->init)(&pLine[j]);
				}
			return(&pLine[iElt*c->iDataSize]);
			}
		if (++iVictim == c->nLines) iVictim = 0;
		}
    /*
     ** Cache Failure!
     */
    sprintf(ach,"MDL CACHE FAILURE: cid == %d, no unlocked lines!\n",cid);
    mdlDiag(mdl,ach);
    exit(1);
    }


void mdlRelease(MDL mdl,int cid,void *p)
{
    CACHE *c = &mdl->cache[cid];
    int iLine,iData;
    
    if (c->iType == MDL_ROCACHE) return;
    iLine = ((char *)p - c->pLine) / c->iLineSize;
    /*
     ** No "local direct" pointers are allowed here!
     */
    if (iLine > 0 && iLine < c->nLines) {
		--c->pTag[iLine].nLock;
		assert(c->pTag[iLine].nLock >= 0);
		}
    else assert(0);
    }


void mdlCacheCheck(MDL mdl)
{
	}


double mdlNumAccess(MDL mdl,int cid)
{
	return(0.0);
	}


double mdlMissRatio(MDL mdl,int cid)
{
	return(0.0);
	}


double mdlCollRatio(MDL mdl,int cid)
{
	return(0.0);
	}


double mdlMinRatio(MDL mdl,int cid)
{
	return(0.0);
	}








