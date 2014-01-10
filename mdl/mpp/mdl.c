/*
 ** T3D pvm/shmem version.  Meant to run in local mode (i.e. no YMP parent)
 */
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <malloc.h>
#include <assert.h>
#include <stdarg.h>
#include <sys/time.h>
#include <sys/resource.h>
#include "mdl.h"
#include <mpp/pvm3.h>
#include <mpp/time.h>
#include <mpp/limits.h>
#include <mpp/shmem.h>

#define MDL_NOCACHE			0
#define MDL_ROCACHE			1
#define MDL_COCACHE			2


#define MDL_DEFAULT_BYTES		4096
#define MDL_DEFAULT_SERVICES	50
#define MDL_DEFAULT_CACHEIDS	5

#define PVM_TRANS_SIZE		50000 
#define MDL_TAG_INIT 		1
#define MDL_TAG_SWAPINIT 	2
#define MDL_TAG_SWAP		3
#define MDL_TAG_REQ	   		4
#define MDL_TAG_RPL			5

#define MDL_MAX_PES                 2048

static long shmem_array[MDL_MAX_PES];
#pragma _CRI cache_align pSync
static long pSync[_SHMEM_COLLECT_SYNC_SIZE];

void _srvNull(void *p1,void *vin,int nIn,void *vout,int *pnOut)
{
	return;
	}


double mdlCpuTimer(MDL mdl)
{
	return( ((double) clock())/CLOCKS_PER_SEC);
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


int mdlInitialize(MDL *pmdl,char **argv,void (*fcnChild)(MDL))
{
	MDL mdl;
	int tidSelf;
	/* int i,nThreads,nArch,bid,bDiag,bThreads,iLen;
	   struct pvmhostinfo *hostp; */
	/* Don't need nArch and hostp anymore.  Use group */
	int i,nThreads,bid,bDiag,bThreads,iLen;
	char *group=NULL;
	char *p,ach[256],achDiag[256],name[256];

	for (i=0;i<_SHMEM_COLLECT_SYNC_SIZE;++i) {
	    pSync[i]=_SHMEM_SYNC_VALUE;
		}
	*pmdl = NULL;
	mdl = malloc(sizeof(struct mdlContext));
	assert(mdl != NULL);
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
	mdl->psrv[0].fcnService = _srvNull;
	/*
	 ** Allocate service buffers.
	 */
	mdl->pszIn = malloc(mdl->nMaxInBytes);
	assert(mdl->pszIn != NULL);
	mdl->pszOut = malloc(mdl->nMaxOutBytes);
	assert(mdl->pszOut != NULL);
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
		tidSelf = pvm_mytid();
		/* pvm_config(&nThreads,&nArch,&hostp); */
		/* pvm_config doesn't work on MPP.  Use pvm_gsize instead */
		nThreads=pvm_gsize(group);
		}
	else if (nThreads > 1) {
		tidSelf = pvm_mytid();
		}
	/* Use PE number instead of pvm_tid's */
	tidSelf = pvm_get_PE(tidSelf);
	mdl->bDiag = bDiag;
	mdl->nThreads = nThreads;
	mdl->atid = malloc(mdl->nThreads*sizeof(int));
	assert(mdl->atid != NULL);
	*pmdl = mdl;
	if (nThreads > 1) {
	    /* "Parent" thread is PE 0 */
		if (tidSelf == 0) {
			/*
			 ** Parent thread, start children.
			 */
			mdl->idSelf = 0;
			mdl->atid[mdl->idSelf] = tidSelf;
			/*
			 ** pvm_spawn don't work on the MPP.  So we make PE 0
			 ** the "parent" thread and just use everyone's PE
			 ** number as task ID's for now.
			 pvm_spawn(argv[0],argv,PvmTaskDefault,NULL,mdl->nThreads-1,
			 &mdl->atid[1]);
			 */
			for (i=1;i<mdl->nThreads;++i) {
			    mdl->atid[i]=i;
				}
			pvm_initsend(PvmDataRaw);
			pvm_pkint(mdl->atid,mdl->nThreads,1);
			iLen = strlen(argv[0])+1;
			pvm_pkint(&iLen,1,1);
			pvm_pkbyte(argv[0],iLen,1);
			pvm_mcast(&mdl->atid[1],nThreads-1,MDL_TAG_INIT);
			if (mdl->bDiag) {
				sprintf(achDiag,"%s/%s.%d",ach,argv[0],mdl->idSelf);
				mdl->fpDiag = fopen(achDiag,"w");
				assert(mdl->fpDiag != NULL);
				}
			}
		else {
			/*
			 ** Child thread, get tid array and determine idSelf.
			 */
			bid = pvm_recv(-1,MDL_TAG_INIT);
			pvm_upkint(mdl->atid,mdl->nThreads,1);
			pvm_upkint(&iLen,1,1);
			assert(iLen >= 1);
			pvm_upkbyte(name,iLen,1);
			pvm_freebuf(bid);
			/* for (mdl->idSelf=1;mdl->idSelf<mdl->nThreads;++mdl->idSelf) {
			   if (mdl->atid[mdl->idSelf] == tidSelf) break;
			   } */
			mdl->idSelf = tidSelf;
			if (mdl->bDiag) {
				sprintf(achDiag,"%s/%s.%d",ach,name,mdl->idSelf);
				mdl->fpDiag = fopen(achDiag,"w");
				assert(mdl->fpDiag != NULL);
				}
			(*fcnChild)(mdl);
			mdlFinish(mdl);
			exit(0);
			}
		}
	else {
		/*
		 ** A unik!
		 */
		mdl->idSelf = 0;
		mdl->atid[mdl->idSelf] = 0;
		if (mdl->bDiag) {
			sprintf(achDiag,"%s.%d",ach,mdl->idSelf);
			mdl->fpDiag = fopen(achDiag,"w");
			assert(mdl->fpDiag != NULL);
			}
		}
	return(nThreads);
	}


void mdlFinish(MDL mdl)
{
	/*
	 ** Close Diagnostic file.
	 */
	if (mdl->bDiag) {
		fclose(mdl->fpDiag);
		}
	/*
	 ** Deregister from PVM and deallocate storage.
	 */
	if (mdl->nThreads > 1) pvm_exit();
	free(mdl->psrv);
	free(mdl->pszIn);
	free(mdl->pszOut);
	free(mdl->cache);
	free(mdl->atid);
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
	int tid,bid,nInBytes,nOutBufBytes,nInMax,nOutMax;
	char *pszBuf = vBuf;
	char *pszIn,*pszOut;
	int cc;

	tid = mdl->atid[id];
	*pnRcvBytes = 0;
	*pnSndBytes = 0;
	/*
	 **	Send number of rejects to target thread amount of free space
	 */ 
	cc = pvm_initsend(PvmDataRaw);
	assert(cc >= 0);
	cc = pvm_pkint(&nOutBytes,1,1);
	assert(cc >= 0);
	cc = pvm_pkint(&nBufBytes,1,1);
	assert(cc >= 0);
	cc = pvm_send(tid,MDL_TAG_SWAPINIT);
	assert(cc >= 0);
	/*
	 ** Receive the number of target thread rejects and target free space
	 */
	bid = pvm_recv(tid,MDL_TAG_SWAPINIT);
	assert(bid >= 0);
	cc = pvm_upkint(&nInBytes,1,1);
	assert(cc >= 0);
	assert(nInBytes >= 0);
	cc = pvm_upkint(&nOutBufBytes,1,1);
	assert(cc >= 0);
	assert(nOutBufBytes >= 0);
	cc = pvm_freebuf(bid);
	assert(cc >= 0);
	/*
	 ** Start bilateral transfers. Note: One processor is GUARANTEED to 
	 ** complete all its transfers.
	 */
	pszOut = &pszBuf[nBufBytes-nOutBytes];
	pszIn = pszBuf;
	while (nOutBytes && nInBytes) {
		/*
		 ** nOutMax is the maximum number of bytes allowed to be sent
		 ** nInMax is the number of bytes which will be received.
		 */
		nOutMax = (nOutBytes < PVM_TRANS_SIZE)?nOutBytes:PVM_TRANS_SIZE;
		nOutMax = (nOutMax < nOutBufBytes)?nOutMax:nOutBufBytes;
		nInMax = (nInBytes < PVM_TRANS_SIZE)?nInBytes:PVM_TRANS_SIZE;
		nInMax = (nInMax < nBufBytes)?nInMax:nBufBytes;
		cc = pvm_initsend(PvmDataRaw);
		assert(cc >= 0);
		cc = pvm_pkbyte(pszOut,nOutMax,1);
		assert(cc >= 0);
		cc = pvm_send(tid,MDL_TAG_SWAP);
		assert(cc >= 0);
		bid = pvm_recv(tid,MDL_TAG_SWAP);
		assert(bid >= 0);
		cc = pvm_upkbyte(pszIn,nInMax,1);
		assert(cc >= 0);
		cc = pvm_freebuf(bid);
		assert(cc >= 0);
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
		nOutMax = (nOutBytes < PVM_TRANS_SIZE)?nOutBytes:PVM_TRANS_SIZE;
		nOutMax = (nOutMax < nOutBufBytes)?nOutMax:nOutBufBytes;
		cc = pvm_initsend(PvmDataRaw);
		assert(cc >= 0);
		cc = pvm_pkbyte(pszOut,nOutMax,1);
		assert(cc >= 0);
		cc = pvm_send(tid,MDL_TAG_SWAP);
		assert(cc >= 0);
		pszOut = &pszOut[nOutMax];
		nOutBytes -= nOutMax;
		nOutBufBytes -= nOutMax;
		*pnSndBytes += nOutMax;
		}
	while (nInBytes && nBufBytes) {
		nInMax = (nInBytes < PVM_TRANS_SIZE)?nInBytes:PVM_TRANS_SIZE;
		nInMax = (nInMax < nBufBytes)?nInMax:nBufBytes;
		bid = pvm_recv(tid,MDL_TAG_SWAP);
		assert(bid >= 0);
		cc = pvm_upkbyte(pszIn,nInMax,1);
		assert(cc >= 0);
		cc = pvm_freebuf(bid);
		assert(cc >= 0);
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
	if (nInBytes > mdl->nMaxInBytes) {
		mdl->pszIn = realloc(mdl->pszIn,nInBytes);
		assert(mdl->pszIn != NULL);
		mdl->nMaxInBytes = nInBytes;
		}
	if (nOutBytes > mdl->nMaxOutBytes) {
		mdl->pszOut = realloc(mdl->pszOut,nOutBytes);
		assert(mdl->pszIn != NULL);
		mdl->nMaxOutBytes = nOutBytes;
		}
	mdl->psrv[sid].p1 = p1;
	mdl->psrv[sid].nInBytes = nInBytes;
	mdl->psrv[sid].nOutBytes = nOutBytes;
	mdl->psrv[sid].fcnService = fcnService;
	}


void mdlReqService(MDL mdl,int id,int sid,void *vin,int nInBytes)
{
	char *pszIn = vin;
	int tid;
	int cc;

	tid = mdl->atid[id];
	cc = pvm_initsend(PvmDataRaw);
	assert(cc >= 0);
	cc = pvm_pkint(&mdl->idSelf,1,1);
	assert(cc >= 0);
	cc = pvm_pkint(&sid,1,1);
	assert(cc >= 0);
	if (!pszIn) nInBytes = 0;
	assert(nInBytes <= mdl->psrv[sid].nInBytes);
	cc = pvm_pkint(&nInBytes,1,1);
	assert(cc >= 0);
	if (nInBytes > 0) cc = pvm_pkbyte(pszIn,nInBytes,1);
	assert(cc >= 0);
	cc = pvm_send(tid,MDL_TAG_REQ);
	assert(cc >= 0);
	}


void mdlGetReply(MDL mdl,int id,void *vout,int *pnOutBytes)
{
	char *pszOut = vout;
	int tid,bid,nOutBytes;
	int cc;

	tid = mdl->atid[id];
	bid = pvm_recv(tid,MDL_TAG_RPL);
	assert(bid >= 0);
	cc = pvm_upkint(&nOutBytes,1,1);
	assert(cc >= 0);
	assert(nOutBytes >= 0);
	if (nOutBytes > 0 && pszOut != NULL) cc = pvm_upkbyte(pszOut,nOutBytes,1);
	assert(cc >= 0);
	cc = pvm_freebuf(bid);
	assert(cc >= 0);
	if (pnOutBytes) *pnOutBytes = nOutBytes;
	}


void mdlHandler(MDL mdl)
{
	int bid,tid;
	int nInBytes,nOutBytes,sid,idFrom;
	int cc;

	sid = 1;
	while (sid != SRV_STOP) {
		bid = pvm_recv(-1,MDL_TAG_REQ);
		assert(bid >= 0);
		cc = pvm_upkint(&idFrom,1,1);
		assert(cc >= 0);
		assert(idFrom >= 0);
		assert(idFrom < mdl->nThreads);
		assert(idFrom != mdl->idSelf);
		cc = pvm_upkint(&sid,1,1);
		assert(cc >= 0);
		assert(sid < mdl->nMaxServices);
		assert(sid >= 0);
		cc = pvm_upkint(&nInBytes,1,1);
		assert(cc >= 0);
		assert(nInBytes <= mdl->psrv[sid].nInBytes);
		assert(nInBytes >= 0);
		if (nInBytes > 0) cc = pvm_upkbyte(mdl->pszIn,nInBytes,1);
		assert(cc >= 0);
		cc = pvm_freebuf(bid);
		assert(cc >= 0);
		nOutBytes = 0;
		assert(mdl->psrv[sid].fcnService != NULL);
		(*mdl->psrv[sid].fcnService)(mdl->psrv[sid].p1,mdl->pszIn,nInBytes,
									 mdl->pszOut,&nOutBytes);
		assert(nOutBytes <= mdl->psrv[sid].nOutBytes);
		tid = mdl->atid[idFrom];
		cc = pvm_initsend(PvmDataRaw);
		assert(cc >= 0);
		cc = pvm_pkint(&nOutBytes,1,1);
		assert(cc >= 0);
		if (nOutBytes > 0) cc = pvm_pkbyte(mdl->pszOut,nOutBytes,1);
		assert(cc >= 0);
		cc = pvm_send(tid,MDL_TAG_RPL);
		assert(cc >= 0);
		}
	}

#define MDL_CACHE_SIZE		2000000
#define MDL_CACHELINE_BITS	0
#define MDL_CACHE_MASK		((1<<MDL_CACHELINE_BITS)-1)

#define MDL_TAG_CACHECOM	10
#define MDL_MID_CACHEIN		1
#define MDL_MID_CACHEREQ	2
#define MDL_MID_CACHERPL	3
#define MDL_MID_CACHEOUT	4
#define MDL_MID_CACHEFLSH	5

#define MDL_CHECK_INTERVAL  100
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
    int i;
    long shmax;

    barrier();
    shmem_fcollect(shmem_array,(long *)&iSize,1,0,0,mdl->nThreads,pSync);
    shmax=0;
    for (i=0;i<mdl->nThreads;++i) {
	if (shmax < shmem_array[i]) shmax=shmem_array[i];
    }
    barrier();
    return(shmalloc(shmax));
}


void mdlFree(MDL mdl,void *p)
{
	shfree(p);
	}



/*
 ** Initialize a caching space.
 */
void mdlROcache(MDL mdl,int cid,void *pData,int iDataSize,int nData)
{
	CACHE *c;
	int i,id,nMaxCacheIds;
	char ach[80];
/*	long shmax; */

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
	c->nLineElts = (1 << MDL_CACHELINE_BITS); 
	c->iLineSize = c->nLineElts*c->iDataSize;
	/*
	 ** Determine the number of cache lines to be allocated.
	 */
	c->nLines = (MDL_CACHE_SIZE/c->iDataSize) >> MDL_CACHELINE_BITS;
	assert(c->nLines < MDL_RANDMOD);
	c->nTrans = 1;
	while(c->nTrans < c->nLines) c->nTrans *= 2;
	c->iTransMask = c->nTrans - 1;
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
		c->pTag[i].iKey = -1;	     
		c->pTag[i].nLock = 0;
		c->pTag[i].nLast = 0;
		c->pTag[i].iLink = 0;
		}
	c->pTag[0].nLock = 1;		/* always locked */
	c->pTag[0].nLast = INT_MAX;	/* always Most Recently Used */
	c->nAccess = 0;
	c->nAccHigh = 0;
	c->nMiss = 0;
	c->nColl = 0;
	c->nMin = 0;
	/*
	 ** Allocate cache data lines.
	 */
	barrier();
	shmem_fcollect(shmem_array,(long *)&nData,1,0,0,mdl->nThreads,pSync);
	c->pDataMax=0;
        for (i=0;i<mdl->nThreads;++i) {
	    if (c->pDataMax < shmem_array[i]) c->pDataMax=shmem_array[i];
	}
	c->pDataMax *= iDataSize;

/*	barrier();
	c->pData=shmalloc(shmax*iDataSize);
	assert(c->pData != NULL);
	memcpy(c->pData,pData,nData*iDataSize);
*/
	
	c->pLine = malloc(c->nLines*c->iLineSize);
	assert(c->pLine != NULL); 
	sprintf(ach,"PE %d: c->pData=%ld\n",mdl->idSelf,(long)c->pData);
	mdlDiag(mdl,ach); 
	barrier();
	}


void mdlFinishCache(MDL mdl,int cid)
{
	CACHE *c = &mdl->cache[cid];
	int i,id,iLine;

	barrier();
	/*
	 ** Free up storage and finish.
	 */
	free(c->pTrans);
	free(c->pTag);
	free(c->pLine);
	c->iType = MDL_NOCACHE;
	AdjustDataSize(mdl);
	}


/* No need to service with shmem's */
void mdlCacheCheck(MDL mdl)
{
	return;
	}


void *mdlAquire(MDL mdl,int cid,int iIndex,int id)
{
	CACHE *c = &mdl->cache[cid];
	char *pLine;
	int iElt,iLine,i;
	int iVictim,iLineVic,*pi;
	int iKey,iKeyVic,idVic;
	char ach[100];
	long iLineSize_64,iLineSize_8;

	++c->nAccess;
	/*
	 ** Determine memory block key value and cache line.
	 */
	iElt = iIndex & MDL_CACHE_MASK;
	iLine = iIndex >> MDL_CACHELINE_BITS;
	iKey = iLine * mdl->nThreads+id;
	i = c->pTrans[iKey & c->iTransMask];
	/*
	 ** Check for a match!
	 */
	if (c->pTag[i].iKey == iKey) {
		++c->pTag[i].nLock;
		c->pTag[i].nLast = c->nAccess;
		pLine = &c->pLine[i*c->iLineSize];
		return(&pLine[iElt*c->iDataSize]);
		}
	i = c->pTag[i].iLink;
	while (i) {
		++c->nColl;
		if (c->pTag[i].iKey == iKey) {
		    ++c->pTag[i].nLock;
		    c->pTag[i].nLast = c->nAccess;
		    pLine = &c->pLine[i*c->iLineSize];
		    return(&pLine[iElt*c->iDataSize]);
			}
		i = c->pTag[i].iLink;
		}
	/*
	 ** Is it a local request?
	 */
	if (id == mdl->idSelf) {
		return(&c->pData[iIndex*c->iDataSize]);
		}
	/*
	 ** Cache Miss.
	 */
	++c->nMiss;
	/*
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
			iKeyVic = c->pTag[iVictim].iKey;
			/*
			 ** 'pLine' will point to the actual data line in the cache.
			 */
			pLine = &c->pLine[iVictim*c->iLineSize];
			if (iKeyVic >= 0) {
				if (c->iType == MDL_COCACHE) {
					}
				/*
				 ** If valid iLine then "unlink" it from the cache.
				 */
				pi = &c->pTrans[iKeyVic & c->iTransMask];
				while (*pi != iVictim) pi = &c->pTag[*pi].iLink;
				*pi = c->pTag[iVictim].iLink;
				}
			c->pTag[iVictim].iKey=iKey;
			c->pTag[iVictim].nLock = 1;
			c->pTag[iVictim].nLast = c->nAccess;
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
	sprintf(ach,"MDL CACHE FAILURE: cid == %d, no unlocked lines!\n",cid);
	mdlDiag(mdl,ach);
	exit(1);
 Await:
	/*
	 ** At this point 'pLine' is the recipient cache line for the 
	 ** data requested from processor 'id'.
	 */	
	if ((iLine+1)*c->iLineSize > c->pDataMax) {
	    iLineSize_8 = c->pDataMax - iLine*c->iLineSize;
	}
	else {
	    iLineSize_8 = c->iLineSize;
	} 

	assert( (iLineSize_8%8) == 0 );
	iLineSize_64=iLineSize_8/8;

/*	sprintf(ach,"PE%d: %ld,%d<-PE%d: %ld,\n",mdl->idSelf,
		(long)pLine,iVictim,id,(long)(c->addr[id] + iLine*c->iLineSize),iLine);
	mdlDiag(mdl,ach); */
	shmem_get((long *)pLine,(long *)(c->pData + iLine*c->iLineSize),
			  iLineSize_64,id);
/*	sprintf(ach,"PE%d: %ld <- PE%d done\n",mdl->idSelf,
		(long)pLine,id);
	mdlDiag(mdl,ach); */

	return(&pLine[iElt*c->iDataSize]);
	}


void mdlRelease(MDL mdl,int cid,void *p)
{
	CACHE *c = &mdl->cache[cid];
	int iLine,iData;
	
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

	return(0.0);
	}






