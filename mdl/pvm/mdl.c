/*
 ** A very low level Machine model. For homogeneous systems only!
 */
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#ifdef sgi
#include <strings.h> /* for rindex() */
#endif
#include <malloc.h>
#include <math.h>
#include <limits.h>
#include <assert.h>
#include <stdarg.h>
#include <sys/time.h>
#include <sys/resource.h>
#include "mdl.h"
#include "pvm3.h"


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


void _srvNull(void *p1,void *vin,int nIn,void *vout,int *pnOut)
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


int mdlInitialize(MDL *pmdl,char **argv,void (*fcnChild)(MDL))
{
	MDL mdl;
	int tidSelf=-1;
	int i,nThreads,nArch,bid,bDiag,bThreads,iLen;
	struct pvmhostinfo *hostp;
	char *p,ach[256],achDiag[256],name[256];

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
		pvm_config(&nThreads,&nArch,&hostp);
		}
    else if (nThreads > 1) {
		tidSelf = pvm_mytid();
		}
	mdl->bDiag = bDiag;
	mdl->nThreads = nThreads;
	mdl->atid = malloc(mdl->nThreads*sizeof(int));
	assert(mdl->atid != NULL);
	*pmdl = mdl;
	if (nThreads > 1) {
		if (pvm_parent() < 0) {
			/*
			 ** Parent thread, start children.
			 */
			mdl->idSelf = 0;
			mdl->atid[mdl->idSelf] = tidSelf;
			pvm_spawn(argv[0],argv,PvmTaskDefault,NULL,mdl->nThreads-1,
					  &mdl->atid[1]);
			pvm_initsend(PvmDataRaw);
			pvm_pkint(mdl->atid,mdl->nThreads,1);
			iLen = strlen(argv[0])+1;
			pvm_pkint(&iLen,1,1);
			pvm_pkbyte(argv[0],iLen,1);
			pvm_mcast(&mdl->atid[1],nThreads-1,MDL_TAG_INIT);
			if (mdl->bDiag) {
			        if(rindex(argv[0], '/'))
				    sprintf(achDiag,"%s/%s.%d",ach,
					    rindex(argv[0],'/')+1,mdl->idSelf);
				else
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
			pvm_upkbyte(name,iLen,1);
			pvm_freebuf(bid);
			for (mdl->idSelf=1;mdl->idSelf<mdl->nThreads;++mdl->idSelf) {
				if (mdl->atid[mdl->idSelf] == tidSelf) break;
				}
			if (mdl->bDiag) {
			        if(rindex(name, '/'))
				    sprintf(achDiag,"%s/%s.%d",ach,
					    rindex(name,'/')+1,mdl->idSelf);
				else
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
			char *tmp = strrchr(argv[0],'/');
			if (!tmp) tmp = argv[0];
			else ++tmp;
			sprintf(achDiag,"%s/%s.%d",ach,tmp,mdl->idSelf);
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

	tid = mdl->atid[id];
	*pnRcvBytes = 0;
	*pnSndBytes = 0;
	/*
	 **	Send number of rejects to target thread amount of free space
	 */ 
	pvm_initsend(PvmDataRaw);
	pvm_pkint(&nOutBytes,1,1);
	pvm_pkint(&nBufBytes,1,1);
	pvm_send(tid,MDL_TAG_SWAPINIT);
	/*
	 ** Receive the number of target thread rejects and target free space
	 */
	bid = pvm_recv(tid,MDL_TAG_SWAPINIT);
	pvm_upkint(&nInBytes,1,1);
	pvm_upkint(&nOutBufBytes,1,1);
	pvm_freebuf(bid);
	/*
	 ** Start bilateral transfers. Note: One processor is GUARANTEED to 
	 ** complete all its transfers.
	 */
	assert(nBufBytes >= nOutBytes);
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
		pvm_initsend(PvmDataRaw);
		pvm_pkbyte(pszOut,nOutMax,1);
		pvm_send(tid,MDL_TAG_SWAP);
		bid = pvm_recv(tid,MDL_TAG_SWAP);
		pvm_upkbyte(pszIn,nInMax,1);
		pvm_freebuf(bid);
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
		pvm_initsend(PvmDataRaw);
		pvm_pkbyte(pszOut,nOutMax,1);
		pvm_send(tid,MDL_TAG_SWAP);
		pszOut = &pszOut[nOutMax];
		nOutBytes -= nOutMax;
		nOutBufBytes -= nOutMax;
		*pnSndBytes += nOutMax;
		}
	while (nInBytes && nBufBytes) {
		nInMax = (nInBytes < PVM_TRANS_SIZE)?nInBytes:PVM_TRANS_SIZE;
		nInMax = (nInMax < nBufBytes)?nInMax:nBufBytes;
		bid = pvm_recv(tid,MDL_TAG_SWAP);
		pvm_upkbyte(pszIn,nInMax,1);
		pvm_freebuf(bid);
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
		assert(mdl->pszOut != NULL);
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

	tid = mdl->atid[id];
	pvm_initsend(PvmDataRaw);
	pvm_pkint(&mdl->idSelf,1,1);
	pvm_pkint(&sid,1,1);
	if (!pszIn) nInBytes = 0;
	pvm_pkint(&nInBytes,1,1);
	if (nInBytes > 0) pvm_pkbyte(pszIn,nInBytes,1);
	pvm_send(tid,MDL_TAG_REQ);
	}


void mdlGetReply(MDL mdl,int id,void *vout,int *pnOutBytes)
{
	char *pszOut = vout;
	int tid,bid,nOutBytes;

	tid = mdl->atid[id];
	bid = pvm_recv(tid,MDL_TAG_RPL);
	pvm_upkint(&nOutBytes,1,1);
	if (nOutBytes > 0 && pszOut != NULL) pvm_upkbyte(pszOut,nOutBytes,1);
	pvm_freebuf(bid);
	if (pnOutBytes) *pnOutBytes = nOutBytes;
	}


void mdlHandler(MDL mdl)
{
	int bid,tid;
	int nInBytes,nOutBytes,sid,idFrom;

	sid = 1;
	while (sid != SRV_STOP) {
		bid = pvm_recv(-1,MDL_TAG_REQ);
		pvm_upkint(&idFrom,1,1);
		pvm_upkint(&sid,1,1);
		assert(sid < mdl->nMaxServices);
		pvm_upkint(&nInBytes,1,1);
		assert(nInBytes <= mdl->psrv[sid].nInBytes);
		if (nInBytes > 0) pvm_upkbyte(mdl->pszIn,nInBytes,1);
		pvm_freebuf(bid);
		nOutBytes = 0;
		assert(mdl->psrv[sid].fcnService != NULL);
		(*mdl->psrv[sid].fcnService)(mdl->psrv[sid].p1,mdl->pszIn,nInBytes,
									 mdl->pszOut,&nOutBytes);
		assert(nOutBytes <= mdl->psrv[sid].nOutBytes);
		tid = mdl->atid[idFrom];
		pvm_initsend(PvmDataRaw);
		pvm_pkint(&nOutBytes,1,1);
		if (nOutBytes > 0) pvm_pkbyte(mdl->pszOut,nOutBytes,1);
		pvm_send(tid,MDL_TAG_RPL);
		}
	}


#define MDL_TAG_CACHECOM	10
#define MDL_MID_CACHEIN		1
#define MDL_MID_CACHEREQ	2
#define MDL_MID_CACHERPL	3
#define MDL_MID_CACHEOUT	4
#define MDL_MID_CACHEFLSH	5

#define MDL_CHECK_MASK  	0x7f
#define BILLION				1000000000

int mdlCacheReceive(MDL mdl,char *pLine)
{
	CACHE *c;
	char *t;
	int bid,msg[4],n,i;

	bid = pvm_recv(-1,MDL_TAG_CACHECOM);
	pvm_upkint(msg,4,1);
	c = &mdl->cache[msg[0]];
	switch (msg[1]) {
	case MDL_MID_CACHEIN:
		pvm_freebuf(bid);
		++c->nCheckIn;
		break;
	case MDL_MID_CACHEOUT:
		pvm_freebuf(bid);
		++c->nCheckOut;
		break;
	case MDL_MID_CACHEREQ:
		pvm_freebuf(bid);
		pvm_initsend(PvmDataRaw);
		c->rpl[3] = msg[3];
		pvm_pkint(c->rpl,4,1);
		/*
		 ** Make sure we don't read beyond the number of data elements!
		 */
		if ((msg[3]+1)*MDL_CACHELINE_ELTS > c->nData) {
			/*
			 ** Copy data to the sentinel line before packing.
			 */
			i = msg[3]*MDL_CACHELINE_ELTS;
			t = &c->pData[i*c->iDataSize];
			n = c->nData - i;
			n *= c->iDataSize;
			for (i=0;i<n;++i) c->pLine[i] = t[i];
			for (;i<c->iLineSize;++i) c->pLine[i] = 0;			
			pvm_pkbyte(c->pLine,c->iLineSize,1);
			}
		else {
			pvm_pkbyte(&c->pData[msg[3]*c->iLineSize],c->iLineSize,1);
			}
		pvm_send(mdl->atid[msg[2]],MDL_TAG_CACHECOM);
		break;
	case MDL_MID_CACHEFLSH:
		assert(c->iType == MDL_COCACHE);
		/*
		 ** Unpack the data into the 'sentinel-line' cache data.
		 */
		pvm_upkbyte(c->pLine,c->iLineSize,1);
		pvm_freebuf(bid);
		i = msg[3]*MDL_CACHELINE_ELTS;
		t = &c->pData[i*c->iDataSize];
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
		break;
	case MDL_MID_CACHERPL:
		/*
		 ** For now assume no prefetching!
		 ** This means that this MUST be the reply to this Aquire
		 ** request.
		 */
		assert(pLine != NULL);
		pvm_upkbyte(pLine,c->iLineSize,1);
		pvm_freebuf(bid);
		if (c->iType == MDL_COCACHE) {
			/*
			 ** Call the initializer function for all elements in 
			 ** the cache line.
			 */
			for (i=0;i<c->iLineSize;i+=c->iDataSize) {
				(*c->init)(&pLine[i]);
				}
			}
		return(1);
		}
	return(0);
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
	c->iLineSize = MDL_CACHELINE_ELTS*c->iDataSize;
	/*
	 ** Determine the number of cache lines to be allocated.
	 */
	c->nLines = (MDL_CACHE_SIZE/c->iDataSize) >> MDL_CACHELINE_BITS;
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
	for (i=0;i<c->nLines;++i) {
		c->pTag[i].iKey = -1;	/* invalid */	
		c->pTag[i].nLock = 0;
		c->pTag[i].nLast = 0;
		c->pTag[i].iLink = 0;
		}
	c->pTag[0].nLock = 1;			/* always locked */
	c->pTag[0].nLast = INT_MAX;  	/* always Most Recently Used */
	c->nAccess = 0;
	c->nAccHigh = 0;
	c->nMiss = 0;
	c->nColl = 0;
	c->nMin = 0;
	c->nKeyMax = 500;
	c->pbKey = malloc(c->nKeyMax);
	assert(c->pbKey != NULL);
	for (i=0;i<c->nKeyMax;++i) c->pbKey[i] = 0;
	/*
	 ** Allocate cache data lines.
	 */
	c->pLine = malloc(c->nLines*c->iLineSize);
	assert(c->pLine != NULL);
	c->req[0] = cid;
	c->req[1] = MDL_MID_CACHEREQ;
	c->req[2] = mdl->idSelf;
	c->rpl[0] = cid;
	c->rpl[1] = MDL_MID_CACHERPL;
	c->rpl[2] = mdl->idSelf;
	c->chko[0] = cid;
	c->chko[1] = MDL_MID_CACHEOUT;
	c->chko[2] = mdl->idSelf;
	c->chki[0] = cid;
	c->chki[1] = MDL_MID_CACHEIN;
	c->chki[2] = mdl->idSelf;
	c->flsh[0] = cid;
	c->flsh[1] = MDL_MID_CACHEFLSH;
	c->flsh[2] = mdl->idSelf;
	c->nCheckIn = 0;
	c->nCheckOut = 0;
	return(c);
	}


/*
 ** Initialize a Read-Only caching space.
 */
void mdlROcache(MDL mdl,int cid,void *pData,int iDataSize,int nData)
{
	CACHE *c;
	int id;

	c = CacheInitialize(mdl,cid,pData,iDataSize,nData);
	c->iType = MDL_ROCACHE;
	c->init = NULL;
	c->combine = NULL;
	/*
	 ** THIS IS A SYNCHRONIZE!!!
	 */
	if(mdl->idSelf == 0) {
	    c->nCheckIn = 1;
	    while(c->nCheckIn < mdl->nThreads)
		mdlCacheReceive(mdl, NULL);
	    }
	else {
	    pvm_initsend(PvmDataRaw);
	    pvm_pkint(c->chki,4,1);
	    pvm_send(mdl->atid[0],MDL_TAG_CACHECOM);
	    }
	if(mdl->idSelf == 0) {
	    for(id = 1; id < mdl->nThreads; id++) {
		pvm_initsend(PvmDataRaw);
		pvm_pkint(c->chki,4,1);
		pvm_send(mdl->atid[id],MDL_TAG_CACHECOM);
		}
	    }
	else {
	    c->nCheckIn = 0;
	    while (c->nCheckIn == 0) {
		    mdlCacheReceive(mdl,NULL);
		    }	
	    }	
	AdjustDataSize(mdl);
	}


/*
 ** Initialize a Combiner caching space.
 */
void mdlCOcache(MDL mdl,int cid,void *pData,int iDataSize,int nData,
				void (*init)(void *),void (*combine)(void *,void *))
{
	CACHE *c;
	int id;

	c = CacheInitialize(mdl,cid,pData,iDataSize,nData);
	c->iType = MDL_COCACHE;
	assert(init);
	c->init = init;
	assert(combine);
	c->combine = combine;
	/*
	 ** THIS IS A SYNCHRONIZE!!!
	 */
	if(mdl->idSelf == 0) {
	    c->nCheckIn = 1;
	    while(c->nCheckIn < mdl->nThreads)
		mdlCacheReceive(mdl, NULL);
	    }
	else {
	    pvm_initsend(PvmDataRaw);
	    pvm_pkint(c->chki,4,1);
	    pvm_send(mdl->atid[0],MDL_TAG_CACHECOM);
	    }
	if(mdl->idSelf == 0) {
	    for(id = 1; id < mdl->nThreads; id++) {
		pvm_initsend(PvmDataRaw);
		pvm_pkint(c->chki,4,1);
		pvm_send(mdl->atid[id],MDL_TAG_CACHECOM);
		}
	    }
	else {
	    c->nCheckIn = 0;
	    while (c->nCheckIn == 0) {
		    mdlCacheReceive(mdl,NULL);
		    }	
	    }	
	AdjustDataSize(mdl);
	}


void mdlFinishCache(MDL mdl,int cid)
{
	CACHE *c = &mdl->cache[cid];
	int id;
	int i,iKey;

	if (c->iType == MDL_COCACHE) {
		/*
		 ** Must flush all valid data elements.
		 */
		for (i=1;i<c->nLines;++i) {
			iKey = c->pTag[i].iKey;
			if (iKey >= 0) {
				/*
				 ** Flush element since it is valid!
				 */
				id = iKey%mdl->nThreads;
				c->flsh[3] = iKey/mdl->nThreads;
				pvm_initsend(PvmDataRaw);
				pvm_pkint(c->flsh,4,1);
				pvm_pkbyte(&c->pLine[i*c->iLineSize],c->iLineSize,1);
				pvm_send(mdl->atid[id],MDL_TAG_CACHECOM);
				mdlCacheCheck(mdl); /* service incoming */
				}
			}
		}
	/*
	 ** THIS IS A SYNCHRONIZE!!!
	 */
	if(mdl->idSelf == 0) {
	    ++c->nCheckOut;
	    while(c->nCheckOut < mdl->nThreads)
		mdlCacheReceive(mdl, NULL);
	    }
	else {
	    pvm_initsend(PvmDataRaw);
	    pvm_pkint(c->chko,4,1);
	    pvm_send(mdl->atid[0],MDL_TAG_CACHECOM);
	    }
	if(mdl->idSelf == 0) {
	    for(id = 1; id < mdl->nThreads; id++) {
		pvm_initsend(PvmDataRaw);
		pvm_pkint(c->chko,4,1);
		pvm_send(mdl->atid[id],MDL_TAG_CACHECOM);
		}
	    }
	else {
	    c->nCheckOut = 0;
	    while (c->nCheckOut == 0) {
		    mdlCacheReceive(mdl,NULL);
		    }	
	    }	
	/*
	 ** Free up storage and finish.
	 */
	free(c->pbKey);
	free(c->pTrans);
	free(c->pTag);
	free(c->pLine);
	c->iType = MDL_NOCACHE;
	AdjustDataSize(mdl);
	}


void mdlCacheCheck(MDL mdl)
{
	while (pvm_probe(-1,MDL_TAG_CACHECOM)) {
		mdlCacheReceive(mdl,NULL);
		}
	}


void *mdlAquire(MDL mdl,int cid,int iIndex,int id)
{
	CACHE *c = &mdl->cache[cid];
	char *pLine;
	int iElt,iLine,i,iKey,iKeyVic,idVic,nKeyNew;
	int iVictim,*pi;
	char ach[80];

	++c->nAccess;
	if (!(c->nAccess & MDL_CHECK_MASK)) {
		while (pvm_probe(-1,MDL_TAG_CACHECOM)) {
			mdlCacheReceive(mdl,NULL);
			}
		}
	/*
	 ** Determine memory block key value and cache line.
	 */
	iLine = iIndex >> MDL_CACHELINE_BITS;
	iKey = iLine*mdl->nThreads + id;
	/*
	 ** Consider the following:
	 ** iKey = (iIndex << c->iKeyShift) | id;
	 */
	i = c->pTrans[iKey & c->iTransMask];
	/*
	 ** Check for a match!
	 */
	if (c->pTag[i].iKey == iKey) {
		++c->pTag[i].nLock;
		c->pTag[i].nLast = c->nAccess;
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
			c->pTag[i].nLast = c->nAccess;
			pLine = &c->pLine[i*c->iLineSize];
			iElt = iIndex & MDL_CACHE_MASK;
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
	c->req[3] = iLine;
	pvm_initsend(PvmDataRaw);
	pvm_pkint(c->req,4,1);
	pvm_send(mdl->atid[id],MDL_TAG_CACHECOM);
	++c->nMiss;
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
	for (i=1;i<c->nLines;++i) {
		if (c->pTag[i].nLast < c->pTag[iVictim].nLast) {
			if (!c->pTag[i].nLock) iVictim = i;
			}
		}
	if (!iVictim) {
		/*
		 ** Cache Failure!
		 */
		sprintf(ach,"MDL CACHE FAILURE: cid == %d, no unlocked lines!\n",cid);
		mdlDiag(mdl,ach);
		exit(1);
		}
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
			idVic = iKeyVic%mdl->nThreads;
			c->flsh[3] = iKeyVic/mdl->nThreads;
			pvm_initsend(PvmDataRaw);
			pvm_pkint(c->flsh,4,1);
			pvm_pkbyte(pLine,c->iLineSize,1);
			pvm_send(mdl->atid[idVic],MDL_TAG_CACHECOM);
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
	c->pTag[iVictim].nLast = c->nAccess;
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
	if (iKey >= c->nKeyMax) {
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
	while (1) {
		if (mdlCacheReceive(mdl,pLine)) {
			return(&pLine[iElt*c->iDataSize]);
			}
		}
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

	if (dAccess > 0.0) return(c->nMin/dAccess);
	else return(0.0);
	}









