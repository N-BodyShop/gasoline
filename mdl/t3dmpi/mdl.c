/*
 ** MPI and shared memory version of MDL for the T3D.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <malloc.h>
#include <math.h>
#include <assert.h>
#include <stdarg.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <mpp/limits.h>
#include <mpp/shmem.h>
#include "mdl.h"


#define MDL_NOCACHE			0
#define MDL_ROCACHE			1
#define MDL_COCACHE			2


#define MDL_DEFAULT_BYTES		80000
#define MDL_DEFAULT_SERVICES	50
#define MDL_DEFAULT_CACHEIDS	5

#define MDL_TRANS_SIZE		50000 
#define MDL_TAG_INIT 		1
#define MDL_TAG_SWAPINIT 	2
#define MDL_TAG_SWAP		3
#define MDL_TAG_REQ	   		4
#define MDL_TAG_RPL			5


#define MDL_MAX_PES                 2048

static long shmem_array[MDL_MAX_PES];
#pragma _CRI cache_align pSync
static long pSync[_SHMEM_COLLECT_SYNC_SIZE];

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
	int i,bDiag,bThreads,nbuf[4];
	char *p,ach[256],achDiag[256];

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
	mdl->nMaxSrvBytes = MDL_DEFAULT_BYTES;
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
	mdl->pszIn = malloc(mdl->nMaxSrvBytes+sizeof(SRVHEAD));
	assert(mdl->pszIn != NULL);
	mdl->pszOut = malloc(mdl->nMaxSrvBytes+sizeof(SRVHEAD));
	assert(mdl->pszOut != NULL);
	mdl->pszBuf = malloc(mdl->nMaxSrvBytes+sizeof(SRVHEAD));
	assert(mdl->pszBuf != NULL);
	/*
	 ** Allocate swapping transfer buffer. This buffer remains fixed.
	 */
	mdl->pszTrans = malloc(MDL_TRANS_SIZE);
	assert(mdl->pszTrans != NULL);
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
			if (argv[i]) bThreads = 1;
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
	if (bThreads) {
		fprintf(stderr,"Warning: -sz parameter ignored, using as many\n");
		fprintf(stderr,"         processors as specified in environment.\n");
		fflush(stderr);
		}
	MPI_Init(0, NULL);

	MPI_Comm_size(MPI_COMM_WORLD, &mdl->nThreads);
	MPI_Comm_rank(MPI_COMM_WORLD, &mdl->idSelf);
	/*
	 ** Allocate caching buffers, with initial data size of 0.
	 ** We need one reply buffer for each thread, to deadlock situations.
	 */
	mdl->iMaxDataSize = 0;
	mdl->iCaBufSize = sizeof(CAHEAD);
	mdl->pszRcv = malloc(mdl->iCaBufSize);
	assert(mdl->pszRcv != NULL);
	mdl->pszFlsh = malloc(mdl->iCaBufSize);
	assert(mdl->pszFlsh != NULL);
	mdl->bDiag = bDiag;
	*pmdl = mdl;
	if (mdl->nThreads > 1) {
		if (mdl->idSelf == 0) {
			/*
			 ** Master thread.
			 */
			if (mdl->bDiag) {
				sprintf(achDiag,"%s/%s.%d",ach,argv[0],mdl->idSelf);
				mdl->fpDiag = fopen(achDiag,"w");
				assert(mdl->fpDiag != NULL);
				}
			}
		else {
			/*
			 ** Child thread.
			 */
			if (mdl->bDiag) {
				sprintf(achDiag,"%s/%s.%d",ach,argv[0],mdl->idSelf);
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
		if (mdl->bDiag) {
			sprintf(achDiag,"%s.%d",ach,mdl->idSelf);
			mdl->fpDiag = fopen(achDiag,"w");
			assert(mdl->fpDiag != NULL);
			}
		}
	return(mdl->nThreads);
	}


void mdlFinish(MDL mdl)
{
  int i;

	MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();
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
	free(mdl->pszIn);
	free(mdl->pszOut);
	free(mdl->pszBuf);
	free(mdl->pszTrans);
	free(mdl->cache);
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
	int nInBytes,nOutBufBytes,nInMax,nOutMax,i;
	int mid,nBytes,iTag,ret,pid;
	char *pszBuf = vBuf;
	char *pszIn,*pszOut;
	struct swapInit {
		int nOutBytes;
		int nBufBytes;
		} swi,swo;
	MPI_Status status;
	MPI_Request request;

	*pnRcvBytes = 0;
	*pnSndBytes = 0;
	/*
	 **	Send number of rejects to target thread amount of free space
	 */ 
	swi.nOutBytes = nOutBytes;
	swi.nBufBytes = nBufBytes;
	MPI_Isend(&swi,sizeof(swi),MPI_CHAR,id,MDL_TAG_SWAPINIT,
		  MPI_COMM_WORLD, &request);
	/*
	 ** Receive the number of target thread rejects and target free space
	 */
	iTag = MDL_TAG_SWAPINIT;
	pid = id;
	MPI_Recv(&swo,sizeof(swo),MPI_CHAR,pid,iTag,MPI_COMM_WORLD, &status);
	MPI_Wait(&request, &status);
	nInBytes = swo.nOutBytes;
	nOutBufBytes = swo.nBufBytes;
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
		nOutMax = (nOutBytes < MDL_TRANS_SIZE)?nOutBytes:MDL_TRANS_SIZE;
		nOutMax = (nOutMax < nOutBufBytes)?nOutMax:nOutBufBytes;
		nInMax = (nInBytes < MDL_TRANS_SIZE)?nInBytes:MDL_TRANS_SIZE;
		nInMax = (nInMax < nBufBytes)?nInMax:nBufBytes;
		/*
		 ** Copy to a temp buffer to be safe.
		 */
		for (i=0;i<nOutMax;++i) mdl->pszTrans[i] = pszOut[i];
		MPI_Isend(mdl->pszTrans,nOutMax,MPI_CHAR,id,MDL_TAG_SWAP,
			 MPI_COMM_WORLD, &request);
		iTag = MDL_TAG_SWAP;
		pid = id;
		MPI_Recv(pszIn,nInMax,MPI_CHAR,pid,iTag,MPI_COMM_WORLD,
			 &status);
		MPI_Get_count(&status, MPI_CHAR, &nBytes);
		assert(nBytes == nInMax);
		MPI_Wait(&request, &status);
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
	 ** Note: use of bsend is mandatory here, also because of this we
	 ** don't need to use the intermediate buffer mdl->pszTrans.
	 */
	while (nOutBytes && nOutBufBytes) {
		nOutMax = (nOutBytes < MDL_TRANS_SIZE)?nOutBytes:MDL_TRANS_SIZE;
		nOutMax = (nOutMax < nOutBufBytes)?nOutMax:nOutBufBytes;
		MPI_Send(pszOut,nOutMax,MPI_CHAR,id,MDL_TAG_SWAP,
			 MPI_COMM_WORLD);
		pszOut = &pszOut[nOutMax];
		nOutBytes -= nOutMax;
		nOutBufBytes -= nOutMax;
		*pnSndBytes += nOutMax;
		}
	while (nInBytes && nBufBytes) {
		nInMax = (nInBytes < MDL_TRANS_SIZE)?nInBytes:MDL_TRANS_SIZE;
		nInMax = (nInMax < nBufBytes)?nInMax:nBufBytes;
		iTag = MDL_TAG_SWAP;
		MPI_Recv(pszIn,nInMax,MPI_CHAR,id,iTag,MPI_COMM_WORLD,
			 &status);
		MPI_Get_count(&status, MPI_CHAR, &nBytes);
		assert(nBytes == nInMax);
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
	int i,nMaxServices,nMaxBytes;

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
	nMaxBytes = (nInBytes > nOutBytes)?nInBytes:nOutBytes;
	if (nMaxBytes > mdl->nMaxSrvBytes) {
		mdl->pszIn = realloc(mdl->pszIn,nMaxBytes+sizeof(SRVHEAD));
		assert(mdl->pszIn != NULL);
		mdl->pszOut = realloc(mdl->pszOut,nMaxBytes+sizeof(SRVHEAD));
		assert(mdl->pszOut != NULL);
		mdl->pszBuf = realloc(mdl->pszBuf,nMaxBytes+sizeof(SRVHEAD));
		assert(mdl->pszBuf != NULL);
		mdl->nMaxSrvBytes = nMaxBytes;
		}
	mdl->psrv[sid].p1 = p1;
	mdl->psrv[sid].nInBytes = nInBytes;
	mdl->psrv[sid].nOutBytes = nOutBytes;
	mdl->psrv[sid].fcnService = fcnService;
	}


void mdlReqService(MDL mdl,int id,int sid,void *vin,int nInBytes)
{
	char *pszIn = vin;
	MPI_Status status;
	/*
	 ** If this looks like dangerous magic, it's because it is!
	 */
	SRVHEAD *ph = (SRVHEAD *)mdl->pszBuf;
	char *pszOut = &mdl->pszBuf[sizeof(SRVHEAD)];
	int i,ret;

	ph->idFrom = mdl->idSelf;
	ph->sid = sid;
	assert(nInBytes <= mdl->psrv[sid].nInBytes);
	if (!pszIn) ph->nInBytes = 0;
	else ph->nInBytes = nInBytes;
	if (nInBytes > 0 && pszIn != NULL) {
		for (i=0;i<nInBytes;++i) pszOut[i] = pszIn[i];
		}
	MPI_Send(mdl->pszBuf,nInBytes+sizeof(SRVHEAD),MPI_CHAR,id,MDL_TAG_REQ,
		 MPI_COMM_WORLD);
	}


void mdlGetReply(MDL mdl,int id,void *vout,int *pnOutBytes)
{
	char *pszOut = vout;
	SRVHEAD *ph = (SRVHEAD *)mdl->pszBuf;
	char *pszIn = &mdl->pszBuf[sizeof(SRVHEAD)];
	int i,ret,iTag,nBytes;
	MPI_Status status;

	iTag = MDL_TAG_RPL;
	MPI_Recv(mdl->pszBuf,mdl->nMaxSrvBytes+sizeof(SRVHEAD),MPI_CHAR,
					id,iTag,MPI_COMM_WORLD, &status);
	MPI_Get_count(&status, MPI_CHAR, &nBytes);
	assert(nBytes == ph->nOutBytes + sizeof(SRVHEAD));
	if (ph->nOutBytes > 0 && pszOut != NULL) {
		for (i=0;i<ph->nOutBytes;++i) pszOut[i] = pszIn[i];
		}
	if (pnOutBytes) *pnOutBytes = ph->nOutBytes;
	}


void mdlHandler(MDL mdl)
{
	SRVHEAD *phi = (SRVHEAD *)mdl->pszIn;
	SRVHEAD *pho = (SRVHEAD *)mdl->pszOut;
	char *pszIn = &mdl->pszIn[sizeof(SRVHEAD)];
	char *pszOut = &mdl->pszOut[sizeof(SRVHEAD)];
	int sid,ret,iTag,id,nOutBytes,nBytes;
	MPI_Status status;

	sid = 1;
	while (sid != SRV_STOP) {
		iTag = MDL_TAG_REQ;
		id = MPI_ANY_SOURCE;
		MPI_Recv(mdl->pszIn,mdl->nMaxSrvBytes+sizeof(SRVHEAD),
			       MPI_CHAR, id,iTag,MPI_COMM_WORLD,&status);
		/*
		 ** Quite a few sanity checks follow.
		 */
		id = status.MPI_SOURCE;
		MPI_Get_count(&status, MPI_CHAR, &nBytes);
		assert(nBytes == phi->nInBytes + sizeof(SRVHEAD));
		assert(id == phi->idFrom);
		sid = phi->sid;
		assert(sid < mdl->nMaxServices);
		assert(phi->nInBytes <= mdl->psrv[sid].nInBytes);
		nOutBytes = 0;
		assert(mdl->psrv[sid].fcnService != NULL);
		(*mdl->psrv[sid].fcnService)(mdl->psrv[sid].p1,pszIn,phi->nInBytes,
									 pszOut,&nOutBytes);
		assert(nOutBytes <= mdl->psrv[sid].nOutBytes);
		pho->idFrom = mdl->idSelf;
		pho->sid = sid;
		pho->nInBytes = phi->nInBytes;
		pho->nOutBytes = nOutBytes;
		MPI_Send(mdl->pszOut,nOutBytes+sizeof(SRVHEAD),
			 MPI_CHAR, id,MDL_TAG_RPL, MPI_COMM_WORLD);
		}
	}

#define MDL_TAG_CACHECOM	10
#define MDL_MID_CACHEIN		1
#define MDL_MID_CACHEREQ	2
#define MDL_MID_CACHERPL	3
#define MDL_MID_CACHEOUT	4
#define MDL_MID_CACHEFLSH	5

#define MDL_RANDMOD		1771875
#define MDL_RAND(mdl) (mdl->uRand = (mdl->uRand*2416+374441)%MDL_RANDMOD)
#define MDL_CHECK_MASK  	0x7f
#define BILLION				1000000000

int mdlCacheReceive(MDL mdl,char *pLine)
{
	CACHE *c;
	CAHEAD *ph = (CAHEAD *)mdl->pszRcv;
	char *pszRcv = &mdl->pszRcv[sizeof(CAHEAD)];
	CAHEAD *phRpl;
	char *pszRpl;
	char *t;
	int n,i,msgid,nBytes;
	int id, iTag;
	MPI_Status status;

	id = MPI_ANY_SOURCE;
	iTag = MDL_TAG_CACHECOM;
	MPI_Recv(mdl->pszRcv,mdl->iCaBufSize, MPI_CHAR, id, iTag,
		 MPI_COMM_WORLD, &status);

	c = &mdl->cache[ph->cid];
	switch (ph->mid) {
	case MDL_MID_CACHEIN:
		++c->nCheckIn;
		return(0);
	case MDL_MID_CACHEOUT:
		++c->nCheckOut;
		return(0);
	case MDL_MID_CACHEFLSH:
		assert(c->iType == MDL_COCACHE);
		/*
		 ** Unpack the data into the 'sentinel-line' cache data.
		 */
		for (i=0;i<c->iLineSize;++i) c->pLine[i] = pszRcv[i];
		i = ph->iLine*MDL_CACHELINE_ELTS;
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
		return(0);
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
		mdl->pszRcv = realloc(mdl->pszRcv,mdl->iCaBufSize);
		assert(mdl->pszRcv != NULL);
		mdl->pszFlsh = realloc(mdl->pszFlsh,mdl->iCaBufSize);
		assert(mdl->pszFlsh != NULL);
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

    MPI_Barrier(MPI_COMM_WORLD);
    shmem_fcollect(shmem_array,(long *)&iSize,1,0,0,mdl->nThreads,pSync);
    shmax=0;
    for (i=0;i<mdl->nThreads;++i) {
	if (shmax < shmem_array[i]) shmax=shmem_array[i];
    }
    MPI_Barrier(MPI_COMM_WORLD);
    return(shmalloc(shmax));
}

void mdlFree(MDL mdl,void *p)
{
	shfree(p);
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
	return(c);
	}

/*
 ** Initialize a caching space.
 */
void mdlROcache(MDL mdl,int cid,void *pData,int iDataSize,int nData)
{
	CACHE *c;
	int i,id;
	CAHEAD caIn;

	c = CacheInitialize(mdl,cid,pData,iDataSize,nData);
	c->iType = MDL_ROCACHE;
	/*
	 ** For an ROcache these two functions are not needed.
	 */
	c->init = NULL;
	c->combine = NULL;
	mdlDiag(mdl, "Before RO sync\n");
	/*
	 ** THIS IS A SYNCHRONIZE!!!
	 */
	caIn.cid = cid;
	caIn.mid = MDL_MID_CACHEIN;
	caIn.id = mdl->idSelf;
	if(mdl->idSelf == 0) {
	    c->nCheckIn = 1;
	    while(c->nCheckIn < mdl->nThreads) {
		mdlCacheReceive(mdl, NULL);
		}
	    }
	else {
		/*
		 ** Must use non-blocking sends here, we will never wait
		 ** for these sends to complete, but will know for sure
		 ** that they have completed.
		 */
		MPI_Send(&caIn,sizeof(CAHEAD),MPI_CHAR, 0,
			       MDL_TAG_CACHECOM, MPI_COMM_WORLD);
		}
	mdlDiag(mdl, "After sends RO sync\n");
	if(mdl->idSelf == 0) {
	    for(id = 1; id < mdl->nThreads; id++) {
		MPI_Send(&caIn,sizeof(CAHEAD),MPI_CHAR, id,
			       MDL_TAG_CACHECOM, MPI_COMM_WORLD);
		}
	    }
	else {
	    c->nCheckIn = 0;
	    while (c->nCheckIn == 0) {
		    mdlCacheReceive(mdl,NULL);
		    }	
	    }	
	mdlDiag(mdl, "After RO sync\n");
	AdjustDataSize(mdl);
	shmem_fcollect(shmem_array,(long *)&nData,1,0,0,mdl->nThreads,pSync);
	c->pDataMax=0;
        for (i=0;i<mdl->nThreads;++i) {
	    if (c->pDataMax < shmem_array[i]) c->pDataMax=shmem_array[i];
	    }
	c->pDataMax *= iDataSize;
	mdlDiag(mdl, "After RO shmem_collect\n");
	}

/*
 ** Initialize a combiner caching space.
 */
void mdlCOcache(MDL mdl,int cid,void *pData,int iDataSize,int nData,
				void (*init)(void *),void (*combine)(void *,void *))
{
	CACHE *c;
	int i,id;
	CAHEAD caIn;

	c = CacheInitialize(mdl,cid,pData,iDataSize,nData);
	c->iType = MDL_COCACHE;
	assert(init);
	c->init = init;
	assert(combine);
	c->combine = combine;
	/*
	 ** THIS IS A SYNCHRONIZE!!!
	 */
	caIn.cid = cid;
	caIn.mid = MDL_MID_CACHEIN;
	caIn.id = mdl->idSelf;
	if(mdl->idSelf == 0) {
	    c->nCheckIn = 1;
	    while(c->nCheckIn < mdl->nThreads) {
		mdlCacheReceive(mdl, NULL);
		}
	    }
	else {
		/*
		 ** Must use non-blocking sends here, we will never wait
		 ** for these sends to complete, but will know for sure
		 ** that they have completed.
		 */
		MPI_Send(&caIn,sizeof(CAHEAD),MPI_CHAR, 0,
			       MDL_TAG_CACHECOM, MPI_COMM_WORLD);
		}
	if(mdl->idSelf == 0) {
	    for(id = 1; id < mdl->nThreads; id++) {
		MPI_Send(&caIn,sizeof(CAHEAD),MPI_CHAR, id,
			       MDL_TAG_CACHECOM, MPI_COMM_WORLD);
		}
	    }
	else {
	    c->nCheckIn = 0;
	    while (c->nCheckIn == 0) {
		    mdlCacheReceive(mdl,NULL);
		    }	
	    }	
	AdjustDataSize(mdl);
	shmem_fcollect(shmem_array,(long *)&nData,1,0,0,mdl->nThreads,pSync);
	c->pDataMax=0;
        for (i=0;i<mdl->nThreads;++i) {
	    if (c->pDataMax < shmem_array[i]) c->pDataMax=shmem_array[i];
	}
	c->pDataMax *= iDataSize;
	}

void mdlFinishCache(MDL mdl,int cid)
{
	CACHE *c = &mdl->cache[cid];
	CAHEAD caOut;
	CAHEAD *caFlsh = (CAHEAD *)mdl->pszFlsh;
	char *pszFlsh = &mdl->pszFlsh[sizeof(CAHEAD)];
	int i,id;
	char *t;
	int j, iKey;

	if (c->iType == MDL_COCACHE) {
		/*
		 ** Must flush all valid data elements.
		 */
		caFlsh->cid = cid;
		caFlsh->mid = MDL_MID_CACHEFLSH;
		caFlsh->id = mdl->idSelf;
		for (i=1;i<c->nLines;++i) {
			iKey = c->pTag[i].iKey;
			if (iKey >= 0) {
				/*
				 ** Flush element since it is valid!
				 */
				id = iKey & c->iIdMask;
				caFlsh->iLine = iKey >> c->iInvKeyShift;
				t = &c->pLine[i*c->iLineSize];
				for(j = 0; j < c->iLineSize; ++j)
				    pszFlsh[j] = t[j];
				MPI_Send(caFlsh, sizeof(CAHEAD)+c->iLineSize,
					 MPI_CHAR, id, MDL_TAG_CACHECOM,
					 MPI_COMM_WORLD); 
				mdlCacheCheck(mdl); /* service incoming */
				}
			}
		}
	/*
	 ** THIS IS A SYNCHRONIZE!!!
	 */
	caOut.cid = cid;
	caOut.mid = MDL_MID_CACHEOUT;
	caOut.id = mdl->idSelf;
	if(mdl->idSelf == 0) {
	    ++c->nCheckOut;
	    while(c->nCheckOut < mdl->nThreads) {
		mdlCacheReceive(mdl, NULL);
		}
	    }
	else {
	    MPI_Send(&caOut,sizeof(CAHEAD),MPI_CHAR, 0,
			       MDL_TAG_CACHECOM, MPI_COMM_WORLD);
	    }
	if(mdl->idSelf == 0) {
	    for(id = 1; id < mdl->nThreads; id++) {
		MPI_Send(&caOut,sizeof(CAHEAD),MPI_CHAR, id,
			       MDL_TAG_CACHECOM, MPI_COMM_WORLD);
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
	free(c->pTrans);
	free(c->pTag);
	free(c->pbKey);
	free(c->pLine);
	c->iType = MDL_NOCACHE;
	AdjustDataSize(mdl);
	}


void mdlCacheCheck(MDL mdl)
{
    int iTag, id,flag;
    MPI_Status status;

    while (1) {
	id = MPI_ANY_SOURCE;
	iTag = MDL_TAG_CACHECOM;
	MPI_Iprobe(id, iTag, MPI_COMM_WORLD, &flag, &status);
	if(flag == 0)
	    break;
	mdlCacheReceive(mdl,NULL);
        }
    }

void *doMiss(MDL mdl, int cid, int iIndex, int id, int iKey, int lock);

void *mdlAquire(MDL mdl,int cid,int iIndex,int id)
{
	CACHE *c = &mdl->cache[cid];
	char *pLine;
	int iElt,i,iKey;
	int iLine;

	++c->nAccess;
	if (!(c->nAccess & MDL_CHECK_MASK))
	        mdlCacheCheck(mdl);
	/*
	 ** Determine memory block key value and cache line.
	iKey = iLine*mdl->nThreads + id;
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
	/*
	 ** Is it a local request?
	 */
	if (id == mdl->idSelf) {
		return(&c->pData[iIndex*c->iDataSize]);
		}
	return(doMiss(mdl, cid, iIndex, id, iKey, 1));
}

void *mdlSteal(MDL mdl,int cid,int iIndex,int id)
{
	CACHE *c = &mdl->cache[cid];
	char *pLine;
	int iElt,i,iKey;
	int iLine;

	++c->nAccess;
	/*
	 ** Determine memory block key value and cache line.
	 */
	iKey = ((iIndex&MDL_INDEX_MASK) << c->iKeyShift)| id;
	/*
	 ** Consider the following:
	 ** iKey = (iIndex << c->iKeyShift) | id;
	 */
	i = c->pTrans[iKey & c->iTransMask];
	/*
	 ** Check for a match!
	 */
	if (c->pTag[i].iKey == iKey) {
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
	return(doMiss(mdl, cid, iIndex, id, iKey, 0));
}

void *doMiss(MDL mdl, int cid, int iIndex, int id, int iKey, int lock)
{
	CACHE *c = &mdl->cache[cid];
	char *pLine;
	int iElt,iLine,i,iKeyVic,nKeyNew;
	int idVic;
	int iVictim,*pi;
	char ach[80];
	unsigned long iLineSize_64,iLineSize_8;
	CAHEAD *caFlsh;
	char *pszFlsh;

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
			    idVic = iKeyVic&c->iIdMask;
			    caFlsh = (CAHEAD *)mdl->pszFlsh;
			    pszFlsh = &mdl->pszFlsh[sizeof(CAHEAD)];
			    caFlsh->cid = cid;
			    caFlsh->mid = MDL_MID_CACHEFLSH;
			    caFlsh->id = mdl->idSelf;
			    caFlsh->iLine = iKeyVic >> c->iInvKeyShift;
			    for(i = 0; i < c->iLineSize; ++i)
				pszFlsh[i] = pLine[i];
			    MPI_Send(caFlsh, sizeof(CAHEAD)+c->iLineSize,
				     MPI_CHAR, idVic,
				     MDL_TAG_CACHECOM, MPI_COMM_WORLD); 
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
	sprintf(ach,"MDL CACHE FAILURE: cid == %d, no unlocked lines!\n",cid);
	mdlDiag(mdl,ach);
	exit(1);
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
		}								/* !!! */
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

	assert( (iLineSize_8%sizeof(long)) == 0 );
	iLineSize_64=iLineSize_8/sizeof(long);

	shmem_get((long *)pLine,(long *)(c->pData + iLine*c->iLineSize),
			  iLineSize_64,id);

	if (c->iType == MDL_COCACHE) {
		/*
		 ** Call the initializer function for all elements in 
		 ** the cache line.
		 */
		for (i=0;i<iLineSize_8;i+=c->iDataSize) {
			(*c->init)(&pLine[i]);
			}
		}

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

	if (dAccess > 0.0) return(c->nMin/dAccess);
	else return(0.0);
	}
