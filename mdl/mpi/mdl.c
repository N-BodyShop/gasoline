/*
 ** MPI version of MDL.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <malloc.h>
#include <math.h>
#include <limits.h>
#include <assert.h>
#include <stdarg.h>
#include <sys/time.h>
#include <sys/resource.h>
#include "mpi.h"
#include "mdl.h"

#define MDL_NOCACHE			0
#define MDL_ROCACHE			1
#define MDL_COCACHE			2
#define MDL_DUMCACHE        3

#define MDL_DEFAULT_BYTES		80000
#define MDL_DEFAULT_SERVICES	50
#define MDL_DEFAULT_CACHEIDS	5

#define MDL_TRANS_SIZE		5000000
#define MDL_TAG_INIT 		1
#define MDL_TAG_SWAPINIT 	2
#define MDL_TAG_SWAP		3
#define MDL_TAG_REQ	   		4
#define MDL_TAG_RPL			5


double mdlVersion(MDL mdl) {
     return MDL_VERSION_NUMBER;
}

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

int mdlInitialize(MDL *pmdl,char **argv,void (*fcnChild)(MDL))
{
	MDL mdl;
	int i,bDiag,bThreads;
	char *p,ach[256],achDiag[256];
	int argc;

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
    ** Initialize work spaces. --JPG
    */
    mdl->work.cWorkList = NULL;
    mdl->work.iNextLocalElt = -1;
    mdl->work.iNextRemoteElt = -1;
    mdl->work.nLocalWorkRemaining = -1;
    mdl->work.idScheduler = -1;
    mdl->work.bDynamic = 0;

	for(argc = 0; argv[argc]; argc++);

	MPI_Init(&argc, &argv);

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
	mdl->ppszRpl = malloc(mdl->nThreads*sizeof(char *));
	assert(mdl->ppszRpl != NULL);
	mdl->pmidRpl = malloc(mdl->nThreads*sizeof(int));
	assert(mdl->pmidRpl != NULL);
	for (i=0;i<mdl->nThreads;++i)
		mdl->pmidRpl[i] = -1;
	mdl->pReqRpl = malloc(mdl->nThreads*sizeof(MPI_Request));
	assert(mdl->pReqRpl != NULL);
	for (i=0;i<mdl->nThreads;++i) {
		mdl->ppszRpl[i] = malloc(mdl->iCaBufSize);
		assert(mdl->ppszRpl[i] != NULL);
		}
	mdl->pszFlsh = malloc(mdl->iCaBufSize);
	assert(mdl->pszFlsh != NULL);
	mdl->bDiag = bDiag;
	*pmdl = mdl;
	if (mdl->bDiag) {
		char *tmp = strrchr(argv[0],'/');
		if (!tmp) tmp = argv[0];
		else ++tmp;
		sprintf(achDiag,"%s/%s.%d",ach,tmp,mdl->idSelf);
		mdl->fpDiag = fopen(achDiag,"w");
		assert(mdl->fpDiag != NULL);
		}
	if (mdl->nThreads > 1 && mdl->idSelf) {
		/*
		 ** Child thread.
		 */
		(*fcnChild)(mdl);
		mdlFinish(mdl);
		exit(0);
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
	free(mdl->pszRcv);
	free(mdl->pszFlsh);
	for (i=0;i<mdl->nThreads;++i) free(mdl->ppszRpl[i]);
	free(mdl->ppszRpl);
	free(mdl->pmidRpl);
	free(mdl->pReqRpl);
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
int mdlSwap(MDL mdl,int id,size_t nBufBytes,void *vBuf,size_t nOutBytes,
			size_t *pnSndBytes,size_t *pnRcvBytes)
{
	size_t nInBytes,nOutBufBytes;
	int nInMax,nOutMax,i;
	int nBytes,iTag,pid;
	char *pszBuf = vBuf;
	char *pszIn,*pszOut;
	struct swapInit {
		size_t nOutBytes;
		size_t nBufBytes;
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
	MPI_Isend(&swi,sizeof(swi),MPI_BYTE,id,MDL_TAG_SWAPINIT,
		  MPI_COMM_WORLD, &request);
	/*
	 ** Receive the number of target thread rejects and target free space
	 */
	iTag = MDL_TAG_SWAPINIT;
	pid = id;
	MPI_Recv(&swo,sizeof(swo),MPI_BYTE,pid,iTag,MPI_COMM_WORLD, &status);
	MPI_Get_count(&status, MPI_BYTE, &nBytes);
	assert(nBytes == sizeof(swo));
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
		MPI_Isend(mdl->pszTrans,nOutMax,MPI_BYTE,id,MDL_TAG_SWAP,
			 MPI_COMM_WORLD, &request);
		iTag = MDL_TAG_SWAP;
		pid = id;
		MPI_Recv(pszIn,nInMax,MPI_BYTE,pid,iTag,MPI_COMM_WORLD,
			 &status);
		MPI_Get_count(&status, MPI_BYTE, &nBytes);
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
	 ** Note: use of Ssend is mandatory here, also because of this we
	 ** don't need to use the intermediate buffer mdl->pszTrans.
	 */
	while (nOutBytes && nOutBufBytes) {
		nOutMax = (nOutBytes < MDL_TRANS_SIZE)?nOutBytes:MDL_TRANS_SIZE;
		nOutMax = (nOutMax < nOutBufBytes)?nOutMax:nOutBufBytes;
		MPI_Ssend(pszOut,nOutMax,MPI_BYTE,id,MDL_TAG_SWAP,
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
		MPI_Recv(pszIn,nInMax,MPI_BYTE,id,iTag,MPI_COMM_WORLD,
			 &status);
		MPI_Get_count(&status, MPI_BYTE, &nBytes);
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
	/*
	 ** If this looks like dangerous magic, it's because it is!
	 */
	SRVHEAD *ph = (SRVHEAD *)mdl->pszBuf;
	char *pszOut = &mdl->pszBuf[sizeof(SRVHEAD)];
	int i;

	ph->idFrom = mdl->idSelf;
	ph->sid = sid;
	if (!pszIn) ph->nInBytes = 0;
	else ph->nInBytes = nInBytes;
	if (nInBytes > 0 && pszIn != NULL) {
		for (i=0;i<nInBytes;++i) pszOut[i] = pszIn[i];
		}
	MPI_Send(mdl->pszBuf,nInBytes+sizeof(SRVHEAD),MPI_BYTE,id,MDL_TAG_REQ,
		 MPI_COMM_WORLD);
	}


void mdlGetReply(MDL mdl,int id,void *vout,int *pnOutBytes)
{
	char *pszOut = vout;
	SRVHEAD *ph = (SRVHEAD *)mdl->pszBuf;
	char *pszIn = &mdl->pszBuf[sizeof(SRVHEAD)];
	int i,iTag,nBytes;
	MPI_Status status;

	iTag = MDL_TAG_RPL;
	MPI_Recv(mdl->pszBuf,mdl->nMaxSrvBytes+sizeof(SRVHEAD),MPI_BYTE,
					id,iTag,MPI_COMM_WORLD, &status);
	MPI_Get_count(&status, MPI_BYTE, &nBytes);
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
	int sid,iTag,id,nOutBytes,nBytes;
	MPI_Status status;

	sid = 1;
	while (sid != SRV_STOP) {
		iTag = MDL_TAG_REQ;
		id = MPI_ANY_SOURCE;
		MPI_Recv(mdl->pszIn,mdl->nMaxSrvBytes+sizeof(SRVHEAD),
			       MPI_BYTE, id,iTag,MPI_COMM_WORLD,&status);
		/*
		 ** Quite a few sanity checks follow.
		 */
		id = status.MPI_SOURCE;
		MPI_Get_count(&status, MPI_BYTE, &nBytes);
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
			 MPI_BYTE, id,MDL_TAG_RPL, MPI_COMM_WORLD);
		}
	}

#define MDL_TAG_CACHECOM	10
#define MDL_MID_CACHEIN		1
#define MDL_MID_CACHEREQ	2
#define MDL_MID_CACHERPL	3
#define MDL_MID_CACHEOUT	4
#define MDL_MID_CACHEFLSH	5
#define MDL_MID_CACHEDONE	6

#define MDL_CHECK_MASK  	0x7f
#define BILLION				1000000000

/*
 * Calling mdlCacheReceive presupposes that you already have posted an MPI_Irecv 
 * with handle mdl->ReqRcv, data going into the buffer mdl->pszRcv (which is
 * of size mdl->iCaBufSize), and tag MDL_TAG_CACHECOM.
 * Initially, mdl->pszRcv is just the size of CAHEAD.  Once a cache is initialized,
 * however, it is realloced (in AdjustDataSize) to the size of 
 * CAHEAD + the size of an MDL cache line. 
 */

int mdlCacheReceive(MDL mdl,char *pLine)
{
	CACHE *c;
	CAHEAD *ph = (CAHEAD *)mdl->pszRcv;
	char *pszRcv = &mdl->pszRcv[sizeof(CAHEAD)];
	CAHEAD *phRpl;
	char *pszRpl;
	char *t;
	int id, iTag;
	int n,i;
	int ret;
	int iLineSize;
	int iDataSize;
	MPI_Status status;
#if (0)
	char achDiag[256];
#endif

	ret = MPI_Wait(&mdl->ReqRcv, &status);
	assert(ret == MPI_SUCCESS);
#if 0
	sprintf(achDiag, "%d: cache %d, message %d, from %d, rec top\n",
		mdl->idSelf, ph->cid, ph->mid, ph->id);
	mdlDiag(mdl, achDiag);
#endif

	c = &mdl->cache[ph->cid];
	assert(c->iType != MDL_NOCACHE);
	
	switch (ph->mid) {
	case MDL_MID_CACHEIN:
		++c->nCheckIn;
		ret = 0;
		break;
	case MDL_MID_CACHEOUT:
		++c->nCheckOut;
		ret = 0;
		break;
	case MDL_MID_CACHEREQ:
        assert(c->iType != MDL_DUMCACHE);
		/*
		 ** This is the tricky part! Here is where the real deadlock
		 ** difficulties surface. Making sure to have one buffer per
		 ** thread solves those problems here.
		 */
		pszRpl = &mdl->ppszRpl[ph->id][sizeof(CAHEAD)];
		phRpl = (CAHEAD *)mdl->ppszRpl[ph->id];
		phRpl->cid = ph->cid;
		phRpl->mid = MDL_MID_CACHERPL;
		phRpl->id = mdl->idSelf;
		t = &c->pData[ph->iLine*c->iLineSize];
		if(t+c->iLineSize > c->pData + c->nData*c->iDataSize)
			iLineSize = c->pData + c->nData*c->iDataSize - t;
		else
			iLineSize = c->iLineSize;
		for (i=0;i<iLineSize;++i) pszRpl[i] = t[i];
		if(mdl->pmidRpl[ph->id] != -1) {
			MPI_Wait(&mdl->pReqRpl[ph->id], &status);
		        }
		mdl->pmidRpl[ph->id] = 0;
		MPI_Isend(phRpl,sizeof(CAHEAD)+iLineSize,MPI_BYTE,
			 ph->id, MDL_TAG_CACHECOM, MPI_COMM_WORLD,
			  &mdl->pReqRpl[ph->id]); 
		ret = 0;
		break;
	case MDL_MID_CACHEFLSH:
		assert(c->iType == MDL_COCACHE);
		i = ph->iLine*MDL_CACHELINE_ELTS;
		t = &c->pData[i*c->iDataSize];
		/*
		 ** Make sure we don't combine beyond the number of data elements!
		 */
		n = i + MDL_CACHELINE_ELTS;
		if (n > c->nData) n = c->nData;
		n -= i;
		n *= c->iDataSize;
		iDataSize = c->iDataSize;
		for (i=0;i<n;i+=iDataSize) {
			(*c->combine)(&t[i],&pszRcv[i]);
			}
		ret = 0;
		break;
	case MDL_MID_CACHERPL:
        assert(c->iType != MDL_DUMCACHE);
		/*
		 ** For now assume no prefetching!
		 ** This means that this WILL be the reply to this Aquire
		 ** request.
		 */
		assert(pLine != NULL);
		iLineSize = c->iLineSize;
		for (i=0;i<iLineSize;++i) pLine[i] = pszRcv[i];
		if (c->iType == MDL_COCACHE && c->init) {
			/*
			 ** Call the initializer function for all elements in 
			 ** the cache line.
			 */
			for (i=0;i<c->iLineSize;i+=c->iDataSize) {
				(*c->init)(&pLine[i]);
				}
			}
		ret = 1;
		break;
	case MDL_MID_CACHEDONE:
	      /*
	       * No more caches, shouldn't get here.
	       */
		assert(0);
		break;
	default:
		assert(0);
		}

#if 0
	sprintf(achDiag, "%d: cache %d, message %d rec bottom\n", mdl->idSelf,
		ph->cid, ph->mid);
	mdlDiag(mdl, achDiag);
#endif
	/*
	 * Fire up next receive
	 */
	id = MPI_ANY_SOURCE;
	iTag = MDL_TAG_CACHECOM;
	MPI_Irecv(mdl->pszRcv,mdl->iCaBufSize, MPI_BYTE, id,
		 iTag, MPI_COMM_WORLD, &mdl->ReqRcv);

	return ret;
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
		MPI_Status status;
		CAHEAD caOut;

		/* cancel outstanding receive by sending a message to
		   myself */

		caOut.cid = 0;
		caOut.mid = MDL_MID_CACHEDONE;
		caOut.id = mdl->idSelf;
		MPI_Send(&caOut,sizeof(CAHEAD),MPI_BYTE, mdl->idSelf,
			 MDL_TAG_CACHECOM, MPI_COMM_WORLD);
		MPI_Wait(&mdl->ReqRcv, &status);

		mdl->iMaxDataSize = iMaxDataSize;
		mdl->iCaBufSize = sizeof(CAHEAD) + 
			iMaxDataSize*(1 << MDL_CACHELINE_BITS);
		mdl->pszRcv = realloc(mdl->pszRcv,mdl->iCaBufSize);
		assert(mdl->pszRcv != NULL);
		for (i=0;i<mdl->nThreads;++i) {
			mdl->ppszRpl[i] = realloc(mdl->ppszRpl[i],mdl->iCaBufSize);
			assert(mdl->ppszRpl[i] != NULL);
			}
		mdl->pszFlsh = realloc(mdl->pszFlsh,mdl->iCaBufSize);
		assert(mdl->pszFlsh != NULL);
		
		/*
		 * Fire up receive again.
		 */
		MPI_Irecv(mdl->pszRcv,mdl->iCaBufSize, MPI_BYTE,
			  MPI_ANY_SOURCE, MDL_TAG_CACHECOM,
			  MPI_COMM_WORLD, &mdl->ReqRcv);
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
     if (iSize > 0)
	  return(malloc(iSize));
     return NULL;
}

/* mdlMallocs and returns the max size malloced on all threads.*/
void *mdlMallocMax(MDL mdl, size_t iSize, int *iMaxSize)
{
     int *iSizes, iMax, i;

     iSizes = (int *)malloc(mdl->nThreads*sizeof(int));
     MPI_Allgather(&iSize, 1, MPI_INT, iSizes, 1, MPI_INT, MPI_COMM_WORLD);
     assert(iSizes[mdl->idSelf] == iSize);
     iMax = iSize;
     for (i=0;i<mdl->nThreads;++i) if (iSizes[i] > iMax) iMax = iSizes[i];
     *iMaxSize = iMax;
     free(iSizes);
     if (iSize > 0)
	  return(malloc(iSize));
     return NULL;
}

/* 
** Requests an mdlMalloc of an array that is at least nThreads*iSize long.
** The *actual* size of the allocation on every PE is the largest iSize called
** by any PE.  This is returned in iEltSize.
*/
void *mdlMallocShared(MDL mdl, size_t stSize, int *iEltSize)
{
     int *iSizes, iSize, iMaxSize, i;

     iSize = (int)stSize;
     iSizes = (int *)malloc(mdl->nThreads*sizeof(int));
     MPI_Allgather(&iSize, 1, MPI_INT, iSizes, 1, MPI_INT, MPI_COMM_WORLD);
     assert(iSizes[mdl->idSelf] == iSize);
     iMaxSize = iSize;
     for (i=0;i<mdl->nThreads;++i) if (iSizes[i] > iMaxSize) iMaxSize = iSizes[i];
     free(iSizes);
     *iEltSize = iMaxSize;
     if (iSize > 0)
	  return(malloc(mdl->nThreads*iMaxSize));
     return NULL;
}

void mdlFree(MDL mdl,void *p)
{
	free(p);
	}


int mdlComputeEltSizeFromType(MDL mdl, int iType)
{
    /* Figure out number of bytes to copy */
    switch (iType) {
    case MDL_TYPE_INT:
        return sizeof(int);
    case MDL_TYPE_LONG:
        return sizeof(long);
    case MDL_TYPE_SHORT:
        return sizeof(short);
    case MDL_TYPE_UNSIGNED_SHORT:
        return sizeof(unsigned short);
    case MDL_TYPE_UNSIGNED:
        return sizeof(unsigned);
    case MDL_TYPE_UNSIGNED_LONG:
        return sizeof(unsigned long);
    case MDL_TYPE_FLOAT:
        return sizeof(float);
    case MDL_TYPE_DOUBLE:
        return sizeof(double);
    case MDL_TYPE_LONG_DOUBLE:
        return sizeof(long double);
    case MDL_TYPE_BYTE:
        return sizeof(char);
    default:
        assert(1);
    }
    return -1;
}


/* Returns the MPI_Datatype of an MDL datatype */
MPI_Datatype mdlMpiDatatypeFromType(MDL mdl, int iType)
{
    switch (iType) {
    case MDL_TYPE_INT:
        return MPI_INT;
    case MDL_TYPE_LONG:
        return MPI_LONG;
    case MDL_TYPE_SHORT:
        return MPI_SHORT;
    case MDL_TYPE_UNSIGNED_SHORT:
        return MPI_UNSIGNED_SHORT;
    case MDL_TYPE_UNSIGNED:
        return MPI_UNSIGNED;
    case MDL_TYPE_UNSIGNED_LONG:
        return MPI_UNSIGNED_LONG;
    case MDL_TYPE_FLOAT:
        return MPI_FLOAT;
    case MDL_TYPE_DOUBLE:
        return MPI_DOUBLE;
    case MDL_TYPE_LONG_DOUBLE:
        return MPI_LONG_DOUBLE;
    case MDL_TYPE_BYTE:
        return MPI_BYTE;
    default:
        assert(1);
    }
    return MPI_BYTE;
}

MPI_Op mdlMpiOpFromReduce(MDL mdl, int iReduce)
{
    switch (iReduce) {
    case MDL_REDUCE_MAX:
        return MPI_MAX;
    case MDL_REDUCE_MIN:
        return MPI_MIN;
    case MDL_REDUCE_SUM:
        return MPI_SUM;
    case MDL_REDUCE_PROD:
        return MPI_PROD;
    case MDL_REDUCE_LAND:
        return MPI_LAND;
    case MDL_REDUCE_BAND:
        return MPI_BAND;
    case MDL_REDUCE_LOR:
        return MPI_LOR;
    case MDL_REDUCE_LXOR:
        return MPI_LXOR;
    case MDL_REDUCE_MAXLOC:
        return MPI_MAXLOC;
    case MDL_REDUCE_MINLOC:
        return MPI_MINLOC;
    default:
        assert(1);
    }
    return MPI_SUM;
}

/*
** Does an MPI_Allgather on array *parray.  iEltSize (which is
** also is the iEltSize that was returned from mdlMallocShared) is
** the size of a single thread's space in bytes.
*/
void mdlCollectShared(MDL mdl, void *parray, int iEltSize)
{
     int *iSizes, i;
     char *array = parray;

     /* Make sure that all iEltSizes are the same on all threads */
     iSizes = (int *)malloc(mdl->nThreads*sizeof(int));
     assert(iSizes != NULL);
     MPI_Allgather(&iEltSize, 1, MPI_INT, iSizes, 1, MPI_INT, MPI_COMM_WORLD);
     for (i=0;i<mdl->nThreads;++i) mdlassert(mdl,iSizes[i]==iEltSize);
     free(iSizes);

     /* Now gather the actual array */
     /* Note to self: I am not 100% positive that the source array can
      * be the same as this PE's piece of the destination array, even though
      * in principle MPI should be smart enough not to touch it. */
     MPI_Allgather(&(array[(mdl->idSelf)*iEltSize]), iEltSize, MPI_BYTE, array,
		   iEltSize, MPI_BYTE, MPI_COMM_WORLD);
     return;
}


void mdlBroadcast(MDL mdl, int iRoot, void *vBuf, int iEltSize)
{
    MPI_Bcast(vBuf, iEltSize, MPI_BYTE, iRoot, MPI_COMM_WORLD);
    return;
}


/*
** Does an MPI_Gather the vOutElts on all threads to vInArray on
** the master thread.  vInArray should be NULL on all other threads.
*/
void mdlGather(MDL mdl, void *vOutElt, void *vInArray, int iEltSize)
{
     int *iSizes, i;
     char *inArray = vInArray;
     char *outElt = vOutElt;

     if (mdl->idSelf == 0) { 
         mdlassert(mdl, vInArray != NULL);
     } else {
         mdlassert(mdl, vInArray == NULL);
     }

     /* Make sure that all iEltSizes are the same on all threads */
     iSizes = (int *)malloc(mdl->nThreads*sizeof(int));
     assert(iSizes != NULL);
     MPI_Allgather(&iEltSize, 1, MPI_INT, iSizes, 1, MPI_INT, MPI_COMM_WORLD);
     for (i=0;i<mdl->nThreads;++i) mdlassert(mdl,iSizes[i]==iEltSize);
     free(iSizes);

     /* Now gather the actual elements into the array */
     MPI_Gather(outElt, iEltSize, MPI_BYTE, inArray,
		   iEltSize, MPI_BYTE, 0, MPI_COMM_WORLD);
     return;
}


/*
** Does an MPI_Gatherv the vOutElts on all threads to an array (returned as void *) on
** the master thread.  On the master PE, iEltInSizes should be an array of size nThreads
** and will be filled with the sizes of each threads output.  iEltOutSize is the size
** of this thread's outgoing data.  piInArrayStride is 0 for all slave threads and 
** is the stride of the data in vInArray on the master thread.
** The return value should be NULL on all slave threads.
*/
void *mdlGatherv(MDL mdl, void *vOutElt, int *iEltInSizes, int iEltOutSize,
                 int *piInArrayStride)
{
     int i, *iStrides, iStride = 0;
     char *inArray;
     char *outElt = vOutElt;

     if (mdl->idSelf == 0) { 
         mdlassert(mdl, iEltInSizes != NULL);
     } else {
         mdlassert(mdl, iEltInSizes == NULL);
     }

     /* Fill iEltInSizes with the values of iEltSizes for all threads */
     MPI_Gather(&iEltOutSize, 1, MPI_INT, iEltInSizes, 1, MPI_INT, 0, MPI_COMM_WORLD);
     if (mdl->idSelf == 0) {
         /* Find the proper stride*/
         for (i=0;i<mdl->nThreads;++i) if (iEltInSizes[i] > iStride) iStride = iEltInSizes[i];
         mdlprintf(mdl,"mdlGatherv found optimal stride of %d\n", iStride);
         /* Malloc and fill array with strides for MPI_Gatherv call */
         iStrides = (int *) malloc(mdl->nThreads * sizeof(int));
         assert(iStrides != NULL);
         for (i=0;i<mdl->nThreads;++i) iStrides[i] = i*iStride;
         /* Malloc inArray */
         inArray = (char *) malloc(mdl->nThreads * iStride);
         assert(inArray != NULL);
     } else {
         iStride = 0;
         iStrides = NULL;
         inArray = NULL;
     }         
     *piInArrayStride = iStride;
         

     /* Now gather the actual elements into the array */
     MPI_Gatherv(outElt, iEltOutSize, MPI_BYTE, inArray,
                 iEltInSizes, iStrides, MPI_BYTE, 0, MPI_COMM_WORLD);

     if (mdl->idSelf == 0) free(iStrides);
     return (void *) inArray;
}

void mdlScatter(MDL mdl, void *vOutArray, void *vInElt, int iEltSize)
{
     int *iSizes, i;
     char *outArray = vOutArray;
     char *inElt = vInElt;

     if (mdl->idSelf == 0) { 
         mdlassert(mdl, vOutArray != NULL);
     } else {
         mdlassert(mdl, vOutArray == NULL);
     }

     /* Make sure that all iEltSizes are the same on all threads */
     iSizes = (int *)malloc(mdl->nThreads*sizeof(int));
     assert(iSizes != NULL);
     MPI_Allgather(&iEltSize, 1, MPI_INT, iSizes, 1, MPI_INT, MPI_COMM_WORLD);
     for (i=0;i<mdl->nThreads;++i) mdlassert(mdl,iSizes[i]==iEltSize);
     free(iSizes);

     /* Now gather the actual elements into the array */
     MPI_Scatter(outArray, iEltSize, MPI_BYTE, inElt,
		   iEltSize, MPI_BYTE, 0, MPI_COMM_WORLD);
     return;
}


void mdlAllReduce(MDL mdl, int iType, int iReduce, void *pSendArray, 
                  void *pReceiveArray, int iEltSize, int nElements)
{
    int *iSizes, i;
    
    mdlassert(mdl, iEltSize == mdlComputeEltSizeFromType(mdl, iType));
    /* Make sure that nElements are the same on all threads */
    iSizes = (int *)malloc(mdl->nThreads*sizeof(int));
    assert(iSizes != NULL);
    MPI_Allgather(&nElements, 1, MPI_INT, iSizes, 1, MPI_INT, MPI_COMM_WORLD);
    for (i=0;i<mdl->nThreads;++i) mdlassert(mdl,iSizes[i]==nElements);
    free(iSizes);

    /* Now do the reduction */
    MPI_Allreduce(pSendArray, pReceiveArray, nElements, mdlMpiDatatypeFromType(mdl, iType),
                  mdlMpiOpFromReduce(mdl, iReduce), MPI_COMM_WORLD);
}


/* Similar to mdlCollectShared except that it does an MPI_Alltoall from
 * *poutarray into *pinarray */
void mdlAllToAll(MDL mdl, void *pScatterArray, void *pGatherArray, int iEltSize)
{
     int *iSizes, i;
     char *outarray = pScatterArray;
     char *inarray = pGatherArray;

     /* Make sure that all iEltSizes are the same on all threads */
     iSizes = (int *)malloc(mdl->nThreads*sizeof(int));
     assert(iSizes != NULL);
     MPI_Allgather(&iEltSize, 1, MPI_INT, iSizes, 1, MPI_INT, MPI_COMM_WORLD);
     for (i=0;i<mdl->nThreads;++i) mdlassert(mdl,iSizes[i]==iEltSize);
     free(iSizes);

     /* Now to the all-to-all */
     MPI_Alltoall(outarray, iEltSize, MPI_BYTE, inarray, iEltSize, MPI_BYTE,
                  MPI_COMM_WORLD);
}



/*
 ** Common initialization for all types of caches.
 */
CACHE *CacheInitialize(MDL mdl,int cid,void *pData,int iDataSize,int nData)
{
	CACHE *c;
	int i,nMaxCacheIds;
	int first;

	/*
	 ** Allocate more cache spaces if required!
	 */
	assert(cid >= 0);
	/*
	 * first cache?
	 */
	first = 1;
	for(i = 0; i < mdl->nMaxCacheIds; ++i) {
	    if(mdl->cache[i].iType != MDL_NOCACHE) {
		first = 0;
		break;
		}
	    }
	if(first) {
	    /*
	     * Fire up first receive
	     */
	    MPI_Irecv(mdl->pszRcv,mdl->iCaBufSize, MPI_BYTE, MPI_ANY_SOURCE,
		      MDL_TAG_CACHECOM, MPI_COMM_WORLD, &mdl->ReqRcv);
	    }
	
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
    /* If this is a "dummy" cache, skip to the end */
    if (nData == -1) goto InitDummyCache;
        
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
	c->iLine = 1;
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
    /* Init cache statistics counters */
	c->nAccess = 0;
	c->nAccHigh = 0;
	c->nMiss = 0;				/* !!!, not NB */
	c->nColl = 0;				/* !!!, not NB */
	c->nMin = 0;				/* !!!, not NB */	
	c->nKeyMax = 500;				/* !!!, not NB */
	c->pbKey = malloc(c->nKeyMax);			/* !!!, not NB */
	assert(c->pbKey != NULL);			/* !!!, not NB */
	for (i=0;i<c->nKeyMax;++i) c->pbKey[i] = 0;	/* !!!, not NB */
    /* Init cache wait timers */
    c->dTimerWaitReplace = 0.;
    c->dTimerWaitFlush = 0.;
	/*
	 ** Allocate cache data lines.
	 */
	c->pLine = malloc(c->nLines*c->iLineSize);
	assert(c->pLine != NULL);

 InitDummyCache:
	c->nCheckOut = 0;
	/*
	 ** Set up the request message as much as possible!
	 */
	c->caReq.cid = cid;
	c->caReq.mid = MDL_MID_CACHEREQ;
	c->caReq.id = mdl->idSelf;
	return(c);
	}

/*
 ** Initialize a Read-Only caching space.
 */
void mdlROcache(MDL mdl,int cid,void *pData,int iDataSize,int nData)
{
	CACHE *c;
	int id;
	CAHEAD caIn;
	char achDiag[256];

	c = CacheInitialize(mdl,cid,pData,iDataSize,nData);
	c->iType = MDL_ROCACHE;
	/*
	 ** For an ROcache these two functions are not needed.
	 */
	c->init = NULL;
	c->combine = NULL;
	sprintf(achDiag, "%d: before CI, cache %d\n", mdl->idSelf, cid);
	mdlDiag(mdl, achDiag);
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
		MPI_Send(&caIn,sizeof(CAHEAD),MPI_BYTE, 0,
			       MDL_TAG_CACHECOM, MPI_COMM_WORLD);
		}
	sprintf(achDiag, "%d: In CI, cache %d\n", mdl->idSelf, cid);
	mdlDiag(mdl, achDiag);
	if(mdl->idSelf == 0) {
	    for(id = 1; id < mdl->nThreads; id++) {
		MPI_Send(&caIn,sizeof(CAHEAD),MPI_BYTE, id,
			       MDL_TAG_CACHECOM, MPI_COMM_WORLD);
		}
	    }
	else {
	    c->nCheckIn = 0;
	    while (c->nCheckIn == 0) {
		  mdlCacheReceive(mdl,NULL);
		}	
	    }	
	sprintf(achDiag, "%d: After CI, cache %d\n", mdl->idSelf, cid);
	mdlDiag(mdl, achDiag);
	AdjustDataSize(mdl);
	MPI_Barrier(MPI_COMM_WORLD);
	}

/*
 ** Initialize a Combiner caching space.
 */
void mdlCOcache(MDL mdl,int cid,void *pData,int iDataSize,int nData,
				void (*init)(void *),void (*combine)(void *,void *))
{
	CACHE *c;
	int id;
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
		MPI_Send(&caIn,sizeof(CAHEAD),MPI_BYTE, 0,
			       MDL_TAG_CACHECOM, MPI_COMM_WORLD);
		}
	if(mdl->idSelf == 0) {
	    for(id = 1; id < mdl->nThreads; id++) {
		MPI_Send(&caIn,sizeof(CAHEAD),MPI_BYTE, id,
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
	MPI_Barrier(MPI_COMM_WORLD);
	}

/*
 ** Initialize a "Dummy" caching space used only for synchronizes.
 */
void mdlDUMcache(MDL mdl,int cid)
{
	CACHE *c;
	int id;
	CAHEAD caIn;
	char achDiag[256];
	MPI_Status status;
	int ret;
    /* Set these to "dummy" values */
    void *pData = NULL;
    int iDataSize = 0;
    int nData = -1;

	c = CacheInitialize(mdl,cid,pData,iDataSize,nData);
	c->iType = MDL_DUMCACHE;
	/*
	 ** For a DUMcache these two functions are not needed.
	 */
	c->init = NULL;
	c->combine = NULL;
	sprintf(achDiag, "%d: before CI, cache %d\n", mdl->idSelf, cid);
	mdlDiag(mdl, achDiag);
	/*
	 ** THIS IS A SYNCHRONIZE!!!
	 */
	caIn.cid = cid;
	caIn.mid = MDL_MID_CACHEIN;
	caIn.id = mdl->idSelf;
	if(mdl->idSelf == 0) {
	    c->nCheckIn = 1;
	    while(c->nCheckIn < mdl->nThreads) {
		  ret = MPI_Wait(&mdl->ReqRcv, &status);
		  assert(ret == MPI_SUCCESS);
		  mdlCacheReceive(mdl, NULL);
	    }
	}
	else {
		/*
		 ** Must use non-blocking sends here, we will never wait
		 ** for these sends to complete, but will know for sure
		 ** that they have completed.
		 */
		MPI_Send(&caIn,sizeof(CAHEAD),MPI_BYTE, 0,
			       MDL_TAG_CACHECOM, MPI_COMM_WORLD);
		}
	sprintf(achDiag, "%d: In CI, cache %d\n", mdl->idSelf, cid);
	mdlDiag(mdl, achDiag);
	if(mdl->idSelf == 0) {
	    for(id = 1; id < mdl->nThreads; id++) {
		MPI_Send(&caIn,sizeof(CAHEAD),MPI_BYTE, id,
			       MDL_TAG_CACHECOM, MPI_COMM_WORLD);
		}
	    }
	else {
	    c->nCheckIn = 0;
	    while (c->nCheckIn == 0) {
		  ret = MPI_Wait(&mdl->ReqRcv, &status);
		  assert(ret == MPI_SUCCESS);
		  mdlCacheReceive(mdl,NULL);
		}	
	    }	
	sprintf(achDiag, "%d: After CI, cache %d\n", mdl->idSelf, cid);
	mdlDiag(mdl, achDiag);
	AdjustDataSize(mdl);
	MPI_Barrier(MPI_COMM_WORLD);
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
	int last;
	MPI_Status status;
	MPI_Request reqFlsh;
	MPI_Request reqBoth[2];
	int index;

	/* Now we can call mdlFinishCache even if we don't know
	 * if the cache is actually active. */
	if (c->iType == MDL_NOCACHE) return;

	if (c->iType == MDL_COCACHE) {
		/*
		 * Extra checkout to let everybody finish before
		 * flushes start.
		* I think this makes for bad synchronizes --trq
		caOut.cid = cid;
		caOut.mid = MDL_MID_CACHEOUT;
		caOut.id = mdl->idSelf;
		for(id = 0; id < mdl->nThreads; id++) {
		    if(id == mdl->idSelf)
			continue;
		    MPI_Send(&caOut,sizeof(CAHEAD),MPI_BYTE, id,
			     MDL_TAG_CACHECOM, MPI_COMM_WORLD);
		    }
		++c->nCheckOut;
		while(c->nCheckOut < mdl->nThreads) {
		  mdlCacheReceive(mdl, NULL);
		}
		c->nCheckOut = 0;
		 */
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
				/*
				 * Use Synchronous send so as not to
				 * overwhelm the receiver.
				 */
				MPI_Issend(caFlsh, sizeof(CAHEAD)+c->iLineSize,
					 MPI_BYTE, id, MDL_TAG_CACHECOM,
					 MPI_COMM_WORLD, &reqFlsh); 
				/*
				 * Wait for the Flush to complete, but
				 * also service any incoming cache requests.
				*/
				reqBoth[0] = mdl->ReqRcv;
				reqBoth[1] = reqFlsh;
				
				while(1) {
				    MPI_Waitany(2, reqBoth, &index, &status);
				    assert(!(index != 0 && reqBoth[0] ==
					   MPI_REQUEST_NULL));
				    mdl->ReqRcv = reqBoth[0];
				    if(index == 1) /* Flush has completed */
					break;
				    else if(index == 0) {
					  mdlCacheReceive(mdl, NULL);
					  reqBoth[0] = mdl->ReqRcv;
					}
				    else
					assert(0);
				    }
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
	    MPI_Send(&caOut,sizeof(CAHEAD),MPI_BYTE, 0,
			       MDL_TAG_CACHECOM, MPI_COMM_WORLD);
	    }
	if(mdl->idSelf == 0) {
	    for(id = 1; id < mdl->nThreads; id++) {
		MPI_Send(&caOut,sizeof(CAHEAD),MPI_BYTE, id,
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
    if (c->iType != MDL_DUMCACHE) {
        free(c->pTrans);
        free(c->pTag);
        free(c->pbKey);
        free(c->pLine);
    }
	c->iType = MDL_NOCACHE;
	  
	AdjustDataSize(mdl);
	/*
	 * last cache?
	 */
	last = 1;
	for(i = 0; i < mdl->nMaxCacheIds; ++i) {
	    if(mdl->cache[i].iType != MDL_NOCACHE) {
		last = 0;
		break;
		}
	    }
	/*
	 * shut down CacheReceive.
	 * Note: I'm sending a message to myself.
	 */
	if(last) {
	    MPI_Status status;
	  
	    caOut.cid = cid;
	    caOut.mid = MDL_MID_CACHEDONE;
	    caOut.id = mdl->idSelf;
	    MPI_Send(&caOut,sizeof(CAHEAD),MPI_BYTE, mdl->idSelf,
		     MDL_TAG_CACHECOM, MPI_COMM_WORLD);
	    MPI_Wait(&mdl->ReqRcv, &status);
	    }
	MPI_Barrier(MPI_COMM_WORLD);
	}


void mdlCacheCheck(MDL mdl)
{
    int flag;
    MPI_Status status;

    while (1) {
        MPI_Test(&mdl->ReqRcv, &flag, &status);
        if(flag == 0)
            break;
        mdlCacheReceive(mdl,NULL);
    }
    /* Also check scheduler requests */
    if (mdl->work.idScheduler >= 0) {
        while (1) {
            MPI_Test(&(mdl->work.handleRcv), &flag, &status);
            if (flag == 0) break;
            mdlWorkReceive(mdl,NULL);
        }
    }
}


void *mdlAquire(MDL mdl,int cid,int iIndex,int id)
{
	CACHE *c = &mdl->cache[cid];
	char *pLine;
	int iElt,iLine,i,iKey,iKeyVic,nKeyNew;
	int idVic;
	int iVictim,*pi;
	/*char ach[80]; unused */
	CAHEAD *caFlsh;
	char *pszFlsh;
	MPI_Status status;
	MPI_Request reqFlsh;
	/*int ret; unused */
    double dStartCpuTimer;

	++c->nAccess;
	if (!(c->nAccess & MDL_CHECK_MASK))
	        mdlCacheCheck(mdl);
	/*
	 ** Is it a local request?
	 */
	if (id == mdl->idSelf) {
		return(&c->pData[iIndex*c->iDataSize]);
		}
	/*
	 ** Determine memory block key value and cache line.
	iLine = iIndex >> MDL_CACHELINE_BITS;
	iKey = iLine*mdl->nThreads + id;
	 */
	iKey = ((iIndex&MDL_INDEX_MASK) << c->iKeyShift)| id;

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
	 ** Cache Miss.
	 */
	iLine = iIndex >> MDL_CACHELINE_BITS;
	c->caReq.cid = cid;
	c->caReq.mid = MDL_MID_CACHEREQ;
	c->caReq.id = mdl->idSelf;
	c->caReq.iLine = iLine;
	MPI_Send(&c->caReq,sizeof(CAHEAD),MPI_BYTE,
		 id,MDL_TAG_CACHECOM, MPI_COMM_WORLD);
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
		mdlprintf(mdl, "MDL CACHE FAILURE: cid == %d, no unlocked lines!\n",cid);
		assert(0);
		}
GotVictim:
	iKeyVic = c->pTag[iVictim].iKey;
	/*
	 ** 'pLine' will point to the actual data line in the cache.
	 */
	pLine = &c->pLine[iVictim*c->iLineSize];
	caFlsh = NULL;
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
			MPI_Isend(caFlsh, sizeof(CAHEAD)+c->iLineSize,
				 MPI_BYTE, idVic,
				 MDL_TAG_CACHECOM, MPI_COMM_WORLD, &reqFlsh); 
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
    dStartCpuTimer = mdlCpuTimer(mdl);
	while (1) {
		if(mdlCacheReceive(mdl,pLine)) {
            c->dTimerWaitReplace += (mdlCpuTimer(mdl) - dStartCpuTimer);
			if(caFlsh) {
                dStartCpuTimer = mdlCpuTimer(mdl);
				MPI_Wait(&reqFlsh, &status);
                c->dTimerWaitFlush += (mdlCpuTimer(mdl) - dStartCpuTimer);
            }
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

/*
** MDL cache stats routines
*/

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

/* Timers */

double mdlWaitReplace(MDL mdl, int cid)
{ 
	CACHE *c = &mdl->cache[cid];
    
    return(c->dTimerWaitReplace);
}

double mdlWaitFlush(MDL mdl, int cid)
{ 
	CACHE *c = &mdl->cache[cid];
    
    return(c->dTimerWaitFlush);
}



/* 
 * Priority queue, using a binary heap with static maximum size.
 * This code was originally posted by Rene Wiermer, rwiermer@googlemail.com,
 * and then modified for use in MDL.
 *
 * The implemementation is based on Heapsort as described in
 * "Introduction to Algorithms" (Cormen, Leiserson, Rivest, 24. printing)
 */

int mdlHeap_isEmpty(PacketHeap *h) {
	return h->size==0;
}

int mdlHeap_isFull(PacketHeap *h) {
	return h->size>=MDL_PQHEAP_SIZE;
}

void heap_heapify(PacketHeap* h,int i) {
	int l,r,smallest;
	Packet tmp;
	l=2*i; /*left child*/
	r=2*i+1; /*right child*/

	if ((l < h->size)&&(h->packets[l].priority < h->packets[i].priority))
		smallest=l;
	else 
		smallest=i;
	if ((r < h->size)&&(h->packets[r].priority < h->packets[smallest].priority))
		smallest=r;
	if (smallest!=i) {
		/*exchange to maintain heap property*/
		tmp=h->packets[smallest];
		h->packets[smallest]=h->packets[i];
		h->packets[i]=tmp;
		heap_heapify(h,smallest);
	}
}


void mdlHeap_init(PacketHeap* h) {
	h->size=0;
}

void mdlHeap_addItem(PacketHeap* h,Packet packet) {
	unsigned int i,parent;	
	h->size=h->size+1;
    assert(h->size <= MDL_PQHEAP_SIZE);
	i=h->size-1;
	parent=i/2;
	/*find the correct place to insert*/
	while ((i > 0)&&(h->packets[parent].priority > packet.priority)) {
		h->packets[i]=h->packets[parent];
		i=parent;
		parent=i/2;
	}
	h->packets[i]=packet;
}

Packet mdlHeap_extractMin(PacketHeap* h) {
	Packet pmax;
	if (mdlHeap_isEmpty(h)) {
        pmax.priority = -9;
        pmax.data = -9;
        return pmax;
    }
	pmax=h->packets[0];
	h->packets[0]=h->packets[h->size-1];
	h->size=h->size-1;
	heap_heapify(h,0);
	return pmax;
}


/*
** New MDL Work management functions
*/

#define MDL_TAG_WORKCOMM     20
#define MDL_TAG_WORKRESPONSE 21
#define MDL_WORK_ERROR    -1
#define MDL_WORK_NULL     0
#define MDL_WORK_CHECKIN  1
#define MDL_WORK_CHECKOUT 2
#define MDL_WORK_REQUEST  3
#define MDL_WORK_ASSIGN   4
#define MDL_WORK_QUERY    5
#define MDL_WORK_NOWORK   6
#define MDL_WORK_UPDATE   7

#define MDL_WORK_UNITS_PER_REQ  3
#define MDL_WORK_UNITS_PREFETCH 1

#define MDL_MIN(a,b)  ((a) < (b)? (a) : (b))


/* Copies work units into buffer wDataOut, adjusts w->nLocalWorkRemaining and iNextRemoteElt.
 * Returns the number of work units actually buffered */
int mdlCreateWorkResponseBuffer(MDL mdl, char *wDataOut)
{
    int iRemoteHead; /* Element that forms the 0th element for outgoing message */
    int nWorkUnits = MDL_MIN(mdl->work.nLocalWorkRemaining, MDL_WORK_UNITS_PER_REQ);

    /* Decrement the local work remaining counter */
    mdl->work.nLocalWorkRemaining -= nWorkUnits;
    assert(mdl->work.nLocalWorkRemaining >= 0);
    /* Move iNextRemoteElt */
    mdl->work.iNextRemoteElt -= MDL_WORK_UNITS_PER_REQ;
    assert(mdl->work.iNextRemoteElt >= 0);
    iRemoteHead = mdl->work.iNextRemoteElt + 1;

    /* Copy work spec to outgoing buffer */
    memcpy(wDataOut, &(mdl->work.cWorkList[(iRemoteHead)*(mdl->work.iWorkEltSize)]),
           mdl->work.iWorkEltSize*nWorkUnits);
    return nWorkUnits;
}
    

/* Returns the ID of the thread that has a valid work unit
 * that can be assigned.  Returns -1 if there are no more
 * work units to be assigned.  If return>=0, the work specification
 * is placed in the location pWorkSpec.  whOrigReq is the WORKHEAD
 * from the message of the original work request.  This routine also
 * modifies the priority queue and keeps it up to date.*/
int mdlFindNextWork(MDL mdl, WORKHEAD *whOrigReq, char *pWorkSpec, int *pnWorkUnits)
{
    WORK *w = &(mdl->work);
    WORKHEAD *whIn, whOut;
    char *pbufRcv; /* We need a local receive buffer*/
    char *wDataIn;
    SCHEDULE *s = &(mdl->sch);
    PacketHeap *pHeap = &(s->pWorkRemainingHeap);
    Packet p;
    int iPE, iResponse, flag;
    int nWorkUnits = 0;
    MPI_Request handleRcv;
    MPI_Status status;

    assert(s->bIAmScheduler);
    assert(w->idScheduler == mdl->idSelf);

    pbufRcv = (char *) malloc(w->iBufSize);
    assert(pbufRcv != NULL);
    whIn = (WORKHEAD *) pbufRcv;
    wDataIn = (pbufRcv) + sizeof(WORKHEAD);
    
    /* This while loops finds the next valid element in the 
     * priority queue, then asks the resulting thread for work.
     * If the thread responds that it has no work remaining,
     * then we repeat.  If we make it all the way through the
     *  priority queue, return -1. */
    iResponse = MDL_WORK_NOWORK;
    while (iResponse == MDL_WORK_NOWORK) {
        /* First, we find a valid priority queue element */
        while( !mdlHeap_isEmpty(pHeap) ) {
            p = mdlHeap_extractMin(pHeap);
            /* Check if the copy in the PQ is the most recent copy, and
             * that it's more than MDL_WORK_UNITS_PER_REQ */
            if (s->iWorkRemainingList[p.data] == (-1 * p.priority) &&
                p.priority < -MDL_WORK_UNITS_PER_REQ) {
                /* This is a valid element */
                if (p.data != mdl->idSelf) /* On another thread */
                    goto FoundValidWork;
                else /* On the scheduler thread */
                    goto FoundValidLocalWork;
            }
        }
        /* If it got here, priority queue is empty */
        iPE = -1;
        goto FinishFindNextWork;
        
    FoundValidWork:
        /* Now we need to query the processor that has this work to see
         * if it still has work availible. */
        iPE = p.data; /* iPE is the thread that we think has work availible */
        whOut.rid = whOrigReq->rid;
        whOut.oid = whOrigReq->oid;
        whOut.id  = mdl->idSelf;
        whOut.nLocalWorkRemaining = w->nLocalWorkRemaining;
        whOut.iMessage = MDL_WORK_QUERY;
        mdlprintf(mdl,"Thread %d should still have work remaining...sending request.\n",
                iPE);
        /*fprintf(stderr,"Thread %d should still have work remaining...sending request.\n",
	          iPE);*/
        MPI_Send(&whOut, sizeof(WORKHEAD), MPI_BYTE, iPE, MDL_TAG_WORKCOMM,
                 MPI_COMM_WORLD);
        /* Post receive for the eventual reply */
        MPI_Irecv(pbufRcv, w->iBufSize, MPI_BYTE, iPE, MDL_TAG_WORKRESPONSE,
                  MPI_COMM_WORLD, &handleRcv);
        /*MPI_Recv(pbufRcv, w->iBufSize, MPI_BYTE, iPE, MDL_TAG_WORKRESPONSE,
   	           MPI_COMM_WORLD, &status); */
        /* Do cache receives until we actually get the reply */
        flag = 0;
        while (!flag) {
            MPI_Test(&mdl->ReqRcv, &flag, &status);
            if (flag != 0) mdlCacheReceive(mdl,NULL);
            MPI_Test(&handleRcv, &flag, &status);
        }
        iResponse = whIn->iMessage;
        assert(whIn->rid == whOrigReq->rid);
        assert(whIn->oid == whOrigReq->oid);
        /* Regardless of response, update iWorkRemainingList */
        s->iWorkRemainingList[iPE] = whIn->nLocalWorkRemaining;
        mdlprintf(mdl,"Received response %d from thread %d.\n",
                iResponse,iPE);
        /*fprintf(stderr,"Received response %d from thread %d.\n",
	          iResponse,iPE);*/
    }

    /* Thread iPE has given us some work! */
    /* Update the priority queue */
    p.data = iPE;
    p.priority = (-1) * whIn->nLocalWorkRemaining;
    mdlHeap_addItem(pHeap, p);

    /*  Copy this into the proper
     * buffer and return */
    nWorkUnits = whIn->nWorkUnits;
    assert(nWorkUnits <= MDL_WORK_UNITS_PER_REQ);
    memcpy(pWorkSpec, wDataIn, w->iWorkEltSize*nWorkUnits);
    goto FinishFindNextWork;

 FoundValidLocalWork:
    assert(p.data == mdl->idSelf);
    assert(mdl->work.nLocalWorkRemaining >= MDL_WORK_UNITS_PER_REQ);
    iPE = mdl->idSelf;
    /* Copy to output buffer */
    nWorkUnits = mdlCreateWorkResponseBuffer(mdl, pWorkSpec);
    /* Update workremaininglist and priority queue */
    s->iWorkRemainingList[iPE] = w->nLocalWorkRemaining;
    /* If there is still work left, add to the priority queue. */
    if (mdl->work.nLocalWorkRemaining) {
      assert(mdl->work.iNextLocalElt <= mdl->work.iNextRemoteElt);
      p.data = iPE;
      p.priority = (-1) * w->nLocalWorkRemaining;
      mdlHeap_addItem(pHeap, p);
    }
    /*fprintf(stderr,"%d: mdlFindNextWork: iNextLocalElt=%d  iNextRemoteElt=%d, nLocalWorkRemaining: %d\n",
    	    mdl->idSelf, mdl->work.iNextLocalElt, mdl->work.iNextRemoteElt,
    	    mdl->work.nLocalWorkRemaining); */
    

 FinishFindNextWork:
    free(pbufRcv);
    *pnWorkUnits = nWorkUnits;
    return iPE;
}


/*
** Like mdlCacheReceive but for work management messages.  This
** routine assumes that there has been an MPI_Irecv posted for handle
** mdl->work.handleRcv with receive buffer mdl->work->pbufRcv,
** and that the presence of a message has been verified (by
** MPI_Wait or similar).
*/
int mdlWorkReceive(MDL mdl, char *pWork)
{
    WORK *w = &(mdl->work);
    WORKHEAD *whIn, *whOut, *whIOut;
    char *wDataIn, *wDataOut, *wIDataOut;
    int nWorkUnitsOut;
    SCHEDULE *s = &(mdl->sch);
    int ret, iPE;
    PacketHeap *pHeap = &(s->pWorkRemainingHeap);
    Packet p;
    MPI_Status status;
    int flag;

    whIn = (WORKHEAD *) w->pbufRcv;
    wDataIn = (w->pbufRcv) + sizeof(WORKHEAD);
    whOut = (WORKHEAD *) w->pbufSnd;
    wDataOut = (w->pbufSnd) + sizeof(WORKHEAD);
    whIOut = (WORKHEAD *) w->pbufISnd;
    wIDataOut = (w->pbufISnd) + sizeof(WORKHEAD);

    switch (whIn->iMessage) {
    case MDL_WORK_CHECKIN:
        assert(s->bIAmScheduler);
        ++(s->nCheckIn);
        s->iWorkRemainingList[whIn->id] = whIn->nLocalWorkRemaining;
        ret = MDL_WORK_CHECKIN;
        break;
    case MDL_WORK_CHECKOUT:
        assert(s->bIAmScheduler);
        ++(s->nCheckOut);
        assert(whIn->nLocalWorkRemaining == 0);
        ret = MDL_WORK_CHECKOUT;
        break;
    case MDL_WORK_QUERY:
        /* A query from the scheduler for work.  Respond using Isend */
        assert(!s->bIAmScheduler);
        MPI_Test(&(w->handleISnd), &flag, &status);
        if (!flag) {
            /* Previous request did not complete */
            fprintf(stderr,"%d: MPI_Isend did not complete, yet!  Waiting to respond to work query.\n",
                    mdl->idSelf);
            MPI_Wait(&(w->handleISnd), &status);
        }
        whIOut->rid = whIn->rid;
        whIOut->oid = whIn->oid;
        whIOut->id  = mdl->idSelf;
        mdlprintf(mdl,"%d: Received work query from master regarding request %d from thread %d...\n",
                  mdl->idSelf, whIn->rid, whIn->oid);
        /*fprintf(stderr,"%d: Received work query from master regarding request %d from thread %d...\n",
	          mdl->idSelf, whIn->rid, whIn->oid);*/
        /* This is the condition that signals there is no work left to do */
        if (mdl->work.iNextLocalElt+(MDL_WORK_UNITS_PER_REQ-1) > mdl->work.iNextRemoteElt) {
            assert(w->nLocalWorkRemaining < MDL_WORK_UNITS_PER_REQ);
            /* Tell scheduler there is not enough local work left */
            whIOut->nLocalWorkRemaining = w->nLocalWorkRemaining;
            whIOut->iMessage = MDL_WORK_NOWORK;
            whIOut->nWorkUnits = 0;
            mdlDiag(mdl,"Sending response NOWORK\n");
            /*fprintf(stderr,"Sending response NOWORK\n");*/
            MPI_Isend(whIOut, sizeof(WORKHEAD), MPI_BYTE, whIn->id,
                     MDL_TAG_WORKRESPONSE, MPI_COMM_WORLD, &(w->handleISnd));
        } else {
            assert(w->nLocalWorkRemaining >= MDL_WORK_UNITS_PER_REQ);
            /* Construct rest of message header */
            whIOut->nLocalWorkRemaining = w->nLocalWorkRemaining;
            whIOut->iMessage = MDL_WORK_ASSIGN;
            /* Copy work spec to outgoing buffer */
            whIOut->nWorkUnits = mdlCreateWorkResponseBuffer(mdl, wIDataOut);
            mdlprintf(mdl,"Sending response ASSIGN with %d work units\n",whIOut->nWorkUnits);
            /*fprintf(stderr,"Sending response ASSIGN\n");*/
            MPI_Isend(whIOut, w->iBufSize, MPI_BYTE, whIn->id,
                      MDL_TAG_WORKRESPONSE, MPI_COMM_WORLD, &(w->handleISnd));
        }
        ret = MDL_WORK_QUERY;
        break;
    case MDL_WORK_REQUEST:
        assert(s->bIAmScheduler);
        assert(whIn->oid == whIn->id);
        assert(whIn->nLocalWorkRemaining <= MDL_WORK_UNITS_PER_REQ);
        /* First, we update iWorkRemainingList */
        s->iWorkRemainingList[whIn->id] = 0;
        s->iWorkRemainingList[mdl->idSelf] = w->nLocalWorkRemaining;
        mdlprintf(mdl,"Received work request ID %d from thread %d\n",whIn->rid, whIn->id);
        /*fprintf(stderr,"Scheduler received work request ID %d from thread %d\n",whIn->rid, whIn->id);*/
        /* Pop the next element off of the priority queue (checking
         * to make sure it has not been updated) */
        MPI_Test(&(w->handleISnd), &flag, &status);
        if (!flag) {
            /* Previous request did not complete */
            fprintf(stderr,"%d: MPI_Isend did not complete, yet!  Waiting to find next work.\n",
                    mdl->idSelf);
            MPI_Wait(&(w->handleISnd), &status);
        }
        /* mdlFindNextWork loads up wIDataOut with nWorkUnitsOut units of for from PE "iPE".
         * If iPE<0, then there is no more work for anyone. */
        iPE = mdlFindNextWork(mdl, whIn, wIDataOut, &nWorkUnitsOut);
        mdlprintf(mdl,"mdlFindNextWork found work from thread %d for thread %d (rid=%d, nWorkUnits=%d)\n",
                iPE, whIn->id, whIn->rid, nWorkUnitsOut);
        /*fprintf(stderr,"mdlFindNextWork found work from thread %d for thread %d (rid=%d)\n",
	          iPE, whIn->id, whIn->rid);*/
        if (iPE < 0) { /* No more work left.  Tell thread to checkout. */
            whIOut->rid = whIn->rid;
            whIOut->oid = whIn->oid;
            whIOut->id  = mdl->idSelf;
            whIOut->nLocalWorkRemaining = 0;
            assert(s->nCheckOut);
            whIOut->iMessage = MDL_WORK_NOWORK;
            whIOut->nWorkUnits = 0;
            mdlprintf(mdl,"Sending NOWORK to thread %d\n",whIn->id);
            /*fprintf(stderr,"Scheduler Sending NOWORK to thread %d\n",whIn->id);*/
            MPI_Isend(whIOut, sizeof(WORKHEAD), MPI_BYTE, whIn->id,
                     MDL_TAG_WORKCOMM, MPI_COMM_WORLD, &(w->handleISnd));
        } else { /* Respond to thread with work assignment. */
            whIOut->rid = whIn->rid;
            whIOut->oid = whIn->oid;
            whIOut->id  = iPE; /* The thread whose work this actually is */
            whIOut->nWorkUnits = nWorkUnitsOut;
            whIOut->nLocalWorkRemaining = 0;
            whIOut->iMessage = MDL_WORK_ASSIGN;
            MPI_Isend(whIOut, w->iBufSize, MPI_BYTE, whIn->id,
                     MDL_TAG_WORKCOMM, MPI_COMM_WORLD, &(w->handleISnd));
        }
        ret = MDL_WORK_REQUEST;
        break;
    case MDL_WORK_UPDATE:
        assert(s->bIAmScheduler);
        /* Update iWorkRemainingList */
        s->iWorkRemainingList[whIn->id] = whIn->nLocalWorkRemaining;
        mdlprintf(mdl,"Received work update from thread %d. Work remaining: %d\n",
                  whIn->id, whIn->nLocalWorkRemaining);
        /* Update the priority queue */
        p.data = whIn->oid;
        p.priority = (-1) * whIn->nLocalWorkRemaining;
        mdlHeap_addItem(pHeap, p);
        ret = MDL_WORK_UPDATE;
        break;
    case MDL_WORK_NOWORK:
        assert(!(s->bIAmScheduler));
        w->retIncomingWork = MDL_WORK_NOWORK;
        w->idIncomingWork = -1;
        ret = MDL_WORK_NOWORK;
        break;
    case MDL_WORK_ASSIGN:
        assert(!(s->bIAmScheduler));
        w->retIncomingWork = MDL_WORK_ASSIGN;
        w->idIncomingWork = whIn->id;
        /* Copy work assignment to incoming work buffer */
        w->nbufIncomingWork = whIn->nWorkUnits;
        memcpy(w->pbufIncomingWork, wDataIn, w->iWorkEltSize * w->nbufIncomingWork);
        ret = MDL_WORK_ASSIGN;
        break;
    default:
        mdlprintf(mdl, "ERROR: received message %d from thread %d (rid=%d, oid=%d)!\n",
                  whIn->iMessage, whIn->id, whIn->rid, whIn->oid);
        mdlassert(mdl, 0);
    }

    /*
    ** Post next receive
    */
    MPI_Irecv(w->pbufRcv, w->iBufSize, MPI_BYTE, MPI_ANY_SOURCE,
              MDL_TAG_WORKCOMM, MPI_COMM_WORLD, &(w->handleRcv));

    return ret;
}

void mdlInitWorkNoDynamic(MDL mdl, void *pWorkList, int iWorkEltSize, int nWorkElts)
{
    mdl->work.cWorkList           = pWorkList;
    mdl->work.iWorkEltSize        = iWorkEltSize;
    mdl->work.nWorkElts           = nWorkElts;
    mdl->work.iNextLocalElt       = 0;
    mdl->work.iNextRemoteElt      = nWorkElts-1;
    mdl->work.nLocalWorkRemaining = nWorkElts;
    mdl->work.idScheduler         = -2;
    return;
}


void mdlInitWork(MDL mdl, void *pWorkList, int iWorkEltSize, int nWorkElts, int bDynamic)
{
    WORK *w = &(mdl->work);
    WORKHEAD wh;
    int ret, i;
    MPI_Status status;

    if (mdl->work.cWorkList != NULL || mdl->work.idScheduler != -1){
        mdlprintf(mdl, "mdlInitWork called after workspace already initialized!\n");
        assert(0);
    }

    mdl->work.bDynamic = bDynamic;
    /* If we are not actually doing dynamic load balancing, call the simple init function
     * and exit. */
    if (!bDynamic) { 
        mdlprintf(mdl, "\nWORK MANAGEMENT requested with no dynamic scheduling.  Workspace initialized for static scheduling.\n");
        mdlInitWorkNoDynamic(mdl, pWorkList, iWorkEltSize, nWorkElts);
        return;
    }

    mdl->work.cWorkList           = (char *) pWorkList;
    mdl->work.iWorkEltSize        = iWorkEltSize;
    mdl->work.nWorkElts           = nWorkElts;
    mdl->work.iNextLocalElt       = 0;
    mdl->work.iNextRemoteElt      = nWorkElts-1;
    mdl->work.nLocalWorkRemaining = nWorkElts;
    mdl->work.idScheduler = (mdl->nThreads)-1; /* Set the scheduler to be the highest
                                                * node on the theory that it probably
                                                * has the least amount of work.*/
    mdl->work.iRidCurrent = 0;
    
    mdlprintf(mdl, "\nDYNAMIC SCHEDULING Requested.  Initializing workspace.  Master scheduler is %d\n",
              w->idScheduler);

    /* Init scheduler space */
    if (mdl->idSelf == w->idScheduler) mdl->sch.bIAmScheduler = 1;
    else mdl->sch.bIAmScheduler = 0;
    mdl->sch.iWorkRemainingList = NULL;
    mdl->sch.nCheckIn = 0;
    mdl->sch.nCheckOut = 0;

    /* Allocate buffer sizes to send/receive work elements */
    w->idRemoteWork = -1;
    w->ibufRemoteWork = 0;
    w->nbufRemoteWork = 0;
    w->pbufRemoteWork = (char *) malloc(MDL_WORK_UNITS_PER_REQ*w->iWorkEltSize);
    assert(w->pbufRemoteWork != NULL);
    w->nbufIncomingWork = 0;
    w->retIncomingWork = MDL_WORK_NULL;
    w->idIncomingWork = -1;
    w->pbufIncomingWork = (char *) malloc(MDL_WORK_UNITS_PER_REQ*w->iWorkEltSize);
    assert(w->pbufIncomingWork != NULL);
    w->iBufSize = sizeof(WORKHEAD) + MDL_WORK_UNITS_PER_REQ*w->iWorkEltSize;
    w->pbufRcv = (char *) malloc(w->iBufSize);
    assert(w->pbufRcv != NULL);
    w->pbufSnd = (char *) malloc(w->iBufSize);
    assert(w->pbufSnd != NULL);
    w->pbufISnd = (char *) malloc(w->iBufSize);
    assert(w->pbufISnd != NULL);
    w->handleISnd = MPI_REQUEST_NULL;
    w->bRequestedMoreWork = 0;
    
    /* Post first work-related receive */
    MPI_Irecv(w->pbufRcv, w->iBufSize, MPI_BYTE, MPI_ANY_SOURCE,
              MDL_TAG_WORKCOMM, MPI_COMM_WORLD, &(w->handleRcv));

    /* Init first message */
    wh.rid = w->iRidCurrent; /* RequestID: 1st message */
    wh.oid = mdl->idSelf; /* Originator ID: This thread */
    wh.id = mdl->idSelf; /* Sender ID: This thread */
    wh.nWorkUnits = 0;
    wh.nLocalWorkRemaining = w->nLocalWorkRemaining;
    wh.iMessage = MDL_WORK_CHECKIN;

    /* Have all threads check in */
    if (mdl->sch.bIAmScheduler) {
        Packet p;
        PacketHeap *pHeap;
        mdlDiag(mdl, "I am the master scheduler!\n");
        /* Init work tracking data structs */
        mdl->sch.iWorkRemainingList = (int *) malloc(mdl->nThreads*sizeof(int));
        assert(mdl->sch.iWorkRemainingList != NULL);
        mdlDiag(mdl, "Waiting to receive check-in from all threads...");
        mdl->sch.nCheckIn = 1;
        mdl->sch.iWorkRemainingList[mdl->idSelf] = w->nLocalWorkRemaining;
        while (mdl->sch.nCheckIn < mdl->nThreads) {
            ret = MPI_Wait(&(w->handleRcv), &status);
            assert(ret == MPI_SUCCESS);
            ret = mdlWorkReceive(mdl, NULL);
            assert(ret != MDL_WORK_ERROR);
        }
        mdlDiag(mdl, "Received check-in from all threads.\n");
        /* Contruct the priority queue */
        pHeap = &(mdl->sch.pWorkRemainingHeap);
        mdlHeap_init(pHeap);
        for (i=0;i<mdl->nThreads;++i) {
            p.data = i; /* Thread ID */
            /* Set to negative since the smallest integer is at the top of the PQ */
            p.priority = (-1) * (mdl->sch.iWorkRemainingList[i]);
            /*fprintf(stderr,"Adding item data: %d  priority: %d\n",
	               p.data,p.priority);*/
            mdlHeap_addItem(pHeap, p);
        }
	/*
        j = 0;
        while (j< pHeap->size) {
            fprintf(stderr,"Heap[%d] data=%d  priority=%d\n", j,
                    pHeap->packets[j].data,
                    pHeap->packets[j].priority);
            ++j;
        }
	*/
        mdl->sch.nCheckOut = 1; /* For the moment, "check-out" the master immediately. */
    } else {
        /* Use blocking sends as this is a check-in anyway */
        MPI_Send(&wh, sizeof(WORKHEAD), MPI_BYTE, w->idScheduler, MDL_TAG_WORKCOMM,
                 MPI_COMM_WORLD);
    }
    mdlDiag(mdl, "Thread check-in complete.  Work management layer has been initialized.\n\n");
    
}

void mdlFinishWork(MDL mdl, void *pWorkList)
{ 
    WORK *w = &(mdl->work);
    SCHEDULE *s = &(mdl->sch);
    /*int ret;*/
    /*MPI_Status status;*/

    if (mdl->work.idScheduler == -1 ) {
        /* Work management was never initialized...just return*/
        assert(mdl->work.cWorkList == NULL);
        return;
    }
    
    /* Always check that the worklist is the same */
    assert(w->cWorkList == (char *)pWorkList);

    /* Receive buffer should have already been deallocated */
    assert(w->pbufRcv == NULL);

    if (mdl->work.bDynamic && s->bIAmScheduler) {
        s->bIAmScheduler = 0;
        /* Work schedule buffer should have already been deallocated */
        assert(s->iWorkRemainingList == NULL);
        s->nCheckIn = -1;
        s->nCheckOut = -1;
    }
    mdl->work.cWorkList = NULL;
    mdl->work.iNextLocalElt = -1;
    mdl->work.iNextRemoteElt = -1;
    mdl->work.nLocalWorkRemaining = -1;
    mdl->work.idScheduler = -1;
    mdl->work.bDynamic = 0;


    /* Cancel the pending MPI_Irecv and we are done with Work management.
     * When this work for real, the master scheduler will initiate the checkout. */
    /*
    ret = MPI_Cancel(&(w->handleRcv));
    assert(ret == MPI_SUCCESS);
    free(w->pbufRcv);
    */

    mdlDiag(mdl, "\nWORKLOAD MANAGEMENT completed.\n\n");
}

void mdlConstructWorkRequest(MDL mdl, WORKHEAD *wh)
{
    WORK *w = &(mdl->work);
    ++(w->iRidCurrent);
    wh->rid = w->iRidCurrent;
    wh->oid = mdl->idSelf;
    wh->id =  mdl->idSelf;
    wh->nWorkUnits = 0;
    wh->nLocalWorkRemaining = w->nLocalWorkRemaining;
    wh->iMessage = MDL_WORK_REQUEST;
}

void mdlSendWorkRequest(MDL mdl)
{
    WORK *w = &(mdl->work);
    WORKHEAD *whiOut = (WORKHEAD *) w->pbufISnd;
    int flag;
    MPI_Status status;

    MPI_Test(&(w->handleISnd), &flag, &status);
    if (!flag) {
        /* Previous request did not complete */
        fprintf(stderr,
                "%d: MPI_Isend did not complete, yet!  Waiting to send work request.\n",
                mdl->idSelf);
        MPI_Wait(&(w->handleISnd), &status);
    }
    mdlConstructWorkRequest(mdl, whiOut);
    MPI_Isend(whiOut, sizeof(WORKHEAD), MPI_BYTE, w->idScheduler,
              MDL_TAG_WORKCOMM, MPI_COMM_WORLD, &(w->handleISnd));
    w->bRequestedMoreWork = 1;
}

void mdlConstructWorkUpdate(MDL mdl, WORKHEAD *wh)
{
    wh->rid = -1;
    wh->oid = mdl->idSelf;
    wh->id = mdl->idSelf;
    wh->nWorkUnits = 0;
    wh->nLocalWorkRemaining = mdl->work.nLocalWorkRemaining;
    wh->iMessage = MDL_WORK_UPDATE;
}

void mdlSendWorkUpdate(MDL mdl)
{
    WORK *w = &(mdl->work);
    WORKHEAD *whiOut = (WORKHEAD *) w->pbufISnd;
    int flag;
    MPI_Status status;

    MPI_Test(&(w->handleISnd), &flag, &status);
    if (!flag) {
        /* Previous request did not complete */
        fprintf(stderr,
                "%d: MPI_Isend did not complete, yet!  Waiting to send work update.\n",
                mdl->idSelf);
        MPI_Wait(&(w->handleISnd), &status);
    }
    mdlConstructWorkUpdate(mdl, whiOut);
    MPI_Isend(whiOut, sizeof(WORKHEAD), MPI_BYTE, w->idScheduler, MDL_TAG_WORKCOMM,
              MPI_COMM_WORLD, &(w->handleISnd));
    /*MPI_Send(&whiOut, sizeof(WORKHEAD), MPI_BYTE, w->idScheduler, MDL_TAG_WORKCOMM,
               MPI_COMM_WORLD);*/
}

void mdlConstructWorkCheckout(MDL mdl, WORKHEAD *wh)
{
    WORK *w = &(mdl->work);

    assert(w->nLocalWorkRemaining == 0);
    ++(w->iRidCurrent);
    wh->rid = w->iRidCurrent;
    wh->oid = mdl->idSelf;
    wh->id = mdl->idSelf;
    wh->nWorkUnits = 0;
    wh->nLocalWorkRemaining = 0;
    wh->iMessage = MDL_WORK_CHECKOUT;
}


void *mdlRequestWorkNoDynamic(MDL mdl, void *pWorkList)
{
    /* This is the condition that signals there is no work left to do */
    if (mdl->work.iNextLocalElt > mdl->work.iNextRemoteElt) return(NULL);
    /* Otherwise, we increment iNextLocalElt to point to the next unassigned
     * element, and return the original iNextLocalElt element. */
    /*fprintf(stderr,"%d: iNextLocalElt=%d  iNextRemoteElt=%d\n", mdl->idSelf,
      mdl->work.iNextLocalElt, mdl->work.iNextRemoteElt);*/
    ++(mdl->work.iNextLocalElt);
    return( &(mdl->work.cWorkList[(mdl->work.iNextLocalElt-1)*mdl->work.iWorkEltSize]) );
}

void *mdlRequestWork(MDL mdl, void *pWorkList, int *pidHome)
{
    WORK *w = &(mdl->work);
    WORKHEAD wh, *whIn, *whiOut;
    char *wDataIn;
    SCHEDULE *s = &(mdl->sch);
    int ret, flag, bFirstTime;
    MPI_Status status;

    /* Always check that the worklist is the same */
    assert(w->cWorkList == (char *)pWorkList);

    if (!(w->bDynamic)) {
        /* Only static scheduling active */
        return mdlRequestWorkNoDynamic(mdl, pWorkList);
    }

    whIn = (WORKHEAD *) w->pbufRcv;
    wDataIn = (w->pbufRcv) + sizeof(WORKHEAD);
    whiOut = (WORKHEAD *) w->pbufISnd;
    *pidHome = -1;
    
    mdlCacheCheck(mdl);

    if (mdl->sch.bIAmScheduler) {
        assert(s->nCheckOut);
        /* This is the condition that signals there is no work left to do */
        if (mdl->work.iNextLocalElt > mdl->work.iNextRemoteElt) {
            mdlDiag(mdl,"Master scheduler is out of work.  Waiting for threads to checkout.\n");
            /* We need to keep servicing the other nodes until everyone checks out */
            while (s->nCheckOut < mdl->nThreads) {
                MPI_Test(&mdl->ReqRcv, &flag, &status);
                if (flag != 0) mdlCacheReceive(mdl,NULL);
                MPI_Test(&(w->handleRcv), &flag, &status);
                if (flag != 0) mdlWorkReceive(mdl, NULL);
            }
            /* All threads are now checked out.  Finish up. */
            mdlDiag(mdl,"All threads have checked out.\n");
            /* Print PQ and work list */
            fprintf(stderr,"iWorkRemainingList:\n");
	    /*
            for (i=0;i<mdl->nThreads;++i) 
                fprintf(stderr,"%d:  %d\n",i,s->iWorkRemainingList[i]);
            fprintf(stderr,"pWorkRemainingHeap:\n");
            while( !mdlHeap_isEmpty(&(s->pWorkRemainingHeap)) ) {
                Packet p;
                p = mdlHeap_extractMin(&(s->pWorkRemainingHeap));
                fprintf(stderr,"PE %d:  %d\n",p.data,p.priority);
	    }
	    */
            /* Reset iWorkRemainingList */
            free(s->iWorkRemainingList);
            s->iWorkRemainingList = NULL;
            /* Cancel the pending MPI_Irecv and free the receive buffer */
            ret = MPI_Cancel(&(w->handleRcv));
            assert(ret == MPI_SUCCESS);
            free(w->pbufRcv);
            free(w->pbufSnd);
            free(w->pbufISnd);
            w->pbufRcv = NULL;
            mdlDiag(mdl,"Pending received cancelled.  Buffers cleared.  Exiting.\n");
            return NULL;
        } else {
            Packet p;
            PacketHeap *pHeap;
            pHeap = &(mdl->sch.pWorkRemainingHeap);
            /* Otherwise, we increment iNextLocalElt to point to the next unassigned
             * element, and return the original iNextLocalElt element. */
            ++(mdl->work.iNextLocalElt);
            --(mdl->work.nLocalWorkRemaining);
            /*fprintf(stderr,"%d: mdlRequestWork: iNextLocalElt=%d  iNextRemoteElt=%d, nLocalWorkRemaining: %d\n",
                      mdl->idSelf, mdl->work.iNextLocalElt, mdl->work.iNextRemoteElt,
                      mdl->work.nLocalWorkRemaining);*/
            /* Update workremaininglist */
            s->iWorkRemainingList[mdl->idSelf] = w->nLocalWorkRemaining;
            /* If there is still work left, add to the priority queue. */
            if (mdl->work.nLocalWorkRemaining) {
                assert(mdl->work.iNextLocalElt <= mdl->work.iNextRemoteElt);
                p.data = mdl->idSelf;
                p.priority = (-1) * w->nLocalWorkRemaining;
                mdlHeap_addItem(pHeap, p);
            }
            *pidHome = mdl->idSelf;
            return( &(mdl->work.cWorkList[(mdl->work.iNextLocalElt-1)*(mdl->work.iWorkEltSize)]) );
        }            

    } else { /* I am not the scheduler */
        
        if (mdl->work.iNextLocalElt > mdl->work.iNextRemoteElt) {
            /* There is no local work left to do... */
            if (w->ibufRemoteWork < w->nbufRemoteWork) {
                /* ...But we have remote work already buffered */
                /* paddr points to the address to return */
                char *paddr = &(w->pbufRemoteWork[(w->ibufRemoteWork) * w->iWorkEltSize]);
                /* The amount of work remaining including the work unit pointed to by paddr */
                int nWorkRemaining = w->nbufRemoteWork - w->ibufRemoteWork ;
                
                ++(w->ibufRemoteWork);
                /* Check to see if we have to send off another work request */
                if (!w->bRequestedMoreWork && nWorkRemaining <= MDL_WORK_UNITS_PREFETCH) {
                    mdlprintf(mdl, "%d units of buffered remote work remaining...", nWorkRemaining);
                    mdlprintf(mdl, "sending work request to scheduler.\n");
                    mdlSendWorkRequest(mdl);
                }
                *pidHome = w->idRemoteWork;
                return(paddr);
            } else {
                /* ...And we need to acquire more remote work */
                assert(w->nLocalWorkRemaining == 0);
                w->ibufRemoteWork = 0;
                w->nbufRemoteWork = 0;
                mdlprintf(mdl,"Out of local work...");
                /* See if we still need to request work from other PEs */
                if (!w->bRequestedMoreWork) {
                    mdlprintf(mdl, "0 units of buffered remote work remaining...");
                    mdlprintf(mdl, "sending work request to scheduler.\n");
                    mdlSendWorkRequest(mdl);
                }

                /* Now do cache receives until we get our work assignment */
                bFirstTime = 1;
                while (w->retIncomingWork == MDL_WORK_NULL) {
                    if (bFirstTime) {
                        fprintf(stderr,"%d: Waiting for work assignment!\n", mdl->idSelf);
                        mdlprintf(mdl, "Waiting for work assignment!...");
                        bFirstTime = 0;
                    }
                    MPI_Test(&mdl->ReqRcv, &flag, &status);
                    if (flag != 0) mdlCacheReceive(mdl,NULL);
                    MPI_Test(&(w->handleRcv), &flag, &status);
                    if (flag != 0) {
                        ret = mdlWorkReceive(mdl, NULL);
                        /* These should be the only three responses possible at this point */
                        assert(ret == MDL_WORK_NOWORK || ret == MDL_WORK_ASSIGN || 
                               ret == MDL_WORK_QUERY);
                    }
                }                            
                mdlprintf(mdl, "Work assignment received.\n");
                w->bRequestedMoreWork = 0;
                /* If our outgoing work request has not completed yet, then
                 * there is something majorly wrong. */
                MPI_Test(&(w->handleISnd), &flag, &status);
                if (!flag) {
                    /* Previous request did not complete */
                    fprintf(stderr,
                            "%d: ERROR! MPI_Isend work request not complete, even though we!\n",
                            mdl->idSelf);
                    fprintf(stderr,"   have just received a work assignment!\n");
                    mdlassert(mdl, 0);
                }
                /* Get the work */
                if (w->retIncomingWork == MDL_WORK_NOWORK) {
                    /* There is no more work.  Send checkout and return NULL. */
                    mdlDiag(mdl,"No more work.  Sending checkout...");
                    mdlConstructWorkCheckout(mdl, &wh);
                    MPI_Send(&wh, sizeof(WORKHEAD), MPI_BYTE, w->idScheduler, MDL_TAG_WORKCOMM,
                             MPI_COMM_WORLD);
                    /* Cancel the pending MPI_Irecv and free the receive buffer */
                    ret = MPI_Cancel(&(w->handleRcv));
                    assert(ret == MPI_SUCCESS);
                    free(w->pbufRcv);
                    w->pbufRcv = NULL;
                    mdlDiag(mdl,"Checkout sent.\n");
                    mdlDiag(mdl,"Pending received cancelled.  Buffers cleared.  Exiting.\n");
                    return NULL;
                } else {
                    assert(w->retIncomingWork == MDL_WORK_ASSIGN);
                    w->idRemoteWork = w->idIncomingWork;
                    /* Copy from incoming work buffer to remote work buffer */
                    assert(w->nbufIncomingWork <= MDL_WORK_UNITS_PER_REQ);
                    w->nbufRemoteWork = w->nbufIncomingWork;
                    memcpy(w->pbufRemoteWork, w->pbufIncomingWork, w->iWorkEltSize * w->nbufRemoteWork);
                    w->ibufRemoteWork = 1; /* Next work unit (not counting the one we are
                                            * returning in the next line) */
                    /* Re-init the incoming work buffer */
                    w->retIncomingWork = MDL_WORK_NULL;
                    w->idIncomingWork = -1;
                    w->nbufIncomingWork = 0;
                    *pidHome = w->idRemoteWork;
                    return(w->pbufRemoteWork); /* Return the 0th element */
                }
            } /* Return remote work buffer or acquire new remote work buffer */

        } /* if no more local work */
        else { /* We still have local work */
            /* Address to of local work unit to return */
            char *paddr = &(mdl->work.cWorkList[(mdl->work.iNextLocalElt) *
                                                (mdl->work.iWorkEltSize)]);
                                                
            mdlCacheCheck(mdl); /*XXX Remove this cache check??*/
            /* Otherwise, we increment iNextLocalElt to point to the next unassigned
             * element, and return the original iNextLocalElt element. */
            mdlprintf(mdl,"%d: iNextLocalElt=%d  iNextRemoteElt=%d\n",
                    mdl->idSelf, mdl->work.iNextLocalElt, mdl->work.iNextRemoteElt);
            /*fprintf(stderr,"%d: iNextLocalElt=%d  iNextRemoteElt=%d\n",
	      mdl->idSelf, mdl->work.iNextLocalElt, mdl->work.iNextRemoteElt);*/
            ++(mdl->work.iNextLocalElt);
            --(mdl->work.nLocalWorkRemaining);
            /* Check to see if we send just and update, or we need to send a work request */
            if (!w->bRequestedMoreWork && 
                mdl->work.nLocalWorkRemaining+1 <= MDL_WORK_UNITS_PREFETCH) {
                /* Isend work request if remaining undone (as opposed to unassigned) work
                 * units == MDL_WORK_UNITS_PREFETCH */
                mdlprintf(mdl,"%d units of local work remaining...sending work request to scheduler.\n",
                          mdl->work.nLocalWorkRemaining+1);
                mdlSendWorkRequest(mdl);
            } else {
                /* Just Send update to scheduler */
                mdlDiag(mdl,"Sending update to scheduler.\n");
                mdlSendWorkUpdate(mdl);
            }
            
            *pidHome = mdl->idSelf;
            return paddr;
        }
    } /* I am not scheduler */

}
