#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <malloc.h>
#include <math.h>
#include <assert.h>
#include "smooth.h"
#include "pkd.h"
#include "smoothfcn.h"


int smInitialize(SMX *psmx,PKD pkd,int nSmooth,int bPeriodic,int bSymmetric,
				 int iSmoothType,int bSmooth)
{
	SMX smx;
	void (*initParticle)(void *);
	void (*init)(void *);
	void (*comb)(void *,void *);
	int pi;
	int nTree;

	smx = (SMX)malloc(sizeof(struct smContext));
	assert(smx != NULL);
	smx->pkd = pkd;
	smx->nSmooth = nSmooth;
	smx->bPeriodic = bPeriodic;

	switch (iSmoothType) {
	case SMX_DENSITY:
		smx->fcnSmooth = bSymmetric?DensitySym:Density;
		initParticle = initDensity;
		init = initDensity;
		comb = combDensity;
		smx->fcnPost = NULL;
		break;
#ifdef SUPERCOOL
	case SMX_MEANVEL:
		smx->fcnSmooth = bSymmetric?MeanVelSym:MeanVel;
		initParticle = initMeanVel;
		init = initMeanVel;
		comb = combMeanVel;
		smx->fcnPost = NULL;
		break;
#endif
#ifdef GASOLINE
	case SMX_HSMDIVV:
		smx->fcnSmooth = bSymmetric?HsmDivvSym:HsmDivv;
		initParticle = initHsmDivv;
		init = initHsmDivv;
		comb = combHsmDivv;
		smx->fcnPost = postHsmDivv;
		break;
	case SMX_GEOMBV:
		assert(bSymmetric != 0);
		smx->fcnSmooth = GeomBVSym;
		initParticle = initGeomBV;
		init = initGeomBV;
		comb = combGeomBV;
		smx->fcnPost = postGeomBV;
		break;
	case SMX_ETHDOTBV:
		assert(bSymmetric != 0);
		smx->fcnSmooth = EthdotBVSym;
		initParticle = initEthdotBV;
		init = initEthdotBV;
		comb = combEthdotBV;
		smx->fcnPost = NULL;
		break;		
	case SMX_ACCSPHBV:
		assert(bSymmetric != 0);
		smx->fcnSmooth = AccsphBVSym;
		initParticle = NULL;		/* Don't zero accelerations */
		init = initAccsph;
		comb = combAccsph;
		smx->fcnPost = NULL;
		break;
#endif
#ifdef PLANETS
	case SMX_TIMESTEP:
		assert(bSymmetric == 0);
		smx->fcnSmooth = SetTimeStep;
		initParticle = NULL;
		init = NULL;
		comb = NULL;
		smx->fcnPost = NULL;
		break;
	case SMX_COLLISION:
		assert(bSymmetric == 0);
		smx->fcnSmooth = CheckForCollision;
		initParticle = NULL;
		init = NULL;
		comb = NULL;
		smx->fcnPost = NULL;
		break;
#endif /* PLANETS */
	default:
		assert(0);
		}
	/*
	 ** Initialize the particles in the tree.
	 */
	nTree = pkd->kdNodes[pkd->iRoot].pUpper + 1;
	for (pi=0;pi<nTree;++pi) {
		if (pkd->pStore[pi].iActive) {
			if (bSmooth) pkd->pStore[pi].fBall2 = -1.0;
			if (initParticle != NULL) {
				initParticle(&pkd->pStore[pi]);
				}
			}
		else {
			assert(bSymmetric == 0);
			if (bSmooth) {
				if (pkd->pStore[pi].fBall2 < 0.0) 
					pkd->pStore[pi].fBall2 = 0.0;
				}
			}
		}
	/*
	 ** Start particle caching space (cell cache is already active).
	 */
	if (bSymmetric) {
		mdlCOcache(pkd->mdl,CID_PARTICLE,pkd->pStore,sizeof(PARTICLE),
				   nTree,init,comb);
		}
	else {
		mdlROcache(pkd->mdl,CID_PARTICLE,pkd->pStore,sizeof(PARTICLE),
				   nTree);
		}
	/*
	 ** Allocate mark array.
	 */
	smx->piMark = (int *)malloc(pkdLocal(pkd)*sizeof(int));
	assert(smx->piMark != NULL);
	/*
	 ** Allocate Nearest-Neighbor List.
	 */
	smx->nListSize = smx->nSmooth+RESMOOTH_SAFE;
	smx->nnList = (NN *)malloc(smx->nListSize*sizeof(NN));
	assert(smx->nnList != NULL);
	smx->pbRelease = (int *)malloc(smx->nListSize*sizeof(int));
	assert(smx->pbRelease != NULL);
	/*
	 ** Allocate priority queue.
	 */
	smx->pq = (PQ *)malloc(nSmooth*sizeof(PQ));
	assert(smx->pq != NULL);
	PQ_INIT(smx->pq,nSmooth);
	/*
	 ** Set up "in priority queue" hash table stuff.
	 */
	smx->nHash = (int)ceil(nSmooth/PQ_LOAD_FACTOR);
	smx->pqHash = (PQ **)malloc(smx->nHash*sizeof(PQ *));
	assert(smx->pqHash != NULL);
	*psmx = smx;	
	return(1);
	}


void smFinish(SMX smx,SMF *smf)
{
	PKD pkd = smx->pkd;
	int pi;
    char achOut[128];

	/*
	 * Output statistics.
	 */
	sprintf(achOut, "Cell Accesses: %g\n",
		mdlNumAccess(smx->pkd->mdl,CID_CELL));
	mdlDiag(smx->pkd->mdl, achOut);
	sprintf(achOut, "    Miss ratio: %g\n",
		mdlMissRatio(smx->pkd->mdl,CID_CELL));
	mdlDiag(smx->pkd->mdl, achOut);
	sprintf(achOut, "    Min ratio: %g\n",
		mdlMinRatio(smx->pkd->mdl,CID_CELL));
	mdlDiag(smx->pkd->mdl, achOut);
	sprintf(achOut, "    Coll ratio: %g\n",
		mdlCollRatio(smx->pkd->mdl,CID_CELL));
	mdlDiag(smx->pkd->mdl, achOut);

	sprintf(achOut, "Particle Accesses: %g\n",
		mdlNumAccess(smx->pkd->mdl,CID_PARTICLE));
	mdlDiag(smx->pkd->mdl, achOut);
	sprintf(achOut, "    Miss ratio: %g\n",
		mdlMissRatio(smx->pkd->mdl,CID_PARTICLE));
	mdlDiag(smx->pkd->mdl, achOut);
	sprintf(achOut, "    Min ratio: %g\n",
		mdlMinRatio(smx->pkd->mdl,CID_PARTICLE));
	mdlDiag(smx->pkd->mdl, achOut);
	sprintf(achOut, "    Coll ratio: %g\n",
		mdlCollRatio(smx->pkd->mdl,CID_PARTICLE));
	mdlDiag(smx->pkd->mdl, achOut);
	/*
	 ** Stop particle caching space.
	 */
	mdlFinishCache(smx->pkd->mdl,CID_PARTICLE);
	/*
	 ** Now do any post calculations, these ususlly involve some sort of
	 ** normalizations of the smoothed quantities, usually division by
	 ** the local density! Do NOT put kernel normalizations in here as
	 ** these do not depend purely on local properties in the case of
	 ** "Gather-Scatter" kernel.
	 */
	if (smx->fcnPost != NULL) {
		for (pi=0;pi<pkd->nActive;++pi) {
			if (pkd->pStore[pi].iActive) {
				smx->fcnPost(&pkd->pStore[pi],smf);
				}
			}
		}
	/*
	 ** Free up context storage.
	 */
	free(smx->pqHash);
	free(smx->pq);
	free(smx->piMark);
	free(smx->nnList);
	free(smx->pbRelease);
	free(smx);
	}


PQ *smBallSearch(SMX smx,PQ *pq,float *ri,int *cpStart)
{
	MDL mdl = smx->pkd->mdl;
	KDN *c = smx->pkd->kdNodes;
	PARTICLE *p = smx->pkd->pStore;
	int cell,cp,ct,pj,pUpper,idSelf,bPeriodic;
	float fBall2,fDist2,dx,dy,dz,lx,ly,lz,sx,sy,sz,x,y,z;

	x = ri[0];
	y = ri[1];
	z = ri[2];
	if(smx->bPeriodic) {
	    lx = smx->pkd->fPeriod[0];
	    ly = smx->pkd->fPeriod[1];
	    lz = smx->pkd->fPeriod[2];
	    }
	else {
	    lx = FLT_MAX;
	    ly = FLT_MAX;
	    lz = FLT_MAX;
	    }
	idSelf = smx->pkd->idSelf;
	fBall2 = pq->fKey;
	cell = ROOT;
	bPeriodic = 0;
	/*
	 ** First find the "local" Bucket.
	 ** This could be the closest to ri[3]! 
	 ** Warning: in parallel ri[3] SHOULD be contained in the LOCAL DOMAIN!
	 */
	while (c[cell].iDim >= 0) {
		if (ri[c[cell].iDim] < c[cell].fSplit) cell = LOWER(cell);
		else cell = UPPER(cell);
		}
	pUpper = c[cell].pUpper;
	for (pj=c[cell].pLower;pj<=pUpper;++pj) {
		dx = x - p[pj].r[0];
		dy = y - p[pj].r[1];
		dz = z - p[pj].r[2];
		fDist2 = dx*dx + dy*dy + dz*dz;
		if (fDist2 < fBall2) {
			*cpStart = cell;
			if (smx->piMark[pj]) continue;
			if (pq->id == idSelf) smx->piMark[pq->p] = 0;
			else {
				mdlRelease(mdl,CID_PARTICLE,pq->pPart);
				PQ_HASHDEL(smx->pqHash,smx->nHash,pq);
				pq->id = idSelf;
				}
			smx->piMark[pj] = 1;
			pq->fKey = fDist2;
			pq->dx = dx;
			pq->dy = dy;
			pq->dz = dz;
			pq->p = pj;
			pq->pPart = &p[pj];
			pq->ax = 0.0;
			pq->ay = 0.0;
			pq->az = 0.0;
			PQ_REPLACE(pq);
			fBall2 = pq->fKey;
			}
		}
	while (cell != ROOT) {
		cp = SIBLING(cell);
		ct = cp;
		SETNEXT(ct);
		while (1) {
			INTERSECT(&c[cp],fBall2,lx,ly,lz,x,y,z,sx,sy,sz,bPeriodic,
					  GetNext_1);
			if (c[cp].iDim >= 0) {
				cp = LOWER(cp);
				continue;
				}
			else {
				pUpper = c[cp].pUpper;
				for (pj=c[cp].pLower;pj<=pUpper;++pj) {
					dx = sx - p[pj].r[0];
					dy = sy - p[pj].r[1];
					dz = sz - p[pj].r[2];
					fDist2 = dx*dx + dy*dy + dz*dz;
					if (fDist2 < fBall2) {
						*cpStart = PARENT(cell);
						if (smx->piMark[pj]) continue;
						if (pq->id == idSelf) smx->piMark[pq->p] = 0;
						else {
							mdlRelease(mdl,CID_PARTICLE,pq->pPart);
							PQ_HASHDEL(smx->pqHash,smx->nHash,pq);
							pq->id = idSelf;
							}
						smx->piMark[pj] = 1;
						pq->fKey = fDist2;
						pq->dx = dx;
						pq->dy = dy;
						pq->dz = dz;
						pq->p = pj;
						pq->pPart = &p[pj];
						pq->ax = sx - x;
						pq->ay = sy - y;
						pq->az = sz - z;
						PQ_REPLACE(pq);
						fBall2 = pq->fKey;
						}
					}
				}
		GetNext_1:
			SETNEXT(cp);
			if (cp == ct) break;
			}
		cell = PARENT(cell);
		}
	if (bPeriodic) *cpStart = 0;
	return(pq);
	}


PQ *smBallSearchNP(SMX smx,PQ *pq,float *ri,int *cpStart)
{
	MDL mdl = smx->pkd->mdl;
	KDN *c = smx->pkd->kdNodes;
	PARTICLE *p = smx->pkd->pStore;
	int cell,cp,ct,pj,pUpper,idSelf;
	float fBall2,fDist2,dx,dy,dz,x,y,z;

	x = ri[0];
	y = ri[1];
	z = ri[2];
	idSelf = smx->pkd->idSelf;
	fBall2 = pq->fKey;
	cell = ROOT;
	/*
	 ** First find the "local" Bucket.
	 ** This could mearly be the closest to ri[3]! 
	 ** Warning: in parallel ri[3] SHOULD be contained in the LOCAL DOMAIN!
	 */
	while (c[cell].iDim >= 0) {
		if (ri[c[cell].iDim] < c[cell].fSplit) cell = LOWER(cell);
		else cell = UPPER(cell);
		}
	pUpper = c[cell].pUpper;
	for (pj=c[cell].pLower;pj<=pUpper;++pj) {
		dx = x - p[pj].r[0];
		dy = y - p[pj].r[1];
		dz = z - p[pj].r[2];
		fDist2 = dx*dx + dy*dy + dz*dz;
		if (fDist2 < fBall2) {
			*cpStart = cell;
			if (smx->piMark[pj]) continue;
			if (pq->id == idSelf) smx->piMark[pq->p] = 0;
			else {
				mdlRelease(mdl,CID_PARTICLE,pq->pPart);
				PQ_HASHDEL(smx->pqHash,smx->nHash,pq);
				pq->id = idSelf;
				}
			smx->piMark[pj] = 1;
			pq->fKey = fDist2;
			pq->dx = dx;
			pq->dy = dy;
			pq->dz = dz;
			pq->p = pj;
			pq->pPart = &p[pj];
			PQ_REPLACE(pq);
			fBall2 = pq->fKey;
			}
		}
	while (cell != ROOT) {
		cp = SIBLING(cell);
		ct = cp;
		SETNEXT(ct);
		while (1) {
			INTERSECTNP(&c[cp],fBall2,x,y,z,GetNext_1);
			if (c[cp].iDim >= 0) {
				cp = LOWER(cp);
				continue;
				}
			else {
				pUpper = c[cp].pUpper;
				for (pj=c[cp].pLower;pj<=pUpper;++pj) {
					dx = x - p[pj].r[0];
					dy = y - p[pj].r[1];
					dz = z - p[pj].r[2];
					fDist2 = dx*dx + dy*dy + dz*dz;
					if (fDist2 < fBall2) {
						*cpStart = PARENT(cell);
						if (smx->piMark[pj]) continue;
						if (pq->id == idSelf) smx->piMark[pq->p] = 0;
						else {
							mdlRelease(mdl,CID_PARTICLE,pq->pPart);
							PQ_HASHDEL(smx->pqHash,smx->nHash,pq);
							pq->id = idSelf;
							}
						smx->piMark[pj] = 1;
						pq->fKey = fDist2;
						pq->dx = dx;
						pq->dy = dy;
						pq->dz = dz;
						pq->p = pj;
						pq->pPart = &p[pj];
						PQ_REPLACE(pq);
						fBall2 = pq->fKey;
						}
					}
				}
		GetNext_1:
			SETNEXT(cp);
			if (cp == ct) break;
			}
		cell = PARENT(cell);
		}
	return(pq);
	}


int smBallGather(SMX smx,float fBall2,float *ri)
{
	KDN *c = smx->pkd->kdNodes;
	PARTICLE *p = smx->pkd->pStore;
	int pj,nCnt,cp,pUpper;
	float dx,dy,dz,x,y,z,lx,ly,lz,sx,sy,sz,fDist2;
	int iDum;

	x = ri[0];
	y = ri[1];
	z = ri[2];
	if (smx->bPeriodic) {
	    lx = smx->pkd->fPeriod[0];
	    ly = smx->pkd->fPeriod[1];
	    lz = smx->pkd->fPeriod[2];
	    }
	else {
	    lx = FLT_MAX;
	    ly = FLT_MAX;
	    lz = FLT_MAX;
	    }
	nCnt = 0;
	cp = ROOT;
	while (1) {
		INTERSECT(&c[cp],fBall2,lx,ly,lz,x,y,z,sx,sy,sz,iDum,GetNextCell);
		/*
		 ** We have an intersection to test.
		 */
		if (c[cp].iDim >= 0) {
			cp = LOWER(cp);
			continue;
			}
		else {
			pUpper = c[cp].pUpper;
			for (pj=c[cp].pLower;pj<=pUpper;++pj) {
				dx = sx - p[pj].r[0];
				dy = sy - p[pj].r[1];
				dz = sz - p[pj].r[2];
				fDist2 = dx*dx + dy*dy + dz*dz;
				if (fDist2 <= fBall2) {
					smx->nnList[nCnt].fDist2 = fDist2;
					smx->nnList[nCnt].dx = dx;
					smx->nnList[nCnt].dy = dy;
					smx->nnList[nCnt].dz = dz;
					smx->nnList[nCnt].pPart = &p[pj];
					smx->pbRelease[nCnt++] = 0;
					}
				}
			}
	GetNextCell:
		SETNEXT(cp);
		if (cp == ROOT) break;
		}
	assert(nCnt <= smx->nListSize);
	return(nCnt);
	}


int smBallGatherNP(SMX smx,float fBall2,float *ri,int cp)
{
	KDN *c = smx->pkd->kdNodes;
	PARTICLE *p = smx->pkd->pStore;
	int pj,nCnt,pUpper;
	float dx,dy,dz,x,y,z,fDist2;

	x = ri[0];
	y = ri[1];
	z = ri[2];
	nCnt = 0;
	while (1) {
		INTERSECTNP(&c[cp],fBall2,x,y,z,GetNextCell);
		/*
		 ** We have an intersection to test.
		 */
		if (c[cp].iDim >= 0) {
			cp = LOWER(cp);
			continue;
			}
		else {
			pUpper = c[cp].pUpper;
			for (pj=c[cp].pLower;pj<=pUpper;++pj) {
				dx = x - p[pj].r[0];
				dy = y - p[pj].r[1];
				dz = z - p[pj].r[2];
				fDist2 = dx*dx + dy*dy + dz*dz;
				if (fDist2 <= fBall2) {
					smx->nnList[nCnt].fDist2 = fDist2;
					smx->nnList[nCnt].dx = dx;
					smx->nnList[nCnt].dy = dy;
					smx->nnList[nCnt].dz = dz;
					smx->nnList[nCnt].pPart = &p[pj];
					smx->pbRelease[nCnt++] = 0;
					}
				}
			}
	GetNextCell:
		SETNEXT(cp);
		if (cp == ROOT) break;
		}
	assert(nCnt <= smx->nListSize);
	return(nCnt);
	}


void smSmooth(SMX smx,SMF *smf)
{
	PKD pkd = smx->pkd;
	MDL mdl = smx->pkd->mdl;
	KDN *c = pkd->kdNodes;
	PARTICLE *p = pkd->pStore;
	KDN *pkdn;
	PARTICLE *pPart;
	int nSmooth,i,j,pi,pj,pNext,nCnt;
	int cell,idcell,cp,id,ct,idct;
	float fBall2,fDist2,x,y,z,dx,dy,dz,lx,ly,lz,sx,sy,sz,h2;
	PQ *pq,*pqi,*pqn;
	int iDum;
	int nTree,nQueue,iLoad;

	nSmooth = smx->nSmooth;
	nTree = c[pkd->iRoot].pUpper + 1;
	if (smx->bPeriodic) {
	    lx = smx->pkd->fPeriod[0];
	    ly = smx->pkd->fPeriod[1];
	    lz = smx->pkd->fPeriod[2];
	    }
	else {
	    lx = FLT_MAX;
	    ly = FLT_MAX;
	    lz = FLT_MAX;
	    }
	/*
	 ** Clear Mark array and pqHash.
	 */
	for (pi=0;pi<pkd->nLocal;++pi) smx->piMark[pi] = 0;
	for (j=0;j<smx->nHash;++j) smx->pqHash[j] = NULL;
	pNext = 0;
 StartParticle:
	/*
	 ** Check if we are finished!
	 */
	if (pNext == nTree) goto DoneSmooth;
	if (pkd->pStore[pNext].fBall2 >= 0.0) {
		++pNext;
		goto StartParticle;
		}
	pi = pNext;
	x = p[pi].r[0];
	y = p[pi].r[1];
	z = p[pi].r[2];
	cell = ROOT;
	while (c[cell].iDim >= 0) {
		if (p[pi].r[c[cell].iDim] < c[cell].fSplit) cell = LOWER(cell);
		else cell = UPPER(cell);
		}
	/*
	 ** Add local stuff to the prioq.
	 */
	iLoad = c[cell].pLower;
	nQueue = 0;
	if (iLoad > nTree - nSmooth) iLoad = nTree - nSmooth;
	/*
	 ** Have to add non-local particles to the queue?
	 */
	for (id=0;id<mdlThreads(pkd->mdl) && (iLoad < 0);++id) {
		if (id == pkd->idSelf) continue;
		pkdn = mdlAquire(mdl,CID_CELL,ROOT,id);
		for (pj=pkdn->pLower;pj<=pkdn->pUpper && (iLoad < 0);++pj) {
			pqi = &smx->pq[nQueue];
			pPart = mdlAquire(mdl,CID_PARTICLE,pj,id);
			dx = x - pPart->r[0];
			dy = y - pPart->r[1];
			dz = z - pPart->r[2];
			fDist2 = dx*dx + dy*dy + dz*dz;
			pqi->fKey = fDist2;
			pqi->dx = dx;
			pqi->dy = dy;
			pqi->dz = dz;
			pqi->p = pj;
			pqi->id = id;
			pqi->pPart = pPart;
			pqi->ax = 0.0;
			pqi->ay = 0.0;
			pqi->az = 0.0;
			PQ_HASHADD(smx->pqHash,smx->nHash,pqi);
			/*
			 ** Note: in this case we DO NOT want to release pPart!
			 */
			++iLoad;
			++nQueue;
			}
		mdlRelease(mdl,CID_CELL,pkdn);
		}
	assert(iLoad >= 0);
	for (;nQueue<nSmooth;++nQueue,++iLoad) {
		pqi = &smx->pq[nQueue];
		smx->piMark[iLoad] = 1;
		dx = x - p[iLoad].r[0];
		dy = y - p[iLoad].r[1];
		dz = z - p[iLoad].r[2];
		pqi->fKey = dx*dx + dy*dy + dz*dz;
		pqi->dx = dx;
		pqi->dy = dy;
		pqi->dz = dz;
		pqi->p = iLoad;
		pqi->id = pkd->idSelf;
		pqi->pPart = &p[iLoad];
		pqi->ax = 0.0;
		pqi->ay = 0.0;
		pqi->az = 0.0;
		}
	PQ_BUILD(smx->pq,nSmooth,pq);
/*
	sprintf(ach,"Start:%d\n",p[pi].iOrder);
	mdlDiag(pkd->mdl,ach);
*/
 PrioqBuilt:
	/*
	 ** Priority Queue must be built. 'pi' must be defined.
	 */
	if (smx->bPeriodic) {
		pq = smBallSearch(smx,pq,p[pi].r,&p[pi].cpStart);
		}
	else {
		pq = smBallSearchNP(smx,pq,p[pi].r,&p[pi].cpStart);
	    }
	/*
	 ** Start non-local search.
	 */
    fBall2 = pq->fKey;
	idcell = -1;	/* We are in the LTT now ! */
	cell = pkd->piLeaf[pkd->idSelf];
	while (!pkdIsRoot(cell,idcell)) {
		cp = SIBLING(cell);
		id = idcell;
		ct = cp;
		idct = id;
		pkdNext(pkd,ct,idct);
		/*
		 ** Check for any cache work.
		 */
		while (1) {
		    if (id >= 0) pkdn = mdlAquire(mdl,CID_CELL,cp,id);
			else pkdn = &pkd->kdTop[cp];
   			INTERSECT(pkdn,fBall2,lx,ly,lz,x,y,z,sx,sy,sz,iDum,GetNext_2);
			if (pkdn->iDim >= 0) {
			    if (id >= 0) mdlRelease(mdl,CID_CELL,pkdn);
				pkdLower(pkd,cp,id);
			    continue;
				}
			for (pj=pkdn->pLower;pj<=pkdn->pUpper;++pj) {
				pPart = mdlAquire(mdl,CID_PARTICLE,pj,id);
				dx = sx - pPart->r[0];
				dy = sy - pPart->r[1];
				dz = sz - pPart->r[2];
				fDist2 = dx*dx + dy*dy + dz*dz;
				if (fDist2 < fBall2) {
					p[pi].cpStart = 0;
					PQ_INQUEUE(smx->pqHash,smx->nHash,pj,id,NextParticle);
					if (pq->id == pkd->idSelf) smx->piMark[pq->p] = 0;
					else {
						mdlRelease(mdl,CID_PARTICLE,pq->pPart);
						PQ_HASHDEL(smx->pqHash,smx->nHash,pq);
						}
					pq->fKey = fDist2;
					pq->dx = dx;
					pq->dy = dy;
					pq->dz = dz;
					pq->p = pj;
					pq->id = id;
					pq->pPart = pPart;
					pq->ax = sx - x;
					pq->ay = sy - y;
					pq->az = sz - z;
					PQ_HASHADD(smx->pqHash,smx->nHash,pq);
					PQ_REPLACE(pq);
					fBall2 = pq->fKey;
					/*
					 ** Note: in this case we DO NOT want to release pPart!
					 */
					/*
					 ** Setting p[pi].cpStart to zero means we have to search 
					 ** non-local tree when doing a BallGather (smReSmooth).
					 */
					continue;
					}
			NextParticle:
				mdlRelease(mdl,CID_PARTICLE,pPart);
				}
		GetNext_2:
			if (id >= 0) mdlRelease(mdl,CID_CELL,pkdn);
			pkdNext(pkd,cp,id);
			if (cp == ct && id == idct) break;
			}
		pkdParent(pkd,cell,idcell);
		}
	/*
	 ** Create Nearest-Neighbor List and try to pick next particle.
	 */
	pqn = NULL;
	nCnt = 0;
	h2 = 2.0*pq->fKey;   /* arbitrarily bigger than pq->fKey! */
	p[pi].fBall2 = h2;
	for (i=0,pqi=smx->pq;i<nSmooth;++i,++pqi) {
/*
 ** There are cases like in collisions where we do want to include the
 ** single nearest neighbor. So include the most distant particle.
 **
		if (pqi == pq) continue;
*/
		/*
		 ** Move relevant data into Nearest Neighbor array.
		 */
		smx->nnList[nCnt].pPart = pqi->pPart;
		smx->nnList[nCnt].fDist2 = pqi->fKey;
		smx->nnList[nCnt].dx = pqi->dx;
		smx->nnList[nCnt].dy = pqi->dy;
		smx->nnList[nCnt].dz = pqi->dz;
		++nCnt;
		if (pqi->id != pkd->idSelf) continue;
		if (pkd->pStore[pqi->p].fBall2 >= 0) continue;
		if (pqi->fKey < h2) {
			pqn = pqi;
			h2 = pqn->fKey;
			}
		}
	smx->fcnSmooth(&p[pi],nCnt,smx->nnList,smf);
	/*
	 ** Try a cache check to improve responsiveness.
	 */
	mdlCacheCheck(mdl);
	if (!pqn) {
		/*
		 ** Clean up piMark array, pqHash, and release all aquired
		 ** pointers!
		 */
		for (i=0,pqi=smx->pq;i<nSmooth;++i,++pqi) {
			if (pqi->id == pkd->idSelf) smx->piMark[pqi->p] = 0;
			else {
				mdlRelease(mdl,CID_PARTICLE,pqi->pPart);
				PQ_HASHDEL(smx->pqHash,smx->nHash,pqi);
				}
			}
		goto StartParticle;
		}
	/*
	 ** Calculate the priority queue using the previous particles.
	 ** "THE SNAKE"
	 */
	pi = pqn->p;
	x = p[pi].r[0];
	y = p[pi].r[1];
	z = p[pi].r[2];
	for (i=0,pqi=smx->pq;i<nSmooth;++i,++pqi) {
		if (pqi == pqn) continue;
		pqi->ax -= pqn->ax;
		pqi->ay -= pqn->ay;
		pqi->az -= pqn->az;
		dx = x - pqi->pPart->r[0] + pqi->ax;
		dy = y - pqi->pPart->r[1] + pqi->ay;
		dz = z - pqi->pPart->r[2] + pqi->az;
		pqi->fKey = dx*dx + dy*dy + dz*dz;
		pqi->dx = dx;
		pqi->dy = dy;
		pqi->dz = dz;
		}
	pqn->fKey = 0.0;
	pqn->ax = 0.0;
	pqn->ay = 0.0;
	pqn->az = 0.0;
	PQ_BUILD(smx->pq,nSmooth,pq);
/*
	sprintf(ach,"Snake:%d\n",p[pi].iOrder);
	mdlDiag(pkd->mdl,ach);
*/
	goto PrioqBuilt;

 DoneSmooth:
	;
	}


void smReSmooth(SMX smx,SMF *smf)
{
	PKD pkd = smx->pkd;
	MDL mdl = smx->pkd->mdl;
	PARTICLE *p = smx->pkd->pStore;
	PARTICLE *pPart;
	KDN *pkdn;
	int pi,pj,nCnt,cp,id,i;
	float x,y,z,lx,ly,lz,sx,sy,sz,dx,dy,dz,fDist2,fBall2;
	int iDum;
	int nTree;
	
	if (smx->bPeriodic) {
	    lx = pkd->fPeriod[0];
	    ly = pkd->fPeriod[1];
	    lz = pkd->fPeriod[2];
	    }
	else {
	    lx = FLT_MAX;
	    ly = FLT_MAX;
	    lz = FLT_MAX;
	    }
	nTree = pkd->kdNodes[pkd->iRoot].pUpper + 1;
	for (pi=0;pi<nTree;++pi) {
		if (p[pi].iActive == 0) continue;
		/*
		 ** Do a Ball Gather at the radius of the most distant particle
		 ** which is smSmooth sets in p[pi].fBall2.
		 */
		fBall2 = p[pi].fBall2;
		cp = p[pi].cpStart;
		if (cp) {
			nCnt = smBallGatherNP(smx,fBall2,p[pi].r,cp);
			smx->fcnSmooth(&p[pi],nCnt,smx->nnList,smf);
			}
		else {
			if (smx->bPeriodic) {
				nCnt = smBallGather(smx,fBall2,p[pi].r);
				}
			else {
				nCnt = smBallGatherNP(smx,fBall2,p[pi].r,ROOT);
				}
			/*
			 ** Start non-local search.
			 */
			x = p[pi].r[0];
			y = p[pi].r[1];
			z = p[pi].r[2];
			cp = ROOT;
			id = -1;	/* We are in the LTT now ! */
			while (1) {
				if (id == pkd->idSelf) goto SkipLocal;
				if (id >= 0) pkdn = mdlAquire(mdl,CID_CELL,cp,id);
				else pkdn = &pkd->kdTop[cp];
				INTERSECT(pkdn,fBall2,lx,ly,lz,x,y,z,sx,sy,sz,iDum,
						  GetNextCell);
				if (pkdn->iDim >= 0) {
					if (id >= 0) mdlRelease(mdl,CID_CELL,pkdn);
					pkdLower(pkd,cp,id);
					continue;
					}
				for (pj=pkdn->pLower;pj<=pkdn->pUpper;++pj) {
					pPart = mdlAquire(mdl,CID_PARTICLE,pj,id);
					dx = sx - pPart->r[0];
					dy = sy - pPart->r[1];
					dz = sz - pPart->r[2];
					fDist2 = dx*dx + dy*dy + dz*dz;
					if (fDist2 <= fBall2) {
						smx->nnList[nCnt].fDist2 = fDist2;
						smx->nnList[nCnt].dx = dx;
						smx->nnList[nCnt].dy = dy;
						smx->nnList[nCnt].dz = dz;
						smx->nnList[nCnt].pPart = pPart;
						smx->pbRelease[nCnt++] = 1;
						continue;
						}
					mdlRelease(mdl,CID_PARTICLE,pPart);
					}
			GetNextCell:
				if (id >= 0) mdlRelease(mdl,CID_CELL,pkdn);
			SkipLocal:
				pkdNext(pkd,cp,id);
				if (pkdIsRoot(cp,id)) break;
				}
			smx->fcnSmooth(&p[pi],nCnt,smx->nnList,smf);
			/*
			 ** Release non-local particle pointers.
			 */
			for (i=0;i<nCnt;++i) {
				if (smx->pbRelease[i]) {
					mdlRelease(mdl,CID_PARTICLE,smx->nnList[i].pPart);
					}
				}
			}
		}
	}




