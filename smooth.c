#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <malloc.h>
#include <math.h>
#include <assert.h>
#include "smooth.h"
#include "pkd.h"



#define INTERSECT(pkdn,fBall2,lx,ly,lz,x,y,z,sx,sy,sz,label)\
{\
	float dx,dy,dz,dx1,dy1,dz1,fDist2;\
	dx = (pkdn)->bnd.fMin[0]-x;\
	dx1 = x-(pkdn)->bnd.fMax[0];\
	dy = (pkdn)->bnd.fMin[1]-y;\
	dy1 = y-(pkdn)->bnd.fMax[1];\
	dz = (pkdn)->bnd.fMin[2]-z;\
	dz1 = z-(pkdn)->bnd.fMax[2];\
	if (dx > 0.0) {\
		dx1 += lx;\
		if (dx1 < dx) {\
			fDist2 = dx1*dx1;\
			sx = x+lx;\
			}\
		else {\
			fDist2 = dx*dx;\
			sx = x;\
			}\
		if (fDist2 > fBall2) goto label;\
		}\
	else if (dx1 > 0.0) {\
		dx += lx;\
		if (dx < dx1) {\
			fDist2 = dx*dx;\
			sx = x-lx;\
			}\
		else {\
			fDist2 = dx1*dx1;\
			sx = x;\
			}\
		if (fDist2 > fBall2) goto label;\
		}\
	else {\
		fDist2 = 0.0;\
		sx = x;\
		}\
	if (dy > 0.0) {\
		dy1 += ly;\
		if (dy1 < dy) {\
			fDist2 += dy1*dy1;\
			sy = y+ly;\
			}\
		else {\
			fDist2 += dy*dy;\
			sy = y;\
			}\
		if (fDist2 > fBall2) goto label;\
		}\
	else if (dy1 > 0.0) {\
		dy += ly;\
		if (dy < dy1) {\
			fDist2 += dy*dy;\
			sy = y-ly;\
			}\
		else {\
			fDist2 += dy1*dy1;\
			sy = y;\
			}\
		if (fDist2 > fBall2) goto label;\
		}\
	else {\
		sy = y;\
		}\
	if (dz > 0.0) {\
		dz1 += lz;\
		if (dz1 < dz) {\
			fDist2 += dz1*dz1;\
			sz = z+lz;\
			}\
		else {\
			fDist2 += dz*dz;\
			sz = z;\
			}\
		if (fDist2 > fBall2) goto label;\
		}\
	else if (dz1 > 0.0) {\
		dz += lz;\
		if (dz < dz1) {\
			fDist2 += dz*dz;\
			sz = z-lz;\
			}\
		else {\
			fDist2 += dz1*dz1;\
			sz = z;\
			}\
		if (fDist2 > fBall2) goto label;\
		}\
	else {\
		sz = z;\
		}\
	}


void initDen(void *p)
{
	((PARTICLE *)p)->fDensity = 0.0;
	}

void combDen(void *p1,void *p2)
{
	((PARTICLE *)p1)->fDensity += ((PARTICLE *)p2)->fDensity;
	}

int smInitialize(SMX *psmx,PKD pkd,int nSmooth,int bCombiner)
{
	SMX smx;
	int pi;

	smx = (SMX)malloc(sizeof(struct smContext));
	assert(smx != NULL);
	smx->pkd = pkd;
	smx->nSmooth = nSmooth;
	/*
	 ** Initialize Densities of the local particles.
	 */
	for (pi=0;pi<pkd->nLocal;++pi) {
		pkd->pStore[pi].fDensity = 0.0;
		pkd->pStore[pi].fBall2 = -1.0;
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
	/*
	 ** Start caching spaces.
	 */
	if (bCombiner) {
		mdlCOcache(pkd->mdl,CID_PARTICLE,pkd->pStore,sizeof(PARTICLE),
				   pkdLocal(pkd),initDen,combDen);
		}
	else {
		mdlROcache(pkd->mdl,CID_PARTICLE,pkd->pStore,sizeof(PARTICLE),
				   pkdLocal(pkd));
		}
	mdlROcache(pkd->mdl,CID_CELL,pkd->kdNodes,sizeof(KDN),pkdNodes(pkd));
	*psmx = smx;	
	return(1);
	}


void smFinish(SMX smx)
{
	/*
	 ** Stop caching spaces.
	 */
	mdlFinishCache(smx->pkd->mdl,CID_CELL);
	mdlFinishCache(smx->pkd->mdl,CID_PARTICLE);
	/*
	 ** Free up context storage.
	 */
	free(smx->pqHash);
	free(smx->pq);
	free(smx->piMark);
	free(smx->nnList);
	free(smx);
	}


PQ *smBallSearch(SMX smx,PQ *pq,float *ri)
{
	MDL mdl = smx->pkd->mdl;
	KDN *c = smx->pkd->kdNodes;
	PARTICLE *p = smx->pkd->pStore;
	int cell,cp,ct,pj,pUpper,idSelf;
	float fBall2,fDist2,dx,dy,dz,lx,ly,lz,sx,sy,sz,x,y,z;

	x = ri[0];
	y = ri[1];
	z = ri[2];
	lx = smx->pkd->fPeriod[0];
	ly = smx->pkd->fPeriod[1];
	lz = smx->pkd->fPeriod[2];
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
			if (smx->piMark[pj]) continue;
			mdlRelease(mdl,CID_PARTICLE,pq->pPart);
			if (pq->id == idSelf) smx->piMark[pq->p] = 0;
			else {
				PQ_HASHDEL(smx->pqHash,smx->nHash,pq);
				pq->id = idSelf;
				}
			smx->piMark[pj] = 1;
			pq->fKey = fDist2;
			pq->p = pj;
			pq->pPart = mdlAquire(mdl,CID_PARTICLE,pj,idSelf);
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
			INTERSECT(&c[cp],fBall2,lx,ly,lz,x,y,z,sx,sy,sz,GetNext_1);
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
						if (smx->piMark[pj]) continue;
						mdlRelease(mdl,CID_PARTICLE,pq->pPart);
						if (pq->id == idSelf) smx->piMark[pq->p] = 0;
						else {
							PQ_HASHDEL(smx->pqHash,smx->nHash,pq);
							pq->id = idSelf;
							}
						smx->piMark[pj] = 1;
						pq->fKey = fDist2;
						pq->p = pj;
						pq->pPart = mdlAquire(mdl,CID_PARTICLE,pj,idSelf);
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
	return(pq);
	}


void smSmooth(SMX smx,void (*fncSmooth)(SMX,int,int,NN *))
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

	nSmooth = smx->nSmooth;
	lx = smx->pkd->fPeriod[0];
	ly = smx->pkd->fPeriod[1];
	lz = smx->pkd->fPeriod[2];
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
	if (pNext == pkd->nLocal) goto DoneSmooth;
	if (pkd->pStore[pNext].fBall2 >= 0) {
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
	pj = c[cell].pLower;
	if (pj > pkd->nLocal - nSmooth) pj = pkd->nLocal - nSmooth;
	for (i=0,pqi=smx->pq;i<nSmooth;++i,++pj,++pqi) {
		smx->piMark[pj] = 1;
		dx = x - p[pj].r[0];
		dy = y - p[pj].r[1];
		dz = z - p[pj].r[2];
		pqi->fKey = dx*dx + dy*dy + dz*dz;
		pqi->p = pj;
		pqi->id = pkd->idSelf;
		pqi->pPart = mdlAquire(mdl,CID_PARTICLE,pj,pkd->idSelf);
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
	pq = smBallSearch(smx,pq,p[pi].r);
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
   			INTERSECT(pkdn,fBall2,lx,ly,lz,x,y,z,sx,sy,sz,GetNext_2);
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
					PQ_INQUEUE(smx->pqHash,smx->nHash,pj,id,NextParticle);
					mdlRelease(mdl,CID_PARTICLE,pq->pPart);
					if (pq->id == pkd->idSelf) smx->piMark[pq->p] = 0;
					else {
						PQ_HASHDEL(smx->pqHash,smx->nHash,pq);
						}
					pq->fKey = fDist2;
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
	h2 = pq->fKey;
	p[pi].fBall2 = h2;
	for (i=0,pqi=smx->pq;i<nSmooth;++i,++pqi) {
		if (pqi == pq) continue;
		/*
		 ** Move relevant data into Nearest Neighbor array.
		 */
		smx->nnList[nCnt].pPart = pqi->pPart;
		smx->nnList[nCnt].fDist2 = pqi->fKey;
		++nCnt;
		if (pqi->id != pkd->idSelf) continue;
		if (pkd->pStore[pqi->p].fBall2 >= 0) continue;
		if (pqi->fKey < h2) {
			pqn = pqi;
			h2 = pqn->fKey;
			}
		}
	(*fncSmooth)(smx,pi,nCnt,smx->nnList);
	if (!pqn) {
		/*
		 ** Clean up piMark array, pqHash, and release all aquired
		 ** pointers!
		 */
		for (i=0,pqi=smx->pq;i<nSmooth;++i,++pqi) {
			mdlRelease(mdl,CID_PARTICLE,pqi->pPart);
			if (pqi->id == pkd->idSelf) smx->piMark[pqi->p] = 0;
			else {
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
		dx = x + pqi->ax - pqi->pPart->r[0];
		dy = y + pqi->ay - pqi->pPart->r[1];
		dz = z + pqi->az - pqi->pPart->r[2];
		pqi->fKey = dx*dx + dy*dy + dz*dz;
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


void smDensity(SMX smx,int pi,int nSmooth,NN *nnList)
{
	PARTICLE *p = &smx->pkd->pStore[pi];
	float ih2,r2,rs,fDensity;
	int i;

	ih2 = 4.0/p->fBall2;
	fDensity = 0.0;
	for (i=0;i<nSmooth;++i) {
		r2 = nnList[i].fDist2*ih2;
		rs = 2.0 - sqrt(r2);
		if (r2 < 1.0) rs = (1.0 - 0.75*rs*r2);
		else rs = 0.25*rs*rs*rs;
		fDensity += rs*nnList[i].pPart->fMass;
		}
	p->fDensity = M_1_PI*sqrt(ih2)*ih2*fDensity; 
	}


void smDensitySym(SMX smx,int pi,int nSmooth,NN *nnList)
{
	PARTICLE *p = &smx->pkd->pStore[pi];
	PARTICLE *q;
	float fNorm,ih2,r2,rs;
	int i;

	ih2 = 4.0/p->fBall2;
	fNorm = 0.5*M_1_PI*sqrt(ih2)*ih2;
	for (i=0;i<nSmooth;++i) {
		r2 = nnList[i].fDist2*ih2;
		rs = 2.0 - sqrt(r2);
		if (r2 < 1.0) rs = (1.0 - 0.75*rs*r2);
		else rs = 0.25*rs*rs*rs;
		rs *= fNorm;
		q = nnList[i].pPart;
		p->fDensity += rs*q->fMass;
		q->fDensity += rs*p->fMass;
		}
	}







