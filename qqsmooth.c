#ifdef OLD_KEPLER /*DEBUG these routines are now broken -- DCR 00-04-11*/
#ifdef COLLISIONS

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <malloc.h>
#include <math.h>
#include <assert.h>
#include "smooth.h"
#include "pkd.h"
#include "smoothfcn.h"

#define CENTRAL_MASS 1.0 /*DEBUG needs to be reconciled with dCentMass*/

FLOAT PARTQQ(PARTICLE *p, int d)
{
    double r2, v2;
    double hx, hy, hz;		/* Angular momentum vector */
    double ex, ey, ez;		/* Runge-Lenz vector */
    double e, a;
    double ir;
    double imu;
    
    FLOAT x = p->r[0];
    FLOAT y = p->r[1];
    FLOAT z = p->r[2];
    FLOAT vx = p->v[0];
    FLOAT vy = p->v[1];
    FLOAT vz = p->v[2];
    
    r2 = x*x + y*y + z*z;
    v2 = vx*vx + vy*vy + vz*vz;
    
    hx = y*vz - vy*z;
    hy = z*vx - vz*x;
    hz = x*vy - vx*y;
    
    ir = 1.0/sqrt(r2);
    imu = 1.0/(CENTRAL_MASS + p->fMass);
    
    ex = (vy*hz - hy*vz)*imu - ir*x;
    ey = (vz*hx - hz*vx)*imu - ir*y;
    ez = (vx*hy - hx*vy)*imu - ir*z;
    e = sqrt(ex*ex + ey*ey + ez*ez);
    
    a = 1.0/(2*ir - v2*imu);
    if(d == 0)
	return a*(1-e);
    if(d == 1)
	return a*(1+e);
    assert(0);
	return 0;
}

void pkdQQCalcBound(PKD pkd,BND *pbnd,BND *pbndActive)
{
	double fq,fQ,fh;
	int i,j;

	/*
	 ** Initialize the bounds to 0 at the beginning
	 */
	for (j=0;j<2;++j) {
		pbnd->fMin[j] = FLOAT_MAXVAL;
		pbnd->fMax[j] = -FLOAT_MAXVAL;
		pbndActive->fMin[j] = FLOAT_MAXVAL;
		pbndActive->fMax[j] = -FLOAT_MAXVAL;
		}
	pbnd->fMin[2] = 0.0;
	pbnd->fMax[2] = 0.0;
	pbndActive->fMin[2] = 0.0;
	pbndActive->fMax[2] = 0.0;
	/*
	 ** Calculate Local Bounds, adjusted for Hill radius.
	 */
	for (i=0;i<pkd->nLocal;++i) {
		fq = PARTQQ(&pkd->pStore[i], 0);
		fQ = PARTQQ(&pkd->pStore[i], 1);
		fh = 0.5*(fq + fQ)*pkd->pStore[i].fHill;
		if (fq - fh < pbnd->fMin[0]) 
			pbnd->fMin[0] = fq - fh;
		if (fq > pbnd->fMax[0])
			pbnd->fMax[0] = fq;
		if (fQ < pbnd->fMin[1]) 
			pbnd->fMin[1] = fQ;
		if (fQ + fh > pbnd->fMax[1])
			pbnd->fMax[1] = fQ + fh;
		}
	/*
	 ** Calculate Active Bounds.
	 */
	for (i=0;i<pkd->nLocal;++i) {
		if (TYPEQueryACTIVE(&pkd->pStore[i])) {
			fq = PARTQQ(&pkd->pStore[i], 0);
			fQ = PARTQQ(&pkd->pStore[i], 1);
			fh = 0.5*(fq + fQ)*pkd->pStore[i].fHill;
			if (fq - fh < pbndActive->fMin[0]) 
				pbndActive->fMin[0] = fq - fh;
			if (fq > pbndActive->fMax[0])
				pbndActive->fMax[0] = fq;
			if (fQ < pbndActive->fMin[1]) 
				pbndActive->fMin[1] = fQ;
			if (fQ + fh > pbndActive->fMax[1])
				pbndActive->fMax[1] = fQ + fh;
			}
	    }
	}

int
pkdLowerQQPart(PKD pkd,int d,FLOAT fSplit,int i,int j)
{
	PARTICLE pTemp;

	if (i > j) goto done;
    while (1) {
        while (PARTQQ(&pkd->pStore[i], d-3) >= fSplit)
            if (++i > j) goto done;
        while (PARTQQ(&pkd->pStore[j], d-3) < fSplit)
            if (i > --j) goto done;
		pTemp = pkd->pStore[i];
		pkd->pStore[i] = pkd->pStore[j];
		pkd->pStore[j] = pTemp;
        }
 done:
    return(i);
	}


int
pkdUpperQQPart(PKD pkd,int d,FLOAT fSplit,int i,int j)
{
	PARTICLE pTemp;

	if (i > j) goto done;
    while (1) {
        while (PARTQQ(&pkd->pStore[i], d-3) < fSplit)
            if (++i > j) goto done;
        while (PARTQQ(&pkd->pStore[j], d-3) >= fSplit)
            if (i > --j) goto done;
		pTemp = pkd->pStore[i];
		pkd->pStore[i] = pkd->pStore[j];
		pkd->pStore[j] = pTemp;
        }
 done:
    return(i);
	}

void
pkdQQUpPass(PKD pkd,int iCell)
{
	KDN *c;
	PARTICLE *p;
	double fq,fQ,fh;
	int l,u,pj,j;

	c = pkd->kdNodes;
	p = pkd->pStore;
	l = c[iCell].pLower;
	u = c[iCell].pUpper;
	if (c[iCell].iDim >= 0) {
		pkdQQUpPass(pkd,LOWER(iCell));
		pkdQQUpPass(pkd,UPPER(iCell));
		pkdCombine(&c[LOWER(iCell)],&c[UPPER(iCell)],&c[iCell]);
		}
	else {
		c[iCell].fMass = 0.0;
		c[iCell].fSoft = 0.0;
		for (j=0;j<2;++j) {
			c[iCell].bnd.fMin[j] = FLOAT_MAXVAL;
			c[iCell].bnd.fMax[j] = -FLOAT_MAXVAL;
			}
		c[iCell].bnd.fMin[2] = 0.0;
		c[iCell].bnd.fMax[2] = 0.0;
		for (pj=l;pj<=u;++pj) {
			fq = PARTQQ(&p[pj], 0);
			fQ = PARTQQ(&p[pj], 1);
			fh = 0.5*(fq + fQ)*p[pj].fHill;
			if (fq - fh < c[iCell].bnd.fMin[0])
				c[iCell].bnd.fMin[0] = fq - fh;
			if (fq > c[iCell].bnd.fMax[0])
				c[iCell].bnd.fMax[0] = fq;
			if (fQ < c[iCell].bnd.fMin[1])
				c[iCell].bnd.fMin[1] = fQ;
			if (fQ + fh > c[iCell].bnd.fMax[1])
				c[iCell].bnd.fMax[1] = fQ + fh;
			/*
			 ** Find center of mass and total mass and mass weighted softening.
			 */
			c[iCell].fMass += p[pj].fMass;
			}
		}
	}

/*
 ** JST's Select
 */
void pkdQQSelect(PKD pkd,int d,int k,int l,int r)
{
	PARTICLE *p,t;
	FLOAT v;
	int i,j;

	p = pkd->pStore;
	while (r > l) {
		v = PARTQQ(&p[k],d);
		t = p[r];
		p[r] = p[k];
		p[k] = t;
		i = l - 1;
		j = r;
		while (1) {
			while (i < j) if (PARTQQ(&p[++i],d) >= v) break;
			while (i < j) if (PARTQQ(&p[--j],d) <= v) break;
			t = p[i];
			p[i] = p[j];
			p[j] = t;
			if (j <= i) break;
			}
		p[j] = p[i];
		p[i] = p[r];
		p[r] = t;
		if (i >= k) r = i - 1;
		if (i <= k) l = i + 1;
		}
	}

void
pkdQQBuild(PKD pkd,int nBucket, int bActiveOnly,KDN *pRoot)
{
	int l,n,i,d,m,j,diff;
	KDN *c;
	char ach[256];
	BND bndDum;

	/*
	 ** Make sure the particles are in Active/Inactive order.
	 */
	pkdActiveOrder(pkd);
	pkd->nBucket = nBucket;
	if (bActiveOnly) n = pkd->nActive;
	else n = pkd->nLocal;
	pkd->nLevels = 1;
	l = 1;
	while (n > nBucket) {
		n = n>>1;
		l = l<<1;
		++pkd->nLevels;
		}
	pkd->nSplit = l;
	pkd->nNodes = l<<1;
	if (pkd->kdNodes) {
		/*
		 ** Close caching space and free up nodes.
		 */
		mdlFinishCache(pkd->mdl,CID_CELL);
		mdlFree(pkd->mdl,pkd->kdNodes);
		}
	if(n == 0) {
	    pkd->kdNodes = NULL;
	    return;
	    }
	/*
	 ** Need to allocate a special extra cell that we will use to calculate
	 ** the acceleration on an arbitrary point in space.
	 */
	pkd->kdNodes = mdlMalloc(pkd->mdl,(pkd->nNodes + 1)*sizeof(KDN));
	assert(pkd->kdNodes != NULL);
	pkd->iFreeCell = pkd->nNodes;
	sprintf(ach,"nNodes:%d nSplit:%d nLevels:%d nBucket:%d\n",
			pkd->nNodes,pkd->nSplit,pkd->nLevels,nBucket);
	mdlDiag(pkd->mdl,ach);
	/*
	 ** Set up ROOT node
	 */
	c = pkd->kdNodes;
	pkd->iRoot = ROOT;
	c[pkd->iRoot].pLower = 0;
	if (bActiveOnly) {
		c[pkd->iRoot].pUpper = pkd->nActive-1;
		}
	else {
		c[pkd->iRoot].pUpper = pkd->nLocal-1;
		}
	/*
	 ** determine the local bound of the particles.
	 */
	if (bActiveOnly) {
		pkdQQCalcBound(pkd,&bndDum,&c[pkd->iRoot].bnd);
		}
	else {
		pkdQQCalcBound(pkd,&c[pkd->iRoot].bnd,&bndDum);
		}
	i = pkd->iRoot;
	while (1) {
		assert(c[i].pUpper - c[i].pLower + 1 > 0);
		if (i < pkd->nSplit && (c[i].pUpper - c[i].pLower) > 0) {
			d = 0;
			for (j=1;j<2;++j) { /* 2D tree */
				if (c[i].bnd.fMax[j]-c[i].bnd.fMin[j] > 
					c[i].bnd.fMax[d]-c[i].bnd.fMin[d]) d = j;
				}
			c[i].iDim = d;
			m = (c[i].pLower + c[i].pUpper)/2;
			pkdQQSelect(pkd,d,m,c[i].pLower,c[i].pUpper);
			c[i].fSplit = PARTQQ(&pkd->pStore[m],d);
			c[i].iLower = LOWER(i);
			c[i].iUpper = UPPER(i);
			c[LOWER(i)].bnd = c[i].bnd;
			c[LOWER(i)].bnd.fMax[d] = c[i].fSplit;
			c[LOWER(i)].pLower = c[i].pLower;
			c[LOWER(i)].pUpper = m;
			c[UPPER(i)].bnd = c[i].bnd;
			c[UPPER(i)].bnd.fMin[d] = c[i].fSplit;
			c[UPPER(i)].pLower = m+1;
			c[UPPER(i)].pUpper = c[i].pUpper;
			diff = (m-c[i].pLower+1)-(c[i].pUpper-m);
			assert(diff == 0 || diff == 1);
			i = LOWER(i);
			}
		else {
			c[i].iDim = -1;		/* to indicate a bucket! */
			c[i].iLower = -1;
			c[i].iUpper = -1;
			SETNEXT(i);
			if (i == pkd->iRoot) break;
			}
		}
	pkdQQUpPass(pkd,pkd->iRoot);
	/*
	 ** Thread the tree.
	 */
	pkdThreadTree(pkd,pkd->iRoot,-1);
	*pRoot = c[pkd->iRoot];
	/*
	 ** Finally activate a read only cache for remote access.
	 */
	mdlROcache(pkd->mdl,CID_CELL,pkd->kdNodes,sizeof(KDN),pkdNodes(pkd));
	}

#define INTERSECTQQ(pkdn,fq,fQ,label)\
{\
	if((pkdn)->bnd.fMin[0] > fQ || (pkdn)->bnd.fMax[1] < fq) goto label;\
	}

int smGatherQQ(SMX smx,FLOAT fq,FLOAT fQ,int cp)
{
	KDN *c = smx->pkd->kdNodes;
	PARTICLE *p = smx->pkd->pStore;
	double fqp,fQp,fhp;
	int pj,nCnt,pUpper;

	nCnt = 0;
	while (1) {
		INTERSECTQQ(&c[cp],fq,fQ,GetNextCell);
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
				fqp = PARTQQ(&p[pj], 0);
				fQp = PARTQQ(&p[pj], 1);
				fhp = 0.5*(fqp + fQp)*p[pj].fHill;
				if (fqp - fhp < fQ && fQp + fhp > fq) {
					if(nCnt >= smx->nListSize)
					    smGrowList(smx);
					smx->nnList[nCnt].pPart = &p[pj];
					smx->nnList[nCnt].iIndex = pj;
					smx->nnList[nCnt].iPid = smx->pkd->idSelf;
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

void smQQSmooth(SMX smx,SMF *smf)
{
	PKD pkd = smx->pkd;
	MDL mdl = smx->pkd->mdl;
	PARTICLE *p = smx->pkd->pStore;
	PARTICLE *pPart;
	KDN *pkdn;
	int pi,pj,nCnt,cp,id,i;
	FLOAT fq,fQ,fh;		/* peri, apo, Hill radius */
	double fqp,fQp,fhp;
	int nTree;
	
	nTree = pkd->kdNodes[pkd->iRoot].pUpper + 1;
	for (pi=0;pi<nTree;++pi) {
		if (!TYPEQueryACTIVE(&(p[pi]))) continue;
		/*
		 ** Do a Gather on q and Q, adjusted for Hill sphere.
		 */
		fq = PARTQQ(&p[pi], 0);
		fQ = PARTQQ(&p[pi], 1);
		fh = 0.5*(fq + fQ)*p[pi].fHill;
		fq -= fh;
		fQ += fh;
		nCnt = smGatherQQ(smx,fq,fQ,ROOT);
		/*
		 ** Start non-local search.
		 */
		cp = ROOT;
		id = -1;	/* We are in the LTT now ! */
		while (1) {
			if (id == pkd->idSelf) goto SkipLocal;
			if (id >= 0) pkdn = mdlAquire(mdl,CID_CELL,cp,id);
			else pkdn = &pkd->kdTop[cp];
			INTERSECTQQ(pkdn,fq,fQ,GetNextCell);
			if (pkdn->iDim >= 0) {
				if (id >= 0) mdlRelease(mdl,CID_CELL,pkdn);
				pkdLower(pkd,cp,id);
				continue;
				}
			for (pj=pkdn->pLower;pj<=pkdn->pUpper;++pj) {
				pPart = mdlAquire(mdl,CID_PARTICLE,pj,id);
				fqp = PARTQQ(pPart, 0);
				fQp = PARTQQ(pPart, 1);
				fhp = 0.5*(fqp + fQp)*pPart->fHill;
				if(fqp - fhp < fQ && fQp + fhp > fq) {
					if(nCnt >= smx->nListSize)
					    smGrowList(smx);
					smx->nnList[nCnt].pPart = pPart;
					smx->nnList[nCnt].iIndex = pj;
					smx->nnList[nCnt].iPid = id;
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

void
CheckForEncounter(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf)
{
	double dInteract(double dTime,double dDelta,double dCentMass,
					 PARTICLE *pi,PARTICLE *pj);

	double min_enc_time = DBL_MAX,enc_time;
	int i;

/*
	p->fDensity = nSmooth;
	return;
*/

/*	(void) printf("CheckForEncounter: p->iOrder = %i, nSmooth = %i\n",
				  p->iOrder,nSmooth);*/
	for (i=0;i<nSmooth;i++) {
		if (p == nnList[i].pPart) continue;
		/* could replace above with "if (p->iOrder >= nnList[i].pPart->iOrder)
		   continue" to save time, assuming multiple encounters are rare */
/*		(void) printf("Checking %i\n",nnList[i].pPart->iOrder);*/
		enc_time = dInteract(smf->dTime,smf->dDelta,smf->dCentMass,p,nnList[i].pPart);
/*		(void) printf("   enc_time = %e\n",enc_time);*/
		if (enc_time == 0) continue;
		if (enc_time < min_enc_time) min_enc_time = enc_time;
		}
    p->fDensity = (min_enc_time == DBL_MAX ? 0 : min_enc_time);
/*	(void) printf("CheckForEncounter: min_enc_time = %e\n",p->fDensity);*/
    }

#endif /* COLLISIONS */
#endif /* OLD_KEPLER */
