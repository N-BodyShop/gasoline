#include <math.h>
#include <stdlib.h>
#include <stddef.h>
#include <malloc.h>
#include <assert.h>
#include "walk.h"
#include "pkd.h"


void pkdLocalWalk(PKD pkd,int iBucket,int bRep,float rOffset[3])
{
	PARTICLE *p;
	KDN *pkdn,*pbuc;
	int iCell,nPart,nCell,n,pj;
	float x,y,z;

	nPart = pkd->nPart;
	nCell = pkd->nCell;
	p = pkd->pStore;
	pbuc = &pkd->kdNodes[iBucket];
	iCell = ROOT;
	while (1) {
		pkdn = &pkd->kdNodes[iCell];
		x = pkdn->r[0] + rOffset[0];
		y = pkdn->r[1] + rOffset[1];
		z = pkdn->r[2] + rOffset[2];
		INTERSECTNP(pbuc,pkdn->fOpen2,x,y,z);
		/*
		 ** Open cell.
		 */
		if (pkdn->iDim >= 0) {
			iCell = LOWER(iCell);
			continue;
			}
		else {
			/*
			 ** Bucket-Bucket Interaction.
			 */
			if (iCell != iBucket || bRep) {
				n = pkdn->pUpper - pkdn->pLower + 1;
				if (nPart + n > pkd->nMaxPart) {
					pkd->nMaxPart += pkd->nMaxPart + n;
					pkd->ilp = realloc(pkd->ilp,pkd->nMaxPart*sizeof(ILP));
					}
				for (pj=pkdn->pLower;pj<=pkdn->pUpper;++pj,++nPart) {
					pkd->ilp[nPart].x = p[pj].r[0] + rOffset[0];
					pkd->ilp[nPart].y = p[pj].r[1] + rOffset[1];
					pkd->ilp[nPart].z = p[pj].r[2] + rOffset[2];
					pkd->ilp[nPart].m = p[pj].fMass;
					pkd->ilp[nPart].h = p[pj].fSoft;
					}
				}
			SETNEXT(iCell);
			if (iCell == ROOT) break;
			continue;
			}
	GetNextCell:
		/*
		 ** Cell-Bucket interaction accepted.
		 */
		if (nCell == pkd->nMaxCell) {
			pkd->nMaxCell *= 2;
			pkd->ilc = realloc(pkd->ilc,pkd->nMaxCell*sizeof(ILC));
			}
		pkd->ilc[nCell].xx = pkdn->fQxx;
		pkd->ilc[nCell].yy = pkdn->fQyy;
		pkd->ilc[nCell].zz = pkdn->fQzz;
		pkd->ilc[nCell].xy = pkdn->fQxy;
		pkd->ilc[nCell].xz = pkdn->fQxz;
		pkd->ilc[nCell].yz = pkdn->fQyz;
		pkd->ilc[nCell].x = x;
		pkd->ilc[nCell].y = y;
		pkd->ilc[nCell].z = z;
		pkd->ilc[nCell].m = pkdn->fMass;
		pkd->ilc[nCell].h = pkdn->fSoft;
		++nCell;
		SETNEXT(iCell);
		if (iCell == ROOT) break;
		}
	pkd->nPart = nPart;
	pkd->nCell = nCell;
	}


void pkdRemoteWalk(PKD pkd,int iBucket,int id,float rOffset[3])
{
	PARTICLE *p;
	KDN *pkdn,*pbuc;
	int iCell,nPart,nCell,n,j;
	float x,y,z;

	assert(id != pkd->idSelf);
	nPart = pkd->nPart;
	nCell = pkd->nCell;
	pbuc = &pkd->kdNodes[iBucket];
	iCell = ROOT;
	while (1) {
		pkdn = mdlAquire(pkd->mdl,CID_CELL,iCell,id);
		x = pkdn->r[0] + rOffset[0];
		y = pkdn->r[1] + rOffset[1];
		z = pkdn->r[2] + rOffset[2];
		INTERSECTNP(pbuc,pkdn->fOpen2,x,y,z);
		/*
		 ** Open cell.
		 */
		if (pkdn->iDim >= 0) {
			mdlRelease(pkd->mdl,CID_CELL,pkdn);
			iCell = LOWER(iCell);
			continue;
			}
		else {
			/*
			 ** Bucket-Bucket Interaction.
			 */
			n = pkdn->pUpper - pkdn->pLower + 1;
			if (nPart + n > pkd->nMaxPart) {
				pkd->nMaxPart += pkd->nMaxPart + n;
				pkd->ilp = realloc(pkd->ilp,pkd->nMaxPart*sizeof(ILP));
				}
			for (j=0;j<n;++j,++nPart) {
				p = mdlAquire(pkd->mdl,CID_PARTICLE,pkdn->pLower+j,id);
				pkd->ilp[nPart].x = p->r[0] + rOffset[0];
				pkd->ilp[nPart].y = p->r[1] + rOffset[1];
				pkd->ilp[nPart].z = p->r[2] + rOffset[2];
				pkd->ilp[nPart].m = p->fMass;
				pkd->ilp[nPart].h = p->fSoft;
				mdlRelease(pkd->mdl,CID_PARTICLE,p);
				}
			mdlRelease(pkd->mdl,CID_CELL,pkdn);
			SETNEXT(iCell);
			if (iCell == ROOT) break;
			continue;
			}
	GetNextCell:
		/*
		 ** Cell-Bucket interaction accepted.
		 */
		if (nCell == pkd->nMaxCell) {
			pkd->nMaxCell *= 2;
			pkd->ilc = realloc(pkd->ilc,pkd->nMaxCell*sizeof(ILC));
			}
		pkd->ilc[nCell].xx = pkdn->fQxx;
		pkd->ilc[nCell].yy = pkdn->fQyy;
		pkd->ilc[nCell].zz = pkdn->fQzz;
		pkd->ilc[nCell].xy = pkdn->fQxy;
		pkd->ilc[nCell].xz = pkdn->fQxz;
		pkd->ilc[nCell].yz = pkdn->fQyz;
		pkd->ilc[nCell].x = x;
		pkd->ilc[nCell].y = y;
		pkd->ilc[nCell].z = z;
		pkd->ilc[nCell].m = pkdn->fMass;
		pkd->ilc[nCell].h = pkdn->fSoft;
	    mdlRelease(pkd->mdl,CID_CELL,pkdn);
		++nCell;
		SETNEXT(iCell);
		if (iCell == ROOT) break;
		}
	pkd->nPart = nPart;
	pkd->nCell = nCell;
	}


void pkdBucketWalk(PKD pkd,int iBucket,int nReps)
{
	KDN *pbuc,*pkdn;
	int iCell,id,ix,iy,iz,bRep;
	float x,y,z,rOffset[3];
	
	pbuc = &pkd->kdNodes[iBucket];
	pkd->nPart = 0;
	pkd->nCell = 0;
	for (ix=-nReps;ix<=nReps;++ix) {
		rOffset[0] = ix*pkd->fPeriod[0];
		for (iy=-nReps;iy<=nReps;++iy) {
			rOffset[1] = iy*pkd->fPeriod[1];
			for (iz=-nReps;iz<=nReps;++iz) {
				rOffset[2] = iz*pkd->fPeriod[2];
				bRep = ix && iy && iz;
				/*
				 ** Walk the top tree first, finding local trees to
				 ** continue walking.
				 */
				iCell = ROOT;
				while (1) {
					id = pkd->kdTop[iCell].pLower;
					if (id == pkd->idSelf) {
						pkdLocalWalk(pkd,iBucket,bRep,rOffset);
						}
					else if (id >= 0) {
						pkdRemoteWalk(pkd,iBucket,id,rOffset);
						}
					else {
						pkdn = &pkd->kdTop[iCell];
						x = pkdn->r[0] + rOffset[0];
						y = pkdn->r[1] + rOffset[1];
						z = pkdn->r[2] + rOffset[2];
						INTERSECTNP(pbuc,pkdn->fOpen2,x,y,z);
						/*
						 ** Open top cell.
						 */
						iCell = LOWER(iCell);
						continue;
					GetNextCell:
						/*
						 ** Cell-Bucket interaction accepted.
						 */
						if (pkd->nCell == pkd->nMaxCell) {
							pkd->nMaxCell *= 2;
							pkd->ilc = realloc(pkd->ilc,pkd->nMaxCell*
											   sizeof(ILC));
							}
						pkd->ilc[pkd->nCell].xx = pkdn->fQxx;
						pkd->ilc[pkd->nCell].yy = pkdn->fQyy;
						pkd->ilc[pkd->nCell].zz = pkdn->fQzz;
						pkd->ilc[pkd->nCell].xy = pkdn->fQxy;
						pkd->ilc[pkd->nCell].xz = pkdn->fQxz;
						pkd->ilc[pkd->nCell].yz = pkdn->fQyz;
						pkd->ilc[pkd->nCell].x = x;
						pkd->ilc[pkd->nCell].y = y;
						pkd->ilc[pkd->nCell].z = z;
						pkd->ilc[pkd->nCell].m = pkdn->fMass;
						pkd->ilc[pkd->nCell].h = pkdn->fSoft;
						++pkd->nCell;
						}
					SETNEXT(iCell);
					if (iCell == ROOT) break;
					} /* of Top tree walk */
				} /* of iz */
			} /* of iy */
		} /* of ix */
	}



