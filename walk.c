#include <math.h>
#include <stdlib.h>
#include <stddef.h>
#include <malloc.h>
#include "walk.h"
#include "pkd.h"


void pkdBucketWalk(PKD pkd,int iBucket,int nReps)
{
	MDL mdl;
	PARTICLE *p;
	KDN *pbuc,*pkdn;
	int nCell,nPart,iCell,id,j,n,ix,iy,iz,bRep;
	float x,y,z,xOffset,yOffset,zOffset;
	
	mdl = pkd->mdl;
	pbuc = &pkd->kdNodes[iBucket];
	nPart = 0;
	nCell = 0;
	for (ix=-nReps;ix<=nReps;++ix) {
		xOffset = ix*pkd->fPeriod[0];
		for (iy=-nReps;iy<=nReps;++iy) {
			yOffset = iy*pkd->fPeriod[1];
			for (iz=-nReps;iz<=nReps;++iz) {
				zOffset = iz*pkd->fPeriod[2];
				bRep = ix && iy && iz;
				iCell = ROOT;
				id = -1;
				while (1) {
					if (id >= 0) pkdn = mdlAquire(mdl,CID_CELL,iCell,id);
					else pkdn = &pkd->kdTop[iCell];
					x = pkdn->r[0] + xOffset;
					y = pkdn->r[1] + yOffset;
					z = pkdn->r[2] + zOffset;
					INTERSECTNP(pbuc,pkdn->fOpen2,x,y,z);
					/*
					 ** Open cell.
					 */
					if (pkdn->iDim >= 0) {
						if (id >= 0) mdlRelease(mdl,CID_CELL,pkdn);
						pkdLower(pkd,iCell,id);
						continue;
						}
					else {
						/*
						 ** Bucket-Bucket Interaction.
						 */
						if (id != pkd->idSelf || iCell != iBucket || bRep) {
							n = pkdn->pUpper - pkdn->pLower + 1;
							if (nPart + n > pkd->nMaxPart) {
								pkd->nMaxPart += pkd->nMaxPart + n;
								pkd->ilp = realloc(pkd->ilp,pkd->nMaxPart*
												   sizeof(ILP));
								}
							for (j=0;j<n;++j) {
								p = mdlAquire(mdl,CID_PARTICLE,
											  pkdn->pLower+j,id);
								pkd->ilp[nPart+j].x = p->r[0] + xOffset;
								pkd->ilp[nPart+j].y = p->r[1] + yOffset;
								pkd->ilp[nPart+j].z = p->r[2] + zOffset;
								pkd->ilp[nPart+j].m = p->fMass;
								pkd->ilp[nPart+j].h = p->fSoft;
								mdlRelease(mdl,CID_PARTICLE,p);
								}
							nPart += n;
							}  
						if (id >= 0) mdlRelease(mdl,CID_CELL,pkdn);
						pkdNext(pkd,iCell,id);
						if (pkdIsRoot(iCell,id)) break;
						continue;
						}
				GetNextCell:
					/*
					 ** Cell-Bucket interaction accepted.
					 */
					if (nCell == pkd->nMaxCell) {
						pkd->nMaxCell *= 2;
						pkd->ilc = realloc(pkd->ilc,pkd->nMaxCell*
										   sizeof(ILC));
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
					if (id >= 0) mdlRelease(mdl,CID_CELL,pkdn);
					++nCell;
					pkdNext(pkd,iCell,id);
					if (pkdIsRoot(iCell,id)) break;
					} /* of walk */
				} /* of iz */
			} /* of iy */
		} /* of ix */
	pkd->nPart = nPart;
	pkd->nCell = nCell;
	}



