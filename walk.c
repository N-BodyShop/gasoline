#include <math.h>
#include <stdlib.h>
#include <stddef.h>
#include <malloc.h>
#include <assert.h>
#include "walk.h"
#include "pkd.h"


void pkdLocalWalk(PKD pkd,int iBucket,float fSoftMax,int bRep,float rOffset[3])
{
	PARTICLE *p;
	KDN *pkdn,*pbuc;
	int iCell,nPart,nCellSoft,nCellNewt,n,pj,bIntersect;
	float x,y,z,twoh2,tr;

	nPart = pkd->nPart;
	nCellSoft = pkd->nCellSoft;
	nCellNewt = pkd->nCellNewt;
	p = pkd->pStore;
	pbuc = &pkd->kdNodes[iBucket];
	iCell = ROOT;
	while (1) {
		pkdn = &pkd->kdNodes[iCell];
		x = pkdn->r[0] + rOffset[0];
		y = pkdn->r[1] + rOffset[1];
		z = pkdn->r[2] + rOffset[2];
		INTERSECTNP(pbuc,pkdn->fOpen2,x,y,z,bIntersect);
		if (bIntersect) {
			/*
			 ** Open cell.
			 */
			if (pkdn->iDim >= 0) {
				iCell = LOWER(iCell);
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
						if (pkd->nMaxPart > pkd->nSqrtTmp) {
						    pkd->nSqrtTmp = pkd->nMaxPart;
						    pkd->sqrttmp = realloc(pkd->sqrttmp,pkd->nSqrtTmp*sizeof(double));
						    pkd->d2a = realloc(pkd->d2a,pkd->nSqrtTmp*sizeof(double));
						}
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
				}
			}
		else {
			/*
			 ** Cell-Bucket interaction accepted.
			 */
			twoh2 = pkdn->fSoft + fSoftMax;
			twoh2 *= twoh2;
			if (twoh2 < pkdn->fOpen2) bIntersect = 0;
			else {
				INTERSECTNP(pbuc,twoh2,x,y,z,bIntersect);
				}
			if (bIntersect) {
				if (nCellSoft == pkd->nMaxCellSoft) {
					pkd->nMaxCellSoft *= 2;
					pkd->ilcs = realloc(pkd->ilcs,pkd->nMaxCellSoft*
										sizeof(ILCS));
					if (pkd->nMaxCellSoft > pkd->nSqrtTmp) {
					    pkd->nSqrtTmp = pkd->nMaxCellSoft;
					    pkd->sqrttmp = realloc(pkd->sqrttmp,pkd->nSqrtTmp*sizeof(double));
					    pkd->d2a = realloc(pkd->d2a,pkd->nSqrtTmp*sizeof(double));
						}
					}
				pkd->ilcs[nCellSoft].m = pkdn->fMass;
				pkd->ilcs[nCellSoft].h = pkdn->fSoft;
				pkd->ilcs[nCellSoft].x = x;
				pkd->ilcs[nCellSoft].y = y;
				pkd->ilcs[nCellSoft].z = z;
				pkd->ilcs[nCellSoft].xx = pkdn->fQxx;
				pkd->ilcs[nCellSoft].yy = pkdn->fQyy;
				pkd->ilcs[nCellSoft].zz = pkdn->fQzz;
				pkd->ilcs[nCellSoft].xy = pkdn->fQxy;
				pkd->ilcs[nCellSoft].xz = pkdn->fQxz;
				pkd->ilcs[nCellSoft].yz = pkdn->fQyz;
				++nCellSoft;
				}
			else {
				if (nCellNewt == pkd->nMaxCellNewt) {
					pkd->nMaxCellNewt *= 2;
					pkd->ilcn = realloc(pkd->ilcn,pkd->nMaxCellNewt*
										sizeof(ILCN));
					if (pkd->nMaxCellNewt > pkd->nSqrtTmp) {
					    pkd->nSqrtTmp = pkd->nMaxCellNewt;
					    pkd->sqrttmp = realloc(pkd->sqrttmp,pkd->nSqrtTmp*sizeof(double));
					    pkd->d2a = realloc(pkd->d2a,pkd->nSqrtTmp*sizeof(double));
						}
					}
				pkd->ilcn[nCellNewt].m = pkdn->fMass;
				pkd->ilcn[nCellNewt].x = x;
				pkd->ilcn[nCellNewt].y = y;
				pkd->ilcn[nCellNewt].z = z;
				tr = pkdn->fQxx + pkdn->fQyy + pkdn->fQzz;
				pkd->ilcn[nCellNewt].xx = 3.0*pkdn->fQxx - tr;
				pkd->ilcn[nCellNewt].yy = 3.0*pkdn->fQyy - tr;
				pkd->ilcn[nCellNewt].xy = 3.0*pkdn->fQxy;
				pkd->ilcn[nCellNewt].xz = 3.0*pkdn->fQxz;
				pkd->ilcn[nCellNewt].yz = 3.0*pkdn->fQyz;
				++nCellNewt;
				}
			SETNEXT(iCell);
			if (iCell == ROOT) break;
			}
		}
	pkd->nPart = nPart;
	pkd->nCellSoft = nCellSoft;
	pkd->nCellNewt = nCellNewt;
	}


void pkdRemoteWalk(PKD pkd,int iBucket,float fSoftMax,int id,float rOffset[3])
{
	PARTICLE *p;
	KDN *pkdn,*pbuc;
	int iCell,nPart,nCellSoft,nCellNewt,n,j,bIntersect;
	float x,y,z,twoh2,tr;

	assert(id != pkd->idSelf);
	nPart = pkd->nPart;
	nCellSoft = pkd->nCellSoft;
	nCellNewt = pkd->nCellNewt;
	pbuc = &pkd->kdNodes[iBucket];
	iCell = ROOT;
	while (1) {
		pkdn = mdlAquire(pkd->mdl,CID_CELL,iCell,id);
		x = pkdn->r[0] + rOffset[0];
		y = pkdn->r[1] + rOffset[1];
		z = pkdn->r[2] + rOffset[2];
		INTERSECTNP(pbuc,pkdn->fOpen2,x,y,z,bIntersect);
		if (bIntersect) {
			/*
			 ** Open cell.
			 */
			if (pkdn->iDim >= 0) {
				mdlRelease(pkd->mdl,CID_CELL,pkdn);
				iCell = LOWER(iCell);
				}
			else {
				/*
				 ** Bucket-Bucket Interaction.
				 */
				n = pkdn->pUpper - pkdn->pLower + 1;
				if (nPart + n > pkd->nMaxPart) {
					pkd->nMaxPart += pkd->nMaxPart + n;
					pkd->ilp = realloc(pkd->ilp,pkd->nMaxPart*sizeof(ILP));
					if (pkd->nMaxPart > pkd->nSqrtTmp) {
					    pkd->nSqrtTmp = pkd->nMaxPart;
					    pkd->sqrttmp = realloc(pkd->sqrttmp,pkd->nSqrtTmp*sizeof(double));
					    pkd->d2a = realloc(pkd->d2a,pkd->nSqrtTmp*sizeof(double));
						}
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
				}
			}
		else {
			/*
			 ** Cell-Bucket interaction accepted.
			 */
			twoh2 = pkdn->fSoft + fSoftMax;
			twoh2 *= twoh2;
			if (twoh2 < pkdn->fOpen2) bIntersect = 0;
			else {
				INTERSECTNP(pbuc,twoh2,x,y,z,bIntersect);
				}
			if (bIntersect) {
				if (nCellSoft == pkd->nMaxCellSoft) {
					pkd->nMaxCellSoft *= 2;
					pkd->ilcs = realloc(pkd->ilcs,pkd->nMaxCellSoft*
										sizeof(ILCS));
					if (pkd->nMaxCellSoft > pkd->nSqrtTmp) {
					    pkd->nSqrtTmp = pkd->nMaxCellSoft;
					    pkd->sqrttmp = realloc(pkd->sqrttmp,pkd->nSqrtTmp*sizeof(double));
					    pkd->d2a = realloc(pkd->d2a,pkd->nSqrtTmp*sizeof(double));
						}
					}
				pkd->ilcs[nCellSoft].m = pkdn->fMass;
				pkd->ilcs[nCellSoft].h = pkdn->fSoft;
				pkd->ilcs[nCellSoft].x = x;
				pkd->ilcs[nCellSoft].y = y;
				pkd->ilcs[nCellSoft].z = z;
				pkd->ilcs[nCellSoft].xx = pkdn->fQxx;
				pkd->ilcs[nCellSoft].yy = pkdn->fQyy;
				pkd->ilcs[nCellSoft].zz = pkdn->fQzz;
				pkd->ilcs[nCellSoft].xy = pkdn->fQxy;
				pkd->ilcs[nCellSoft].xz = pkdn->fQxz;
				pkd->ilcs[nCellSoft].yz = pkdn->fQyz;
				++nCellSoft;
				}
			else {
				if (nCellNewt == pkd->nMaxCellNewt) {
					pkd->nMaxCellNewt *= 2;
					pkd->ilcn = realloc(pkd->ilcn,pkd->nMaxCellNewt*
										sizeof(ILCN));
					}
				pkd->ilcn[nCellNewt].m = pkdn->fMass;
				pkd->ilcn[nCellNewt].x = x;
				pkd->ilcn[nCellNewt].y = y;
				pkd->ilcn[nCellNewt].z = z;
				tr = pkdn->fQxx + pkdn->fQyy + pkdn->fQzz;
				pkd->ilcn[nCellNewt].xx = 3.0*pkdn->fQxx - tr;
				pkd->ilcn[nCellNewt].yy = 3.0*pkdn->fQyy - tr;
				pkd->ilcn[nCellNewt].xy = 3.0*pkdn->fQxy;
				pkd->ilcn[nCellNewt].xz = 3.0*pkdn->fQxz;
				pkd->ilcn[nCellNewt].yz = 3.0*pkdn->fQyz;
				++nCellNewt;
				}
			mdlRelease(pkd->mdl,CID_CELL,pkdn);
			SETNEXT(iCell);
			if (iCell == ROOT) break;
			}
		}
	pkd->nPart = nPart;
	pkd->nCellSoft = nCellSoft;
	pkd->nCellNewt = nCellNewt;
	}


void pkdBucketWalk(PKD pkd,int iBucket,int nReps)
{
	KDN *pbuc,*pkdn;
	int iCell,id,ix,iy,iz,bRep,bIntersect,pj;
	float x,y,z,rOffset[3],fSoftMax,twoh2,tr;
	
	pbuc = &pkd->kdNodes[iBucket];
	/*
	 ** Find the maximum softening for particles in the bucket.
	 */
	fSoftMax = 0.0;
	for (pj=pbuc->pLower;pj<=pbuc->pUpper;++pj) {
		if (pkd->pStore[pj].fSoft > fSoftMax) {
			fSoftMax = pkd->pStore[pj].fSoft;
			}
		}
	pkd->nPart = 0;
	pkd->nCellSoft = 0;
	pkd->nCellNewt = 0;
	for (ix=-nReps;ix<=nReps;++ix) {
		rOffset[0] = ix*pkd->fPeriod[0];
		for (iy=-nReps;iy<=nReps;++iy) {
			rOffset[1] = iy*pkd->fPeriod[1];
			for (iz=-nReps;iz<=nReps;++iz) {
				rOffset[2] = iz*pkd->fPeriod[2];
				bRep = ix || iy || iz;
				/*
				 ** Walk the top tree first, finding local trees to
				 ** continue walking.
				 */
				iCell = ROOT;
				while (1) {
					id = pkd->kdTop[iCell].pLower;
					if (id == pkd->idSelf) {
						pkdLocalWalk(pkd,iBucket,fSoftMax,bRep,rOffset);
						SETNEXT(iCell);
						if (iCell == ROOT) break;
						}
					else if (id >= 0) {
						pkdRemoteWalk(pkd,iBucket,fSoftMax,id,rOffset);
						SETNEXT(iCell);
						if (iCell == ROOT) break;
						}
					else {
						pkdn = &pkd->kdTop[iCell];
						x = pkdn->r[0] + rOffset[0];
						y = pkdn->r[1] + rOffset[1];
						z = pkdn->r[2] + rOffset[2];
						INTERSECTNP(pbuc,pkdn->fOpen2,x,y,z,bIntersect);
						if (bIntersect) {
							/*
							 ** Open top cell.
							 */
							iCell = LOWER(iCell);
							}
						else {
							/*
							 ** Cell-Bucket interaction accepted.
							 ** Decide whether it is safe to use
							 ** Newtonian form.
							 */
							twoh2 = pkdn->fSoft + fSoftMax;
							twoh2 *= twoh2;
							if (twoh2 < pkdn->fOpen2) bIntersect = 0;
							else {
								INTERSECTNP(pbuc,twoh2,x,y,z,bIntersect);
								}
							if (bIntersect) {
								/*
								 ** May need to use softened Cell!
								 */
								if (pkd->nCellSoft == pkd->nMaxCellSoft) {
									pkd->nMaxCellSoft *= 2;
									pkd->ilcs = realloc(pkd->ilcs,
														pkd->nMaxCellSoft*
														sizeof(ILCS));
									if (pkd->nMaxCellSoft > pkd->nSqrtTmp) {
									    pkd->nSqrtTmp = pkd->nMaxCellSoft;
									    pkd->sqrttmp = realloc(pkd->sqrttmp,pkd->nSqrtTmp*sizeof(double));
									    pkd->d2a = realloc(pkd->d2a,pkd->nSqrtTmp*sizeof(double));
										}
									}
								pkd->ilcs[pkd->nCellSoft].m = pkdn->fMass;
								pkd->ilcs[pkd->nCellSoft].h = pkdn->fSoft;
								pkd->ilcs[pkd->nCellSoft].x = x;
								pkd->ilcs[pkd->nCellSoft].y = y;
								pkd->ilcs[pkd->nCellSoft].z = z;
								pkd->ilcs[pkd->nCellSoft].xx = pkdn->fQxx;
								pkd->ilcs[pkd->nCellSoft].yy = pkdn->fQyy;
								pkd->ilcs[pkd->nCellSoft].zz = pkdn->fQzz;
								pkd->ilcs[pkd->nCellSoft].xy = pkdn->fQxy;
								pkd->ilcs[pkd->nCellSoft].xz = pkdn->fQxz;
								pkd->ilcs[pkd->nCellSoft].yz = pkdn->fQyz;
								++pkd->nCellSoft;
								}
							else {
								/*
								 ** Can use Newtonian Cell interaction.
								 */
								if (pkd->nCellNewt == pkd->nMaxCellNewt) {
									pkd->nMaxCellNewt *= 2;
									pkd->ilcn = realloc(pkd->ilcn,
														pkd->nMaxCellNewt*
														sizeof(ILCN));
									if (pkd->nMaxCellNewt > pkd->nSqrtTmp) {
									    pkd->nSqrtTmp = pkd->nMaxCellNewt;
									    pkd->sqrttmp = realloc(pkd->sqrttmp,pkd->nSqrtTmp*sizeof(double));
									    pkd->d2a = realloc(pkd->d2a,pkd->nSqrtTmp*sizeof(double));
										}
									}
								pkd->ilcn[pkd->nCellNewt].m = pkdn->fMass;
								pkd->ilcn[pkd->nCellNewt].x = x;
								pkd->ilcn[pkd->nCellNewt].y = y;
								pkd->ilcn[pkd->nCellNewt].z = z;
								tr = pkdn->fQxx + pkdn->fQyy + pkdn->fQzz;
								pkd->ilcn[pkd->nCellNewt].xx = 
									3.0*pkdn->fQxx - tr;
								pkd->ilcn[pkd->nCellNewt].yy = 
									3.0*pkdn->fQyy - tr;
								pkd->ilcn[pkd->nCellNewt].xy = 3.0*pkdn->fQxy;
								pkd->ilcn[pkd->nCellNewt].xz = 3.0*pkdn->fQxz;
								pkd->ilcn[pkd->nCellNewt].yz = 3.0*pkdn->fQyz;
								++pkd->nCellNewt;
								}
							SETNEXT(iCell);
							if (iCell == ROOT) break;
							}
						}
					} /* of Top tree walk */
				} /* of iz */
			} /* of iy */
		} /* of ix */
	}



