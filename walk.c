#include <math.h>
#include <stdlib.h>
#include <stddef.h>
#include <malloc.h>
#include <assert.h>
#include "walk.h"
#include "pkd.h"


#ifdef COMPLETE_LOCAL
#define SETILIST(iOrder,ilcn,pkdn,x,y,z)\
{\
	switch (iOrder) {\
	case 4:\
		ilcn.xxxx = pkdn->mom.Hxxxx;\
		ilcn.xyyy = pkdn->mom.Hxyyy;\
		ilcn.xxxy = pkdn->mom.Hxxxy;\
		ilcn.yyyy = pkdn->mom.Hyyyy;\
		ilcn.xxxz = pkdn->mom.Hxxxz;\
		ilcn.yyyz = pkdn->mom.Hyyyz;\
		ilcn.xxyy = pkdn->mom.Hxxyy;\
		ilcn.xxyz = pkdn->mom.Hxxyz;\
		ilcn.xyyz = pkdn->mom.Hxyyz;\
		ilcn.xxzz = pkdn->mom.Hxxzz;\
		ilcn.xyzz = pkdn->mom.Hxyzz;\
		ilcn.xzzz = pkdn->mom.Hxzzz;\
		ilcn.yyzz = pkdn->mom.Hyyzz;\
		ilcn.yzzz = pkdn->mom.Hyzzz;\
		ilcn.zzzz = pkdn->mom.Hzzzz;\
	case 3:\
		ilcn.xxx = pkdn->mom.Oxxx;\
		ilcn.xyy = pkdn->mom.Oxyy;\
		ilcn.xxy = pkdn->mom.Oxxy;\
		ilcn.yyy = pkdn->mom.Oyyy;\
		ilcn.xxz = pkdn->mom.Oxxz;\
		ilcn.yyz = pkdn->mom.Oyyz;\
		ilcn.xyz = pkdn->mom.Oxyz;\
		ilcn.xzz = pkdn->mom.Oxzz;\
		ilcn.yzz = pkdn->mom.Oyzz;\
		ilcn.zzz = pkdn->mom.Ozzz;\
	case 2:\
		ilcn.xx = pkdn->mom.Qxx;\
		ilcn.yy = pkdn->mom.Qyy;\
                ilcn.zz = pkdn->mom.Qzz;\
		ilcn.xy = pkdn->mom.Qxy;\
		ilcn.xz = pkdn->mom.Qxz;\
		ilcn.yz = pkdn->mom.Qyz;\
	case 1:\
	default:\
		ilcn.m = pkdn->fMass;\
		ilcn.x = x;\
		ilcn.y = y;\
		ilcn.z = z;\
		}\
	}
#else
#define SETILIST(iOrder,ilcn,pkdn,x,y,z)\
{\
	switch (iOrder) {\
        double tr;\
	case 4:\
		ilcn.xxxx = pkdn->mom.Hxxxx;\
		ilcn.xyyy = pkdn->mom.Hxyyy;\
		ilcn.xxxy = pkdn->mom.Hxxxy;\
		ilcn.yyyy = pkdn->mom.Hyyyy;\
		ilcn.xxxz = pkdn->mom.Hxxxz;\
		ilcn.yyyz = pkdn->mom.Hyyyz;\
		ilcn.xxyy = pkdn->mom.Hxxyy;\
		ilcn.xxyz = pkdn->mom.Hxxyz;\
		ilcn.xyyz = pkdn->mom.Hxyyz;\
		ilcn.xxzz = pkdn->mom.Hxxzz;\
		ilcn.xyzz = pkdn->mom.Hxyzz;\
		ilcn.xzzz = pkdn->mom.Hxzzz;\
		ilcn.yyzz = pkdn->mom.Hyyzz;\
		ilcn.yzzz = pkdn->mom.Hyzzz;\
		ilcn.zzzz = pkdn->mom.Hzzzz;\
	case 3:\
		ilcn.xxx = pkdn->mom.Oxxx;\
		ilcn.xyy = pkdn->mom.Oxyy;\
		ilcn.xxy = pkdn->mom.Oxxy;\
		ilcn.yyy = pkdn->mom.Oyyy;\
		ilcn.xxz = pkdn->mom.Oxxz;\
		ilcn.yyz = pkdn->mom.Oyyz;\
		ilcn.xyz = pkdn->mom.Oxyz;\
		ilcn.xzz = pkdn->mom.Oxzz;\
		ilcn.yzz = pkdn->mom.Oyzz;\
		ilcn.zzz = pkdn->mom.Ozzz;\
	case 2:\
		tr = pkdn->mom.Qxx + pkdn->mom.Qyy + pkdn->mom.Qzz;\
		ilcn.xx = pkdn->mom.Qxx - tr/3.0;\
		ilcn.yy = pkdn->mom.Qyy - tr/3.0;\
                ilcn.zz = pkdn->mom.Qzz - tr/3.0;\
		ilcn.xy = pkdn->mom.Qxy;\
		ilcn.xz = pkdn->mom.Qxz;\
		ilcn.yz = pkdn->mom.Qyz;\
	case 1:\
	default:\
		ilcn.m = pkdn->fMass;\
		ilcn.x = x;\
		ilcn.y = y;\
		ilcn.z = z;\
		}\
	}
#endif

void pkdLocalWalk(PKD pkd,int iBucket,FLOAT fSoftMax,int bRep,FLOAT rOffset[3],
				  int iOrder)
{
	PARTICLE *p;
	KDN *pkdn,*pbuc;
	int iCell,nPart,nCellSoft,nCellNewt,n,pj,bIntersect;
	FLOAT x,y,z,twoh2;

	nPart = pkd->nPart;
	nCellSoft = pkd->nCellSoft;
	nCellNewt = pkd->nCellNewt;
	p = pkd->pStore;
	pbuc = &pkd->kdNodes[iBucket];
	iCell = pkd->iRoot;
	while (iCell != -1) {
		pkdn = &pkd->kdNodes[iCell];
		x = pkdn->r[0] + rOffset[0];
		y = pkdn->r[1] + rOffset[1];
		z = pkdn->r[2] + rOffset[2];
		INTERSECTNP(pbuc,pkdn->fOpen2,x,y,z,bIntersect);
		/*
		 ** If the cell has less than 4 particles in it then open it.
		 */
		if (pkdn->pUpper - pkdn->pLower < 3) bIntersect = 1;
		if (bIntersect) {
			if (pkdn->iLower != -1) {
				/*
				 ** Open cell.
				 */
				iCell = pkdn->iLower;
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
						assert(pkd->ilp != NULL);
						if (pkd->nMaxPart > pkd->nSqrtTmp) {
						    pkd->nSqrtTmp = pkd->nMaxPart;
						    pkd->sqrttmp = realloc(pkd->sqrttmp,pkd->nSqrtTmp*sizeof(double));
							assert(pkd->sqrttmp != NULL);
						    pkd->d2a = realloc(pkd->d2a,pkd->nSqrtTmp*sizeof(double));
							assert(pkd->d2a != NULL);
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
				iCell = pkdn->iUpper;
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
					assert(pkd->ilcs != NULL);
					if (pkd->nMaxCellSoft > pkd->nSqrtTmp) {
					    pkd->nSqrtTmp = pkd->nMaxCellSoft;
					    pkd->sqrttmp = realloc(pkd->sqrttmp,pkd->nSqrtTmp*sizeof(double));
						assert(pkd->sqrttmp != NULL);
					    pkd->d2a = realloc(pkd->d2a,pkd->nSqrtTmp*sizeof(double));
						assert(pkd->d2a != NULL);
						}
					}
				pkd->ilcs[nCellSoft].m = pkdn->fMass;
				pkd->ilcs[nCellSoft].h = pkdn->fSoft;
				pkd->ilcs[nCellSoft].x = x;
				pkd->ilcs[nCellSoft].y = y;
				pkd->ilcs[nCellSoft].z = z;
				pkd->ilcs[nCellSoft].xx = pkdn->mom.Qxx;
				pkd->ilcs[nCellSoft].yy = pkdn->mom.Qyy;
				pkd->ilcs[nCellSoft].zz = pkdn->mom.Qzz;
				pkd->ilcs[nCellSoft].xy = pkdn->mom.Qxy;
				pkd->ilcs[nCellSoft].xz = pkdn->mom.Qxz;
				pkd->ilcs[nCellSoft].yz = pkdn->mom.Qyz;
				++nCellSoft;
				}
			else {
				if (nCellNewt == pkd->nMaxCellNewt) {
					pkd->nMaxCellNewt *= 2;
					pkd->ilcn = realloc(pkd->ilcn,pkd->nMaxCellNewt*
										sizeof(ILCN));
					assert(pkd->ilcn != NULL);
					if (pkd->nMaxCellNewt > pkd->nSqrtTmp) {
					    pkd->nSqrtTmp = pkd->nMaxCellNewt;
					    pkd->sqrttmp = realloc(pkd->sqrttmp,pkd->nSqrtTmp*sizeof(double));
						assert(pkd->sqrttmp != NULL);
					    pkd->d2a = realloc(pkd->d2a,pkd->nSqrtTmp*sizeof(double));
						assert(pkd->d2a != NULL);
						}
					}
				SETILIST(iOrder,pkd->ilcn[nCellNewt],pkdn,x,y,z);
				++nCellNewt;
				}
			iCell = pkdn->iUpper;
			}
		}
	pkd->nPart = nPart;
	pkd->nCellSoft = nCellSoft;
	pkd->nCellNewt = nCellNewt;
	}


void pkdRemoteWalk(PKD pkd,int iBucket,FLOAT fSoftMax,int id,FLOAT rOffset[3],
				   int iOrder)
{
	PARTICLE *p;
	KDN *pkdn,*pbuc;
	int iCell,nPart,nCellSoft,nCellNewt,n,j,bIntersect;
	FLOAT x,y,z,twoh2;

	assert(id != pkd->idSelf);
	nPart = pkd->nPart;
	nCellSoft = pkd->nCellSoft;
	nCellNewt = pkd->nCellNewt;
	pbuc = &pkd->kdNodes[iBucket];
	iCell = pkd->iRoot;
	while (iCell != -1) {
		pkdn = mdlAquire(pkd->mdl,CID_CELL,iCell,id);
		x = pkdn->r[0] + rOffset[0];
		y = pkdn->r[1] + rOffset[1];
		z = pkdn->r[2] + rOffset[2];
		INTERSECTNP(pbuc,pkdn->fOpen2,x,y,z,bIntersect);
		/*
		 ** If the cell has less than 4 particles in it then open it.
		 */
		if (pkdn->pUpper - pkdn->pLower < 3) bIntersect = 1;
		if (bIntersect) {
			if (pkdn->iLower != -1) {
				/*
				 ** Open cell.
				 */
				iCell = pkdn->iLower;
				mdlRelease(pkd->mdl,CID_CELL,pkdn);
				}
			else {
				/*
				 ** Bucket-Bucket Interaction.
				 */
				n = pkdn->pUpper - pkdn->pLower + 1;
				if (nPart + n > pkd->nMaxPart) {
					pkd->nMaxPart += pkd->nMaxPart + n;
					pkd->ilp = realloc(pkd->ilp,pkd->nMaxPart*sizeof(ILP));
					assert(pkd->ilp != NULL);
					if (pkd->nMaxPart > pkd->nSqrtTmp) {
					    pkd->nSqrtTmp = pkd->nMaxPart;
					    pkd->sqrttmp = realloc(pkd->sqrttmp,pkd->nSqrtTmp*sizeof(double));
						assert(pkd->sqrttmp != NULL);
					    pkd->d2a = realloc(pkd->d2a,pkd->nSqrtTmp*sizeof(double));
						assert(pkd->d2a != NULL);
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
				iCell = pkdn->iUpper;
				mdlRelease(pkd->mdl,CID_CELL,pkdn);
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
					assert(pkd->ilcs != NULL);
					if (pkd->nMaxCellSoft > pkd->nSqrtTmp) {
					    pkd->nSqrtTmp = pkd->nMaxCellSoft;
					    pkd->sqrttmp = realloc(pkd->sqrttmp,pkd->nSqrtTmp*sizeof(double));
						assert(pkd->sqrttmp != NULL);
					    pkd->d2a = realloc(pkd->d2a,pkd->nSqrtTmp*sizeof(double));
						assert(pkd->d2a != NULL);
						}
					}
				pkd->ilcs[nCellSoft].m = pkdn->fMass;
				pkd->ilcs[nCellSoft].h = pkdn->fSoft;
				pkd->ilcs[nCellSoft].x = x;
				pkd->ilcs[nCellSoft].y = y;
				pkd->ilcs[nCellSoft].z = z;
				pkd->ilcs[nCellSoft].xx = pkdn->mom.Qxx;
				pkd->ilcs[nCellSoft].yy = pkdn->mom.Qyy;
				pkd->ilcs[nCellSoft].zz = pkdn->mom.Qzz;
				pkd->ilcs[nCellSoft].xy = pkdn->mom.Qxy;
				pkd->ilcs[nCellSoft].xz = pkdn->mom.Qxz;
				pkd->ilcs[nCellSoft].yz = pkdn->mom.Qyz;
				++nCellSoft;
				}
			else {
				if (nCellNewt == pkd->nMaxCellNewt) {
					pkd->nMaxCellNewt *= 2;
					pkd->ilcn = realloc(pkd->ilcn,pkd->nMaxCellNewt*
										sizeof(ILCN));
					assert(pkd->ilcn != NULL);
					if (pkd->nMaxCellNewt > pkd->nSqrtTmp) {
					    pkd->nSqrtTmp = pkd->nMaxCellNewt;
					    pkd->sqrttmp = realloc(pkd->sqrttmp,pkd->nSqrtTmp*sizeof(double));
						assert(pkd->sqrttmp != NULL);
					    pkd->d2a = realloc(pkd->d2a,pkd->nSqrtTmp*sizeof(double));
						assert(pkd->d2a != NULL);
						}
					}
				SETILIST(iOrder,pkd->ilcn[nCellNewt],pkdn,x,y,z);
				++nCellNewt;
				}
			iCell = pkdn->iUpper;
			mdlRelease(pkd->mdl,CID_CELL,pkdn);
			}
		}
	pkd->nPart = nPart;
	pkd->nCellSoft = nCellSoft;
	pkd->nCellNewt = nCellNewt;
	}


void pkdBucketWalk(PKD pkd,int iBucket,int nReps,int iOrder)
{
	KDN *pbuc,*pkdn;
	int iCell,id,ix,iy,iz,bRep,bIntersect,pj;
	FLOAT x,y,z,rOffset[3],fSoftMax,twoh2;
	
	pkd->nPart = 0;
	pkd->nCellSoft = 0;
	pkd->nCellNewt = 0;
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
						pkdLocalWalk(pkd,iBucket,fSoftMax,bRep,rOffset,iOrder);
						SETNEXT(iCell);
						if (iCell == ROOT) break;
						}
					else if (id >= 0) {
						pkdRemoteWalk(pkd,iBucket,fSoftMax,id,rOffset,iOrder);
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
									assert(pkd->ilcs != NULL);
									if (pkd->nMaxCellSoft > pkd->nSqrtTmp) {
									    pkd->nSqrtTmp = pkd->nMaxCellSoft;
									    pkd->sqrttmp = realloc(pkd->sqrttmp,pkd->nSqrtTmp*sizeof(double));
										assert(pkd->sqrttmp != NULL);
									    pkd->d2a = realloc(pkd->d2a,pkd->nSqrtTmp*sizeof(double));
										assert(pkd->d2a != NULL);
										}
									}
								pkd->ilcs[pkd->nCellSoft].m = pkdn->fMass;
								pkd->ilcs[pkd->nCellSoft].h = pkdn->fSoft;
								pkd->ilcs[pkd->nCellSoft].x = x;
								pkd->ilcs[pkd->nCellSoft].y = y;
								pkd->ilcs[pkd->nCellSoft].z = z;
								pkd->ilcs[pkd->nCellSoft].xx = pkdn->mom.Qxx;
								pkd->ilcs[pkd->nCellSoft].yy = pkdn->mom.Qyy;
								pkd->ilcs[pkd->nCellSoft].zz = pkdn->mom.Qzz;
								pkd->ilcs[pkd->nCellSoft].xy = pkdn->mom.Qxy;
								pkd->ilcs[pkd->nCellSoft].xz = pkdn->mom.Qxz;
								pkd->ilcs[pkd->nCellSoft].yz = pkdn->mom.Qyz;
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
									assert(pkd->ilcn != NULL);
									if (pkd->nMaxCellNewt > pkd->nSqrtTmp) {
									    pkd->nSqrtTmp = pkd->nMaxCellNewt;
									    pkd->sqrttmp = realloc(pkd->sqrttmp,pkd->nSqrtTmp*sizeof(double));
										assert(pkd->sqrttmp != NULL);
									    pkd->d2a = realloc(pkd->d2a,pkd->nSqrtTmp*sizeof(double));
										assert(pkd->d2a != NULL);
										}
									}
								SETILIST(iOrder,pkd->ilcn[pkd->nCellNewt],pkdn,x,y,z);
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



