#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <malloc.h>
#include <math.h>
#include <assert.h>
#include <sys/time.h>

#ifndef CRAY_T3D
#include <rpc/types.h>
#include <rpc/xdr.h>
#endif

#include "pkd.h"
#include "ewald.h"
#include "grav.h"
#include "walk.h"
#include "opentype.h"
#include "mdl.h"
#include "tipsydefs.h"
#include "checkdefs.h"


double pkdGetTimer(PKD pkd,int iTimer)
{
	return(pkd->ti[iTimer].sec);
	}


void pkdClearTimer(PKD pkd,int iTimer)
{
	int i;

	if (iTimer >= 0) {
		pkd->ti[iTimer].sec = 0.0;
		}
	else {
		for (i=0;i<MAX_TIMERS;++i) {
			pkd->ti[i].sec = 0.0;
			}
		}
	}


void pkdStartTimer(PKD pkd,int iTimer)
{
	pkd->ti[iTimer].stamp = mdlCpuTimer(pkd->mdl);
	}


void pkdStopTimer(PKD pkd,int iTimer)
{
	double sec;
	sec = mdlCpuTimer(pkd->mdl) - pkd->ti[iTimer].stamp;
	if (sec < 0.0) sec = 0.0;
	pkd->ti[iTimer].sec += sec;
	}


void pkdInitialize(PKD *ppkd,MDL mdl,int iOrder,int nStore,int nLvl,
				   float *fPeriod)
{
	PKD pkd;
	int j;
	
	pkd = (PKD)malloc(sizeof(struct pkdContext));
	assert(pkd != NULL);
	pkd->mdl = mdl;
	pkd->iOrder = iOrder;
	pkd->idSelf = mdlSelf(mdl);
	pkd->nThreads = mdlThreads(mdl);
	pkd->nStore = nStore;
	pkd->nRejects = 0;
	for (j=0;j<3;++j) {
		pkd->fPeriod[j] = fPeriod[j];
		}
	/*
	 ** Allocate the main particle store.
	 */
	pkd->pStore = (PARTICLE *)malloc(nStore*sizeof(PARTICLE));
	assert(pkd->pStore != NULL);
	pkd->kdNodes = NULL;
	pkd->piLeaf = NULL;
	pkd->kdTop = NULL;
	/*
	 ** Allocate initial interaction lists
	 */
	pkd->nMaxPart = 500;
	pkd->nMaxCellSoft = 500;
	pkd->nMaxCellNewt = 500;
	pkd->nSqrtTmp = 500;
	pkd->ilp = malloc(pkd->nMaxPart*sizeof(ILP));
	assert(pkd->ilp != NULL);
	pkd->ilcs = malloc(pkd->nMaxCellSoft*sizeof(ILCS));
	assert(pkd->ilcs != NULL);
	pkd->ilcn = malloc(pkd->nMaxCellNewt*sizeof(ILCN));
	assert(pkd->ilcn != NULL);
	pkd->sqrttmp = malloc(pkd->nSqrtTmp*sizeof(double));
	assert(pkd->sqrttmp != NULL);
	pkd->d2a = malloc(pkd->nSqrtTmp*sizeof(double));
	assert(pkd->d2a != NULL);
	/*
	 ** Ewald stuff!
	 */
	pkd->nMaxEwhLoop = 100;
	pkd->ewt = malloc(pkd->nMaxEwhLoop*sizeof(EWT));
	assert(pkd->ewt != NULL);
	*ppkd = pkd;
	}


void pkdFinish(PKD pkd)
{
	if (pkd->kdNodes) free(pkd->kdNodes);
	if (pkd->kdTop) free(pkd->kdTop);
	if (pkd->piLeaf) free(pkd->piLeaf);
	free(pkd->ilp);
	free(pkd->ilcs);
	free(pkd->ilcn);
	free(pkd->sqrttmp);
	free(pkd->d2a);
	free(pkd->ewt);
	free(pkd->pStore);
	free(pkd);
	}


void pkdReadTipsy(PKD pkd,char *pszFileName,int nStart,int nLocal)
{
	FILE *fp;
	int i,j;
	struct dark_particle dp;
	long lStart;

	pkd->nLocal = nLocal;
	/*
	 ** Seek past the header and up to nStart.
	 */
	fp = fopen(pszFileName,"r");
	assert(fp != NULL);
	lStart = sizeof(struct dump)+nStart*sizeof(struct dark_particle);
	fseek(fp,lStart,0);
	/*
	 ** Read Stuff!
	 */
	for (i=0;i<nLocal;++i) {
		fread(&dp,sizeof(struct dark_particle),1,fp);
		for (j=0;j<3;++j) {
			pkd->pStore[i].r[j] = dp.pos[j];
			pkd->pStore[i].v[j] = dp.vel[j];
			}
		pkd->pStore[i].fMass = dp.mass;
		pkd->pStore[i].fSoft = dp.eps;
		pkd->pStore[i].iOrder = nStart + i;
		pkd->pStore[i].fWeight = 1.0;
		}
	fclose(fp);
	}


void pkdCalcBound(PKD pkd,BND *pbnd)
{
	int i,j;

	/*
	 ** Calculate Local Bounds.
	 */
	for (j=0;j<3;++j) {
		pbnd->fMin[j] = pkd->pStore[0].r[j];
		pbnd->fMax[j] = pkd->pStore[0].r[j];
		}
	for (i=1;i<pkd->nLocal;++i) {
		for (j=0;j<3;++j) {
			if (pkd->pStore[i].r[j] < pbnd->fMin[j]) 
				pbnd->fMin[j] = pkd->pStore[i].r[j];
			else if (pkd->pStore[i].r[j] > pbnd->fMax[j])
				pbnd->fMax[j] = pkd->pStore[i].r[j];
			}
		}
	}


int pkdWeight(PKD pkd,int d,float fSplit,int iSplitSide,int iFrom,int iTo,
			  int *pnLow,int *pnHigh,float *pfLow,float *pfHigh)
{
	int i,iPart;
	float fLower,fUpper;
	
	/*
	 ** First partition the memory about fSplit for particles iFrom to iTo.
	 */
	if (iSplitSide) {
		iPart = pkdLowerPart(pkd,d,fSplit,iFrom,iTo);
		*pnLow = pkd->nLocal-iPart;
		*pnHigh = iPart;
		}
	else {
		iPart = pkdUpperPart(pkd,d,fSplit,iFrom,iTo);
		*pnLow = iPart;
		*pnHigh = pkd->nLocal-iPart;
		}
	/*
	 ** Calculate the lower weight and upper weight BETWEEN the particles
	 ** iFrom to iTo!
	 */
	fLower = 0.0;
	for (i=iFrom;i<iPart;++i) {
		fLower += pkd->pStore[i].fWeight;
		}
	fUpper = 0.0;
	for (i=iPart;i<=iTo;++i) {
		fUpper += pkd->pStore[i].fWeight;
		}
	if (iSplitSide) {
		*pfLow = fUpper;
		*pfHigh = fLower;
		}
	else {
		*pfLow = fLower;
		*pfHigh = fUpper;
		}
	return(iPart);
	}


int pkdLowerPart(PKD pkd,int d,float fSplit,int i,int j)
{
	PARTICLE pTemp;

	if (i > j) goto done;
    while (1) {
        while (pkd->pStore[i].r[d] >= fSplit)
            if (++i > j) goto done;
        while (pkd->pStore[j].r[d] < fSplit)
            if (i > --j) goto done;
		pTemp = pkd->pStore[i];
		pkd->pStore[i] = pkd->pStore[j];
		pkd->pStore[j] = pTemp;
        }
 done:
    return(i);
	}


int pkdUpperPart(PKD pkd,int d,float fSplit,int i,int j)
{
	PARTICLE pTemp;

	if (i > j) goto done;
    while (1) {
        while (pkd->pStore[i].r[d] < fSplit)
            if (++i > j) goto done;
        while (pkd->pStore[j].r[d] >= fSplit)
            if (i > --j) goto done;
		pTemp = pkd->pStore[i];
		pkd->pStore[i] = pkd->pStore[j];
		pkd->pStore[j] = pTemp;
        }
 done:
    return(i);
	}


int pkdLowerOrdPart(PKD pkd,int nOrdSplit,int i,int j)
{
	PARTICLE pTemp;

	if (i > j) goto done;
    while (1) {
        while (pkd->pStore[i].iOrder >= nOrdSplit)
            if (++i > j) goto done;
        while (pkd->pStore[j].iOrder < nOrdSplit)
            if (i > --j) goto done;
		pTemp = pkd->pStore[i];
		pkd->pStore[i] = pkd->pStore[j];
		pkd->pStore[j] = pTemp;
        }
 done:
    return(i);
	}


int pkdUpperOrdPart(PKD pkd,int nOrdSplit,int i,int j)
{
	PARTICLE pTemp;

	if (i > j) goto done;
    while (1) {
        while (pkd->pStore[i].iOrder < nOrdSplit)
            if (++i > j) goto done;
        while (pkd->pStore[j].iOrder >= nOrdSplit)
            if (i > --j) goto done;
		pTemp = pkd->pStore[i];
		pkd->pStore[i] = pkd->pStore[j];
		pkd->pStore[j] = pTemp;
        }
 done:
    return(i);
	}


int pkdColRejects(PKD pkd,int d,float fSplit,int iSplitSide)
{
	int nSplit,iRejects,i;

	if (iSplitSide) nSplit = pkdLowerPart(pkd,d,fSplit,0,pkd->nLocal-1);
	else nSplit = pkdUpperPart(pkd,d,fSplit,0,pkd->nLocal-1);
	pkd->nRejects = pkd->nLocal - nSplit;
	iRejects = pkd->nStore - pkd->nRejects;
	pkd->nLocal = nSplit;
	/*
	 ** Move rejects to High memory.
	 */
	for (i=pkd->nRejects-1;i>=0;--i)
		pkd->pStore[iRejects+i] = pkd->pStore[pkd->nLocal+i];
	return(pkd->nRejects);
	}


int pkdSwapRejects(PKD pkd,int idSwap)
{
	int nBuf;
	int nOutBytes,nSndBytes,nRcvBytes;

	if (idSwap != -1) {
		nBuf = (pkd->nStore - pkd->nLocal)*sizeof(PARTICLE);
		nOutBytes = pkd->nRejects*sizeof(PARTICLE);
		mdlSwap(pkd->mdl,idSwap,nBuf,&pkd->pStore[pkd->nLocal],
				nOutBytes,&nSndBytes,&nRcvBytes);
		pkd->nLocal += nRcvBytes/sizeof(PARTICLE);
		pkd->nRejects -= nSndBytes/sizeof(PARTICLE);
		}
	return(pkd->nRejects);
	}


int pkdSwapSpace(PKD pkd)
{
	return(pkdFreeStore(pkd) - pkdLocal(pkd));
	}


int pkdFreeStore(PKD pkd)
{
	return(pkd->nStore);
	}


int pkdLocal(PKD pkd)
{
	return(pkd->nLocal);
	}

int pkdNodes(PKD pkd)
{
	return(pkd->nNodes);
	}


void pkdDomainColor(PKD pkd)
{
	assert(0);
	/*
	int i;
	
	for (i=0;i<pkd->nLocal;++i) {
		pkd->pStore[i].fColor = (float)pkd->idSelf;
		}
		*/
	}


int pkdColOrdRejects(PKD pkd,int nOrdSplit,int iSplitSide)
{
	int nSplit,iRejects,i;

	if (iSplitSide) nSplit = pkdLowerOrdPart(pkd,nOrdSplit,0,pkd->nLocal-1);
	else nSplit = pkdUpperOrdPart(pkd,nOrdSplit,0,pkd->nLocal-1);
	pkd->nRejects = pkd->nLocal - nSplit;
	iRejects = pkd->nStore - pkd->nRejects;
	pkd->nLocal = nSplit;
	/*
	 ** Move rejects to High memory.
	 */
	for (i=pkd->nRejects-1;i>=0;--i)
		pkd->pStore[iRejects+i] = pkd->pStore[pkd->nLocal+i];
	return(pkd->nRejects);
	}


int cmpParticles(const void *pva,const void *pvb)
{
	PARTICLE *pa = (PARTICLE *)pva;
	PARTICLE *pb = (PARTICLE *)pvb;

	return(pa->iOrder - pb->iOrder);
	}


void _pkdOrder(PKD pkd,PARTICLE *pTemp,int nStart,int i,int j)
{
	int nSplit,k,iOrd;

	if (j-i+1 > PKD_ORDERTEMP) {
		nSplit = (i+j+1)/2 + nStart;
		nSplit = pkdUpperOrdPart(pkd,nSplit,i,j);
		_pkdOrder(pkd,pTemp,nStart,i,nSplit-1);
		_pkdOrder(pkd,pTemp,nStart,nSplit,j);
		}
	else {
		/*
		 ** Reorder using the temporary storage!
		 */
		for (k=i;k<=j;++k) {
			iOrd = pkd->pStore[k].iOrder - nStart - i;
			pTemp[iOrd] = pkd->pStore[k];
			}
		for (k=i;k<=j;++k) {
			pkd->pStore[k] = pTemp[k-i];
			}
		}
	}


void pkdLocalOrder(PKD pkd,int nStart)
{
	/*
	 ** My ordering function seems to break on some machines, but
	 ** I still don't understand why. Use qsort instead.
	 **	PARTICLE pTemp[PKD_ORDERTEMP];
	 ** _pkdOrder(pkd,pTemp,nStart,0,pkd->nLocal-1);
	 */
	qsort(pkd->pStore,pkdLocal(pkd),sizeof(PARTICLE),cmpParticles);
	}


void pkdWriteTipsy(PKD pkd,char *pszFileName,int nStart,int nEnd,
				   int bStandard)
{
	FILE *fp;
	int i,j;
	struct dark_particle dp;
	long lStart;

	/*
	 ** First verify order of the particles!
	 */
	for (i=0;i<pkd->nLocal;++i) {
		if (pkd->pStore[i].iOrder != i+nStart) {
			mdlDiag(pkd->mdl,"Order error\n");
			break;
			}
		}
	if (nEnd - nStart + 1 != pkd->nLocal) mdlDiag(pkd->mdl,"Number error\n");
	/*
	 ** Seek past the header and up to nStart.
	 */
	fp = fopen(pszFileName,"r+");
	assert(fp != NULL);
	if (bStandard) {
#ifndef CRAY_T3D
		XDR xdrs;
		/*
		 ** Seek according to true XDR size structures!
		 ** This may be a bit dicey, but it should work as long
		 ** as no one changes the tipsy binary format!
		 */
		assert(sizeof(struct dark_particle) == 9*sizeof(float));
		assert(sizeof(struct dump) == sizeof(double)+6*sizeof(int));
		lStart = 32 + nStart*36;
		fseek(fp,lStart,0);
		/* 
		 ** Write Stuff!
		 */
		xdrstdio_create(&xdrs,fp,XDR_ENCODE);
		for (i=0;i<pkd->nLocal;++i) {
			xdr_float(&xdrs,&pkd->pStore[i].fMass);
			for (j=0;j<3;++j) {
				xdr_float(&xdrs,&pkd->pStore[i].r[j]);
				}
			for (j=0;j<3;++j) {
				xdr_float(&xdrs,&pkd->pStore[i].v[j]);
				}
			xdr_float(&xdrs,&pkd->pStore[i].fSoft);
			xdr_float(&xdrs,&pkd->pStore[i].fPot);
			}
		xdr_destroy(&xdrs);
#endif
		}
	else {
		lStart = sizeof(struct dump)+nStart*sizeof(struct dark_particle);
		fseek(fp,lStart,0);
		/* 
		 ** Write Stuff!
		 */
		for (i=0;i<pkd->nLocal;++i) {
			for (j=0;j<3;++j) {
				dp.pos[j] = pkd->pStore[i].r[j];
				dp.vel[j] = pkd->pStore[i].v[j];
				}
			dp.mass = pkd->pStore[i].fMass;
			dp.eps = pkd->pStore[i].fSoft;
			dp.phi = pkd->pStore[i].fPot;
			fwrite(&dp,sizeof(struct dark_particle),1,fp);
			}
		}
	fclose(fp);
	}


void pkdSetSoft(PKD pkd,double dSoft)
{
	PARTICLE *p;
	int i,n;

	p = pkd->pStore;
	n = pkdLocal(pkd);
	for (i=0;i<n;++i) {
             p[i].fSoft = dSoft;
             }
        }


void pkdCombine(KDN *p1,KDN *p2,KDN *pOut)
{
	KDN *p;
	int j;
	float dx,dy,dz;	

	/*
	 ** Combine the bounds.
	 */
	for (j=0;j<3;++j) {
		if (p2->bnd.fMin[j] < p1->bnd.fMin[j])
			pOut->bnd.fMin[j] = p2->bnd.fMin[j];
		else
			pOut->bnd.fMin[j] = p1->bnd.fMin[j];
		if (p2->bnd.fMax[j] > p1->bnd.fMax[j])
			pOut->bnd.fMax[j] = p2->bnd.fMax[j];
		else
			pOut->bnd.fMax[j] = p1->bnd.fMax[j];
		}
	/*
	 ** Find the center of mass and mass weighted softening.
	 */
    pOut->fMass = p1->fMass + p2->fMass;
	pOut->fSoft = (p1->fMass*p1->fSoft + p2->fMass*p2->fSoft)/pOut->fMass;
	for (j=0;j<3;++j) {
		pOut->r[j] = (p1->fMass*p1->r[j] + p2->fMass*p2->r[j])/pOut->fMass;
		}
	}


void pkdCalcCell(PKD pkd,KDN *pkdn,float *rcm,int iOrder,
				 struct pkdCalcCellStruct *pcc)
{
	int pj;
	double m,dx,dy,dz,d2,d1;

	/*
	 ** Initialize moments.
	 ** Initialize various B numbers.
	 */
	switch (iOrder) {	
	case 4:
		pcc->Hxxxx = 0.0;
		pcc->Hxyyy = 0.0;
		pcc->Hxxxy = 0.0;
		pcc->Hyyyy = 0.0;
		pcc->Hxxxz = 0.0;
		pcc->Hyyyz = 0.0;
		pcc->Hxxyy = 0.0;
		pcc->Hxxyz = 0.0;
		pcc->Hxyyz = 0.0;
		pcc->B6 = 0.0;
	case 3:
		pcc->Oxxx = 0.0;
		pcc->Oxyy = 0.0;
		pcc->Oxxy = 0.0;
		pcc->Oyyy = 0.0;
		pcc->Oxxz = 0.0;
		pcc->Oyyz = 0.0;
		pcc->Oxyz = 0.0;
		pcc->B5 = 0.0;
	default:
		pcc->Qxx = 0.0;
		pcc->Qyy = 0.0;
		pcc->Qzz = 0.0;
		pcc->Qxy = 0.0;
		pcc->Qxz = 0.0;
		pcc->Qyz = 0.0;
		pcc->B2 = 0.0;
		pcc->B3 = 0.0;
		pcc->B4 = 0.0;
		}
	pcc->Bmax = 0.0;
	/*
	 ** Calculate moments and B numbers about center-of-mass.
	 */
	if (pkdn == NULL) pkdn = &pkd->kdNodes[ROOT];
	for (pj=pkdn->pLower;pj<=pkdn->pUpper;++pj) {
		m = pkd->pStore[pj].fMass;
		dx = pkd->pStore[pj].r[0] - rcm[0];
		dy = pkd->pStore[pj].r[1] - rcm[1];
		dz = pkd->pStore[pj].r[2] - rcm[2];
		d2 = dx*dx + dy*dy + dz*dz;
		d1 = sqrt(d2);
		if (d1 > pcc->Bmax) pcc->Bmax = d1;
		switch (iOrder) {
		case 4:
			/*
			 ** Calculate reduced hexadecapole moment...
			 */
			pcc->Hxxxx += m*(dx*dx*dx*dx - 6.0/7.0*d2*(dx*dx - 0.1*d2));
			pcc->Hxyyy += m*(dx*dy*dy*dy - 3.0/7.0*d2*dx*dy);
			pcc->Hxxxy += m*(dx*dx*dx*dy - 3.0/7.0*d2*dx*dy);
			pcc->Hyyyy += m*(dy*dy*dy*dy - 6.0/7.0*d2*(dy*dy - 0.1*d2));
			pcc->Hxxxz += m*(dx*dx*dx*dz - 3.0/7.0*d2*dx*dz);
			pcc->Hyyyz += m*(dy*dy*dy*dz - 3.0/7.0*d2*dy*dz);
			pcc->Hxxyy += m*(dx*dx*dy*dy - 1.0/7.0*d2*(dx*dx + dy*dy - 0.2*d2));
			pcc->Hxxyz += m*(dx*dx*dy*dz - 1.0/7.0*d2*dy*dz);
			pcc->Hxyyz += m*(dx*dy*dy*dz - 1.0/7.0*d2*dx*dz);
			pcc->B6 += m*d2*d2*d2;
		case 3:
			/*
			 ** Calculate reduced octopole moment...
			 */
			pcc->Oxxx += m*(dx*dx*dx - 0.6*d2*dx);
			pcc->Oxyy += m*(dx*dy*dy - 0.2*d2*dx);
			pcc->Oxxy += m*(dx*dx*dy - 0.2*d2*dy);
			pcc->Oyyy += m*(dy*dy*dy - 0.6*d2*dy);
			pcc->Oxxz += m*(dx*dx*dz - 0.2*d2*dz);
			pcc->Oyyz += m*(dy*dy*dz - 0.2*d2*dz);
			pcc->Oxyz += m*dx*dy*dz;
			pcc->B5 += m*d2*d2*d1;
		default:
			/*
			 ** Calculate quadrupole moment...
			 */
			pcc->Qxx += m*dx*dx;
			pcc->Qyy += m*dy*dy;
			pcc->Qzz += m*dz*dz;
			pcc->Qxy += m*dx*dy;
			pcc->Qxz += m*dx*dz;
			pcc->Qyz += m*dy*dz;
			pcc->B2 += m*d2;
			pcc->B3 += m*d2*d1;
			pcc->B4 += m*d2*d2;
			}
		}
	}


double fcnAbsMono(KDN *pkdn,double r)
{
	double t;

	t = r*(r - pkdn->mom.Bmax);
	t *= t;
	t = (3.0*pkdn->mom.B2 - 2.0*pkdn->mom.B3/r)/t;
	return t;
	}


double fcnRelMono(KDN *pkdn,double r)
{
	double t;

	t = (r - pkdn->mom.Bmax);
	t *= pkdn->fMass*t;
	t = (3.0*pkdn->mom.B2 - 2.0*pkdn->mom.B3/r)/t;
	return t;
	}


double fcnAbsQuad(KDN *pkdn,double r)
{
	double t;

	t = r*(r - pkdn->mom.Bmax);
	t *= r*t;
	t = (4.0*pkdn->mom.B3 - 3.0*pkdn->mom.B4/r)/t;
	return t;
	}


double fcnRelQuad(KDN *pkdn,double r)
{
	double t;

	t = (r - pkdn->mom.Bmax);
	t *= pkdn->fMass*r*t;
	t = (4.0*pkdn->mom.B3 - 3.0*pkdn->mom.B4/r)/t;
	return t;
	}


double fcnAbsOct(KDN *pkdn,double r)
{
	double t;

	t = r*r*(r - pkdn->mom.Bmax);
	t *= t;
	t = (5.0*pkdn->mom.B4 - 4.0*pkdn->mom.B5/r)/t;
	return t;
	}


double fcnRelOct(KDN *pkdn,double r)
{
	double t;

	t = r*(r - pkdn->mom.Bmax);
	t *= pkdn->fMass*t;
	t = (5.0*pkdn->mom.B4 - 4.0*pkdn->mom.B5/r)/t;
	return t;
	}


double fcnAbsHex(KDN *pkdn,double r)
{
	double t;

	t = r*r*(r - pkdn->mom.Bmax);
	t *= r*t;
	t = (6.0*pkdn->mom.B5 - 5.0*pkdn->mom.B6/r)/t;
	return t;
	}


double fcnRelHex(KDN *pkdn,double r)
{
	double t;

	t = r*(r - pkdn->mom.Bmax);
	t *= pkdn->fMass*r*t;
	t = (6.0*pkdn->mom.B5 - 5.0*pkdn->mom.B6/r)/t;
	return t;
	}


double dRootBracket(KDN *pkdn,double dErrBnd,double (*fcn)(KDN *,double))
{
	double dLower,dUpper,dMid;
	double dLowerErr,dErrDif,dErrCrit;
	int iter = 0;

	dErrCrit = 1e-6*dErrBnd;
	dLower = (1.0 + 1e-6)*pkdn->mom.Bmax;
	dLowerErr = (*fcn)(pkdn,dLower);
	if (dLowerErr < dErrBnd) {
/*
		printf("iter:%d %g\n",iter,dLower);
*/
		return dLower;
		}
	/*
	 ** Hunt phase for upper bound.
	 */
	dUpper = 2*dLower;
	while (1) {
		dErrDif = dErrBnd - (*fcn)(pkdn,dUpper);
		if (dErrDif > dErrCrit) break;
		dUpper = 2*dUpper;
		}
	/*
	 ** Bracket root.
	 */
	while (1) {
		dMid = 0.5*(dLower + dUpper);
		dErrDif = dErrBnd - (*fcn)(pkdn,dMid);
		if (dErrDif < 0) {
			dLower = dMid;
			}
		else {
			dUpper = dMid;
			if (dErrDif < dErrCrit) break;
			}
		if (++iter > 32) break;
		}
/*
	printf("iter:%d %g\n",iter,dMid);
*/
	return dMid;
	}


double pkdCalcOpen(KDN *pkdn,int iOpenType,double dCrit,int iOrder)
{
	double dOpen;

	if (iOpenType == OPEN_ABSPAR) {
		switch (iOrder) {
		case 1:
			dOpen = dRootBracket(pkdn,dCrit,fcnAbsMono);
			break;
		case 2:
			dOpen = dRootBracket(pkdn,dCrit,fcnAbsQuad);
			break;
		case 3:
			dOpen = dRootBracket(pkdn,dCrit,fcnAbsOct);
			break;
		case 4:
			dOpen = dRootBracket(pkdn,dCrit,fcnAbsHex);
			break;
			}
		}
	else if (iOpenType == OPEN_RELPAR) {
		switch (iOrder) {
		case 1:
			dOpen = dRootBracket(pkdn,dCrit,fcnRelMono);
			break;
		case 2:
			dOpen = dRootBracket(pkdn,dCrit,fcnRelQuad);
			break;
		case 3:
			dOpen = dRootBracket(pkdn,dCrit,fcnRelOct);
			break;
		case 4:
			dOpen = dRootBracket(pkdn,dCrit,fcnRelHex);
			break;
			}
		}
	else if (iOpenType == OPEN_JOSH) {
		/*
		 ** Set openening criterion to an approximation of Josh's theta.
		 ** Will be equivalent to Josh's theta for sufficiently cubical
		 ** cells. The longest side will be chosen.
		 */
		int j;
		double d,Bdel,dx,dy,dz,d2;
		float rc[3];

		Bdel = 0.0;
		for (j=0;j<3;++j) {
			rc[j] = 0.5*(pkdn->bnd.fMin[j] + pkdn->bnd.fMax[j]);
			d = pkdn->bnd.fMax[j] - pkdn->bnd.fMin[j];
			if (d > Bdel) Bdel = d;
			}
		dx = pkdn->r[0] - rc[0];
		dy = pkdn->r[1] - rc[1];
		dz = pkdn->r[2] - rc[2];
		d2 = dx*dx + dy*dy + dz*dz;
		dOpen = 2/sqrt(3.0)*pkdn->mom.Bmax/dCrit;
		if (dOpen < pkdn->mom.Bmax) dOpen = pkdn->mom.Bmax;
		}
	else {
		/*
		 ** Set opening criterion to the minimal, i.e., Bmax.
		 */
		dOpen = pkdn->mom.Bmax;
		}
	return(dOpen);
	}


void pkdUpPass(PKD pkd,int iCell,int iOpenType,double dCrit,int iOrder)
{
	KDN *c;
	PARTICLE *p;
	int l,u,pj,j;
	double dx,dy,dz;
	double dOpen;
	float rx,ry,rz,d2,d1;

	c = pkd->kdNodes;
	p = pkd->pStore;
	l = c[iCell].pLower;
	u = c[iCell].pUpper;
	if (c[iCell].iDim >= 0) {
		pkdUpPass(pkd,LOWER(iCell),iOpenType,dCrit,iOrder);
		pkdUpPass(pkd,UPPER(iCell),iOpenType,dCrit,iOrder);
		pkdCombine(&c[LOWER(iCell)],&c[UPPER(iCell)],&c[iCell]);
		}
	else {
		c[iCell].fMass = 0.0;
		c[iCell].fSoft = 0.0;
		for (j=0;j<3;++j) {
		    	c[iCell].bnd.fMin[j] = p[u].r[j];
			c[iCell].bnd.fMax[j] = p[u].r[j];
			c[iCell].r[j] = 0.0;
			}
		for (pj=l;pj<=u;++pj) {
			for (j=0;j<3;++j) {
				if (p[pj].r[j] < c[iCell].bnd.fMin[j])
					c[iCell].bnd.fMin[j] = p[pj].r[j];
				if (p[pj].r[j] > c[iCell].bnd.fMax[j])
					c[iCell].bnd.fMax[j] = p[pj].r[j];
				}
			/*
			 ** Find center of mass and total mass and mass weighted softening.
			 */
			c[iCell].fMass += p[pj].fMass;
			c[iCell].fSoft += p[pj].fMass*p[pj].fSoft;
			for (j=0;j<3;++j) {
				c[iCell].r[j] += p[pj].fMass*p[pj].r[j];
				}
			}
		for (j=0;j<3;++j) c[iCell].r[j] /= c[iCell].fMass;
		c[iCell].fSoft /= c[iCell].fMass;
		}
	/*
	 ** Calculate multipole moments.
	 */
	pkdCalcCell(pkd,&c[iCell],c[iCell].r,iOrder,&c[iCell].mom);
	dOpen = pkdCalcOpen(&c[iCell],iOpenType,dCrit,iOrder);
	c[iCell].fOpen2 = dOpen*dOpen;
	}


/*
 ** JST's Select
 */
void pkdSelect(PKD pkd,int d,int k,int l,int r)
{
	PARTICLE *p,t;
	float v;
	int i,j;

	p = pkd->pStore;
	while (r > l) {
		v = p[k].r[d];
		t = p[r];
		p[r] = p[k];
		p[k] = t;
		i = l - 1;
		j = r;
		while (1) {
			while (i < j) if (p[++i].r[d] >= v) break;
			while (i < j) if (p[--j].r[d] <= v) break;
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


void pkdBuildLocal(PKD pkd,int nBucket,int iOpenType,double dCrit,
				   int iOrder,KDN *pRoot)
{
	int l,n,i,d,m,j,diff;
	KDN *c;
	char ach[256];

	pkd->nBucket = nBucket;
	n = pkd->nLocal;
	pkd->nLevels = 1;
	l = 1;
	while (n > nBucket) {
		n = n>>1;
		l = l<<1;
		++pkd->nLevels;
		}
	pkd->nSplit = l;
	pkd->nNodes = l<<1;
	if (pkd->kdNodes) free(pkd->kdNodes);
	pkd->kdNodes = (KDN *)malloc(pkd->nNodes*sizeof(KDN));
	assert(pkd->kdNodes != NULL);
	sprintf(ach,"nNodes:%d nSplit:%d nLevels:%d nBucket:%d\n",
			pkd->nNodes,pkd->nSplit,pkd->nLevels,nBucket);
	mdlDiag(pkd->mdl,ach);
	/*
	 ** Set up ROOT node
	 */
	c = pkd->kdNodes;
	c[ROOT].pLower = 0;
	c[ROOT].pUpper = pkd->nLocal-1;
	/*
	 ** determine the local bound of the particles.
	 */
	pkdCalcBound(pkd,&c[ROOT].bnd);
	i = ROOT;
	while (1) {
		assert(c[i].pUpper - c[i].pLower + 1 > 0);
		if (i < pkd->nSplit && (c[i].pUpper - c[i].pLower) > 0) {
			d = 0;
			for (j=1;j<3;++j) {
				if (c[i].bnd.fMax[j]-c[i].bnd.fMin[j] > 
					c[i].bnd.fMax[d]-c[i].bnd.fMin[d]) d = j;
				}
			c[i].iDim = d;

			m = (c[i].pLower + c[i].pUpper)/2;
			pkdSelect(pkd,d,m,c[i].pLower,c[i].pUpper);

			c[i].fSplit = pkd->pStore[m].r[d];
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
			SETNEXT(i);
			if (i == ROOT) break;
			}
		}
	pkdUpPass(pkd,ROOT,iOpenType,dCrit,iOrder);
	*pRoot = c[ROOT];
	}


void pkdBucketWeight(PKD pkd,int iBucket,float fWeight)
{
	KDN *pbuc;
	int pj;
	
	pbuc = &pkd->kdNodes[iBucket];
	for (pj=pbuc->pLower;pj<=pbuc->pUpper;++pj) {
		pkd->pStore[pj].fWeight = fWeight;
		}
	}


void pkdGravAll(PKD pkd,int nReps,int bPeriodic,int iOrder,int iEwOrder,
				double fEwCut,double fEwhCut,
				double *pdPartSum,double *pdCellSum,CASTAT *pcsPart,
				CASTAT *pcsCell)
{
	KDN *c;
	int iCell,n;
	float fWeight;

	pkdClearTimer(pkd,1);
	pkdClearTimer(pkd,2);
	pkdClearTimer(pkd,3);
	*pdPartSum = 0.0;
	*pdCellSum = 0.0;
	if (bPeriodic) {
		pkdStartTimer(pkd,3);
		pkdEwaldInit(pkd,fEwhCut,iEwOrder);
		pkdStopTimer(pkd,3);
		}
	/*
	 ** Start caching spaces.
	 */
	mdlROcache(pkd->mdl,CID_PARTICLE,pkd->pStore,sizeof(PARTICLE),
			   pkdLocal(pkd));
	mdlROcache(pkd->mdl,CID_CELL,pkd->kdNodes,sizeof(KDN),pkdNodes(pkd));
	/*
	 ** Walk over the local buckets!
	 */
	c = pkd->kdNodes;
	iCell = ROOT;
	while (1) {
		if (c[iCell].iDim >= 0) iCell = LOWER(iCell);
		else {
			n = c[iCell].pUpper - c[iCell].pLower + 1;
			pkdStartTimer(pkd,1);
			pkdBucketWalk(pkd,iCell,nReps,iOrder);
			pkdStopTimer(pkd,1);
			*pdPartSum += n*pkd->nPart + n*(n-1)/2;
			*pdCellSum += n*(pkd->nCellSoft + pkd->nCellNewt);
			pkdStartTimer(pkd,2);
			pkdBucketInteract(pkd,iCell,iOrder);
			pkdStopTimer(pkd,2);
			if (bPeriodic) {
				pkdStartTimer(pkd,3);
				pkdBucketEwald(pkd,iCell,nReps,fEwCut,iEwOrder);
				pkdStopTimer(pkd,3);
				}
			fWeight = 2.0*(pkd->nCellSoft + pkd->nCellNewt) + 
			    1.0*(pkd->nPart + (n-1)/2.0);
			pkdBucketWeight(pkd,iCell,fWeight);
			SETNEXT(iCell);
			if (iCell == ROOT) break;
			}
		}
	/*
	 ** Get caching statistics.
	 */
	mdlCacheStat(pkd->mdl,CID_CELL,pcsCell);
	mdlCacheStat(pkd->mdl,CID_PARTICLE,pcsPart);
	/*
	 ** Stop caching spaces.
	 */
	mdlFinishCache(pkd->mdl,CID_CELL);
	mdlFinishCache(pkd->mdl,CID_PARTICLE);
	}


void pkdCalcE(PKD pkd,double *T,double *U)
{
	PARTICLE *p;
	int i,n;
	float vx,vy,vz;

	p = pkd->pStore;
	n = pkdLocal(pkd);
	*T = 0.0;
	*U = 0.0;
	for (i=0;i<n;++i) {
		*U += 0.5*p[i].fMass*p[i].fPot;
		vx = p[i].v[0];
		vy = p[i].v[1];
		vz = p[i].v[2];
		*T += 0.5*p[i].fMass*(vx*vx + vy*vy + vz*vz);
		}
	}


void pkdDrift(PKD pkd,double dDelta,float fCenter[3])
{
	PARTICLE *p;
	int i,j,n;

	p = pkd->pStore;
	n = pkdLocal(pkd);
	for (i=0;i<n;++i) {
		for (j=0;j<3;++j) {
			p[i].r[j] += dDelta*p[i].v[j];
			if (p[i].r[j] >= fCenter[j]+0.5*pkd->fPeriod[j])
				p[i].r[j] -= pkd->fPeriod[j];
			if (p[i].r[j] < fCenter[j]-0.5*pkd->fPeriod[j])
				p[i].r[j] += pkd->fPeriod[j];
			}
		}
	}


void pkdKick(PKD pkd,double dvFacOne,double dvFacTwo)
{
	PARTICLE *p;
	int i,j,n;

	p = pkd->pStore;
	n = pkdLocal(pkd);
	for (i=0;i<n;++i) {
		for (j=0;j<3;++j) {
			p[i].v[j] = p[i].v[j]*dvFacOne + p[i].a[j]*dvFacTwo;
			}
		}
	}


void pkdReadCheckNew(PKD pkd,char *pszFileName,int nStart,int nLocal)
{
	FILE *fp;
	CHKPART cp;
	long lStart;
	int i,j;

	pkd->nLocal = nLocal;
	/*
	 ** Seek past the header and up to nStart.
	 */
	fp = fopen(pszFileName,"r");
	assert(fp != NULL);
	lStart = sizeof(struct msrCheckPointHeader)+nStart*sizeof(CHKPART);
	fseek(fp,lStart,0);
	/*
	 ** Read Stuff!
	 */
	for (i=0;i<nLocal;++i) {
		fread(&cp,sizeof(CHKPART),1,fp);
		pkd->pStore[i].iOrder = cp.iOrder;
		pkd->pStore[i].fMass = cp.fMass;
		pkd->pStore[i].fSoft = cp.fSoft;
		for (j=0;j<3;++j) {
			pkd->pStore[i].r[j] = cp.r[j];
			pkd->pStore[i].v[j] = cp.v[j];
			}
		pkd->pStore[i].fWeight = 1.0;	/* set the initial weight to 1.0 */
		}
	fclose(fp);	
	}


void pkdWriteCheckNew(PKD pkd,char *pszFileName,int nStart)
{
	FILE *fp;
	CHKPART cp;
	long lStart;
	int i,j,nLocal;

	/*
	 ** Seek past the header and up to nStart.
	 */
	fp = fopen(pszFileName,"r+");
	assert(fp != NULL);
	lStart = sizeof(struct msrCheckPointHeader)+nStart*sizeof(CHKPART);
	fseek(fp,lStart,0);
	/* 
	 ** Write Stuff!
	 */
	nLocal = pkdLocal(pkd);
	for (i=0;i<nLocal;++i) {
		cp.iOrder = pkd->pStore[i].iOrder;
		cp.fMass = pkd->pStore[i].fMass;
		cp.fSoft = pkd->pStore[i].fSoft;
		for (j=0;j<3;++j) {
			cp.r[j] = pkd->pStore[i].r[j];
			cp.v[j] = pkd->pStore[i].v[j];
			}
		fwrite(&cp,sizeof(CHKPART),1,fp);
		}
	fclose(fp);
	}


void pkdReadCheckOld(PKD pkd,char *pszFileName,int nStart,int nLocal)
{
	FILE *fp;
	OCHKPART cp;
	long lStart;
	int i,j;

	pkd->nLocal = nLocal;
	/*
	 ** Seek past the header and up to nStart.
	 */
	fp = fopen(pszFileName,"r");
	assert(fp != NULL);
	lStart = sizeof(struct msrCheckPointHeader)+nStart*sizeof(OCHKPART);
	fseek(fp,lStart,0);
	/*
	 ** Read Stuff!
	 */
	for (i=0;i<nLocal;++i) {
		fread(&cp,sizeof(OCHKPART),1,fp);
		pkd->pStore[i].iOrder = cp.iOrder;
		pkd->pStore[i].fMass = cp.fMass;
		pkd->pStore[i].fSoft = cp.fSoft;
		for (j=0;j<3;++j) {
			pkd->pStore[i].r[j] = cp.r[j];
			pkd->pStore[i].v[j] = cp.v[j];
			}
		pkd->pStore[i].fWeight = 1.0;	/* set the initial weight to 1.0 */
		}
	fclose(fp);	
	}


void pkdDistribCells(PKD pkd,int nCell,KDN *pkdn)
{
	int i;

	if (pkd->kdTop != NULL) free(pkd->kdTop);
	if (pkd->piLeaf != NULL) free(pkd->piLeaf);
	pkd->kdTop = malloc(nCell*sizeof(KDN));
	assert(pkd->kdTop != NULL);
	pkd->piLeaf = malloc(pkd->nThreads*sizeof(int));
	assert(pkd->piLeaf != NULL);
	for (i=1;i<nCell;++i) {
		pkd->kdTop[i] = pkdn[i];
		if (pkdn[i].pLower >= 0) pkd->piLeaf[pkdn[i].pLower] = i;
		}
	}




