#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <malloc.h>
#include <math.h>
#include <assert.h>
#include <sys/time.h>

#include <rpc/types.h>
#include <rpc/xdr.h>

#include "pkd.h"
#include "ewald.h"
#include "grav.h"
#include "walk.h"
#include "opentype.h"
#include "mdl.h"
#include "tipsydefs.h"


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
	pkd->nLocal = 0;
	pkd->nRejects = 0;
	for (j=0;j<3;++j) {
		pkd->fPeriod[j] = fPeriod[j];
		}
	/*
	 ** Allocate the main particle store.
	 ** Need to use mdlMalloc() since the particles will need to be
	 ** visible to all other processors thru mdlAquire() later on.
	 */
	pkd->pStore = mdlMalloc(pkd->mdl,nStore*sizeof(PARTICLE));
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
	if (pkd->kdNodes) mdlFree(pkd->mdl,pkd->kdNodes);
	if (pkd->kdTop) free(pkd->kdTop);
	if (pkd->piLeaf) free(pkd->piLeaf);
	free(pkd->ilp);
	free(pkd->ilcs);
	free(pkd->ilcn);
	free(pkd->sqrttmp);
	free(pkd->d2a);
	free(pkd->ewt);
	mdlFree(pkd->mdl,pkd->pStore);
	free(pkd);
	}


void pkdReadTipsy(PKD pkd,char *pszFileName,int nStart,int nLocal,double dvFac)
{
	FILE *fp;
	int i,j;
	struct dark_particle dp;
	long lStart;

	pkd->nLocal = nLocal;
	pkd->nActive = nLocal;
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
			pkd->pStore[i].v[j] = dvFac*dp.vel[j];
			}
		pkd->pStore[i].fMass = dp.mass;
		pkd->pStore[i].fSoft = dp.eps;
		pkd->pStore[i].iOrder = nStart + i;
		pkd->pStore[i].iActive = 1;
		pkd->pStore[i].iRung = 0;
		pkd->pStore[i].fWeight = 1.0;
		pkd->pStore[i].fDensity = 0.0;
		pkd->pStore[i].fBall2 = 0.0;
		}
	fclose(fp);
	}


void pkdCalcBound(PKD pkd,BND *pbnd,BND *pbndActive)
{
	int i,j,bFirst;

	/*
	 ** Initialize the bounds to 0 at the beginning
	 */
	for (j=0;j<3;++j) {
		pbnd->fMin[j] = 0.0;
		pbnd->fMax[j] = 0.0;
		pbndActive->fMin[j] = 0.0;
		pbndActive->fMax[j] = 0.0;
		}
	/*
	 ** Calculate Local Bounds.
	 */
	bFirst = 1;
	for (i=0;i<pkd->nLocal;++i) {
		if (bFirst) {
			bFirst = 0;
			for (j=0;j<3;++j) {
				pbnd->fMin[j] = pkd->pStore[i].r[j];
				pbnd->fMax[j] = pkd->pStore[i].r[j];
				}
			}
		else {
			for (j=0;j<3;++j) {
				if (pkd->pStore[i].r[j] < pbnd->fMin[j]) 
					pbnd->fMin[j] = pkd->pStore[i].r[j];
				else if (pkd->pStore[i].r[j] > pbnd->fMax[j])
					pbnd->fMax[j] = pkd->pStore[i].r[j];
				}
			}			
		}
	/*
	 ** Calculate Active Bounds.
	 */
	bFirst = 1;
	for (i=0;i<pkd->nLocal;++i) {
		if (pkd->pStore[i].iActive) {
			if (bFirst) {
				bFirst = 0;
				for (j=0;j<3;++j) {
					pbndActive->fMin[j] = pkd->pStore[i].r[j];
					pbndActive->fMax[j] = pkd->pStore[i].r[j];
					}
				}
			else {
				for (j=0;j<3;++j) {
					if (pkd->pStore[i].r[j] < pbndActive->fMin[j]) 
						pbndActive->fMin[j] = pkd->pStore[i].r[j];
					else if (pkd->pStore[i].r[j] > pbndActive->fMax[j])
						pbndActive->fMax[j] = pkd->pStore[i].r[j];
					}
				}			
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
		*pnLow = pkdLocal(pkd)-iPart;
		*pnHigh = iPart;
		}
	else {
		iPart = pkdUpperPart(pkd,d,fSplit,iFrom,iTo);
		*pnLow = iPart;
		*pnHigh = pkdLocal(pkd)-iPart;
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


void pkdActiveOrder(PKD pkd)
{
	PARTICLE pTemp;
	int i=0;
	int j=pkdLocal(pkd)-1;

	if (i > j) goto done;
    while (1) {
        while (pkd->pStore[i].iActive)
            if (++i > j) goto done;
        while (!pkd->pStore[j].iActive)
            if (i > --j) goto done;
		pTemp = pkd->pStore[i];
		pkd->pStore[i] = pkd->pStore[j];
		pkd->pStore[j] = pTemp;
        }
 done:
	pkd->nActive = i;
	}


int pkdColRejects(PKD pkd,int d,float fSplit,float fSplitInactive,
				  int iSplitSide)
{
	PARTICLE pTemp;
	int nSplit,nSplitInactive,iRejects,i,j;

	assert(pkd->nRejects == 0);
	if (iSplitSide) {
		nSplit = pkdLowerPart(pkd,d,fSplit,0,pkdActive(pkd)-1);
		}
	else {
		nSplit = pkdUpperPart(pkd,d,fSplit,0,pkdActive(pkd)-1);
		}
	if (iSplitSide) {
		nSplitInactive = pkdLowerPart(pkd,d,fSplitInactive,
					      pkdActive(pkd),pkdLocal(pkd)-1);
		}
	else {
		nSplitInactive = pkdUpperPart(pkd,d,fSplitInactive,
					      pkdActive(pkd),pkdLocal(pkd)-1);
		}
	for(i = 0; i < nSplit; ++i)
	    assert(pkd->pStore[i].iActive == 1);
	for(i = pkdActive(pkd); i < nSplitInactive; ++i)
	    assert(pkd->pStore[i].iActive == 0);

	nSplitInactive -= pkdActive(pkd);
	/*
	 ** Now do some fancy rearrangement.
	 */
	i = nSplit;
	j = pkdActive(pkd);
	while (j < pkdActive(pkd) + nSplitInactive) {
		pTemp = pkd->pStore[i];
		pkd->pStore[i] = pkd->pStore[j];
		pkd->pStore[j] = pTemp;
		++i; ++j;
		}
	pkd->nRejects = pkdLocal(pkd) - nSplit - nSplitInactive;
	iRejects = pkdFreeStore(pkd) - pkd->nRejects;
	pkd->nActive = nSplit;
	pkd->nLocal = nSplit + nSplitInactive;
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
		nBuf = (pkdSwapSpace(pkd))*sizeof(PARTICLE);
		nOutBytes = pkd->nRejects*sizeof(PARTICLE);
		assert(pkdLocal(pkd) + pkd->nRejects <= pkdFreeStore(pkd));
		mdlSwap(pkd->mdl,idSwap,nBuf,&pkd->pStore[pkdLocal(pkd)],
				nOutBytes,&nSndBytes,&nRcvBytes);
		pkd->nLocal += nRcvBytes/sizeof(PARTICLE);
		pkd->nRejects -= nSndBytes/sizeof(PARTICLE);
		}
	return(pkd->nRejects);
	}

void pkdSwapAll(PKD pkd, int idSwap)
{
    int nBuf;
    int nOutBytes,nSndBytes,nRcvBytes;
    int i;
    int iBuf;
    
    /*
     ** Move particles to High memory.
     */
    iBuf = pkdSwapSpace(pkd);
    for (i=pkdLocal(pkd)-1;i>=0;--i)
	pkd->pStore[iBuf+i] = pkd->pStore[i];

    nBuf = pkdFreeStore(pkd)*sizeof(PARTICLE);
    nOutBytes = pkdLocal(pkd)*sizeof(PARTICLE);
    mdlSwap(pkd->mdl,idSwap,nBuf,&pkd->pStore[0], nOutBytes,
	    &nSndBytes, &nRcvBytes);
    assert(nSndBytes/sizeof(PARTICLE) == pkdLocal(pkd));
    pkd->nLocal = nRcvBytes/sizeof(PARTICLE);
    }

int pkdSwapSpace(PKD pkd)
{
	return(pkdFreeStore(pkd) - pkdLocal(pkd));
	}


int pkdFreeStore(PKD pkd)
{
	return(pkd->nStore);
	}

int pkdActive(PKD pkd)
{
	return(pkd->nActive);
	}

int pkdInactive(PKD pkd)
{
	return(pkd->nLocal - pkd->nActive);
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
#if (0)
	int i;
	
	for (i=0;i<pkd->nLocal;++i) {
		pkd->pStore[i].fColor = (float)pkd->idSelf;
		}
#endif
	}


int pkdColOrdRejects(PKD pkd,int nOrdSplit,int iSplitSide)
{
	int nSplit,iRejects,i;

	if (iSplitSide) nSplit = pkdLowerOrdPart(pkd,nOrdSplit,0,pkdLocal(pkd)-1);
	else nSplit = pkdUpperOrdPart(pkd,nOrdSplit,0,pkdLocal(pkd)-1);
	pkd->nRejects = pkdLocal(pkd) - nSplit;
	iRejects = pkdFreeStore(pkd) - pkd->nRejects;
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
				   int bStandard,double dvFac)
{
	FILE *fp;
	int i,j;
	struct dark_particle dp;
	long lStart;
	int nout;

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
		float vTemp;
		XDR xdrs;
		/*
		 ** Seek according to true XDR size structures!
		 ** This may be a bit dicey, but it should work as long
		 ** as no one changes the tipsy binary format!
		 */
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
				vTemp = dvFac*pkd->pStore[i].v[j];			
				xdr_float(&xdrs,&vTemp);
				}
			xdr_float(&xdrs,&pkd->pStore[i].fSoft);
			xdr_float(&xdrs,&pkd->pStore[i].fPot);
			}
		xdr_destroy(&xdrs);
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
				dp.vel[j] = dvFac*pkd->pStore[i].v[j];
				}
			dp.mass = pkd->pStore[i].fMass;
			dp.eps = pkd->pStore[i].fSoft;
			dp.phi = pkd->pStore[i].fPot;
			nout = fwrite(&dp,sizeof(struct dark_particle),1,fp);
			assert(nout == 1);
			}
		}
	nout = fclose(fp);
	assert(nout == 0);
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
	int j;

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
	