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

#ifdef COLLISIONS
#include "ssdefs.h"
#include "collision.h"
#endif /* COLLISIONS */

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
				   FLOAT *fPeriod,int nDark,int nGas,int nStar)
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
	pkd->nDark = nDark;
	pkd->nGas = nGas;
	pkd->nStar = nStar;
	pkd->nMaxOrderGas = nGas - 1;
	pkd->nMaxOrderDark = nGas + nDark - 1;
	pkd->nRejects = 0;
	for (j=0;j<3;++j) {
		pkd->fPeriod[j] = fPeriod[j];
		}
	/*
	 ** Allocate the main particle store.
	 ** Need to use mdlMalloc() since the particles will need to be
	 ** visible to all other processors thru mdlAquire() later on.
	 **
	 ** We need one EXTRA storage location at the very end to use for 
	 ** calculating acceleration on arbitrary positions in space, for example
	 ** determining the force on the sun. The easiest way to do this is to
	 ** allocate one hidden particle, which won't interfere with the rest of
	 ** the code (hopefully). pkd->pStore[pkd->nStore] is this particle.
	 */
	pkd->pStore = mdlMalloc(pkd->mdl,(nStore+1)*sizeof(PARTICLE));
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
	if (pkd->kdNodes) {
		/*
		 ** Close caching space and free up nodes.
		 */
		mdlFinishCache(pkd->mdl,CID_CELL);
		mdlFree(pkd->mdl,pkd->kdNodes);
		}
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


void pkdSeek(PKD pkd,FILE *fp,int nStart,int bStandard)
{
	long lStart;

	/*
	 ** Seek according to true XDR size structures when bStandard is true.
	 ** This may be a bit dicey, but it should work as long
	 ** as no one changes the tipsy binary format!
	 */
	if (bStandard) lStart = 32;
	else lStart = sizeof(struct dump);
	if (nStart > pkd->nGas) {
		if (bStandard) lStart += pkd->nGas*48;
		else lStart += pkd->nGas*sizeof(struct gas_particle);
		nStart -= pkd->nGas;
		if (nStart > pkd->nDark) {
			if (bStandard) lStart += pkd->nDark*36;
			else lStart += pkd->nDark*sizeof(struct dark_particle);
			nStart -= pkd->nDark;
			if (bStandard) lStart += nStart*44;
			lStart += nStart*sizeof(struct star_particle);
			}
		else {
			if (bStandard) lStart += nStart*36;
			else lStart += nStart*sizeof(struct dark_particle);
			}
		}
	else {
		if (bStandard) lStart += nStart*48;
		else lStart += nStart*sizeof(struct gas_particle);
		}
	fseek(fp,lStart,0);
	}


void pkdReadTipsy(PKD pkd,char *pszFileName,int nStart,int nLocal,
				  int bStandard,double dvFac,double dTuFac)
{
	FILE *fp;
	int i,j;
	PARTICLE *p;
	struct dark_particle dp;
	struct gas_particle gp;
	struct star_particle sp;
	float fTmp;

	pkd->nLocal = nLocal;
	pkd->nActive = nLocal;
	/*
	 ** General initialization.
	 */
	for (i=0;i<nLocal;++i) {
		p = &pkd->pStore[i];
		TYPEClear(p);
		p->iTouchRung = 0;
		p->iRung = 0;
		p->fWeight = 1.0;
		p->fDensity = 0.0;
		p->fBall2 = 0.0;
#ifdef GASOLINE
		for (j=0;j<3;++j) {
			p->vPred[j] = 0.0;
			}
		p->u = 0.0;
		p->uPred = 0.0;
		p->Y_HI = (pkd->cl).Y_H;
		p->Y_HeI = (pkd->cl).Y_He;
		p->Y_HeII = 0.0;
		
		p->c = 0.0;
		p->fMetals = 0.0;
		p->fTimeForm = 0.0;
#endif
		}
	/*
	 ** Seek past the header and up to nStart.
	 */
	fp = fopen(pszFileName,"r");
	assert(fp != NULL);
	/*
	 ** Seek to right place in file.
	 */
	pkdSeek(pkd,fp,nStart,bStandard);
	/*
	 ** Read Stuff!
	 */
	if (bStandard) {
		FLOAT vTemp;
		XDR xdrs;
		xdrstdio_create(&xdrs,fp,XDR_DECODE);
		for (i=0;i<nLocal;++i) {
			p = &pkd->pStore[i];
			p->iOrder = nStart + i;
			if (pkdIsDark(pkd,p)) {
				xdr_float(&xdrs,&fTmp);
				p->fMass = fTmp;
				for (j=0;j<3;++j) {
					xdr_float(&xdrs,&fTmp);
					p->r[j] = fTmp;
					}
				for (j=0;j<3;++j) {
					xdr_float(&xdrs,&fTmp);
					vTemp = fTmp;
					p->v[j] = dvFac*vTemp;			
					}
				xdr_float(&xdrs,&fTmp);
				p->fSoft = fTmp;
				xdr_float(&xdrs,&fTmp);
				p->fPot = fTmp;
				}
			else if (pkdIsGas(pkd,p)) {
				xdr_float(&xdrs,&fTmp);
				p->fMass = fTmp;
				for (j=0;j<3;++j) {
					xdr_float(&xdrs,&fTmp);
					p->r[j] = fTmp;
					}
				for (j=0;j<3;++j) {
					xdr_float(&xdrs,&fTmp);
					vTemp = fTmp;
					p->v[j] = dvFac*vTemp;			
#ifdef GASOLINE
					p->vPred[j] = dvFac*vTemp;
#endif
					}
#ifdef GASOLINE
				xdr_float(&xdrs,&fTmp);
				p->fDensity = fTmp;
				/*
				 ** Convert Temperature to Thermal energy.
				 */
				xdr_float(&xdrs,&fTmp);
				vTemp = fTmp;
				p->u = dTuFac*vTemp;
				p->uPred = dTuFac*vTemp;
				xdr_float(&xdrs,&fTmp);
				p->fSoft = fTmp;
				xdr_float(&xdrs,&fTmp);
				p->fMetals = fTmp;
#else
				xdr_float(&xdrs,&fTmp);
				xdr_float(&xdrs,&fTmp);
				xdr_float(&xdrs,&fTmp);
				p->fSoft = fTmp;
				xdr_float(&xdrs,&fTmp);
#endif
				xdr_float(&xdrs,&fTmp);
				p->fPot = fTmp;
				}
			else if (pkdIsStar(pkd,p)) {
				xdr_float(&xdrs,&fTmp);
				p->fMass = fTmp;
				for (j=0;j<3;++j) {
					xdr_float(&xdrs,&fTmp);
					p->r[j] = fTmp;
					}
				for (j=0;j<3;++j) {
					xdr_float(&xdrs,&fTmp);
					vTemp = fTmp;
					p->v[j] = dvFac*vTemp;			
					}
#ifdef GASOLINE
				xdr_float(&xdrs,&fTmp);
				p->fMetals = fTmp;
				xdr_float(&xdrs,&fTmp);
				p->fTimeForm = fTmp;
#else
				xdr_float(&xdrs,&fTmp);
				xdr_float(&xdrs,&fTmp);
#endif
				xdr_float(&xdrs,&fTmp);
				p->fSoft = fTmp;
				xdr_float(&xdrs,&fTmp);
				p->fPot = fTmp;
				}
			else assert(0);
			}
		xdr_destroy(&xdrs);
		}
	else {
		for (i=0;i<nLocal;++i) {
			p = &pkd->pStore[i];
			p->iOrder = nStart + i;
			if (pkdIsDark(pkd,p)) {
				fread(&dp,sizeof(struct dark_particle),1,fp);
				for (j=0;j<3;++j) {
					p->r[j] = dp.pos[j];
					p->v[j] = dvFac*dp.vel[j];
					}
				p->fMass = dp.mass;
				p->fSoft = dp.eps;
				p->fPot = dp.phi;
				}
			else if (pkdIsGas(pkd,p)) {
				fread(&gp,sizeof(struct gas_particle),1,fp);
				for (j=0;j<3;++j) {
					p->r[j] = gp.pos[j];
					p->v[j] = dvFac*gp.vel[j];
#ifdef GASOLINE
					p->vPred[j] = dvFac*gp.vel[j];
#endif
					}
				p->fMass = gp.mass;
				p->fSoft = gp.hsmooth;
				p->fPot = gp.phi;
#ifdef GASOLINE
				p->fDensity = gp.rho;
				p->u = dTuFac*gp.temp;
				p->uPred = dTuFac*gp.temp;
				p->fMetals = gp.metals;
#endif
				}
			else if (pkdIsStar(pkd,p)) {
				fread(&sp,sizeof(struct star_particle),1,fp);
				for (j=0;j<3;++j) {
					p->r[j] = sp.pos[j];
					p->v[j] = dvFac*sp.vel[j];
					}
				p->fMass = sp.mass;
				p->fSoft = sp.eps;
				p->fPot = sp.phi;
#ifdef GASOLINE
				p->fMetals = sp.metals;
				p->fTimeForm = sp.tform;		
#endif
				}
			else assert(0);
			}
		}
	fclose(fp);
	}


void pkdCalcBound(PKD pkd,BND *pbnd,BND *pbndActive,BND *pbndTreeActive)
{
	int i,j;

	/*
	 ** Initialize the bounds to 0 at the beginning
	 */
	for (j=0;j<3;++j) {
		pbnd->fMin[j] = FLOAT_MAXVAL;
		pbnd->fMax[j] = -FLOAT_MAXVAL;
		pbndActive->fMin[j] = FLOAT_MAXVAL;
		pbndActive->fMax[j] = -FLOAT_MAXVAL;
		pbndTreeActive->fMin[j] = FLOAT_MAXVAL;
		pbndTreeActive->fMax[j] = -FLOAT_MAXVAL;
		}
	/*
	 ** Calculate Local Bounds.
	 */
	for (i=0;i<pkd->nLocal;++i) {
		for (j=0;j<3;++j) {
			if (pkd->pStore[i].r[j] < pbnd->fMin[j]) 
				pbnd->fMin[j] = pkd->pStore[i].r[j];
			if (pkd->pStore[i].r[j] > pbnd->fMax[j])
				pbnd->fMax[j] = pkd->pStore[i].r[j];
			}
		}
	/*
	 ** Calculate Active Bounds.
	 */
	for (i=0;i<pkd->nLocal;++i) {
		if (TYPEQueryACTIVE(&(pkd->pStore[i]))) {
			for (j=0;j<3;++j) {
				if (pkd->pStore[i].r[j] < pbndActive->fMin[j]) 
					pbndActive->fMin[j] = pkd->pStore[i].r[j];
				if (pkd->pStore[i].r[j] > pbndActive->fMax[j])
					pbndActive->fMax[j] = pkd->pStore[i].r[j];
				}
			}
		}
	/*
	 ** Calculate TreeActive Bounds.
	 */
	for (i=0;i<pkd->nLocal;++i) {
		if (TYPEQueryTREEACTIVE(&(pkd->pStore[i]))) {
			for (j=0;j<3;++j) {
				if (pkd->pStore[i].r[j] < pbndTreeActive->fMin[j]) 
					pbndTreeActive->fMin[j] = pkd->pStore[i].r[j];
				if (pkd->pStore[i].r[j] > pbndTreeActive->fMax[j])
					pbndTreeActive->fMax[j] = pkd->pStore[i].r[j];
				}
			}
		}
	}


void pkdGasWeight(PKD pkd)
{
    PARTICLE *p;
    int i, nActive = 0;

    for(i=0;i<pkdLocal(pkd);++i) {
                p = &pkd->pStore[i];
		if (TYPETest( p, TYPE_GAS )) p->fWeight=1.0;
		else p->fWeight=0.0;
		}
    }

/*
 ** Partition particles between iFrom and iTo into those < fSplit and
 ** those >= to fSplit.  Find number and weight in each partition.
 */
int pkdWeight(PKD pkd,int d,FLOAT fSplit,int iSplitSide,int iFrom,int iTo,
			  int *pnLow,int *pnHigh,FLOAT *pfLow,FLOAT *pfHigh)
{
	int i,iPart;
	FLOAT fLower,fUpper;

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


int pkdOrdWeight(PKD pkd,int iOrdSplit,int iSplitSide,int iFrom,int iTo,
				 int *pnLow,int *pnHigh)
{
	int iPart;
	
	/*
	 ** First partition the memory about fSplit for particles iFrom to iTo.
	 */
	if (iSplitSide) {
		iPart = pkdLowerOrdPart(pkd,iOrdSplit,iFrom,iTo);
		*pnLow = pkdLocal(pkd)-iPart;
		*pnHigh = iPart;
		}
	else {
		iPart = pkdUpperOrdPart(pkd,iOrdSplit,iFrom,iTo);
		*pnLow = iPart;
		*pnHigh = pkdLocal(pkd)-iPart;
		}
	return(iPart);
	}


int pkdLowerPart(PKD pkd,int d,FLOAT fSplit,int i,int j)
{
	PARTICLE pTemp;

#ifdef COLLISIONS
	assert(d < 5);
	if(d > 2)
	    return pkdLowerQQPart(pkd, d, fSplit, i, j);
#else /* COLLISIONS */
	assert(d < 3);
#endif /* !COLLISIONS */
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


int pkdUpperPart(PKD pkd,int d,FLOAT fSplit,int i,int j)
{
	PARTICLE pTemp;
	int i0,j0,k;

	i0=i; j0=j;
        if (i0==0 && j0==3557) {
	  k=0;
	  }

#ifdef COLLISIONS
	assert(d < 5);
	if(d > 2)
	    return pkdUpperQQPart(pkd, d, fSplit, i, j);
#else /* COLLISIONS */
	assert(d < 3);
#endif /* !COLLISIONS */
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


int pkdActiveTypeOrder(PKD pkd, unsigned int iTestMask)
{
	PARTICLE pTemp;
	int i=0;
	int j=pkdLocal(pkd)-1;

	if (i > j) goto done;
    while (1) {
        while (TYPETest(&(pkd->pStore[i]), iTestMask ))
            if (++i > j) goto done;
        while (!TYPETest(&(pkd->pStore[j]), iTestMask ))
            if (i > --j) goto done;
		pTemp = pkd->pStore[i];
		pkd->pStore[i] = pkd->pStore[j];
		pkd->pStore[j] = pTemp;
        }
 done:
        if ( iTestMask & TYPE_ACTIVE )       pkd->nActive = i;
        if ( iTestMask & TYPE_TREEACTIVE )   pkd->nTreeActive = i;
        if ( iTestMask & TYPE_SMOOTHACTIVE ) pkd->nSmoothActive = i;
	return (i);
        }


int pkdActiveOrder(PKD pkd)
{
	PARTICLE pTemp;
	int i=0;
	int j=pkdLocal(pkd)-1;

	if (i > j) goto done;
    while (1) {
        while (TYPEQueryACTIVE(&(pkd->pStore[i])))
            if (++i > j) goto done;
        while (!TYPEQueryACTIVE(&(pkd->pStore[j])))
            if (i > --j) goto done;
		pTemp = pkd->pStore[i];
		pkd->pStore[i] = pkd->pStore[j];
		pkd->pStore[j] = pTemp;
        }
 done:
	return (pkd->nActive = i);
	}


int pkdColRejects(PKD pkd,int d,FLOAT fSplit,FLOAT fSplitInactive,
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
	    assert(TYPEQueryACTIVE(&(pkd->pStore[i])));
	for(i = pkdActive(pkd); i < nSplitInactive; ++i)
	    assert(!TYPEQueryACTIVE(&(pkd->pStore[i])));

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

int pkdTreeActive(PKD pkd)
{
	return(pkd->nTreeActive);
	}

int pkdSmoothActive(PKD pkd)
{
	return(pkd->nSmoothActive);
	}

int pkdInactive(PKD pkd)
{
	return(pkd->nLocal - pkd->nActive);
	}

int pkdTreeInactive(PKD pkd)
{
	return(pkd->nLocal - pkd->nTreeActive);
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
#ifdef COLORCODE
	int i;
	
	for (i=0;i<pkd->nLocal;++i) {
		pkd->pStore[i].fColor = (FLOAT)pkd->idSelf;
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


void pkdLocalOrder(PKD pkd)
{
	qsort(pkd->pStore,pkdLocal(pkd),sizeof(PARTICLE),cmpParticles);
	}


void pkdWriteTipsy(PKD pkd,char *pszFileName,int nStart,
				   int bStandard,double dvFac,double duTFac,int iGasModel)
{
	PARTICLE *p;
	FILE *fp;
	int i,j;
	struct dark_particle dp;
	struct gas_particle gp;
	struct star_particle sp;
	int nout;
	float fTmp;

	/*
	 ** Seek past the header and up to nStart.
	 */
	fp = fopen(pszFileName,"r+");
	assert(fp != NULL);
	pkdSeek(pkd,fp,nStart,bStandard);
	if (bStandard) {
		FLOAT vTemp;
		XDR xdrs;
		/* 
		 ** Write Stuff!
		 */
		xdrstdio_create(&xdrs,fp,XDR_ENCODE);
		for (i=0;i<pkdLocal(pkd);++i) {
			p = &pkd->pStore[i];
			if (pkdIsDark(pkd,p)) {
				fTmp = p->fMass;
				xdr_float(&xdrs,&fTmp);
				for (j=0;j<3;++j) {
					fTmp = p->r[j];
					xdr_float(&xdrs,&fTmp);
					}
				for (j=0;j<3;++j) {
					vTemp = dvFac*p->v[j];			
					fTmp = vTemp;
					xdr_float(&xdrs,&fTmp);
					}
				fTmp = p->fSoft;
				xdr_float(&xdrs,&fTmp);
				fTmp = p->fPot;
				xdr_float(&xdrs,&fTmp);
				}
			else if (pkdIsGas(pkd,p)) {
				fTmp = p->fMass;
				xdr_float(&xdrs,&fTmp);
				for (j=0;j<3;++j) {
					fTmp = p->r[j];
					xdr_float(&xdrs,&fTmp);
					}
				for (j=0;j<3;++j) {
					vTemp = dvFac*p->v[j];			
					fTmp = vTemp;
					xdr_float(&xdrs,&fTmp);
					}
				fTmp = p->fDensity;
				xdr_float(&xdrs,&fTmp);
#ifdef GASOLINE
				/*
				 ** Convert thermal energy to tempertature.
				 */
			        switch (iGasModel) {
				case GASMODEL_COOLING:
				case GASMODEL_COOLING_NONEQM:
				        vTemp = clTemperature( 2*(pkd->cl).Y_H - p->Y_HI + 
							       3*(pkd->cl).Y_He - 2*p->Y_HeI - p->Y_HeII,
							       p->u * (pkd->cl).dErgPerGmUnit );
					break;
				default:
				        vTemp = duTFac*p->u;
				        }
				fTmp = vTemp;
				xdr_float(&xdrs,&fTmp);
				/* fTmp = sqrt(0.25*p->fBall2);  Write softening in tipsy outputs */
                                fTmp = p->fSoft;
				xdr_float(&xdrs,&fTmp);
#ifdef DEBUG
				/* Store divv in metals for now */
				fTmp = p->divv;
#else
				fTmp = p->fMetals;
#endif
				xdr_float(&xdrs,&fTmp);
#else
				fTmp = 0.0;
				xdr_float(&xdrs,&fTmp);
				fTmp = p->fSoft;
				xdr_float(&xdrs,&fTmp);
				fTmp = 0.0;
				xdr_float(&xdrs,&fTmp);
#endif
				fTmp = p->fPot;
				xdr_float(&xdrs,&fTmp);
				}
			else if (pkdIsStar(pkd,p)) {
				fTmp = p->fMass;
				xdr_float(&xdrs,&fTmp);
				for (j=0;j<3;++j) {
					fTmp = p->r[j];
					xdr_float(&xdrs,&fTmp);
					}
				for (j=0;j<3;++j) {
					vTemp = dvFac*p->v[j];			
					fTmp = vTemp;
					xdr_float(&xdrs,&fTmp);
					}
#ifdef GASOLINE
				fTmp = p->fMetals;
				xdr_float(&xdrs,&fTmp);
				fTmp = p->fTimeForm;
				xdr_float(&xdrs,&fTmp);
#else
				fTmp = 0.0;
				xdr_float(&xdrs,&fTmp);
				xdr_float(&xdrs,&fTmp);			
#endif
				fTmp = p->fSoft;
				xdr_float(&xdrs,&fTmp);
				fTmp = p->fPot;
				xdr_float(&xdrs,&fTmp);
				}
			else assert(0);
			}
		xdr_destroy(&xdrs);
		}
	else {
		/* 
		 ** Write Stuff!
		 */
		for (i=0;i<pkdLocal(pkd);++i) {
			p = &pkd->pStore[i];
			if (pkdIsDark(pkd,p)) {
				for (j=0;j<3;++j) {
					dp.pos[j] = p->r[j];
					dp.vel[j] = dvFac*p->v[j];
					}
				dp.mass = p->fMass;
				dp.eps = p->fSoft;
				dp.phi = p->fPot;
				nout = fwrite(&dp,sizeof(struct dark_particle),1,fp);
				assert(nout == 1);
				}
			else if (pkdIsGas(pkd,p)) {
				for (j=0;j<3;++j) {
					gp.pos[j] = p->r[j];
					gp.vel[j] = dvFac*p->v[j];
					}
				gp.mass = p->fMass;
				gp.hsmooth = p->fSoft;
				gp.phi = p->fPot;
				gp.rho = p->fDensity;
#ifdef GASOLINE
				gp.temp = duTFac*p->u;
				gp.metals = p->fMetals;
#else
				gp.temp = 0.0;
				gp.metals = 0.0;
#endif
				nout = fwrite(&gp,sizeof(struct gas_particle),1,fp);
				assert(nout == 1);
				}
			else if (pkdIsStar(pkd,p)) {
				for (j=0;j<3;++j) {
					sp.pos[j] = p->r[j];
					sp.vel[j] = dvFac*p->v[j];
					}
				sp.mass = p->fMass;
				sp.eps = p->fSoft;
				sp.phi = p->fPot;
#ifdef GASOLINE
				sp.metals = p->fMetals;
				sp.tform = p->fTimeForm;
#else
				sp.metals = 0.0;
				sp.tform = 0.0;
#endif
				nout = fwrite(&sp,sizeof(struct dark_particle),1,fp);
				assert(nout == 1);
				}
			else assert(0);
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
	pOut->fSoft = p1->fMass*p1->fSoft + p2->fMass*p2->fSoft;
	for (j=0;j<3;++j) {
		pOut->r[j] = p1->fMass*p1->r[j] + p2->fMass*p2->r[j];
		}
	if (pOut->fMass > 0) {
		pOut->fSoft /= pOut->fMass;
		for (j=0;j<3;++j) {
			pOut->r[j] /= pOut->fMass;
			}
		}
	}


void pkdCalcCell(PKD pkd,KDN *pkdn,FLOAT *rcm,int iOrder,
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
		pcc->Hxxzz = 0.0;
		pcc->Hxyzz = 0.0;
		pcc->Hxzzz = 0.0;
		pcc->Hyyzz = 0.0;
		pcc->Hyzzz = 0.0;
		pcc->Hzzzz = 0.0;
		pcc->B6 = 0.0;
	case 3:
		pcc->Oxxx = 0.0;
		pcc->Oxyy = 0.0;
		pcc->Oxxy = 0.0;
		pcc->Oyyy = 0.0;
		pcc->Oxxz = 0.0;
		pcc->Oyyz = 0.0;
		pcc->Oxyz = 0.0;
		pcc->Oxzz = 0.0;
		pcc->Oyzz = 0.0;
		pcc->Ozzz = 0.0;
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
	if (pkdn == NULL) pkdn = &pkd->kdNodes[pkd->iRoot];
	for (pj=pkdn->pLower;pj<=pkdn->pUpper;++pj) {
		m = pkd->pStore[pj].fMass;
		dx = pkd->pStore[pj].r[0] - rcm[0];
		dy = pkd->pStore[pj].r[1] - rcm[1];
		dz = pkd->pStore[pj].r[2] - rcm[2];
		d2 = dx*dx + dy*dy + dz*dz;
		d1 = sqrt(d2);
		if (d1 > pcc->Bmax) pcc->Bmax = d1;
		pcc->B2 += m*d2;
		pcc->B3 += m*d2*d1;
		pcc->B4 += m*d2*d2;
		pcc->B5 += m*d2*d2*d1;
		pcc->B6 += m*d2*d2*d2;
#ifdef COMPLETE_LOCAL
		d2 = 0.0;
#endif
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
			pcc->Hxxzz += m*(dx*dx*dz*dz - 1.0/7.0*d2*(dx*dx + dz*dz - 0.2*d2));
			pcc->Hxyzz += m*(dx*dy*dz*dz - 1.0/7.0*d2*dx*dy);
			pcc->Hxzzz += m*(dx*dz*dz*dz - 3.0/7.0*d2*dx*dz);
			pcc->Hyyzz += m*(dy*dy*dz*dz - 1.0/7.0*d2*(dy*dy + dz*dz - 0.2*d2));
			pcc->Hyzzz += m*(dy*dz*dz*dz - 3.0/7.0*d2*dy*dz);
			pcc->Hzzzz += m*(dz*dz*dz*dz - 6.0/7.0*d2*(dz*dz - 0.1*d2));
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
			pcc->Oxzz += m*(dx*dz*dz - 0.2*d2*dx);
			pcc->Oyzz += m*(dy*dz*dz - 0.2*d2*dy);
			pcc->Ozzz += m*(dz*dz*dz - 0.6*d2*dz);;
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


double fcnAbsQuad(KDN *pkdn,double r)
{
	double t;

	t = r*(r - pkdn->mom.Bmax);
	t *= r*t;
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


double fcnAbsHex(KDN *pkdn,double r)
{
	double t;

	t = r*r*(r - pkdn->mom.Bmax);
	t *= r*t;
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
		default:
		        assert(0);
			}
		}
	else if (iOpenType == OPEN_JOSH) {
		/*
		 ** Set openening criterion to an approximation of Josh's theta.
		 */
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


void pkdUpPass(PKD pkd,int iCell,int iOpenType,double dCrit,
	       int iOrder, int bGravity)
{
	KDN *c;
	PARTICLE *p;
	int l,u,pj,j;
	double dOpen;

	c = pkd->kdNodes;
	p = pkd->pStore;
	l = c[iCell].pLower;
	u = c[iCell].pUpper;
	if (c[iCell].iDim >= 0) {
		pkdUpPass(pkd,LOWER(iCell),iOpenType,dCrit,iOrder, bGravity);
		pkdUpPass(pkd,UPPER(iCell),iOpenType,dCrit,iOrder, bGravity);
		pkdCombine(&c[LOWER(iCell)],&c[UPPER(iCell)],&c[iCell]);
		}
	else {
		c[iCell].fMass = 0.0;
		c[iCell].fSoft = 0.0;
		for (j=0;j<3;++j) {
			/*
			 ** We initialize the bounds to these extreme crazy values so
			 ** that any particle will set the bounds correctly, if there 
			 ** are no particles in the loop then these remain the bounds.
			 */
			c[iCell].bnd.fMin[j] = FLOAT_MAXVAL;
			c[iCell].bnd.fMax[j] = -FLOAT_MAXVAL;
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
		if (c[iCell].fMass > 0) {
			for (j=0;j<3;++j) {
				c[iCell].r[j] /= c[iCell].fMass;
				}
			c[iCell].fSoft /= c[iCell].fMass;
			}
		}
	/*
	 ** Calculate multipole moments.
	 */
	if(bGravity) {
	    pkdCalcCell(pkd,&c[iCell],c[iCell].r,iOrder,&c[iCell].mom);
	    dOpen = pkdCalcOpen(&c[iCell],iOpenType,dCrit,iOrder);
	    c[iCell].fOpen2 = dOpen*dOpen;
	    }
	}


/*
 ** JST's Select
 */
void pkdSelect(PKD pkd,int d,int k,int l,int r)
{
	PARTICLE *p,t;
	FLOAT v;
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


int NumBinaryNodes(PKD pkd,int nBucket,int pLower,int pUpper)
{
	BND bnd;
	int i,j,m,d,l,u;
	FLOAT fSplit;
	int bGoodBounds;	/* is the cell a finite size? */

	if (pLower > pUpper) return(0);
	else {
		/*
		 ** We need to find the bounding box.
		 */
		for (j=0;j<3;++j) {
			bnd.fMin[j] = pkd->pStore[pLower].r[j];
			bnd.fMax[j] = pkd->pStore[pLower].r[j];
			}
		for (i=pLower+1;i<=pUpper;++i) {
			for (j=0;j<3;++j) {
				if (pkd->pStore[i].r[j] < bnd.fMin[j]) 
					bnd.fMin[j] = pkd->pStore[i].r[j];
				else if (pkd->pStore[i].r[j] > bnd.fMax[j])
					bnd.fMax[j] = pkd->pStore[i].r[j];
				}
			}	
		bGoodBounds = 0;
		for (j=0;j<3;++j) {
			if (bnd.fMax[j] > bnd.fMin[j])
			    bGoodBounds = 1;
			}
		if ((pUpper-pLower+1 > nBucket) && bGoodBounds) {
			/*
			 ** Now we need to determine the longest axis.
			 */
			d = 0;
			for (j=1;j<3;++j) {
				if (bnd.fMax[j]-bnd.fMin[j] > bnd.fMax[d]-bnd.fMin[d]) {
					d = j;
					}
				}
			/*
			 ** Now we do the split.
			 */
			fSplit = 0.5*(bnd.fMin[d]+bnd.fMax[d]);
			m = pkdUpperPart(pkd,d,fSplit,pLower,pUpper);
			/*
			 ** Recursive node count.
			 */
			l = NumBinaryNodes(pkd,nBucket,pLower,m-1);
			u = NumBinaryNodes(pkd,nBucket,m,pUpper);
			/*
			 ** Careful, this assert only applies when we are doing the
			 ** squeezing!
			 */
			assert(l > 0 && u > 0);
			return(l+u+1);
			}
		else {
			return(1);
			}
		}
	}


int BuildBinary(PKD pkd,int nBucket,int pLower,int pUpper,int iOpenType,
		double dCrit,int iOrder, int bGravity)
{
	KDN *pkdn;
	int i,j,m,d,c;
	FLOAT fm;
	double dOpen;
	int bGoodBounds;	/* Is the cell a finite size? */

	if (pLower > pUpper) return(-1);
	else {
		/*
		 ** Grab a cell from the cell storage!
		 */
		assert(pkd->iFreeCell < pkd->nNodes);
		c = pkd->iFreeCell++;
		pkdn = &pkd->kdNodes[c];
		pkdn->pLower = pLower;
		pkdn->pUpper = pUpper;
		/*
		 ** We need to find the bounding box.
		 */
		for (j=0;j<3;++j) {
			pkdn->bnd.fMin[j] = pkd->pStore[pLower].r[j];
			pkdn->bnd.fMax[j] = pkd->pStore[pLower].r[j];
			}
		for (i=pLower+1;i<=pUpper;++i) {
			for (j=0;j<3;++j) {
				if (pkd->pStore[i].r[j] < pkdn->bnd.fMin[j]) 
					pkdn->bnd.fMin[j] = pkd->pStore[i].r[j];
				else if (pkd->pStore[i].r[j] > pkdn->bnd.fMax[j])
					pkdn->bnd.fMax[j] = pkd->pStore[i].r[j];
				}
			}	
		bGoodBounds = 0;
		for (j=0;j<3;++j) {
			if (pkdn->bnd.fMax[j] > pkdn->bnd.fMin[j])
			    bGoodBounds = 1;
			}
		if ((pUpper-pLower+1 > nBucket) && bGoodBounds) {
			/*
			 ** Now we need to determine the longest axis.
			 */
			d = 0;
			for (j=1;j<3;++j) {
				if (pkdn->bnd.fMax[j]-pkdn->bnd.fMin[j] > 
					pkdn->bnd.fMax[d]-pkdn->bnd.fMin[d]) {
					d = j;
					}
				}
			pkdn->iDim = d;
			/*
			 ** Now we do the split.
			 */
			pkdn->fSplit = 0.5*(pkdn->bnd.fMin[d]+pkdn->bnd.fMax[d]);
			m = pkdUpperPart(pkd,d,pkdn->fSplit,pLower,pUpper);
			pkdn->iLower = BuildBinary(pkd,nBucket,pLower,m-1,
						   iOpenType, dCrit,iOrder,
						   bGravity);
			pkdn->iUpper = BuildBinary(pkd,nBucket,m,pUpper,
						   iOpenType, dCrit,iOrder,
						   bGravity);
			/*
			 ** Careful, this assert only applies when we are doing the
			 ** squeezing!
			 */
			assert(pkdn->iLower != -1 && pkdn->iUpper != -1);
			/*
			 ** Now calculate the mass, center of mass and mass weighted
			 ** softening radius.
			 */
			pkdn->fMass = 0.0;
			pkdn->fSoft = 0.0;
			for (j=0;j<3;++j) {
				pkdn->r[j] = 0.0;
				}
			if (pkdn->iLower != -1) {
				fm = pkd->kdNodes[pkdn->iLower].fMass;
				pkdn->fMass += fm;
				pkdn->fSoft += fm*pkd->kdNodes[pkdn->iLower].fSoft;
				for (j=0;j<3;++j) {
					pkdn->r[j] += fm*pkd->kdNodes[pkdn->iLower].r[j];
					}
				}
			if (pkdn->iUpper != -1) {
				fm = pkd->kdNodes[pkdn->iUpper].fMass;
				pkdn->fMass += fm;
				pkdn->fSoft += fm*pkd->kdNodes[pkdn->iUpper].fSoft;
				for (j=0;j<3;++j) {
					pkdn->r[j] += fm*pkd->kdNodes[pkdn->iUpper].r[j];
					}
				}
			if (pkdn->fMass > 0) {
				pkdn->fSoft /= pkdn->fMass;
				for (j=0;j<3;++j) {
					pkdn->r[j] /= pkdn->fMass;
					}
				}
			}
		else {
			pkdn->iDim = -1; /* it is a bucket! */
			pkdn->fSplit = 0.0;
			pkdn->iLower = -1;
			pkdn->iUpper = -1;
			/*
			 ** Calculate the bucket quantities.
			 */
			pkdn->fMass = 0.0;
			pkdn->fSoft = 0.0;
			for (j=0;j<3;++j) {
				pkdn->r[j] = 0.0;
				}
			for (i=pkdn->pLower;i<=pkdn->pUpper;++i) {
				fm = pkd->pStore[i].fMass;
				pkdn->fMass += fm;
				pkdn->fSoft += fm*pkd->pStore[i].fSoft;
				for (j=0;j<3;++j) {
					pkdn->r[j] += fm*pkd->pStore[i].r[j];
					}
				}
			if (pkdn->fMass > 0) {
				pkdn->fSoft /= pkdn->fMass;
				for (j=0;j<3;++j) {
					pkdn->r[j] /= pkdn->fMass;
					}
				}
			}
		/*
		 ** Calculate multipole moments.
		 */
		if(bGravity) {
		    pkdCalcCell(pkd,pkdn,pkdn->r,iOrder,&pkdn->mom);
		    dOpen = pkdCalcOpen(pkdn,iOpenType,dCrit,iOrder);
		    pkdn->fOpen2 = dOpen*dOpen;
		    }
		return(c);
		}
	}


void pkdThreadTree(PKD pkd,int iCell,int iNext)
{
	int l,u;

	if (iCell == -1) return;
	else if (pkd->kdNodes[iCell].iDim != -1) {
		l = pkd->kdNodes[iCell].iLower;
		u = pkd->kdNodes[iCell].iUpper;
		if (u == -1) {
			assert(l != -1);
			pkdThreadTree(pkd,l,iNext);
			}
		else if (l == -1) {
			assert(u != -1);
			pkdThreadTree(pkd,u,iNext);
			/* 
			 ** It is convenient to change the "down" pointer in this case.
			 */
			pkd->kdNodes[iCell].iLower = u;
			}
		else {
			pkdThreadTree(pkd,l,u);
			pkdThreadTree(pkd,u,iNext);
			}
		pkd->kdNodes[iCell].iUpper = iNext;
		}
	else {
		pkd->kdNodes[iCell].iLower = -1;	/* Just make sure! */
		pkd->kdNodes[iCell].iUpper = iNext;
		}
	}


void pkdBuildBinary(PKD pkd,int nBucket,int iOpenType,double dCrit,
		    int iOrder,int bTreeActiveOnly,int bGravity, KDN *pRoot)
{
	/*
	 ** Make sure the particles are in Active/Inactive order.
	 */
	pkdActiveTypeOrder(pkd, TYPE_ACTIVE|TYPE_TREEACTIVE );
	if (pkd->kdNodes) {
		/*
		 ** Close caching space and free up nodes.
		 */
		mdlFinishCache(pkd->mdl,CID_CELL);
		mdlFree(pkd->mdl,pkd->kdNodes);
		}
	/*
	 ** First problem is to figure out how many cells we need to
	 ** allocate. We need to do this in advance because we need to
	 ** synchronize to allocate the memory with mdlMalloc().
	 */
	if (bTreeActiveOnly) {
		pkd->nNodes = NumBinaryNodes(pkd,nBucket,0,pkd->nTreeActive-1);
		}
	else {
		pkd->nNodes = NumBinaryNodes(pkd,nBucket,0,pkd->nLocal-1);
		}
	/*
	 ** We need at least one particle per processor.
	 */
	assert(pkd->nNodes > 0);

	/*
	 ** Need to allocate a special extra cell that we will use to calculate
	 ** the acceleration on an arbitrary point in space.
	 */
	pkd->kdNodes = mdlMalloc(pkd->mdl,(pkd->nNodes + 1)*sizeof(KDN));
	assert(pkd->kdNodes != NULL);
	/*
	 ** Now we really build the tree.
	 */
	pkd->iFreeCell = 0;
	if (bTreeActiveOnly) {
		pkd->iRoot = BuildBinary(pkd,nBucket,0,pkd->nTreeActive-1,
					 iOpenType,dCrit,iOrder, bGravity);
		}
	else {
		pkd->iRoot = BuildBinary(pkd,nBucket,0,pkd->nLocal-1,
					 iOpenType,dCrit,iOrder, bGravity);
		}
	assert(pkd->iFreeCell == pkd->nNodes);
	/*
	 ** Thread the tree.
	 */
	pkdThreadTree(pkd,pkd->iRoot,-1);
	*pRoot = pkd->kdNodes[pkd->iRoot];
	/*
	 ** Finally activate a read only cache for remote access.
	 */
	mdlROcache(pkd->mdl,CID_CELL,pkd->kdNodes,sizeof(KDN),pkdNodes(pkd));
	}


void pkdBuildLocal(PKD pkd,int nBucket,int iOpenType,double dCrit,
		   int iOrder,int bTreeActiveOnly,int bGravity, KDN *pRoot)
{
	int l,n,i,d,m,j,diff;
	KDN *c;
	char ach[256];
	BND bndDum;
	int bGoodBounds;

	/*
	 ** Make sure the particles are in Active/Inactive order.
	 */
	pkdActiveTypeOrder(pkd, TYPE_ACTIVE|TYPE_TREEACTIVE );
	pkd->nBucket = nBucket;
	if (bTreeActiveOnly) n = pkd->nTreeActive;
	else n = pkd->nLocal;
	if (n == 0) {
		printf("id:%d has an empty local tree\n",pkd->idSelf);
		}
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
	if (bTreeActiveOnly) {
		c[pkd->iRoot].pUpper = pkd->nTreeActive-1;
		}
	else {
		c[pkd->iRoot].pUpper = pkd->nLocal-1;
		}
	/*
	 ** determine the local bound of the particles.
	 */
	if (bTreeActiveOnly) {
		pkdCalcBound(pkd,&bndDum,&bndDum,&c[pkd->iRoot].bnd);
		}
	else {
		pkdCalcBound(pkd,&c[pkd->iRoot].bnd,&bndDum,&bndDum);
		}
	i = pkd->iRoot;
	while (1) {
		bGoodBounds = 0;
		for (j=0;j<3;++j) {
			if (c[i].bnd.fMax[j] > c[i].bnd.fMin[j])
				bGoodBounds = 1;
			}
		if (i < pkd->nSplit && (c[i].pUpper - c[i].pLower) > 0 && bGoodBounds) {
			d = 0;
			for (j=1;j<3;++j) {
				if (c[i].bnd.fMax[j]-c[i].bnd.fMin[j] > 
					c[i].bnd.fMax[d]-c[i].bnd.fMin[d]) d = j;
				}
			c[i].iDim = d;
			m = (c[i].pLower + c[i].pUpper)/2;
			pkdSelect(pkd,d,m,c[i].pLower,c[i].pUpper);
			c[i].fSplit = pkd->pStore[m].r[d];
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
	pkdUpPass(pkd,pkd->iRoot,iOpenType,dCrit,iOrder, bGravity);
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


void pkdBucketWeight(PKD pkd,int iBucket,FLOAT fWeight)
{
	KDN *pbuc;
	int pj;
	
	pbuc = &pkd->kdNodes[iBucket];
	for (pj=pbuc->pLower;pj<=pbuc->pUpper;++pj) {
		if (TYPEQueryACTIVE(&(pkd->pStore[pj])))
			pkd->pStore[pj].fWeight = fWeight;
		}
	}


void pkdColorCell(PKD pkd,int iCell,FLOAT fColor)
{
#ifdef COLORCODE
	KDN *pkdn;
	int pj;
	
	pkdn = &pkd->kdNodes[iCell];
	for (pj=pkdn->pLower;pj<=pkdn->pUpper;++pj) {
		pkd->pStore[pj].fColor = fColor;
		}
#endif	
	}


void pkdGravAll(PKD pkd,int nReps,int bPeriodic,int iOrder,int iEwOrder,
				double fEwCut,double fEwhCut,int bDoSun,
				double *aSun,int *nActive,
				double *pdPartSum,double *pdCellSum,double *pdSoftSum,
				CASTAT *pcs,double *pdFlop)
{
	KDN *c = pkd->kdNodes;
	int iCell,n;
	FLOAT fWeight;
	double dFlopI, dFlopE;
	int i,j;
	BND bndActive;
	BND bndTmp;
	char achDiag[256];

	*pdFlop = 0.0;
	pkdClearTimer(pkd,1);
	pkdClearTimer(pkd,2);
	pkdClearTimer(pkd,3);
	/*
	 ** Set up Ewald tables and stuff.
	 */
	if (bPeriodic) {
	    pkdEwaldInit(pkd,fEwhCut,iEwOrder);	/* ignored in Flop count! */
		}
	/*
	 ** Start particle caching space (cell cache already active).
	 */
	mdlROcache(pkd->mdl,CID_PARTICLE,pkd->pStore,sizeof(PARTICLE),
			   pkdLocal(pkd));
	/*
	 ** Walk over the local buckets!
	 */
	*nActive = 0;
	*pdPartSum = 0.0;
	*pdCellSum = 0.0;
	*pdSoftSum = 0.0;
	iCell = pkd->iRoot;
	while (iCell != -1) {
		if (c[iCell].iLower != -1) {
			iCell = c[iCell].iLower;
			continue;
			}
		n = 0;
		for(j = 0; j < 3; j++) {
		    bndActive.fMin[j] = FLOAT_MAXVAL;
		    bndActive.fMax[j] = -FLOAT_MAXVAL;
		    }
		for (i=c[iCell].pLower;i<=c[iCell].pUpper;++i) {
			if (TYPEQueryACTIVE(&(pkd->pStore[i]))) {
			    ++n;
			    for (j=0;j<3;++j) {
				if (pkd->pStore[i].r[j] < bndActive.fMin[j]) 
				    bndActive.fMin[j] = pkd->pStore[i].r[j];
				if (pkd->pStore[i].r[j] > bndActive.fMax[j])
				    bndActive.fMax[j] = pkd->pStore[i].r[j];
				}
			    }
			}
		if (n > 0) {
			/* 
			 * set bounds to bounds of active particles.
			 */
			bndTmp = c[iCell].bnd;
			c[iCell].bnd = bndActive;
			/*
			 ** Calculate gravity on this bucket.
			 */
			pkdStartTimer(pkd,1);
			pkdBucketWalk(pkd,iCell,nReps,iOrder); /* ignored in Flop count! */
			pkdStopTimer(pkd,1);
			c[iCell].bnd = bndTmp;
			*nActive += n;
			*pdPartSum += n*pkd->nPart + 
				n*(2*(c[iCell].pUpper-c[iCell].pLower) - n + 1)/2;
			*pdCellSum += n*pkd->nCellNewt;
			*pdSoftSum += n*pkd->nCellSoft;
			pkdStartTimer(pkd,2);
			dFlopI = pkdBucketInteract(pkd,iCell,iOrder);
			*pdFlop += dFlopI;
			pkdStopTimer(pkd,2);
			/*
			 * Now do Ewald part.
			 */
			if (bPeriodic) {
				pkdStartTimer(pkd,3);
				dFlopE = pkdBucketEwald(pkd,iCell,nReps,fEwCut,iEwOrder);
				*pdFlop += dFlopE;
				pkdStopTimer(pkd,3);
				}
			else {
			    dFlopE = 0.0;
			    }
			fWeight = dFlopI + dFlopE;
			/*
			 ** pkdBucketWeight, only updates the weights of the active
			 ** particles. Although this isn't really a requirement it
			 ** might be a good idea, if weights correspond to different 
			 ** tasks at different times.
			 */
			pkdBucketWeight(pkd,iCell,fWeight);
			}
		iCell = c[iCell].iUpper;
		}
	if (bDoSun) {
		const double dTinyBox = 1e-14;
		int iDummy = pkd->nStore;
		/*
		 ** Calculate the indirect interaction for solar system problems.
		 ** we need a "dummy" particle that is pointed to by the cell,
		 ** pkd->iFreeCell which also needs its bounds set correctly for 
		 ** this to work.
		 ** Don't allow periodic BCs at the moment.
		 */
		assert(nReps == 0);
		assert(bPeriodic == 0);
		for (j=0;j<3;++j) {
			pkd->pStore[iDummy].r[j] = 0.0;
			pkd->pStore[iDummy].a[j] = 0.0;
			c[pkd->iFreeCell].bnd.fMin[j] = -dTinyBox;
			c[pkd->iFreeCell].bnd.fMax[j] = dTinyBox;
			}
		TYPESet(&(pkd->pStore[iDummy]),TYPE_ACTIVE);
		pkd->pStore[iDummy].fPot = 0;
		c[pkd->iFreeCell].pLower = iDummy;
		c[pkd->iFreeCell].pUpper = iDummy;
		c[pkd->iFreeCell].iDim = -1;	/* It is really a bucket! */

		pkdBucketWalk(pkd,pkd->iFreeCell,nReps,iOrder);
		pkdBucketInteract(pkd,pkd->iFreeCell,iOrder);

		/*
		 ** Now we should have the indirect acceleration contained in the 
		 ** particle pointed to by iDummy. Now we have to bring this acceleration
		 ** back up to the master level. We can then take care of it with an
		 ** external potential call.
		 */
		for (j=0;j<3;++j) {
			aSun[j] = pkd->pStore[iDummy].a[j];
			}
		}
	/*
	 ** Get caching statistics.
	 */
	pcs->dcNumAccess = mdlNumAccess(pkd->mdl,CID_CELL);
	pcs->dcMissRatio = mdlMissRatio(pkd->mdl,CID_CELL);
	pcs->dcCollRatio = mdlCollRatio(pkd->mdl,CID_CELL);
	pcs->dcMinRatio = mdlMinRatio(pkd->mdl,CID_CELL);
	pcs->dpNumAccess = mdlNumAccess(pkd->mdl,CID_PARTICLE);
	pcs->dpMissRatio = mdlMissRatio(pkd->mdl,CID_PARTICLE);
	pcs->dpCollRatio = mdlCollRatio(pkd->mdl,CID_PARTICLE);
	pcs->dpMinRatio = mdlMinRatio(pkd->mdl,CID_PARTICLE);
	/*
	 ** Stop particle caching space.
	 */
	mdlFinishCache(pkd->mdl,CID_PARTICLE);
	sprintf(achDiag, "nMaxPart: %d, nMaxSoftCell: %d, nMaxNewtCell: %d\n",
		pkd->nMaxPart, pkd->nMaxCellSoft, pkd->nMaxCellNewt);
	mdlDiag(pkd->mdl, achDiag);

	}


void pkdSunIndirect(PKD pkd,double *aSun,int bDoSun,double dSunMass)
{
	PARTICLE *p;
	double t;
	int i,j,n;

	p = pkd->pStore;
	n = pkdLocal(pkd);
	for (i=0;i<n;++i) {
		if (TYPEQueryACTIVE(&(p[i]))) {
			t = 0;
			for (j=0;j<3;++j) t += p[i].r[j]*p[i].r[j];
			t = 1.0/sqrt(t);
			t = t*t*t;
			/*
			 ** The way that aSun[j] gets added in here is confusing, but this
			 ** is the correct way!
			 */
			if (bDoSun) {
				t *= dSunMass;
				for (j=0;j<3;++j) {					
					p[i].a[j] -= aSun[j] + p[i].r[j]*t;
					}				
				}
			else {
				t *= p[i].fMass;
				for (j=0;j<3;++j) {
					p[i].a[j] -= aSun[j] - p[i].r[j]*t;
					}
				}
			}
		}
	}


void pkdLogHalo(PKD pkd)
{
	PARTICLE *p;
	int i,n;

	/* v in (kpc)/(4.691822922e16s) 128 km/s x 1.5233 */
	const double v = 194.9848486;
	const double d = 12.0;	/* in kpc */
	const double fl = 1.0;	/* flattening */
	p = pkd->pStore;
	n = pkdLocal(pkd);
	for (i=0;i<n;++i) {
		if (TYPEQueryACTIVE(&(p[i]))) {
			double x = p[i].r[0];
			double y = p[i].r[1];
			double z = p[i].r[2];
			/*
			 **	Do the logarithmic halo potential.
			 */
			double r = sqrt(x*x + y*y + z*z/(fl*fl));
			double C = 1.0/(r*r + d*d);	
			p[i].a[0] -= 2*v*v*x*C;
			p[i].a[1] -= 2*v*v*y*C;
			p[i].a[2] -= 2*v*v*z*C/(fl*fl);
			p[i].fPot += v*v*log(r*r + d*d);
			}
		}
	}


void pkdHernquistSpheroid(PKD pkd)
{
	PARTICLE *p;
	int i,n;

	const double M_s = 3.4e5;	/* in 10^5 M_sun */
	const double c = 0.7;		/* in kpc */
	p = pkd->pStore;
	n = pkdLocal(pkd);
	for (i=0;i<n;++i) {
		if (TYPEQueryACTIVE(&(p[i]))) {
			double x = p[i].r[0];
			double y = p[i].r[1];
			double z = p[i].r[2];
			/*
			 **	Do the spheroid potential
			 */
			double r = sqrt(x*x + y*y + z*z);
			double A = 1.0/(r + c);	
			p[i].a[0] -= M_s*A*A*x/r;
			p[i].a[1] -= M_s*A*A*y/r;
			p[i].a[2] -= M_s*A*A*z/r;
			p[i].fPot -= M_s/(r + c);
			}
		}
	}


void pkdMiyamotoDisk(PKD pkd)
{
	PARTICLE *p;
	int i,n;

	const double M_d = 1.0e6;	/* in 10^5 M_sun */
	const double aa = 6.5;		/* in kpc */
	const double b = 0.26;		/* in kpc */
	p = pkd->pStore;
	n = pkdLocal(pkd);
	for (i=0;i<n;++i) {
		if (TYPEQueryACTIVE(&(p[i]))) {
			double x = p[i].r[0];
			double y = p[i].r[1];
			double z = p[i].r[2];
			/*
			 **	Do the Miyamoto Disk potential
			 */
			double U = sqrt(z*z + b*b);
			double W = aa + U;
			double A = sqrt(x*x + y*y + W*W);
			p[i].a[0] -= M_d*x/(A*A*A);
			p[i].a[1] -= M_d*y/(A*A*A);
			p[i].a[2] -= M_d*z*W/(A*A*A*U);
			p[i].fPot -= M_d/A;
			}
		}
	}


void pkdCalcE(PKD pkd,double *T,double *U,double *Eth)
{
	PARTICLE *p;
	int i,n;
	FLOAT vx,vy,vz;

	p = pkd->pStore;
	n = pkdLocal(pkd);
	*T = 0.0;
	*U = 0.0;
	*Eth = 0.0;
	for (i=0;i<n;++i) {
		*U += 0.5*p[i].fMass*p[i].fPot;
		vx = p[i].v[0];
		vy = p[i].v[1];
		vz = p[i].v[2];
		*T += 0.5*p[i].fMass*(vx*vx + vy*vy + vz*vz);
#ifdef GASOLINE
		if (pkdIsGas(pkd,&p[i]))
		         *Eth += p[i].fMass*p[i].u;
#endif
		}
	}


/* 
 * Use the f and g functions to advance an unperturbed orbit.
 */
void fg(MDL mdl,double mu,FLOAT *x,FLOAT *v,double dt) {
	double f,g,fd,gd;			/* Gauss's f, g, fdot and gdot */
	double r,vsq;
	double u;					/* r v cos(phi) */
	double a;					/* semi-major axis */
	double e;					/* eccentricity */
	double ec,es;				/* e cos(E), e sin(E) */
	double en;					/* mean motion */
	double nf;					/* orbital frequency */
	double dec;					/* delta E */
	double dm;					/* delta mean anomoly */
	double lo = -2*M_PI;
	double up = 2*M_PI;
	double w;					/* function to zero */
	double wp;					/* first derivative */
	double wpp;					/* second derivative */
	double wppp;				/* third derivative */
	double dx,s,c,sh;
	double next;
	double converge;			/* converge criterion */
	int iter,i,j;
	const double DOUBLE_EPS = 1.2e-16;

	/* 
	 * Evaluate some orbital quantites.
	 */
	r = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
	vsq = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
	u = x[0]*v[0] + x[1]*v[1] + x[2]*v[2];
	a = 1/(2/r-vsq/mu);
	en = sqrt(mu/(a*a*a));
	ec = 1-r/a;
	es = u/(en*a*a);
	e = sqrt(ec*ec + es*es);
#if (0)
	/*
	 ** Only for GR!
	 */
	dt *= 1 - 1.5*mu/(csq*a);
#endif
	nf = en/(2*M_PI);
	j = (int)(nf*dt);
	dt -= j/nf;					/* reduce to single orbit */
	dm = en*dt - es;
	if ((es*cos(dm)+ec*sin(dm)) > 0) dec = dm + 0.85*e;
	else dec = dm - 0.85*e;
	dm = en*dt;					/* reset mean anomoly */
	converge = fabs(dm*DOUBLE_EPS);
	/*
	 * First double precision iterations for dec solution.
	 * This solves Kepler's equation in difference form:
	 * dM = dE - e cos(E_0) sin(dE) + e sin(E_0)(1 - cos(dE))
	 * Note that (1 - cos(dE)) == 2*sin(dE/2)*sin(dE/2).  The latter is
	 * probably more accurate for small dE.
	 */
	for(iter=1;iter<=32;iter++) {
		s = sin(dec);
		c = cos(dec);
		sh = sin(0.5*dec);
		w = dec - ec*s + es*2*sh*sh - dm;
		if (w > 0) up = dec;
		else lo = dec;
		wp = 1 - ec*c + es*s;
		wpp = ec*s + es*c;
		wppp = ec*c - es*s;
		dx = -w/wp;
		dx = -w/(wp + 0.5*dx*wpp);
		dx = -w/(wp + 0.5*dx*wpp + dx*dx*wppp/6);
		next = dec + dx;
		if (fabs(dx) <= converge) break;
		if (next > lo && next < up) dec = next;
		else dec = 0.5*(lo + up);
		if (dec==lo || dec==up) break;
		}
	if (iter>32) {
		char ach[256];

		sprintf(ach,"dec soln failed, dconverge: %g dx: %g dec: %g\n",
				converge, dx, dec);
		mdlDiag(mdl,ach);
		exit(1);
		}
	/*
	 * Update the orbit.
	 * See Danby, eq. 6.8.11 and 6.8.12, for expressions for f, g, fdot,
	 * and gdot.
	 * JS: changed f and gd expressions to take advantage of the half angle 
	 * formula given above.
	 */
	f = 1 - (a/r)*2*sh*sh;
	g = dt + (s-dec)/en;
	fd = -(a/(r*wp))*en*s;
	gd = 1 - 2*sh*sh/wp;
	for (i=0;i<3;i++) {
		s = f*x[i]+g*v[i];
		v[i] = fd*x[i]+gd*v[i];
		x[i] = s;
		}
	}


void pkdDrift(PKD pkd,double dDelta,FLOAT fCenter[3],int bPeriodic,int bFandG,
			  FLOAT fCentMass)
{
	PARTICLE *p;
	int i,j,n;

	p = pkd->pStore;
	n = pkdLocal(pkd);
	for (i=0;i<n;++i) {
#ifdef COLLISIONS
		if (p[i].iDriftType == KEPLER) /*DEBUG a bit ugly...*/
#else /* COLLISIONS */
		if (bFandG)
#endif /* !COLLISIONS */
			fg(pkd->mdl,fCentMass + p[i].fMass,p[i].r,p[i].v,dDelta);
		else {
			for (j=0;j<3;++j) {
				p[i].r[j] += dDelta*p[i].v[j];
				if (bPeriodic) {
					if (p[i].r[j] >= fCenter[j]+0.5*pkd->fPeriod[j])
						p[i].r[j] -= pkd->fPeriod[j];
					if (p[i].r[j] < fCenter[j]-0.5*pkd->fPeriod[j])
						p[i].r[j] += pkd->fPeriod[j];
					}
				}
			}
		}
	}

void pkdKick(PKD pkd, double dvFacOne, double dvFacTwo, double dvPredFacOne,
	     double dvPredFacTwo, double duDelta, double duPredDelta, int iGasModel,
	     double z, double duDotLimit )
{
	PARTICLE *p;
	int i,j,n;
#ifdef GASOLINE
	int bCool = 0;
	CL *cl;
	PERBARYON Y;
	double E,dt,dtPred;

	pkdClearTimer(pkd,1);
	pkdStartTimer(pkd,1);

	switch (iGasModel) {
	case GASMODEL_COOLING:
	case GASMODEL_COOLING_NONEQM:
	         bCool = 1;
	         cl = &(pkd->cl);
                 clRatesRedshift( cl, z );
		 dt = duDelta * cl->dSecUnit;
		 dtPred = duPredDelta * cl->dSecUnit;
		 break;
		 }
#endif

#ifdef DEBUG
	printf("pkdKick: %f %f %f %f %f %f\n",dvFacOne,dvFacTwo,
	       dvPredFacOne,dvPredFacTwo,duDelta,duPredDelta );
#endif
	p = pkd->pStore;
	n = pkdLocal(pkd);
	for (i=0;i<n;++i,++p) {
#ifdef GASOLINE
#ifdef DEBUG
		if (pkdIsGas(pkd,p) && (p->iOrder % 300)==0) {
		        printf("Kicking %i: %i %i %i %f %f %f %f   %f %f %f %f %f %f\n",
			       p->iOrder,TYPEQueryACTIVE(p),
                               TYPEQueryTREEACTIVE(p),p->iRung,
			       p->r[0],p->fDensity,p->u,p->uPred,
			       p->a[0],p->a[1],p->a[2],p->PdV,
			       p->PdVpres,p->PdVvisc);
		}
#endif
#endif
	    if(TYPEQueryACTIVE(p)) {
#ifdef GASOLINE
	                if(pkdIsGas(pkd, p)) {
				for (j=0;j<3;++j) {
					p->vPred[j] = p->v[j]*dvPredFacOne + 
						p->a[j]*dvPredFacTwo;
					}
				if (bCool) {
				       pkdPARTICLE2PERBARYON(&Y, p, cl->Y_H, cl->Y_He);
				       E = p->u * cl->dErgPerGmUnit;
#ifdef COOLDEBUG
				       if ((p->iOrder % 300)==0) {
					    printf("COOL: %i ",p->iOrder);
					    clIntegrateEnergyDEBUG(cl, &Y, &E, p->PdV*cl->dErgPerGmPerSecUnit, 
							   p->fDensity*cl->dComovingGmPerCcUnit, dtPred);
					    printf("\n");
				            }
				       else {
					    clIntegrateEnergy(cl, &Y, &E, p->PdV*cl->dErgPerGmPerSecUnit, 
							   p->fDensity*cl->dComovingGmPerCcUnit, dtPred);
				            }
#else				       
				       clIntegrateEnergy(cl, &Y, &E, p->PdV*cl->dErgPerGmPerSecUnit, 
							 p->fDensity*cl->dComovingGmPerCcUnit, dtPred);
#endif
				       p->uPred = E * cl->diErgPerGmUnit;

				       pkdPARTICLE2PERBARYON(&Y, p, cl->Y_H, cl->Y_He);
				       E = p->u * cl->dErgPerGmUnit;
#ifdef COOLDEBUG
				       if ((p->iOrder % 300)==0) {
					    printf("COOL: %i ",p->iOrder);
					    clIntegrateEnergyDEBUG(cl, &Y, &E, p->PdV*cl->dErgPerGmPerSecUnit, 
							   p->fDensity*cl->dComovingGmPerCcUnit, dt);
					    printf("\n");
				            }
				       else {
					    clIntegrateEnergy(cl, &Y, &E, p->PdV*cl->dErgPerGmPerSecUnit, 
							   p->fDensity*cl->dComovingGmPerCcUnit, dt);
				            }
#else				       
				       clIntegrateEnergy(cl, &Y, &E, p->PdV*cl->dErgPerGmPerSecUnit, 
							 p->fDensity*cl->dComovingGmPerCcUnit, dt);
#endif
				       pkdPERBARYON2PARTICLE(&Y, p);
				       p->u = E * cl->diErgPerGmUnit;
				       }
				else { 
 				       if (p->PdV*duPredDelta < duDotLimit*p->u) 
					      p->uPred = p->u * exp(p->PdV*duPredDelta/p->u);
				       else 
					      p->uPred = p->u + p->PdV*duPredDelta;

 				       if (p->PdV*duDelta < duDotLimit*p->u) 
					      p->u = p->u * exp(p->PdV*duDelta/p->u);
				       else 
				              p->u = p->u + p->PdV*duDelta;
				       }
				}
			
#endif
			for (j=0;j<3;++j) {
				p->v[j] = p->v[j]*dvFacOne + p->a[j]*dvFacTwo;
#ifdef SAND_PILE
				if (p->r[2] > 2*p->fSoft)
					p->v[2] -= 0.25*dvFacTwo; /* uniform 0.25 -z accel. */
				if (p->r[2] > 10)
					p->v[2] = -1; /* uniform velocity above z = 10 */
#endif /* SAND_PILE */
				}
			}
	    }

	pkdStopTimer(pkd,1);
	}


void pkdReadCheck(PKD pkd,char *pszFileName,int iVersion,int iOffset,
				  int nStart,int nLocal)
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
	lStart = iOffset+nStart*sizeof(CHKPART);
	fseek(fp,lStart,SEEK_SET);
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
#ifdef GASOLINE
		pkd->pStore[i].u = cp.u;
		pkd->pStore[i].uPred = cp.u;
#endif
#ifdef COLLISIONS
		for (j=0;j<3;++j)
			pkd->pStore[i].w[j] = cp.w[j];
		pkd->pStore[i].iColor = cp.iColor;
#endif /* COLLISIONS */
		TYPESet(&(pkd->pStore[i]),TYPE_ACTIVE);
		pkd->pStore[i].iRung = 0;
		pkd->pStore[i].fWeight = 1.0;	/* set the initial weight to 1.0 */
		pkd->pStore[i].fDensity = 0.0;
		pkd->pStore[i].fBall2 = 0.0;
		}
	fclose(fp);	
	}


void pkdWriteCheck(PKD pkd,char *pszFileName,int iOffset,int nStart)
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
	lStart = iOffset+nStart*sizeof(CHKPART);
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
#ifdef GASOLINE
		cp.u = pkd->pStore[i].u;
#endif
#ifdef COLLISIONS
		for (j=0;j<3;++j)
			cp.w[j] = pkd->pStore[i].w[j];
		cp.iColor = pkd->pStore[i].iColor;
#endif /* COLLISIONS */
		fwrite(&cp,sizeof(CHKPART),1,fp);
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
		if (pkdn[i].pUpper) {
			pkd->kdTop[i] = pkdn[i];
			if (pkdn[i].pLower >= 0) pkd->piLeaf[pkdn[i].pLower] = i;
			}
		}
	}


void pkdCalcRoot(PKD pkd,struct ilCellNewt *pcc)
{
	KDN *pkdn;
	int pj;
	double m,dx,dy,dz;
	double d2;

	/*
	 ** Initialize moments.
	 */
	pcc->xxxx = 0.0;
	pcc->xyyy = 0.0;
	pcc->xxxy = 0.0;
	pcc->yyyy = 0.0;
	pcc->xxxz = 0.0;
	pcc->yyyz = 0.0;
	pcc->xxyy = 0.0;
	pcc->xxyz = 0.0;
	pcc->xyyz = 0.0;
	pcc->xxzz = 0.0;
	pcc->xyzz = 0.0;
	pcc->xzzz = 0.0;
	pcc->yyzz = 0.0;
	pcc->yzzz = 0.0;
	pcc->zzzz = 0.0;
	pcc->xxx = 0.0;
	pcc->xyy = 0.0;
	pcc->xxy = 0.0;
	pcc->yyy = 0.0;
	pcc->xxz = 0.0;
	pcc->yyz = 0.0;
	pcc->xyz = 0.0;
	pcc->xzz = 0.0;
	pcc->yzz = 0.0;
	pcc->zzz = 0.0;
	/*
	 ** Calculate moments.
	 */
	pkdn = &pkd->kdNodes[pkd->iRoot];
	for (pj=pkdn->pLower;pj<=pkdn->pUpper;++pj) {
		m = pkd->pStore[pj].fMass;
		dx = pkd->pStore[pj].r[0] - pkdn->r[0];
		dy = pkd->pStore[pj].r[1] - pkdn->r[1];
		dz = pkd->pStore[pj].r[2] - pkdn->r[2];
#ifdef REDUCED_EWALD
		d2 = dx*dx + dy*dy + dz*dz;
#else
		d2 = 0.0;
#endif
		pcc->xxxx += m*(dx*dx*dx*dx - 6.0/7.0*d2*(dx*dx - 0.1*d2));
		pcc->xyyy += m*(dx*dy*dy*dy - 3.0/7.0*d2*dx*dy);
		pcc->xxxy += m*(dx*dx*dx*dy - 3.0/7.0*d2*dx*dy);
		pcc->yyyy += m*(dy*dy*dy*dy - 6.0/7.0*d2*(dy*dy - 0.1*d2));
		pcc->xxxz += m*(dx*dx*dx*dz - 3.0/7.0*d2*dx*dz);
		pcc->yyyz += m*(dy*dy*dy*dz - 3.0/7.0*d2*dy*dz);
		pcc->xxyy += m*(dx*dx*dy*dy - 1.0/7.0*d2*(dx*dx + dy*dy - 0.2*d2));
		pcc->xxyz += m*(dx*dx*dy*dz - 1.0/7.0*d2*dy*dz);
		pcc->xyyz += m*(dx*dy*dy*dz - 1.0/7.0*d2*dx*dz);
		pcc->xxzz += m*(dx*dx*dz*dz - 1.0/7.0*d2*(dx*dx + dz*dz - 0.2*d2));
		pcc->xyzz += m*(dx*dy*dz*dz - 1.0/7.0*d2*dx*dy);
		pcc->xzzz += m*(dx*dz*dz*dz - 3.0/7.0*d2*dx*dz);
		pcc->yyzz += m*(dy*dy*dz*dz - 1.0/7.0*d2*(dy*dy + dz*dz - 0.2*d2));
		pcc->yzzz += m*(dy*dz*dz*dz - 3.0/7.0*d2*dy*dz);
		pcc->zzzz += m*(dz*dz*dz*dz - 6.0/7.0*d2*(dz*dz - 0.1*d2));
		/*
		 ** Calculate octopole moment...
		 */
		pcc->xxx += m*(dx*dx*dx - 0.6*d2*dx);
		pcc->xyy += m*(dx*dy*dy - 0.2*d2*dx);
		pcc->xxy += m*(dx*dx*dy - 0.2*d2*dy);
		pcc->yyy += m*(dy*dy*dy - 0.6*d2*dy);
		pcc->xxz += m*(dx*dx*dz - 0.2*d2*dz);
		pcc->yyz += m*(dy*dy*dz - 0.2*d2*dz);
		pcc->xyz += m*dx*dy*dz;
		pcc->xzz += m*(dx*dz*dz - 0.2*d2*dx);
		pcc->yzz += m*(dy*dz*dz - 0.2*d2*dy);
		pcc->zzz += m*(dz*dz*dz - 0.6*d2*dz);
		}
	}


void pkdDistribRoot(PKD pkd,struct ilCellNewt *pcc)
{
	KDN *pkdn;
	double tr;

	pkd->ilcnRoot = *pcc;
	/*
	 ** Must set the quadrupole, mass and cm.
	 */
	pkdn = &pkd->kdTop[ROOT];
#ifdef REDUCED_EWALD
	tr = pkdn->mom.Qxx + pkdn->mom.Qyy + pkdn->mom.Qzz;
#else
	tr = 0.0;
#endif
	pkd->ilcnRoot.m = pkdn->fMass;
	pkd->ilcnRoot.x = pkdn->r[0];
	pkd->ilcnRoot.y = pkdn->r[1];
	pkd->ilcnRoot.z = pkdn->r[2];
	pkd->ilcnRoot.xx = pkdn->mom.Qxx - tr/3.0;
	pkd->ilcnRoot.xy = pkdn->mom.Qxy;
	pkd->ilcnRoot.xz = pkdn->mom.Qxz;
	pkd->ilcnRoot.yy = pkdn->mom.Qyy - tr/3.0;
	pkd->ilcnRoot.yz = pkdn->mom.Qyz;
	pkd->ilcnRoot.zz = pkdn->mom.Qzz - tr/3.0;
	}


double pkdMassCheck(PKD pkd) 
{
	double dMass=0.0;
	int i,iRej;

	for (i=0;i<pkdLocal(pkd);++i) {
		dMass += pkd->pStore[i].fMass;
		}
	iRej = pkdFreeStore(pkd) - pkd->nRejects;
	for (i=0;i<pkd->nRejects;++i) {
		dMass += pkd->pStore[iRej+i].fMass;
		}
	return(dMass);
	}

void
pkdSetRung(PKD pkd, int iRung)
{
    int i;
    
    for(i = 0; i < pkdLocal(pkd); ++i) {
	pkd->pStore[i].iRung = iRung;
	}
    }

int
pkdActiveRung(PKD pkd, int iRung, int bGreater)
{
    int i;
    int nActive;
    char out[128];
    
    nActive = 0;
    for(i = 0; i < pkdLocal(pkd); ++i) {
	if(pkd->pStore[i].iRung == iRung
	   || (bGreater && pkd->pStore[i].iRung > iRung)) {
	    TYPESet(&(pkd->pStore[i]),TYPE_ACTIVE);
	    ++nActive;
	    }
	else
	    TYPEReset(&(pkd->pStore[i]),TYPE_ACTIVE);
	}
    sprintf(out, "nActive: %d\n", nActive);
    mdlDiag(pkd->mdl, out);
    pkd->nActive = nActive;
    return nActive;
    }

int
pkdCurrRung(PKD pkd, int iRung)
{
    int i;
    int iCurrent;
    
    iCurrent = 0;
    for(i = 0; i < pkdLocal(pkd); ++i) {
	if(pkd->pStore[i].iRung == iRung) {
	    iCurrent = 1;
	    break;
	    }
	}
    return iCurrent;
    }

void
pkdDensityStep(PKD pkd, double dEta, double dRhoFac)
{
    int i;
    double dT;
    
    for(i = 0; i < pkdLocal(pkd); ++i) {
		if(TYPEQueryACTIVE(&(pkd->pStore[i]))) {
			dT = dEta/sqrt(pkd->pStore[i].fDensity*dRhoFac);
			if(dT < pkd->pStore[i].dt)
				pkd->pStore[i].dt = dT;
			}
		}
    }


void
pkdAccelStep(PKD pkd, double dEta, double dVelFac, double dAccFac, int bDoGravity,
             int bEpsVel, int bSqrtPhi)
{
    int i;
    double vel;
    double acc;
    int j;
    double dT;
    
    for(i = 0; i < pkdLocal(pkd); ++i) {
		if(TYPEQueryACTIVE(&(pkd->pStore[i]))) {
			vel = 0;
			acc = 0;
			for(j = 0; j < 3; j++) {
				vel += pkd->pStore[i].v[j]*pkd->pStore[i].v[j];
				acc += pkd->pStore[i].a[j]*pkd->pStore[i].a[j];
				}
			assert(vel >= 0);
			vel = sqrt(vel)*dVelFac;
			assert(acc > 0);
			acc = sqrt(acc)*dAccFac;
			dT = FLOAT_MAXVAL;
			if(bEpsVel)
#ifdef GASOLINE
				if (pkdIsGas(pkd, &(pkd->pStore[i])) &&
					pkd->pStore[i].fBall2*0.25 < pkd->pStore[i].fSoft*pkd->pStore[i].fSoft) 
					dT = dEta*sqrt(sqrt(0.25*pkd->pStore[i].fBall2)/acc);
                else
					dT = dEta*sqrt(pkd->pStore[i].fSoft/acc);
#else
			dT = dEta*sqrt(pkd->pStore[i].fSoft/acc);
#endif
#ifdef DEBUG
			if ((pkd->pStore[i].iOrder % 300)==0 || dT<1e-6) {
				printf("dt_a %i: %i %i %f %g %g %g  %f %f %f %f   dt_a %g\n",
				       pkd->pStore[i].iOrder,TYPEQueryACTIVE(&(pkd->pStore[i])),
					   TYPEQueryTREEACTIVE(&(pkd->pStore[i])),
				       sqrt(0.25*pkd->pStore[i].fBall2),pkd->pStore[i].PdV,pkd->pStore[i].u,pkd->pStore[i].PoverRho2,pkd->pStore[i].a[0],
				       pkd->pStore[i].a[1],pkd->pStore[i].a[2],pkd->pStore[i].fSoft,dT);
				}
#endif
			if(bSqrtPhi) {
				double dtemp =
					dEta*3.5*sqrt(dAccFac*fabs(pkd->pStore[i].fPot))/acc;
				if(dtemp < dT)
					dT = dtemp;
				}
			if(dT < pkd->pStore[i].dt)
				pkd->pStore[i].dt = dT;
			}
		}
    }

int
pkdDtToRung(PKD pkd, int iRung, double dDelta, int iMaxRung, int bAll)
{
    int i;
    int iMaxRungOut;
    int iTempRung;
    int iSteps;
    
    iMaxRungOut = 0;
    for(i = 0; i < pkdLocal(pkd); ++i) {
		if(pkd->pStore[i].iRung >= iRung) {
			assert(TYPEQueryACTIVE(&(pkd->pStore[i])));
			if(bAll) {          /* Assign all rungs at iRung and above */
				iSteps = dDelta/pkd->pStore[i].dt;
				/* insure that integer boundary goes
				   to the lower rung. */
				if(fmod(dDelta,pkd->pStore[i].dt) == 0.0)
				    iSteps--;
				iTempRung = iRung;
				if(iSteps < 0)
				    iSteps = 0;
				while(iSteps) {
					++iTempRung;
					iSteps >>= 1;
					}
				if(iTempRung >= iMaxRung)
					iTempRung = iMaxRung-1;
				pkd->pStore[i].iRung = iTempRung;
				}
			else {
				if(dDelta <= pkd->pStore[i].dt) {
					pkd->pStore[i].iRung = iRung;
					}
				else {
					pkd->pStore[i].iRung = iRung+1;
					}
                }
			}
		if(pkd->pStore[i].iRung > iMaxRungOut)
			iMaxRungOut = pkd->pStore[i].iRung;
		}
    return iMaxRungOut;
    }

void
pkdInitDt(PKD pkd, double dDelta)
{
    int i;
    
    for(i = 0; i < pkdLocal(pkd); ++i) {
		if(TYPEQueryACTIVE(&(pkd->pStore[i])))
			pkd->pStore[i].dt = dDelta;
		}
    }

    
int pkdRungParticles(PKD pkd,int iRung)
{
	int i,n;

	n = 0;
	for (i=0;i<pkdLocal(pkd);++i) {
		if (pkd->pStore[i].iRung == iRung) ++n;
		}
	return n;
	}

void
pkdDeleteParticle(PKD pkd, int i)
{
    pkd->pStore[i].iOrder = -2 - pkd->pStore[i].iOrder;
	TYPEReset(&(pkd->pStore[i]),TYPE_ACTIVE);
    }

void
pkdNewParticle(PKD pkd, PARTICLE p)
{
    assert(pkd->nLocal < pkd->nStore);
    pkd->pStore[pkd->nLocal] = p;
    pkd->pStore[pkd->nLocal].iOrder = -1;
    pkd->nLocal++;
    }

void
pkdColNParts(PKD pkd, int *pnNew, int *nDeltaGas, int *nDeltaDark,
	     int *nDeltaStar)
{
    int pi, pj;
    int nNew;
    int ndGas;
    int ndDark;
    int ndStar;
    int newnLocal;
    PARTICLE *p;
    
    nNew = 0;
    ndGas = 0;
    ndDark = 0;
    ndStar = 0;
    newnLocal = pkdLocal(pkd);
    for(pi = 0, pj = 0; pi < pkdLocal(pkd); pi++) {
	if(pj < pi)
	    pkd->pStore[pj] = pkd->pStore[pi];
	p = &pkd->pStore[pi];
	if(p->iOrder == -1) {
	    ++pj;
	    ++nNew;
#ifdef GASOLINE
	    ++ndStar;
#else
	    ++ndDark;
#endif
	    if(TYPEQueryACTIVE(p))
		++pkd->nActive;
	    continue;
	    }
	else if(p->iOrder < -1){
	    --newnLocal;
	    p->iOrder = -2 - p->iOrder;
	    if(pkdIsGas(pkd, p))
		--ndGas;
	    else if(pkdIsDark(pkd, p))
		--ndDark;
	    else if(pkdIsStar(pkd, p))
		--ndStar;
	    else
		assert(0);
	    if(TYPEQueryACTIVE(p))
		--pkd->nActive;
	    }
	else {
	    ++pj;
	    }
	}

    *pnNew = nNew;
    *nDeltaGas = ndGas;
    *nDeltaDark = ndDark;
    *nDeltaStar = ndStar;
    pkd->nLocal = newnLocal;
    }

void
pkdNewOrder(PKD pkd, int nStart)
{
    int pi;
    
    for(pi = 0; pi < pkdLocal(pkd); pi++) {
	if(pkd->pStore[pi].iOrder == -1) {
	    pkd->pStore[pi].iOrder = nStart++;
	    }
	}
    }

void
pkdSetNParts(PKD pkd, int nGas, int nDark, int nStar, int nMaxOrderGas,
	     int nMaxOrderDark)
{
    pkd->nGas = nGas;
    pkd->nDark = nDark;
    pkd->nStar = nStar;
    pkd->nMaxOrderGas = nMaxOrderGas;
    pkd->nMaxOrderDark = nMaxOrderDark;
}

void pkdCoolVelocity(PKD pkd,int nSuperCool,double dCoolFac,
					 double dCoolDens,double dCoolMaxDens)
{
	PARTICLE *p;
	int i,j,n;

	p = pkd->pStore;
	n = pkdLocal(pkd);
	for (i=0;i<n;++i) {
		if (p[i].iOrder < nSuperCool && p[i].fDensity >= dCoolDens &&
			p[i].fDensity < dCoolMaxDens) {
			for (j=0;j<3;++j) {
#ifdef SUPERCOOL
				p[i].v[j] = p[i].vMean[j] + (p[i].v[j]-p[i].vMean[j])*dCoolFac;
#else
				p[i].v[j] *= dCoolFac;
#endif
				}
			}
		}
	}

int pkdResetTouchRung(PKD pkd, unsigned int iTestMask, unsigned int iSetMask)
{
    PARTICLE *p;
    int i, nActive = 0;

    for(i=0;i<pkdLocal(pkd);++i) { 
                p = &pkd->pStore[i];
		if (TYPETest( p, iTestMask )) {
		        if (p->iTouchRung & iSetMask) {
			       p->iTouchRung &= (~iSetMask);
			       TYPESet( p, TYPE_SMOOTHACTIVE );
			       nActive++;
			       }
		        }
		}
    return nActive;
    }

int pkdActiveExactType(PKD pkd, unsigned int iFilterMask, unsigned int iTestMask, unsigned int iSetMask)
{
    PARTICLE *p;
    int i, nActive = 0;

    for(i=0;i<pkdLocal(pkd);++i) { 
                p = &pkd->pStore[i];
                /* DEBUG: Paranoia check */
                assert( TYPETest( p, TYPE_ALL ) );
		if (TYPEFilter( p, iFilterMask, iTestMask )) {
			TYPESet( p, iSetMask );
                        nActive++;
			}
		else {
			TYPEReset( p, iSetMask );
			}
		}
    if (iSetMask & TYPE_ACTIVE) pkd->nActive = nActive;
    if (iSetMask & TYPE_TREEACTIVE) pkd->nTreeActive = nActive;
    if (iSetMask & TYPE_SMOOTHACTIVE) pkd->nSmoothActive = nActive;
    return nActive;
    }

int pkdSetType(PKD pkd, unsigned int iTestMask, unsigned int iSetMask)
{
    PARTICLE *p;
    int i, nActive = 0;

    for(i=0;i<pkdLocal(pkd);++i) {
                p = &pkd->pStore[i];
                /* DEBUG: Paranoia check */
                assert( TYPETest( p, TYPE_ALL ) );
		if (TYPETest( p, iTestMask )) {
			TYPESet( p, iSetMask );
                        nActive++;
			}
		}
    /*
      Need to fix this up:
    if (iSetMask & TYPE_ACTIVE) pkd->nActive = nActive;
    if (iSetMask & TYPE_TREEACTIVE) pkd->nTreeActive = nActive;
    if (iSetMask & TYPE_SMOOTHACTIVE) pkd->nSmoothActive = nActive;
    */
    return nActive;
    }

int pkdResetType(PKD pkd, unsigned int iTestMask, unsigned int iSetMask)
{
    PARTICLE *p;
    int i, nActive = 0;

    for(i=0;i<pkdLocal(pkd);++i) {
                p = &pkd->pStore[i];
                /* DEBUG: Paranoia check */
                assert( TYPETest( p, TYPE_ALL ) );
		if (TYPETest( p, iTestMask )) {
			TYPEReset( p, iSetMask );
                        nActive++;
			}
		}
    /*
    if (iSetMask & TYPE_ACTIVE) pkd->nActive = nActive;
    if (iSetMask & TYPE_TREEACTIVE) pkd->nTreeActive = nActive;
    if (iSetMask & TYPE_SMOOTHACTIVE) pkd->nSmoothActive = nActive;
    */
    return nActive;
    }

int pkdCountType(PKD pkd, unsigned int iFilterMask, unsigned int iTestMask)
{
    PARTICLE *p;
    int i, nActive = 0;

    for(i=0;i<pkdLocal(pkd);++i) {
                p = &pkd->pStore[i];
		if (TYPEFilter( p, iFilterMask, iTestMask )) {
                        nActive++;
			}
		}
    {
       char debug[100];
       sprintf(debug, "Filter %d:%d, Counted: %d\n",iFilterMask,iTestMask,nActive);
       mdlDiag(pkd->mdl,debug);
       }
    return nActive;
    }

int pkdActiveType(PKD pkd, unsigned int iTestMask, unsigned int iSetMask)
{
    PARTICLE *p;
    int i, nActive = 0;

    for(i=0;i<pkdLocal(pkd);++i) {
                p = &pkd->pStore[i];
                /* DEBUG: Paranoia check */
                assert( TYPETest( p, TYPE_ALL ) );
		if (TYPETest( p, iTestMask )) {
			TYPESet( p, iSetMask );
                        nActive++;
			}
		else {
			TYPEReset( p, iSetMask );
			}
		}
    if (iSetMask & TYPE_ACTIVE      ) pkd->nActive       = nActive;
    if (iSetMask & TYPE_TREEACTIVE  ) pkd->nTreeActive   = nActive;
    if (iSetMask & TYPE_SMOOTHACTIVE) pkd->nSmoothActive = nActive;
    return nActive;
    }

int
pkdActiveMaskRung(PKD pkd, unsigned iSetMask, int iRung, int bGreater)
{
    PARTICLE *p;
    int i;
    int nActive;
    char out[128];
    
    nActive = 0;
    for(i = 0; i < pkdLocal(pkd); ++i) {
        p = &pkd->pStore[i];
	if(p->iRung == iRung
	   || (bGreater && p->iRung > iRung)) {
	    TYPESet( p, iSetMask );
	    ++nActive;
	    }
	else
	    TYPEReset( p, iSetMask );
	}
    sprintf(out, "nActive: %d\n", nActive);
    mdlDiag(pkd->mdl, out);

    if ( iSetMask & TYPE_ACTIVE      ) pkd->nActive       = nActive;
    if ( iSetMask & TYPE_TREEACTIVE  ) pkd->nTreeActive   = nActive;
    if ( iSetMask & TYPE_SMOOTHACTIVE) pkd->nSmoothActive = nActive;
    return nActive;
    }

int
pkdActiveTypeRung(PKD pkd, unsigned iTestMask, unsigned iSetMask, int iRung, int bGreater)
{
    PARTICLE *p;
    int i;
    int nActive;
    char out[128];
    
    nActive = 0;
    for(i = 0; i < pkdLocal(pkd); ++i) {
        p = &pkd->pStore[i];
        /* DEBUG: Paranoia check */
        assert( TYPETest( p, TYPE_ALL ) );
	if( TYPETest( p, iTestMask ) && 
           ( p->iRung == iRung || (bGreater && p->iRung > iRung) )) {
	    TYPESet( p, iSetMask );
	    ++nActive;
	    }
	else
	    TYPEReset( p, iSetMask );
	}
    sprintf(out, "nActive: %d\n", nActive);
    mdlDiag(pkd->mdl, out);

    if ( iSetMask & TYPE_ACTIVE      ) pkd->nActive       = nActive;
    if ( iSetMask & TYPE_TREEACTIVE  ) pkd->nTreeActive   = nActive;
    if ( iSetMask & TYPE_SMOOTHACTIVE) pkd->nSmoothActive = nActive;
    return nActive;
    }

void pkdSetParticleTypes(PKD pkd, int nSuperCool) {
    PARTICLE *p;
    int i, iSetMask;
    
    for(i=0;i<pkdLocal(pkd);++i) {
                p = &pkd->pStore[i];
                iSetMask = 0;
                if (pkdIsGas (pkd,p)) iSetMask |= TYPE_GAS;
                if (pkdIsDark(pkd,p)) iSetMask |= TYPE_DARK;
                if (pkdIsStar(pkd,p)) iSetMask |= TYPE_STAR;
                if (p->iOrder < nSuperCool) iSetMask |= TYPE_SUPERCOOL;

                TYPESet( p, iSetMask );
                }
    }

void pkdGrowMass(PKD pkd,int nGrowMass, double dDeltaM)
{
    int i;
    
#ifdef GROWMASS
    for(i=0;i<pkdLocal(pkd);++i) {
		if (pkd->pStore[i].iOrder < nGrowMass) {
			pkd->pStore[i].fMass += dDeltaM;
			}
		}
#endif
    }

void pkdInitAccel(PKD pkd)
{
    int i,j;
    
    for(i=0;i<pkdLocal(pkd);++i) {
		if (TYPEQueryACTIVE(&(pkd->pStore[i]))) {
			for (j=0;j<3;++j) {
				pkd->pStore[i].a[j] = 0;
				}
			}
		}
    }

int pkdIsGas(PKD pkd,PARTICLE *p) {
	if (p->iOrder <= pkd->nMaxOrderGas) return 1;
	else return 0;
	}

int pkdIsDark(PKD pkd,PARTICLE *p) {
	if (p->iOrder > pkd->nMaxOrderGas && p->iOrder <= pkd->nMaxOrderDark)
	    return 1;
	else return 0;
	}

int pkdIsStar(PKD pkd,PARTICLE *p) {
	if (p->iOrder > pkd->nMaxOrderDark) return 1;
	else return 0;
	}

#ifdef GASOLINE

/* Note: Uses uPred */
void pkdAdiabaticGasPressure(PKD pkd, double gammam1, double gamma)
{
    PARTICLE *p;
    double PoverRho;
    int i;

    p = pkd->pStore;
    for(i=0;i<pkdLocal(pkd);++i,++p) {
                if (TYPEQueryTREEACTIVE(p)) {
		        PoverRho = gammam1*p->uPred;
			p->PoverRho2 = PoverRho/p->fDensity;
   			p->c = sqrt(gamma*PoverRho);
		        }
#ifdef DEBUG
		if (pkdIsGas(pkd,p) && (p->iOrder % 1000)==0) {
		        printf("Pressure %i: %i %i %f %f %f  %f %f %f %f %f\n",
			       p->iOrder,TYPEQueryACTIVE(p),TYPEQueryTREEACTIVE(p),
			       p->r[0],p->r[1],p->r[2],sqrt(0.25*p->fBall2),p->fDensity,p->uPred,
			       p->PoverRho2*p->fDensity*p->fDensity,p->c);
		        } 
#endif            
                }
    }

void pkdInitEnergy(PKD pkd, double dTuFac, double z)
{
    PARTICLE *p;
    int i;
    CL *cl;
    PERBARYON Y;
    RATE r;
    double T;

    cl = &(pkd->cl);
    clRatesRedshift( cl, z );

    p = pkd->pStore;
    for(i=0;i<pkdLocal(pkd);++i,++p) {
                if (TYPEQueryTREEACTIVE(p) && pkdIsGas(pkd,p)) {
		        T = p->u / dTuFac;
			pkdPARTICLE2PERBARYON(&Y, p, cl->Y_H, cl->Y_He);
			clRates( cl, &r, T );
			clAbunds( cl, &Y, &r, p->fDensity * cl->dComovingGmPerCcUnit );
			pkdPERBARYON2PARTICLE(&Y, p);
                        p->u = clThermalEnergy( Y.Total, T ) * cl->diErgPerGmUnit;
			p->uPred = p->u;
#ifdef DEBUG
			if ((p->iOrder % 1000)==0) {
			         printf("InitEnergy %i: %f %g   %f %f %f %g\n",
				        p->iOrder,T,p->u * cl->dErgPerGmUnit,
					Y.HI,Y.HeI,Y.HeII,p->fDensity * cl->dComovingGmPerCcUnit);
		                 } 
#endif            
		        }
                }
    }

#ifdef GLASS
/* Currently wired to have no more than two regions with different
   Pressures (densities) split by x=0 with a linear connection */
void pkdGlassGasPressure(PKD pkd, void *vin)
{
    PARTICLE *p;
    double PoverRho,xx,nsp=2.5;
    int i;
    struct inGetGasPressure *in = vin;

    for(i=0;i<pkdLocal(pkd);++i) {
                p = &pkd->pStore[i];
                if (TYPEQueryTREEACTIVE(p)) {
    		        if (p->r[0] < -nsp*in->dGlassxL) {
			  if (p->r[0] > in->dxBoundL + nsp*in->dGlassxL)
			       PoverRho=in->dGlassPoverRhoL;  
			  else {
			       xx = ( p->r[0] - in->dxBoundL + nsp*in->dGlassxR )
				 / ( nsp*in->dGlassxL+nsp*in->dGlassxR );
			       xx = xx*xx* ( -2*xx + 3 );
		               PoverRho = in->dGlassPoverRhoR + xx
				 *(in->dGlassPoverRhoL - in->dGlassPoverRhoR);
              		       }
			  }
                        else if (p->r[0] > nsp*in->dGlassxR) {
			  if (p->r[0] < in->dxBoundR - nsp*in->dGlassxR)
			       PoverRho=in->dGlassPoverRhoR;  
			  else {
			       xx = ( p->r[0] - in->dxBoundR + nsp*in->dGlassxR )
				 / ( nsp*in->dGlassxL+nsp*in->dGlassxR );
			       xx = xx*xx* ( -2*xx + 3 );
		               PoverRho = in->dGlassPoverRhoR + xx
				 *(in->dGlassPoverRhoL - in->dGlassPoverRhoR);
			       }
			  }
			else {
			       xx = ( p->r[0] + nsp*in->dGlassxL )
				 / ( nsp*in->dGlassxL+nsp*in->dGlassxR );
			       xx = xx*xx* ( -2*xx + 3 );
		               PoverRho = in->dGlassPoverRhoL + xx
				 *(in->dGlassPoverRhoR - in->dGlassPoverRhoL);
			       }

			p->u = PoverRho;
			p->uPred = PoverRho;
			p->PoverRho2 = PoverRho/p->fDensity;
   			p->c = sqrt(in->gamma*PoverRho);
		        }
#ifdef DEBUG
		if (pkdIsGas(pkd,p) && (p->iOrder % 1000)==0) 
		        printf("Glass P %i: %i %i %f %f %f  %f %f %f %f %f\n",
			       p->iOrder,TYPEQueryACTIVE(p),TYPEQueryTREEACTIVE(p),
			       p->r[0],p->r[1],p->r[2],sqrt(0.25*p->fBall2),p->fDensity,p->uPred,p->PoverRho2*p->fDensity*p->fDensity,p->c);
#endif
                }
    }

#endif

void pkdKickVpred(PKD pkd, double dvFacOne, double dvFacTwo, double duDelta,int iGasModel,
		  double z, double duDotLimit)
{
	PARTICLE *p;
	int i,j,n;
#ifdef GASOLINE
	int bCool = 0;
	CL *cl;
	PERBARYON Y;
	double E,dt;

	pkdClearTimer(pkd,1);
	pkdStartTimer(pkd,1);

	switch (iGasModel) {
	case GASMODEL_COOLING:
	case GASMODEL_COOLING_NONEQM:
	         bCool = 1;
		 cl = &(pkd->cl);
		 clRatesRedshift( cl, z );
		 dt = duDelta * cl->dSecUnit;
		 break;
		 }
#endif

#ifdef DEBUG
	printf("pkdKickVpred: %f %f %f\n",dvFacOne,dvFacTwo,duDelta);
#endif
        p = pkd->pStore;
	n = pkdLocal(pkd);
	for (i=0;i<n;++i,++p) {
#ifdef DEBUG
		if (pkdIsGas(pkd,p) && (p->iOrder % 3000)==0) {
		        printf("Vpreding %i: %i %i %f %f %f %f   %f %f %f %f\n",
			       p->iOrder,TYPEQueryACTIVE(p),
                               TYPEQueryTREEACTIVE(p),
			       sqrt(0.25*p->fBall2),p->fDensity,p->u,p->uPred,
			       p->a[0],p->a[1],p->a[2],p->PdV);
		        }
#endif
#ifdef GASOLINE
 	        if (pkdIsGas(pkd, p)) {
			for (j=0;j<3;++j) {
				p->vPred[j] = p->vPred[j]*dvFacOne + p->a[j]*dvFacTwo;
				}
			if (bCool) {
                                pkdPARTICLE2PERBARYON(&Y, p, cl->Y_H, cl->Y_He);
				E = p->uPred * cl->dErgPerGmUnit;
#ifdef COOLDEBUG
				if ((p->iOrder % 300)==0) {
					printf("COOL: %i ",p->iOrder);
					clIntegrateEnergyDEBUG(cl, &Y, &E, p->PdV*cl->dErgPerGmPerSecUnit, 
							       p->fDensity*cl->dComovingGmPerCcUnit, dt);
					printf("\n");
				        }
				else {
					clIntegrateEnergy(cl, &Y, &E, p->PdV*cl->dErgPerGmPerSecUnit, 
							  p->fDensity*cl->dComovingGmPerCcUnit, dt);
				        }
#else				       
				clIntegrateEnergy(cl, &Y, &E, p->PdV*cl->dErgPerGmPerSecUnit, 
						  p->fDensity*cl->dComovingGmPerCcUnit, dt);
#endif
				p->uPred = E * cl->diErgPerGmUnit;
			        }
			else {
			        if (p->PdV*duDelta < duDotLimit*p->uPred) 
				         p->uPred = p->uPred * exp(p->PdV*duDelta/p->uPred);
				else 
			                 p->uPred = p->uPred + p->PdV*duDelta;
			        }
		        }
#endif
	        }

        pkdStopTimer(pkd,1);
	}

void pkdKickRhopred(PKD pkd, double dHubbFac, double dDelta)
{
	PARTICLE *p;
	int i,j,n;

#ifdef DEBUG
	printf("pkdKickRhopred: %g %g\n",dHubbFac,dDelta);
#endif
	p = pkd->pStore;
	n = pkdLocal(pkd);
	for (i=0;i<n;++i) {
#ifdef DEBUG
		if (pkdIsGas(pkd,(p+i)) && ((p+i)->iOrder % 3000)==0) {
		        printf("Rhopreding %i: %i %i %f %f %f %f   %f %f %f %f\n",
			       (p+i)->iOrder,TYPEQueryACTIVE(p+i),
                               TYPEQueryTREEACTIVE(p+i),
			       sqrt(0.25*(p+i)->fBall2),(p+i)->fDensity,(p+i)->u,(p+i)->uPred,
			       (p+i)->a[0],(p+i)->a[1],(p+i)->a[2],(p+i)->PdV);
		        }
#endif
 	        if(TYPEFilter( &p[i], TYPE_GAS|TYPE_ACTIVE, TYPE_GAS )) {
			p[i].fDensity = p[i].fDensity*( 1 + dDelta*( dHubbFac - p[i].divv ));
			}
	        }
	}

int pkdSphCurrRung(PKD pkd, int iRung, int bGreater)
{
    int i;
    int iCurrent;
    
    iCurrent = 0;
    for(i = 0; i < pkdLocal(pkd); ++i) {
		if(pkdIsGas(pkd, &pkd->pStore[i]) && (pkd->pStore[i].iRung == iRung ||
                    (pkd->pStore[i].iRung > iRung && bGreater))) {
			iCurrent = 1;
			break;
			}
		}
    return iCurrent;
    }

void
pkdSphStep(PKD pkd, double dCosmoFac, double dEtaCourant, double dEtauDot)
{
    int i;
    PARTICLE *p;    
    double dT,dTu;

    for(i = 0; i < pkdLocal(pkd); ++i) {
        p = &pkd->pStore[i];
        if(pkdIsGas(pkd, p) && TYPEQueryACTIVE(p)) {
	       /*
		* Courant condition goes here.
		*/
  	       if (p->mumax>0.0) 
	                dT = dEtaCourant*dCosmoFac*( sqrt(0.25*p->fBall2)/(p->c + 0.6*(p->c + 2*p->mumax)) );
	       else
                        dT = dEtaCourant*dCosmoFac*( sqrt(0.25*p->fBall2)/(1.6*p->c) );

	       if (dEtauDot > 0.0 && p->PdV != 0.0) {
		 dTu = dEtauDot*p->u/fabs(p->PdV);
		 if (dTu < dT) 
		        dT = dTu;
	       }
#ifdef DEBUG	      
	       if ((p->iOrder % 300) == 0 || dT<1e-6) 
		        printf("dt_C %i: %i %i %i %f %f %f %f %f %f   dT_c %g\n",
			       p->iOrder,TYPEQueryACTIVE(p),TYPEQueryTREEACTIVE(p),p->iRung,
			       sqrt(0.25*p->fBall2),p->fDensity,p->u,p->uPred,p->c,p->mumax,
			       dT);
#endif
	       if(dT < p->dt)
			p->dt = dT;
	       }
	}
    }

void 
pkdSphViscosityLimiter(PKD pkd, int bOn)
{
    int i;
    PARTICLE *p;    
    double frog;

    if (bOn) {
        for(i = 0; i < pkdLocal(pkd); ++i) {
               p = &pkd->pStore[i];
               if(pkdIsGas(pkd, p)) {
#ifdef DEBUG
	                if ((p->iOrder % 300) == 0) 
		               printf("Particle %i: %i %i %i %f %f %f %f %f %f   Curl\n",
				 p->iOrder,TYPEQueryACTIVE(p),TYPEQueryTREEACTIVE(p),p->iRung,
				 sqrt(0.25*p->fBall2),p->fDensity,p->divv,p->curlv[0],p->curlv[1],p->curlv[2]
				 );
#endif		
    	                if (p->divv!=0.0) {         	 
                               p->BalsaraSwitch = fabs(p->divv)/(fabs(p->divv)+sqrt(p->curlv[0]*p->curlv[0]+
			         p->curlv[1]*p->curlv[1]+p->curlv[2]*p->curlv[2]));
                               }
                        else { 
                               p->BalsaraSwitch = 0.0;
                               }
	                }
	       }
        }
    else {
        for(i = 0; i < pkdLocal(pkd); ++i) {
               p = &pkd->pStore[i];
               if(pkdIsGas(pkd, p)) {
                        p->BalsaraSwitch = 1.0;
                        }
	       }
        }
    }

#endif

#ifdef GLASS

void
pkdRandomVelocities(PKD pkd, double dMaxVL, double dMaxVR)
{
    int i,j;
    PARTICLE *p;  
    double v;

    for(i = 0; i < pkdLocal(pkd); ++i) {
        p = &pkd->pStore[i];
        if (p->r[0]<0.0) {
                for (j=0;j<3;j++) {
                         p->v[j]+=(v=((random()%20001)/10000.0-1.0)*dMaxVL);
                         p->vPred[j]+=v;
                         }
                }
        else {
                for (j=0;j<3;j++) {
                         p->v[j]+=(v=((random()%20001)/10000.0-1.0)*dMaxVR);
                         p->vPred[j]+=v;
                         }
                }
	}
    }

#endif

#ifdef COLLISIONS

int
pkdNumRejects(PKD pkd)
{
	int i,nRej = 0;

	for (i=0;i<pkd->nLocal;++i)
		if (REJECT(&pkd->pStore[i])) ++nRej;

	return nRej;
	}

void
pkdReadSS(PKD pkd,char *pszFileName,int nStart,int nLocal)
{
	FILE *fp;
	XDR xdrs;
	PARTICLE *p;
	double dDum;
	int i,j,iDum;

	pkd->nLocal = nLocal;
	pkd->nActive = nLocal;
	/*
	 ** General initialization.
	 */
	for (i=0;i<nLocal;++i) {
		p = &pkd->pStore[i];
		TYPEClear(p);
                TYPESet(p,TYPE_ACTIVE);
		p->iTouchRung = 0;
		p->iRung = 0;
		p->fWeight = 1.0;
		p->fDensity = 0.0;
		p->fBall2 = 0.0;
		p->fPot = 0.0;
		p->iDriftType = NORMAL; /*DEBUG must initialize b/c pkdDrift()...*/
		p->dTEnc = 0.0; /*DEBUG must be initialized >= 0 (reject ==> -1) */
#ifdef SAND_PILE
		p->bStuck = 0;
#endif /* SAND_PILE */
		}
	/*
	 ** Seek past the header and up to nStart.
	 */
	fp = fopen(pszFileName,"r");
	assert(fp != NULL);
	/*
	 ** Seek to right place in file.
	 */
	fseek(fp,SS_HEAD_SIZE + nStart*SS_DATA_SIZE,SEEK_SET);
	/*
	 ** Read Stuff!
	 */
	xdrstdio_create(&xdrs,fp,XDR_DECODE);
	for (i=0;i<nLocal;++i) {
		p = &pkd->pStore[i];
		p->iOrder = nStart + i;
		if (!pkdIsDark(pkd,p)) assert(0);
		xdr_double(&xdrs,&dDum); p->fMass = dDum; /* SS format always double */
		xdr_double(&xdrs,&dDum); p->fSoft = 0.5*dDum;
		for (j=0;j<3;++j)
			{xdr_double(&xdrs,&dDum); p->r[j] = dDum;}
		for (j=0;j<3;++j)
			{xdr_double(&xdrs,&dDum); p->v[j] = dDum;}
		for (j=0;j<3;++j)
			{xdr_double(&xdrs,&dDum); p->w[j] = dDum;}
		xdr_int(&xdrs,&p->iColor);
		xdr_int(&xdrs,&iDum);
		}
	xdr_destroy(&xdrs);
	fclose(fp);
	}

void
pkdWriteSS(PKD pkd,char *pszFileName,int nStart)
{
	FILE *fp;
	XDR xdrs;
	PARTICLE *p;
	double dDum;
	int i,j,iDum,nout;

	/*
	 ** Seek past the header and up to nStart.
	 */
	fp = fopen(pszFileName,"r+");
	assert(fp != NULL);
	fseek(fp,SS_HEAD_SIZE + nStart*SS_DATA_SIZE,SEEK_SET);
	/* 
	 ** Write Stuff!
	 */
	xdrstdio_create(&xdrs,fp,XDR_ENCODE);
	for (i=0;i<pkdLocal(pkd);++i) {
		p = &pkd->pStore[i];
		if (!pkdIsDark(pkd,p)) assert(0);
		dDum = p->fMass; xdr_double(&xdrs,&dDum); /* SS format always double */
		dDum = 2*p->fSoft; xdr_double(&xdrs,&dDum);
		for (j=0;j<3;++j)
			{dDum = p->r[j]; xdr_double(&xdrs,&dDum);}
		for (j=0;j<3;++j)
			{dDum = p->v[j]; xdr_double(&xdrs,&dDum);}
		for (j=0;j<3;++j)
			{dDum = p->w[j]; xdr_double(&xdrs,&dDum);}
		xdr_int(&xdrs,&p->iColor);
		iDum = -1; xdr_int(&xdrs,&iDum);
		}
	xdr_destroy(&xdrs);
	nout = fclose(fp);
	assert(nout == 0);
	}

void
pkdCalcHill(PKD pkd,double dCentMass)
{
	/*
	 ** Note: normally the reduced Hill sphere contains the SUM of the masses
	 ** of the two interacting bodies, relative to the central mass. Here only
	 ** one mass is used. This is equivalent to assuming the mass of the third
	 ** body is negligible. Alternatively, if the masses are comparable, the
	 ** factor HILL_SCALE could be multiplied by 2^1/3.
	 */

	const double OneThird = 1.0/3.0;

	PARTICLE *p;
	int i;

	for (i=0;i<pkdLocal(pkd);i++) {
		p = &pkd->pStore[i];
		if (dCentMass == 0)
			p->fRedHill = HUGE_VAL;
		else
			p->fRedHill = HILL_SCALE*pow(OneThird*p->fMass/dCentMass,OneThird);
		}
	}

void
pkdHillStep(PKD pkd,double dEta)
{
	/*DEBUG assumes heliocentric coords...still good approx otherwise*/

	PARTICLE *p;
	double r2,a2,dt;
	int i,j;

	for (i=0;i<pkdLocal(pkd);i++) {
		p = &pkd->pStore[i];
		if (TYPEQueryACTIVE(p)) {
			r2 = a2 = 0;
			/*DEBUG could use Manhattan metric here...*/
			for (j=0;j<3;j++) {
				r2 += p->r[j]*p->r[j];
				a2 += p->a[j]*p->a[j];
				}
			dt = dEta*sqrt(p->fRedHill*sqrt(r2/a2));
			if (dt < p->dt)
				p->dt = dt;
			}
		}
	}

void
pkdFindEncounter(PKD pkd,double *dNext)
{
	PARTICLE *p;
	int i;

	for (i=0;i<pkdLocal(pkd);++i) {
		p = &pkd->pStore[i];
		if (p->dTEnc < *dNext) *dNext = p->dTEnc;
		}
	}

void
pkdMarkEncounters(PKD pkd,double dt)
{
	PARTICLE *p;
	int i;

	if (dt <= 0) { /* initialize */
		for (i=0;i<pkdLocal(pkd);++i) {
			p = &pkd->pStore[i];
			TYPESet(p,TYPE_ACTIVE);/*DEBUG redundant?*/
			p->iDriftType = KEPLER;
			p->dTEnc = HUGE_VAL;
/*DEBUG			p->idNbr.iPid = p->idNbr.iIndex =
				p->idNbr.iOrder = -1; */
			}
		}
	else {
		for (i=0;i<pkdLocal(pkd);++i) {
			if (p->dTEnc < dt) {
				TYPESet(p,TYPE_ACTIVE);
				p->iDriftType = NORMAL;
				}
			else {
				TYPEReset(p,TYPE_ACTIVE);
				p->iDriftType = KEPLER;
				}
			}
		}
	}

#endif /* COLLISIONS */
