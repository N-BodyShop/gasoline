#include <math.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <malloc.h>
#include <assert.h>
#include "mdl.h"
#include "pst.h"
#include "pkd.h"	
#include "outtype.h"
#include "smooth.h"


void pstInitialize(PST *ppst,MDL mdl,LCL *plcl)
{
	PST pst;
	int i;

	pst = (PST)malloc(sizeof(struct pstContext));
	assert(pst != NULL);
	*ppst = pst;
	pst->plcl = plcl;
	pst->mdl = mdl;
	pst->idSelf = mdlSelf(mdl);
	pst->pstLower = NULL;
	pst->idUpper = -1;	/* invalidate upper 'id' */
	pst->nLeaves = 1;
	pst->nLower = 0;
	pst->nUpper = 0;
	}


void pstFinish(PST pst)
{
	PST pstKill;

	while (pst) {
		pstKill = pst;
		pst = pst->pstLower;
		free(pstKill);
		}
	}


void pstSetAdd(PST pst,char *in,int nIn,char *out,int *pnOut)
{
	PST pstNew;

	if (pst->nLeaves > 1) {
		++pst->nLeaves;
		if (pst->nLower > pst->nUpper) {
			++pst->nUpper;
			mdlReqService(pst->mdl,pst->idUpper,PST_SETADD,in,nIn);
			mdlGetReply(pst->mdl,pst->idUpper,out,pnOut);
			}
		else {
			++pst->nLower;
			pstSetAdd(pst->pstLower,in,nIn,out,pnOut);
			}
		}
	else {
		++pst->nLeaves;
		++pst->nLower;
		++pst->nUpper;
		pstInitialize(&pstNew,pst->mdl,pst->plcl);
		pst->pstLower = pstNew;
		pst->idUpper = DATA(in,inSetAdd)->id;
		}
	if (pnOut) *pnOut = 0;
	}


void pstLevelize(PST pst,char *in,int nIn,char *out,int *pnOut)
{
	LCL *plcl = pst->plcl;
	int l,n;

	pst->iLvl = DATA(in,inLevelize)->iLvl;
	if (pst->nLeaves > 1) {
		++DATA(in,inLevelize)->iLvl;
		mdlReqService(pst->mdl,pst->idUpper,PST_LEVELIZE,in,nIn);
		pstLevelize(pst->pstLower,in,nIn,out,pnOut);
		mdlGetReply(pst->mdl,pst->idUpper,out,pnOut);
		}
	else {
		plcl->nPstLvl = pst->iLvl+1;

		/*
		 ** Placing this here is quite a HACK!!!
		 ** Allocate storage for kdTop.
		 */
		l = 1;
		n = mdlThreads(pst->mdl);
		while (n > l) {
			l = l << 1;
			}
		plcl->nTopNodes = l << 1;
		plcl->kdTop = malloc(plcl->nTopNodes*sizeof(KDN));	
		assert(plcl->kdTop != NULL);
		plcl->piLeaf = malloc(mdlThreads(pst->mdl)*sizeof(int));
		assert(plcl->piLeaf != NULL);

		}
	if (pnOut) *pnOut = 0;
	}


void pstShowPst(PST pst,char *in,int nIn,char *out,int *pnOut)
{
	char buf[100000];
	int i,iLvl,nBytes;
	
	iLvl = pst->iLvl;
	*out = 0;
	*pnOut = 1;
	strcpy(buf,"   ");
	for (i=0;i<iLvl;++i) {
		strcat(out,buf);
		*pnOut += strlen(buf);
		}
	sprintf(buf,"%d nLeaves:%d iLvl:%d\n",pst->idSelf,pst->nLeaves,pst->iLvl);
	strcat(out,buf);
	*pnOut += strlen(buf);
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_SHOWPST,in,nIn);
		pstShowPst(pst->pstLower,in,nIn,buf,&nBytes);
		strcat(out,buf);
		*pnOut += strlen(buf);
		mdlGetReply(pst->mdl,pst->idUpper,buf,&nBytes);
		strcat(out,buf);
		*pnOut += strlen(buf);
		}
	}


void pstReadTipsy(PST pst,char *in,int nIn,char *out,int *pnOut)
{
	int nTotal,nStore;
	char achInFile[PST_FILENAME_SIZE];
	LCL *plcl = pst->plcl;

	pst->nStart = DATA(in,inReadTipsy)->nStart;
	pst->nEnd = DATA(in,inReadTipsy)->nEnd;
	nTotal = pst->nEnd - pst->nStart + 1;
	if (pst->nLeaves > 1) {
		pst->nOrdSplit = pst->nStart + pst->nLower*nTotal/pst->nLeaves;
		DATA(in,inReadTipsy)->nStart = pst->nOrdSplit;
		mdlReqService(pst->mdl,pst->idUpper,PST_READTIPSY,in,nIn);
		DATA(in,inReadTipsy)->nStart = pst->nStart;
		DATA(in,inReadTipsy)->nEnd = pst->nOrdSplit - 1;
		pstReadTipsy(pst->pstLower,in,nIn,out,pnOut);
		DATA(in,inReadTipsy)->nEnd = pst->nEnd;
		mdlGetReply(pst->mdl,pst->idUpper,NULL,0);
		}
	else {
		/*
		 ** Add the local Data Path to the provided filename.
		 */
		achInFile[0] = 0;
		if (plcl->pszDataPath) {
			strcat(achInFile,plcl->pszDataPath);
			strcat(achInFile,"/");
			}
		strcat(achInFile,DATA(in,inReadTipsy)->achInFile);
		/*
		 ** Determine the size of the local particle store.
		 */
		nStore = nTotal + (int)ceil(nTotal*DATA(in,inReadTipsy)->fExtraStore);
		pkdInitialize(&plcl->pkd,pst->mdl,DATA(in,inReadTipsy)->iOrder,
					  nStore,plcl->nPstLvl,DATA(in,inReadTipsy)->fPeriod);
		pkdReadTipsy(plcl->pkd,achInFile,pst->nStart,nTotal);
		}
	if (pnOut) *pnOut = 0;
	}


int _pstRejMatch(PST pst,int n1,OREJ *p1,int n2,OREJ *p2,int *pidSwap)
{
	int id,i,i1,i2,nLarge,id1,id2;
	
	/*
	 ** First invalidate the pidSwap array.
	 */
	for (id=0;id<mdlThreads(pst->mdl);++id) pidSwap[id] = -1;
	/*
	 ** Now map largest nReject of p1 to largest nSpace of p2.
	 */
	while (1) {
		nLarge = 0;
		for (i=0;i<n1;++i) {
			if (p1[i].nRejects > nLarge) {
				nLarge = p1[i].nRejects;
				i1 = i;
				}
			}
		if (nLarge == 0) break;
		nLarge = 0;
		for (i=0;i<n2;++i) {
			if (p2[i].nSpace > nLarge) {
				nLarge = p2[i].nSpace;
				i2 = i;
				}
			}
		if (nLarge == 0) break;
		p1[i1].nRejects = 0;
		p1[i1].nSpace = 0;
		p2[i2].nRejects = 0;
		p2[i2].nSpace = 0;
		id1 = p1[i1].id;
		id2 = p2[i2].id;
		pidSwap[id1] = id2;
		pidSwap[id2] = id1;
		}
	/*
	 ** Now map largest nReject of p2 to largest nSpace of p1.
	 ** However, already mapped stuff is ignored, by the above!
	 */
	while (1) {
		nLarge = 0;
		for (i=0;i<n2;++i) {
			if (p2[i].nRejects > nLarge) {
				nLarge = p2[i].nRejects;
				i2 = i;
				}
			}
		if (nLarge == 0) break;
		nLarge = 0;
		for (i=0;i<n1;++i) {
			if (p1[i].nSpace > nLarge) {
				nLarge = p1[i].nSpace;
				i1 = i;
				}
			}
		if (nLarge == 0) break;
		p1[i1].nRejects = 0;
		p1[i1].nSpace = 0;
		p2[i2].nRejects = 0;
		p2[i2].nSpace = 0;
		id1 = p1[i1].id;
		id2 = p2[i2].id;
		pidSwap[id1] = id2;
		pidSwap[id2] = id1;
		}
	for (i=0;i<mdlThreads(pst->mdl);++i)
		if (pidSwap[i] != -1) return(1);
	return(0);
	}


void _pstRootSplit(PST pst,int iSplitDim)
{
	int d,ittr,iDum;
	int nLow,nHigh,nLowerStore,nUpperStore;
	float fLow,fHigh;
	float fl,fu,fm;
	char outFree[SIZE(outFreeStore)];
	char inWt[SIZE(inWeight)];
	char outWtLow[SIZE(outWeight)];
	char outWtHigh[SIZE(outWeight)];
	char inCol[SIZE(inColRejects)];
	OREJ *pLowerRej,*pUpperRej;
	int *pidSwap,iRet;

	/*
	 ** First find out how much free storage there is available for particles
	 ** on the lower and upper subset of processors.
	 */
	mdlReqService(pst->mdl,pst->idUpper,PST_FREESTORE,NULL,0);
	pstFreeStore(pst->pstLower,NULL,0,outFree,&iDum);
	nLowerStore = DATA(outFree,outFreeStore)->nFreeStore;
	mdlGetReply(pst->mdl,pst->idUpper,outFree,&iDum);
	nUpperStore = DATA(outFree,outFreeStore)->nFreeStore;
	/*
	 ** Now start the ROOT finder based on balancing weight ALONE!
	 */
	d = iSplitDim;
	fl = pst->bnd.fMin[d];
	fu = pst->bnd.fMax[d];
	fm = (fl + fu)/2;
	ittr = 0;
	while (fl < fm && fm < fu) {
		DATA(inWt,inWeight)->iSplitDim = d;
		DATA(inWt,inWeight)->fSplit = fm;
		DATA(inWt,inWeight)->ittr = ittr;
		DATA(inWt,inWeight)->iSplitSide = 1;
		mdlReqService(pst->mdl,pst->idUpper,PST_WEIGHT,inWt,SIZE(inWeight));
		DATA(inWt,inWeight)->iSplitSide = 0;
		pstWeight(pst->pstLower,inWt,SIZE(inWeight),outWtLow,&iDum);
		mdlGetReply(pst->mdl,pst->idUpper,outWtHigh,&iDum);
		/*
		 ** Add lower and Upper subsets weights and numbers
		 */
		nLow = DATA(outWtLow,outWeight)->nLow +
			DATA(outWtHigh,outWeight)->nLow;
		nHigh = DATA(outWtLow,outWeight)->nHigh +
			DATA(outWtHigh,outWeight)->nHigh;
		fLow = DATA(outWtLow,outWeight)->fLow +
			DATA(outWtHigh,outWeight)->fLow;
		fHigh = DATA(outWtLow,outWeight)->fHigh +
			DATA(outWtHigh,outWeight)->fHigh;
/*
		printf("ittr:%d l:%d u:%d lw:%f uw:%f\n",ittr,nLow,nHigh,fLow,fHigh);
*/
		if (fLow/pst->nLower > fHigh/pst->nUpper) fu = fm;
		else if (fLow/pst->nLower < fHigh/pst->nUpper) fl = fm;
		else break;
		fm = (fl + fu)/2;
		++ittr;
		}
	/*
	 ** If we exceed the local pStore in the lower or upper subsets then
	 ** we have to relax the condition that work be balanced and try to
	 ** max out the number of particles in the subset which had too many.
	 */
	if (nLow > nLowerStore) {
		fl = pst->bnd.fMin[d];
		fu = fm;
		fm = (fl + fu)/2;
		ittr = 0;
	    while (fl < fm && fm < fu) {
			DATA(inWt,inWeight)->iSplitDim = d;
			DATA(inWt,inWeight)->fSplit = fm;
			DATA(inWt,inWeight)->ittr = ittr;
			DATA(inWt,inWeight)->iSplitSide = 1;
			mdlReqService(pst->mdl,pst->idUpper,PST_WEIGHT,inWt,
						  SIZE(inWeight));
			DATA(inWt,inWeight)->iSplitSide = 0;
			pstWeight(pst->pstLower,inWt,SIZE(inWeight),outWtLow,&iDum);
			mdlGetReply(pst->mdl,pst->idUpper,outWtHigh,&iDum);
			/*
			 ** Add lower and Upper subsets weights and numbers
			 */
			nLow = DATA(outWtLow,outWeight)->nLow +
			DATA(outWtHigh,outWeight)->nLow;
/*
			printf("Fit ittr:%d l:%d\n",ittr,nLow);
*/
			if (nLow > nLowerStore) fu = fm;
			else if (nLow < nLowerStore) fl = fm;
			else {
				fl = fm;
				break;
				}
			fm = (fl + fu)/2;
			++ittr;
			}
		fm = fl;
		}
	else if (nHigh > nUpperStore) {
		fl = fm;
		fu = pst->bnd.fMax[d];
		fm = (fl + fu)/2;
		ittr = 0;
	    while (fl < fm && fm < fu) {
			DATA(inWt,inWeight)->iSplitDim = d;
			DATA(inWt,inWeight)->fSplit = fm;
			DATA(inWt,inWeight)->ittr = ittr;
			DATA(inWt,inWeight)->iSplitSide = 1;
			mdlReqService(pst->mdl,pst->idUpper,PST_WEIGHT,inWt,
						  SIZE(inWeight));
			DATA(inWt,inWeight)->iSplitSide = 0;
			pstWeight(pst->pstLower,inWt,SIZE(inWeight),outWtLow,&iDum);
			mdlGetReply(pst->mdl,pst->idUpper,outWtHigh,&iDum);
			/*
			 ** Add lower and Upper subsets weights and numbers
			 */
			nHigh = DATA(outWtLow,outWeight)->nHigh +
			DATA(outWtHigh,outWeight)->nHigh;
/*
			printf("Fit ittr:%d u:%d\n",ittr,nHigh);
*/
			if (nHigh > nUpperStore) fl = fm;
			else if (nHigh < nUpperStore) fu = fm;
			else {
				fu = fm;
				break;
				}
			fm = (fl + fu)/2;
			++ittr;
			}
		fm = fu;
		}
	pst->fSplit = fm;
	pst->iSplitDim = d;
	/*
	 ** First Collect rejects.
	 */
	pLowerRej = (OREJ *)malloc(pst->nLower*sizeof(OREJ));
	assert(pLowerRej != NULL);
	pUpperRej = (OREJ *)malloc(pst->nLower*sizeof(OREJ));
	assert(pUpperRej != NULL);
	pidSwap = (int *)malloc(mdlThreads(pst->mdl)*sizeof(int));
	assert(pidSwap != NULL);

	DATA(inCol,inColRejects)->fSplit = pst->fSplit;
	DATA(inCol,inColRejects)->iSplitDim = pst->iSplitDim;
	DATA(inCol,inColRejects)->iSplitSide = 1;
	mdlReqService(pst->mdl,pst->idUpper,PST_COLREJECTS,inCol,
				  SIZE(inColRejects));

	DATA(inCol,inColRejects)->iSplitSide = 0;
	pstColRejects(pst->pstLower,inCol,SIZE(inColRejects),
				  (char *)pLowerRej,&iDum);
	assert(iDum/sizeof(OREJ) == pst->nLower);

	mdlGetReply(pst->mdl,pst->idUpper,(char *)pUpperRej,&iDum);
	assert(iDum/sizeof(OREJ) == pst->nUpper);

	while (1) {
		iRet = _pstRejMatch(pst,pst->nLower,pLowerRej,
							pst->nUpper,pUpperRej,pidSwap);
		if (!iRet) break;
		mdlReqService(pst->mdl,pst->idUpper,PST_SWAPREJECTS,(char *)pidSwap,
					  mdlThreads(pst->mdl)*sizeof(int));
		pstSwapRejects(pst->pstLower,(char *)pidSwap,
					   mdlThreads(pst->mdl)*sizeof(int),
					   (char *)pLowerRej,&iDum);
		assert(iDum/sizeof(OREJ) == pst->nLower);

		mdlGetReply(pst->mdl,pst->idUpper,(char *)pUpperRej,&iDum);
		assert(iDum/sizeof(OREJ) == pst->nUpper);
		}
	free(pLowerRej);
	free(pUpperRej);
	free(pidSwap);
	}


void pstDomainDecomp(PST pst,char *in,int nIn,char *out,int *pnOut)
{
	int nOut,d,j;
	char outBnd[SIZE(outCalcBound)];
	
	if (pst->nLeaves > 1) {
		/*
		 ** First calculate the Bound for the set.
		 */
		pstCalcBound(pst,NULL,0,outBnd,&nOut);
		pst->bnd = DATA(outBnd,outCalcBound)->bnd;
		/*
		 ** Next determine the longest axis based on the bounds.
		 */
		d = 0;
		for (j=1;j<3;++j) {
			if (pst->bnd.fMax[j]-pst->bnd.fMin[j] > 
				pst->bnd.fMax[d]-pst->bnd.fMin[d]) d = j;
			}
		pst->iSplitDim = d;
		_pstRootSplit(pst,d);
		/*
		 ** Now go on to DD of next levels.
		 */
		if (pst->nUpper > 1) 
			mdlReqService(pst->mdl,pst->idUpper,PST_DOMAINDECOMP,in,nIn);
		if (pst->nLower > 1) 
			pstDomainDecomp(pst->pstLower,in,nIn,out,pnOut);
		if (pst->nUpper > 1) 
			mdlGetReply(pst->mdl,pst->idUpper,out,pnOut);
		}
	if (pnOut) *pnOut = 0;
	}


void pstCalcBound(PST pst,char *in,int nIn,char *out,int *pnOut)
{
	LCL *plcl = pst->plcl;
	int j;
	char outBnd[SIZE(outCalcBound)];
	BND bnd;
	
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_CALCBOUND,in,nIn);
		pstCalcBound(pst->pstLower,in,nIn,out,pnOut);
		mdlGetReply(pst->mdl,pst->idUpper,outBnd,pnOut);
		bnd = DATA(outBnd,outCalcBound)->bnd;
		for (j=0;j<3;++j) {
			if (bnd.fMin[j] < DATA(out,outCalcBound)->bnd.fMin[j]) {
				DATA(out,outCalcBound)->bnd.fMin[j] = bnd.fMin[j];
				}
			if (bnd.fMax[j] > DATA(out,outCalcBound)->bnd.fMax[j]) {
				DATA(out,outCalcBound)->bnd.fMax[j] = bnd.fMax[j];
				}
			}
		}
	else {
		pkdCalcBound(plcl->pkd,&DATA(out,outCalcBound)->bnd);
		*pnOut = SIZE(outCalcBound);
		}
	}


void pstWeight(PST pst,char *in,int nIn,char *out,int *pnOut)
{
	LCL *plcl = pst->plcl;
	char outWt[SIZE(outWeight)];
	float fSplit,fLow,fHigh;
	int iSplitSide;

	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_WEIGHT,in,nIn);
		pstWeight(pst->pstLower,in,nIn,out,pnOut);
		mdlGetReply(pst->mdl,pst->idUpper,outWt,pnOut);
		DATA(out,outWeight)->nLow += DATA(outWt,outWeight)->nLow;
		DATA(out,outWeight)->nHigh += DATA(outWt,outWeight)->nHigh;
		DATA(out,outWeight)->fLow += DATA(outWt,outWeight)->fLow;
		DATA(out,outWeight)->fHigh += DATA(outWt,outWeight)->fHigh;
		}
	else {
		fSplit = DATA(in,inWeight)->fSplit;
		iSplitSide = DATA(in,inWeight)->iSplitSide;
		if (DATA(in,inWeight)->ittr == 0) {
			/*
			 ** Initialize.
			 */
			plcl->fSplit = fSplit;
			plcl->iWtFrom = 0;
			plcl->iWtTo = pkdLocal(plcl->pkd)-1;
			plcl->fWtLow = 0.0;
			plcl->fWtHigh = 0.0;
			}
		else {
			/*
			 ** Update the Weight Sums and use smaller weight region.
			 */
			if (fSplit < plcl->fSplit) {
				plcl->fWtHigh += plcl->fHigh;
				if (iSplitSide) plcl->iWtFrom = plcl->iPart;
				else plcl->iWtTo = plcl->iPart-1;
				}
			else {
				plcl->fWtLow += plcl->fLow;
				if (iSplitSide) plcl->iWtTo = plcl->iPart-1;
				else plcl->iWtFrom = plcl->iPart;
				}
			plcl->fSplit = fSplit;
			}
		plcl->iPart = pkdWeight(plcl->pkd,DATA(in,inWeight)->iSplitDim,
								fSplit,iSplitSide,
								plcl->iWtFrom,plcl->iWtTo,
								&DATA(out,outWeight)->nLow,
								&DATA(out,outWeight)->nHigh,
								&fLow,&fHigh);
		DATA(out,outWeight)->fLow = fLow + plcl->fWtLow;
		DATA(out,outWeight)->fHigh = fHigh + plcl->fWtHigh;
		plcl->fLow = fLow;
		plcl->fHigh = fHigh;
		*pnOut = SIZE(outWeight);
		}
	}


void pstFreeStore(PST pst,char *in,int nIn,char *out,int *pnOut)
{
	LCL *plcl = pst->plcl;
	int nLowerStore,nUpperStore;

	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_FREESTORE,in,nIn);
		pstFreeStore(pst->pstLower,in,nIn,out,pnOut);
		nLowerStore = DATA(out,outFreeStore)->nFreeStore;
		mdlGetReply(pst->mdl,pst->idUpper,out,pnOut);
		nUpperStore = DATA(out,outFreeStore)->nFreeStore;
		DATA(out,outFreeStore)->nFreeStore = nLowerStore + nUpperStore;
		}
	else {
		DATA(out,outFreeStore)->nFreeStore = pkdFreeStore(plcl->pkd);
		*pnOut = SIZE(outFreeStore);
		}
	}


void pstColRejects(PST pst,char *in,int nIn,char *out,int *pnOut)
{
	LCL *plcl = pst->plcl;
	int nLower,nUpper,iUpper;
	OREJ *pOutRej = (OREJ *)out;
	
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_COLREJECTS,in,nIn);
		pstColRejects(pst->pstLower,in,nIn,(char *)&pOutRej[0],&nLower);
		iUpper = nLower/sizeof(OREJ);
		mdlGetReply(pst->mdl,pst->idUpper,(char *)&pOutRej[iUpper],&nUpper);
		*pnOut = nLower + nUpper;
		}
	else {
	    pOutRej->nRejects = pkdColRejects(plcl->pkd,
										  DATA(in,inColRejects)->iSplitDim,
										  DATA(in,inColRejects)->fSplit,
										  DATA(in,inColRejects)->iSplitSide);
		pOutRej->nSpace = pkdSwapSpace(plcl->pkd);
		pOutRej->id = pst->idSelf;
		*pnOut = sizeof(OREJ);
		}
	}


void pstColOrdRejects(PST pst,char *in,int nIn,char *out,int *pnOut)
{
	LCL *plcl = pst->plcl;
	int nLower,nUpper,iUpper;
	OREJ *pOutRej = (OREJ *)out;
	
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_COLORDREJECTS,in,nIn);
		pstColOrdRejects(pst->pstLower,in,nIn,(char *)&pOutRej[0],&nLower);
		iUpper = nLower/sizeof(OREJ);
		mdlGetReply(pst->mdl,pst->idUpper,(char *)&pOutRej[iUpper],&nUpper);
		*pnOut = nLower + nUpper;
		}
	else {
	    pOutRej->nRejects =
			pkdColOrdRejects(plcl->pkd,
							 DATA(in,inColOrdRejects)->nOrdSplit,
							 DATA(in,inColOrdRejects)->iSplitSide);
		pOutRej->nSpace = pkdSwapSpace(plcl->pkd);
		pOutRej->id = pst->idSelf;
		*pnOut = sizeof(OREJ);
		}
	}


void pstSwapRejects(PST pst,char *in,int nIn,char *out,int *pnOut)
{
	LCL *plcl = pst->plcl;
	int nLower,nUpper,iUpper,idSwap;
	OREJ *pOutRej = (OREJ *)out;
	int *pidSwap = (int *)in;

	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_SWAPREJECTS,in,nIn);
		pstSwapRejects(pst->pstLower,in,nIn,(char *)&pOutRej[0],&nLower);
		iUpper = nLower/sizeof(OREJ);
		mdlGetReply(pst->mdl,pst->idUpper,(char *)&pOutRej[iUpper],&nUpper);
		*pnOut = nLower + nUpper;
		}
	else {
		idSwap = pidSwap[pst->idSelf];
		pOutRej->nRejects = pkdSwapRejects(plcl->pkd,idSwap);
		pOutRej->nSpace = pkdSwapSpace(plcl->pkd);
		pOutRej->id = pst->idSelf;
		*pnOut = sizeof(OREJ);
		}
	}


void pstDomainColor(PST pst,char *in,int nIn,char *out,int *pnOut)
{
	LCL *plcl = pst->plcl;

	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_DOMAINCOLOR,in,nIn);
		pstDomainColor(pst->pstLower,in,nIn,out,pnOut);
		mdlGetReply(pst->mdl,pst->idUpper,out,pnOut);
		}
	else {
		pkdDomainColor(plcl->pkd);
		}
	if (pnOut) *pnOut = 0;
	}


void _pstOrdSplit(PST pst)
{
	char inCol[SIZE(inColOrdRejects)];
	int iDum;
	OREJ *pLowerRej,*pUpperRej;
	int *pidSwap,iRet;

	/*
	 ** First Collect rejects.
	 */
	pLowerRej = (OREJ *)malloc(pst->nLower*sizeof(OREJ));
	assert(pLowerRej != NULL);
	pUpperRej = (OREJ *)malloc(pst->nLower*sizeof(OREJ));
	assert(pUpperRej != NULL);
	pidSwap = (int *)malloc(mdlThreads(pst->mdl)*sizeof(int));
	assert(pidSwap != NULL);

	DATA(inCol,inColOrdRejects)->nOrdSplit = pst->nOrdSplit;
	DATA(inCol,inColOrdRejects)->iSplitSide = 1;
	mdlReqService(pst->mdl,pst->idUpper,PST_COLORDREJECTS,inCol,
				  SIZE(inColOrdRejects));

	DATA(inCol,inColOrdRejects)->iSplitSide = 0;
	pstColOrdRejects(pst->pstLower,inCol,SIZE(inColOrdRejects),
					 (char *)pLowerRej,&iDum);
	assert(iDum/sizeof(OREJ) == pst->nLower);

	mdlGetReply(pst->mdl,pst->idUpper,(char *)pUpperRej,&iDum);
	assert(iDum/sizeof(OREJ) == pst->nUpper);

	while (1) {
		iRet = _pstRejMatch(pst,pst->nLower,pLowerRej,
							pst->nUpper,pUpperRej,pidSwap);
		if (!iRet) break;
		mdlReqService(pst->mdl,pst->idUpper,PST_SWAPREJECTS,(char *)pidSwap,
					  mdlThreads(pst->mdl)*sizeof(int));
		pstSwapRejects(pst->pstLower,(char *)pidSwap,
					   mdlThreads(pst->mdl)*sizeof(int),
					   (char *)pLowerRej,&iDum);
		assert(iDum/sizeof(OREJ) == pst->nLower);

		mdlGetReply(pst->mdl,pst->idUpper,(char *)pUpperRej,&iDum);
		assert(iDum/sizeof(OREJ) == pst->nUpper);
		}
	free(pLowerRej);
	free(pUpperRej);
	free(pidSwap);
	}


void pstDomainOrder(PST pst,char *in,int nIn,char *out,int *pnOut)
{
	if (pst->nLeaves > 1) {
		_pstOrdSplit(pst);
		/*
		 ** Now go on to Domain Order of next levels.
		 */
		if (pst->nUpper > 1) 
			mdlReqService(pst->mdl,pst->idUpper,PST_DOMAINORDER,in,nIn);
		if (pst->nLower > 1) 
			pstDomainOrder(pst->pstLower,in,nIn,out,pnOut);
		if (pst->nUpper > 1) 
			mdlGetReply(pst->mdl,pst->idUpper,out,pnOut);
		}
	if (pnOut) *pnOut = 0;
	}


void pstLocalOrder(PST pst,char *in,int nIn,char *out,int *pnOut)
{
	LCL *plcl = pst->plcl;

	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_LOCALORDER,in,nIn);
		pstLocalOrder(pst->pstLower,in,nIn,out,pnOut);
		mdlGetReply(pst->mdl,pst->idUpper,out,pnOut);
		}
	else {
		pkdLocalOrder(plcl->pkd,pst->nStart);
		}
	if (pnOut) *pnOut = 0;
	}


void pstOutArray(PST pst,char *in,int nIn,char *out,int *pnOut)
{
	char achOutFile[PST_FILENAME_SIZE];
	LCL *plcl = pst->plcl;

	if (pst->nLeaves > 1) {
		/*
		 ** Non-Recursive Text output.
		 */
		pstOutArray(pst->pstLower,in,nIn,out,pnOut);
		mdlReqService(pst->mdl,pst->idUpper,PST_OUTARRAY,in,nIn);
		mdlGetReply(pst->mdl,pst->idUpper,NULL,0);
		}
	else {
		/*
		 ** Add the local Data Path to the provided filename.
		 */
		achOutFile[0] = 0;
		if (plcl->pszDataPath) {
			strcat(achOutFile,plcl->pszDataPath);
			strcat(achOutFile,"/");
			}
		strcat(achOutFile,DATA(in,inOutArray)->achOutFile);
		pkdOutArray(plcl->pkd,achOutFile,DATA(in,inOutArray)->iType);
		}
	if (pnOut) *pnOut = 0;
	}


void pstOutVector(PST pst,char *in,int nIn,char *out,int *pnOut)
{
	char achOutFile[PST_FILENAME_SIZE];
	LCL *plcl = pst->plcl;

	if (pst->nLeaves > 1) {
		/*
		 ** Non-Recursive Text output.
		 */
		pstOutVector(pst->pstLower,in,nIn,out,pnOut);
		mdlReqService(pst->mdl,pst->idUpper,PST_OUTVECTOR,in,nIn);
		mdlGetReply(pst->mdl,pst->idUpper,NULL,0);
		}
	else {
		/*
		 ** Add the local Data Path to the provided filename.
		 */
		achOutFile[0] = 0;
		if (plcl->pszDataPath) {
			strcat(achOutFile,plcl->pszDataPath);
			strcat(achOutFile,"/");
			}
		strcat(achOutFile,DATA(in,inOutVector)->achOutFile);
		pkdOutVector(plcl->pkd,achOutFile,DATA(in,inOutVector)->iDim,
					 DATA(in,inOutVector)->iType);
		}
	if (pnOut) *pnOut = 0;
	}


void pstWriteTipsy(PST pst,char *in,int nIn,char *out,int *pnOut)
{
	char achOutFile[PST_FILENAME_SIZE];
	LCL *plcl = pst->plcl;

	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_WRITETIPSY,in,nIn);
		pstWriteTipsy(pst->pstLower,in,nIn,out,pnOut);
		mdlGetReply(pst->mdl,pst->idUpper,NULL,0);
		}
	else {
		/*
		 ** Add the local Data Path to the provided filename.
		 */
		achOutFile[0] = 0;
		if (plcl->pszDataPath) {
			strcat(achOutFile,plcl->pszDataPath);
			strcat(achOutFile,"/");
			}
		strcat(achOutFile,DATA(in,inWriteTipsy)->achOutFile);
		pkdWriteTipsy(plcl->pkd,achOutFile,pst->nStart,pst->nEnd);
		}
	if (pnOut) *pnOut = 0;
	}


void pstSetSoft(PST pst,char *in,int nIn,char *out,int *pnOut)
{
	LCL *plcl = pst->plcl;

	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_SETSOFT,in,nIn);
		pstSetSoft(pst->pstLower,in,nIn,out,pnOut);
		mdlGetReply(pst->mdl,pst->idUpper,out,pnOut);
		}
	else {
                pkdSetSoft(plcl->pkd,DATA(in,inSetSoft)->dSoft);
		}
	if (pnOut) *pnOut = 0;
	}

void pstBuildTree(PST pst,char *in,int nIn,char *out,int *pnOut)
{
	LCL *plcl = pst->plcl;
	int iCell;
	
	iCell = DATA(in,inBuildTree)->iCell;
	if (pst->nLeaves > 1) {
		plcl->kdTop[iCell].iDim = pst->iSplitDim;
		plcl->kdTop[iCell].fSplit = pst->fSplit;
		plcl->kdTop[iCell].bnd = pst->bnd;
		plcl->kdTop[iCell].pLower = pst->idSelf;
		plcl->kdTop[iCell].pUpper = pst->idUpper;
		DATA(in,inBuildTree)->iCell = UPPER(iCell);
		mdlReqService(pst->mdl,pst->idUpper,PST_BUILDTREE,in,nIn);
		DATA(in,inBuildTree)->iCell = LOWER(iCell);
		pstBuildTree(pst->pstLower,in,nIn,NULL,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
		}
	else {
		pkdBuildLocal(plcl->pkd,DATA(in,inBuildTree)->nBucket,
					  DATA(in,inBuildTree)->iOpenType,
					  DATA(in,inBuildTree)->dCrit,
					  &plcl->kdTop[iCell]);
		plcl->kdTop[iCell].pLower = pst->idSelf;
		plcl->kdTop[iCell].pUpper = -1;
		pkdBuildTop(plcl->pkd,DATA(in,inBuildTree)->iOpenType,
					DATA(in,inBuildTree)->dCrit,
					plcl->kdTop,plcl->nTopNodes,plcl->piLeaf);
		}
	if (pnOut) *pnOut = 0;
	}


void pstDensity(PST pst,char *in,int nIn,char *out,int *pnOut)
{
	LCL *plcl = pst->plcl;
	SMX smx;
	int nSmooth,bGatherScatter;

	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_DENSITY,in,nIn);
		pstDensity(pst->pstLower,in,nIn,NULL,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
		}
	else {
		nSmooth = DATA(in,inDensity)->nSmooth;
		bGatherScatter = DATA(in,inDensity)->bGatherScatter;
		smInitialize(&smx,plcl->pkd,nSmooth,bGatherScatter);
		if (bGatherScatter) {
			smSmooth(smx,smDensitySym);
			}
		else {
			smSmooth(smx,smDensity);
			}
		smFinish(smx);
		}
	if (pnOut) *pnOut = 0;
	}


void pstGravity(PST pst,char *in,int nIn,char *out,int *pnOut)
{
	LCL *plcl = pst->plcl;
	SMX smx;
	double dPartSum,dCellSum;
	char outUp[SIZE(outGravity)];
	int iDum;

	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_GRAVITY,in,nIn);
		pstGravity(pst->pstLower,in,nIn,out,&iDum);
		assert(iDum == SIZE(outGravity));
		mdlGetReply(pst->mdl,pst->idUpper,outUp,&iDum);
		assert(iDum == SIZE(outGravity));
		DATA(out,outGravity)->dPartSum += DATA(outUp,outGravity)->dPartSum;
		DATA(out,outGravity)->dCellSum += DATA(outUp,outGravity)->dCellSum;
		DATA(out,outGravity)->dWSum += DATA(outUp,outGravity)->dWSum;
		DATA(out,outGravity)->dISum += DATA(outUp,outGravity)->dISum;
		DATA(out,outGravity)->dESum += DATA(outUp,outGravity)->dESum;
		if (DATA(outUp,outGravity)->dWMax > DATA(out,outGravity)->dWMax)
			DATA(out,outGravity)->dWMax = DATA(outUp,outGravity)->dWMax;
		if (DATA(outUp,outGravity)->dIMax > DATA(out,outGravity)->dIMax)
			DATA(out,outGravity)->dIMax = DATA(outUp,outGravity)->dIMax;
		if (DATA(outUp,outGravity)->dEMax > DATA(out,outGravity)->dEMax)
			DATA(out,outGravity)->dEMax = DATA(outUp,outGravity)->dEMax;
		if (DATA(outUp,outGravity)->dWMin < DATA(out,outGravity)->dWMin)
			DATA(out,outGravity)->dWMin = DATA(outUp,outGravity)->dWMin;
		if (DATA(outUp,outGravity)->dIMin < DATA(out,outGravity)->dIMin)
			DATA(out,outGravity)->dIMin = DATA(outUp,outGravity)->dIMin;
		if (DATA(outUp,outGravity)->dEMin < DATA(out,outGravity)->dEMin)
			DATA(out,outGravity)->dEMin = DATA(outUp,outGravity)->dEMin;
		}
	else {
		pkdGravAll(plcl->pkd,DATA(in,inGravity)->nReps,
				   DATA(in,inGravity)->bPeriodic,DATA(in,inGravity)->dEwCut,
				   DATA(in,inGravity)->dEwhCut,
				   &DATA(out,outGravity)->dPartSum,
				   &DATA(out,outGravity)->dCellSum);
		DATA(out,outGravity)->dWSum = pkdGetTimer(plcl->pkd,1);
		DATA(out,outGravity)->dISum = pkdGetTimer(plcl->pkd,2);
		DATA(out,outGravity)->dESum = pkdGetTimer(plcl->pkd,3);
		DATA(out,outGravity)->dWMax = DATA(out,outGravity)->dWSum;
		DATA(out,outGravity)->dIMax = DATA(out,outGravity)->dISum;
		DATA(out,outGravity)->dEMax = DATA(out,outGravity)->dESum;
		DATA(out,outGravity)->dWMin = DATA(out,outGravity)->dWSum;
		DATA(out,outGravity)->dIMin = DATA(out,outGravity)->dISum;
		DATA(out,outGravity)->dEMin = DATA(out,outGravity)->dESum;
		}
	*pnOut = SIZE(outGravity);
	}


void pstCalcE(PST pst,char *in,int nIn,char *out,int *pnOut)
{
	char outE[SIZE(outCalcE)];
	LCL *plcl = pst->plcl;
	int iDum;

	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_CALCE,in,nIn);
		pstCalcE(pst->pstLower,in,nIn,out,&iDum);
		mdlGetReply(pst->mdl,pst->idUpper,outE,&iDum);
		DATA(out,outCalcE)->T += DATA(outE,outCalcE)->T;
		DATA(out,outCalcE)->U += DATA(outE,outCalcE)->U;
		}
	else {
		pkdCalcE(plcl->pkd,&DATA(out,outCalcE)->T,&DATA(out,outCalcE)->U);
		}
	*pnOut = SIZE(outCalcE);
	}


void pstDrift(PST pst,char *in,int nIn,char *out,int *pnOut)
{
	LCL *plcl = pst->plcl;

	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_DRIFT,in,nIn);
		pstDrift(pst->pstLower,in,nIn,out,pnOut);
		mdlGetReply(pst->mdl,pst->idUpper,out,pnOut);
		}
	else {
		pkdDrift(plcl->pkd,DATA(in,inDrift)->dDelta,
				 DATA(in,inDrift)->fCenter);
		}
	if (pnOut) *pnOut = 0;
	}


void pstKick(PST pst,char *in,int nIn,char *out,int *pnOut)
{
	LCL *plcl = pst->plcl;

	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_KICK,in,nIn);
		pstKick(pst->pstLower,in,nIn,out,pnOut);
		mdlGetReply(pst->mdl,pst->idUpper,out,pnOut);
		}
	else {
		pkdKick(plcl->pkd,DATA(in,inKick)->dvFacOne,
				DATA(in,inKick)->dvFacTwo);
		}
	if (pnOut) *pnOut = 0;
	}


void pstReadCheck(PST pst,char *in,int nIn,char *out,int *pnOut)
{
	int nTotal,nStore;
	char achInFile[PST_FILENAME_SIZE];
	LCL *plcl = pst->plcl;

	pst->nStart = DATA(in,inReadCheck)->nStart;
	pst->nEnd = DATA(in,inReadCheck)->nEnd;
	nTotal = pst->nEnd - pst->nStart + 1;
	if (pst->nLeaves > 1) {
		pst->nOrdSplit = pst->nStart + pst->nLower*nTotal/pst->nLeaves;
		DATA(in,inReadCheck)->nStart = pst->nOrdSplit;
		mdlReqService(pst->mdl,pst->idUpper,PST_READCHECK,in,nIn);
		DATA(in,inReadCheck)->nStart = pst->nStart;
		DATA(in,inReadCheck)->nEnd = pst->nOrdSplit - 1;
		pstReadCheck(pst->pstLower,in,nIn,out,pnOut);
		DATA(in,inReadCheck)->nEnd = pst->nEnd;
		mdlGetReply(pst->mdl,pst->idUpper,NULL,0);
		}
	else {
		/*
		 ** Add the local Data Path to the provided filename.
		 */
		achInFile[0] = 0;
		if (plcl->pszDataPath) {
			strcat(achInFile,plcl->pszDataPath);
			strcat(achInFile,"/");
			}
		strcat(achInFile,DATA(in,inReadCheck)->achInFile);
		/*
		 ** Determine the size of the local particle store.
		 */
		nStore = nTotal + (int)ceil(nTotal*DATA(in,inReadCheck)->fExtraStore);
		pkdInitialize(&plcl->pkd,pst->mdl,DATA(in,inReadCheck)->iOrder,
					  nStore,plcl->nPstLvl,DATA(in,inReadCheck)->fPeriod);
		pkdReadCheck(plcl->pkd,achInFile,pst->nStart,nTotal);
		}
	if (pnOut) *pnOut = 0;
	}


void pstSetTotal(PST pst,char *in,int nIn,char *out,int *pnOut)
{	
	char oute[SIZE(outSetTotal)];
	int iDum;
	LCL *plcl = pst->plcl;

	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_SETTOTAL,in,nIn);
		pstSetTotal(pst->pstLower,in,nIn,out,pnOut);
		mdlGetReply(pst->mdl,pst->idUpper,oute,&iDum);
		DATA(out,outSetTotal)->nTotal += DATA(oute,outSetTotal)->nTotal;
		pst->nTotal = DATA(out,outSetTotal)->nTotal;
		}
	else {
		pst->nTotal = pkdLocal(plcl->pkd);
		DATA(out,outSetTotal)->nTotal = pst->nTotal;
		}
	*pnOut = SIZE(outSetTotal);
	}


void pstWriteCheck(PST pst,char *in,int nIn,char *out,int *pnOut)
{
	int nStart;
	char achOutFile[PST_FILENAME_SIZE];
	LCL *plcl = pst->plcl;

	if (pst->nLeaves > 1) {
		nStart = DATA(in,inWriteCheck)->nStart;
		DATA(in,inWriteCheck)->nStart += pst->pstLower->nTotal;
		mdlReqService(pst->mdl,pst->idUpper,PST_WRITECHECK,in,nIn);
		DATA(in,inWriteCheck)->nStart = nStart;
		pstWriteCheck(pst->pstLower,in,nIn,out,pnOut);
		mdlGetReply(pst->mdl,pst->idUpper,NULL,0);
		}
	else {
		/*
		 ** Add the local Data Path to the provided filename.
		 */
		achOutFile[0] = 0;
		if (plcl->pszDataPath) {
			strcat(achOutFile,plcl->pszDataPath);
			strcat(achOutFile,"/");
			}
		strcat(achOutFile,DATA(in,inWriteCheck)->achOutFile);
		pkdWriteCheck(plcl->pkd,achOutFile,DATA(in,inWriteCheck)->nStart);
		}
	if (pnOut) *pnOut = 0;
	}



