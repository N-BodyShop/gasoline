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


void pstSetAdd(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	PST pstNew;
	struct inSetAdd *in = vin;

	assert(nIn == sizeof(struct inSetAdd));
	if (pst->nLeaves > 1) {
		++pst->nLeaves;
		if (pst->nLower > pst->nUpper) {
			++pst->nUpper;
			mdlReqService(pst->mdl,pst->idUpper,PST_SETADD,in,nIn);
			mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
			}
		else {
			++pst->nLower;
			pstSetAdd(pst->pstLower,in,nIn,NULL,NULL);
			}
		}
	else {
		++pst->nLeaves;
		++pst->nLower;
		++pst->nUpper;
		pstInitialize(&pstNew,pst->mdl,pst->plcl);
		pst->pstLower = pstNew;
		pst->idUpper = in->id;
		}
	if (pnOut) *pnOut = 0;
	}


void pstLevelize(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inLevelize *in = vin;
	int l,n;

	assert(nIn == sizeof(struct inLevelize));
	pst->iLvl = in->iLvl;
	if (pst->nLeaves > 1) {
		++in->iLvl;
		mdlReqService(pst->mdl,pst->idUpper,PST_LEVELIZE,in,nIn);
		pstLevelize(pst->pstLower,in,nIn,NULL,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
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


#if (0)
/*
 ** This code segment is no longer considered a proper use of the
 ** services primitives of MDL!
 */
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
#endif


void pstReadTipsy(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inReadTipsy *in = vin;
	int nTotal,nStore;
	char achInFile[PST_FILENAME_SIZE];

	assert(nIn == sizeof(struct inReadTipsy));
	pst->nStart = in->nStart;
	pst->nEnd = in->nEnd;
	nTotal = pst->nEnd - pst->nStart + 1;
	if (pst->nLeaves > 1) {
		pst->nOrdSplit = pst->nStart + pst->nLower*nTotal/pst->nLeaves;
		in->nStart = pst->nOrdSplit;
		mdlReqService(pst->mdl,pst->idUpper,PST_READTIPSY,in,nIn);
		in->nStart = pst->nStart;
		in->nEnd = pst->nOrdSplit - 1;
		pstReadTipsy(pst->pstLower,in,nIn,NULL,NULL);
		in->nEnd = pst->nEnd;
		mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
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
		strcat(achInFile,in->achInFile);
		/*
		 ** Determine the size of the local particle store.
		 */
		nStore = nTotal + (int)ceil(nTotal*in->fExtraStore);
		pkdInitialize(&plcl->pkd,pst->mdl,in->iOrder,nStore,plcl->nPstLvl,
					  in->fPeriod);
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
	int d,ittr,nOut;
	int nLow,nHigh,nLowerStore,nUpperStore;
	float fLow,fHigh;
	float fl,fu,fm;
	struct outFreeStore outFree;
	struct inWeight inWt;
	struct outWeight outWtLow;
	struct outWeight outWtHigh;
	struct inColRejects inCol;
	OREJ *pLowerRej,*pUpperRej;
	int *pidSwap,iRet;

	/*
	 ** First find out how much free storage there is available for particles
	 ** on the lower and upper subset of processors.
	 */
	mdlReqService(pst->mdl,pst->idUpper,PST_FREESTORE,NULL,0);
	pstFreeStore(pst->pstLower,NULL,0,&outFree,NULL);
	nLowerStore = outFree.nFreeStore;
	mdlGetReply(pst->mdl,pst->idUpper,&outFree,NULL);
	nUpperStore = outFree.nFreeStore;
	/*
	 ** Now start the ROOT finder based on balancing weight ALONE!
	 */
	d = iSplitDim;
	fl = pst->bnd.fMin[d];
	fu = pst->bnd.fMax[d];
	fm = (fl + fu)/2;
	ittr = 0;
	while (fl < fm && fm < fu) {
		inWt.iSplitDim = d;
		inWt.fSplit = fm;
		inWt.ittr = ittr;
		inWt.iSplitSide = 1;
		mdlReqService(pst->mdl,pst->idUpper,PST_WEIGHT,&inWt,sizeof(inWt));
		inWt.iSplitSide = 0;
		pstWeight(pst->pstLower,&inWt,sizeof(inWt),&outWtLow,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,&outWtHigh,NULL);
		/*
		 ** Add lower and Upper subsets weights and numbers
		 */
		nLow = outWtLow.nLow + outWtHigh.nLow;
		nHigh = outWtLow.nHigh + outWtHigh.nHigh;
		fLow = outWtLow.fLow + outWtHigh.fLow;
		fHigh = outWtLow.fHigh + outWtHigh.fHigh;
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
			inWt.iSplitDim = d;
			inWt.fSplit = fm;
			inWt.ittr = ittr;
			inWt.iSplitSide = 1;
			mdlReqService(pst->mdl,pst->idUpper,PST_WEIGHT,&inWt,sizeof(inWt));
			inWt.iSplitSide = 0;
			pstWeight(pst->pstLower,&inWt,sizeof(inWt),&outWtLow,NULL);
			mdlGetReply(pst->mdl,pst->idUpper,&outWtHigh,NULL);
			/*
			 ** Add lower and Upper subsets weights and numbers
			 */
			nLow = outWtLow.nLow + outWtHigh.nLow;
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
			inWt.iSplitDim = d;
			inWt.fSplit = fm;
			inWt.ittr = ittr;
			inWt.iSplitSide = 1;
			mdlReqService(pst->mdl,pst->idUpper,PST_WEIGHT,&inWt,sizeof(inWt));
			inWt.iSplitSide = 0;
			pstWeight(pst->pstLower,&inWt,sizeof(inWt),&outWtLow,NULL);
			mdlGetReply(pst->mdl,pst->idUpper,&outWtHigh,NULL);
			/*
			 ** Add lower and Upper subsets weights and numbers
			 */
			nHigh = outWtLow.nHigh + outWtHigh.nHigh;
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
	 **
	 ** Careful, SERVICE PST_COLREJECTS does NOT conform strictly to
	 ** the proper use of MDL. This should be fixed in the future.
	 */
	pLowerRej = malloc(pst->nLower*sizeof(OREJ));
	assert(pLowerRej != NULL);
	pUpperRej = malloc(pst->nUpper*sizeof(OREJ));
	assert(pUpperRej != NULL);
	pidSwap = malloc(mdlThreads(pst->mdl)*sizeof(int));
	assert(pidSwap != NULL);
	inCol.fSplit = pst->fSplit;
	inCol.iSplitDim = pst->iSplitDim;
	inCol.iSplitSide = 1;
	mdlReqService(pst->mdl,pst->idUpper,PST_COLREJECTS,&inCol,sizeof(inCol));
	inCol.iSplitSide = 0;
	pstColRejects(pst->pstLower,&inCol,sizeof(inCol),pLowerRej,&nOut);
	assert(nOut/sizeof(OREJ) == pst->nLower);
	mdlGetReply(pst->mdl,pst->idUpper,pUpperRej,&nOut);
	assert(nOut/sizeof(OREJ) == pst->nUpper);
	while (1) {
		iRet = _pstRejMatch(pst,pst->nLower,pLowerRej,
							pst->nUpper,pUpperRej,pidSwap);
		if (!iRet) break;
		mdlReqService(pst->mdl,pst->idUpper,PST_SWAPREJECTS,pidSwap,
					  mdlThreads(pst->mdl)*sizeof(int));
		pstSwapRejects(pst->pstLower,pidSwap,
					   mdlThreads(pst->mdl)*sizeof(int),pLowerRej,&nOut);
		assert(nOut/sizeof(OREJ) == pst->nLower);
		mdlGetReply(pst->mdl,pst->idUpper,pUpperRej,&nOut);
		assert(nOut/sizeof(OREJ) == pst->nUpper);
		}
	free(pLowerRej);
	free(pUpperRej);
	free(pidSwap);
	}


void pstDomainDecomp(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	int nOut,d,j;
	struct outCalcBound outBnd;
	
	assert(nIn == 0);
	if (pst->nLeaves > 1) {
		/*
		 ** First calculate the Bound for the set.
		 */
		pstCalcBound(pst,NULL,0,&outBnd,&nOut);
		pst->bnd = outBnd.bnd;
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
			mdlReqService(pst->mdl,pst->idUpper,PST_DOMAINDECOMP,NULL,0);
		if (pst->nLower > 1) 
			pstDomainDecomp(pst->pstLower,NULL,0,NULL,NULL);
		if (pst->nUpper > 1) 
			mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
		}
	if (pnOut) *pnOut = 0;
	}


void pstCalcBound(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct outCalcBound *out = vout;
	struct outCalcBound outBnd;
	int j;
	
	assert(nIn == 0);
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_CALCBOUND,NULL,0);
		pstCalcBound(pst->pstLower,NULL,0,out,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,&outBnd,NULL);
		for (j=0;j<3;++j) {
			if (outBnd.bnd.fMin[j] < out->bnd.fMin[j]) {
				out->bnd.fMin[j] = outBnd.bnd.fMin[j];
				}
			if (outBnd.bnd.fMax[j] > out->bnd.fMax[j]) {
				out->bnd.fMax[j] = outBnd.bnd.fMax[j];
				}
			}
		}
	else {
		pkdCalcBound(plcl->pkd,&out->bnd);
		}
	if (pnOut) *pnOut = sizeof(struct outCalcBound); 
	}


void pstWeight(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inWeight *in = vin;
	struct outWeight *out = vout;
	struct outWeight outWt;
	float fSplit,fLow,fHigh;
	int iSplitSide;

	assert(nIn == sizeof(struct inWeight));
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_WEIGHT,in,nIn);
		pstWeight(pst->pstLower,in,nIn,out,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,&outWt,NULL);
		out->nLow += outWt.nLow;
		out->nHigh += outWt.nHigh;
		out->fLow += outWt.fLow;
		out->fHigh += outWt.fHigh;
		}
	else {
		fSplit = in->fSplit;
		iSplitSide = in->iSplitSide;
		if (in->ittr == 0) {
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
		plcl->iPart = pkdWeight(plcl->pkd,in->iSplitDim,fSplit,iSplitSide,
								plcl->iWtFrom,plcl->iWtTo,
								&out->nLow,&out->nHigh,&fLow,&fHigh);
		out->fLow = fLow + plcl->fWtLow;
		out->fHigh = fHigh + plcl->fWtHigh;
		plcl->fLow = fLow;
		plcl->fHigh = fHigh;
		}
	if (pnOut) *pnOut = sizeof(struct outWeight); 
	}


void pstFreeStore(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct outFreeStore *out = vout;
	int nLowerStore,nUpperStore;

	assert(nIn == 0);
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_FREESTORE,NULL,0);
		pstFreeStore(pst->pstLower,NULL,0,out,NULL);
		nLowerStore = out->nFreeStore;
		mdlGetReply(pst->mdl,pst->idUpper,out,NULL);
		nUpperStore = out->nFreeStore;
		out->nFreeStore = nLowerStore + nUpperStore;
		}
	else {
		out->nFreeStore = pkdFreeStore(plcl->pkd);
		}
	if (pnOut) *pnOut = sizeof(struct outFreeStore);
	}


/*
 ** This function does NOT conform to proper use of the MDL services
 ** primitives, as it returns a variable length arguement. This 
 ** cannot be fully supported by MDL. Have to fix this sometime!
 */
void pstColRejects(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inColRejects *in = vin;
	OREJ *pOutRej = vout;
	int nLower,nUpper,iUpper;
	
	assert(nIn == sizeof(struct inColRejects));
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_COLREJECTS,in,nIn);
		pstColRejects(pst->pstLower,in,nIn,&pOutRej[0],&nLower);
		iUpper = nLower/sizeof(OREJ);
		mdlGetReply(pst->mdl,pst->idUpper,&pOutRej[iUpper],&nUpper);
		if (pnOut) *pnOut = nLower + nUpper;
		}
	else {
	    pOutRej->nRejects = pkdColRejects(plcl->pkd,in->iSplitDim,
										  in->fSplit,in->iSplitSide);
		pOutRej->nSpace = pkdSwapSpace(plcl->pkd);
		pOutRej->id = pst->idSelf;
		if (pnOut) *pnOut = sizeof(OREJ);
		}
	}


/*
 ** This function does NOT conform to proper use of the MDL services
 ** primitives, as it returns a variable length arguement. This 
 ** cannot be fully supported by MDL. Have to fix this sometime!
 */
void pstColOrdRejects(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inColOrdRejects *in = vin;
	OREJ *pOutRej = vout;
	int nLower,nUpper,iUpper;
	
	assert(nIn == sizeof(struct inColOrdRejects));
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_COLORDREJECTS,in,nIn);
		pstColOrdRejects(pst->pstLower,in,nIn,&pOutRej[0],&nLower);
		iUpper = nLower/sizeof(OREJ);
		mdlGetReply(pst->mdl,pst->idUpper,&pOutRej[iUpper],&nUpper);
		if (pnOut) *pnOut = nLower + nUpper;
		}
	else {
	    pOutRej->nRejects = pkdColOrdRejects(plcl->pkd,in->nOrdSplit,
											 in->iSplitSide);
		pOutRej->nSpace = pkdSwapSpace(plcl->pkd);
		pOutRej->id = pst->idSelf;
		if (pnOut) *pnOut = sizeof(OREJ);
		}
	}


/*
 ** This function does NOT conform to proper use of the MDL services
 ** primitives, as it uses variable length arguements. This 
 ** cannot be fully supported by MDL. Have to fix this sometime!
 */
void pstSwapRejects(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	int *pidSwap = vin;
	OREJ *pOutRej = vout;
	int nLower,nUpper,iUpper,idSwap;

	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_SWAPREJECTS,vin,nIn);
		pstSwapRejects(pst->pstLower,vin,nIn,&pOutRej[0],&nLower);
		iUpper = nLower/sizeof(OREJ);
		mdlGetReply(pst->mdl,pst->idUpper,&pOutRej[iUpper],&nUpper);
		if (pnOut) *pnOut = nLower + nUpper;
		}
	else {
		idSwap = pidSwap[pst->idSelf];
		pOutRej->nRejects = pkdSwapRejects(plcl->pkd,idSwap);
		pOutRej->nSpace = pkdSwapSpace(plcl->pkd);
		pOutRej->id = pst->idSelf;
		if (pnOut) *pnOut = sizeof(OREJ);
		}
	}


void pstDomainColor(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;

	assert(nIn == 0);
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_DOMAINCOLOR,NULL,0);
		pstDomainColor(pst->pstLower,NULL,0,NULL,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
		}
	else {
		pkdDomainColor(plcl->pkd);
		}
	if (pnOut) *pnOut = 0;
	}


void _pstOrdSplit(PST pst)
{
	struct inColOrdRejects inCol;
	OREJ *pLowerRej,*pUpperRej;
	int *pidSwap,iRet,nOut;

	/*
	 ** First Collect rejects.
	 ** Again NON-SAFE use of MDL here!
	 */
	pLowerRej = malloc(pst->nLower*sizeof(OREJ));
	assert(pLowerRej != NULL);
	pUpperRej = malloc(pst->nUpper*sizeof(OREJ));
	assert(pUpperRej != NULL);
	pidSwap = malloc(mdlThreads(pst->mdl)*sizeof(int));
	assert(pidSwap != NULL);
	inCol.nOrdSplit = pst->nOrdSplit;
	inCol.iSplitSide = 1;
	mdlReqService(pst->mdl,pst->idUpper,PST_COLORDREJECTS,&inCol,
				  sizeof(inCol));
	inCol.iSplitSide = 0;
	pstColOrdRejects(pst->pstLower,&inCol,sizeof(inCol),pLowerRej,&nOut);
	assert(nOut/sizeof(OREJ) == pst->nLower);
	mdlGetReply(pst->mdl,pst->idUpper,pUpperRej,&nOut);
	assert(nOut/sizeof(OREJ) == pst->nUpper);
	while (1) {
		iRet = _pstRejMatch(pst,pst->nLower,pLowerRej,pst->nUpper,
							pUpperRej,pidSwap);
		if (!iRet) break;
		mdlReqService(pst->mdl,pst->idUpper,PST_SWAPREJECTS,pidSwap,
					  mdlThreads(pst->mdl)*sizeof(int));
		pstSwapRejects(pst->pstLower,pidSwap,mdlThreads(pst->mdl)*sizeof(int),
					   pLowerRej,&nOut);
		assert(nOut/sizeof(OREJ) == pst->nLower);
		mdlGetReply(pst->mdl,pst->idUpper,pUpperRej,&nOut);
		assert(nOut/sizeof(OREJ) == pst->nUpper);
		}
	free(pLowerRej);
	free(pUpperRej);
	free(pidSwap);
	}


void pstDomainOrder(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	assert(nIn == 0);
	if (pst->nLeaves > 1) {
		_pstOrdSplit(pst);
		/*
		 ** Now go on to Domain Order of next levels.
		 */
		if (pst->nUpper > 1) 
			mdlReqService(pst->mdl,pst->idUpper,PST_DOMAINORDER,NULL,0);
		if (pst->nLower > 1) 
			pstDomainOrder(pst->pstLower,NULL,0,NULL,NULL);
		if (pst->nUpper > 1) 
			mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
		}
	if (pnOut) *pnOut = 0;
	}


void pstLocalOrder(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;

	assert(nIn == 0);
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_LOCALORDER,NULL,0);
		pstLocalOrder(pst->pstLower,NULL,0,NULL,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
		}
	else {
		pkdLocalOrder(plcl->pkd,pst->nStart);
		}
	if (pnOut) *pnOut = 0;
	}


void pstOutArray(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inOutArray *in = vin;
	char achOutFile[PST_FILENAME_SIZE];

	assert(nIn == sizeof(struct inOutArray));
	if (pst->nLeaves > 1) {
		/*
		 ** Non-Recursive Text output.
		 */
		pstOutArray(pst->pstLower,in,nIn,NULL,NULL);
		mdlReqService(pst->mdl,pst->idUpper,PST_OUTARRAY,in,nIn);
		mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
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
		strcat(achOutFile,in->achOutFile);
		pkdOutArray(plcl->pkd,achOutFile,in->iType);
		}
	if (pnOut) *pnOut = 0;
	}


void pstOutVector(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inOutVector *in = vin;
	char achOutFile[PST_FILENAME_SIZE];

	assert(nIn == sizeof(struct inOutVector));
	if (pst->nLeaves > 1) {
		/*
		 ** Non-Recursive Text output.
		 */
		pstOutVector(pst->pstLower,in,nIn,NULL,NULL);
		mdlReqService(pst->mdl,pst->idUpper,PST_OUTVECTOR,in,nIn);
		mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
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
		strcat(achOutFile,in->achOutFile);
		pkdOutVector(plcl->pkd,achOutFile,in->iDim,in->iType);
		}
	if (pnOut) *pnOut = 0;
	}


void pstWriteTipsy(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inWriteTipsy *in = vin;
	char achOutFile[PST_FILENAME_SIZE];

	assert(nIn == sizeof(struct inWriteTipsy));
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_WRITETIPSY,in,nIn);
		pstWriteTipsy(pst->pstLower,in,nIn,NULL,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
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
		strcat(achOutFile,in->achOutFile);
		pkdWriteTipsy(plcl->pkd,achOutFile,pst->nStart,pst->nEnd,
					  in->bStandard);
		}
	if (pnOut) *pnOut = 0;
	}


void pstSetSoft(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inSetSoft *in = vin;

	assert(nIn == sizeof(struct inSetSoft));
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_SETSOFT,in,nIn);
		pstSetSoft(pst->pstLower,in,nIn,NULL,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
		}
	else {
		pkdSetSoft(plcl->pkd,in->dSoft);
		}
	if (pnOut) *pnOut = 0;
	}


void pstBuildTree(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inBuildTree *in = vin;
	int iCell;
	
	assert(nIn == sizeof(struct inBuildTree));
	iCell = in->iCell;
	if (pst->nLeaves > 1) {
		plcl->kdTop[iCell].iDim = pst->iSplitDim;
		plcl->kdTop[iCell].fSplit = pst->fSplit;
		plcl->kdTop[iCell].bnd = pst->bnd;
		plcl->kdTop[iCell].pLower = pst->idSelf;
		plcl->kdTop[iCell].pUpper = pst->idUpper;
		in->iCell = UPPER(iCell);
		mdlReqService(pst->mdl,pst->idUpper,PST_BUILDTREE,in,nIn);
		in->iCell = LOWER(iCell);
		pstBuildTree(pst->pstLower,in,nIn,NULL,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
		}
	else {
		pkdBuildLocal(plcl->pkd,in->nBucket,in->iOpenType,in->dCrit,
					  &plcl->kdTop[iCell]);
		plcl->kdTop[iCell].pLower = pst->idSelf;
		plcl->kdTop[iCell].pUpper = -1;
		pkdBuildTop(plcl->pkd,in->iOpenType,in->dCrit,plcl->kdTop,
					plcl->nTopNodes,plcl->piLeaf);
		}
	if (pnOut) *pnOut = 0;
	}


void pstDensity(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	struct inDensity *in = vin;

	assert(nIn == sizeof(struct inDensity));
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_DENSITY,in,nIn);
		pstDensity(pst->pstLower,in,nIn,NULL,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
		}
	else {
#ifdef SMOOTH_CODE
		LCL *plcl = pst->plcl;
		SMX smx;

		smInitialize(&smx,plcl->pkd,in->nSmooth,in->bGatherScatter);
		if (in->bGatherScatter) {
			smSmooth(smx,smDensitySym);
			}
		else {
			smSmooth(smx,smDensity);
			}
		smFinish(smx);
#endif
		}
	if (pnOut) *pnOut = 0;
	}


void pstGravity(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	CASTAT csPart,csCell;
	struct inGravity *in = vin;
	struct outGravity *out = vout;
	struct outGravity outUp;

	assert(nIn == sizeof(struct inGravity));
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_GRAVITY,in,nIn);
		pstGravity(pst->pstLower,in,nIn,out,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,&outUp,NULL);
		out->dPartSum += outUp.dPartSum;
		out->dCellSum += outUp.dCellSum;
		out->dWSum += outUp.dWSum;
		out->dISum += outUp.dISum;
		out->dESum += outUp.dESum;
		if (outUp.dWMax > out->dWMax) out->dWMax = outUp.dWMax;
		if (outUp.dIMax > out->dIMax) out->dIMax = outUp.dIMax;
		if (outUp.dEMax > out->dEMax) out->dEMax = outUp.dEMax;
		if (outUp.dWMin < out->dWMin) out->dWMin = outUp.dWMin;
		if (outUp.dIMin < out->dIMin) out->dIMin = outUp.dIMin;
		if (outUp.dEMin < out->dEMin) out->dEMin = outUp.dEMin;
		/*
		 ** Cache statistics sums.
		 */
		out->dpASum += outUp.dpASum;
		out->dpMSum += outUp.dpMSum;
		out->dpCSum += outUp.dpCSum;
		out->dpTSum += outUp.dpTSum;
		out->dcASum += outUp.dcASum;
		out->dcMSum += outUp.dcMSum;
		out->dcCSum += outUp.dcCSum;
		out->dcTSum += outUp.dcTSum;
		}
	else {
		pkdGravAll(plcl->pkd,in->nReps,in->bPeriodic,in->dEwCut,in->dEwhCut,
				   &out->dPartSum,&out->dCellSum,&csPart,&csCell);
		out->dWSum = pkdGetTimer(plcl->pkd,1);
		out->dISum = pkdGetTimer(plcl->pkd,2);
		out->dESum = pkdGetTimer(plcl->pkd,3);
		out->dWMax = out->dWSum;
		out->dIMax = out->dISum;
		out->dEMax = out->dESum;
		out->dWMin = out->dWSum;
		out->dIMin = out->dISum;
		out->dEMin = out->dESum;
		/*
		 ** Cache statistics
		 */
		out->dpASum = csPart.dAccess;
		out->dpMSum = csPart.dMissRatio;
		out->dpCSum = csPart.dCollRatio;
		out->dpTSum = csPart.dMinRatio;
		out->dcASum = csCell.dAccess;
		out->dcMSum = csCell.dMissRatio;
		out->dcCSum = csCell.dCollRatio;
		out->dcTSum = csCell.dMinRatio;
		}
	if (pnOut) *pnOut = sizeof(struct outGravity);
	}


void pstCalcE(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct outCalcE *out = vout;
	struct outCalcE outE;

	assert(nIn == 0);
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_CALCE,NULL,0);
		pstCalcE(pst->pstLower,NULL,0,out,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,&outE,NULL);
		out->T += outE.T;
		out->U += outE.U;
		}
	else {
		pkdCalcE(plcl->pkd,&out->T,&out->U);
		}
	if (pnOut) *pnOut = sizeof(struct outCalcE);
	}


void pstDrift(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inDrift *in = vin;

	assert(nIn == sizeof(struct inDrift));
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_DRIFT,in,nIn);
		pstDrift(pst->pstLower,in,nIn,NULL,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
		}
	else {
		pkdDrift(plcl->pkd,in->dDelta,in->fCenter);
		}
	if (pnOut) *pnOut = 0;
	}


void pstKick(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inKick *in = vin;

	assert(nIn == sizeof(struct inKick));
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_KICK,in,nIn);
		pstKick(pst->pstLower,in,nIn,NULL,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
		}
	else {
		pkdKick(plcl->pkd,in->dvFacOne,in->dvFacTwo);
		}
	if (pnOut) *pnOut = 0;
	}


void pstReadCheck(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inReadCheck *in = vin;
	char achInFile[PST_FILENAME_SIZE];
	int nTotal,nStore;

	assert(nIn == sizeof(struct inReadCheck));
	pst->nStart = in->nStart;
	pst->nEnd = in->nEnd;
	nTotal = pst->nEnd - pst->nStart + 1;
	if (pst->nLeaves > 1) {
		pst->nOrdSplit = pst->nStart + pst->nLower*nTotal/pst->nLeaves;
		in->nStart = pst->nOrdSplit;
		mdlReqService(pst->mdl,pst->idUpper,PST_READCHECK,in,nIn);
		in->nStart = pst->nStart;
		in->nEnd = pst->nOrdSplit - 1;
		pstReadCheck(pst->pstLower,in,nIn,NULL,NULL);
		in->nEnd = pst->nEnd;
		mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
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
		strcat(achInFile,in->achInFile);
		/*
		 ** Determine the size of the local particle store.
		 */
		nStore = nTotal + (int)ceil(nTotal*in->fExtraStore);
		pkdInitialize(&plcl->pkd,pst->mdl,in->iOrder,nStore,plcl->nPstLvl,
					  in->fPeriod);
		if (in->bNewCheck) {
			pkdReadCheckNew(plcl->pkd,achInFile,pst->nStart,nTotal);
			}	
		else {
			pkdReadCheckOld(plcl->pkd,achInFile,pst->nStart,nTotal);
			}
		}
	if (pnOut) *pnOut = 0;
	}


void pstSetTotal(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{	
	LCL *plcl = pst->plcl;
	struct outSetTotal *out = vout;
	struct outSetTotal oute;

	assert(nIn == 0);
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_SETTOTAL,NULL,0);
		pstSetTotal(pst->pstLower,NULL,0,out,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,&oute,NULL);
		out->nTotal += oute.nTotal;
		pst->nTotal = out->nTotal;
		}
	else {
		pst->nTotal = pkdLocal(plcl->pkd);
		out->nTotal = pst->nTotal;
		}
	if (pnOut) *pnOut = sizeof(struct outSetTotal);
	}


void pstWriteCheck(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inWriteCheck *in = vin;
	char achOutFile[PST_FILENAME_SIZE];
	int nStart;

	assert(nIn == sizeof(struct inWriteCheck));
	if (pst->nLeaves > 1) {
		nStart = in->nStart;
		in->nStart += pst->pstLower->nTotal;
		mdlReqService(pst->mdl,pst->idUpper,PST_WRITECHECK,in,nIn);
		in->nStart = nStart;
		pstWriteCheck(pst->pstLower,in,nIn,NULL,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
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
		strcat(achOutFile,in->achOutFile);
		pkdWriteCheckNew(plcl->pkd,achOutFile,in->nStart);
		}
	if (pnOut) *pnOut = 0;
	}






