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


void pstAddServices(PST pst,MDL mdl)
{
	int nThreads,nCell;

	nThreads = mdlThreads(mdl);
	mdlAddService(mdl,PST_SETADD,pst,pstSetAdd,
				  sizeof(struct inSetAdd),0);
	mdlAddService(mdl,PST_LEVELIZE,pst,pstLevelize,
				  sizeof(struct inLevelize),0);
	mdlAddService(mdl,PST_READTIPSY,pst,pstReadTipsy,
				  sizeof(struct inReadTipsy),0);
	mdlAddService(mdl,PST_DOMAINDECOMP,pst,pstDomainDecomp,
				  0,0);
	mdlAddService(mdl,PST_CALCBOUND,pst,pstCalcBound,
				  0,sizeof(struct outCalcBound));
	mdlAddService(mdl,PST_WEIGHT,pst,pstWeight,
				  sizeof(struct inWeight),sizeof(struct outWeight));
	mdlAddService(mdl,PST_FREESTORE,pst,pstFreeStore,
				  0,sizeof(struct outFreeStore));
	mdlAddService(mdl,PST_COLREJECTS,pst,pstColRejects,
				  sizeof(struct inColRejects),nThreads*sizeof(OREJ));
	mdlAddService(mdl,PST_SWAPREJECTS,pst,pstSwapRejects,
				  nThreads*sizeof(int),nThreads*sizeof(OREJ));
	mdlAddService(mdl,PST_DOMAINCOLOR,pst,pstDomainColor,
				  0,0);
	mdlAddService(mdl,PST_COLORDREJECTS,pst,pstColOrdRejects,
				  sizeof(struct inColOrdRejects),nThreads*sizeof(OREJ));
	mdlAddService(mdl,PST_DOMAINORDER,pst,pstDomainOrder,
				  0,0);
	mdlAddService(mdl,PST_LOCALORDER,pst,pstLocalOrder,
				  0,0);
	mdlAddService(mdl,PST_OUTARRAY,pst,pstOutArray,
				  sizeof(struct inOutArray),0);
	mdlAddService(mdl,PST_OUTVECTOR,pst,pstOutVector,
				  sizeof(struct inOutVector),0);
	mdlAddService(mdl,PST_WRITETIPSY,pst,pstWriteTipsy,
				  sizeof(struct inWriteTipsy),0);
	mdlAddService(mdl,PST_BUILDTREE,pst,pstBuildTree,
				  sizeof(struct inBuildTree),sizeof(struct outBuildTree));
	mdlAddService(mdl,PST_DENSITY,pst,pstDensity,
				  sizeof(struct inDensity),0);
	mdlAddService(mdl,PST_GRAVITY,pst,pstGravity,
				  sizeof(struct inGravity),sizeof(struct outGravity));
	mdlAddService(mdl,PST_CALCE,pst,pstCalcE,
				  0,sizeof(struct outCalcE));
	mdlAddService(mdl,PST_DRIFT,pst,pstDrift,
				  sizeof(struct inDrift),0);
	mdlAddService(mdl,PST_KICK,pst,pstKick,
				  sizeof(struct inKick),0);
	mdlAddService(mdl,PST_READCHECK,pst,pstReadCheck,
				  sizeof(struct inReadCheck),0);
	mdlAddService(mdl,PST_WRITECHECK,pst,pstWriteCheck,
				  sizeof(struct inWriteCheck),0);
	mdlAddService(mdl,PST_SETSOFT,pst,pstSetSoft,
				  sizeof(struct inSetSoft),0);
	mdlAddService(mdl,PST_SETTOTAL,pst,pstSetTotal,
				  0,sizeof(struct outSetTotal));
	mdlAddService(mdl,PST_CALCCELL,pst,pstCalcCell,
				  sizeof(struct inCalcCell),sizeof(struct outCalcCell));
	/*
	 ** Calculate the number of levels in the top tree and use it to 
	 ** define the size of the messages.
	 */
	nCell = 1<<(1+(int)ceil(log((double)nThreads)/log(2.0)));
	mdlAddService(mdl,PST_COLCELLS,pst,pstColCells,
				  sizeof(struct inColCells),nCell*sizeof(KDN));
	mdlAddService(mdl,PST_DISTRIBCELLS,pst,pstDistribCells,
				  nCell*sizeof(KDN),0);
	mdlAddService(mdl,PST_CALCROOT,pst,pstCalcRoot,
				  0,sizeof(struct ioCalcRoot));
	mdlAddService(mdl,PST_DISTRIBROOT,pst,pstDistribRoot,
				  sizeof(struct ioCalcRoot),0);
	mdlAddService(mdl,PST_ONENODEREADINIT,pst,pstOneNodeReadInit,
		      sizeof(struct inReadTipsy), nThreads*sizeof(int));
	mdlAddService(mdl,PST_SWAPALL,pst,pstSwapAll,
				  sizeof(int),0);
	mdlAddService(mdl,PST_MASSCHECK,pst,
		      (void (*)(void *,void *,int,void *,int *))pstMassCheck,
				  0,sizeof(struct outMassCheck));
	mdlAddService(mdl,PST_ACTIVEORDER,pst,
		      (void (*)(void *,void *,int,void *,int *))pstActiveOrder,
		      0,0);
	mdlAddService(mdl,PST_SETRUNG,pst,
		      (void (*)(void *,void *,int,void *,int *))pstSetRung,
		      sizeof(struct inSetRung), 0);
	mdlAddService(mdl,PST_ACTIVERUNG,pst,
		      (void (*)(void *,void *,int,void *,int *))pstActiveRung,
		      sizeof(struct inActiveRung), 0);
	mdlAddService(mdl,PST_CURRRUNG,pst,
		      (void (*)(void *,void *,int,void *,int *))pstCurrRung,
		      sizeof(struct inCurrRung), sizeof(struct outCurrRung));
	mdlAddService(mdl,PST_DENSITYRUNG,pst,
		      (void (*)(void *,void *,int,void *,int *))pstDensityRung,
		      sizeof(struct inDensityRung),
		      sizeof(struct outDensityRung));
	mdlAddService(mdl,PST_RUNGSTATS,pst,pstRungStats,
				  sizeof(struct inRungStats),
				  sizeof(struct outRungStats));
	mdlAddService(mdl,PST_GETMAP,pst,pstGetMap,
				  sizeof(struct inGetMap),nThreads*sizeof(int));
	mdlAddService(mdl,PST_VELOCITYRUNG,pst,
		      (void (*)(void *,void *,int,void *,int *))pstVelocityRung,
		      sizeof(struct inVelocityRung),
		      sizeof(struct outVelocityRung));
	}


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
		}
	if (pnOut) *pnOut = 0;
	}

void pstGetMap(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	struct inGetMap *in = vin;
	int *out = vout;
	int *tmp,i;

	assert(nIn == sizeof(struct inGetMap));
	if (pst->nLeaves > 1) {
		tmp = malloc(mdlThreads(pst->mdl)*sizeof(int));
		assert(tmp != NULL);
		pstGetMap(pst->pstLower,in,nIn,vout,pnOut);
		in->nStart += pst->nLower;
		mdlReqService(pst->mdl,pst->idUpper,PST_GETMAP,in,nIn);
		mdlGetReply(pst->mdl,pst->idUpper,tmp,pnOut);
		for (i=0;i<pst->nUpper;++i) {
			out[in->nStart+i] = tmp[in->nStart+i];
			}
		in->nStart -= pst->nLower;
		free(tmp);
		}
	else {
		out[in->nStart] = pst->idSelf;
		}
	if (pnOut) *pnOut = mdlThreads(pst->mdl)*sizeof(int);
	}

void pstOneNodeReadInit(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inReadTipsy *in = vin;
	int *pout = vout;
	int nTotal,nStore;
	int *ptmp;
	int i;

	assert(nIn == sizeof(struct inReadTipsy));
	pst->nStart = in->nStart;
	pst->nEnd = in->nEnd;
	nTotal = pst->nEnd - pst->nStart + 1;
	if (pst->nLeaves > 1) {
		pst->nOrdSplit = pst->nStart
		    + pst->nLower*(nTotal/pst->nLeaves);
		in->nStart = pst->nOrdSplit;
		mdlReqService(pst->mdl,pst->idUpper,PST_ONENODEREADINIT,
			      in, nIn);
		in->nStart = pst->nStart;
		in->nEnd = pst->nOrdSplit - 1;
		pstOneNodeReadInit(pst->pstLower,in,nIn,vout,pnOut);
		in->nEnd = pst->nEnd;
		ptmp = malloc(pst->mdl->nThreads*sizeof(*ptmp));
		assert(ptmp != NULL);
		mdlGetReply(pst->mdl,pst->idUpper,ptmp,pnOut);
		for(i = 0; i < pst->mdl->nThreads; i++) {
		    if(ptmp[i] != -1)
			pout[i] = ptmp[i];
		    }
		free(ptmp);
		}
	else {
		for(i = 0; i < pst->mdl->nThreads; i++)
		    pout[i] = -1;
		/*
		 ** Determine the size of the local particle store.
		 */
		nStore = nTotal + (int)ceil(nTotal*in->fExtraStore);
		
		pkdInitialize(&plcl->pkd,pst->mdl,in->iOrder,nStore,
		    plcl->nPstLvl, in->fPeriod);
		pout[pst->idSelf] = nTotal;
		}
	if (pnOut) *pnOut = pst->mdl->nThreads*sizeof(*pout);
	}

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
		pst->nOrdSplit = pst->nStart
		    + pst->nLower*(nTotal/pst->nLeaves);
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
		pkdReadTipsy(plcl->pkd,achInFile,pst->nStart,nTotal,in->dvFac);
		}
	if (pnOut) *pnOut = 0;
	}


int _pstRejMatch(PST pst,int n1,OREJ *p1,int n2,OREJ *p2,int *pidSwap)
{
	int id,i,i1,i2,nLarge,id1,id2;
	int s1,s2,r1,r2;

	/*
	 ** Check to see if there is enough space...
	 */
	s1 = 0;
    r1 = 0;
	for (i=0;i<n1;++i) {
		s1 += p1[i].nSpace;
		r1 += p1[i].nRejects;
		}
	s2 = 0;
	r2 = 0;
	for (i=0;i<n2;++i) {
		s2 += p2[i].nSpace;
		r2 += p2[i].nRejects;
		}
	assert(r1 <= s2);
	assert(r2 <= s1);
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


#define NUM_SAFETY	2

void _pstRootSplit(PST pst,int iSplitDim,double dMass)
{
	int d,ittr,nOut;
	int nLow,nHigh,nLowerStore,nUpperStore;
	float fLow,fHigh;
	float fl,fu,fm,fmm;
	struct outFreeStore outFree;
	struct inWeight inWt;
	struct outWeight outWtLow;
	struct outWeight outWtHigh;
	struct inColRejects inCol;
	OREJ *pLowerRej,*pUpperRej;
	int *pidSwap,iRet;
	int nLowTot,nHighTot;

	struct outMassCheck outMass;
#ifdef PARANOID_CHECK
	int i,iLowSum,iHighSum;
#endif

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
	 ** Make sure that the particles are ordered into active-inactive order.
	 */
	mdlReqService(pst->mdl,pst->idUpper,PST_ACTIVEORDER,NULL,0);
	pstActiveOrder(pst->pstLower,NULL,0,NULL,NULL);
	mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
	/*
	 ** Now start the ROOT finder based on balancing active weight ALONE!
	 */
	d = iSplitDim;
	fl = pst->bnd.fMin[d];
	fu = pst->bnd.fMax[d];
	fmm = (fl + fu)/2;
	ittr = 0;
	while (fl < fmm && fmm < fu) {
		fm = fmm;
		inWt.iSplitDim = d;
		inWt.fSplit = fm;
		inWt.ittr = ittr;
		inWt.iSplitSide = 1;
		inWt.pFlag = 1;
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

		fmm = (fl + fu)/2;
		++ittr;
		}
	/*
	 ** If we exceed the local pStore in the lower or upper subsets then
	 ** we have to relax the condition that work be balanced and try to
	 ** max out the number of particles in the subset which had too many.
	 */
	if (nLow > nLowerStore-NUM_SAFETY) {
		fl = pst->bnd.fMin[d];
		fu = fm;
		fmm = (fl + fu)/2;
		ittr = 0;
	    while (fl < fmm && fmm < fu) {
			fm = fmm;
			inWt.iSplitDim = d;
			inWt.fSplit = fm;
			inWt.ittr = ittr;
			inWt.iSplitSide = 1;
			inWt.pFlag = 1;
			mdlReqService(pst->mdl,pst->idUpper,PST_WEIGHT,&inWt,sizeof(inWt));
			inWt.iSplitSide = 0;
			pstWeight(pst->pstLower,&inWt,sizeof(inWt),&outWtLow,NULL);
			mdlGetReply(pst->mdl,pst->idUpper,&outWtHigh,NULL);
			/*
			 ** Add lower and Upper subsets number of active particles
			 */
			nLow = outWtLow.nLow + outWtHigh.nLow;
			nHigh = outWtLow.nHigh + outWtHigh.nHigh;
			/*
			printf("Fit ittr:%d l:%d\n",ittr,nLow);
			*/
			if (nLow > nLowerStore) fu = fm;
			else if (nLow < nLowerStore) fl = fm;
			else {
				fl = fm;
				break;
				}
			fmm = (fl + fu)/2;
			++ittr;
			}
		/*
		printf("Fit ittr:%d l:%d <= %d\n",ittr,nLow,nLowerStore);
		*/
		assert(nLow <= nLowerStore);
		}
	else if (nHigh > nUpperStore-NUM_SAFETY) {
		fl = fm;
		fu = pst->bnd.fMax[d];
		fmm = (fl + fu)/2;
		ittr = 0;
	    while (fl < fmm && fmm < fu) {
			fm = fmm;
			inWt.iSplitDim = d;
			inWt.fSplit = fm;
			inWt.ittr = ittr;
			inWt.iSplitSide = 1;
			inWt.pFlag = 1;
			mdlReqService(pst->mdl,pst->idUpper,PST_WEIGHT,&inWt,sizeof(inWt));
			inWt.iSplitSide = 0;
			pstWeight(pst->pstLower,&inWt,sizeof(inWt),&outWtLow,NULL);
			mdlGetReply(pst->mdl,pst->idUpper,&outWtHigh,NULL);
			/*
			 ** Add lower and Upper subsets number of active particles
			 */
			nLow = outWtLow.nLow + outWtHigh.nLow;
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
			fmm = (fl + fu)/2;
			++ittr;
			}
		/*
		printf("Fit ittr:%d u:%d <= %d\n",ittr,nHigh,nUpperStore);
		*/
		assert(nHigh <= nUpperStore);
		}
	pst->fSplit = fm;
	pst->iSplitDim = d;
	/*
	 ** Now we see if the TOTAL number of particles in the lower and upper
	 ** subsets exceeds the local pStores. If so then we need to find a new
	 ** boundary to distribute the INACTIVE particles so that everything 
	 ** fits.
	 */
	inWt.iSplitDim = d;
	inWt.fSplit = fm;
	inWt.ittr = 0;
	inWt.iSplitSide = 1;
	inWt.pFlag = -1;
	mdlReqService(pst->mdl,pst->idUpper,PST_WEIGHT,&inWt,sizeof(inWt));
	inWt.iSplitSide = 0;
	pstWeight(pst->pstLower,&inWt,sizeof(inWt),&outWtLow,NULL);
	mdlGetReply(pst->mdl,pst->idUpper,&outWtHigh,NULL);
	/*
	 ** Add lower and Upper subsets numbers of particles
	 */
	nLowTot = nLow + outWtLow.nLow + outWtHigh.nLow;
	nHighTot = nHigh + outWtLow.nHigh + outWtHigh.nHigh;
	printf("Fit stats: Active: %d %d Inactive: %d %d Sum: %d %d\n", nLow, nHigh, outWtLow.nLow + outWtHigh.nLow, outWtLow.nHigh +
	       outWtHigh.nHigh, nLowTot, nHighTot);
	if (nLowTot > nLowerStore-NUM_SAFETY) {
		fl = pst->bnd.fMin[d];
		fu = fm;
		fmm = (fl + fu)/2;
		ittr = 1;
	    while (fl < fmm && fmm < fu) {
			fm = fmm;
			inWt.iSplitDim = d;
			inWt.fSplit = fm;
			inWt.ittr = ittr;
			inWt.iSplitSide = 1;
			inWt.pFlag = -1;
			mdlReqService(pst->mdl,pst->idUpper,PST_WEIGHT,&inWt,sizeof(inWt));
			inWt.iSplitSide = 0;
			pstWeight(pst->pstLower,&inWt,sizeof(inWt),&outWtLow,NULL);
			mdlGetReply(pst->mdl,pst->idUpper,&outWtHigh,NULL);
			/*
			 ** Add lower and Upper subsets numbers of particles
			 */
			nLowTot = nLow + outWtLow.nLow + outWtHigh.nLow;
			nHighTot = nHigh + outWtLow.nHigh + outWtHigh.nHigh;
			/*
			printf("Inactive Fit ittr:%d l:%d\n",ittr,nLowTot);
			*/
			if (nLowTot > nLowerStore) fu = fm;
			else if (nLowTot < nLowerStore) fl = fm;
			else {
				fl = fm;
				break;
				}
			fmm = (fl + fu)/2;
			++ittr;
			}
		/*
		printf("Inactive Fit ittr:%d l:%d <= %d\n",ittr,nLowTot,nLowerStore);
		*/
		assert(nLowTot <= nLowerStore);
		}
	else if (nHighTot > nUpperStore-NUM_SAFETY) {
		fl = fm;
		fu = pst->bnd.fMax[d];
		fmm = (fl + fu)/2;
		ittr = 1;
	    while (fl < fmm && fmm < fu) {
			fm = fmm;
			inWt.iSplitDim = d;
			inWt.fSplit = fm;
			inWt.ittr = ittr;
			inWt.iSplitSide = 1;
			inWt.pFlag = -1;
			mdlReqService(pst->mdl,pst->idUpper,PST_WEIGHT,&inWt,sizeof(inWt));
			inWt.iSplitSide = 0;
			pstWeight(pst->pstLower,&inWt,sizeof(inWt),&outWtLow,NULL);
			mdlGetReply(pst->mdl,pst->idUpper,&outWtHigh,NULL);
			/*
			 ** Add lower and Upper subsets numbers of particles
			 */
			nLowTot = nLow + outWtLow.nLow + outWtHigh.nLow;
			nHighTot = nHigh + outWtLow.nHigh + outWtHigh.nHigh;
			/*
			printf("Inactive Fit ittr:%d u:%d\n",ittr,nHighTot);
			*/
			if (nHighTot > nUpperStore) fl = fm;
			else if (nHighTot < nUpperStore) fu = fm;
			else {
				fu = fm;
				break;
				}
			fmm = (fl + fu)/2;
			++ittr;
			}
		/*
		printf("Inactive Fit ittr:%d u:%d <= %d\n",ittr,nHighTot,nUpperStore);
		*/
		assert(nHighTot <= nUpperStore);
		}
	pst->fSplitInactive = fm;
	/*
	 ** First Collect rejects.
	 **
	 ** Careful, SERVICE PST_COLREJECTS does NOT conform strictly to
	 ** the proper use of MDL. This should be fixed in the future.
	 ** FIXED -- MDL modified.
	 */
	pLowerRej = malloc(pst->nLower*sizeof(OREJ));
	assert(pLowerRej != NULL);
	pUpperRej = malloc(pst->nUpper*sizeof(OREJ));
	assert(pUpperRej != NULL);
	pidSwap = malloc(mdlThreads(pst->mdl)*sizeof(int));
	assert(pidSwap != NULL);

	inCol.fSplit = pst->fSplit;
	inCol.fSplitInactive = pst->fSplitInactive;
	inCol.iSplitDim = pst->iSplitDim;
	inCol.iSplitSide = 1;
	mdlReqService(pst->mdl,pst->idUpper,PST_COLREJECTS,&inCol,sizeof(inCol));
	inCol.iSplitSide = 0;
	pstColRejects(pst->pstLower,&inCol,sizeof(inCol),pLowerRej,&nOut);
	assert(nOut/sizeof(OREJ) == pst->nLower);
	mdlGetReply(pst->mdl,pst->idUpper,pUpperRej,&nOut);
	assert(nOut/sizeof(OREJ) == pst->nUpper);

#ifdef PARANOID_CHECK
	/*
	 ** This paranoid check no longer works for Active particles, need
	 ** to modify or remove this!
	 */
	iLowSum = 0;
	iHighSum = 0;
	for (i=0;i<pst->nLower;++i) {
		iLowSum += pLowerRej[i].nLocal;
		iHighSum += pLowerRej[i].nRejects;
		}
	for (i=0;i<pst->nUpper;++i) {
		iHighSum += pUpperRej[i].nLocal;
		iLowSum += pUpperRej[i].nRejects;
		}
	printf("%d l:%d nLow:%d == %d, nHigh:%d == %d\n",pst->idSelf,pst->iLvl,
		   nLow,iLowSum,nHigh,iHighSum);
	iLowSum = 0;
	iHighSum = 0;
	for (i=0;i<pst->nLower;++i) {
		iLowSum += pLowerRej[i].nLocal;
		iHighSum += pLowerRej[i].nRejects;
		}
	printf("%d l:%d nLow:%d == %d, nHigh:%d == %d\n",pst->idSelf,pst->iLvl,
		   outWtLow.nLow,iLowSum,outWtLow.nHigh,iHighSum);
	iLowSum = 0;
	iHighSum = 0;
	for (i=0;i<pst->nUpper;++i) {
		iHighSum += pUpperRej[i].nLocal;
		iLowSum += pUpperRej[i].nRejects;
		}
	printf("%d l:%d nLow:%d == %d, nHigh:%d == %d\n",pst->idSelf,pst->iLvl,
		   outWtHigh.nLow,iLowSum,outWtHigh.nHigh,iHighSum);
#endif
	
	pstMassCheck(pst,NULL,0,&outMass,NULL);
	if (dMass != outMass.dMass) {
		printf("ERROR id:%d lvl:%d:in _pstRootSplit after ColRejects\n",
			   pst->idSelf,pst->iLvl);
		}

	ittr = 0;
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

		pstMassCheck(pst,NULL,0,&outMass,NULL);
		if (dMass != outMass.dMass) {
			printf("ERROR id:%d lvl:%d iter:%d in _pstRootSplit after Swap\n",
				   pst->idSelf,pst->iLvl,ittr);
			}
		++ittr;
		}
	free(pLowerRej);
	free(pUpperRej);
	free(pidSwap);
	}


void pstDomainDecomp(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	int nOut,d,j;
	struct outCalcBound outBnd;
	struct outMassCheck outMass;
	double dMass;

	assert(nIn == 0);
	if (pst->nLeaves > 1) {
		/*
		 ** First calculate the Bounds for the set.
		 */
		pstCalcBound(pst,NULL,0,&outBnd,&nOut);
		pst->bndActive = outBnd.bndActive;
		pst->bnd = outBnd.bnd;
		/*
		 ** Next determine the longest axis based on the active bounds.
		 */
		d = 0;
		for (j=1;j<3;++j) {
			if (pst->bndActive.fMax[j]-pst->bndActive.fMin[j] > 
				pst->bndActive.fMax[d]-pst->bndActive.fMin[d]) d = j;
			}
		pst->iSplitDim = d;
		
		pstMassCheck(pst,NULL,0,&outMass,NULL);
		dMass = outMass.dMass;
		_pstRootSplit(pst,d,dMass);
		pstMassCheck(pst,NULL,0,&outMass,NULL);
		if (dMass != outMass.dMass) {
			printf("ERROR id:%d lvl:%d:_pstRootSplit mass not cons.\n",
				   pst->idSelf,pst->iLvl);
			}
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
			if (outBnd.bndActive.fMin[j] < out->bndActive.fMin[j]) {
				out->bndActive.fMin[j] = outBnd.bndActive.fMin[j];
				}
			if (outBnd.bndActive.fMax[j] > out->bndActive.fMax[j]) {
				out->bndActive.fMax[j] = outBnd.bndActive.fMax[j];
				}
			}
		}
	else {
		pkdCalcBound(plcl->pkd,&out->bnd,&out->bndActive);
		}
	if (pnOut) *pnOut = sizeof(struct outCalcBound); 
	}


/*
 ** Make sure that the local particles are split into active and inactive
 ** when passing pFlag != 0.
 */
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
			if (in->pFlag == 0) {
				plcl->iWtFrom = 0;
				plcl->iWtTo = pkdLocal(plcl->pkd)-1;
				}
			else if (in->pFlag > 0) {
				/*
				 ** Particles must be in the active-inactive order here!
				 */
				plcl->iWtFrom = 0;
				plcl->iWtTo = pkdActive(plcl->pkd)-1;
				}
			else {
				/*
				 ** Particles must be in the active-inactive order here!
				 */
				plcl->iWtFrom = pkdActive(plcl->pkd);
				plcl->iWtTo = pkdLocal(plcl->pkd)-1;
				}
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
		if(in->pFlag > 0) {
		    if(iSplitSide) out->nLow -= pkdInactive(plcl->pkd);
		    else out->nHigh -= pkdInactive(plcl->pkd);
		    }
		if(in->pFlag < 0) {
		    if(iSplitSide) out->nHigh -= pkdActive(plcl->pkd);
		    else out->nLow -= pkdActive(plcl->pkd);
		    }
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
										  in->fSplit,in->fSplitInactive,
										  in->iSplitSide);
		pOutRej->nSpace = pkdSwapSpace(plcl->pkd);
		pOutRej->id = pst->idSelf;
		pOutRej->nLocal = pkdLocal(plcl->pkd);

		if (pnOut) *pnOut = sizeof(OREJ);
		}
	}


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

/*
 * Routine to swap all particles.  Note that this does not walk the
 * but simply works with one other processor.
 */
void pstSwapAll(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl;
	int *pidSwap = vin;
	PST lpst;

	assert(nIn == sizeof(*pidSwap));
	lpst = pst;
	while(lpst->nLeaves > 1)
	    lpst = lpst->pstLower;

	plcl = lpst->plcl;
	pkdSwapAll(plcl->pkd, *pidSwap);
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


void pstActiveOrder(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;

	assert(nIn == 0);
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_ACTIVEORDER,NULL,0);
		pstActiveOrder(pst->pstLower,NULL,0,NULL,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
		}
	else {
		pkdActiveOrder(plcl->pkd);
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
					  in->bStandard,in->dvFac);
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
	struct outBuildTree *out = vout;
	struct outBuildTree out1,out2;
	struct inCalcCell inc;
	struct outCalcCell outc;
	int j;
	float dOpen;
	
	assert(nIn == sizeof(struct inBuildTree));
	if (pst->nLeaves > 1) {
		pst->kdn.iDim = pst->iSplitDim;
		pst->kdn.fSplit = pst->fSplit;
		pst->kdn.pLower = -1;
		pst->kdn.pUpper = 1;
		mdlReqService(pst->mdl,pst->idUpper,PST_BUILDTREE,in,nIn);
		pstBuildTree(pst->pstLower,in,nIn,&out1,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,&out2,NULL);
		/*
		 ** Combine to find cell mass, cm, softening.
		 ** This also combine the bounds to set pst->kdn.bnd!
		 */
		pkdCombine(&out1.kdn,&out2.kdn,&pst->kdn);
		/*
		 ** Calculate moments and other center-of-mass related quantities.
		 */
		for (j=0;j<3;++j) inc.rcm[j] = pst->kdn.r[j];
		inc.iOrder = in->iOrder;
		pstCalcCell(pst,&inc,sizeof(struct inCalcCell),&outc,NULL);
		pst->kdn.mom = outc.mom;
		/*
		 ** Calculate an opening radius.
		 */
		dOpen = pkdCalcOpen(&pst->kdn,in->iOpenType,in->dCrit,in->iOrder);
		pst->kdn.fOpen2 = dOpen*dOpen;
		}
	else {
		if (in->bBinary) {
			pkdBuildBinary(plcl->pkd,in->nBucket,in->iOpenType,in->dCrit,
						   in->iOrder,in->bActiveOnly,&pst->kdn);
			}
		else {
			pkdBuildLocal(plcl->pkd,in->nBucket,in->iOpenType,in->dCrit,
						  in->iOrder,in->bActiveOnly,&pst->kdn);
			}
		pst->kdn.pLower = pst->idSelf;
		pst->kdn.pUpper = 1;
		}
	/*
	 ** Calculated all cell properties, now pass up this cell info.
	 */
	out->kdn = pst->kdn;
	if (pnOut) *pnOut = sizeof(struct outBuildTree);
	}


void pstDensity(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	struct inDensity *in = vin;

	fprintf(stderr,"In pstDensity...\n");
	assert(nIn == sizeof(struct inDensity));
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_DENSITY,in,nIn);
		pstDensity(pst->pstLower,in,nIn,NULL,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
		}
	else {
		LCL *plcl = pst->plcl;
		SMX smx;

		smInitialize(&smx,plcl->pkd,in->nSmooth,in->bGatherScatter,
			     in->bPeriodic);
		if (in->bGatherScatter) {
                        fprintf(stderr,"In GatherScatter\n");
			smSmooth(smx,smDensitySym);
			}
		else {
                        fprintf(stderr,"In no GatherScatter\n");
			smSmooth(smx,smDensity);
			}
                fprintf(stderr,"Before smFinish\n");
		smFinish(smx);
                fprintf(stderr,"After smFinish\n");
		}
	if (pnOut) *pnOut = 0;
	}


void pstGravity(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	CASTAT cs;
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
		out->nActive += outUp.nActive;
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
		pkdGravAll(plcl->pkd,in->nReps,in->bPeriodic,in->iOrder,in->iEwOrder,
				   in->dEwCut,in->dEwhCut,&out->nActive,
				   &out->dPartSum,&out->dCellSum,&cs);
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
		out->dpASum = cs.dpNumAccess;
		out->dpMSum = cs.dpMissRatio;
		out->dpCSum = cs.dpCollRatio;
		out->dpTSum = cs.dpMinRatio;
		out->dcASum = cs.dcNumAccess;
		out->dcMSum = cs.dcMissRatio;
		out->dcCSum = cs.dcCollRatio;
		out->dcTSum = cs.dcMinRatio;
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
		pkdDrift(plcl->pkd,in->dDelta,in->fCenter,in->bPeriodic);
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
		pst->nOrdSplit = pst->nStart
		    + pst->nLower*(nTotal/pst->nLeaves);
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
		pkdReadCheck(plcl->pkd,achInFile,in->iVersion,in->iOffset,
					 pst->nStart,nTotal);
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
		pkdWriteCheck(plcl->pkd,achOutFile,in->iOffset,in->nStart);
		}
	if (pnOut) *pnOut = 0;
	}


void pstCalcCell(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inCalcCell *in = vin;
	struct outCalcCell *out = vout;
	struct outCalcCell occ;
	struct pkdCalcCellStruct pcc;

	assert(nIn == sizeof(struct inCalcCell));
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_CALCCELL,in,nIn);
		pstCalcCell(pst->pstLower,in,nIn,out,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,&occ,NULL);
		switch (in->iOrder) {
		case 4:
			out->mom.Hxxxx += occ.mom.Hxxxx;
			out->mom.Hxyyy += occ.mom.Hxyyy;
			out->mom.Hxxxy += occ.mom.Hxxxy;
			out->mom.Hyyyy += occ.mom.Hyyyy;
			out->mom.Hxxxz += occ.mom.Hxxxz;
			out->mom.Hyyyz += occ.mom.Hyyyz;
			out->mom.Hxxyy += occ.mom.Hxxyy;
			out->mom.Hxxyz += occ.mom.Hxxyz;
			out->mom.Hxyyz += occ.mom.Hxyyz;
			out->mom.Hxxzz += occ.mom.Hxxzz;
			out->mom.Hxyzz += occ.mom.Hxyzz;
			out->mom.Hxzzz += occ.mom.Hxzzz;
			out->mom.Hyyzz += occ.mom.Hyyzz;
			out->mom.Hyzzz += occ.mom.Hyzzz;
			out->mom.Hzzzz += occ.mom.Hzzzz;
			out->mom.B6 += occ.mom.B6;
		case 3:
			out->mom.Oxxx += occ.mom.Oxxx;
			out->mom.Oxyy += occ.mom.Oxyy;
			out->mom.Oxxy += occ.mom.Oxxy;
			out->mom.Oyyy += occ.mom.Oyyy;
			out->mom.Oxxz += occ.mom.Oxxz;
			out->mom.Oyyz += occ.mom.Oyyz;
			out->mom.Oxyz += occ.mom.Oxyz;
			out->mom.Oxzz += occ.mom.Oxzz;
			out->mom.Oyzz += occ.mom.Oyzz;
			out->mom.Ozzz += occ.mom.Ozzz;
			out->mom.B5 += occ.mom.B5;
		default:
			out->mom.Qxx += occ.mom.Qxx;
			out->mom.Qyy += occ.mom.Qyy;
			out->mom.Qzz += occ.mom.Qzz;
			out->mom.Qxy += occ.mom.Qxy;
			out->mom.Qxz += occ.mom.Qxz;
			out->mom.Qyz += occ.mom.Qyz;
			out->mom.B2 += occ.mom.B2;
			out->mom.B3 += occ.mom.B3;
			out->mom.B4 += occ.mom.B4;
			}
		if (occ.mom.Bmax > out->mom.Bmax) out->mom.Bmax = occ.mom.Bmax;
		}
	else {
		pkdCalcCell(plcl->pkd,NULL,in->rcm,in->iOrder,&pcc);
		switch (in->iOrder) {
		case 4:
			out->mom.Hxxxx = pcc.Hxxxx;
			out->mom.Hxyyy = pcc.Hxyyy;
			out->mom.Hxxxy = pcc.Hxxxy;
			out->mom.Hyyyy = pcc.Hyyyy;
			out->mom.Hxxxz = pcc.Hxxxz;
			out->mom.Hyyyz = pcc.Hyyyz;
			out->mom.Hxxyy = pcc.Hxxyy;
			out->mom.Hxxyz = pcc.Hxxyz;
			out->mom.Hxyyz = pcc.Hxyyz;
			out->mom.Hxxzz = pcc.Hxxzz;
			out->mom.Hxyzz = pcc.Hxyzz;
			out->mom.Hxzzz = pcc.Hxzzz;
			out->mom.Hyyzz = pcc.Hyyzz;
			out->mom.Hyzzz = pcc.Hyzzz;
			out->mom.Hzzzz = pcc.Hzzzz;
			out->mom.B6 = pcc.B6;
		case 3:	
			out->mom.Oxxx = pcc.Oxxx;
			out->mom.Oxyy = pcc.Oxyy;
			out->mom.Oxxy = pcc.Oxxy;
			out->mom.Oyyy = pcc.Oyyy;
			out->mom.Oxxz = pcc.Oxxz;
			out->mom.Oyyz = pcc.Oyyz;
			out->mom.Oxyz = pcc.Oxyz;
			out->mom.Oxzz = pcc.Oxzz;
			out->mom.Oyzz = pcc.Oyzz;
			out->mom.Ozzz = pcc.Ozzz;
			out->mom.B5 = pcc.B5;
		default:
			out->mom.Qxx = pcc.Qxx;
			out->mom.Qyy = pcc.Qyy;
			out->mom.Qzz = pcc.Qzz;
			out->mom.Qxy = pcc.Qxy;
			out->mom.Qxz = pcc.Qxz;
			out->mom.Qyz = pcc.Qyz;
			out->mom.B2 = pcc.B2;
			out->mom.B3 = pcc.B3;
			out->mom.B4 = pcc.B4;
			}
		out->mom.Bmax = pcc.Bmax;
		}
	if (pnOut) *pnOut = sizeof(struct outCalcCell);
	}


void pstColCells(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	struct inColCells *in = vin;
	KDN *pkdn = vout;
	KDN *ptmp;
	int i,iCell;

	assert(nIn == sizeof(struct inColCells));
	iCell = in->iCell;
	if (pst->nLeaves > 1) {
		in->iCell = UPPER(iCell);
		mdlReqService(pst->mdl,pst->idUpper,PST_COLCELLS,in,nIn);
		in->iCell = LOWER(iCell);
		pstColCells(pst->pstLower,in,nIn,vout,pnOut);
		in->iCell = iCell;
		ptmp = malloc(in->nCell*sizeof(KDN));
		assert(ptmp != NULL);
		mdlGetReply(pst->mdl,pst->idUpper,ptmp,pnOut);
		for (i=1;i<in->nCell;++i) {
			if (ptmp[i].pUpper) {
				pkdn[i] = ptmp[i];
				}
			}
		free(ptmp);
		pst->kdn.iLower = LOWER(iCell);
		pst->kdn.iUpper = UPPER(iCell);
		pkdn[iCell] = pst->kdn;
		assert(pkdn[iCell].pUpper != 0);
		}
	else {
		for (i=1;i<in->nCell;++i) pkdn[i].pUpper = 0; /* used flag = unused */
		pst->kdn.iLower = -1;
		pst->kdn.iUpper = -1;
		pkdn[iCell] = pst->kdn;
		assert(pkdn[iCell].pUpper != 0);
		}
	if (pnOut) *pnOut = in->nCell*sizeof(KDN);
	}


void pstDistribCells(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	KDN *pkdn = vin;
	int nCell;

	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_DISTRIBCELLS,vin,nIn);
		pstDistribCells(pst->pstLower,vin,nIn,NULL,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
		}
	else {
		nCell = nIn/sizeof(KDN);
		pkdDistribCells(plcl->pkd,nCell,pkdn);
		}
	if (pnOut) *pnOut = 0;
	}


void pstCalcRoot(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct ioCalcRoot *out = vout;
	struct ioCalcRoot occ;

	assert(nIn == 0);
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_CALCROOT,vin,nIn);
		pstCalcRoot(pst->pstLower,vin,nIn,out,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,&occ,NULL);
		out->ilcn.xxxx += occ.ilcn.xxxx;
		out->ilcn.xyyy += occ.ilcn.xyyy;
		out->ilcn.xxxy += occ.ilcn.xxxy;
		out->ilcn.yyyy += occ.ilcn.yyyy;
		out->ilcn.xxxz += occ.ilcn.xxxz;
		out->ilcn.yyyz += occ.ilcn.yyyz;
		out->ilcn.xxyy += occ.ilcn.xxyy;
		out->ilcn.xxyz += occ.ilcn.xxyz;
		out->ilcn.xyyz += occ.ilcn.xyyz;
		out->ilcn.xxzz += occ.ilcn.xxzz;
		out->ilcn.xyzz += occ.ilcn.xyzz;
		out->ilcn.xzzz += occ.ilcn.xzzz;
		out->ilcn.yyzz += occ.ilcn.yyzz;
		out->ilcn.yzzz += occ.ilcn.yzzz;
		out->ilcn.zzzz += occ.ilcn.zzzz;
		out->ilcn.xxx += occ.ilcn.xxx;
		out->ilcn.xyy += occ.ilcn.xyy;
		out->ilcn.xxy += occ.ilcn.xxy;
		out->ilcn.yyy += occ.ilcn.yyy;
		out->ilcn.xxz += occ.ilcn.xxz;
		out->ilcn.yyz += occ.ilcn.yyz;
		out->ilcn.xyz += occ.ilcn.xyz;
		out->ilcn.xzz += occ.ilcn.xzz;
		out->ilcn.yzz += occ.ilcn.yzz;
		out->ilcn.zzz += occ.ilcn.zzz;
		}
	else {
		pkdCalcRoot(plcl->pkd,&occ.ilcn);
		out->ilcn.xxxx = occ.ilcn.xxxx;
		out->ilcn.xyyy = occ.ilcn.xyyy;
		out->ilcn.xxxy = occ.ilcn.xxxy;
		out->ilcn.yyyy = occ.ilcn.yyyy;
		out->ilcn.xxxz = occ.ilcn.xxxz;
		out->ilcn.yyyz = occ.ilcn.yyyz;
		out->ilcn.xxyy = occ.ilcn.xxyy;
		out->ilcn.xxyz = occ.ilcn.xxyz;
		out->ilcn.xyyz = occ.ilcn.xyyz;
		out->ilcn.xxzz = occ.ilcn.xxzz;
		out->ilcn.xyzz = occ.ilcn.xyzz;
		out->ilcn.xzzz = occ.ilcn.xzzz;
		out->ilcn.yyzz = occ.ilcn.yyzz;
		out->ilcn.yzzz = occ.ilcn.yzzz;
		out->ilcn.zzzz = occ.ilcn.zzzz;
		out->ilcn.xxx = occ.ilcn.xxx;
		out->ilcn.xyy = occ.ilcn.xyy;
		out->ilcn.xxy = occ.ilcn.xxy;
		out->ilcn.yyy = occ.ilcn.yyy;
		out->ilcn.xxz = occ.ilcn.xxz;
		out->ilcn.yyz = occ.ilcn.yyz;
		out->ilcn.xyz = occ.ilcn.xyz;
		out->ilcn.xzz = occ.ilcn.xzz;
		out->ilcn.yzz = occ.ilcn.yzz;
		out->ilcn.zzz = occ.ilcn.zzz;
		}
	if (pnOut) *pnOut = sizeof(struct ioCalcRoot);
	}


void pstDistribRoot(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct ioCalcRoot *in = vin;

	assert(nIn == sizeof(struct ioCalcRoot));
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_DISTRIBROOT,vin,nIn);
		pstDistribRoot(pst->pstLower,vin,nIn,NULL,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
		}
	else {
		pkdDistribRoot(plcl->pkd,&in->ilcn);
		}
	if (pnOut) *pnOut = 0;
	}


void pstMassCheck(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct outMassCheck *out = vout;
	double dMass;
	
	assert(nIn == 0);
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_MASSCHECK,NULL,0);
		pstMassCheck(pst->pstLower,NULL,0,vout,pnOut);
		dMass = out->dMass;
		mdlGetReply(pst->mdl,pst->idUpper,vout,pnOut);
		out->dMass += dMass;
		}
	else {
		out->dMass = pkdMassCheck(plcl->pkd);
		}
	if (pnOut) *pnOut = sizeof(struct outMassCheck);
	}

void
pstSetRung(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inSetRung *in = vin;
	
	assert(nIn == sizeof(*in));
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_SETRUNG,vin,nIn);
		pstSetRung(pst->pstLower,vin,nIn,NULL,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
		}
	else {
		pkdSetRung(plcl->pkd, in->iRung);
		}
	if (pnOut) *pnOut = 0;
	}

void
pstActiveRung(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inActiveRung *in = vin;
	
	assert(nIn == sizeof(*in));
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_ACTIVERUNG,vin,nIn);
		pstActiveRung(pst->pstLower,vin,nIn,NULL,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
		}
	else {
		pkdActiveRung(plcl->pkd, in->iRung, in->bGreater);
		}
	if (pnOut) *pnOut = 0;
	}

void
pstCurrRung(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inCurrRung *in = vin;
	struct outCurrRung *out = vout;
	int iCurrent;
	
	assert(nIn == sizeof(*in));
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_CURRRUNG,vin,nIn);
		pstCurrRung(pst->pstLower,vin,nIn,vout,pnOut);
		iCurrent = out->iCurrent;
		mdlGetReply(pst->mdl,pst->idUpper,vout,pnOut);
		if(iCurrent)
		    out->iCurrent = iCurrent;
		}
	else {
		out->iCurrent = pkdCurrRung(plcl->pkd, in->iRung);
		}
	if (pnOut) *pnOut = sizeof(*out);
	}

void
pstDensityRung(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inDensityRung *in = vin;
	struct outDensityRung *out = vout;
	int iMaxRung;
	
	assert(nIn == sizeof(*in));
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_DENSITYRUNG,vin,nIn);
		pstDensityRung(pst->pstLower,vin,nIn,vout,pnOut);
		iMaxRung = out->iMaxRung;
		mdlGetReply(pst->mdl,pst->idUpper,vout,pnOut);
		if(iMaxRung > out->iMaxRung)
		    out->iMaxRung = iMaxRung;
		}
	else {
		out->iMaxRung = pkdDensityRung(plcl->pkd, in->iRung,
					       in->dDelta, in->dEta,
					       in->dRhoFac, in->bAll);
		}
	if (pnOut) *pnOut = sizeof(*out);
	}

void
pstVelocityRung(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
      LCL *plcl = pst->plcl;
      struct inVelocityRung *in = vin;
      struct outVelocityRung *out = vout;
      int iMaxRung;

      assert(nIn == sizeof(*in));
      if (pst->nLeaves > 1) {
              mdlReqService(pst->mdl,pst->idUpper,PST_VELOCITYRUNG,vin,nIn);
              pstVelocityRung(pst->pstLower,vin,nIn,vout,pnOut);
              iMaxRung = out->iMaxRung;
              mdlGetReply(pst->mdl,pst->idUpper,vout,pnOut);
              if(iMaxRung > out->iMaxRung)
                  out->iMaxRung = iMaxRung;
              }
      else {
              out->iMaxRung = pkdVelocityRung(plcl->pkd, in->iRung,
                                             in->dDelta, in->dEta,
                                             in->iMaxRung, in->dVelFac,
					      in->dAccFac, in->bAll);
              }
      if (pnOut) *pnOut = sizeof(*out);
      }

void pstRungStats(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inRungStats *in = vin;
	struct outRungStats *out = vout;
	int nParticles;

	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_RUNGSTATS,vin,nIn);
		pstRungStats(pst->pstLower,vin,nIn,out,NULL);
		nParticles = out->nParticles;
		mdlGetReply(pst->mdl,pst->idUpper,out,NULL);
		out->nParticles += nParticles;
		}
	else {
		out->nParticles = pkdRungParticles(plcl->pkd,in->iRung);
		}
	if (pnOut) *pnOut = sizeof(struct outRungStats);
	}


