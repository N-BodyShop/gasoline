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
	mdlAddService(mdl,PST_SETADD,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstSetAdd,
				  sizeof(struct inSetAdd),0);
	mdlAddService(mdl,PST_LEVELIZE,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstLevelize,
				  sizeof(struct inLevelize),0);
	mdlAddService(mdl,PST_READTIPSY,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstReadTipsy,
				  sizeof(struct inReadTipsy),0);
	mdlAddService(mdl,PST_DOMAINDECOMP,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstDomainDecomp,
				  sizeof(struct inDomainDecomp),0);
	mdlAddService(mdl,PST_CALCBOUND,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstCalcBound,
				  0,sizeof(struct outCalcBound));
	mdlAddService(mdl,PST_GASWEIGHT,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstGasWeight,
				  0,0);
	mdlAddService(mdl,PST_WEIGHT,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstWeight,
				  sizeof(struct inWeight),sizeof(struct outWeight));
	mdlAddService(mdl,PST_ORDWEIGHT,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstOrdWeight,
				  sizeof(struct inOrdWeight),sizeof(struct outOrdWeight));
	mdlAddService(mdl,PST_FREESTORE,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstFreeStore,
				  0,sizeof(struct outFreeStore));
	mdlAddService(mdl,PST_COLREJECTS,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstColRejects,
				  sizeof(struct inColRejects),nThreads*sizeof(OREJ));
	mdlAddService(mdl,PST_SWAPREJECTS,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstSwapRejects,
				  nThreads*sizeof(int),nThreads*sizeof(OREJ));
	mdlAddService(mdl,PST_DOMAINCOLOR,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstDomainColor,
				  0,0);
	mdlAddService(mdl,PST_COLORDREJECTS,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstColOrdRejects,
				  sizeof(struct inColOrdRejects),nThreads*sizeof(OREJ));
	mdlAddService(mdl,PST_DOMAINORDER,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstDomainOrder,
				  sizeof(struct inDomainOrder),0);
	mdlAddService(mdl,PST_LOCALORDER,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstLocalOrder,
				  0,0);
	mdlAddService(mdl,PST_OUTARRAY,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstOutArray,
				  sizeof(struct inOutArray),0);
	mdlAddService(mdl,PST_OUTVECTOR,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstOutVector,
				  sizeof(struct inOutVector),0);
	mdlAddService(mdl,PST_WRITETIPSY,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstWriteTipsy,
				  sizeof(struct inWriteTipsy),0);
	mdlAddService(mdl,PST_BUILDTREE,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstBuildTree,
				  sizeof(struct inBuildTree),sizeof(struct outBuildTree));
	mdlAddService(mdl,PST_SMOOTH,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstSmooth,
				  sizeof(struct inSmooth),0);
	mdlAddService(mdl,PST_GRAVITY,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstGravity,
				  sizeof(struct inGravity),sizeof(struct outGravity));
	mdlAddService(mdl,PST_CALCE,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstCalcE,
				  0,sizeof(struct outCalcE));
	mdlAddService(mdl,PST_DRIFT,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstDrift,
				  sizeof(struct inDrift),0);
	mdlAddService(mdl,PST_KICK,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstKick,
				  sizeof(struct inKick),sizeof(struct outKick));
	mdlAddService(mdl,PST_READCHECK,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstReadCheck,
				  sizeof(struct inReadCheck), 0);
	mdlAddService(mdl,PST_WRITECHECK,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstWriteCheck,
				  sizeof(struct inWriteCheck),0);
	mdlAddService(mdl,PST_SETSOFT,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstSetSoft,
				  sizeof(struct inSetSoft),0);
	mdlAddService(mdl,PST_SETTOTAL,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstSetTotal,
				  0,sizeof(struct outSetTotal));
	mdlAddService(mdl,PST_SETWRITESTART,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstSetWriteStart,
				  sizeof(struct inSetWriteStart),0);
	mdlAddService(mdl,PST_CALCCELL,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstCalcCell,
				  sizeof(struct inCalcCell),sizeof(struct outCalcCell));
	/*
	 ** Calculate the number of levels in the top tree and use it to 
	 ** define the size of the messages.
	 */
	nCell = 1<<(1+(int)ceil(log((double)nThreads)/log(2.0)));
	mdlAddService(mdl,PST_COLCELLS,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstColCells,
				  sizeof(struct inColCells),nCell*sizeof(KDN));
	mdlAddService(mdl,PST_DISTRIBCELLS,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstDistribCells,
				  nCell*sizeof(KDN),0);
	mdlAddService(mdl,PST_CALCROOT,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstCalcRoot,
				  0,sizeof(struct ioCalcRoot));
	mdlAddService(mdl,PST_DISTRIBROOT,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstDistribRoot,
				  sizeof(struct ioCalcRoot),0);
	mdlAddService(mdl,PST_ONENODEREADINIT,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstOneNodeReadInit,
				  sizeof(struct inReadTipsy), nThreads*sizeof(int));
	mdlAddService(mdl,PST_SWAPALL,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstSwapAll,
				  sizeof(int),0);
	mdlAddService(mdl,PST_MASSCHECK,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstMassCheck,
				  0,sizeof(struct outMassCheck));
	mdlAddService(mdl,PST_ACTIVETYPEORDER,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstActiveTypeOrder,
				  sizeof(struct inActiveTypeOrder),sizeof(int));
	mdlAddService(mdl,PST_ACTIVEORDER,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstActiveOrder,
				  0,sizeof(int));
	mdlAddService(mdl,PST_SETRUNG,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstSetRung,
				  sizeof(struct inSetRung),0);
	mdlAddService(mdl,PST_BALLMAX,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstBallMax,
				  sizeof(struct inBallMax),0);
	mdlAddService(mdl,PST_ACTIVERUNG,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstActiveRung,
				  sizeof(struct inActiveRung),sizeof(int));
	mdlAddService(mdl,PST_CURRRUNG,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstCurrRung,
				  sizeof(struct inCurrRung),sizeof(struct outCurrRung));
	mdlAddService(mdl,PST_DENSITYSTEP,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstDensityStep,
				  sizeof(struct inDensityStep),0);
	mdlAddService(mdl,PST_RUNGSTATS,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstRungStats,
				  sizeof(struct inRungStats),sizeof(struct outRungStats));
	mdlAddService(mdl,PST_GETMAP,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstGetMap,
				  sizeof(struct inGetMap),nThreads*sizeof(int));
	mdlAddService(mdl,PST_ACCELSTEP,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstAccelStep,
				  sizeof(struct inAccelStep), 0);
	mdlAddService(mdl,PST_COOLVELOCITY,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstCoolVelocity,
				  sizeof(struct inCoolVelocity),0);
	mdlAddService(mdl,PST_RESETTOUCHRUNG,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstResetTouchRung,
				  sizeof(struct inResetTouchRung),sizeof(int));
	mdlAddService(mdl,PST_ACTIVEEXACTTYPE,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstActiveExactType,
				  sizeof(struct inActiveType),sizeof(int));
	mdlAddService(mdl,PST_ACTIVETYPE,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstActiveType,
				  sizeof(struct inActiveType),sizeof(int));
	mdlAddService(mdl,PST_SETTYPE,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstSetType,
				  sizeof(struct inActiveType),sizeof(int));
	mdlAddService(mdl,PST_RESETTYPE,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstResetType,
				  sizeof(struct inActiveType),sizeof(int));
	mdlAddService(mdl,PST_COUNTTYPE,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstCountType,
				  sizeof(struct inActiveType),sizeof(int));
	mdlAddService(mdl,PST_ACTIVEMASKRUNG,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstActiveMaskRung,
				  sizeof(struct inActiveType),sizeof(int));
	mdlAddService(mdl,PST_ACTIVETYPERUNG,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstActiveTypeRung,
				  sizeof(struct inActiveType),sizeof(int));
	mdlAddService(mdl,PST_SETPARTICLETYPES,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstSetParticleTypes,
				  sizeof(struct inSetParticleTypes),0);
	mdlAddService(mdl,PST_GROWMASS,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstGrowMass,
				  sizeof(struct inGrowMass),0);
	mdlAddService(mdl,PST_MARKSMOOTH,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstMarkSmooth,
				  sizeof(struct inMarkSmooth),0);
	mdlAddService(mdl,PST_RESMOOTH,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstReSmooth,
				  sizeof(struct inReSmooth),0);
	mdlAddService(mdl,PST_INITACCEL,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstInitAccel,0,0);
	mdlAddService(mdl,PST_DTTORUNG,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstDtToRung,
				  sizeof(struct inDtToRung),sizeof(struct outDtToRung));
	mdlAddService(mdl,PST_INITDT,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstInitDt,
				  sizeof(struct inInitDt),0);
#ifdef GASOLINE
	mdlAddService(mdl,PST_UPDATEUDOT,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstUpdateuDot,
				  sizeof(struct inUpdateuDot),sizeof(struct outUpdateuDot));
	mdlAddService(mdl,PST_SPHCURRRUNG,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstSphCurrRung,
				  sizeof(struct inSphCurrRung),sizeof(struct outSphCurrRung));
	mdlAddService(mdl,PST_GETGASPRESSURE,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstGetGasPressure,
				  sizeof(struct inGetGasPressure),0);
	mdlAddService(mdl,PST_INITENERGY,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstInitEnergy,
				  sizeof(struct inInitEnergy),0);
	mdlAddService(mdl,PST_KICKVPRED,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstKickVpred, 
				  sizeof(struct inKickVpred),sizeof(struct outKick));
	mdlAddService(mdl,PST_KICKRHOPRED,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstKickRhopred, 
				  sizeof(struct inKickRhopred),0);
	mdlAddService(mdl,PST_SPHSTEP,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstSphStep,
				  sizeof(struct inSphStep),0);
	mdlAddService(mdl,PST_SPHVISCOSITYLIMITER,pst,
				  (void (*)(void *,void *,int,void *,int *)) 
				  pstSphViscosityLimiter, 
				  sizeof(struct inSphViscosityLimiter),0);
	mdlAddService(mdl,PST_INITCOOLING,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstInitCooling,
				  sizeof(struct inInitCooling),0);
	mdlAddService(mdl,PST_INITUV,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstInitUV,
				  CL_NMAXBYTETABLE,0);
	mdlAddService(mdl,PST_DENSCHECK,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstDensCheck,
				  sizeof(struct inDensCheck),sizeof(struct outDensCheck));
#endif /* GASOLINE */
#ifdef GLASS
	mdlAddService(mdl,PST_RANDOMVELOCITIES,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstRandomVelocities, 
				  sizeof(struct inRandomVelocities),0);
#endif
	mdlAddService(mdl,PST_COLNPARTS,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstColNParts,
				  0,nThreads*sizeof(struct outColNParts));
	mdlAddService(mdl,PST_NEWORDER,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstNewOrder,
				  nThreads*sizeof(int),0);
	mdlAddService(mdl,PST_SETNPARTS,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstSetNParts,
				  sizeof(struct inSetNParts),0);
	mdlAddService(mdl,PST_GRAVEXTERNAL,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstGravExternal,
				  sizeof(struct inGravExternal),0);
#ifdef COLLISIONS
	mdlAddService(mdl,PST_NUMREJECTS,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstNumRejects,
				  0,sizeof(struct outNumRejects));
	mdlAddService(mdl,PST_READSS,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstReadSS,
				  sizeof(struct inReadSS),0);
	mdlAddService(mdl,PST_WRITESS,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstWriteSS,
				  sizeof(struct inWriteSS),0);
	mdlAddService(mdl,PST_CALCHILL,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstCalcHill,
				  sizeof(struct inCalcHill),0);
	mdlAddService(mdl,PST_HELIOSTEP,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstHelioStep,
				  sizeof(struct inHelioStep),0);
	mdlAddService(mdl,PST_KICKUNIFGRAV,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstKickUnifGrav,
				  sizeof(struct inKickUnifGrav),0);
	mdlAddService(mdl,PST_NEXTENCOUNTER,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstNextEncounter,
				  0,sizeof(struct outNextEncounter));
	mdlAddService(mdl,PST_MARKENCOUNTERS,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstMarkEncounters,
				  sizeof(struct inMarkEncounters),0);
	mdlAddService(mdl,PST_NEXTCOLLISION,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstNextCollision,
				  0,sizeof(struct outNextCollision));
	mdlAddService(mdl,PST_GETCOLLIDERINFO,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstGetColliderInfo,
				  sizeof(struct inGetColliderInfo),
				  sizeof(struct outGetColliderInfo));
	mdlAddService(mdl,PST_DOCOLLISION,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstDoCollision,
				  sizeof(struct inDoCollision),sizeof(struct outDoCollision));
	mdlAddService(mdl,PST_RESETCOLLIDERS,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstResetColliders,
				  sizeof(struct inResetColliders),0);
	mdlAddService(mdl,PST_QQCALCBOUND,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstQQCalcBound,
				  0,sizeof(struct outCalcBound));
	mdlAddService(mdl,PST_QQDOMAINDECOMP,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstQQDomainDecomp,
				  0,0);
	mdlAddService(mdl,PST_QQBUILDTREE,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstQQBuildTree,
				  sizeof(struct inBuildTree),sizeof(struct outBuildTree));
	mdlAddService(mdl,PST_QQSMOOTH,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstQQSmooth,
				  sizeof(struct inSmooth),0);
#endif /* COLLISIONS */
#ifdef SLIDING_PATCH
	mdlAddService(mdl,PST_KICKVPRED,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstKickVpred, 
				  sizeof(struct inKickVpred),0);
#endif
	mdlAddService(mdl,PST_CLEARTIMER,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstClearTimer,
				  sizeof(struct inClearTimer),0);
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
		if(pst->nLeaves == 1 && pst->plcl->pkd)
			pkdFinish(pst->plcl->pkd);
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
	int nFileStart,nFileEnd,nFileTotal,nFileSplit,nStore;
	int *ptmp;
	int i;

	assert(nIn == sizeof(struct inReadTipsy));
	nFileStart = in->nFileStart;
	nFileEnd = in->nFileEnd;
	nFileTotal = nFileEnd - nFileStart + 1;
	if (pst->nLeaves > 1) {
		nFileSplit = nFileStart + pst->nLower*(nFileTotal/pst->nLeaves);
		in->nFileStart = nFileSplit;
		mdlReqService(pst->mdl,pst->idUpper,PST_ONENODEREADINIT,in,nIn);
		in->nFileStart = nFileStart;
		in->nFileEnd = nFileSplit - 1;
		pstOneNodeReadInit(pst->pstLower,in,nIn,vout,pnOut);
		in->nFileEnd = nFileEnd;
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
		nStore = nFileTotal + (int)ceil(nFileTotal*in->fExtraStore);
		
		pkdInitialize(&plcl->pkd,pst->mdl,in->iOrder,nStore,
					  plcl->nPstLvl,in->fPeriod,in->nDark,in->nGas,in->nStar);
		pout[pst->idSelf] = nFileTotal;
		}
	if (pnOut) *pnOut = pst->mdl->nThreads*sizeof(*pout);
	}

void pstReadTipsy(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inReadTipsy *in = vin;
	int nFileStart,nFileEnd,nFileTotal,nFileSplit,nStore;
	char achInFile[PST_FILENAME_SIZE];

	assert(nIn == sizeof(struct inReadTipsy));
	nFileStart = in->nFileStart;
	nFileEnd = in->nFileEnd;
	nFileTotal = nFileEnd - nFileStart + 1;
	if (pst->nLeaves > 1) {
		nFileSplit = nFileStart + pst->nLower*(nFileTotal/pst->nLeaves);
		in->nFileStart = nFileSplit;
		mdlReqService(pst->mdl,pst->idUpper,PST_READTIPSY,in,nIn);
		in->nFileStart = nFileStart;
		in->nFileEnd = nFileSplit - 1;
		pstReadTipsy(pst->pstLower,in,nIn,NULL,NULL);
		in->nFileEnd = nFileEnd;
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
		nStore = nFileTotal + (int)ceil(nFileTotal*in->fExtraStore);
		pkdInitialize(&plcl->pkd,pst->mdl,in->iOrder,nStore,plcl->nPstLvl,
					  in->fPeriod,in->nDark,in->nGas,in->nStar);
		pkdReadTipsy(plcl->pkd,achInFile,nFileStart,nFileTotal,in->bStandard,
					 in->dvFac,in->dTuFac);
		}
	if (pnOut) *pnOut = 0;
	}


int _pstRejMatch(PST pst,int n1,OREJ *p1,int n2,OREJ *p2,int *pidSwap)
{
	int id,i,i1=-1,i2=-1,nLarge,id1,id2;
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


#define MAX_ITTR	64
#define EPS_BOUND	0.01
#define MASS_EPS	1e-11

void _pstRootSplit(PST pst,int iSplitDim,double dMass, int bDoRootFind)
{
	int NUM_SAFETY = 4; 	/* slop space when filling up memory */
	int nSafeTot;		/* total slop space we have to play with */
	int d,ittr,nOut;
	int nLow=-1,nHigh=-1,nLowerStore,nUpperStore;
	FLOAT fLow,fHigh;
	FLOAT fl,fu,fm=-1,fmm;
	struct outFreeStore outFree;
	struct inWeight inWt;
	struct outWeight outWtLow;
	struct outWeight outWtHigh;
	struct inColRejects inCol;
	OREJ *pLowerRej,*pUpperRej;
	int *pidSwap,iRet;
	int nLowTot,nHighTot;
	int nActive;
	int nLast;		/* number of particles at the last split
					   iteration */
	int nDiff=0;	/* Difference between one iteration and the next */
	char ach[256];	/* Debug */

	struct outMassCheck outMass;
#ifdef PARANOID_CHECK
	int i,iLowSum,iHighSum;
#endif
	int pFlag,	/* 0 => we are splitting all particles at once. 
				   1 => we first split active, and then inactive. */
	    nTotalActive;
	int dBnd;

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
	pstActiveOrder(pst->pstLower,NULL,0,&nActive,NULL);
	mdlGetReply(pst->mdl,pst->idUpper,&nActive,NULL);

	/* Debug */
	/*
	sprintf(ach,"id: %d _pstRootSplit\n", pst->mdl->idSelf );
	mdlDiag(pst->mdl,ach);
	*/

	if(!bDoRootFind) {
	    d = dBnd = pst->iSplitDim;
#ifdef COLLISIONS
	    if(d > 2)
		dBnd -= 3;
#else
	    assert(d < 3);
#endif
	    fm = pst->fSplit;
	    pFlag = 0;
	    nLow = 0;
	    nHigh = 0;
	    goto DoneRootFind;
	    }
	    
	d = dBnd = iSplitDim;
#ifdef COLLISIONS
	if(d > 2)
	    dBnd -= 3;
#else
	assert(d < 3);
#endif
	fl = pst->bnd.fMin[dBnd];
	fu = pst->bnd.fMax[dBnd];
	fmm = (fl + fu)/2;
	/*
	 * First find total number of active particles.
	 */
	inWt.iSplitDim = d;
	inWt.fSplit = fmm;
	inWt.ittr = 0;
	inWt.iSplitSide = 1;
	inWt.pFlag = 1;
	mdlReqService(pst->mdl,pst->idUpper,PST_WEIGHT,&inWt,sizeof(inWt));
	inWt.iSplitSide = 0;
	pstWeight(pst->pstLower,&inWt,sizeof(inWt),&outWtLow,NULL);
	mdlGetReply(pst->mdl,pst->idUpper,&outWtHigh,NULL);
	nTotalActive = outWtLow.nLow + outWtHigh.nLow
	    + outWtLow.nHigh + outWtHigh.nHigh;
	/*
	 * If we only have a few active particles, just split
	 * all the particles.
	 */
	if(nTotalActive < 16) pFlag = 0;
	else pFlag = 1;
	/*
	 ** Now start the ROOT finder based on balancing active weight ALONE!
	 ** (unless pFlag == 0)
	 */
	ittr = 0;
	while (fl < fmm && fmm < fu && ittr < MAX_ITTR) {
		fm = fmm;
		inWt.iSplitDim = d;
		inWt.fSplit = fm;
		inWt.ittr = ittr;
		inWt.iSplitSide = 1;
		inWt.pFlag = pFlag;
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
		if(nLow == 1 && nHigh == 1) /* break on trivial case */
		    break;
		if(pFlag) {	/* split on work */
		    if (fLow/pst->nLower > fHigh/pst->nUpper) fu = fm;
		    else if (fLow/pst->nLower < fHigh/pst->nUpper) fl = fm;
		    else break;
		    }
		else {		/* split on number */
		    if (nLow/(double)pst->nLower >
			nHigh/(double)pst->nUpper) fu = fm;
		    else if (nLow/(double)pst->nLower <
			     nHigh/(double)pst->nUpper) fl = fm;
		    else break;
		    }
		fmm = (fl + fu)/2;
		++ittr;
		}
	/*
	 ** If we exceed the local pStore in the lower or upper subsets then
	 ** we have to relax the condition that work be balanced and try to
	 ** max out the number of particles in the subset which had too many.
	 */
	nSafeTot = nLowerStore + nUpperStore - nTotalActive;
	if(nSafeTot/pst->nLeaves < NUM_SAFETY) {
		NUM_SAFETY = nSafeTot/pst->nLeaves;
		fprintf(stderr, "id: %d tripped active NUM_SAFETY %d\n",
		pst->mdl->idSelf, NUM_SAFETY);
		}
	if (nLow > nLowerStore-NUM_SAFETY*pst->nLower) {
		fl = pst->bnd.fMin[dBnd];
		fu = fm;
		fmm = (fl + fu)/2;
		ittr = 0;
	    while (fl < fmm && fmm < fu && ittr < MAX_ITTR) {
			fm = fmm;
			inWt.iSplitDim = d;
			inWt.fSplit = fm;
			inWt.ittr = ittr;
			inWt.iSplitSide = 1;
			inWt.pFlag = pFlag;
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
			if (nLow > nLowerStore-NUM_SAFETY*pst->nLower) fu = fm;
			else if (nLow < nLowerStore-NUM_SAFETY*pst->nLower) fl = fm;
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
	else if (nHigh > nUpperStore-NUM_SAFETY*pst->nUpper) {
		fl = fm;
		fu = pst->bnd.fMax[dBnd];
		fmm = (fl + fu)/2;
		ittr = 0;
	    while (fl < fmm && fmm < fu && ittr < MAX_ITTR) {
			fm = fmm;
			inWt.iSplitDim = d;
			inWt.fSplit = fm;
			inWt.ittr = ittr;
			inWt.iSplitSide = 1;
			inWt.pFlag = pFlag;
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
			if (nHigh > nUpperStore-NUM_SAFETY*pst->nUpper) fl = fm;
			else if (nHigh < nUpperStore-NUM_SAFETY*pst->nUpper) fu = fm;
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

 DoneRootFind:
	if(pFlag || !bDoRootFind) {
	    /*
	     ** Now we see if the TOTAL number of particles in the lower and upper
	     ** subsets exceeds the local pStores. If so then we need to find a new
	     ** boundary to distribute the INACTIVE particles so that everything 
	     ** fits.  However, if we did not do a rootfind, we need
	     ** to split ALL particles.
	     */
	    inWt.iSplitDim = d;
	    inWt.fSplit = fm;
	    inWt.ittr = 0;
	    inWt.iSplitSide = 1;
	    inWt.pFlag = (!bDoRootFind ? 0 : -1);
	    mdlReqService(pst->mdl,pst->idUpper,PST_WEIGHT,&inWt,sizeof(inWt));
	    inWt.iSplitSide = 0;
	    pstWeight(pst->pstLower,&inWt,sizeof(inWt),&outWtLow,NULL);
	    mdlGetReply(pst->mdl,pst->idUpper,&outWtHigh,NULL);
	    /*
	     ** Add lower and Upper subsets numbers of particles
	     */ 
	    nLowTot = nLow + outWtLow.nLow + outWtHigh.nLow;
	    nHighTot = nHigh + outWtLow.nHigh + outWtHigh.nHigh;

	    nSafeTot = nLowerStore + nUpperStore - (nLowTot + nHighTot);
	    if(nSafeTot/pst->nLeaves < NUM_SAFETY) {
		NUM_SAFETY = nSafeTot/pst->nLeaves;
		sprintf(ach,"id: %d tripped inactive NUM_SAFETY %d\n",
				pst->mdl->idSelf, NUM_SAFETY);
		mdlDiag(pst->mdl,ach);
		fprintf(stderr, "id: %d tripped inactive NUM_SAFETY %d\n",
		pst->mdl->idSelf, NUM_SAFETY);
		}
	    if (nLowTot > nLowerStore-NUM_SAFETY*pst->nLower) {
			sprintf(ach,"id: %d: nLowTot > nLowerStore-NUM_SAFETY*pst->nLower %d %d %d %d\n",
					pst->mdl->idSelf, nLowTot, nLowerStore, NUM_SAFETY, pst->nLower);
		    mdlDiag(pst->mdl,ach);
		    fl = pst->bnd.fMin[dBnd];
		    fu = fm;
		    fmm = (fl + fu)/2;
		    ittr = 1;
		    nLast = nLowTot;
			while (fl < fmm && fmm < fu && ittr < MAX_ITTR) {
				fm = fmm;
				inWt.iSplitDim = d;
				inWt.fSplit = fm;
				inWt.ittr = ittr;
				inWt.iSplitSide = 1;
				inWt.pFlag = (!bDoRootFind ? 0 : -1);
				mdlReqService(pst->mdl,pst->idUpper,PST_WEIGHT,&inWt,sizeof(inWt));
				inWt.iSplitSide = 0;
				pstWeight(pst->pstLower,&inWt,sizeof(inWt),&outWtLow,NULL);
				mdlGetReply(pst->mdl,pst->idUpper,&outWtHigh,NULL);
				/*
				 ** Add lower and Upper subsets numbers of particles
				 */
				nLowTot = nLow + outWtLow.nLow + outWtHigh.nLow;
				nHighTot = nHigh + outWtLow.nHigh + outWtHigh.nHigh;
				if(nLowTot != nLast)
					nDiff = nLowTot - nLast;
				nLast = nLowTot;
				if (nLowTot > nLowerStore-NUM_SAFETY*pst->nLower) fu = fm;
				else if (nLowTot < nLowerStore-NUM_SAFETY*pst->nLower) fl = fm;
				else {
					fl = fm;
					break;
					}
				fmm = (fl + fu)/2;
				++ittr;
				}
		    if(nLowTot != nLowerStore-NUM_SAFETY*pst->nLower) {
				if(abs(nDiff) > 1)
					fprintf(stderr, "id: %d delta of %d, check NUM_SAFTEY\n",
							pst->mdl->idSelf, nDiff);
				}
		    assert(nLowTot <= nLowerStore);
		    }
	    else if (nHighTot > nUpperStore-NUM_SAFETY*pst->nUpper) {
			sprintf(ach,"id: %d: nHighTot > nUpperStore-NUM_SAFETY*pst->nUpper %d %d %d %d\n",
					pst->mdl->idSelf, nHighTot, nUpperStore, NUM_SAFETY, pst->nUpper);
		    mdlDiag(pst->mdl,ach);
		    fl = fm;
		    fu = pst->bnd.fMax[dBnd];
		    fmm = (fl + fu)/2;
		    ittr = 1;
		    nLast = nLowTot;
			while (fl < fmm && fmm < fu && ittr < MAX_ITTR) {
				fm = fmm;
				inWt.iSplitDim = d;
				inWt.fSplit = fm;
				inWt.ittr = ittr;
				inWt.iSplitSide = 1;
				inWt.pFlag = (!bDoRootFind ? 0 : -1);
				mdlReqService(pst->mdl,pst->idUpper,PST_WEIGHT,&inWt,sizeof(inWt));
				inWt.iSplitSide = 0;
				pstWeight(pst->pstLower,&inWt,sizeof(inWt),&outWtLow,NULL);
				mdlGetReply(pst->mdl,pst->idUpper,&outWtHigh,NULL);
				/*
				 ** Add lower and Upper subsets numbers of particles
				 */
				nLowTot = nLow + outWtLow.nLow + outWtHigh.nLow;
				nHighTot = nHigh + outWtLow.nHigh + outWtHigh.nHigh;
				if(nLowTot != nLast)
					nDiff = nLowTot - nLast;
				nLast = nLowTot;
				if (nHighTot > nUpperStore-NUM_SAFETY*pst->nUpper) fl = fm;
				else if (nHighTot < nUpperStore-NUM_SAFETY*pst->nUpper) fu = fm;
				else {
					fu = fm;
					break;
					}
				fmm = (fl + fu)/2;
				++ittr;
				}
		    if(nHighTot != nUpperStore-NUM_SAFETY*pst->nUpper) {
				if(abs(nDiff) > 1)
					fprintf(stderr, "id: %d delta of %d, check NUM_SAFTEY\n",
							pst->mdl->idSelf, nDiff);
				}
		    assert(nHighTot <= nUpperStore);
		    }
	    /*
	     ** Now we make sure that there is at least one particle per
	     ** processor.  If not, redistribute the INACTIVE particles so
	     ** every processor has at least one.
	     */
 	    if (nLowTot < NUM_SAFETY*pst->nLower) {
			sprintf(ach,"id: %d: nLowTot < NUM_SAFETY*pst->nLower %d %d %d\n",
					pst->mdl->idSelf, nLowTot, NUM_SAFETY, pst->nLower);
		    mdlDiag(pst->mdl,ach);
		    fl = fm;
			/* Try to catch highest particle if needed. */
		    fu = pst->bnd.fMax[dBnd]
				+ EPS_BOUND*(pst->bnd.fMax[dBnd] - pst->bnd.fMin[dBnd]);
		    fmm = (fl + fu)/2;
		    ittr = 1;
			while (fl < fmm && fmm < fu && ittr < MAX_ITTR) {
				fm = fmm;
				inWt.iSplitDim = d;
				inWt.fSplit = fm;
				inWt.ittr = ittr;
				inWt.iSplitSide = 1;
				inWt.pFlag = (!bDoRootFind ? 0 : -1);
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
				if (nLowTot > NUM_SAFETY*pst->nLower) fu = fm;
				else if (nLowTot < NUM_SAFETY*pst->nLower) fl = fm;
				else {
					fu = fm;
					break;
					}
				fmm = (fl + fu)/2;
				++ittr;
				}
		    }
	    else if (nHighTot < NUM_SAFETY*pst->nUpper) {
			sprintf(ach,"id: %d: nHighTot < NUM_SAFETY*pst->nUpper %d %d %d\n",
					pst->mdl->idSelf, nHighTot, NUM_SAFETY, pst->nUpper);
		    mdlDiag(pst->mdl,ach);
			/* Try to catch lowest particle if needed. */
		    fl = pst->bnd.fMin[dBnd]
				- EPS_BOUND*(pst->bnd.fMax[dBnd] - pst->bnd.fMin[dBnd]);
		    fu = fm;
		    fmm = (fl + fu)/2;
		    ittr = 1;
			while (fl < fmm && fmm < fu && ittr < MAX_ITTR) {
				fm = fmm;
				inWt.iSplitDim = d;
			    inWt.fSplit = fm;
			    inWt.ittr = ittr;
			    inWt.iSplitSide = 1;
			    inWt.pFlag = (!bDoRootFind ? 0 : -1);
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
			    if (nHighTot > NUM_SAFETY*pst->nUpper) fl = fm;
			    else if (nHighTot < NUM_SAFETY*pst->nUpper) fu = fm;
			    else {
				    fl = fm;
				    break;
				    }
			    fmm = (fl + fu)/2;
			    ++ittr;
			    }
		    }
		/*
		   fprintf(stderr, "id: %dFit stats: Active: %d %d Inactive: %d %d Sum: %d %d Space: %d %d\n",
		   pst->mdl->idSelf, nLow, nHigh,
		   outWtLow.nLow + outWtHigh.nLow, outWtLow.nHigh
		   + outWtHigh.nHigh,
		   nLowTot, nHighTot, nLowerStore, nUpperStore);
		   */
	    /*
	     ** Make sure everything is OK.
	     */
	    assert(nLowTot >= pst->nLower);
	    assert(nHighTot >= pst->nUpper);
	    assert(nLowTot <= nLowerStore);
	    assert(nHighTot <= nUpperStore);
	    }
	else {
	    /*
	     ** Make sure that the particles are back in
	     ** active-inactive order after we've split all the particles.
	     */
	    mdlReqService(pst->mdl,pst->idUpper,PST_ACTIVEORDER,NULL,0);
	    pstActiveOrder(pst->pstLower,NULL,0,&nActive,NULL);
	    mdlGetReply(pst->mdl,pst->idUpper,&nActive,NULL);
	    /*
		   printf("Fit stats: Sum: %d %d\n", nLow, nHigh);
		   */
	    /*
	     ** Make sure everything is OK.
	     */
	    assert(nLow >= pst->nLower);
	    assert(nHigh >= pst->nUpper);
	    assert(nLow <= nLowerStore);
	    assert(nHigh <= nUpperStore);
	    }
    
	pst->fSplitInactive = fm;
	if(!bDoRootFind)
	    pst->fSplit = fm;
	
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
#if 0
	if (dMass != outMass.dMass) {
#else
	if (fabs(dMass - outMass.dMass) > MASS_EPS*dMass) {
#endif
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
#if 0
		if (dMass != outMass.dMass) {
#else
		if (fabs(dMass - outMass.dMass) > MASS_EPS*dMass) {
#endif
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
	struct inDomainDecomp *in = vin;
	
	char ach[256]; /* debug */

	pkdStartTimer(pst->plcl->pkd,6);
	assert(nIn == sizeof(struct inDomainDecomp));

	if (pst->nLeaves > 1) {
		/*
		 ** First calculate the Bounds for the set.
		 */
		pstCalcBound(pst,NULL,0,&outBnd,&nOut);
		pst->bnd = outBnd.bnd;
		pst->bndActive = outBnd.bndActive;
		pst->bndTreeActive = outBnd.bndTreeActive;
		/*
		 ** Next determine the longest axis based on the active bounds.
		 */
		d = 0;
		if(pst->bndActive.fMax[0] > pst->bndActive.fMin[0]) {
		    for (j=1;j<3;++j) {
				if (pst->bndActive.fMax[j]-pst->bndActive.fMin[j] > 
					pst->bndActive.fMax[d]-pst->bndActive.fMin[d]) d = j;
			    }
			}
		else {		/* no active particles */
		    for (j=1;j<3;++j) {
				if (pst->bnd.fMax[j]-pst->bnd.fMin[j] > 
					pst->bnd.fMax[d]-pst->bnd.fMin[d]) d = j;
			    }
			}

		/* pst->iSplitDim = d; */
		
		pstMassCheck(pst,NULL,0,&outMass,NULL);
		dMass = outMass.dMass;
		_pstRootSplit(pst,d,dMass, in->bDoRootFind);
		pstMassCheck(pst,NULL,0,&outMass,NULL);
#if 0
		if (dMass != outMass.dMass) {
#else
  if (fabs(dMass - outMass.dMass) > MASS_EPS*dMass) {
#endif
			printf("ERROR id:%d lvl:%d:_pstRootSplit mass not cons.\n",
				   pst->idSelf,pst->iLvl);
			}
		/*
		 ** Now go on to DD of next levels.
		 */
		if (pst->nUpper > 1) 
			mdlReqService(pst->mdl,pst->idUpper,PST_DOMAINDECOMP,
						  vin,sizeof(*in));
		if (pst->nLower > 1) 
			pstDomainDecomp(pst->pstLower,vin,sizeof(*in),NULL,NULL);
		if (pst->nUpper > 1) 
			mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
		}

	if (pnOut) *pnOut = 0;

	pkdStopTimer(pst->plcl->pkd,6);
	sprintf(ach,"id: %d DD: %g %g %g  Weight: %g %g %g  SwapRej: %g %g %g  ColRej: %g %g %g  AO: %g %g %g  FS: %g %g %g\n",pst->mdl->idSelf,
			pkdGetTimer(pst->plcl->pkd,6),pkdGetSystemTimer(pst->plcl->pkd,6),pkdGetWallClockTimer(pst->plcl->pkd,6),
			pkdGetTimer(pst->plcl->pkd,7),pkdGetSystemTimer(pst->plcl->pkd,7),pkdGetWallClockTimer(pst->plcl->pkd,7),
			pkdGetTimer(pst->plcl->pkd,8),pkdGetSystemTimer(pst->plcl->pkd,8),pkdGetWallClockTimer(pst->plcl->pkd,8), 
			pkdGetTimer(pst->plcl->pkd,9),pkdGetSystemTimer(pst->plcl->pkd,9),pkdGetWallClockTimer(pst->plcl->pkd,9),
			pkdGetTimer(pst->plcl->pkd,4),pkdGetSystemTimer(pst->plcl->pkd,4),pkdGetWallClockTimer(pst->plcl->pkd,4), 
			pkdGetTimer(pst->plcl->pkd,5),pkdGetSystemTimer(pst->plcl->pkd,5),pkdGetWallClockTimer(pst->plcl->pkd,5) );
	mdlDiag(pst->mdl,ach);
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
			if (outBnd.bndTreeActive.fMin[j] < out->bndTreeActive.fMin[j]) {
				out->bndTreeActive.fMin[j] = outBnd.bndTreeActive.fMin[j];
				}
			if (outBnd.bndTreeActive.fMax[j] > out->bndTreeActive.fMax[j]) {
				out->bndTreeActive.fMax[j] = outBnd.bndTreeActive.fMax[j];
				}
			if (outBnd.bndBall.fMin[j] < out->bndBall.fMin[j]) {
				out->bndBall.fMin[j] = outBnd.bndBall.fMin[j];
				}
			if (outBnd.bndBall.fMax[j] > out->bndBall.fMax[j]) {
				out->bndBall.fMax[j] = outBnd.bndBall.fMax[j];
				}
			}
		}
	else {
		pkdCalcBound(plcl->pkd,&out->bnd,&out->bndActive,&out->bndTreeActive,&out->bndBall);
		}
	if (pnOut) *pnOut = sizeof(struct outCalcBound); 
	}

void pstGasWeight(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	
	assert(nIn == 0);
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_GASWEIGHT,vin,nIn);
		pstGasWeight(pst->pstLower,vin,nIn,NULL,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
		}
	else {
		pkdGasWeight(plcl->pkd);
		}
	if (pnOut) *pnOut = 0;
	}


/*
 ** Make sure that the local particles are split into active and inactive
 ** when passing pFlag != 0.
 ** pFlag == 0 => weight all particles.
 ** pFlag > 0 => weight active particles.
 ** pFlag < 0 => weight inactive particles.
 */
void pstWeight(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inWeight *in = vin;
	struct outWeight *out = vout;
	struct outWeight outWt;
	FLOAT fSplit,fLow,fHigh;
	int iSplitSide;

	assert(nIn == sizeof(struct inWeight));
	pkdStartTimer(plcl->pkd,7);
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
	pkdStopTimer(plcl->pkd,7);
	}


/*
 ** Weight request for splitting into iOrder order.
 */
void pstOrdWeight(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inOrdWeight *in = vin;
	struct outOrdWeight *out = vout;
	struct outOrdWeight outWt;

	assert(nIn == sizeof(struct inOrdWeight));
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_ORDWEIGHT,in,nIn);
		pstOrdWeight(pst->pstLower,in,nIn,out,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,&outWt,NULL);
		out->nLow += outWt.nLow;
		out->nHigh += outWt.nHigh;
		}
	else {
		if (in->ittr == 0) {
			/*
			 ** Initialize.
			 */
			plcl->iOrdSplit = in->iOrdSplit;
			plcl->iWtFrom = 0;
			plcl->iWtTo = pkdLocal(plcl->pkd)-1;
			}
		else {
			/*
			 ** Update the Weight Sums and use smaller weight region.
			 */
			if (in->iOrdSplit < plcl->iOrdSplit) {
				if (in->iSplitSide) plcl->iWtFrom = plcl->iPart;
				else plcl->iWtTo = plcl->iPart-1;
				}
			else {
				if (in->iSplitSide) plcl->iWtTo = plcl->iPart-1;
				else plcl->iWtFrom = plcl->iPart;
				}
			plcl->iOrdSplit = in->iOrdSplit;
			}
		plcl->iPart = pkdOrdWeight(plcl->pkd,in->iOrdSplit,in->iSplitSide,
								   plcl->iWtFrom,plcl->iWtTo,
								   &out->nLow,&out->nHigh);
		}
	if (pnOut) *pnOut = sizeof(struct outOrdWeight); 
	}


void pstFreeStore(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct outFreeStore *out = vout;
	int nLowerStore,nUpperStore;

	assert(nIn == 0);
	pkdStartTimer(plcl->pkd,4);
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
	pkdStopTimer(plcl->pkd,4);
	}


void pstColRejects(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inColRejects *in = vin;
	OREJ *pOutRej = vout;
	int nLower,nUpper,iUpper;
	
	assert(nIn == sizeof(struct inColRejects));
	pkdStartTimer(plcl->pkd,9);
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
	pkdStopTimer(plcl->pkd,9);
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
	    pOutRej->nRejects = pkdColOrdRejects(plcl->pkd,in->iOrdSplit,
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
	
	pkdStartTimer(plcl->pkd,8);
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
	pkdStopTimer(plcl->pkd,8);
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


void _pstOrdSplit(PST pst,int iMaxOrder)
{
	struct outFreeStore outFree;
	struct inOrdWeight inWt;
	struct outOrdWeight outWtLow,outWtHigh;
	int im=-1,imm,il,iu;
	int nLowerStore,nUpperStore,nLow=-1,nHigh=-1;
	struct inColOrdRejects inCol;
	OREJ *pLowerRej,*pUpperRej;
	int *pidSwap,iRet,nOut,ittr;

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
	 ** Find the correct iOrdSplit, such that all processors will
	 ** have close to the same number of particles in the end.
	 ** Start the ROOT finder based on balancing number of particles.
	 */
	il = 0;
	iu = iMaxOrder;
	imm = (il + iu)/2;
	ittr = 0;
	while (il < imm && imm < iu && ittr < MAX_ITTR) {
		im = imm;
		inWt.iOrdSplit = im;
		inWt.ittr = ittr;
		inWt.iSplitSide = 1;
		mdlReqService(pst->mdl,pst->idUpper,PST_ORDWEIGHT,&inWt,sizeof(inWt));
		inWt.iSplitSide = 0;
		pstOrdWeight(pst->pstLower,&inWt,sizeof(inWt),&outWtLow,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,&outWtHigh,NULL);
		/*
		 ** Add lower and Upper subsets weights and numbers
		 */
		nLow = outWtLow.nLow + outWtHigh.nLow;
		nHigh = outWtLow.nHigh + outWtHigh.nHigh;
		/*
		printf("ittr:%d l:%d u:%d\n",ittr,nLow,nHigh);
		*/
		if(nLow == 1 && nHigh == 1) /* break on trivial case */
		    break;
		else {		/* split on number */
		    if (nLow/(double)pst->nLower >
			nHigh/(double)pst->nUpper) iu = im;
		    else if (nLow/(double)pst->nLower <
			     nHigh/(double)pst->nUpper) il = im;
		    else break;
		    }
		imm = (il + iu)/2;
		++ittr;
		}
	assert(nLow <= nLowerStore);
	assert(nHigh <= nUpperStore);
	pst->iOrdSplit = im;
	/*
	 ** Collect rejects.
	 */
	pLowerRej = malloc(pst->nLower*sizeof(OREJ));
	assert(pLowerRej != NULL);
	pUpperRej = malloc(pst->nUpper*sizeof(OREJ));
	assert(pUpperRej != NULL);
	pidSwap = malloc(mdlThreads(pst->mdl)*sizeof(int));
	assert(pidSwap != NULL);
	inCol.iOrdSplit = pst->iOrdSplit;
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
	struct inDomainOrder *in = vin;

	assert(nIn == sizeof(struct inDomainOrder));
	if (pst->nLeaves > 1) {
		_pstOrdSplit(pst,in->iMaxOrder);
		/*
		 ** Now go on to Domain Order of next levels.
		 */
		if (pst->nUpper > 1) 
			mdlReqService(pst->mdl,pst->idUpper,PST_DOMAINORDER,in,nIn);
		if (pst->nLower > 1) 
			pstDomainOrder(pst->pstLower,in,nIn,NULL,NULL);
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
		pkdLocalOrder(plcl->pkd);
		}
	if (pnOut) *pnOut = 0;
	}


void pstActiveTypeOrder(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	int *pnActive = vout;
	struct inActiveTypeOrder *in = vin;

	assert(nIn == sizeof(*in));
	if (pst->nLeaves > 1) {
		int nActiveLeaf;
		mdlReqService(pst->mdl,pst->idUpper,PST_ACTIVETYPEORDER,vin,nIn);
		pstActiveTypeOrder(pst->pstLower,vin,nIn,vout,pnOut);
		mdlGetReply(pst->mdl,pst->idUpper,&nActiveLeaf,pnOut);
		*pnActive += nActiveLeaf;
		}
	else {
		*pnActive = pkdActiveTypeOrder(plcl->pkd, in->iTestMask );
		}
	if (pnOut) *pnOut = sizeof(int);
	}


void pstActiveOrder(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	int *pnActive = vout;

	assert(nIn == 0);
	pkdStartTimer(pst->plcl->pkd,5);
	if (pst->nLeaves > 1) {
		int nActiveLeaf;
		mdlReqService(pst->mdl,pst->idUpper,PST_ACTIVEORDER,NULL,0);
		pstActiveOrder(pst->pstLower,NULL,0,vout,pnOut);
		mdlGetReply(pst->mdl,pst->idUpper,&nActiveLeaf,pnOut);
		*pnActive += nActiveLeaf;
		}
	else {
		*pnActive = pkdActiveOrder(plcl->pkd);
		}
	if (pnOut) *pnOut = sizeof(int);
	pkdStopTimer(pst->plcl->pkd,5);
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
		pkdWriteTipsy(plcl->pkd,achOutFile,plcl->nWriteStart,
					  in->bStandard,in->dvFac,in->duTFac,in->iGasModel);
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
	double dOpen;
        
	assert(nIn == sizeof(struct inBuildTree));
	if (pst->nLeaves > 1) {
		pst->kdn.iDim = pst->iSplitDim;
		pst->kdn.fSplit = pst->fSplit;
		pst->kdn.pLower = -1;
		mdlReqService(pst->mdl,pst->idUpper,PST_BUILDTREE,in,nIn);
		pstBuildTree(pst->pstLower,in,nIn,&out1,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,&out2,NULL);
		/*
		 ** Combine to find cell mass, cm, softening.
		 ** This also combine the bounds to set pst->kdn.bnd!
		 */
		pkdCombine(&out1.kdn,&out2.kdn,&pst->kdn);
		/*
		 ** Check if both subtrees are null, i.e. pUpper < 0
		 ** if so set pUpper to -1;
		 */
		if (out1.kdn.pUpper < 0 && out2.kdn.pUpper < 0) {
			pst->kdn.pUpper = -1;
			}
		else {
			pst->kdn.pUpper = 1;
			}			
		if(in->bGravity) {
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
		}
	else {
		if (in->bBinary) {
			pkdBuildBinary(plcl->pkd,in->nBucket,in->iOpenType,in->dCrit,
				       in->iOrder,in->bTreeActiveOnly,
				       in->bGravity, &pst->kdn);
			}
		else {
			pkdBuildLocal(plcl->pkd,in->nBucket,in->iOpenType,in->dCrit,
				      in->iOrder,in->bTreeActiveOnly,
				      in->bGravity, &pst->kdn);
			}
		pst->kdn.pLower = pst->idSelf;
		if (pst->kdn.pUpper < 0) {
			pst->kdn.pUpper = -1;
			}
		else {
			pst->kdn.pUpper = 1;
			}	
		}
	/*
	 ** Calculated all cell properties, now pass up this cell info.
	 */
	out->kdn = pst->kdn;
	if (pnOut) *pnOut = sizeof(struct outBuildTree);
	}


void pstSmooth(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	struct inSmooth *in = vin;

	assert(nIn == sizeof(struct inSmooth));
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_SMOOTH,in,nIn);
		pstSmooth(pst->pstLower,in,nIn,NULL,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
		}
	else {
		LCL *plcl = pst->plcl;
		SMX smx;

		(&in->smf)->pkd = pst->plcl->pkd;
		smInitialize(&smx,plcl->pkd,&in->smf,in->nSmooth,in->bPeriodic,
					 in->bSymmetric,in->iSmoothType,1);
		smSmooth(smx,&in->smf);
		smFinish(smx,&in->smf);
		}
	if (pnOut) *pnOut = 0;
	}


void pstMarkSmooth(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	struct inMarkSmooth *in = vin;

	assert(nIn == sizeof(struct inMarkSmooth));
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_MARKSMOOTH,in,nIn);
		pstMarkSmooth(pst->pstLower,in,nIn,NULL,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
		}
	else {
		LCL *plcl = pst->plcl;
		SMX smx;

		(&in->smf)->pkd = pst->plcl->pkd;
		smInitialize(&smx,plcl->pkd,&in->smf,in->nSmooth,in->bPeriodic,
					 in->bSymmetric,in->iSmoothType,0);
		smMarkSmooth(smx,&in->smf,in->iMarkType);
		smFinish(smx,&in->smf);
		}
	if (pnOut) *pnOut = 0;
	}


void pstReSmooth(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	struct inReSmooth *in = vin;

	assert(nIn == sizeof(struct inReSmooth));
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_RESMOOTH,in,nIn);
		pstReSmooth(pst->pstLower,in,nIn,NULL,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
		}
	else {
		LCL *plcl = pst->plcl;
		SMX smx;

		(&in->smf)->pkd = pst->plcl->pkd;
		smInitialize(&smx,plcl->pkd,&in->smf,in->nSmooth,in->bPeriodic,
					 in->bSymmetric,in->iSmoothType,0);
		smReSmooth(smx,&in->smf);
		smFinish(smx,&in->smf);
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
	int j,bSunLower,bSunUpper;

	assert(nIn == sizeof(struct inGravity));
	if (pst->nLeaves > 1) {
		/*
		 ** Here we need to determine which domain is closest, or containing,
		 ** coordinate (0,0,0), the location of the Sun!
		 */
		if (pst->fSplit < 0) { 
			bSunLower = 0;
			bSunUpper = in->bDoSun;
			}
		else {
			bSunLower = in->bDoSun;
			bSunUpper = 0;
			}
		in->bDoSun = bSunUpper;
		mdlReqService(pst->mdl,pst->idUpper,PST_GRAVITY,in,nIn);
		in->bDoSun = bSunLower;
		pstGravity(pst->pstLower,in,nIn,out,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,&outUp,NULL);
		if (bSunUpper) {
			for (j=0;j<3;++j) out->aSun[j] = outUp.aSun[j];
			}
		out->dPartSum += outUp.dPartSum;
		out->dCellSum += outUp.dCellSum;
		out->dSoftSum += outUp.dSoftSum;
		out->dFlop += outUp.dFlop;
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
#ifdef SLIDING_PATCH
		pkdGravAll(plcl->pkd,in->nReps,in->bPeriodic,in->iOrder,in->bEwald,
				   in->iEwOrder,in->dEwCut,in->dEwhCut,in->bDoSun,out->aSun,
				   &out->nActive,&out->dPartSum,&out->dCellSum,&out->dSoftSum,
				   &cs,&out->dFlop,in->dOrbFreq,in->dTime);
#else
		pkdGravAll(plcl->pkd,in->nReps,in->bPeriodic,in->iOrder,in->bEwald,
				   in->iEwOrder,in->dEwCut,in->dEwhCut,in->bDoSun,out->aSun,
				   &out->nActive,&out->dPartSum,&out->dCellSum,&out->dSoftSum,
				   &cs,&out->dFlop);
#endif
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
		out->Eth += outE.Eth;
		}
	else {
		pkdCalcE(plcl->pkd,&out->T,&out->U,&out->Eth);
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
#ifdef SLIDING_PATCH
		pkdDrift(plcl->pkd,in->dDelta,in->fCenter,in->bPeriodic,in->bFandG,
				 in->fCentMass,in->dOrbFreq,in->dTime);
#else
		pkdDrift(plcl->pkd,in->dDelta,in->fCenter,in->bPeriodic,in->bFandG,
				 in->fCentMass);
#endif
		}
	if (pnOut) *pnOut = 0;
	}


void pstKick(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inKick *in = vin;
	struct outKick *out = vout;
	struct outKick outUp;

	assert(nIn == sizeof(struct inKick));

	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_KICK,in,nIn);
		pstKick(pst->pstLower,in,nIn,out,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,&outUp,NULL);

		out->SumTime += outUp.SumTime;
		out->nSum += outUp.nSum;
		if (outUp.MaxTime > out->MaxTime) out->MaxTime = outUp.MaxTime;
		}
	else {
		pkdKick(plcl->pkd,in->dvFacOne,in->dvFacTwo,
				in->dvPredFacOne,in->dvPredFacTwo,in->duDelta,in->duPredDelta,
				in->iGasModel,in->z,in->duDotLimit);
		out->Time = pkdGetTimer(plcl->pkd,1);
		out->MaxTime = out->Time;
		out->SumTime = out->Time;
		out->nSum = 1;
		}
	if (pnOut) *pnOut = sizeof(struct outKick);
	}


void pstReadCheck(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inReadCheck *in = vin;
	char achInFile[PST_FILENAME_SIZE];
	int nFileStart,nFileEnd,nFileTotal,nFileSplit,nStore;

	assert(nIn == sizeof(struct inReadCheck));
	nFileStart = in->nFileStart;
	nFileEnd = in->nFileEnd;
	nFileTotal = nFileEnd - nFileStart + 1;
	if (pst->nLeaves > 1) {
		nFileSplit = nFileStart + pst->nLower*(nFileTotal/pst->nLeaves);
		in->nFileStart = nFileSplit;
		mdlReqService(pst->mdl,pst->idUpper,PST_READCHECK,in,nIn);
		in->nFileStart = nFileStart;
		in->nFileEnd = nFileSplit - 1;
		pstReadCheck(pst->pstLower,in,nIn,vout,pnOut);
		in->nFileEnd = nFileEnd;
		mdlGetReply(pst->mdl,pst->idUpper,vout,pnOut);
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
		nStore = nFileTotal + (int)ceil(nFileTotal*in->fExtraStore);
		pkdInitialize(&plcl->pkd,pst->mdl,in->iOrder,nStore,plcl->nPstLvl,
					  in->fPeriod,in->nDark,in->nGas,in->nStar);
		pkdReadCheck(plcl->pkd,achInFile, in->iVersion,in->iOffset,
					 nFileStart,nFileTotal);
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


void pstSetWriteStart(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{	
	LCL *plcl = pst->plcl;
	struct inSetWriteStart *in = vin;
	int nWriteStart;

	assert(nIn == sizeof(struct inSetWriteStart));
	nWriteStart = in->nWriteStart;
	if (pst->nLeaves > 1) {
		in->nWriteStart = nWriteStart + pst->pstLower->nTotal;
		mdlReqService(pst->mdl,pst->idUpper,PST_SETWRITESTART,in,nIn);
		in->nWriteStart = nWriteStart;
		pstSetWriteStart(pst->pstLower,in,nIn,NULL,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
		}
	else {
		plcl->nWriteStart = nWriteStart;
		}
	if (pnOut) *pnOut = 0;
	}


void pstWriteCheck(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inWriteCheck *in = vin;
	char achOutFile[PST_FILENAME_SIZE];

	assert(nIn == sizeof(struct inWriteCheck));
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_WRITECHECK,in,nIn);
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
		pkdWriteCheck(plcl->pkd,achOutFile,in->iOffset,plcl->nWriteStart);
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
pstBallMax(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inBallMax *in = vin;
	
	assert(nIn == sizeof(*in));
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_BALLMAX,vin,nIn);
		pstBallMax(pst->pstLower,vin,nIn,NULL,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
		}
	else {
		pkdBallMax(plcl->pkd, in->iRung, in->bGreater, in->dhFac);
		}
	if (pnOut) *pnOut = 0;
	}

void
pstActiveRung(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inActiveRung *in = vin;
        int *pnActive = vout;
	
	assert(nIn == sizeof(*in));
	if (pst->nLeaves > 1) {
		int nActiveLeaf;
		mdlReqService(pst->mdl,pst->idUpper,PST_ACTIVERUNG,vin,nIn);
		pstActiveRung(pst->pstLower,vin,nIn,vout,pnOut);
		mdlGetReply(pst->mdl,pst->idUpper,&nActiveLeaf,pnOut);
		*pnActive += nActiveLeaf;
		}
	else {
		*pnActive = pkdActiveRung(plcl->pkd, in->iRung, in->bGreater);
		}
	if (pnOut) *pnOut = sizeof(int);
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


void pstDensityStep(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inDensityStep *in = vin;
	
	assert(nIn == sizeof(*in));
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_DENSITYSTEP,vin,nIn);
		pstDensityStep(pst->pstLower,vin,nIn,vout,pnOut);
		mdlGetReply(pst->mdl,pst->idUpper,vout,pnOut);
		}
	else {
		pkdDensityStep(plcl->pkd, in->dEta, in->dRhoFac);
		}
	if (pnOut) *pnOut = 0;
	}


void pstAccelStep(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inAccelStep *in = vin;

	assert(nIn == sizeof(struct inAccelStep));
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_ACCELSTEP,vin,nIn);
		pstAccelStep(pst->pstLower,vin,nIn,vout,pnOut);
		mdlGetReply(pst->mdl,pst->idUpper,vout,pnOut);
		}
	else {
		pkdAccelStep(plcl->pkd, in->dEta, in->dVelFac, in->dAccFac, in->bDoGravity,
			     in->bEpsVel, in->bSqrtPhi);
		}
	if (pnOut) *pnOut = 0;
	}

void pstDtToRung(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inDtToRung *in = vin;
	struct outDtToRung *out = vout;
	int iMaxRung;

	assert(nIn == sizeof(*in));
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_DTTORUNG,vin,nIn);
		pstDtToRung(pst->pstLower,vin,nIn,vout,pnOut);
		iMaxRung = out->iMaxRung;
		mdlGetReply(pst->mdl,pst->idUpper,vout,pnOut);
		if(iMaxRung > out->iMaxRung)
			out->iMaxRung = iMaxRung;
		}
	else {
		out->iMaxRung = pkdDtToRung(plcl->pkd, in->iRung,
									in->dDelta, in->iMaxRung, in->bAll);
		}
	if (pnOut) *pnOut = sizeof(*out);
	}

void pstInitDt(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inInitDt *in = vin;

	assert(nIn == sizeof(*in));
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_INITDT,vin,nIn);
		pstInitDt(pst->pstLower,vin,nIn,vout,pnOut);
		mdlGetReply(pst->mdl,pst->idUpper,vout,pnOut);
		}
	else {
		pkdInitDt(plcl->pkd, in->dDelta);
		}
	if (pnOut) *pnOut = 0;
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


void pstCoolVelocity(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inCoolVelocity *in = vin;

	assert(nIn == sizeof(struct inCoolVelocity));
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_COOLVELOCITY,in,nIn);
		pstCoolVelocity(pst->pstLower,in,nIn,NULL,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
		}
	else {
		pkdCoolVelocity(plcl->pkd,in->nSuperCool,in->dCoolFac,
						in->dCoolDens,in->dCoolMaxDens);
		}
	if (pnOut) *pnOut = 0;
	}


void pstResetTouchRung(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
        struct inResetTouchRung *in = vin;
        int *pnActive = vout;
	
	assert(nIn == sizeof(*in));

	if (pst->nLeaves > 1) {
		int nActiveLeaf;
		mdlReqService(pst->mdl,pst->idUpper,PST_RESETTOUCHRUNG,vin,nIn);
		pstResetTouchRung(pst->pstLower,vin,nIn,vout,pnOut);
		mdlGetReply(pst->mdl,pst->idUpper,&nActiveLeaf,pnOut);
		*pnActive += nActiveLeaf;
		}
	else {
		*pnActive = pkdResetTouchRung(plcl->pkd,in->iTestMask,in->iSetMask);
		}
	if (pnOut) *pnOut = sizeof(int);
	}


void pstActiveExactType(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inActiveType *in = vin;
	int *pnActive = vout;
	
	assert(nIn == sizeof(*in));

	if (pst->nLeaves > 1) {
		int nActiveLeaf;
		mdlReqService(pst->mdl,pst->idUpper,PST_ACTIVEEXACTTYPE,vin,nIn);
		pstActiveExactType(pst->pstLower,vin,nIn,vout,pnOut);
		mdlGetReply(pst->mdl,pst->idUpper,&nActiveLeaf,pnOut);
		*pnActive += nActiveLeaf;
		}
	else {
		*pnActive = pkdActiveExactType(plcl->pkd,in->iFilterMask,in->iTestMask,in->iSetMask);
		}
	if (pnOut) *pnOut = sizeof(int);
	}


void pstSetType(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inActiveType *in = vin;
	int *pnActive = vout;
	
	assert(nIn == sizeof(struct inActiveType));

	if (pst->nLeaves > 1) {
		int nActiveLeaf;
		mdlReqService(pst->mdl,pst->idUpper,PST_SETTYPE,vin,nIn);
		pstSetType(pst->pstLower,vin,nIn,vout,pnOut);
		mdlGetReply(pst->mdl,pst->idUpper,&nActiveLeaf,pnOut);
		*pnActive += nActiveLeaf;
		}
	else {
		*pnActive = pkdSetType(plcl->pkd,in->iTestMask,in->iSetMask);
		}
	if (pnOut) *pnOut = sizeof(int);
	}


void pstResetType(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inActiveType *in = vin;
	int *pnActive = vout;
	
	assert(nIn == sizeof(struct inActiveType));

	if (pst->nLeaves > 1) {
		int nActiveLeaf;
		mdlReqService(pst->mdl,pst->idUpper,PST_RESETTYPE,vin,nIn);
		pstResetType(pst->pstLower,vin,nIn,vout,pnOut);
		mdlGetReply(pst->mdl,pst->idUpper,&nActiveLeaf,pnOut);
		*pnActive += nActiveLeaf;
		}
	else {
		*pnActive = pkdResetType(plcl->pkd,in->iTestMask,in->iSetMask);
		}
	if (pnOut) *pnOut = sizeof(int);
	}

void pstCountType(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inActiveType *in = vin;
	int *pnActive = vout;
	
	assert(nIn == sizeof(*in));

	if (pst->nLeaves > 1) {
		int nActiveLeaf;
		mdlReqService(pst->mdl,pst->idUpper,PST_COUNTTYPE,vin,nIn);
		pstCountType(pst->pstLower,vin,nIn,vout,pnOut);
		mdlGetReply(pst->mdl,pst->idUpper,&nActiveLeaf,pnOut);
		*pnActive += nActiveLeaf;
		}
	else {
		*pnActive = pkdCountType(plcl->pkd,in->iFilterMask,in->iTestMask);
		}
	if (pnOut) *pnOut = sizeof(int);
	}


void pstActiveType(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inActiveType *in = vin;
	int *pnActive = vout;
	
	assert(nIn == sizeof(*in));

	if (pst->nLeaves > 1) {
		int nActiveLeaf;
		mdlReqService(pst->mdl,pst->idUpper,PST_ACTIVETYPE,vin,nIn);
		pstActiveType(pst->pstLower,vin,nIn,vout,pnOut);
		mdlGetReply(pst->mdl,pst->idUpper,&nActiveLeaf,pnOut);
		*pnActive += nActiveLeaf;
		}
	else {
		*pnActive = pkdActiveType(plcl->pkd,in->iTestMask,in->iSetMask);
		}
	if (pnOut) *pnOut = sizeof(int);
	}


void pstActiveMaskRung(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inActiveType *in = vin;
	int *pnActive = vout;
	
	assert(nIn == sizeof(*in));
	if (pst->nLeaves > 1) {
		int nActiveLeaf;
		mdlReqService(pst->mdl,pst->idUpper,PST_ACTIVEMASKRUNG,vin,nIn);
		pstActiveMaskRung(pst->pstLower,vin,nIn,vout,pnOut);
		mdlGetReply(pst->mdl,pst->idUpper,&nActiveLeaf,pnOut);
		*pnActive += nActiveLeaf;
		}
	else {
		*pnActive = pkdActiveMaskRung(plcl->pkd,in->iSetMask,in->iRung,in->bGreater);
		}
	if (pnOut) *pnOut = sizeof(int);
	}


void pstActiveTypeRung(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inActiveType *in = vin;
	int *pnActive = vout;
	
	assert(nIn == sizeof(*in));
	if (pst->nLeaves > 1) {
		int nActiveLeaf;
		mdlReqService(pst->mdl,pst->idUpper,PST_ACTIVETYPERUNG,vin,nIn);
		pstActiveTypeRung(pst->pstLower,vin,nIn,vout,pnOut);
		mdlGetReply(pst->mdl,pst->idUpper,&nActiveLeaf,pnOut);
		*pnActive += nActiveLeaf;
		}
	else {
		*pnActive = pkdActiveTypeRung(plcl->pkd,in->iTestMask,in->iSetMask,
									  in->iRung,in->bGreater);
		}
	if (pnOut) *pnOut = sizeof(int);
	}


void pstSetParticleTypes(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inSetParticleTypes *in = vin;
	
	assert(nIn == sizeof(*in));
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_SETPARTICLETYPES,vin,nIn);
		pstSetParticleTypes(pst->pstLower,vin,nIn,NULL,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
		}
	else {
		pkdSetParticleTypes(plcl->pkd,in->nSuperCool);
		}
	if (pnOut) *pnOut = 0;
	}

void pstGrowMass(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inGrowMass *in = vin;

	assert(nIn == sizeof(struct inGrowMass));
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_GROWMASS,in,nIn);
		pstGrowMass(pst->pstLower,in,nIn,NULL,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
		}
	else {
		pkdGrowMass(plcl->pkd,in->nGrowMass,in->dDeltaM);
		}
	if (pnOut) *pnOut = 0;
	}

void pstInitAccel(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	
	assert(nIn == 0);
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_INITACCEL,vin,nIn);
		pstInitAccel(pst->pstLower,vin,nIn,NULL,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
		}
	else {
		pkdInitAccel(plcl->pkd);
		}
	if (pnOut) *pnOut = 0;
	}


void pstGravExternal(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inGravExternal *in = vin;

	assert(nIn == sizeof(struct inGravExternal));
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_GRAVEXTERNAL,in,nIn);
		pstGravExternal(pst->pstLower,in,nIn,NULL,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
		}
	else {
		if (in->bIndirect) {
			pkdSunIndirect(plcl->pkd,in->aSun,in->bDoSun,in->dSunMass);
			}
		if (in->bLogHalo) {
			pkdLogHalo(plcl->pkd);
			}
		if (in->bHernquistSpheroid) {
			pkdHernquistSpheroid(plcl->pkd);
			}
		if (in->bMiyamotoDisk) {
			pkdMiyamotoDisk(plcl->pkd);
			}
#ifdef SLIDING_PATCH
		if (in->bPatch) {
			pkdPatch(plcl->pkd,in->dOrbFreq,in->dOrbFreqZ2);
			}
#endif
		}
	if (pnOut) *pnOut = 0;
	}


#ifdef GASOLINE

void pstUpdateuDot(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inUpdateuDot *in = vin;
	struct outUpdateuDot *out = vout;
	struct outUpdateuDot outUp;

	assert(nIn == sizeof(struct inUpdateuDot));

	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_UPDATEUDOT,in,nIn);
		pstUpdateuDot(pst->pstLower,in,nIn,out,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,&outUp,NULL);

		out->SumTime += outUp.SumTime;
		out->nSum += outUp.nSum;
		if (outUp.MaxTime > out->MaxTime) out->MaxTime = outUp.MaxTime;
		}
	else {
		pkdUpdateuDot(plcl->pkd,in->duDelta,in->z,in->iGasModel,in->bUpdateY);
		out->Time = pkdGetTimer(plcl->pkd,1);
		out->MaxTime = out->Time;
		out->SumTime = out->Time;
		out->nSum = 1;
		}
	if (pnOut) *pnOut = sizeof(struct outUpdateuDot);
	}

void pstGetGasPressure(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inGetGasPressure *in = vin;
	
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_GETGASPRESSURE,vin,nIn);
		pstGetGasPressure(pst->pstLower,vin,nIn,NULL,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
		}
	else {
		switch( in->iGasModel ) {
		case GASMODEL_ADIABATIC: 
		case GASMODEL_COOLING:
		case GASMODEL_COOLING_NONEQM:
			pkdAdiabaticGasPressure(plcl->pkd, in->gammam1,in->gamma);
			break;
		case GASMODEL_ISOTHERMAL:
			assert(0);
			break;
		case GASMODEL_GLASS:
#ifdef GLASS		  
			pkdGlassGasPressure(plcl->pkd, in);
#else
			assert(0);
#endif
			break;
			}
		}
	if (pnOut) *pnOut = 0;
	}


void pstInitEnergy(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inInitEnergy *in = vin;
	
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_INITENERGY,vin,nIn);
		pstInitEnergy(pst->pstLower,vin,nIn,NULL,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
		}
	else {
		pkdInitEnergy(plcl->pkd,in->dTuFac,in->z);
		}
	if (pnOut) *pnOut = 0;
	}


void pstKickVpred(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inKickVpred *in = vin;
	struct outKick *out = vout;
	struct outKick outUp;

	assert(nIn == sizeof(struct inKickVpred));
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_KICKVPRED,in,nIn);
		pstKickVpred(pst->pstLower,in,nIn,out,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,&outUp,NULL);
		}
	else {
		pkdKickVpred(plcl->pkd,in->dvFacOne,in->dvFacTwo,in->duDelta,
					 in->iGasModel,in->z,in->duDotLimit);
		out->Time = pkdGetTimer(plcl->pkd,1);
		out->MaxTime = out->Time;
		out->SumTime = out->Time;
		out->nSum = 1;
		}
	if (pnOut) *pnOut = sizeof(struct outKick);
	}

void pstKickRhopred(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inKickRhopred *in = vin;

	assert(nIn == sizeof(struct inKickRhopred));
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_KICKRHOPRED,in,nIn);
		pstKickRhopred(pst->pstLower,in,nIn,NULL,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
		}
	else {
		pkdKickRhopred(plcl->pkd,in->dHubbFac,in->dDelta);
		}
	if (pnOut) *pnOut = 0;
	}

void
pstSphCurrRung(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inSphCurrRung *in = vin;
	struct outSphCurrRung *out = vout;
	int iCurrent;
	
	assert(nIn == sizeof(*in));
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_SPHCURRRUNG,vin,nIn);
		pstSphCurrRung(pst->pstLower,vin,nIn,vout,pnOut);
		iCurrent = out->iCurrent;
		mdlGetReply(pst->mdl,pst->idUpper,vout,pnOut);
		if(iCurrent)
		    out->iCurrent = iCurrent;
		}
	else {
		out->iCurrent = pkdSphCurrRung(plcl->pkd, in->iRung, in->bGreater);
		}
	if (pnOut) *pnOut = sizeof(*out);
	}

void pstSphStep(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inSphStep *in = vin;

	assert(nIn == sizeof(struct inSphStep));
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_SPHSTEP,in,nIn);
		pstSphStep(pst->pstLower,in,nIn,NULL,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
		}
	else {
		pkdSphStep(plcl->pkd,in->dCosmoFac,in->dEtaCourant,in->dEtauDot);
		}
	if (pnOut) *pnOut = 0;
	}

void pstSphViscosityLimiter(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inSphViscosityLimiter *in = vin;

	assert(nIn == sizeof(struct inSphViscosityLimiter));
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_SPHVISCOSITYLIMITER,in,nIn);
		pstSphViscosityLimiter(pst->pstLower,in,nIn,NULL,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
		}
	else {
		pkdSphViscosityLimiter(plcl->pkd,in->bOn);
		}
	if (pnOut) *pnOut = 0;
	}

void pstInitCooling(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inInitCooling *in = vin;

	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_INITCOOLING,in,nIn);
		pstInitCooling(pst->pstLower,in,nIn,NULL,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
		}
	else {
		(plcl->pkd->cl).mdl = plcl->pkd->mdl;
		clInitConstants(&(plcl->pkd->cl),in->dGmPerCcUnit,in->dComovingGmPerCcUnit,
						in->dErgPerGmUnit,in->dSecUnit,in->dMassFracHelium);
		clInitRatesTable(&(plcl->pkd->cl),in->Tmin,in->Tmax,in->nTable);
		clRatesRedshift(&(plcl->pkd->cl),in->z);
		}
	if (pnOut) *pnOut = 0;
	}

void pstInitUV(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;

	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_INITUV,vin,nIn);
		pstInitUV(pst->pstLower,vin,nIn,NULL,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
		}
	else {
		clInitUV(&(plcl->pkd->cl),vin,nIn/sizeof(UVSPECTRUM));
		}
	if (pnOut) *pnOut = 0;
	}

void
pstDensCheck(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inDensCheck *in = vin;
	struct outDensCheck *out = vout, tmp;
	
	assert(nIn == sizeof(*in));
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_DENSCHECK,vin,nIn);
		pstDensCheck(pst->pstLower,vin,nIn,vout,pnOut);
		mdlGetReply(pst->mdl,pst->idUpper,&tmp,pnOut);
		if (tmp.dMaxDensError > out->dMaxDensError)
			out->dMaxDensError = tmp.dMaxDensError;
		out->dAvgDensError = (tmp.dAvgDensError*tmp.nTotal + out->dAvgDensError*out->nTotal)/(tmp.nTotal+out->nTotal);
		out->nTotal += tmp.nTotal;
		out->nError += tmp.nError;
		}
	else {
		pkdDensCheck(plcl->pkd,in->iRung,in->bGreater,in->iMeasure,out);
		}
	if (pnOut) *pnOut = sizeof(struct outDensCheck);
	}

#endif

#ifdef GLASS
void pstRandomVelocities(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inRandomVelocities *in = vin;

	assert(nIn == sizeof(struct inRandomVelocities));
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_RANDOMVELOCITIES,in,nIn);
		pstRandomVelocities(pst->pstLower,in,nIn,NULL,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
		}
	else {
		pkdRandomVelocities(plcl->pkd,in->dMaxVelocityL,in->dMaxVelocityR);
		}
	if (pnOut) *pnOut = 0;
	}
#endif

void
pstColNParts(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
    LCL *plcl = pst->plcl;
    struct outColNParts *out = vout;
    struct outColNParts *ptmp;
    int i;
    
    for(i=0;i<pst->mdl->nThreads;i++)
		out[i].nNew = -1;
    if(pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_COLNPARTS,vin,nIn);
		pstColNParts(pst->pstLower,vin,nIn,vout,pnOut);
		ptmp = malloc(pst->mdl->nThreads*sizeof(*ptmp));
		assert(ptmp != NULL);
		mdlGetReply(pst->mdl,pst->idUpper,ptmp,pnOut);
		for(i = 0; i < pst->mdl->nThreads; i++) {
			if(ptmp[i].nNew != -1)
				out[i] = ptmp[i];
			}
		free(ptmp);
		}
    else {
		pkdColNParts(plcl->pkd, &out[pst->idSelf].nNew,
					 &out[pst->idSelf].nDeltaGas,
					 &out[pst->idSelf].nDeltaDark,
					 &out[pst->idSelf].nDeltaStar);
		}
    if(pnOut) *pnOut = pst->mdl->nThreads*sizeof(*out);
    }

void
pstNewOrder(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
    LCL *plcl = pst->plcl;
    int *in = vin;
    
    if(pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_NEWORDER,vin,nIn);
		pstNewOrder(pst->pstLower, vin, nIn, vout, pnOut);
		mdlGetReply(pst->mdl, pst->idUpper, vout, pnOut);
		}
    else {
		pkdNewOrder(plcl->pkd, in[pst->idSelf]);
		}
    if(pnOut) *pnOut = 0;
    }

void
pstSetNParts(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
    struct inSetNParts *in = vin;
    
    if(pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_SETNPARTS,vin,nIn);
		pstSetNParts(pst->pstLower, vin, nIn, vout, pnOut);
		mdlGetReply(pst->mdl, pst->idUpper, vout, pnOut);
		}
    else {
		pkdSetNParts(pst->plcl->pkd, in->nGas, in->nDark, in->nStar,
					 in->nMaxOrderGas, in->nMaxOrderDark);
		}
    if(pnOut) *pnOut = 0;
    }

#ifdef COLLISIONS

void
pstNumRejects(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct outNumRejects *out = vout;
	int nRej;

	assert(nIn == 0);
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_NUMREJECTS,NULL,0);
		pstNumRejects(pst->pstLower,NULL,0,vout,pnOut);
		nRej = out->nRej;
		mdlGetReply(pst->mdl,pst->idUpper,vout,pnOut);
		out->nRej += nRej;
		}
	else {
		out->nRej = pkdNumRejects(plcl->pkd);
		}
	if (pnOut) *pnOut = sizeof(*out);
	}

void
pstReadSS(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inReadSS *in = vin;
	int nFileStart,nFileEnd,nFileTotal,nFileSplit,nStore;
	char achInFile[PST_FILENAME_SIZE];

	assert(nIn == sizeof(struct inReadSS));
	nFileStart = in->nFileStart;
	nFileEnd = in->nFileEnd;
	nFileTotal = nFileEnd - nFileStart + 1;
	if (pst->nLeaves > 1) {
		nFileSplit = nFileStart + pst->nLower*(nFileTotal/pst->nLeaves);
		in->nFileStart = nFileSplit;
		mdlReqService(pst->mdl,pst->idUpper,PST_READSS,vin,nIn);
		in->nFileStart = nFileStart;
		in->nFileEnd = nFileSplit - 1;
		pstReadSS(pst->pstLower,vin,nIn,NULL,NULL);
		in->nFileEnd = nFileEnd;
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
		nStore = nFileTotal + (int)ceil(nFileTotal*in->fExtraStore);
		pkdInitialize(&plcl->pkd,pst->mdl,in->iOrder,nStore,plcl->nPstLvl,
					  in->fPeriod,in->nDark,in->nGas,in->nStar);
		pkdReadSS(plcl->pkd,achInFile,nFileStart,nFileTotal);
		}
	if (pnOut) *pnOut = 0;
	}

void
pstWriteSS(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inWriteSS *in = vin;
	char achOutFile[PST_FILENAME_SIZE];

	assert(nIn == sizeof(struct inWriteSS));
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_WRITESS,vin,nIn);
		pstWriteSS(pst->pstLower,vin,nIn,NULL,NULL);
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
		pkdWriteSS(plcl->pkd,achOutFile,plcl->nWriteStart);
		}
	if (pnOut) *pnOut = 0;
	}

void
pstCalcHill(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inCalcHill *in = vin;

	assert(nIn == sizeof(struct inCalcHill));
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_CALCHILL,vin,nIn);
		pstCalcHill(pst->pstLower,vin,nIn,NULL,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
		}
	else {
		pkdCalcHill(plcl->pkd,in->dCentMass);
		}
	if (pnOut) *pnOut = 0;
	}

void
pstHelioStep(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inHelioStep *in = vin;

	assert(nIn == sizeof(struct inHelioStep));
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_HELIOSTEP,vin,nIn);
		pstHelioStep(pst->pstLower,vin,nIn,NULL,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
		}
	else {
		pkdHelioStep(plcl->pkd,in->dEta);
		}
	if (pnOut) *pnOut = 0;
	}

void
pstKickUnifGrav(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inKickUnifGrav *in = vin;

	assert(nIn == sizeof(struct inKickUnifGrav));
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_KICKUNIFGRAV,in,nIn);
		pstKickUnifGrav(pst->pstLower,in,nIn,NULL,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
		}
	else {
		pkdKickUnifGrav(plcl->pkd,in->dvx,in->dvy,in->dvz);
		}
	if (pnOut) *pnOut = 0;
	}

void
pstNextEncounter(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct outNextEncounter local,*out = vout;

	assert(nIn == 0);
	/*DEBUG MUST INITIALIZE HERE! (cf. pstNextCollision())*/
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_NEXTENCOUNTER,NULL,0);
		pstNextEncounter(pst->pstLower,NULL,0,vout,NULL);
		local = *out;
		mdlGetReply(pst->mdl,pst->idUpper,vout,NULL);
		if (local.dt < out->dt) out->dt = local.dt;
		}
	else {
		pkdNextEncounter(plcl->pkd,&out->dt);
		}
	if (pnOut) *pnOut = sizeof(*out);
	}

void
pstMarkEncounters(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inMarkEncounters *in = vin;
	
	assert(nIn == sizeof(*in));
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_MARKENCOUNTERS,vin,nIn);
		pstMarkEncounters(pst->pstLower,vin,nIn,NULL,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
		}
	else {
		pkdMarkEncounters(plcl->pkd,in->dt);
		}
	if (pnOut) *pnOut = 0;
	}

void
pstNextCollision(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct outNextCollision local,*out = vout;

	assert(nIn == 0);
	local.dt = out->dt = DBL_MAX; /* trust me */
	local.iOrder1 = local.iOrder2 = out->iOrder1 = out->iOrder2 = -1;
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_NEXTCOLLISION,NULL,0);
		pstNextCollision(pst->pstLower,NULL,0,vout,NULL);
		local = *out;
		mdlGetReply(pst->mdl,pst->idUpper,vout,NULL);
		if (local.dt < out->dt) *out = local;
		}
	else {
		pkdNextCollision(plcl->pkd,&out->dt,&out->iOrder1,&out->iOrder2);
		}
	if (pnOut) *pnOut = sizeof(*out);
	}

void
pstGetColliderInfo(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inGetColliderInfo *in = vin;
	struct outGetColliderInfo local,*out = vout;

	assert(nIn == sizeof(*in));
	local.Collider.id.iOrder = out->Collider.id.iOrder = -1;
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_GETCOLLIDERINFO,vin,nIn);
		pstGetColliderInfo(pst->pstLower,vin,nIn,vout,NULL);
		if (out->Collider.id.iOrder == in->iOrder) local = *out;
		mdlGetReply(pst->mdl,pst->idUpper,vout,NULL);
		if (local.Collider.id.iOrder == in->iOrder) *out = local;
		}
	else {
		pkdGetColliderInfo(plcl->pkd,in->iOrder,&out->Collider);
		}
	if (pnOut) *pnOut = sizeof(*out);
	}

void
pstDoCollision(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inDoCollision *in = vin;
	struct outDoCollision local,*out = vout;

	assert(nIn == sizeof(*in));
	local.nOut = out->nOut = 0;
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_DOCOLLISION,vin,nIn);
		pstDoCollision(pst->pstLower,vin,nIn,vout,NULL);
		if (out->nOut) local = *out;
		mdlGetReply(pst->mdl,pst->idUpper,vout,NULL);
		if (local.nOut) *out = local;
		}
	else {
		PKD pkd = plcl->pkd;
		if (in->Collider1.id.iPid == pkd->idSelf ||
			in->Collider2.id.iPid == pkd->idSelf) {
#ifdef SLIDING_PATCH
			pkdDoCollision(pkd,in->dTime,in->dt,&in->Collider1,&in->Collider2,
						   in->bPeriodic,&in->CP,in->dOrbFreq,
						   &out->dImpactEnergy,&out->iOutcome,out->Out,
						   &out->nOut);
#else
			pkdDoCollision(pkd,in->dTime,in->dt,&in->Collider1,&in->Collider2,
						   in->bPeriodic,&in->CP,&out->dImpactEnergy,
						   &out->iOutcome,out->Out,&out->nOut);
#endif
			}
		}
	if (pnOut) *pnOut = sizeof(*out);
	}

void
pstResetColliders(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inResetColliders *in = vin;
	
	assert(nIn == sizeof(*in));
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_RESETCOLLIDERS,vin,nIn);
		pstResetColliders(pst->pstLower,vin,nIn,NULL,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
		}
	else {
		pkdResetColliders(plcl->pkd,in->iOrder1,in->iOrder2);
		}
	if (pnOut) *pnOut = 0;
	}

void
pstQQCalcBound(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct outCalcBound *out = vout;
	struct outCalcBound outBnd;
	int j;
	
	assert(nIn == 0);
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_QQCALCBOUND,NULL,0);
		pstQQCalcBound(pst->pstLower,NULL,0,out,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,&outBnd,NULL);
		for (j=0;j<2;++j) { /* 2D tree */
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
		pkdQQCalcBound(plcl->pkd,&out->bnd,&out->bndActive);
		}
	if (pnOut) *pnOut = sizeof(struct outCalcBound); 
	}

void
pstQQDomainDecomp(PST pst,void *vin,int nIn,void *vout,int *pnOut)
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
		pstQQCalcBound(pst,NULL,0,&outBnd,&nOut);
		pst->bndActive = outBnd.bndActive;
		pst->bnd = outBnd.bnd;
		/*
		 ** Next determine the longest axis based on the active bounds.
		 */
		d = 0;
		if(pst->bndActive.fMax[0] > pst->bndActive.fMin[0]) {
		    for (j=1;j<2;++j) {	/* 2D tree */
			if (pst->bndActive.fMax[j]-pst->bndActive.fMin[j] > 
				pst->bndActive.fMax[d]-pst->bndActive.fMin[d]) d = j;
			    }
		    	}
		else {		/* no active particles */
		    for (j=1;j<2;++j) {	/* 2D tree */
			if (pst->bnd.fMax[j]-pst->bnd.fMin[j] > 
				pst->bnd.fMax[d]-pst->bnd.fMin[d]) d = j;
			    }
		    	}

		pst->iSplitDim = d+3; /* Scary! */
		
		pstMassCheck(pst,NULL,0,&outMass,NULL);
		dMass = outMass.dMass;
		_pstRootSplit(pst,d+3,dMass,1); /* Scary! */
		pstMassCheck(pst,NULL,0,&outMass,NULL);
		if (fabs(dMass - outMass.dMass) > MASS_EPS*dMass) {
			printf("ERROR id:%d lvl:%d:_pstQQRootSplit mass not cons.\n",
				   pst->idSelf,pst->iLvl);
			}
		/*
		 ** Now go on to DD of next levels.
		 */
		if (pst->nUpper > 1) 
			mdlReqService(pst->mdl,pst->idUpper,PST_QQDOMAINDECOMP,NULL,0);
		if (pst->nLower > 1) 
			pstQQDomainDecomp(pst->pstLower,NULL,0,NULL,NULL);
		if (pst->nUpper > 1) 
			mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
		}
	if (pnOut) *pnOut = 0;
	}

void
pstQQBuildTree(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inBuildTree *in = vin;
	struct outBuildTree *out = vout;
	struct outBuildTree out1,out2;
	
	assert(nIn == sizeof(struct inBuildTree));
	if (pst->nLeaves > 1) {
		pst->kdn.iDim = pst->iSplitDim;
		pst->kdn.fSplit = pst->fSplit;
		pst->kdn.pLower = -1;
		pst->kdn.pUpper = 1;
		mdlReqService(pst->mdl,pst->idUpper,PST_QQBUILDTREE,in,nIn);
		pstQQBuildTree(pst->pstLower,in,nIn,&out1,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,&out2,NULL);
		/*
		 ** Combine to find cell mass, cm, softening.
		 ** This also combine the bounds to set pst->kdn.bnd!
		 */
		pkdCombine(&out1.kdn,&out2.kdn,&pst->kdn);
		}
	else {
		pkdQQBuild(plcl->pkd,in->nBucket,in->bActiveOnly,&pst->kdn);
		pst->kdn.pLower = pst->idSelf;
		pst->kdn.pUpper = 1;
		}
	/*
	 ** Calculated all cell properties, now pass up this cell info.
	 */
	out->kdn = pst->kdn;
	if (pnOut) *pnOut = sizeof(struct outBuildTree);
	}

void
pstQQSmooth(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	struct inSmooth *in = vin;

	assert(nIn == sizeof(struct inSmooth));
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_QQSMOOTH,in,nIn);
		pstQQSmooth(pst->pstLower,in,nIn,NULL,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
		}
	else {
		LCL *plcl = pst->plcl;
		SMX smx;

		smInitialize(&smx,plcl->pkd,&in->smf,in->nSmooth,in->bPeriodic,
					 in->bSymmetric,in->iSmoothType,0);
		smQQSmooth(smx,&in->smf);
		smFinish(smx,&in->smf);
		}
	if (pnOut) *pnOut = 0;
	}

#endif /* COLLISIONS */

#ifdef SLIDING_PATCH

void pstKickVpred(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inKickVpred *in = vin;

	assert(nIn == sizeof(struct inKickVpred));
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_KICKVPRED,in,nIn);
		pstKickVpred(pst->pstLower,in,nIn,NULL,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
		}
	else {
		pkdKickVpred(plcl->pkd,in->dvFacOne,in->dvFacTwo);
		}
	if (pnOut) *pnOut = 0;
	}

#endif

void 
pstClearTimer(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	struct inClearTimer *in = vin;

	assert(nIn == sizeof(struct inClearTimer));
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_CLEARTIMER,in,nIn);
		pstClearTimer(pst->pstLower,in,nIn,NULL,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
		}
	else {
		pkdClearTimer(pst->plcl->pkd,in->iTimer);
		}
	if (pnOut) *pnOut = 0;
	}
