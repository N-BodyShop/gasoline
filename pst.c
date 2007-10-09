#include <math.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <assert.h>
#include "mdl.h"
#include "pst.h"
#include "pkd.h"	
#include "outtype.h"
#include "smooth.h"

void
pstAddServices(PST pst,MDL mdl)
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
	mdlAddService(mdl,PST_RUNGDDWEIGHT,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstRungDDWeight,
				  sizeof(struct inRungDDWeight),0);
	mdlAddService(mdl,PST_WEIGHT,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstWeight,
				  sizeof(struct inWeight),sizeof(struct outWeight));
	mdlAddService(mdl,PST_WEIGHTWRAP,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstWeightWrap,
				  sizeof(struct inWeightWrap),sizeof(struct outWeight));
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
	mdlAddService(mdl,PST_MOVEPART,pst,
		      (void (*)(void *,void *,int,void *,int *)) pstMoveParticle,
		      sizeof(struct inMoveParticle),0);
	
	mdlAddService(mdl,PST_OUTARRAY,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstOutArray,
				  sizeof(struct inOutArray),0);
	mdlAddService(mdl,PST_OUTNCVECTOR,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstOutNCVector,
				  sizeof(struct inOutput),sizeof(struct outNC));
	mdlAddService(mdl,PST_OUTVECTOR,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstOutVector,
				  sizeof(struct inOutput),0);
	mdlAddService(mdl,PST_WRITETIPSY,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstWriteTipsy,
				  sizeof(struct inWriteTipsy),0);
	mdlAddService(mdl,PST_BUILDTREE,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstBuildTree,
				  sizeof(struct inBuildTree),sizeof(struct outBuildTree));
	mdlAddService(mdl,PST_SMOOTH,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstSmooth,
				  sizeof(struct inSmooth),sizeof(struct outSmooth));
	mdlAddService(mdl,PST_GRAVITY,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstGravity,
				  sizeof(struct inGravity),sizeof(struct outGravity));
	mdlAddService(mdl,PST_GRAVEXTERNAL,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstGravExternal,
				  sizeof(struct inGravExternal),sizeof(struct outGravExternal));
	mdlAddService(mdl,PST_CALCEANDL,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstCalcEandL,
				  0,sizeof(struct outCalcEandL));
	mdlAddService(mdl,PST_CALCEANDLEXT,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstCalcEandLExt,
				  sizeof(struct inCalcEandLExt),
				  sizeof(struct outCalcEandLExt));
	mdlAddService(mdl,PST_DRIFT,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstDrift,
				  sizeof(struct inDrift),0);
	mdlAddService(mdl,PST_KICK,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstKick,
				  sizeof(struct inKick),sizeof(struct outKick));
	mdlAddService(mdl,PST_KICKPATCH,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstKickPatch,
				  sizeof(struct inKickPatch),sizeof(struct outKick));
	mdlAddService(mdl,PST_READCHECK,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstReadCheck,
				  sizeof(struct inReadCheck), 0);
	mdlAddService(mdl,PST_WRITECHECK,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstWriteCheck,
				  sizeof(struct inWriteCheck),0);
	mdlAddService(mdl,PST_SETSOFT,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstSetSoft,
				  sizeof(struct inSetSoft),0);
#ifdef CHANGESOFT
	mdlAddService(mdl,PST_PHYSICALSOFT,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstPhysicalSoft,
				  sizeof(struct inPhysicalSoft),0);
	mdlAddService(mdl,PST_PREVARIABLESOFT,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstPreVariableSoft,
				  sizeof(struct inPreVariableSoft),0);
	mdlAddService(mdl,PST_POSTVARIABLESOFT,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstPostVariableSoft,
				  sizeof(struct inPostVariableSoft),0);
#endif
	mdlAddService(mdl,PST_SETTOTAL,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstSetTotal,
				  0,sizeof(struct outSetTotal));
	mdlAddService(mdl,PST_SETTOTALS,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstSetTotals,
				  0,sizeof(struct outSetTotals));
	mdlAddService(mdl,PST_SETWRITESTART,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstSetWriteStart,
				  sizeof(struct inSetWriteStart),0);
	mdlAddService(mdl,PST_SETNCWRITESTART,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstSetNCWriteStart,
				  sizeof(struct inSetNCWriteStart),0);
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
	mdlAddService(mdl,PST_MASSMETALSENERGYCHECK,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstMassMetalsEnergyCheck,
				  0,sizeof(struct outMassMetalsEnergyCheck));
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
	mdlAddService(mdl,PST_GRAVSTEP,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstGravStep,
				  sizeof(struct inGravStep), 0);
	mdlAddService(mdl,PST_ACCELSTEP,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstAccelStep,
				  sizeof(struct inAccelStep), 0);
	mdlAddService(mdl,PST_COOLVELOCITY,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstCoolVelocity,
				  sizeof(struct inCoolVelocity),0);
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
	mdlAddService(mdl,PST_SETTYPEFROMFILE,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstSetTypeFromFile,
				  sizeof(struct inSetTypeFromFile),sizeof(struct outSetTypeFromFile));
	mdlAddService(mdl,PST_SETPARTICLETYPES,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstSetParticleTypes,
				  sizeof(struct inSetParticleTypes),0);
	mdlAddService(mdl,PST_SOUGHTPARTICLELIST,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstSoughtParticleList,
				  sizeof(struct inSoughtParticleList),sizeof(struct inoutParticleList));
	mdlAddService(mdl,PST_COOLUSINGPARTICLELIST,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstCoolUsingParticleList,
				  sizeof(struct inoutParticleList),0);
	mdlAddService(mdl,PST_GROWMASS,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstGrowMass,
				  sizeof(struct inGrowMass),0);
	mdlAddService(mdl,PST_MARKSMOOTH,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstMarkSmooth,
				  sizeof(struct inMarkSmooth),0);
	mdlAddService(mdl,PST_RESMOOTH,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstReSmooth,
				  sizeof(struct inReSmooth),sizeof(struct outReSmooth));
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
	mdlAddService(mdl,PST_UPDATESHOCKTRACKER,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstUpdateShockTracker,
				  sizeof(struct inUpdateShockTracker),0);
	mdlAddService(mdl,PST_SPHCURRRUNG,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstSphCurrRung,
				  sizeof(struct inSphCurrRung),sizeof(struct outSphCurrRung));
	mdlAddService(mdl,PST_GETGASPRESSURE,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstGetGasPressure,
				  sizeof(struct inGetGasPressure),0);
	mdlAddService(mdl,PST_GETDENSITYU,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstGetDensityU,
				  sizeof(double),0);
	mdlAddService(mdl,PST_LOWERSOUNDSPEED,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstLowerSoundSpeed,
				  sizeof(struct inLowerSoundSpeed),0);
	mdlAddService(mdl,PST_KICKRHOPRED,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstKickRhopred, 
				  sizeof(struct inKickRhopred),0);
	mdlAddService(mdl,PST_SPHSTEP,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstSphStep,
				  sizeof(struct inSphStep),sizeof(double));
	mdlAddService(mdl,PST_SINKSTEP,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstSinkStep,
				  sizeof(struct inSinkStep),0);
	mdlAddService(mdl,PST_SPHVISCOSITYLIMITER,pst,
				  (void (*)(void *,void *,int,void *,int *)) 
				  pstSphViscosityLimiter, 
				  sizeof(struct inSphViscosityLimiter),0);
#ifndef NOCOOLING
	mdlAddService(mdl,PST_INITENERGY,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstInitEnergy,
				  sizeof(struct inInitEnergy),0);
	mdlAddService(mdl,PST_INITCOOLING,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstInitCooling,
				  sizeof(struct inInitCooling),0);
	mdlAddService(mdl,PST_COOLTABLEREAD,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstCoolTableRead,
				  CL_NMAXBYTETABLE,0);
#endif
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
	mdlAddService(mdl,PST_GETNPARTS,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstGetNParts,
				  0,sizeof(struct outGetNParts));
	mdlAddService(mdl,PST_SETNPARTS,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstSetNParts,
				  sizeof(struct inSetNParts),0);
#ifdef SPECIAL_PARTICLES
	mdlAddService(mdl,PST_GETSPECIALPARTICLES,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstGetSpecialParticles,
				  sizeof(struct inGetSpecial),sizeof(struct outGetSpecial));
	mdlAddService(mdl,PST_DOSPECIALPARTICLES,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstDoSpecialParticles,
				  sizeof(struct inDoSpecial),sizeof(struct outDoSpecial));
#endif
#ifdef COLLISIONS
	mdlAddService(mdl,PST_SETBALL,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstSetBall,
				  sizeof(struct inSetBall),0);
	mdlAddService(mdl,PST_NUMREJECTS,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstNumRejects,
				  0,sizeof(struct outNumRejects));
	mdlAddService(mdl,PST_READSS,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstReadSS,
				  sizeof(struct inReadSS),0);
	mdlAddService(mdl,PST_WRITESS,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstWriteSS,
				  sizeof(struct inWriteSS),0);
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
	mdlAddService(mdl,PST_FINDBINARY,pst,
		                  (void (*)(void *,void *,int,void *,int *)) pstFindTightestBinary,
		      sizeof(struct inBinary),sizeof(struct outBinary));
	mdlAddService(mdl,PST_MERGEBINARY,pst,
		                  (void (*)(void *,void *,int,void *,int *)) pstMergeBinary,
		      sizeof(struct inMrgBnry),sizeof(struct outMrgBnry));
#ifdef OLD_KEPLER
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
#endif /* OLD_KEPLER */

#ifdef SLIDING_PATCH
	mdlAddService(mdl,PST_FINDLM,pst,
		                  (void (*)(void *,void *,int,void *,int *)) pstFindLargeMasses,
				  sizeof(struct inLargeMass),sizeof(struct outLargeMass));
	mdlAddService(mdl,PST_GETNEIGHBORS,pst,
		                  (void (*)(void *,void *,int,void *,int *)) pstGetNeighborParticles,
				  sizeof(struct inGetNeighbors),sizeof(struct outGetNeighbors));
#endif /* SLIDING_PATCH */
#endif /* COLLISIONS */
#ifdef AGGS
	mdlAddService(mdl,PST_AGGSFIND,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstAggsFind,
				  0,sizeof(struct outAggsFind));
	mdlAddService(mdl,PST_AGGSCONFIRM,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstAggsConfirm,
				  sizeof(struct inAggsConfirm),sizeof(struct outAggsConfirm));
	mdlAddService(mdl,PST_AGGSMERGE,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstAggsMerge,
				  sizeof(struct inAggsMerge),0);
	mdlAddService(mdl,PST_AGGSBACKDRIFT,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstAggsBackDrift,
				  sizeof(struct inAggsBackDrift),0);
	mdlAddService(mdl,PST_AGGSGETCOM,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstAggsGetCOM,
				  sizeof(struct inAggsGetCOM),sizeof(struct outAggsGetCOM));
	mdlAddService(mdl,PST_AGGSGETAXESANDSPIN,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstAggsGetAxesAndSpin,
				  sizeof(struct inAggsGetAxesAndSpin),sizeof(struct outAggsGetAxesAndSpin));
	mdlAddService(mdl,PST_AGGSSETBODYPOS,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstAggsSetBodyPos,
				  sizeof(struct inAggsSetBodyPos),0);
	mdlAddService(mdl,PST_AGGSSETSPACEPOS,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstAggsSetSpacePos,
				  sizeof(struct inAggsSetSpacePos),0);
	mdlAddService(mdl,PST_AGGSSETSPACEVEL,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstAggsSetSpaceVel,
				  sizeof(struct inAggsSetSpaceVel),0);
	mdlAddService(mdl,PST_AGGSSETSPACESPINS,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstAggsSetSpaceSpins,
				  sizeof(struct inAggsSetSpaceSpins),0);
	mdlAddService(mdl,PST_AGGSDELETE,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstAggsDelete,
				  sizeof(struct inAggsDelete),sizeof(struct outAggsDelete));
	mdlAddService(mdl,PST_AGGSGETACCEL,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstAggsGetAccel,
				  sizeof(struct inAggsGetAccel),sizeof(struct outAggsGetAccel));
	mdlAddService(mdl,PST_AGGSCHECKSTRESS,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstAggsCheckStress,
				  sizeof(struct inAggsCheckStress),sizeof(struct outAggsCheckStress));
	mdlAddService(mdl,PST_AGGSGETTORQUE,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstAggsGetTorque,
				  sizeof(struct inAggsGetTorque),sizeof(struct outAggsGetTorque));
#endif /* AGGS */
#ifdef SLIDING_PATCH
	mdlAddService(mdl,PST_RANDAZWRAP,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstRandAzWrap,
				  sizeof(struct inRandAzWrap),sizeof(struct outRandAzWrap));
#endif /* SLIDING_PATCH */
#ifdef RUBBLE_ZML
	mdlAddService(mdl,PST_DUSTBINSGETMASS,pst,
				  (void (*)(void*,void*,int,void*,int*)) pstDustBinsGetMass,
				  sizeof(struct inDustBinsGetMass),sizeof(struct outDustBinsGetMass));
	mdlAddService(mdl,PST_DUSTBINSAPPLY,pst,
				  (void (*)(void*,void*,int,void*,int*)) pstDustBinsApply,
				  sizeof(struct inDustBinsApply),0);
	mdlAddService(mdl,PST_RUBBLERESETCOLFLAG,pst,
				  (void (*)(void*,void*,int,void*,int*)) pstRubbleResetColFlag,0,0);
	mdlAddService(mdl,PST_RUBBLECHECKFORKDKRESTART,pst,
				  (void (*)(void*,void*,int,void*,int*)) pstRubbleCheckForKDKRestart,
				  0,sizeof(struct outRubbleCheckForKDKRestart));
	mdlAddService(mdl,PST_RUBBLESTEP,pst,
				  (void (*)(void*,void*,int,void*,int*)) pstRubbleStep,
				  sizeof(struct inRubbleStep),0);
	mdlAddService(mdl,PST_RUBCLEANUP,pst,
				  (void (*)(void*,void*,int,void*,int*)) pstRubCleanup,
				  sizeof(struct inRubCleanup),sizeof(struct outRubCleanup));
	mdlAddService(mdl,PST_RUBINTERPCLEANUP,pst,
				  (void (*)(void*,void*,int,void*,int*)) pstRubInterpCleanup,
				  sizeof(struct inRubInterpCleanup),
				  sizeof(struct outRubInterpCleanup));
#endif /* RUBBLE_ZML */
	mdlAddService(mdl,PST_FORMSINKS,pst,
		      (void (*)(void *,void *,int,void *,int *)) pstFormSinks,
		      sizeof(struct inFormSinks),sizeof(struct outFormSinks));
#ifdef STARFORM
	mdlAddService(mdl,PST_FORMSTARS,pst,
		      (void (*)(void *,void *,int,void *,int *)) pstFormStars,
		      sizeof(struct inFormStars),sizeof(struct outFormStars));
	
	mdlAddService(mdl,PST_FEEDBACK,pst,
		      (void (*)(void *,void *,int,void *,int *))
		      pstFeedback, sizeof(struct inFeedback),
		      sizeof(struct outFeedback));
	mdlAddService(mdl,PST_INITSTARLOG,pst,
		      (void (*)(void *,void *,int,void *,int *))
		      pstInitStarLog, 0, 0);
	mdlAddService(mdl,PST_FLUSHSTARLOG,pst,
		      (void (*)(void *,void *,int,void *,int *))
		      pstFlushStarLog, sizeof(struct inFlushStarLog), 0);
#endif
#ifdef SIMPLESF
	mdlAddService(mdl,PST_SIMPLESTARFORM,pst,
		      (void (*)(void *,void *,int,void *,int *)) pstSimpleStarForm,
		      sizeof(struct inSimpleStarForm),sizeof(struct outSimpleStarForm));
#endif
	mdlAddService(mdl,PST_CLEARTIMER,pst,
		      (void (*)(void *,void *,int,void *,int *)) pstClearTimer,
		      sizeof(struct inClearTimer),0);
	mdlAddService(mdl,PST_MASSINR,pst,
		      (void (*)(void *,void *,int,void *,int *)) pstMassInR,
		      sizeof(struct inMassInR), sizeof(struct outMassInR));
	mdlAddService(mdl,PST_ROTBARINIT,pst,
		      (void (*)(void *,void *,int,void *,int *)) pstInitRotBar,
		      sizeof(struct inRotBar), 0);
#ifdef NEED_VPRED
	mdlAddService(mdl,PST_KICKVPRED,pst,
				  (void (*)(void *,void *,int,void *,int *)) pstKickVpred, 
				  sizeof(struct inKickVpred),sizeof(struct outKick));
#endif
	mdlAddService(mdl,PST_DUMPFRAME,pst,
		      (void (*)(void *,void *,int,void *,int *))
		      pstDumpFrame, sizeof(struct inDumpFrame),
		      DF_NBYTEDUMPFRAME );
	mdlAddService(mdl,PST_COM,pst,
		      (void (*)(void *,void *,int,void *,int *))
		      pstCOM, 0,
		      sizeof(double)*12 );
	mdlAddService(mdl,PST_COMBYTYPE,pst,
		      (void (*)(void *,void *,int,void *,int *))
		      pstCOMByType, sizeof(int),
		      sizeof(double)*4 );
	mdlAddService(mdl,PST_OLDESTSTAR,pst,
		      (void (*)(void *,void *,int,void *,int *))
		      pstOldestStar, 0,
		      sizeof(double)*4 );
	mdlAddService(mdl,PST_SETSINK,pst,
		      (void (*)(void *,void *,int,void *,int *))
		      pstSetSink, sizeof(struct inSetSink),
		      sizeof(struct outSetSink) );
#ifdef VOXEL
	mdlAddService(mdl,PST_DUMPVOXEL,pst,
		      (void (*)(void *,void *,int,void *,int *))
		      pstDumpVoxel, sizeof(struct inDumpVoxel), 0);
#endif
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
	pst->iSplitDim = -1;
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

	mdlassert(pst->mdl,nIn == sizeof(struct inSetAdd));
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

	mdlassert(pst->mdl,nIn == sizeof(struct inLevelize));
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

	mdlassert(pst->mdl,nIn == sizeof(struct inGetMap));
	if (pst->nLeaves > 1) {
		tmp = malloc(mdlThreads(pst->mdl)*sizeof(int));
		mdlassert(pst->mdl,tmp != NULL);
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
#ifdef COLLISIONS
	struct inReadSS *in = vin;
#else
	struct inReadTipsy *in = vin;
#endif
	int *pout = vout;
	int nFileStart,nFileEnd,nFileTotal,nFileSplit,nStore;
	int *ptmp;
	int i;

#ifdef COLLISIONS
	mdlassert(pst->mdl,nIn == sizeof(struct inReadSS));
#else
	mdlassert(pst->mdl,nIn == sizeof(struct inReadTipsy));
#endif
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
		mdlassert(pst->mdl,ptmp != NULL);
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

	mdlassert(pst->mdl,nIn == sizeof(struct inReadTipsy));
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
		pkdReadTipsy(plcl->pkd,achInFile,nFileStart,nFileTotal,in->bStandard,in->iReadIOrder,
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
	mdlassert(pst->mdl,r1 <= s2);
	mdlassert(pst->mdl,r2 <= s1);
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
/* #define PARANOID_CHECK */

void _pstRootSplit(PST pst,int iSplitDim,double dMass, int bDoRootFind,
		   int bDoSplitDimFind, int bSplitWork)
{
#ifdef STARFORM
	int NUM_SAFETY = 64; 	/* Larger margin for extra particles */
#else
	int NUM_SAFETY = 4; 	/* minimum margin space per processor when
				   filling up memory */
#endif
	int nSafeTot;		/* total slop space we have to play with */
	int margin;             /* margin of accuracy for filling up memory */
	int iSloppy = 2;	/* controls whether we are sloppy in filling
				   up memory: 2=>sloppy, 1=> not */
	int d,ittr,nOut;
	int nLow=-1,nHigh=-1,nLowerStore,nUpperStore;
	FLOAT fLow,fHigh;
	FLOAT fl,fu,fm=-1,fmm;
	struct outFreeStore outFree;
	struct inWeight inWt;
	struct inWeightWrap inWtWrap;
	struct outWeight outWtLow;
	struct outWeight outWtHigh;
	struct inColRejects inCol;
	OREJ *pLowerRej,*pUpperRej;
	int *pidSwap,iRet;
	int nLowTot,nHighTot;
	int nLast;		/* number of particles at the last split
					   iteration */
	int nDiff=0;	/* Difference between one iteration and the next */
	char ach[256];	/* Debug */

	mdlTimer t;

	struct outMassCheck outMass;
#ifdef PARANOID_CHECK
	int i,iLowSum,iHighSum;
#endif
	int pFlag;	/* 0 => we are splitting all particles at once. 
			   1 => we first split active, and then inactive. */
	int nTotalActive;
	int dBnd;

	mdlZeroTimer(pst->mdl,&t);
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
	 ** Particles are ordered into active-inactive order in pstDomainDecomp.
	 */

	mdlprintf(pst->mdl,"_pstRootSplit: id %d Level %d\n",pst->mdl->idSelf,pst->iLvl);
	mdlPrintTimer(pst->mdl,"TIME START _pstRootSplit ",&t);
	/* Debug */
	/*
	sprintf(ach,"id: %d _pstRootSplit\n", pst->mdl->idSelf );
	mdlDiag(pst->mdl,ach);
	*/

	if(bDoSplitDimFind || pst->iSplitDim == -1) {
	    pst->iSplitDim = iSplitDim;
	}

	d = dBnd = pst->iSplitDim;

#ifdef COLLISIONS
	if(d > 2)
	  dBnd -= 3;
#else
	mdlassert(pst->mdl,d < 3);
#endif

	/*
	 * If we only have a few active particles do nothin'
	 */

	fl = pst->bnd.fMin[dBnd];
	fu = pst->bnd.fMax[dBnd];
	fm = pst->fSplit; /* This is from the previous timestep */
	ittr = -1;
	
        if (bDoRootFind || fm<fl || fm>fu ) {
	     fmm = (fl + fu)/2;
	     /*
	      * First find total number of active particles.
	      */
	     inWt.iSplitDim = d;
	     inWt.fSplit = fmm;
	     inWt.ittr = 0;
	     inWt.iSplitSide = 1;
	     if(bSplitWork)
		 pFlag = 1; /* Sometimes parallel runs decomposed better
				    with this set to 0 */
	     else
		 pFlag = 0;
	     inWt.pFlag = pFlag;

	     mdlReqService(pst->mdl,pst->idUpper,PST_WEIGHT,&inWt,sizeof(inWt));
	     inWt.iSplitSide = 0;
	     pstWeight(pst->pstLower,&inWt,sizeof(inWt),&outWtLow,NULL);
	     mdlGetReply(pst->mdl,pst->idUpper,&outWtHigh,NULL);
	     nTotalActive = outWtLow.nLow + outWtHigh.nLow
	       + outWtLow.nHigh + outWtHigh.nHigh;
	     if (!bDoRootFind) {
	          pFlag = 0; /* Divide them all */
		  mdlReqService(pst->mdl,pst->idUpper,PST_WEIGHT,&inWt,sizeof(inWt));
		  inWt.iSplitSide = 0;
		  pstWeight(pst->pstLower,&inWt,sizeof(inWt),&outWtLow,NULL);
		  mdlGetReply(pst->mdl,pst->idUpper,&outWtHigh,NULL);
		  nTotalActive = outWtLow.nLow + outWtHigh.nLow
		    + outWtLow.nHigh + outWtHigh.nHigh;
	          }

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
	     
	     mdlPrintTimer(pst->mdl,"TIME active split _pstRootSplit ",&t);
	     }

	pst->fSplit = fm;

	mdlprintf(pst->mdl, "id: %d (%d) Chose split: %f (%f,%f) %d %d\n",
		  pst->mdl->idSelf, pst->iLvl, fm, pst->bnd.fMin[dBnd],
		  pst->bnd.fMax[dBnd], pst->nLower, pst->nUpper);
	if(ittr != -1)
	    mdlprintf(pst->mdl, "  Low %d %f,  High %d %f\n",
		      nLow,outWtLow.fLow + outWtHigh.fLow, nHigh,
		      outWtLow.fHigh + outWtHigh.fHigh);
	nLow = 0;
	nHigh = 0;
	fLow = 0.0;
	fHigh = 0.0;

	{
	    /*
	     ** Now we see if the TOTAL number of particles in the lower and upper
	     ** subsets exceeds the local pStores. If so then we need to find a new
	     ** boundary to distribute the INACTIVE particles so that everything 
	     ** fits.  
	     */
	    inWtWrap.iSplitDim = d;
	    fl = pst->fSplit + 1e-6*(pst->bnd.fMax[dBnd]-pst->bnd.fMin[dBnd]);
	    fu = pst->fSplit - 1e-6*(pst->bnd.fMax[dBnd]-pst->bnd.fMin[dBnd]);
	    
	    if (!bDoSplitDimFind) fm = pst->fSplitInactive;
	    else {
	      if (fu > fl) fm = 0.5*(fl+fu);
	      else {
		fm = 0.5*(fl+fu+pst->bnd.fMax[dBnd]-pst->bnd.fMin[dBnd]);
		if (fm > pst->bnd.fMax[dBnd]) fm = 0.5*(fl+fu-(pst->bnd.fMax[dBnd]-pst->bnd.fMin[dBnd]));
	      }
	      if (fabs(fm-pst->bnd.fMin[dBnd]) < fabs(fm-pst->bnd.fMax[dBnd])) {
		fm = pst->bnd.fMin[dBnd]- 1e-6*(pst->bnd.fMax[dBnd]-pst->bnd.fMin[dBnd]);
	      }
	      else {
		fm = pst->bnd.fMax[dBnd]+ 1e-6*(pst->bnd.fMax[dBnd]-pst->bnd.fMin[dBnd]);
	      }
	    }
	    mdlprintf(pst->mdl, "id: %d (%d) Zeroeth guess reverse split: %f (%f,%f)\n",
		      pst->mdl->idSelf, pst->iLvl, fm, pst->bnd.fMin[dBnd], pst->bnd.fMax[dBnd]);
            inWtWrap.fSplit = fm;	
	    inWtWrap.fSplit2 = pst->fSplit;
	    inWtWrap.ittr = 0;
	    inWtWrap.iSplitSide = 1;
	    inWtWrap.pFlag = 0;
	    mdlReqService(pst->mdl,pst->idUpper,PST_WEIGHTWRAP,&inWtWrap,sizeof(inWtWrap));
	    inWtWrap.iSplitSide = 0;
	    pstWeightWrap(pst->pstLower,&inWtWrap,sizeof(inWtWrap),&outWtLow,NULL);
	    mdlGetReply(pst->mdl,pst->idUpper,&outWtHigh,NULL);
	    /*
	    ** Add lower and Upper subsets numbers of particles
	    */ 
	    nLowTot = outWtLow.nLow + outWtHigh.nLow;
	    nHighTot = outWtLow.nHigh + outWtHigh.nHigh;
	    
	    nSafeTot = nLowerStore + nUpperStore - (nLowTot + nHighTot);
	    if(nSafeTot <= iSloppy*NUM_SAFETY*pst->nLeaves)
		iSloppy = 1; /* NO slop to play with */

	    if(nSafeTot/pst->nLeaves < NUM_SAFETY) {
	           NUM_SAFETY = nSafeTot/pst->nLeaves;
		   mdlprintf(pst->mdl,"id: %d tripped inactive NUM_SAFETY %d  Low %d/%d  High %d/%d\n",
			     pst->mdl->idSelf, NUM_SAFETY, nLowTot,
			     nLowerStore, nHighTot, nUpperStore);
	    }
	    
	    /*
	     ** Only be accurate to 5% of available space
	     */
	    margin = 0.05*nSafeTot/pst->nLeaves;
	    if (margin < NUM_SAFETY) margin = NUM_SAFETY;

	    mdlprintf(pst->mdl,"id: %d  %d Low %d/%d   %d High %d/%d  NUM_SAFETY %d margin %d\n",
		      pst->mdl->idSelf, pst->nLower,nLowTot, nLowerStore, pst->nUpper,nHighTot, nUpperStore,NUM_SAFETY,margin);
	    
	    
	    if (nLowTot > nLowerStore-NUM_SAFETY*pst->nLower) {
	            sprintf(ach,"id: %d: nLowTot > nLowerStore-NUM_SAFETY*pst->nLower %d %d %d %d\n",
			    pst->mdl->idSelf, nLowTot, nLowerStore, NUM_SAFETY, pst->nLower);
		    mdlDiag(pst->mdl,ach);
		    if (fm > pst->bnd.fMax[dBnd]) fm=pst->bnd.fMax[dBnd];
		    if (fm < pst->bnd.fMin[dBnd]) fm=pst->bnd.fMin[dBnd];
		    fl = fm;
		    if (fu > fl) fmm = 0.5*(fl+fu);
		    else {
		      fmm = 0.5*(fl+fu+pst->bnd.fMax[dBnd]-pst->bnd.fMin[dBnd]);
		      if (fmm > pst->bnd.fMax[dBnd]) fmm = 0.5*(fl+fu-(pst->bnd.fMax[dBnd]-pst->bnd.fMin[dBnd]));
		      mdlassert(pst->mdl, fmm >= pst->bnd.fMin[dBnd] && fmm <= pst->bnd.fMax[dBnd]);
		    }
		    ittr = 1;
		    nLast = nLowTot;
		    while (ittr < MAX_ITTR) {
		          fm = fmm;
			  inWtWrap.iSplitDim = d;
			  inWtWrap.fSplit = fm;	
			  inWtWrap.fSplit2 = pst->fSplit;
			  inWtWrap.ittr = ittr;
			  inWtWrap.iSplitSide = 1;
			  inWtWrap.pFlag = 0;
			  mdlReqService(pst->mdl,pst->idUpper,PST_WEIGHTWRAP,&inWtWrap,sizeof(inWtWrap));
			  inWtWrap.iSplitSide = 0;
			  pstWeightWrap(pst->pstLower,&inWtWrap,sizeof(inWtWrap),&outWtLow,NULL);
			  mdlGetReply(pst->mdl,pst->idUpper,&outWtHigh,NULL);
				/*
				** Add lower and Upper subsets numbers of particles
				*/
			  nLowTot = outWtLow.nLow + outWtHigh.nLow;
			  nHighTot = outWtLow.nHigh + outWtHigh.nHigh;
				/*
				  mdlprintf(pst->mdl, "id: %d (%d) %d th guess reverse split: %f (%f,%f) (%f,%f) Low %d High %d\n",
				  pst->mdl->idSelf, pst->iLvl, ittr, fm, fl, fu, pst->bnd.fMin[dBnd], pst->bnd.fMax[dBnd],nLowTot,nHighTot);
				*/
			  if(nLowTot != nLast)
			    nDiff = nLowTot - nLast;
			  nLast = nLowTot;
			  if (nLowTot > nLowerStore-margin*pst->nLower) fl = fm;
			  else if (nLowTot < nLowerStore-iSloppy*margin*pst->nLower) fu = fm;
			  else {
			    fl = fm;
			    break;
			  }
			  if (fu == fl)
			      break;
			  else if (fu > fl) fmm = 0.5*(fl+fu);
			  else {
			    fmm = 0.5*(fl+fu+(pst->bnd.fMax[dBnd]-pst->bnd.fMin[dBnd]));
			    if (fmm > pst->bnd.fMax[dBnd])
				fmm = 0.5*(fl+fu-(pst->bnd.fMax[dBnd]-pst->bnd.fMin[dBnd]));
			    /*
		 	     * Make sure split is inside bounds (within
			     * roundoff)
			     */
			    mdlassert(pst->mdl,
				      fmm - pst->bnd.fMin[dBnd] >= -1e-14*fabs(fmm)
				      && fmm <= pst->bnd.fMax[dBnd]);
			  }
			  ++ittr;
		    }
		    mdlprintf(pst->mdl, "id: %d (%d) Fix Low %d th guess reverse split: %f (%f,%f) (%f,%f) Low %d High %d\n",
			      pst->mdl->idSelf, pst->iLvl, ittr, fm, fl, fu, pst->bnd.fMin[dBnd], pst->bnd.fMax[dBnd],nLowTot,nHighTot);
		    if(nLowTot != nLowerStore-NUM_SAFETY*pst->nLower) {
		      if(abs(nDiff) > 1)
			mdlprintf(pst->mdl, "id: %d delta of %d, check NUM_SAFTEY\n",
				  pst->mdl->idSelf, nDiff);
		    }
		    mdlassert(pst->mdl,nLowTot <= nLowerStore);
		    mdlPrintTimer(pst->mdl,"TIME fix lower II _pstRootSplit ",&t);
	    }
	    else if (nHighTot > nUpperStore-NUM_SAFETY*pst->nUpper) {
	            sprintf(ach,"id: %d: nHighTot > nUpperStore-NUM_SAFETY*pst->nUpper %d %d %d %d\n",
		      pst->mdl->idSelf, nHighTot, nUpperStore, NUM_SAFETY, pst->nUpper);
		    mdlDiag(pst->mdl,ach);
		    if (fm > pst->bnd.fMax[dBnd]) fm=pst->bnd.fMax[dBnd];
		    if (fm < pst->bnd.fMin[dBnd]) fm=pst->bnd.fMin[dBnd];
		    fu = fm;
		    if (fu > fl) fmm = 0.5*(fl+fu);
		    else {
		      fmm = 0.5*(fl+fu+pst->bnd.fMax[dBnd]-pst->bnd.fMin[dBnd]);
		      if (fmm > pst->bnd.fMax[dBnd]) fmm = 0.5*(fl+fu-(pst->bnd.fMax[dBnd]-pst->bnd.fMin[dBnd]));
		      mdlassert(pst->mdl, fmm >= pst->bnd.fMin[dBnd] && fmm <= pst->bnd.fMax[dBnd]);
		    }
		    ittr = 1; 
		    nLast = nLowTot;
		    while (ittr < MAX_ITTR) {
		          fm = fmm;
			  inWtWrap.iSplitDim = d;
			  inWtWrap.fSplit = fm;
			  inWtWrap.fSplit2 = pst->fSplit;
			  inWtWrap.ittr = ittr;
			  inWtWrap.iSplitSide = 1;
			  inWtWrap.pFlag = 0;
			  mdlReqService(pst->mdl,pst->idUpper,PST_WEIGHTWRAP,&inWtWrap,sizeof(inWtWrap));
			  inWtWrap.iSplitSide = 0;
			  pstWeightWrap(pst->pstLower,&inWtWrap,sizeof(inWtWrap),&outWtLow,NULL);
			  mdlGetReply(pst->mdl,pst->idUpper,&outWtHigh,NULL);
			  /*
			  ** Add lower and Upper subsets numbers of particles
			  */
			  nLowTot = outWtLow.nLow + outWtHigh.nLow;
			  nHighTot = outWtLow.nHigh + outWtHigh.nHigh;
				/*				
								mdlprintf(pst->mdl, "id: %d (%d) %d th guess reverse split: %f (%f,%f) (%f,%f) Low %d High %d\n",
								pst->mdl->idSelf, pst->iLvl, ittr, fm, fl, fu, pst->bnd.fMin[dBnd], pst->bnd.fMax[dBnd],nLowTot,nHighTot);
				*/
			  if(nLowTot != nLast)
			    nDiff = nLowTot - nLast;
			  nLast = nLowTot;
			  if (nHighTot > nUpperStore-margin*pst->nUpper) fu = fm;
			  else if (nHighTot < nUpperStore-iSloppy*margin*pst->nUpper) fl = fm;
			  else {
			    fu = fm;
			    break;
			  }
			  if (fu == fl)
			      break;
			  else if (fu > fl) fmm = 0.5*(fl+fu);
			  else {
			    fmm = 0.5*(fl+fu+(pst->bnd.fMax[dBnd]-pst->bnd.fMin[dBnd]));
			    if (fmm > pst->bnd.fMax[dBnd])
				fmm = 0.5*(fl+fu-(pst->bnd.fMax[dBnd]-pst->bnd.fMin[dBnd]));
			    /*
		 	     * Make sure split is inside bounds (within
			     * roundoff)
			     */
			    mdlassert(pst->mdl,
				      fmm - pst->bnd.fMin[dBnd] >= -1e-14*fabs(fmm)
				      && fmm <= pst->bnd.fMax[dBnd]);
			  }
			  ++ittr;
		    }
		    mdlprintf(pst->mdl, "id: %d (%d) Fix High %d th guess reverse split: %f (%f,%f) (%f,%f) Low %d High %d\n",
			      pst->mdl->idSelf, pst->iLvl, ittr, fm, fl, fu, pst->bnd.fMin[dBnd], pst->bnd.fMax[dBnd],nLowTot,nHighTot);
		    if(nHighTot != nUpperStore-NUM_SAFETY*pst->nUpper) {
		      if(abs(nDiff) > 1)
			mdlprintf(pst->mdl, "id: %d delta of %d, check NUM_SAFETY\n",
				  pst->mdl->idSelf, nDiff);
		    }
		    mdlassert(pst->mdl,nHighTot <= nUpperStore);
		    mdlPrintTimer(pst->mdl,"TIME fix upper II _pstRootSplit ",&t);
	    }

	    mdlassert(pst->mdl, nLowTot >= 0);
	    mdlassert(pst->mdl, nHighTot >= 0);
	    mdlassert(pst->mdl, nLowTot <= nLowerStore);
	    mdlassert(pst->mdl, nHighTot <= nUpperStore);

	}
	    
	mdlPrintTimer(pst->mdl,"TIME Total Split _pstRootSplit ",&t);
	
	mdlprintf(pst->mdl, "id: %d (%d) Chose reverse split: %f (%f,%f)\n",
		pst->mdl->idSelf, pst->iLvl, fm, pst->bnd.fMin[dBnd], pst->bnd.fMax[dBnd]);
	pst->fSplitInactive = fm;
	
	/*
	 ** First Collect rejects.
	 **
	 ** Careful, SERVICE PST_COLREJECTS does NOT conform strictly to
	 ** the proper use of MDL. This should be fixed in the future.
	 ** FIXED -- MDL modified.
	 */
	pLowerRej = malloc(pst->nLower*sizeof(OREJ));
	mdlassert(pst->mdl,pLowerRej != NULL);
	pUpperRej = malloc(pst->nUpper*sizeof(OREJ));
	mdlassert(pst->mdl,pUpperRej != NULL);
	pidSwap = malloc(mdlThreads(pst->mdl)*sizeof(int));
	mdlassert(pst->mdl,pidSwap != NULL);

	inCol.fSplit = pst->fSplit;
	inCol.fSplitInactive = pst->fSplitInactive;
	inCol.iSplitDim = pst->iSplitDim;
	inCol.iSplitSide = 1;
	mdlReqService(pst->mdl,pst->idUpper,PST_COLREJECTS,&inCol,sizeof(inCol));
	inCol.iSplitSide = 0;
	pstColRejects(pst->pstLower,&inCol,sizeof(inCol),pLowerRej,&nOut);
	mdlassert(pst->mdl,nOut/sizeof(OREJ) == pst->nLower);
	mdlGetReply(pst->mdl,pst->idUpper,pUpperRej,&nOut);
	mdlassert(pst->mdl,nOut/sizeof(OREJ) == pst->nUpper);

	mdlPrintTimer(pst->mdl,"TIME Collected Rejects _pstRootSplit ",&t);

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
	mdlprintf(pst->mdl,"%d l:%d Paranoid Check Final nLow:%d == %d, nHigh:%d == %d\n",pst->idSelf,pst->iLvl,
		   nLowTot,iLowSum,nHighTot,iHighSum);
	iLowSum = 0;
	iHighSum = 0;
	for (i=0;i<pst->nLower;++i) {
		iLowSum += pLowerRej[i].nLocal;
		iHighSum += pLowerRej[i].nRejects;
		}
	mdlprintf(pst->mdl,"%d l:%d Paranoid Check Low: nLow:%d == %d, nHigh(rejects):%d == %d\n",pst->idSelf,pst->iLvl,
		   outWtLow.nLow,iLowSum,outWtLow.nHigh,iHighSum);
	iLowSum = 0;
	iHighSum = 0;
	for (i=0;i<pst->nUpper;++i) {
		iHighSum += pUpperRej[i].nLocal;
		iLowSum += pUpperRej[i].nRejects;
		}
	mdlprintf(pst->mdl,"%d l:%d Paranoid Check High: nLow(rejects):%d == %d, nHigh:%d == %d\n",pst->idSelf,pst->iLvl,
		   outWtHigh.nLow,iLowSum,outWtHigh.nHigh,iHighSum);
#endif
	
	pstMassCheck(pst,NULL,0,&outMass,NULL);
#if 0
	if (dMass != outMass.dMass) 
#else
	if (fabs(dMass - outMass.dMass) > MASS_EPS*dMass) 
#endif
	        {
		mdlprintf(pst->mdl,"ERROR id:%d lvl:%d:in _pstRootSplit after ColRejects\n",
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
		mdlassert(pst->mdl,nOut/sizeof(OREJ) == pst->nLower);
		mdlGetReply(pst->mdl,pst->idUpper,pUpperRej,&nOut);
		mdlassert(pst->mdl,nOut/sizeof(OREJ) == pst->nUpper);

		pstMassCheck(pst,NULL,0,&outMass,NULL);
#if 0
		if (dMass != outMass.dMass) 
#else
		if (fabs(dMass - outMass.dMass) > MASS_EPS*dMass) 
#endif
		        {
			printf("ERROR id:%d lvl:%d iter:%d in _pstRootSplit after Swap\n",
				   pst->idSelf,pst->iLvl,ittr);
		        }
		++ittr;
	        }

	free(pLowerRej);
	free(pUpperRej);
	free(pidSwap);

	mdlPrintTimer(pst->mdl,"TIME (FINISH) Swapped Rejects _pstRootSplit ",&t);
	
	}


void _pstRootSplit_Active_Inactive(PST pst,int iSplitDim,double dMass, int bDoRootFind)
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

	if(bDoRootFind) {
	    pst->iSplitDim = iSplitDim;
	}

	d = dBnd = pst->iSplitDim;

#ifdef COLLISIONS
	if(d > 2)
	  dBnd -= 3;
#else
	mdlassert(pst->mdl,d < 3);
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
		mdlprintf(pst->mdl, "id: %d tripped active NUM_SAFETY %d\n",
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
		mdlassert(pst->mdl,nLow <= nLowerStore);
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
		mdlassert(pst->mdl,nHigh <= nUpperStore);
		}
	pst->fSplit = fm;

	if(pFlag) {
	    if (!bDoRootFind) fm = pst->fSplitInactive;
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

	    nSafeTot = nLowerStore + nUpperStore - (nLowTot + nHighTot);
	    if(nSafeTot/pst->nLeaves < NUM_SAFETY) {
		NUM_SAFETY = nSafeTot/pst->nLeaves;
		sprintf(ach,"id: %d tripped inactive NUM_SAFETY %d\n",
				pst->mdl->idSelf, NUM_SAFETY);
		mdlDiag(pst->mdl,ach);
		mdlprintf(pst->mdl, "id: %d tripped inactive NUM_SAFETY %d\n",
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
					mdlprintf(pst->mdl, "id: %d delta of %d, check NUM_SAFTEY\n",
							pst->mdl->idSelf, nDiff);
				}
		    mdlassert(pst->mdl,nLowTot <= nLowerStore);
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
					mdlprintf(pst->mdl, "id: %d delta of %d, check NUM_SAFTEY\n",
							pst->mdl->idSelf, nDiff);
				}
		    mdlassert(pst->mdl,nHighTot <= nUpperStore);
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
	     ** Make sure everything is OK.
	     */
	    mdlassert(pst->mdl,nLowTot >= pst->nLower);
	    mdlassert(pst->mdl,nHighTot >= pst->nUpper);
	    mdlassert(pst->mdl,nLowTot <= nLowerStore);
	    mdlassert(pst->mdl,nHighTot <= nUpperStore);
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
	    mdlassert(pst->mdl,nLow >= pst->nLower);
	    mdlassert(pst->mdl,nHigh >= pst->nUpper);
	    mdlassert(pst->mdl,nLow <= nLowerStore);
	    mdlassert(pst->mdl,nHigh <= nUpperStore);
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
	mdlassert(pst->mdl,pLowerRej != NULL);
	pUpperRej = malloc(pst->nUpper*sizeof(OREJ));
	mdlassert(pst->mdl,pUpperRej != NULL);
	pidSwap = malloc(mdlThreads(pst->mdl)*sizeof(int));
	mdlassert(pst->mdl,pidSwap != NULL);

	inCol.fSplit = pst->fSplit;
	inCol.fSplitInactive = pst->fSplitInactive;
	inCol.iSplitDim = pst->iSplitDim;
	inCol.iSplitSide = 1;
	mdlReqService(pst->mdl,pst->idUpper,PST_COLREJECTS,&inCol,sizeof(inCol));
	inCol.iSplitSide = 0;
	pstColRejects(pst->pstLower,&inCol,sizeof(inCol),pLowerRej,&nOut);
	mdlassert(pst->mdl,nOut/sizeof(OREJ) == pst->nLower);
	mdlGetReply(pst->mdl,pst->idUpper,pUpperRej,&nOut);
	mdlassert(pst->mdl,nOut/sizeof(OREJ) == pst->nUpper);

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
		mdlassert(pst->mdl,nOut/sizeof(OREJ) == pst->nLower);
		mdlGetReply(pst->mdl,pst->idUpper,pUpperRej,&nOut);
		mdlassert(pst->mdl,nOut/sizeof(OREJ) == pst->nUpper);

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

#define NEWSPLITDIMCUT 0.707
#define NMINFORROOTFIND 16

void pstDomainDecomp(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	int nOut,d=0,j;
	double dimsize;
	struct outCalcBound outBnd;
	struct outMassCheck outMass;
	double dMass;
	struct inDomainDecomp *in = vin;
	int nActive;
	
	mdlTimer t;

	mdlassert(pst->mdl,nIn == sizeof(struct inDomainDecomp));

	mdlZeroTimer(pst->mdl,&t);
	mdlprintf(pst->mdl,"Starting pstDomainDecomp\n");

	if (pst->nLeaves > 1) {
		/*
		 ** First calculate the Bounds for the set.
		 */
	        pstActiveOrder(pst, NULL, 0, &nActive, NULL);
		if (pst->iSplitDim != -1 && nActive < NMINFORROOTFIND) {
		  mdlprintf(pst->mdl,"Aborting RootFind -- too few actives.\n");
		  in->bDoSplitDimFind = 0;
		  in->bDoRootFind = 0;
		}

		mdlPrintTimer(pst->mdl,"TIME active order done in pstDomainDecomp",&t);

		pstCalcBound(pst,NULL,0,&outBnd,&nOut);
		pst->bnd = outBnd.bnd;
		pst->bndActive = outBnd.bndActive;
		pst->bndTreeActive = outBnd.bndTreeActive;
		if(pst->bnd.fMin[0] > pst->bnd.fMax[0]) {
		/*
		 ** No particles in this branch!
		 */
		    pstMassCheck(pst,NULL,0,&outMass,NULL); /* sanity check */
		    assert(outMass.dMass == 0.0);
		    if (pnOut) *pnOut = 0;
		    return;
		    }
		/*
		 ** Next determine the longest axis 
		 */
		if (in->bDoSplitDimFind) {
			d = pst->iSplitDim;
			if (d==-1) dimsize = -1;
			else dimsize = (pst->bnd.fMax[d]-pst->bnd.fMin[d])*NEWSPLITDIMCUT;
			for (j=0;j<3;++j) {
				if ((pst->bnd.fMax[j]-pst->bnd.fMin[j]) > dimsize) {
					d=j;
					dimsize = (pst->bnd.fMax[d]-pst->bnd.fMin[d]);
					}
				}
			}

		mdlPrintTimer(pst->mdl,"TIME CalcBound done in pstDomainDecomp",&t);
		pstMassCheck(pst,NULL,0,&outMass,NULL);
		dMass = outMass.dMass;
		mdlPrintTimer(pst->mdl,"TIME Mass Check done in pstDomainDecomp",&t);
		_pstRootSplit(pst,d,dMass, in->bDoRootFind,
			      in->bDoSplitDimFind, in->bSplitWork );
		mdlPrintTimer(pst->mdl,"TIME RootSplit done in pstDomainDecomp",&t);
		pstMassCheck(pst,NULL,0,&outMass,NULL);
  		if (fabs(dMass - outMass.dMass) > MASS_EPS*dMass) {
			printf("ERROR id:%d lvl:%d:_pstRootSplit mass not cons.\n",
				   pst->idSelf,pst->iLvl);
			}
		mdlPrintTimer(pst->mdl,"TIME Mass Check done in pstDomainDecomp",&t);
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
	mdlPrintTimer(pst->mdl,"TIME Sublevels done in pstDomainDecomp ",&t);

	if (pnOut) *pnOut = 0;

	mdlPrintTimer(pst->mdl,"TIME DONE pstDomainDecomp ",&t);
	}


void pstCalcBound(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct outCalcBound *out = vout;
	struct outCalcBound outBnd;
	int j;
	
	mdlassert(pst->mdl,nIn == 0);
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
	
	mdlassert(pst->mdl,nIn == 0);
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


void pstRungDDWeight(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inRungDDWeight *in = vin;
	
	mdlassert(pst->mdl,nIn == 0);
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_RUNGDDWEIGHT,vin,nIn);
		pstRungDDWeight(pst->pstLower,vin,nIn,NULL,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
		}
	else {
		pkdRungDDWeight(plcl->pkd,in->iMaxRung,in->dWeight);
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

	mdlassert(pst->mdl,nIn == sizeof(struct inWeight));
	/*
	pkdStartTimer(plcl->pkd,7);
	*/
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
	/*
	pkdStopTimer(plcl->pkd,7);
	*/
	}


/*
 ** Make sure that the local particles are split into active and inactive
 ** when passing pFlag != 0.
 ** pFlag == 0 => weight all particles.
 ** pFlag > 0 => weight active particles.
 ** pFlag < 0 => weight inactive particles.
 */
void pstWeightWrap(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inWeightWrap *in = vin;
	struct outWeight *out = vout;
	struct outWeight outWt;
	FLOAT fSplit,fLow,fHigh;
	int iSplitSide;

	mdlassert(pst->mdl,nIn == sizeof(struct inWeightWrap));
	/*
	pkdStartTimer(plcl->pkd,7);
	*/
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_WEIGHTWRAP,in,nIn);
		pstWeightWrap(pst->pstLower,in,nIn,out,NULL);
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
			if ((fSplit < plcl->fSplit && fSplit<in->fSplit2 && plcl->fSplit<in->fSplit2) ||
			    (fSplit < plcl->fSplit && fSplit>in->fSplit2 && plcl->fSplit>in->fSplit2) ||
			    (fSplit > plcl->fSplit && fSplit>in->fSplit2 && plcl->fSplit<in->fSplit2)) {
				plcl->fWtHigh += plcl->fHigh;
				if (!iSplitSide) plcl->iWtFrom = plcl->iPart;
				else plcl->iWtTo = plcl->iPart-1;
				}
			else {
				plcl->fWtLow += plcl->fLow;
				if (!iSplitSide) plcl->iWtTo = plcl->iPart-1;
				else plcl->iWtFrom = plcl->iPart;
				}
			plcl->fSplit = fSplit;
			}
		plcl->iPart = pkdWeightWrap(plcl->pkd,in->iSplitDim,fSplit,in->fSplit2,iSplitSide,
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
	/*
	pkdStopTimer(plcl->pkd,7);
	*/
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

	mdlassert(pst->mdl,nIn == sizeof(struct inOrdWeight));
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

	mdlassert(pst->mdl,nIn == 0);
	/*
	pkdStartTimer(plcl->pkd,4);
	*/
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
	/*
	pkdStopTimer(plcl->pkd,4);
	*/
	}


void pstColRejects(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inColRejects *in = vin;
	OREJ *pOutRej = vout;
	int nLower,nUpper,iUpper;
	
	mdlassert(pst->mdl,nIn == sizeof(struct inColRejects));
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
	
	mdlassert(pst->mdl,nIn == sizeof(struct inColOrdRejects));
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
	
	/*	pkdStartTimer(plcl->pkd,8);*/
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
	/*	pkdStopTimer(plcl->pkd,8); */
	}

/*
 * Routine to swap all particles.  Note that this does not walk the pst
 * but simply works with one other processor.
 */
void pstSwapAll(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl;
	int *pidSwap = vin;
	PST lpst;

	mdlassert(pst->mdl,nIn == sizeof(*pidSwap));
	lpst = pst;
	while(lpst->nLeaves > 1)
	    lpst = lpst->pstLower;

	plcl = lpst->plcl;
	pkdSwapAll(plcl->pkd, *pidSwap);
	}


void pstDomainColor(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;

	mdlassert(pst->mdl,nIn == 0);
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
	mdlassert(pst->mdl,nLow <= nLowerStore);
	mdlassert(pst->mdl,nHigh <= nUpperStore);
	pst->iOrdSplit = im;
	/*
	 ** Collect rejects.
	 */
	pLowerRej = malloc(pst->nLower*sizeof(OREJ));
	mdlassert(pst->mdl,pLowerRej != NULL);
	pUpperRej = malloc(pst->nUpper*sizeof(OREJ));
	mdlassert(pst->mdl,pUpperRej != NULL);
	pidSwap = malloc(mdlThreads(pst->mdl)*sizeof(int));
	mdlassert(pst->mdl,pidSwap != NULL);
	inCol.iOrdSplit = pst->iOrdSplit;
	inCol.iSplitSide = 1;
	mdlReqService(pst->mdl,pst->idUpper,PST_COLORDREJECTS,&inCol,
				  sizeof(inCol));
	inCol.iSplitSide = 0;
	pstColOrdRejects(pst->pstLower,&inCol,sizeof(inCol),pLowerRej,&nOut);
	mdlassert(pst->mdl,nOut/sizeof(OREJ) == pst->nLower);
	mdlGetReply(pst->mdl,pst->idUpper,pUpperRej,&nOut);
	mdlassert(pst->mdl,nOut/sizeof(OREJ) == pst->nUpper);
	while (1) {
		iRet = _pstRejMatch(pst,pst->nLower,pLowerRej,pst->nUpper,
							pUpperRej,pidSwap);
		if (!iRet) break;
		mdlReqService(pst->mdl,pst->idUpper,PST_SWAPREJECTS,pidSwap,
					  mdlThreads(pst->mdl)*sizeof(int));
		pstSwapRejects(pst->pstLower,pidSwap,mdlThreads(pst->mdl)*sizeof(int),
					   pLowerRej,&nOut);
		mdlassert(pst->mdl,nOut/sizeof(OREJ) == pst->nLower);
		mdlGetReply(pst->mdl,pst->idUpper,pUpperRej,&nOut);
		mdlassert(pst->mdl,nOut/sizeof(OREJ) == pst->nUpper);
		}
	free(pLowerRej);
	free(pUpperRej);
	free(pidSwap);
	}


void pstDomainOrder(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	struct inDomainOrder *in = vin;

	mdlassert(pst->mdl,nIn == sizeof(struct inDomainOrder));
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

	mdlassert(pst->mdl,nIn == 0);
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

	mdlassert(pst->mdl,nIn == sizeof(*in));
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

	mdlassert(pst->mdl,nIn == 0);
	/*
	pkdStartTimer(pst->plcl->pkd,5);
	*/
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
	/*
	pkdStopTimer(pst->plcl->pkd,5);
	*/
	}


void pstOutArray(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inOutput *in = vin;
	char achOutFile[PST_FILENAME_SIZE];

	mdlassert(pst->mdl,nIn == sizeof(struct inOutArray));
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_OUTARRAY,in,nIn);
		pstOutArray(pst->pstLower,in,nIn,NULL,NULL);
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
		pkdOutVector(plcl->pkd,achOutFile,plcl->nWriteStart, 1,in->iType, in->iBinaryOutput, in->N,in->bStandard);
		}
	if (pnOut) *pnOut = 0;
	}


void pstOutNCVector(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inOutput *in = vin;
        struct outNC *out = vout;
        struct outNC outUp;
        int i,iDim;
	char achOutFile[PST_FILENAME_SIZE];

	mdlassert(pst->mdl,nIn == sizeof(struct inOutput));
	if (pst->nLeaves > 1) {
		/*
		 ** Non-Recursive Text output.
		 */
		mdlReqService(pst->mdl,pst->idUpper,PST_OUTNCVECTOR,in,nIn);
		pstOutNCVector(pst->pstLower,in,nIn,out,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,&outUp,NULL);
                if (out != NULL){
                    for(i=0; i<3; i++){
                        for(iDim=0;iDim<3;iDim++){
                            out->min[i][iDim]=(out->min[i][iDim] < outUp.min[i][iDim]) ? out->min[i][iDim]: outUp.min[i][iDim];
                            out->max[i][iDim]=(out->max[i][iDim] > outUp.max[i][iDim]) ? out->max[i][iDim]: outUp.max[i][iDim];
                            }
                        }
                    }
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
		pkdOutNChilada(plcl->pkd,achOutFile,plcl->nGasWriteStart,
			       plcl->nDarkWriteStart, plcl->nStarWriteStart,
			       in->iType,out->min, out->max, in->duTFac,
			       in->dvFac);
		}
	if (pnOut) *pnOut = sizeof(struct outNC);
	}


void pstOutVector(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inOutput *in = vin;
	char achOutFile[PST_FILENAME_SIZE];

	mdlassert(pst->mdl,nIn == sizeof(struct inOutput));
	if (pst->nLeaves > 1) {
		/*
		 ** Non-Recursive Text output.
		 */
		mdlReqService(pst->mdl,pst->idUpper,PST_OUTVECTOR,in,nIn);
		pstOutVector(pst->pstLower,in,nIn,NULL,NULL);
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
		pkdOutVector(plcl->pkd,achOutFile,plcl->nWriteStart, in->iDim,in->iType,in->iBinaryOutput,in->N,in->bStandard);
		}
	if (pnOut) *pnOut = 0;
	}


void pstWriteTipsy(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inWriteTipsy *in = vin;
	char achOutFile[PST_FILENAME_SIZE];

	mdlassert(pst->mdl,nIn == sizeof(struct inWriteTipsy));
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

	mdlassert(pst->mdl,nIn == sizeof(struct inSetSoft));
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

#ifdef CHANGESOFT
void pstPhysicalSoft(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inPhysicalSoft *in = vin;

	mdlassert(pst->mdl,nIn == sizeof(struct inPhysicalSoft));
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_PHYSICALSOFT,in,nIn);
		pstPhysicalSoft(pst->pstLower,in,nIn,NULL,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
		}
	else {
		pkdPhysicalSoft(plcl->pkd,in->dSoftMax,in->dFac,in->bSoftMaxMul);
		}
	if (pnOut) *pnOut = 0;
	}


void pstPreVariableSoft(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inPreVariableSoft *in = vin;

	mdlassert(pst->mdl,nIn == sizeof(struct inPreVariableSoft));
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_PREVARIABLESOFT,in,nIn);
		pstPreVariableSoft(pst->pstLower,in,nIn,NULL,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
		}
	else {
		pkdPreVariableSoft(plcl->pkd,in->iVariableSoftType);
		}
	if (pnOut) *pnOut = 0;
	}

void pstPostVariableSoft(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inPostVariableSoft *in = vin;

	mdlassert(pst->mdl,nIn == sizeof(struct inPostVariableSoft));
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_POSTVARIABLESOFT,in,nIn);
		pstPostVariableSoft(pst->pstLower,in,nIn,NULL,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
		}
	else {
		pkdPostVariableSoft(plcl->pkd,in->dSoftMax,in->bSoftMaxMul,in->iVariableSoftType);
		}
	if (pnOut) *pnOut = 0;
	}
#endif

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
        
	mdlassert(pst->mdl,nIn == sizeof(struct inBuildTree));
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
	struct outSmooth *out = vout;
	struct outSmooth outUp;
	CASTAT cs;

	mdlassert(pst->mdl,nIn == sizeof(struct inSmooth));
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_SMOOTH,in,nIn);
		pstSmooth(pst->pstLower,in,nIn,out,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,&outUp,NULL);
		if (out != NULL) {
		    out->iSmoothFlags |= outUp.iSmoothFlags;
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
		}
	else {
		LCL *plcl = pst->plcl;
		SMX smx;

		(&in->smf)->pkd = pst->plcl->pkd;
		smInitialize(&smx,plcl->pkd,&in->smf,in->nSmooth,in->bPeriodic,
					 in->bSymmetric,in->iSmoothType,1,in->dfBall2OverSoft2);
		smSmooth(smx,&in->smf);
		smFinish(smx,&in->smf, &cs);
		
		if (out != NULL) {
		    out->iSmoothFlags |= in->smf.iSmoothFlags;
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
		}
	if (pnOut) *pnOut = sizeof(struct outSmooth);
	}


void pstMarkSmooth(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	struct inMarkSmooth *in = vin;
	CASTAT cs;

	mdlassert(pst->mdl,nIn == sizeof(struct inMarkSmooth));
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
					 in->bSymmetric,in->iSmoothType,0,0.0);
		smMarkSmooth(smx,&in->smf,in->iMarkType);
		smFinish(smx,&in->smf, &cs);
		}
	if (pnOut) *pnOut = 0;
	}


void pstReSmooth(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	struct inReSmooth *in = vin;
	struct outReSmooth *out = vout;
	struct outReSmooth outUp;
	CASTAT cs;

	mdlassert(pst->mdl,nIn == sizeof(struct inReSmooth));
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_RESMOOTH,in,nIn);
		pstReSmooth(pst->pstLower,in,nIn,out,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,&outUp,NULL);

		if (out != NULL) {
		    out->iSmoothFlags |= outUp.iSmoothFlags;
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
		}
	else {
		LCL *plcl = pst->plcl;
		SMX smx;

		(&in->smf)->pkd = pst->plcl->pkd;
		smInitialize(&smx,plcl->pkd,&in->smf,in->nSmooth,in->bPeriodic,
					 in->bSymmetric,in->iSmoothType,0,in->dfBall2OverSoft2);
		smReSmooth(smx,&in->smf);
		smFinish(smx,&in->smf, &cs);

		if (out != NULL) {
		    out->iSmoothFlags |= in->smf.iSmoothFlags;
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
		}
	if (pnOut) *pnOut = sizeof(struct outReSmooth);
	}


void pstSoughtParticleList(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	struct inSoughtParticleList *in = vin;
	struct inoutParticleList *out = vout;
	struct inoutParticleList outUp;
	int i,n;

	mdlassert(pst->mdl,nIn == sizeof(struct inSoughtParticleList));
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_SOUGHTPARTICLELIST,in,nIn);
		pstSoughtParticleList(pst->pstLower,in,nIn,out,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,&outUp,NULL);

		n = out->n+outUp.n;
		if (n > in->nMax) n = in->nMax;
		for (i=out->n;i<n;i++) {
		  out->p[i] = outUp.p[i-out->n];
		}
		out->n += outUp.n; /* If this exceeds the max it quietly returns but stops
							  copying data: Caller should test  */
	    }
	else {
		pkdSoughtParticleList(pst->plcl->pkd, in->iTypeSought, in->nMax, &out->n, &(out->p[0]));
		}

	if (pnOut) *pnOut = sizeof(struct inoutParticleList);
	}


void pstCoolUsingParticleList(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	struct inoutParticleList *list = vin;

	mdlassert(pst->mdl,nIn == sizeof(struct inoutParticleList));
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_COOLUSINGPARTICLELIST,list,nIn);
		pstCoolUsingParticleList(pst->pstLower,list,nIn,vout,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,vout,NULL);
	    }
	else {
		pkdCoolUsingParticleList(pst->plcl->pkd, list->n, &(list->p[0]));
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

	mdlassert(pst->mdl,nIn == sizeof(struct inGravity));
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
		plcl->pkd->dTime = in->dTime;
		plcl->pkd->PP = &in->PP;
#endif
		pkdGravAll(plcl->pkd,in->nReps,in->bPeriodic,in->iOrder,in->bEwald,
				   in->iEwOrder,in->dEwCut,in->dEwhCut,in->bDoSun,in->dSunSoft,out->aSun,
				   &out->nActive,&out->dPartSum,&out->dCellSum,&out->dSoftSum,
				   &cs,&out->dFlop);
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


void pstGravExternal(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inGravExternal *in = vin;
	struct outGravExternal *out = vout;
	struct outGravExternal outLcl;

	mdlassert(pst->mdl,nIn == sizeof(struct inGravExternal));
	if (pst->nLeaves > 1) {
	        int j;
	    
		mdlReqService(pst->mdl,pst->idUpper,PST_GRAVEXTERNAL,in,nIn);
		pstGravExternal(pst->pstLower,in,nIn,out,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,&outLcl,NULL);
		if (out != NULL) {
		  for(j = 0; j < 3; j++) {
		    out->dAcc[j] += outLcl.dAcc[j];
		    out->dTorque[j] += outLcl.dTorque[j];
		  }
		}
		}
	else {
	    int j;

		if (out != NULL) {
		  for(j = 0; j < 3; j++) {
			out->dAcc[j] = 0.0;
			out->dTorque[j] = 0.0;
		  }
		}
		if (in->bIndirect) {
			pkdSunIndirect(plcl->pkd,in->aSun,in->bDoSun,in->dSunMass,in->dSunSoft);
			}
		if (in->bLogHalo) {
			pkdLogHalo(plcl->pkd);
			}
		if (in->bHernquistSpheroid) {
			pkdHernquistSpheroid(plcl->pkd);
			}
		if (in->bNFWSpheroid) {
			pkdNFWSpheroid(plcl->pkd,in->dNFWm200,in->dNFWr200,in->dNFWconc,in->dNFWsoft);
			}
		if (in->bElliptical) {
		    pkdElliptical(plcl->pkd,in->bEllipticalDarkNFW);
		    }
		if (in->bHomogSpheroid) {
			pkdHomogSpheroid(plcl->pkd);
			}
		if (in->bBodyForce) {
			pkdBodyForce(plcl->pkd, in->dBodyForceConst);
			}
		if (in->bMiyamotoDisk) {
			pkdMiyamotoDisk(plcl->pkd);
			}
		if (in->bTimeVarying) {
			pkdTimeVarying(plcl->pkd,in->dTime);
			}
		if (in->bRotatingBar) {
		  assert(out != NULL);
			pkdRotatingBar(plcl->pkd, in->dRotBarAmp,
				       in->dRotBarPosAng,
				       in->dRotBarB5, in->aCom, out->dAcc,
				       out->dTorque);
			}
#ifdef ROT_FRAME
		if (in->bRotFrame) {
			pkdRotFrame(plcl->pkd,in->dOmega,in->dOmegaDot);
			}
#endif
#ifdef SLIDING_PATCH
		if (in->bPatch) {
		  /*
		  ** Could just pass these to pkdPatch directly, but since the
		  ** PKD struct already has space for them, we'll just use that.
		  */
		  plcl->pkd->dTime = in->dTime;
		  plcl->pkd->PP = &in->PP;
		  pkdPatch(plcl->pkd);
		}
#endif
#ifdef SIMPLE_GAS_DRAG
		if (in->bSimpleGasDrag) {
			pkdSimpleGasDrag(plcl->pkd,in->iFlowOpt,in->bEpstein,in->dGamma,
							 in->dTime);
			}
#endif
		}
	if (pnOut) *pnOut = sizeof(struct outGravExternal);
	}


void pstCalcEandL(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct outCalcEandL *out = vout;
	struct outCalcEandL outLcl;
	int k;

	mdlassert(pst->mdl,nIn == 0);
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_CALCEANDL,NULL,0);
		pstCalcEandL(pst->pstLower,NULL,0,out,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,&outLcl,NULL);
		out->T += outLcl.T;
		out->U += outLcl.U;
		out->Eth += outLcl.Eth;
		for (k=0;k<3;k++) out->L[k] += outLcl.L[k];
		}
	else {
		pkdCalcEandL(plcl->pkd,&out->T,&out->U,&out->Eth,out->L);
		}
	if (pnOut) *pnOut = sizeof(struct outCalcEandL);
	}


void pstCalcEandLExt(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inCalcEandLExt *in = vin;
	struct outCalcEandLExt *out = vout,outLcl;

	mdlassert(pst->mdl,nIn == sizeof(struct inCalcEandLExt));
	mdlassert(pst->mdl,in->bHeliocentric); /* only one option supported currently */
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_CALCEANDLEXT,in,nIn);
		pstCalcEandLExt(pst->pstLower,in,nIn,out,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,&outLcl,NULL);
		}
	else {
		pkdCalcEandLExt(plcl->pkd,&out->dMass,out->dSumMR,out->dSumMV,
						&out->dPot);
		}
	if (pnOut) *pnOut = sizeof(struct outCalcEandLExt);
	}


void pstDrift(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inDrift *in = vin;

	mdlassert(pst->mdl,nIn == sizeof(struct inDrift));
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_DRIFT,in,nIn);
		pstDrift(pst->pstLower,in,nIn,NULL,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
		}
	else {
#ifdef SLIDING_PATCH
		plcl->pkd->dTime = in->dTime;
		plcl->pkd->PP = &in->PP;
#endif
		pkdDrift(plcl->pkd,in->dDelta,in->fCenter,in->bPeriodic,in->bFandG,
				 in->fCentMass);
		}
	if (pnOut) *pnOut = 0;
	}


void pstKick(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inKick *in = vin;
	struct outKick *out = vout;
	struct outKick outUp;

	mdlassert(pst->mdl,nIn == sizeof(struct inKick));

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

void pstKickPatch(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inKickPatch *in = vin;
	struct outKick *out = vout;
	struct outKick outUp;

	mdlassert(pst->mdl,nIn == sizeof(struct inKickPatch));

	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_KICKPATCH,in,nIn);
		pstKickPatch(pst->pstLower,in,nIn,out,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,&outUp,NULL);

		out->SumTime += outUp.SumTime;
		out->nSum += outUp.nSum;
		if (outUp.MaxTime > out->MaxTime) out->MaxTime = outUp.MaxTime;
		}
	else {
		pkdKickPatch(plcl->pkd,in->dvFacOne,in->dvFacTwo,
			     in->dOrbFreq, in->bOpen);
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

	mdlassert(pst->mdl,nIn == sizeof(struct inReadCheck));
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

	mdlassert(pst->mdl,nIn == 0);
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


void pstSetTotals(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{	
	LCL *plcl = pst->plcl;
	struct outSetTotals *out = vout;
	struct outSetTotals oute;
        int nDark, nGas, nStar;

	mdlassert(pst->mdl,nIn == 0);
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_SETTOTALS,NULL,0);
		pstSetTotals(pst->pstLower,NULL,0,out,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,&oute,NULL);
		out->nDark += oute.nDark;
		pst->nDark = out->nDark;
		out->nGas += oute.nGas;
		pst->nGas = out->nGas;
		out->nStar += oute.nStar;
		pst->nStar = out->nStar;
		}
	else {
		pkdTotals(plcl->pkd,&nDark,&nGas,&nStar);
                pst->nDark = nDark;
		out->nDark = pst->nDark;
		pst->nGas = nGas;
		out->nGas = pst->nGas;
		pst->nStar = nStar;
		out->nStar = pst->nStar;
		}
	if (pnOut) *pnOut = sizeof(struct outSetTotals);
	}


void pstSetWriteStart(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{	
	LCL *plcl = pst->plcl;
	struct inSetWriteStart *in = vin;
	int nWriteStart;

	mdlassert(pst->mdl,nIn == sizeof(struct inSetWriteStart));
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


void pstSetNCWriteStart(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{	
	LCL *plcl = pst->plcl;
	struct inSetNCWriteStart *in = vin;
	int nDarkWriteStart,nGasWriteStart,nStarWriteStart;

	mdlassert(pst->mdl,nIn == sizeof(struct inSetNCWriteStart));
	nDarkWriteStart = in->nDarkWriteStart;
	nGasWriteStart = in->nGasWriteStart;
	nStarWriteStart = in->nStarWriteStart;
	if (pst->nLeaves > 1) {
		in->nDarkWriteStart = nDarkWriteStart + pst->pstLower->nDark;
		in->nGasWriteStart = nGasWriteStart + pst->pstLower->nGas;
		in->nStarWriteStart = nStarWriteStart + pst->pstLower->nStar;
		mdlReqService(pst->mdl,pst->idUpper,PST_SETNCWRITESTART,in,nIn);
		in->nDarkWriteStart = nDarkWriteStart;
		in->nGasWriteStart = nGasWriteStart;
		in->nStarWriteStart = nStarWriteStart;
		pstSetNCWriteStart(pst->pstLower,in,nIn,NULL,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
		}
	else {
		plcl->nDarkWriteStart = nDarkWriteStart;
		plcl->nGasWriteStart = nGasWriteStart;
		plcl->nStarWriteStart = nStarWriteStart;
		}
	if (pnOut) *pnOut = 0;
	}


void pstWriteCheck(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inWriteCheck *in = vin;
	char achOutFile[PST_FILENAME_SIZE];

	mdlassert(pst->mdl,nIn == sizeof(struct inWriteCheck));
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

	mdlassert(pst->mdl,nIn == sizeof(struct inCalcCell));
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

	mdlassert(pst->mdl,nIn == sizeof(struct inColCells));
	iCell = in->iCell;
	if (pst->nLeaves > 1) {
		in->iCell = UPPER(iCell);
		mdlReqService(pst->mdl,pst->idUpper,PST_COLCELLS,in,nIn);
		in->iCell = LOWER(iCell);
		pstColCells(pst->pstLower,in,nIn,vout,pnOut);
		in->iCell = iCell;
		ptmp = malloc(in->nCell*sizeof(KDN));
		mdlassert(pst->mdl,ptmp != NULL);
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
		mdlassert(pst->mdl,pkdn[iCell].pUpper != 0);
		}
	else {
		for (i=1;i<in->nCell;++i) pkdn[i].pUpper = 0; /* used flag = unused */
		pst->kdn.iLower = -1;
		pst->kdn.iUpper = -1;
		pkdn[iCell] = pst->kdn;
		mdlassert(pst->mdl,pkdn[iCell].pUpper != 0);
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

	mdlassert(pst->mdl,nIn == 0);
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

	mdlassert(pst->mdl,nIn == sizeof(struct ioCalcRoot));
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
	
	mdlassert(pst->mdl,nIn == 0);
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
        
void pstMassMetalsEnergyCheck(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct outMassMetalsEnergyCheck *out = vout;
	struct outMassMetalsEnergyCheck outUp;
	
	mdlassert(pst->mdl,nIn == 0);
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_MASSMETALSENERGYCHECK,NULL,0);
		pstMassMetalsEnergyCheck(pst->pstLower,NULL,0,vout,pnOut);
		mdlGetReply(pst->mdl,pst->idUpper,&outUp,pnOut);
		out->dTotMass += outUp.dTotMass;
		out->dTotMetals += outUp.dTotMetals;
		out->dTotFe += outUp.dTotFe;
		out->dTotOx += outUp.dTotOx;
		out->dTotEnergy += outUp.dTotEnergy;
		}
	else {
		pkdMassMetalsEnergyCheck(plcl->pkd,&out->dTotMass,&out->dTotMetals,
                            &out->dTotOx,&out->dTotFe,&out->dTotEnergy);
		}
	if (pnOut) *pnOut = sizeof(struct outMassMetalsEnergyCheck);
	}

void
pstSetRung(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inSetRung *in = vin;
	
	mdlassert(pst->mdl,nIn == sizeof(*in));
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
	
	mdlassert(pst->mdl,nIn == sizeof(*in));
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
	
	mdlassert(pst->mdl,nIn == sizeof(*in));
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
	
	mdlassert(pst->mdl,nIn == sizeof(*in));
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
pstGravStep(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inGravStep *in = vin;

	mdlassert(pst->mdl,nIn == sizeof(struct inGravStep));
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_GRAVSTEP,vin,nIn);
		pstGravStep(pst->pstLower,vin,nIn,vout,pnOut);
		mdlGetReply(pst->mdl,pst->idUpper,vout,pnOut);
		}
	else {
		pkdGravStep(plcl->pkd,in->dEta);
		}
	if (pnOut) *pnOut = 0;
	}

void
pstAccelStep(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inAccelStep *in = vin;

	mdlassert(pst->mdl,nIn == sizeof(struct inAccelStep));
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_ACCELSTEP,vin,nIn);
		pstAccelStep(pst->pstLower,vin,nIn,vout,pnOut);
		mdlGetReply(pst->mdl,pst->idUpper,vout,pnOut);
		}
	else {
		pkdAccelStep(plcl->pkd,in->dEta,in->dVelFac,in->dAccFac,
					 in->bDoGravity,in->bEpsAcc,in->bSqrtPhi,in->dhMinOverSoft);
		}
	if (pnOut) *pnOut = 0;
	}

void
pstDensityStep(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inDensityStep *in = vin;
	
	mdlassert(pst->mdl,nIn == sizeof(*in));
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_DENSITYSTEP,vin,nIn);
		pstDensityStep(pst->pstLower,vin,nIn,vout,pnOut);
		mdlGetReply(pst->mdl,pst->idUpper,vout,pnOut);
		}
	else {
		pkdDensityStep(plcl->pkd,in->dEta,in->dRhoFac);
		}
	if (pnOut) *pnOut = 0;
	}

void pstDtToRung(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inDtToRung *in = vin;
	struct outDtToRung *out = vout;
	int iMaxRung, nMaxRung;
	int iMaxRungIdeal;

	mdlassert(pst->mdl,nIn == sizeof(*in));
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_DTTORUNG,vin,nIn);
		pstDtToRung(pst->pstLower,vin,nIn,vout,pnOut);
		iMaxRung = out->iMaxRung;
		nMaxRung = out->nMaxRung;
		iMaxRungIdeal = out->iMaxRungIdeal;
       
		mdlGetReply(pst->mdl,pst->idUpper,vout,pnOut);
		if(iMaxRung > out->iMaxRung) {
			out->iMaxRung = iMaxRung;
			out->nMaxRung = nMaxRung;
		        }
		else if (iMaxRung == out->iMaxRung) 
		        out->nMaxRung += nMaxRung;
		        
		if(iMaxRungIdeal > out->iMaxRungIdeal)
		    out->iMaxRungIdeal = iMaxRungIdeal;

	    }
	else {
		out->iMaxRung = pkdDtToRung(plcl->pkd, in->iRung,
					    in->dDelta, in->iMaxRung, in->bAll,
					    &(out->nMaxRung),
					    &(out->iMaxRungIdeal) );
		}
	if (pnOut) *pnOut = sizeof(*out);
	}

void pstInitDt(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inInitDt *in = vin;

	mdlassert(pst->mdl,nIn == sizeof(*in));
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_INITDT,vin,nIn);
		pstInitDt(pst->pstLower,vin,nIn,vout,pnOut);
		mdlGetReply(pst->mdl,pst->idUpper,vout,pnOut);
		}
	else {
		pkdInitDt(plcl->pkd,in->dDelta);
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

	mdlassert(pst->mdl,nIn == sizeof(struct inCoolVelocity));
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


void pstActiveExactType(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inActiveType *in = vin;
	int *pnActive = vout;
	
	mdlassert(pst->mdl,nIn == sizeof(*in));

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
	
	mdlassert(pst->mdl,nIn == sizeof(struct inActiveType));

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
	
	mdlassert(pst->mdl,nIn == sizeof(struct inActiveType));

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
	
	mdlassert(pst->mdl,nIn == sizeof(*in));

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
	
	mdlassert(pst->mdl,nIn == sizeof(*in));

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
	
	mdlassert(pst->mdl,nIn == sizeof(*in));
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
	
	mdlassert(pst->mdl,nIn == sizeof(*in));
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

void pstSetTypeFromFile(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inSetTypeFromFile *in = vin;
	struct outSetTypeFromFile *out = vout, outtmp;
	
	mdlassert(pst->mdl,nIn == sizeof(struct inSetTypeFromFile));

	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_SETTYPEFROMFILE,vin,nIn);
		pstSetTypeFromFile(pst->pstLower,vin,nIn,vout,pnOut);
		mdlGetReply(pst->mdl,pst->idUpper,&outtmp,pnOut);
		assert(out->niOrder == outtmp.niOrder);
		out->nSet += outtmp.nSet;
		out->nSetiGasOrder += outtmp.nSetiGasOrder;
		}
	else {
		pkdSetTypeFromFile(plcl->pkd,in->iSetMask,in->biGasOrder,in->file,&out->niOrder,&out->nSet,&out->nSetiGasOrder);
		}

	if (pnOut) *pnOut = sizeof(out);
	}

void pstSetParticleTypes(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inSetParticleTypes *in = vin;
	
	mdlassert(pst->mdl,nIn == sizeof(*in));
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

	mdlassert(pst->mdl,nIn == sizeof(struct inGrowMass));
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_GROWMASS,in,nIn);
		pstGrowMass(pst->pstLower,in,nIn,NULL,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
		}
	else {
		pkdGrowMass(plcl->pkd,in->nGrowMass,in->iGrowType,in->dDeltaM,in->dMinM,in->dMaxM);
		}
	if (pnOut) *pnOut = 0;
	}


void pstInitAccel(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	
	mdlassert(pst->mdl,nIn == 0);
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


#ifdef GASOLINE

void pstUpdateuDot(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inUpdateuDot *in = vin;
	struct outUpdateuDot *out = vout;
	struct outUpdateuDot outUp;

	mdlassert(pst->mdl,nIn == sizeof(struct inUpdateuDot));

	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_UPDATEUDOT,in,nIn);
		pstUpdateuDot(pst->pstLower,in,nIn,out,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,&outUp,NULL);

		out->SumTime += outUp.SumTime;
		out->nSum += outUp.nSum;
		if (outUp.MaxTime > out->MaxTime) out->MaxTime = outUp.MaxTime;
		}
	else {
		pkdUpdateuDot(plcl->pkd,in->duDelta,in->dTime,in->z,in->iGasModel,in->bUpdateState);
		out->Time = pkdGetTimer(plcl->pkd,1);
		out->MaxTime = out->Time;
		out->SumTime = out->Time;
		out->nSum = 1;
		}
	if (pnOut) *pnOut = sizeof(struct outUpdateuDot);
	}

void pstUpdateShockTracker(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inUpdateShockTracker *in = vin;

	mdlassert(pst->mdl,nIn == sizeof(struct inUpdateShockTracker));

	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_UPDATESHOCKTRACKER,vin,nIn);
		pstUpdateShockTracker(pst->pstLower,vin,nIn,vout,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,vout,NULL);
		}
	else {
		pkdUpdateShockTracker(plcl->pkd,in->dDelta,in->dShockTrackerA,in->dShockTrackerB);
		}
	if (pnOut) *pnOut = 0;
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
		case GASMODEL_ISOTHERMAL:
			pkdAdiabaticGasPressure(plcl->pkd, in->gammam1,in->gamma);
			break;
		case GASMODEL_COOLING:
#ifndef NOCOOLING
			pkdCoolingGasPressure(plcl->pkd, in->gammam1,in->gamma);
#endif
			break;
		case GASMODEL_GLASS:
#ifdef GLASS		  
			pkdGlassGasPressure(plcl->pkd, in);
#else
			mdlassert(pst->mdl,0);
#endif
			break;
			}
		}
	if (pnOut) *pnOut = 0;
	}


void pstGetDensityU(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	double *uMin = vin;
	
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_GETDENSITYU,vin,nIn);
		pstGetDensityU(pst->pstLower,vin,nIn,NULL,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
		}
	else {
  	        pkdGetDensityU(plcl->pkd, *uMin);
		}
	if (pnOut) *pnOut = 0;
}


void pstLowerSoundSpeed(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inLowerSoundSpeed *in = vin;
	
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_LOWERSOUNDSPEED,vin,nIn);
		pstLowerSoundSpeed(pst->pstLower,vin,nIn,NULL,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
		}
	else {
		pkdLowerSoundSpeed(plcl->pkd,in->dhMinOverSoft);
		}
	if (pnOut) *pnOut = 0;
	}


#ifndef NOCOOLING
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
		pkdInitEnergy(plcl->pkd,in->dTuFac,in->z,in->dTime);
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
#if defined(COOLDEBUG)
		(plcl->pkd->Cool)->mdl = plcl->pkd->mdl;
#endif
		clInitConstants((plcl->pkd->Cool),in->dGmPerCcUnit,in->dComovingGmPerCcUnit,
						in->dErgPerGmUnit,in->dSecUnit,in->dKpcUnit,in->CoolParam);
		CoolInitRatesTable((plcl->pkd->Cool),in->CoolParam);
		/* NOT DONE HERE: must be done before use: */
        /* CoolSetTime( (plcl->pkd->Cool), in->dTime, in->z  );*/
		}
	if (pnOut) *pnOut = 0;
	}

void pstCoolTableRead(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;

	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_COOLTABLEREAD,vin,nIn);
		pstCoolTableRead(pst->pstLower,vin,nIn,NULL,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
		}
	else {
		CoolTableRead( (plcl->pkd->Cool),nIn,vin);
		}
	if (pnOut) *pnOut = 0;
	}
#endif


void pstKickRhopred(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inKickRhopred *in = vin;

	mdlassert(pst->mdl,nIn == sizeof(struct inKickRhopred));
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
	
	mdlassert(pst->mdl,nIn == sizeof(*in));
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
	double *pdtMinGas = vout;

	mdlassert(pst->mdl,nIn == sizeof(struct inSphStep));
	if (pst->nLeaves > 1) {
   	        double dtMinGas;
		mdlReqService(pst->mdl,pst->idUpper,PST_SPHSTEP,in,nIn);
		pstSphStep(pst->pstLower,in,nIn,vout,pnOut);
		mdlGetReply(pst->mdl,pst->idUpper,&dtMinGas,pnOut);
		if (dtMinGas < *pdtMinGas) *pdtMinGas = dtMinGas;
		}
	else {
		pkdSphStep(plcl->pkd,in->dCosmoFac,in->dEtaCourant,in->dEtauDot,in->bViscosityLimitdt,pdtMinGas);
		}
	if (pnOut) *pnOut = sizeof(double);
	}

void pstSinkStep(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inSinkStep *in = vin;

	mdlassert(pst->mdl,nIn == sizeof(struct inSinkStep));
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_SINKSTEP,in,nIn);
		pstSinkStep(pst->pstLower,in,nIn,NULL,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
		}
	else {
		pkdSinkStep(plcl->pkd,in->dtMax);
		}
	if (pnOut) *pnOut = 0;
	}

void pstSphViscosityLimiter(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inSphViscosityLimiter *in = vin;

	mdlassert(pst->mdl,nIn == sizeof(struct inSphViscosityLimiter));
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_SPHVISCOSITYLIMITER,in,nIn);
		pstSphViscosityLimiter(pst->pstLower,in,nIn,NULL,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
		}
	else {
		pkdSphViscosityLimiter(plcl->pkd,in->bOn,in->bShockTracker);
		}
	if (pnOut) *pnOut = 0;
	}

void
pstDensCheck(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inDensCheck *in = vin;
	struct outDensCheck *out = vout, tmp;
	
	mdlassert(pst->mdl,nIn == sizeof(*in));
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

	mdlassert(pst->mdl,nIn == sizeof(struct inRandomVelocities));
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
		mdlassert(pst->mdl,ptmp != NULL);
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
pstGetNParts(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
    /*LCL *plcl = pst->plcl; -- not used: DCR 12/19/02*/
    struct outGetNParts *out = vout;
    /*int i; -- not used: DCR 12/19/02*/
    
    if(pst->nLeaves > 1) {
		struct outGetNParts outtmp;
		mdlReqService(pst->mdl,pst->idUpper,PST_GETNPARTS,vin,nIn);
		pstGetNParts(pst->pstLower,vin,nIn,vout,pnOut);
		mdlGetReply(pst->mdl,pst->idUpper,(void *) &outtmp,pnOut);

		out->n += outtmp.n;
		out->nGas += outtmp.nGas;
		out->nDark += outtmp.nDark;
		out->nStar += outtmp.nStar;
		if (outtmp.iMaxOrderGas > out->iMaxOrderGas) out->iMaxOrderGas = outtmp.iMaxOrderGas;
		if (outtmp.iMaxOrderDark > out->iMaxOrderDark) out->iMaxOrderDark = outtmp.iMaxOrderDark;
		if (outtmp.iMaxOrderStar > out->iMaxOrderStar) out->iMaxOrderStar = outtmp.iMaxOrderStar;
		}
    else {
		pkdGetNParts(pst->plcl->pkd, out);
		}
    if(pnOut) *pnOut = sizeof(struct outGetNParts);
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
					 in->nMaxOrder, in->nMaxOrderGas, in->nMaxOrderDark);
		}
    if(pnOut) *pnOut = 0;
    }

#ifdef SPECIAL_PARTICLES

void
pstGetSpecialParticles(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inGetSpecial *in = vin;
	struct outGetSpecial local,*out = vout;
	int i;

	mdlassert(pst->mdl,nIn == sizeof(*in));
	for (i=0;i<in->nSpecial;i++)
		local.sInfo[i].iOrder = out->sInfo[i].iOrder = -1;
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_GETSPECIALPARTICLES,vin,nIn);
		pstGetSpecialParticles(pst->pstLower,vin,nIn,vout,NULL);
		local = *out;
		mdlGetReply(pst->mdl,pst->idUpper,vout,NULL);
		for (i=0;i<in->nSpecial;i++)
			if (local.sInfo[i].iOrder != -1)
				out->sInfo[i] = local.sInfo[i];
		}
	else {
		pkdGetSpecialParticles(plcl->pkd,in->nSpecial,in->iId,&in->mInfo,
							   out->sInfo);
		}
	if (pnOut) *pnOut = sizeof(*out);
	}

void
pstDoSpecialParticles(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inDoSpecial *in = vin;
	struct outDoSpecial *out = vout;

	mdlassert(pst->mdl,nIn == sizeof(*in));
	if (pst->nLeaves > 1) {
		double aFrame[3];
		int k;
		mdlReqService(pst->mdl,pst->idUpper,PST_DOSPECIALPARTICLES,
			      vin,nIn);
		pstDoSpecialParticles(pst->pstLower,vin,nIn,vout,pnOut);
		if (in->mInfo.bNonInertial) for (k=0;k<3;k++) aFrame[k] = out->aFrame[k];
		mdlGetReply(pst->mdl,pst->idUpper,vout,pnOut);
		if (in->mInfo.bNonInertial) for (k=0;k<3;k++) out->aFrame[k] += aFrame[k];
		}
	else {
		pkdDoSpecialParticles(plcl->pkd,in->nSpecial,&in->mInfo,
							  in->sData,in->sInfo,out->aFrame);
		}
	if (pnOut) *pnOut = sizeof(*out);
	}

#endif /* SPECIAL_PARTICLES */

#ifdef COLLISIONS

void pstSetBall(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inSetBall *in = vin;

	mdlassert(pst->mdl,nIn == sizeof(struct inSetBall));
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_SETBALL,in,nIn);
		pstSetBall(pst->pstLower,in,nIn,NULL,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
		}
	else {
		pkdSetBall(plcl->pkd,in->dDelta,in->dBallVelFact);
		}
	if (pnOut) *pnOut = 0;
	}
  

void
pstNumRejects(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct outNumRejects *out = vout;
	int nRej;

	mdlassert(pst->mdl,nIn == 0);
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

	mdlassert(pst->mdl,nIn == sizeof(struct inReadSS));
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

	mdlassert(pst->mdl,nIn == sizeof(struct inWriteSS));
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
pstKickUnifGrav(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inKickUnifGrav *in = vin;

	mdlassert(pst->mdl,nIn == sizeof(struct inKickUnifGrav));
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

	mdlassert(pst->mdl,nIn == 0);
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

	mdlassert(pst->mdl,nIn == sizeof(*in));
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

	mdlassert(pst->mdl,nIn == 0);
	out->dt = DBL_MAX;
	out->iOrder1 = out->iOrder2 = -1;
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

	mdlassert(pst->mdl,nIn == sizeof(*in));
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
pstFindTightestBinary(PST pst,void *vin,int nIn,void *vout,int *pnOut) 
{
  LCL *plcl = pst->plcl;
  struct inBinary *in = vin;
  struct outBinary local,*out = vout;
  
  mdlassert(pst->mdl,nIn == sizeof(*in));
  local.dBindEn = out->dBindEn = 0;
  local.n = out->n = 0;
  if (pst->nLeaves > 1) {
    mdlReqService(pst->mdl,pst->idUpper,PST_FINDBINARY,vin,nIn);
    pstFindTightestBinary(pst->pstLower,vin,nIn,&local,NULL);
    mdlGetReply(pst->mdl,pst->idUpper,vout,NULL);
    if (local.dBindEn < out->dBindEn) *out = local;
  } else {
    pkdFindTightestBinary(plcl->pkd,&out->dBindEn,&out->iOrder1,
			  &out->iOrder2,&out->n);
  }
  if (pnOut) *pnOut = sizeof(*out);
}

void
pstMergeBinary(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
  LCL *plcl = pst->plcl;
  struct inMrgBnry *in = vin;
  struct outMrgBnry local,*out = vout;

  mdlassert(pst->mdl,nIn == sizeof(*in));
  local.n=0;
  if (pst->nLeaves > 1) {
    mdlReqService(pst->mdl,pst->idUpper,PST_MERGEBINARY,vin,nIn);
    pstMergeBinary(pst->pstLower,vin,nIn,&local,NULL);
    mdlGetReply(pst->mdl,pst->idUpper,vout,NULL);
    if (local.n) *out = local;
  } else {
    PKD pkd = plcl->pkd;
    if (in->c1.id.iPid == pkd->idSelf || in->c2.id.iPid == pkd->idSelf) {
#ifdef SLIDING_PATCH
      pkd->PP = &in->PP;
      pkd->dTime = in->dTime;
#endif
      pkdMergeBinary(pkd,(const COLLIDER *) &in->c1,(const COLLIDER *) &in->c2,&out->cOut,in->bPeriodic,in->dBaseStep,in->dTimeNow,in->iTime0,in->dDensity,&out->n);
    }
    if (pnOut) *pnOut = sizeof(*out);
  }
}

void
pstDoCollision(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inDoCollision *in = vin;
	struct outDoCollision local,*out = vout;

	mdlassert(pst->mdl,nIn == sizeof(*in));
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
			pkd->dTime = in->dTime;
			pkd->PP = &in->PP;
#endif
#ifdef AGGS
			pkdAggsDoCollision(pkd,in->dt,&in->Collider1,&in->Collider2,
							   in->bPeriodic,&in->CP,in->iAggNewIdx,
							   &out->iOutcome,&out->dT,out->Out,&out->nOut);
#else
			pkdDoCollision(pkd,in->dt,&in->Collider1,&in->Collider2,
						   in->bPeriodic,in->iTime0,in->dBaseTime,
					in->dTime,&in->CP,&out->iOutcome,&out->dT,
						   out->Out,&out->nOut);
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
	
	mdlassert(pst->mdl,nIn == sizeof(*in));
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

#ifdef OLD_KEPLER
void
pstQQCalcBound(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct outCalcBound *out = vout;
	struct outCalcBound outBnd;
	int j;
	
	mdlassert(pst->mdl,nIn == 0);
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

	mdlassert(pst->mdl,nIn == 0);
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
		_pstRootSplit(pst,d+3,dMass,1,1,1); /* Scary! */
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
	
	mdlassert(pst->mdl,nIn == sizeof(struct inBuildTree));
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

	mdlassert(pst->mdl,nIn == sizeof(struct inSmooth));
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_QQSMOOTH,in,nIn);
		pstQQSmooth(pst->pstLower,in,nIn,NULL,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
		}
	else {
		LCL *plcl = pst->plcl;
		SMX smx;

		smInitialize(&smx,plcl->pkd,&in->smf,in->nSmooth,in->bPeriodic,
					 in->bSymmetric,in->iSmoothType,0,0.0);
		smQQSmooth(smx,&in->smf);
		smFinish(smx,&in->smf);
		}
	if (pnOut) *pnOut = 0;
	}
#endif /* OLD_KEPLER */

#ifdef SLIDING_PATCH
/* These are the subroutines for randomizing masses which are so large
   they will violate the shearing patch Hamiltonian. */ 
void 
pstFindLargeMasses(PST pst,void *vin,int nIn,void *vOut,int *pnOut)
{
    struct inLargeMass *in = vin;
    struct outLargeMass local,*out = vOut;
    int i;

	mdlassert(pst->mdl,nIn == sizeof(struct inLargeMass));
	local.n = out->n = 0;
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_FINDLM,in,nIn);
		pstFindLargeMasses(pst->pstLower,vin,nIn,vOut,NULL);
		local = *out;		
		mdlGetReply(pst->mdl,pst->idUpper,vOut,NULL);
		mdlprintf(pst->mdl,"After GetReply\n");
		if (local.n > 0) {
		    /* assert((local.n+out->n) < MAXLARGEMASS);	*/	    
		    for (i=0;i<local.n;i++) {
			out->dRadius[out->n + i]=local.dRadius[i];
			out->p[out->n + i]=local.p[i];		    
			}
		    out->n += local.n;
		    }		
	    }	
	else {
	    pkdFindLargeMasses(pst->plcl->pkd,in->dMass,in->dCentMass,in->dOrbRad,in->fNumHillSphere,out->p,out->dRadius,&out->n);
	    }
	
	if (pnOut) *pnOut = sizeof(struct outLargeMass);
	

    }

void
pstGetNeighborParticles(PST pst,void *vin,int nIn,void *vOut,int *pnOut)
{
    struct inGetNeighbors *in = vin;
    struct outGetNeighbors local,*out = vOut;
    int i;    
    
        mdlassert(pst->mdl,nIn == sizeof(struct inGetNeighbors));
	local.n = out->n = 0;
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_GETNEIGHBORS,in,nIn);
		pstGetNeighborParticles(pst->pstLower,in,nIn,&local,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,vOut,NULL);
		if (local.n > 0) {
		    assert((local.n + out->n) < MAXNEIGHBORS);		    
		    for (i=0;i<local.n;i++) {
			out->dSep2[i+out->n]=local.dSep2[i];
			out->p[i+out->n]=local.p[i];
			}
		    out->n += local.n;
		    }
		
	    }
	else {
	  pst->plcl->pkd->PP = &in->PP; /* PKD struct already has space for patch params pointer, so use that */
	    pkdGetNeighborParticles(pst->plcl->pkd,in->x,in->dDist*in->dDist,in->id,in->dTime,out->p,out->dSep2,&out->n);
	    }
	if (pnOut) *pnOut = sizeof(struct outGetNeighbors);

    }


#endif /* SLIDING_PATCH */

#endif /* COLLISIONS */

void 
pstMoveParticle(PST pst,void *vin,int nIn,void *vOut,int *pnOut)
{
    struct inMoveParticle *in = vin;

        mdlassert(pst->mdl,nIn == sizeof(struct inMoveParticle));
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_MOVEPART,in,nIn);
		pstMoveParticle(pst->pstLower,in,nIn,NULL,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
		}
	else {
#ifdef SLIDING_PATCH
	  pst->plcl->pkd->PP = &in->PP;
#endif
	    pkdMoveParticle(pst->plcl->pkd,in->dOrigin,in->dRelx,in->p.iOrder);    
	    }
	
    }
 


#ifdef AGGS

void pstAggsFind(PST pst,void *vIn,int nIn,void *vOut,int *pnOut)
{
	struct outAggsFind local,*out = vOut;

	mdlassert(pst->mdl,nIn == 0);
	out->iMaxIdx = -1;
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_AGGSFIND,NULL,0);
		pstAggsFind(pst->pstLower,NULL,0,vOut,NULL);
		local = *out;
		mdlGetReply(pst->mdl,pst->idUpper,vOut,NULL);
		if (local.iMaxIdx > out->iMaxIdx) *out = local;
		}
	else
		pkdAggsFind(pst->plcl->pkd,&out->iMaxIdx);
	if (pnOut) *pnOut = sizeof(*out);
	}

void pstAggsConfirm(PST pst,void *vIn,int nIn,void *vOut,int *pnOut)
{
	/* called by msrAggsFind() */

	struct inAggsConfirm *in = vIn;
	struct outAggsConfirm local,*out = vOut;

	mdlassert(pst->mdl,nIn == sizeof(*in));
	out->bAssigned = 0;
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_AGGSCONFIRM,vIn,nIn);
		pstAggsConfirm(pst->pstLower,vIn,nIn,vOut,NULL);
		local = *out;
		mdlGetReply(pst->mdl,pst->idUpper,vOut,NULL);
		if (local.bAssigned) *out = local;
		}
	else
		pkdAggsConfirm(pst->plcl->pkd,in->iAggIdx,&out->bAssigned);
	if (pnOut) *pnOut = sizeof(*out);
	}

void pstAggsMerge(PST pst,void *vIn,int nIn,void *vOut,int *pnOut)
{
	/* called by msrAggsMerge() */

	struct inAggsMerge *in = vIn;

	mdlassert(pst->mdl,nIn == sizeof(*in));
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_AGGSMERGE,vIn,nIn);
		pstAggsMerge(pst->pstLower,vIn,nIn,NULL,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
		}
	else
		pkdAggsMerge(pst->plcl->pkd,in->iOldIdx,in->iNewIdx);
	if (pnOut) *pnOut = 0;
	}

void pstAggsBackDrift(PST pst,void *vIn,int nIn,void *vOut,int *pnOut)
{
	struct inAggsBackDrift *in = vIn;

	mdlassert(pst->mdl,nIn == sizeof(*in));
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_AGGSBACKDRIFT,vIn,nIn);
		pstAggsBackDrift(pst->pstLower,vIn,nIn,NULL,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
		}
	else
		pkdAggsBackDrift(pst->plcl->pkd,in->iAggIdx,in->dt);
	if (pnOut) *pnOut = 0;
	}

void pstAggsGetCOM(PST pst,void* vIn,int nIn,void* vOut,int* pnOut)
{
	struct inAggsGetCOM *in = vIn;
	struct outAggsGetCOM local,*out = vOut;
	int k;

	mdlassert(pst->mdl,nIn == sizeof(*in));
	out->m = 0.0;
	for (k=0;k<3;k++) out->mr[k] = out->mv[k] = 0.0;
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_AGGSGETCOM,vIn,nIn);
		pstAggsGetCOM(pst->pstLower,vIn,nIn,vOut,NULL);
		local = *out;
		mdlGetReply(pst->mdl,pst->idUpper,vOut,NULL);
		out->m += local.m;
		for (k=0;k<3;k++) {
			out->mr[k] += local.mr[k];
			out->mv[k] += local.mv[k];
			}
		}
	else
		pkdAggsGetCOM(pst->plcl->pkd,in->iAggIdx,&out->m,out->mr,out->mv);
	if (pnOut) *pnOut = sizeof(*out);
	}

void pstAggsGetAxesAndSpin(PST pst,void* vIn,int nIn,void* vOut,int* pnOut)
{
	struct inAggsGetAxesAndSpin *in = vIn;
	struct outAggsGetAxesAndSpin local,*out = vOut;

	mdlassert(pst->mdl,nIn == sizeof(*in));
	out->I[0][0] = out->I[0][1] = out->I[0][2] = out->I[1][1] =
	out->I[1][2] = out->I[2][2] = out->L[0] = out->L[1] = out->L[2] = 0.0;
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_AGGSGETAXESANDSPIN,vIn,nIn);
		pstAggsGetAxesAndSpin(pst->pstLower,vIn,nIn,vOut,NULL);
		local = *out;
		mdlGetReply(pst->mdl,pst->idUpper,vOut,NULL);
		out->I[0][0] += local.I[0][0];
		out->I[0][1] += local.I[0][1];
		out->I[0][2] += local.I[0][2];
		out->I[1][1] += local.I[1][1];
		out->I[1][2] += local.I[1][2];
		out->I[2][2] += local.I[2][2];
		out->L[0] += local.L[0];
		out->L[1] += local.L[1];
		out->L[2] += local.L[2];
		}
	else
		pkdAggsGetAxesAndSpin(pst->plcl->pkd,in->iAggIdx,in->r_com,in->v_com,
							  out->I,out->L);
	if (pnOut) *pnOut = sizeof(*out);
	}

void pstAggsSetBodyPos(PST pst,void* vIn,int nIn,void* vOut,int* pnOut)
{
	struct inAggsSetBodyPos *in = vIn;

	mdlassert(pst->mdl,nIn == sizeof(*in));
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_AGGSSETBODYPOS,vIn,nIn);
		pstAggsSetBodyPos(pst->pstLower,vIn,nIn,NULL,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
		}
	else
		pkdAggsSetBodyPos(pst->plcl->pkd,in->iAggIdx,in->spaceToBody);
	if (pnOut) *pnOut = 0;
	}

void pstAggsSetSpacePos(PST pst,void* vIn,int nIn,void* vOut,int* pnOut)
{
	struct inAggsSetSpacePos *in = vIn;

	mdlassert(pst->mdl,nIn == sizeof(*in));
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_AGGSSETSPACEPOS,vIn,nIn);
		pstAggsSetSpacePos(pst->pstLower,vIn,nIn,NULL,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
		}
	else
		pkdAggsSetSpacePos(pst->plcl->pkd,in->iAggIdx,in->r_com,in->lambda);
	if (pnOut) *pnOut = 0;
	}

void pstAggsSetSpaceVel(PST pst,void* vIn,int nIn,void* vOut,int* pnOut)
{
	struct inAggsSetSpaceVel *in = vIn;

	mdlassert(pst->mdl,nIn == sizeof(*in));
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_AGGSSETSPACEVEL,vIn,nIn);
		pstAggsSetSpaceVel(pst->pstLower,vIn,nIn,NULL,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
		}
	else
		pkdAggsSetSpaceVel(pst->plcl->pkd,in->iAggIdx,in->v_com,in->omega,
						   in->lambda);
	if (pnOut) *pnOut = 0;
	}

void pstAggsSetSpaceSpins(PST pst,void* vIn,int nIn,void* vOut,int* pnOut)
{
	struct inAggsSetSpaceSpins *in = vIn;

	mdlassert(pst->mdl,nIn == sizeof(*in));
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_AGGSSETSPACESPINS,vIn,nIn);
		pstAggsSetSpaceSpins(pst->pstLower,vIn,nIn,NULL,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
		}
	else
		pkdAggsSetSpaceSpins(pst->plcl->pkd,in->iAggIdx,in->omega);
	if (pnOut) *pnOut = 0;
	}

void pstAggsDelete(PST pst,void* vIn,int nIn,void* vOut,int* pnOut)
{
	struct inAggsDelete *in = vIn;
	struct outAggsDelete local,*out = vOut;

	mdlassert(pst->mdl,nIn == sizeof(*in));
	out->bFound = 0;
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_AGGSDELETE,vIn,nIn);
		pstAggsDelete(pst->pstLower,vIn,nIn,NULL,NULL);
		local = *out;
		mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
		if (local.bFound == 1) out->bFound = 1;
		}
	else
		pkdAggsDelete(pst->plcl->pkd,in->iAggIdx,&out->bFound);
	if (pnOut) *pnOut = sizeof(*out);
	}

void pstAggsGetAccel(PST pst,void* vIn,int nIn,void* vOut,int* pnOut)
{
	struct inAggsGetAccel *in = vIn;
	struct outAggsGetAccel local,*out = vOut;
	int k;

	mdlassert(pst->mdl,nIn == sizeof(*in));
	out->m = 0.0;
	for (k=0;k<3;k++) out->ma[k] = 0.0;
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_AGGSGETACCEL,vIn,nIn);
		pstAggsGetAccel(pst->pstLower,vIn,nIn,vOut,NULL);
		local = *out;
		mdlGetReply(pst->mdl,pst->idUpper,vOut,NULL);
		out->m += local.m;
		for (k=0;k<3;k++) out->ma[k] += local.ma[k];
		}
	else
		pkdAggsGetAccel(pst->plcl->pkd,in->iAggIdx,&out->m,out->ma);
	if (pnOut) *pnOut = sizeof(*out);
	}

void pstAggsCheckStress(PST pst,void* vIn,int nIn,void* vOut,int* pnOut)
{
	struct inAggsCheckStress *in = vIn;
	struct outAggsCheckStress local,*out = vOut;

	mdlassert(pst->mdl,nIn == sizeof(*in));
	out->nLost = out->nLeft = 0;
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_AGGSCHECKSTRESS,vIn,nIn);
		pstAggsCheckStress(pst->pstLower,vIn,nIn,vOut,NULL);
		local = *out;
		mdlGetReply(pst->mdl,pst->idUpper,vOut,NULL);
		out->nLost += local.nLost;
		out->nLeft += local.nLeft;
		}
	else
		pkdAggsCheckStress(pst->plcl->pkd,in->iAggIdx,in->r_com,in->a_com,
						   in->omega,in->fTensileStrength,in->fShearStrength,
						   &out->nLost,&out->nLeft);
	if (pnOut) *pnOut = sizeof(*out);
	}

void pstAggsGetTorque(PST pst,void* vIn,int nIn,void* vOut,int* pnOut)
{
	struct inAggsGetTorque *in = vIn;
	struct outAggsGetTorque local,*out = vOut;
	int k;

	mdlassert(pst->mdl,nIn == sizeof(*in));
	for (k=0;k<3;k++) out->torque[k] = 0.0;
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_AGGSGETTORQUE,vIn,nIn);
		pstAggsGetTorque(pst->pstLower,vIn,nIn,vOut,NULL);
		local = *out;
		mdlGetReply(pst->mdl,pst->idUpper,vOut,NULL);
		for (k=0;k<3;k++) out->torque[k] += local.torque[k];
		}
	else
		pkdAggsGetTorque(pst->plcl->pkd,in->iAggIdx,in->r_com,in->a_com,
						 out->torque);
	if (pnOut) *pnOut = sizeof(*out);
	}

#endif /* AGGS */

#ifdef SLIDING_PATCH

void pstRandAzWrap(PST pst,void *vIn,int nIn,void *vOut,int *pnOut)
{
  LCL *plcl = pst->plcl;
  struct inRandAzWrap *in = vIn;
  struct outRandAzWrap *out = vOut;
  int nRand;

  mdlassert(pst->mdl,nIn == sizeof(*in));
  if (pst->nLeaves > 1) {
	mdlReqService(pst->mdl,pst->idUpper,PST_RANDAZWRAP,vIn,nIn);
	pstRandAzWrap(pst->pstLower,vIn,nIn,vOut,pnOut);
	nRand = out->nRandomized;
	mdlGetReply(pst->mdl,pst->idUpper,vOut,pnOut);
	out->nRandomized += nRand;
  }
  else {
	/*
	** Like pkdPatch(), could just pass this directly, but since the
	** PKD struct already has space for them, we'll just use that.
	*/
	plcl->pkd->PP = &in->PP;
	out->nRandomized = pkdRandAzWrap(plcl->pkd);
  }
  if (pnOut) *pnOut = sizeof(*out);
}

#endif /* SLIDING_PATCH */

#ifdef RUBBLE_ZML

void 
pstDustBinsGetMass(PST pst,void* vIn,int nIn,void* vOut,int* pnOut)
{
	struct inDustBinsGetMass *in = vIn;
	struct outDustBinsGetMass local,*out = vOut;
	int i;

	mdlassert(pst->mdl,nIn == sizeof(*in));
	for (i=0;i<in->DB.nDustBins;i++) out->aDustBinsMassLoss[i] = 0.0;
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_DUSTBINSGETMASS,vIn,nIn);
		pstDustBinsGetMass(pst->pstLower,vIn,nIn,vOut,NULL);
		local = *out;
		mdlGetReply(pst->mdl,pst->idUpper,vOut,NULL);
		for (i=0;i<in->DB.nDustBins;i++)
			out->aDustBinsMassLoss[i] += local.aDustBinsMassLoss[i];
		}
	else
		pkdDustBinsGetMass(pst->plcl->pkd,&in->DB,in->aDustBins,in->dTimeInt,in->dCentMass,
						   out->aDustBinsMassLoss);
	if (pnOut) *pnOut = sizeof(*out);
	}				   

void 
pstDustBinsApply(PST pst,void* vIn,int nIn,void* vOut,int* pnOut)
{
	struct inDustBinsApply *in = vIn;

	mdlassert(pst->mdl,nIn == sizeof(*in));
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_DUSTBINSAPPLY,vIn,nIn);
		pstDustBinsApply(pst->pstLower,vIn,nIn,NULL,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
		}
	else
		pkdDustBinsApply(pst->plcl->pkd,in->dCentMass,in->aMassIncrFrac,in->nBins);
	}				   

void 
pstRubbleResetColFlag(PST pst,void* vIn,int nIn,void* vOut,int* pnOut)
{
	mdlassert(pst->mdl,nIn == 0);
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_RUBBLERESETCOLFLAG,NULL,0);
		pstRubbleResetColFlag(pst->pstLower,NULL,0,NULL,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
		}
	else
		pkdRubbleResetColFlag(pst->plcl->pkd);
}

void 
pstRubbleCheckForKDKRestart(PST pst,void* vIn,int nIn,void* vOut,int* pnOut)
{
	struct outRubbleCheckForKDKRestart local,*out = vOut;

	mdlassert(pst->mdl,nIn == 0);
	out->bRestart = 0;
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_RUBBLECHECKFORKDKRESTART,NULL,0);
		pstRubbleCheckForKDKRestart(pst->pstLower,NULL,0,vOut,NULL);
		local = *out;
		mdlGetReply(pst->mdl,pst->idUpper,vOut,NULL);
		if (local.bRestart) out->bRestart = 1;
		}
	else
		out->bRestart = pkdRubbleCheckForKDKRestart(pst->plcl->pkd);
}

void 
pstRubbleStep(PST pst,void* vIn,int nIn,void* vOut,int* pnOut)
{
	struct inRubbleStep *in = vIn;

	mdlassert(pst->mdl,nIn == sizeof(*in));
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_RUBBLESTEP,vIn,nIn);
		pstRubbleStep(pst->pstLower,vIn,nIn,NULL,NULL);
		mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
		}
	else
		pkdRubbleStep(pst->plcl->pkd,in->dMaxStep,in->dMinStep);
}

void 
pstRubCleanup(PST pst,void* vIn,int nIn,void* vOut,int* pnOut)
{
	struct inRubCleanup *in = vIn;
	struct outRubCleanup local,*out = vOut;
	int i;

	mdlassert(pst->mdl,nIn == sizeof(*in));

	for (i=0;i<in->DB.nDustBins;i++) out->aDustBinsMassGain[i] = 0.0;
	out->dDustBinsRubTrash = 0.0;
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_RUBCLEANUP,vIn,nIn);
		pstRubCleanup(pst->pstLower,vIn,nIn,vOut,NULL);
		local = *out;
		mdlGetReply(pst->mdl,pst->idUpper,vOut,NULL);
		for (i=0;i<in->DB.nDustBins;i++) 
			out->aDustBinsMassGain[i] += local.aDustBinsMassGain[i];
		out->dDustBinsRubTrash += local.dDustBinsRubTrash;
	}
	else
		pkdRubCleanup(pst->plcl->pkd,in->iColor,&in->DB,in->dCentMass,
					  out->aDustBinsMassGain,&out->dDustBinsRubTrash);
	if (pnOut) *pnOut = sizeof(*out);
}

void 
pstRubInterpCleanup(PST pst,void* vIn,int nIn,void* vOut,int* pnOut)
{
	struct inRubInterpCleanup *in = vIn;
	struct outRubInterpCleanup local,*out = vOut;

	mdlassert(pst->mdl,nIn == sizeof(*in));
	if (pst->nLeaves > 1) {
		mdlReqService(pst->mdl,pst->idUpper,PST_RUBINTERPCLEANUP,vIn,nIn);
		pstRubInterpCleanup(pst->pstLower,vIn,nIn,vOut,NULL);
		local = *out;
		mdlGetReply(pst->mdl,pst->idUpper,vOut,NULL);
		if (local.dDustBinsInterpMass > 0.0) {
			assert(out->dDustBinsInterpMass == 0.0); /* sanity check: only one processor can report dust mass */
			*out = local;
			}
		}
	else
		pkdRubInterpCleanup(pst->plcl->pkd,&in->DB,in->dCentMass,in->iOrder,
							&out->iBin,&out->dDustBinsInterpMass);
	if (pnOut) *pnOut = sizeof(*out);
}

#endif /* RUBBLE_ZML */

void 
pstClearTimer(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	struct inClearTimer *in = vin;

	mdlassert(pst->mdl,nIn == sizeof(struct inClearTimer));
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

void pstMassInR(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
        LCL *plcl = pst->plcl;
        struct inMassInR *in = vin;
        struct outMassInR *out = vout;
        struct outMassInR outLcl;

        mdlassert(pst->mdl,nIn == sizeof(struct inMassInR));
        if (pst->nLeaves > 1) {
                int k;
                double dTMass;

                mdlReqService(pst->mdl,pst->idUpper,PST_MASSINR,in,nIn);
                pstMassInR(pst->pstLower,in,nIn,out,NULL);
                mdlGetReply(pst->mdl,pst->idUpper,&outLcl,NULL);
                dTMass = out->dMass + outLcl.dMass;
                if(dTMass == 0.0) {
                    for(k = 0; k < 3; k++)
                        out->com[k] = 0.0;
                    }
                else {
                    for(k = 0; k < 3; k++)
                        out->com[k] = (out->com[k]*out->dMass +
                            outLcl.com[k]*outLcl.dMass)/dTMass;
                    }
                out->dMass = dTMass;
                }
        else {
                pkdMassInR(plcl->pkd, in->R, &out->dMass, out->com);
                }
        if (pnOut) *pnOut = sizeof(struct outMassInR);
        }

void pstInitRotBar(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
        LCL *plcl = pst->plcl;
        struct inRotBar *in = vin;

        mdlassert(pst->mdl,nIn == sizeof(struct inRotBar));
        if (pst->nLeaves > 1) {
                mdlReqService(pst->mdl,pst->idUpper,PST_ROTBARINIT,in,nIn);
                pstInitRotBar(pst->pstLower,in,nIn,NULL,NULL);
                mdlGetReply(pst->mdl,pst->idUpper,NULL,NULL);
                }
        else {
            pkdInitRotBar(plcl->pkd, &(in->rotbar));
            }
        if (pnOut) *pnOut = 0;
        }

#ifdef NEED_VPRED
#ifdef GASOLINE
void pstKickVpred(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inKickVpred *in = vin;
	struct outKick *out = vout;
	struct outKick outUp;

	mdlassert(pst->mdl,nIn == sizeof(struct inKickVpred));
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
#else
void pstKickVpred(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	LCL *plcl = pst->plcl;
	struct inKickVpred *in = vin;

	mdlassert(pst->mdl,nIn == sizeof(struct inKickVpred));
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
#endif /* NEED_VPRED */

void
pstFormSinks(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	struct inFormSinks *in = vin;
	struct outFormSinks *out = vout;

	mdlassert(pst->mdl,nIn == sizeof(struct inFormSinks));
	if (pst->nLeaves > 1) {
	    struct outFormSinks fsStats;
	    
		mdlReqService(pst->mdl,pst->idUpper,PST_FORMSINKS,in,nIn);
		pstFormSinks(pst->pstLower,in,nIn,vout,pnOut);
		mdlGetReply(pst->mdl,pst->idUpper,&fsStats,NULL);
		out->nCandidates += fsStats.nCandidates;
		if (fsStats.Jvalmin < out->Jvalmin) out->Jvalmin = fsStats.Jvalmin;
		}
	else {
	    pkdFormSinks(pst->plcl->pkd,in->bJeans,in->dJConst2,in->bDensity,in->dDensityCut,
			 in->dTime,in->iKickRung, in->bSimple, &out->nCandidates, &out->Jvalmin);
		}
	if (pnOut) *pnOut = sizeof(struct outFormSinks);
	}

#ifdef STARFORM
void
pstFormStars(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	struct inFormStars *in = vin;
	struct outFormStars *out = vout;

	mdlassert(pst->mdl,nIn == sizeof(struct inFormStars));
	if (pst->nLeaves > 1) {
	    struct outFormStars fsStats;
	    
		mdlReqService(pst->mdl,pst->idUpper,PST_FORMSTARS,in,nIn);
		pstFormStars(pst->pstLower,in,nIn,vout,pnOut);
		mdlGetReply(pst->mdl,pst->idUpper,&fsStats,NULL);
		out->nFormed += fsStats.nFormed;
		out->nDeleted += fsStats.nDeleted;
		out->dMassFormed += fsStats.dMassFormed;
		}
	else {
		pkdFormStars(pst->plcl->pkd,&in->stfm, in->dTime,
			     &out->nFormed, &out->dMassFormed, &out->nDeleted);
		}
	if (pnOut) *pnOut = sizeof(struct outFormStars);
	}

void
pstFeedback(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	struct inFeedback *in = vin;
	struct outFeedback *out = vout;

	mdlassert(pst->mdl,nIn == sizeof(struct inFeedback));
	if (pst->nLeaves > 1) {
	    FBEffects fbTotals[FB_NFEEDBACKS];
	    int i;

	    mdlReqService(pst->mdl,pst->idUpper,PST_FEEDBACK,in,nIn);
	    pstFeedback(pst->pstLower,in,nIn, vout, pnOut);
	    mdlGetReply(pst->mdl,pst->idUpper,fbTotals, NULL);
	    for(i = 0; i < FB_NFEEDBACKS; i++){
		out->fbTotals[i].dMassLoss += fbTotals[i].dMassLoss;
		out->fbTotals[i].dEnergy += fbTotals[i].dEnergy;
		out->fbTotals[i].dMetals += fbTotals[i].dMetals;
		out->fbTotals[i].dMIron += fbTotals[i].dMIron;
		out->fbTotals[i].dMOxygen += fbTotals[i].dMOxygen;
		}
	    }
	else {
		pkdFeedback(pst->plcl->pkd,&in->fb, &in->sn, in->dTime,
			    in->dDelta, out->fbTotals);
		}
	if (pnOut) *pnOut = FB_NFEEDBACKS*sizeof(FBEffects);
	}

void
pstInitStarLog(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
    if (pst->nLeaves > 1) {
	mdlReqService(pst->mdl,pst->idUpper,PST_INITSTARLOG,NULL,0);
	pstInitStarLog(pst->pstLower,NULL,0, NULL, NULL);
	mdlGetReply(pst->mdl,pst->idUpper,NULL, NULL);
    }
    else {
	pkdStarLogInit(pst->plcl->pkd);
	}
    }

void
pstFlushStarLog(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
    struct inFlushStarLog *in = vin;

    mdlassert(pst->mdl,nIn == sizeof(struct inFlushStarLog));
    if (pst->nLeaves > 1) {
	/*
	 * N.B. NOT parallel
	 */
	pstFlushStarLog(pst->pstLower, in, nIn, NULL, NULL);
	mdlReqService(pst->mdl,pst->idUpper,PST_FLUSHSTARLOG, in, nIn);
	mdlGetReply(pst->mdl,pst->idUpper,NULL, NULL);
    }
    else {
	pkdStarLogFlush(pst->plcl->pkd, in->achStarLogFile);
	}
    }

#endif

#ifdef SIMPLESF
void
pstSimpleStarForm(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	struct inSimpleStarForm *in = vin;
	struct outSimpleStarForm *out = vout;

	mdlassert(pst->mdl,nIn == sizeof(struct inSimpleStarForm));
	if (pst->nLeaves > 1) {
	    struct outSimpleStarForm fsStats;
	    
		mdlReqService(pst->mdl,pst->idUpper,PST_SIMPLESTARFORM,in,nIn);
		pstSimpleStarForm(pst->pstLower,in,nIn,vout,pnOut);
		mdlGetReply(pst->mdl,pst->idUpper,&fsStats,NULL);
		out->nFormed += fsStats.nFormed;
		out->nDeleted += fsStats.nDeleted;
		out->dMassFormed += fsStats.dMassFormed;
		}
	else {
		pkdSimpleStarForm(pst->plcl->pkd, in->dRateCoeff, in->dTMax, in->dDenMin, in->dDelta, in->dTime,
						  in->dInitStarMass, in->dESNPerStarMass, in->dtCoolingShutoff,in->bdivv,
						  &out->nFormed, &out->dMassFormed, &out->nDeleted);
		}

	if (pnOut) *pnOut = sizeof(struct outSimpleStarForm);
	}
#endif

void
pstDumpFrame(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	struct inDumpFrame *in = vin;

	mdlassert(pst->mdl,nIn == sizeof(struct inDumpFrame));
	mdlassert(pst->mdl,pnOut != NULL );
	if (pst->nLeaves > 1) {
		void * Image2;
		int nImage2;
		
		Image2 = malloc( DF_NBYTEDUMPFRAME );
		mdlassert(pst->mdl, Image2 != NULL );

		in->bNonLocal = 1;
	    mdlReqService(pst->mdl,pst->idUpper,PST_DUMPFRAME,in,nIn);
		in->bNonLocal = 0;
        pstDumpFrame(pst->pstLower,in,nIn, vout, pnOut);
	    mdlGetReply(pst->mdl,pst->idUpper, Image2, &nImage2);
		dfMergeImage( in, vout, pnOut, Image2, &nImage2 );

		free( Image2 );
		}
	else {
		PARTICLE *p = pst->plcl->pkd->pStore;
		dfClearImage( in, vout, pnOut );
		dfRenderParticlesInit( in, TYPE_GAS, TYPE_DARK, TYPE_STAR,
							   &p->r[0], &p->fMass, &p->fSoft, &p->fBall2, &p->iActive, 
#ifdef GASOLINE 
							   &p->fTimeForm,
#else
		/* N.B. This is just a place holder when we don't have stars */
							   &p->fMass,
#endif
							   p, sizeof(p[0]) );
		dfRenderParticles( in, vout, p, pst->plcl->pkd->nLocal );
		}
	}

void
pstCOM(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
  double *com = vout;
  int nOut = 0;

        if (pst->nLeaves > 1) {
		  double comtmp[12]; 
		  int i;
		  
		  mdlReqService(pst->mdl,pst->idUpper,PST_COM,vin,nIn);
		  
		  pstCOM(pst->pstLower, vin,nIn, vout, pnOut);
		  
		  mdlGetReply(pst->mdl,pst->idUpper, &comtmp[0], &nOut);
		  
		  for (i=0;i<12;i++) com[i] += comtmp[i];
		  }
        else {
                pkdCOM( pst->plcl->pkd, &com[0] );
                if (pnOut) *pnOut = sizeof(double)*12;
                }
        }


void
pstCOMByType(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
  double *com = vout;
  int *type = vin;
  int nOut = 0;

        if (pst->nLeaves > 1) {
		  double comtmp[4]; 
		  int i;
		  
		  mdlReqService(pst->mdl,pst->idUpper,PST_COMBYTYPE,vin,nIn);
		  
		  pstCOMByType(pst->pstLower, vin,nIn, vout, pnOut);
		  
		  mdlGetReply(pst->mdl,pst->idUpper, &comtmp[0], &nOut);
		  
		  for (i=0;i<4;i++) com[i] += comtmp[i];
		  }
        else {
                pkdCOMByType( pst->plcl->pkd, *type, &com[0] );
                if (pnOut) *pnOut = sizeof(double)*4;
                }
        }


void
pstOldestStar(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
  double *com = vout;
  int nOut = 0;

        if (pst->nLeaves > 1) {
		  double comtmp[12]; 
		  int i;
		  
		  mdlReqService(pst->mdl,pst->idUpper,PST_OLDESTSTAR,vin,nIn);
		  
		  pstOldestStar(pst->pstLower, vin,nIn, vout, pnOut);
		  
		  mdlGetReply(pst->mdl,pst->idUpper, &comtmp[0], &nOut);
		  
		  if (comtmp[3] < com[3]) for (i=0;i<4;i++) com[i] = comtmp[i];
		  }
        else {
                pkdOldestStar( pst->plcl->pkd, &com[0] );
                if (pnOut) *pnOut = sizeof(double)*4;
                }
        }


void
pstSetSink(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
  struct inSetSink *in = vin;
  struct outSetSink *out = vout, outtmp;

  assert(nIn == sizeof(struct inSetSink));
   
        if (pst->nLeaves > 1) {
		  mdlReqService(pst->mdl,pst->idUpper,PST_SETSINK,vin,nIn);
		  pstSetSink(pst->pstLower, vin,nIn, out, pnOut);
		  mdlGetReply(pst->mdl,pst->idUpper, &outtmp, pnOut);
		  out->nSink += outtmp.nSink;
		  }
        else {
                out->nSink = pkdSetSink( pst->plcl->pkd, in->dSinkMassMin );
                }
		if (pnOut) *pnOut = sizeof(struct outSetSink);
        }


#ifdef VOXEL
void
pstDumpVoxel(PST pst,void *vin,int nIn,void *vout,int *pnOut)
{
	struct inDumpVoxel *in = vin;

	mdlassert(pst->mdl,nIn == sizeof(struct inDumpVoxel));
	if (pst->nLeaves > 1) {
	    mdlReqService(pst->mdl,pst->idUpper,PST_DUMPFRAME,in,nIn);
        pstDumpFrame(pst->pstLower,in,nIn, NULL, NULL);
	    mdlGetReply(pst->mdl,pst->idUpper, NULL, NULL);
		}
	else {
		dfRenderVoxel(pst->plcl->pkd, in );
		}
	}
#endif
