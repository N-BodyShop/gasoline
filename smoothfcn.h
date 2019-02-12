
#ifndef SMOOTHFCN_INCLUDED
#define SMOOTHFCN_INCLUDED

#include "pkd.h"
#include "floattype.h"

#ifdef SAND_PILE
#include "collision.h" /* for definition of WALLS */
#endif

#ifdef STARFORM
#include "supernova.h"
#endif

#define VISCOSITYLIMITER_NONE      0
#define VISCOSITYLIMITER_BALSARA   1
#define VISCOSITYLIMITER_JW        2

#define SMF_SMOOTHAGAIN    1

typedef struct smfParameters {
    double H;
    double a;
    double dDeltaAccelFac;
    double dSinkRadius;
    double dSinkBoundOrbitRadius;
    double dSinkMustAccreteRadius;
    double dBHSinkEddFactor;
    double dBHSinkFeedbackEff;
    double dBHSinkAlphaFactor;
    double dBHSinkFeedbackFactor;
    int bSmallBHSmooth;
    int bBHTurnOffCooling;
    int bBHMindv;
    int bBHAccreteAll;
    int bDoBHKick;
    int bStepZero;
    double dSinkCurrentDelta;
    double dDeltaStarForm;
#ifdef GASOLINE
    double dTauAlpha;
    double dAlphaMax;
    double dAlphaMin;
    double dNAlphaNoise;
    double dAFac;
    double alpha;
    double beta;
    double gamma;
    double algam;
    double Pext;
    double uMin;
    double dtMin; /* Read/Write */
    double dtFacCourant;
    double dtFacDiffusion;
    double dEtaCourantLong;
    double dDelta;
    int bGeometric;
    int bCannonical;
    int bGrowSmoothList;
    int iViscosityLimiter;
#endif
    int bSinkThermal;
    int bSinkFormDivV;
    int bSinkFormDivAcc;
    int bSinkFormDV;
    int bSinkFormPotMin;
    int nSinkFormMin;
    int iSinkCurrentRung;
    int iSmoothFlags; /* Read/Write locally.  Master sets initial value. */
#if defined(WENDLAND) || defined(WENDLANDC4)
    double Wzero;
#endif
    int nSmoothed;
    int nSmoothedInner;
    int nSmoothedFixh;
    double dSinkFormDivVCoeff;
    double dSinkFormDivAccCoeff;
    double dSinkTimeEligible;
    double dTime;
    double dEvapCoeffCode;
    double dEvapMinTemp;
    double dMetalDiffusionCoeff;
    double dThermalDiffusionCoeff;
    int bConstantDiffusion;
#ifdef TWOPHASE
    double dFBInitialMassLoad;
    double dMultiPhaseMinTemp;
    double dHotInitTemp;
    double dHotInitCodeDensity;
#endif
#ifdef STARFORM
    double dMinMassFrac;
    double dMaxGasMass;
    double dMinGasMass;
    int bShortCoolShutoff;
    int bSNTurnOffCooling;
    int bSmallSNSmooth;
    double dSecUnit;
    double dErgUnit;
    double dKmPerSecUnit;
    double dGmUnit;
    struct snContext sn;
    int bIonize;
    double dIonizeTime;
    double dIonizeMultiple;
    double dIonizeTMin;
    double dIonizeT;
    int bESF;
    double dESFTime;
    double dESFEnergy; 
    double dStarClusterRatio;
#endif    
#ifdef COLLISIONS
    double dDelta;
    double dCentMass;
    double dStart; /* collision search time interval */
    double dEnd;
    double dCollapseLimit; /* limit for inelastic collapse checks */
    int bFixCollapse;
    double dMaxBinaryEcc;
#endif
#ifdef SLIDING_PATCH
    PATCH_PARAMS PP;
#endif
#ifdef SAND_PILE
    WALLS walls;
#endif
    PKD pkd; /* useful for diagnostics, etc. */
    } SMF;


typedef struct nNeighbor {
	int iPid;
	int iIndex;
	PARTICLE *pPart;
	FLOAT fDist2;
	FLOAT dx;
	FLOAT dy;
	FLOAT dz;
	} NN;

typedef struct { double r2; NN *pNN; } ISORT;

int CompISORT(const void * a, const void * b);

enum smx_smoothtype {
  SMX_NULL,
  SMX_DENSITY,
  SMX_DENSITYTMP,
  SMX_MARKDENSITY,
  SMX_MARKIIDENSITY,
  SMX_MARK,
  SMX_DT,
  SMX_MEANVEL,
  SMX_DELTAACCEL,
  SMX_SINKACCRETETEST,
  SMX_SINKACCRETE,
  SMX_SINKINGAVERAGE,
  SMX_SINKINGFORCESHARE,
  SMX_BHDENSITY,
  SMX_BHSINKACCRETE,
  SMX_BHSINKIDENTIFY,
  SMX_BHSINKMERGE,
  SMX_SINKFORMTEST,
  SMX_SINKFORM,
  SMX_SINKMERGETEST,
  SMX_SINKMERGE,
#ifdef GASOLINE
  SMX_SPHPRESSURETERMS,
  SMX_DENDVDX,
  SMX_SURFACENORMAL,
  SMX_SURFACEAREA,
  SMX_SMOOTHBSW,
  SMX_DIVVORT,
  SMX_SHOCKTRACK,
  SMX_HKPRESSURETERMS,
  SMX_SPHPRESSURE,
  SMX_SPHVISCOSITY,
  SMX_HKVISCOSITY,
  SMX_DIST_DELETED_GAS,
#ifdef STARFORM
  SMX_STARCLUSTERFORM,
  SMX_PROMOTE_TO_HOT_GAS,
  SMX_EVAPORATE_TO_HOT_GAS,
  SMX_SHARE_WITH_HOT_GAS,
  SMX_DELETE_GAS,
  SMX_DIST_FB_ENERGY,
#endif
#endif
#ifdef COLLISIONS
  SMX_REJECTS,
  SMX_COLLISION,
  SMX_FINDBINARY,
#endif /* COLLISIONS */
#ifdef SLIDING_PATCH
  SMX_FIND_OVERLAPS,
#endif /* SLIDING_PATCH */
  };


/*  SMX_NULL */
void NullSmooth(PARTICLE *,int,NN *,SMF *);

/* SMX_DENSITY */
void initDensity(void *);
void combDensity(void *,void *);
void Density(PARTICLE *,int,NN *,SMF *);
void DensitySym(PARTICLE *,int,NN *,SMF *);

/* SMX_DENSITYTMP */
void initDensityTmp(void *);
void combDensityTmp(void *,void *);
void DensityTmp(PARTICLE *,int,NN *,SMF *);
void DensityTmpSym(PARTICLE *,int,NN *,SMF *);

/* SMX_MARKDENSITY */
void initParticleMarkDensity(void *);
void initMarkDensity(void *);
void combMarkDensity(void *,void *);
void MarkDensity(PARTICLE *,int,NN *,SMF *);
void MarkDensitySym(PARTICLE *,int,NN *,SMF *);

/* SMX_MARKIIDENSITY */
void initParticleMarkIIDensity(void *);
void initMarkIIDensity(void *);
void combMarkIIDensity(void *,void *);
void MarkIIDensity(PARTICLE *,int,NN *,SMF *);
void MarkIIDensitySym(PARTICLE *,int,NN *,SMF *);

/* SMX_MARK */
void initMark(void *);
void combMark(void *,void *);

/* SMX_MEANVEL */
void initMeanVel(void *);
void combMeanVel(void *,void *);
void MeanVel(PARTICLE *,int,NN *,SMF *);
void MeanVelSym(PARTICLE *,int,NN *,SMF *);

/* SMX_DELTAACCEL */
void DeltaAccel(PARTICLE *,int,NN *,SMF *);
void initDeltaAccel(void *);
void combDeltaAccel(void *,void *);

/* SMX_SINKACCRETETEST */
void SinkAccreteTest(PARTICLE *,int,NN *,SMF *);
void initSinkAccreteTest(void *);
void combSinkAccreteTest(void *,void *);

/* SMX_SINKACCRETE */
void SinkAccrete(PARTICLE *,int,NN *,SMF *);
void initSinkAccrete(void *);
void combSinkAccrete(void *,void *);

/* SMX_SINKINGAVERAGE */
void SinkingAverage(PARTICLE *,int,NN *,SMF *);

/* SMX_SINKINGFORCESHARE */
void SinkingForceShare(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf);
void initSinkingForceShare(void *);
void combSinkingForceShare(void *,void *);

/* SMX_BHDENSITY */
void initTreeParticleBHSinkDensity(void *);
void BHSinkDensity(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf);
void initBHSinkDensity(void *);
void combBHSinkDensity(void *,void *);

/* SMX_BHSINKACCRETE */
void initTreeParticleBHSinkAccrete(void *p1);
void BHSinkAccrete(PARTICLE *,int,NN *,SMF *);
void initBHSinkAccrete(void *);
void combBHSinkAccrete(void *,void *);
void postBHSinkAccrete(PARTICLE *p1, SMF *smf);

/* SMX_BHSINKIDENTIFY */
void BHSinkIdentify(PARTICLE *,int,NN *,SMF *);
void initBHSinkIdentify(void *);
void combBHSinkIdentify(void *,void *);

/* SMX_BHSINKMERGE */
void BHSinkMerge(PARTICLE *,int,NN *,SMF *);
void initBHSinkMerge(void *);
void combBHSinkMerge(void *,void *);

/* SMX_SINKMERGETEST */
void SinkMergeTest(PARTICLE *,int,NN *,SMF *);
void initSinkMergeTest(void *);
void combSinkMergeTest(void *,void *);

/* SMX_SINKMERGE */
void SinkMerge(PARTICLE *,int,NN *,SMF *);
void initSinkMerge(void *);
void combSinkMerge(void *,void *);

/* SMX_SINKFORMTEST */
void SinkFormTest(PARTICLE *,int,NN *,SMF *);
void initSinkFormTest(void *);
void combSinkFormTest(void *,void *);

/* SMX_SINKFORM */
void SinkForm(PARTICLE *,int,NN *,SMF *);
void initSinkForm(void *);
void combSinkForm(void *,void *);

#ifdef GASOLINE

/* SMX_SPHPRESSURETERMS */
void initSphPressureTermsParticle(void *);
void initSphPressureTerms(void *);
void combSphPressureTerms(void *,void *);
void SphPressureTerms(PARTICLE *,int,NN *,SMF *);
void SphPressureTermsSym(PARTICLE *,int,NN *,SMF *);

/* SMX_DENDVDX */
void initDenDVDX(void *);
void combDenDVDX(void *,void *);
void DenDVDX(PARTICLE *,int,NN *,SMF *);
void postDenDVDX(PARTICLE *, SMF *);

/* SMX_SURFACENORMAL */
void initSurfaceNormal(void *);
void combSurfaceNormal(void *,void *);
void SurfaceNormal(PARTICLE *,int,NN *,SMF *);

/* SMX_SURFACEAREA */
void initSurfaceArea(void *);
void combSurfaceArea(void *,void *);
void SurfaceArea(PARTICLE *,int,NN *,SMF *);

/* SMX_DENBSW */
void initSmoothBSw(void *);
void combSmoothBSw(void *,void *);
void SmoothBSw(PARTICLE *,int,NN *,SMF *);

/* SMX_DIVVORT */
void initDivVort(void *);
void combDivVort(void *,void *);
void DivVort(PARTICLE *,int,NN *,SMF *);
void DivVortSym(PARTICLE *,int,NN *,SMF *);

/* SMX_SHOCKTRACK */
void initShockTrack(void *);
void combShockTrack(void *,void *);
void ShockTrack(PARTICLE *,int,NN *,SMF *);
void ShockTrackSym(PARTICLE *,int,NN *,SMF *);

/* SMX_HKPRESSURETERMS */
void initHKPressureTermsParticle(void *);
void initHKPressureTerms(void *);
void combHKPressureTerms(void *,void *);
void HKPressureTerms(PARTICLE *,int,NN *,SMF *);
void HKPressureTermsSym(PARTICLE *,int,NN *,SMF *);

/* SMX_SPHPRESSURE */
void initSphPressureParticle(void *);
void initSphPressure(void *);
void combSphPressure(void *,void *);
void postSphPressure(PARTICLE *,SMF *);
void SphPressure(PARTICLE *,int,NN *,SMF *);
void SphPressureSym(PARTICLE *,int,NN *,SMF *);

/* SMX_SPHVISCOSITY */
void initSphViscosityParticle(void *);
void initSphViscosity(void *);
void combSphViscosity(void *,void *);
void SphViscosity(PARTICLE *,int,NN *,SMF *);
void SphViscositySym(PARTICLE *,int,NN *,SMF *);

/* SMX_HKVISCOSITY */
void initHKViscosityParticle(void *);
void initHKViscosity(void *);
void combHKViscosity(void *,void *);
void HKViscosity(PARTICLE *,int,NN *,SMF *);
void HKViscositySym(PARTICLE *,int,NN *,SMF *);

#ifdef STARFORM

/* SMX_STARCLUSTERFORM */
void initStarClusterForm(void *);
void combStarClusterForm(void *,void *);
void StarClusterForm(PARTICLE *,int,NN *,SMF *);

/* SMX_DIST_DELETED_GAS */
void initDistDeletedGas(void *p1);
void combDistDeletedGas(void *p1,void *p2);
void DistDeletedGas(PARTICLE *, int, NN *, SMF *);

/* SMX_EVAPORATE_TO_HOT_GAS */
void initEvaporateToHotGas(void *p1);
void combEvaporateToHotGas(void *p1,void *p2);
void EvaporateToHotGas(PARTICLE *, int, NN *, SMF *);

/* SMX_PROMOTE_TO_HOT_GAS */
void initPromoteToHotGas(void *p1);
void combPromoteToHotGas(void *p1,void *p2);
void PromoteToHotGas(PARTICLE *, int, NN *, SMF *);

/* SMX_SHARE_WITH_HOT_GAS */
void initShareWithHotGas(void *p1);
void combShareWithHotGas(void *p1,void *p2);
void ShareWithHotGas(PARTICLE *, int, NN *, SMF *);

/* SMX_DELETE_GAS */
void DeleteGas(PARTICLE *, int, NN *, SMF *);

/* SMX_DIST_FB_ENERGY */
void initTreeParticleDistFBEnergy(void *p1);
void initDistFBEnergy(void *p1);
void combDistFBEnergy(void *p1,void *p2);
void DistFBEnergy(PARTICLE *p, int, NN *, SMF *);
void postDistFBEnergy(PARTICLE *p1, SMF *smf);

#endif


#endif

#ifdef COLLISIONS

/* SMX_REJECTS */
void initFindRejects(void *p);
void combFindRejects(void *p1,void *p2);
void FindRejects(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf);


/* SMX_COLLISION */
void CheckForCollision(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf);

void FindBinary(PARTICLE *p, int nSmooth, NN *nnList, SMF *smf);

#endif /* COLLISIONS */

#ifdef SLIDING_PATCH

/* SMX_FIND_OVERLAPS */
void initFindOverlaps(void *p);
void combFindOverlaps(void *p1,void *p2);
void FindOverlaps(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf);

#endif /* SLIDING_PATCH */

#endif
