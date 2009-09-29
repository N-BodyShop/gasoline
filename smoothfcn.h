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
    int bDoBHKick;
    double dSinkCurrentDelta;
  double dDeltaStarForm;
#ifdef GASOLINE
	double alpha;
	double beta;
	double gamma;
        double algam;
        double Pext;
        double uMin;
	int bGeometric;
	int bCannonical;
	int bGrowSmoothList;
        int iViscosityLimiter;
#endif
    int bSinkThermal;
    int nSinkFormMin;
    int iSinkCurrentRung;
    int iSmoothFlags; /* Read/Write locally.  Master sets initial value. */
    double dTime;
#ifdef DIFFUSION
    double dMetalDiffusionCoeff;
    double dThermalDiffusionCoeff;
    int bConstantDiffusion;
#endif
#ifdef STARFORM
        double dMinMassFrac;
	int bShortCoolShutoff;
        int bSNTurnOffCooling;
	int bSmallSNSmooth;
        double dSecUnit;
        double dErgUnit;
        double dKmPerSecUnit;
        double dGmUnit;
        struct snContext sn;
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


enum smx_smoothtype {
  SMX_NULL,
  SMX_DENSITY,
  SMX_MARKDENSITY,
  SMX_MARKIIDENSITY,
  SMX_MARK,
  SMX_MEANVEL,
  SMX_DELTAACCEL,
  SMX_SINKTEST,
  SMX_SINKACCRETE,
  SMX_SINKINGAVERAGE,
  SMX_SINKINGFORCESHARE,
  SMX_BHDENSITY,
  SMX_BHSINKACCRETE,
  SMX_BHSINKIDENTIFY,
  SMX_BHSINKMERGE,
  SMX_SINKFORMTEST,
  SMX_SINKFORM,
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
#ifdef STARFORM
  SMX_DIST_DELETED_GAS,
  SMX_DELETE_GAS,
  SMX_DIST_SN_ENERGY,
#endif
#ifdef SIMPLESF
  SMX_SIMPLESF_FEEDBACK,
#endif
#endif
#ifdef COLLISIONS
  SMX_REJECTS,
#ifdef OLD_KEPLER
  SMX_ENCOUNTER,
#endif
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

/* SMX_SINKTEST */
void SinkTest(PARTICLE *,int,NN *,SMF *);
void initSinkTest(void *);
void combSinkTest(void *,void *);

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
void initBHSinkMerge(void *,SMF *);

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

/* SMX_DIST_DELETED_GAS */
void initDistDeletedGas(void *p1);
void combDistDeletedGas(void *p1,void *p2);
void DistDeletedGas(PARTICLE *, int, NN *, SMF *);

/* SMX_DELETE_GAS */
void DeleteGas(PARTICLE *, int, NN *, SMF *);

/* SMX_DIST_SN_ENERGY */
void initTreeParticleDistSNEnergy(void *p1);
void initDistSNEnergy(void *p1);
void combDistSNEnergy(void *p1,void *p2);
void DistSNEnergy(PARTICLE *p, int, NN *, SMF *);
void postDistSNEnergy(PARTICLE *p1, SMF *smf);

#endif

#ifdef SIMPLESF
/* SMX_SIMPLESF_FEEDBACK */
void initSimpleSF_Feedback(void *p1);
void combSimpleSF_Feedback(void *p1,void *p2);
void SimpleSF_Feedback(PARTICLE *, int, NN *, SMF *);

#endif

#endif

#ifdef COLLISIONS

/* SMX_REJECTS */
void initFindRejects(void *p);
void combFindRejects(void *p1,void *p2);
void FindRejects(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf);

#ifdef OLD_KEPLER
/* SMX_ENCOUNTER */
void CheckForEncounter(PARTICLE *p,int nSmooth,NN *nnList,SMF *smf);
#endif

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
