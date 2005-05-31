#ifndef SMOOTHFCN_INCLUDED
#define SMOOTHFCN_INCLUDED

#include "pkd.h"
#include "floattype.h"

#ifdef SAND_PILE
#include "collision.h" /* for definition of WALLS */
#endif


/* Benz Method (Default) */
#if !defined(PRES_MONAGHAN) && !defined(PRES_HK)
#define PRES_PDV(a,b) (a)
#define PRES_ACC(a,b) (a+b)
#endif

/* Monaghan Method */
#ifdef PRES_MONAGHAN
#define PRES_PDV(a,b) ((a+b)*0.5)
#define PRES_ACC(a,b) (a+b)
#endif

/* HK */
#ifdef PRES_HK
#define PRES_PDV(a,b) sqrt(a*b)
#define PRES_ACC(a,b) (sqrt(a*b)*2)
#endif

typedef struct smfParameters {
	double H;
	double a;
    double dDeltaAccelFac;
    double dSinkRadius;
    double dSinkBoundOrbitRadius;
#ifdef GASOLINE
	double alpha;
	double beta;
	double gamma;
	double algam;
	int bGeometric;
	int bCannonical;
	int bGrowSmoothList;
#endif
    int bSinkThermal;
#ifdef STARFORM
        double dMinMassFrac;
        double dRadPreFactor;
        double dTimePreFactor;
	double dTime;
	int bShortCoolShutoff;
	int bSmallSNSmooth;
#endif    
#ifdef COLLISIONS
	double dTime;
	double dDelta;
	double dCentMass;
	double dStart; /* collision search time interval */
	double dEnd;
	double dCollapseLimit; /* limit for inelastic collapse checks */
#endif
#ifdef SLIDING_PATCH
	double dOrbFreq;
	FLOAT fLx;
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
  SMX_SINKACCRETE,
#ifdef GASOLINE
  SMX_SPHPRESSURETERMS,
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
  SMX_COLLISION
#endif /* COLLISIONS */
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

/* SMX_SINKACCRETE */
void SinkAccrete(PARTICLE *,int,NN *,SMF *);
void initSinkAccrete(void *);
void combSinkAccrete(void *,void *);

#ifdef GASOLINE

/* SMX_SPHPRESSURETERMS */
void initSphPressureTermsParticle(void *);
void initSphPressureTerms(void *);
void combSphPressureTerms(void *,void *);
void SphPressureTerms(PARTICLE *,int,NN *,SMF *);
void SphPressureTermsSym(PARTICLE *,int,NN *,SMF *);

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
void combFindRejects(void *p1, void *p2);
void FindRejects(PARTICLE *p, int nSmooth, NN *nnList, SMF *smf);

#ifdef OLD_KEPLER
/* SMX_ENCOUNTER */
void CheckForEncounter(PARTICLE *p, int nSmooth, NN *nnList, SMF *smf);
#endif

/* SMX_COLLISION */
void CheckForCollision(PARTICLE *p, int nSmooth, NN *nnList, SMF *smf);

#endif /* COLLISIONS */

#endif
