#ifndef SMOOTHFCN_INCLUDED
#define SMOOTHFCN_INCLUDED

#include "pkd.h"
#include "floattype.h"

typedef struct smfParameters {
	double H;
	double a;
	double alpha;
	double beta;
	double gamma;
	double algam;
	int bGeometric;
        int bCannonical;
        int bGrowSmoothList;
#ifdef COLLISIONS
	double dTime;
	double dDelta;
	double dCentMass;
	double dStart;	/* collision search time interval */
	double dEnd;
	PKD pkd;		/* pointer to processor's PKD structure */
#endif /* COLLISIONS */
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

#define SMX_DENSITY		1
void initDensity(void *);
void combDensity(void *,void *);
void Density(PARTICLE *,int,NN *,SMF *);
void DensitySym(PARTICLE *,int,NN *,SMF *);

#define SMX_MARKDENSITY		6
void initMarkDensity(void *);
void combMarkDensity(void *,void *);
void MarkDensity(PARTICLE *,int,NN *,SMF *);
void MarkDensitySym(PARTICLE *,int,NN *,SMF *);

#define SMX_MEANVEL		2
void initMeanVel(void *);
void combMeanVel(void *,void *);
void MeanVel(PARTICLE *,int,NN *,SMF *);
void MeanVelSym(PARTICLE *,int,NN *,SMF *);

#ifdef GASOLINE

#define SMX_SPHPRESSURETERMS    3
void initSphPressureTermsParticle(void *);
void initSphPressureTerms(void *);
void combSphPressureTerms(void *,void *);
void SphPressureTerms(PARTICLE *,int,NN *,SMF *);
void SphPressureTermsSym(PARTICLE *,int,NN *,SMF *);

#define SMX_DIVVORT             4
void initDivVort(void *);
void combDivVort(void *,void *);
void DivVort(PARTICLE *,int,NN *,SMF *);
void DivVortSym(PARTICLE *,int,NN *,SMF *);

#define SMX_HKPRESSURETERMS    5
void initHKPressureTermsParticle(void *);
void initHKPressureTerms(void *);
void combHKPressureTerms(void *,void *);
void HKPressureTerms(PARTICLE *,int,NN *,SMF *);
void HKPressureTermsSym(PARTICLE *,int,NN *,SMF *);

#endif

#ifdef OLDGASOLINE
#define SMX_HSMDIVV		3
void initHsmDivv(void *);
void combHsmDivv(void *,void *);
void postHsmDivv(PARTICLE *,SMF *);
void HsmDivv(PARTICLE *,int,NN *,SMF *);
void HsmDivvSym(PARTICLE *,int,NN *,SMF *);

#define SMX_GEOMBV		4
void initGeomBV(void *);
void combGeomBV(void *,void *);
void postGeomBV(PARTICLE *,SMF *);
void GeomBVSym(PARTICLE *,int,NN *,SMF *);

#define SMX_ETHDOTBV	5
void initEthdotBV(void *);
void combEthdotBV(void *,void *);
void EthdotBVSym(PARTICLE *,int,NN *,SMF *);

void initAccsph(void *);
void combAccsph(void *,void *);

#define SMX_ACCSPHBV	6
void AccsphBVSym(PARTICLE *,int,NN *,SMF *);

#endif

#ifdef COLLISIONS

#define SMX_REJECTS		7
void FindRejects(PARTICLE *p, int nSmooth, NN *nnList, SMF *smf);

#define SMX_TIMESTEP	8
void SetTimeStep(PARTICLE *p, int nSmooth, NN *nnList, SMF *smf);

#define SMX_ENCOUNTER	9
void CheckForEncounter(PARTICLE *p, int nSmooth, NN *nnList, SMF *smf);

#define SMX_COLLISION	10
void CheckForCollision(PARTICLE *p, int nSmooth, NN *nnList, SMF *smf);

#endif /* COLLISIONS */

#endif
