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
#ifdef PLANETS
	double dStart; /* collision search time interval */
	double dEnd;
#endif /* PLANETS */
	} SMF;


typedef struct nNeighbor {
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

#define SMX_MEANVEL		2
void initMeanVel(void *);
void combMeanVel(void *,void *);
void MeanVel(PARTICLE *,int,NN *,SMF *);
void MeanVelSym(PARTICLE *,int,NN *,SMF *);

#ifdef GASOLINE

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

#ifdef PLANETS

#define SMX_TIMESTEP	7
void SetTimeStep(PARTICLE *,int,NN *,SMF *);

#define SMX_COLLISION	8
void CheckForCollision(PARTICLE *,int,NN *,SMF *);

#endif /* PLANETS */

#endif













