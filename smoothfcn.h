#ifndef SMOOTHFCN_INCLUDED
#define SMOOTHFCN_INCLUDED

#include "pkd.h"

typedef struct smfParameters {
	double H;
	double a;
	} SMF;


typedef struct nNeighbor {
	PARTICLE *pPart;
	float fDist2;
	float dx;
	float dy;
	float dz;
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

#define SMX_DENDIVV		3
void initDendivv(void *);
void combDendivv(void *,void *);
void postDendivv(PARTICLE *);
void Dendivv(PARTICLE *,int,NN *,SMF *);
void DendivvSym(PARTICLE *,int,NN *,SMF *);

#endif

#endif


