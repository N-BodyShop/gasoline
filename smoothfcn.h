#ifndef SMOOTHFCN_INCLUDED
#define SMOOTHFCN_INCLUDED

#include "pkd.h"
#include "smooth.h"


#define SMX_DENSITY		1
void initDensity(void *);
void combDensity(void *,void *);
void Density(PARTICLE *,int,NN *);
void DensitySym(PARTICLE *,int,NN *);

#define SMX_MEANVEL		2
void initMeanVel(void *);
void combMeanVel(void *,void *);
void MeanVel(PARTICLE *,int,NN *);
void MeanVelSym(PARTICLE *,int,NN *);

#endif
