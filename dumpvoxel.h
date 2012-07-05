#ifndef DUMPVOXEL_HINCLUDED
#define DUMPVOXEL_HINCLUDED

/* Rewrite to use particles Interface in dumpframe! so pkd can be avoided */
#include "pkd.h"


#include "dumpframe.h"
/* PST */

struct inDumpVoxel {
	double dTime;
	double dStep;

	/* Other info */
	int iRender;       /* Rendering */
	int iEncode;       /* Encoding */
	
	/* Particle Filters */
	double dMassStarMin,dMassStarMax;
	double dMassGasMin, dMassGasMax;
	double dMassDarkMin,dMassDarkMax;

	int bNonLocal; /* Is this data going non-local? */
	int bVDetails;
	char fileout[160];
    };


void dfSetupVoxel( struct DumpFrameContext *df, double dTime, double dStep, struct inDumpVoxel *vin );

void dfAddParticlesVoxel( PKD pkd, int type, double dMassMin, double dMassMax);

void dfRenderVoxel( PKD pkd, struct inDumpVoxel *in );

void dfFinishVoxel( struct DumpFrameContext *df, double dTime, double dStep, struct inDumpVoxel *in );

#endif

