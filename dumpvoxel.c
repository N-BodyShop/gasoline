#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>

#include "pkd.h"
#include "dumpvoxel.h"

void dfSetupVoxel( struct DumpFrameContext *df, double dTime, double dStep, struct inDumpVoxel *vin ) {

	vin->dTime = dTime;
	vin->dStep = dStep;
	
	vin->dMassStarMin = df->dMassStarMin;
	vin->dMassStarMax = df->dMassStarMax;
	vin->dMassGasMin  = df->dMassGasMin;
	vin->dMassGasMax  = df->dMassGasMax;
	vin->dMassDarkMin = df->dMassDarkMin;
	vin->dMassDarkMax = df->dMassDarkMax;

	vin->bVDetails = df->bVDetails;

	switch( df->iNumbering ) {
	case DF_NUMBERING_FRAME:
		sprintf(vin->fileout, df->FileName, df->nFrame );
		break;
	case DF_NUMBERING_STEP:
		sprintf(vin->fileout, df->FileName, dStep );
		break;
	case DF_NUMBERING_TIME:
		sprintf(vin->fileout, df->FileName, dTime );
		break;
		}

    }


void dfAddParticlesVoxel( PKD pkd, int type, double dMassMin, double dMassMax) {
	PARTICLE *p,*pend;

	pend = pkd->pStore + pkd->nLocal;

    if (dMassMax == DBL_MAX && dMassMin == 0) {
		for (p = pkd->pStore;p<pend;p++) {
			if (TYPETest( p, type)) {
				/* Add Particle */
				}
			}
		}
	else {
		for (p = pkd->pStore;p<pend;p++) {
			if (TYPETest( p, type )) {
				if (p->fMass <= dMassMax &&
					p->fMass >= dMassMin) {
					/* Add Particle */
					}
				}
			}
		}
	}


void dfRenderVoxel( PKD pkd, struct inDumpVoxel *in ) {

	if (in->dMassGasMin < DBL_MAX) {
		dfAddParticlesVoxel( pkd, TYPE_GAS, in->dMassGasMin, in->dMassGasMax);
		}
	if (in->dMassDarkMin < DBL_MAX) {
		dfAddParticlesVoxel( pkd, TYPE_DARK, in->dMassDarkMin, in->dMassDarkMax);
		}
	if (in->dMassStarMin < DBL_MAX) {
		dfAddParticlesVoxel( pkd, TYPE_STAR, in->dMassStarMin, in->dMassStarMax);
		}

	}


void dfFinishVoxel( struct DumpFrameContext *df, double dTime, double dStep, struct inDumpVoxel *in ) {

	df->nFrame++; /* NB: need to sort out something for restarts */
	
    if (df->dDumpFrameTime > 0 && dTime >= df->dTime)
        df->dTime = df->dTime + df->dDumpFrameTime;

	if (df->dDumpFrameStep > 0 && dStep >= df->dStep) 
        df->dStep = df->dStep + df->dDumpFrameStep;

	}
