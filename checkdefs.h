#ifndef CHECKDEFS_HINCLUDED
#define CHECKDEFS_HINCLUDED

#include "parameters.h"

#define MAX_REDSHIFTS	1000

struct msrCheckPointHeader {
	double dTime;
	int iStep;
	int N;
	struct parameters param;
	int iOpenType;
	double dCrit;
	/*
	 ** Save comoving coordinartes stuff.
	 */
	double dRedshift;
	double dHubble;
	double dCosmoFac;
	double dEcosmo;
	double dHubbleOld,dUOld,dTimeOld;
	/*
	 ** Save output redshifts.
	 */
	int nRed;
	double dRedOut[MAX_REDSHIFTS];
	};

#endif
