#ifndef WALK_HINCLUDED
#define WALK_HINCLUDED

#include "pkd.h"
#include "floattype.h"

#define CID_PARTICLE	0
#define CID_CELL		1


#define INTERSECTNP(pkdn,fBall2,x,y,z,bIntersect)\
{\
	FLOAT INTRSCT_dx,INTRSCT_dy,INTRSCT_dz;\
	FLOAT INTRSCT_dx1,INTRSCT_dy1,INTRSCT_dz1,INTRSCT_fDist2;\
	INTRSCT_dx = (pkdn)->bnd.fMin[0] - x;\
	INTRSCT_dx1 = x - (pkdn)->bnd.fMax[0];\
	INTRSCT_dy = (pkdn)->bnd.fMin[1] - y;\
	INTRSCT_dy1 = y - (pkdn)->bnd.fMax[1];\
	INTRSCT_dz = (pkdn)->bnd.fMin[2] - z;\
	INTRSCT_dz1 = z - (pkdn)->bnd.fMax[2];\
	if (INTRSCT_dx > 0.0) INTRSCT_fDist2 = INTRSCT_dx*INTRSCT_dx;\
	else if (INTRSCT_dx1 > 0.0) INTRSCT_fDist2 = INTRSCT_dx1*INTRSCT_dx1;\
	else INTRSCT_fDist2 = 0.0;\
	if (INTRSCT_dy > 0.0) INTRSCT_fDist2 += INTRSCT_dy*INTRSCT_dy;\
	else if (INTRSCT_dy1 > 0.0) INTRSCT_fDist2 += INTRSCT_dy1*INTRSCT_dy1;\
	if (INTRSCT_dz > 0.0) INTRSCT_fDist2 += INTRSCT_dz*INTRSCT_dz;\
	else if (INTRSCT_dz1 > 0.0) INTRSCT_fDist2 += INTRSCT_dz1*INTRSCT_dz1;\
	bIntersect = (INTRSCT_fDist2 <= fBall2);\
	}

void pkdBucketWalk(PKD,int,int,int);

#endif
