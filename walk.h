#ifndef WALK_HINCLUDED
#define WALK_HINCLUDED

#include "pkd.h"

#define CID_PARTICLE	0
#define CID_CELL		1


#define INTERSECTNP(pkdn,fBall2,x,y,z)\
{\
	float dx,dy,dz,dx1,dy1,dz1,fDist2;\
	dx = pkdn->bnd.fMin[0] - x;\
	dx1 = x - pkdn->bnd.fMax[0];\
	dy = pkdn->bnd.fMin[1] - y;\
	dy1 = y - pkdn->bnd.fMax[1];\
	dz = pkdn->bnd.fMin[2] - z;\
	dz1 = z - pkdn->bnd.fMax[2];\
	if (dx > 0.0) fDist2 = dx*dx;\
	else if (dx1 > 0.0) fDist2 = dx1*dx1;\
	else fDist2 = 0.0;\
	if (dy > 0.0) fDist2 += dy*dy;\
	else if (dy1 > 0.0) fDist2 += dy1*dy1;\
	if (dz > 0.0) fDist2 += dz*dz;\
	else if (dz1 > 0.0) fDist2 += dz1*dz1;\
	if (fDist2 > fBall2) goto GetNextCell;\
	}


void pkdBucketWalk(PKD,int,int);

#endif
