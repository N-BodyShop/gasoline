#ifndef WALK_HINCLUDED
#define WALK_HINCLUDED

#include "pkd.h"

#define CID_PARTICLE	0
#define CID_CELL		1


#define INTERSECTNP(pkdn,fBall2,x,y,z,bIntersect)\
{\
	float INTRSCT_dx,INTRSCT_dy,INTRSCT_dz;\
	float INTRSCT_dx1,INTRSCT_dy1,INTRSCT_dz1,INTRSCT_fDist2;\
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
	bIntersect = (INTRSCT_fDist2 < fBall2);\
	}


#ifdef COMPLETE_LOCAL
#define SETILIST(iOrder,ilcn,pkdn,x,y,z)\
{\
	switch (iOrder) {\
	case 4:\
		ilcn.xxxx = pkdn->mom.Hxxxx;\
		ilcn.xyyy = pkdn->mom.Hxyyy;\
		ilcn.xxxy = pkdn->mom.Hxxxy;\
		ilcn.yyyy = pkdn->mom.Hyyyy;\
		ilcn.xxxz = pkdn->mom.Hxxxz;\
		ilcn.yyyz = pkdn->mom.Hyyyz;\
		ilcn.xxyy = pkdn->mom.Hxxyy;\
		ilcn.xxyz = pkdn->mom.Hxxyz;\
		ilcn.xyyz = pkdn->mom.Hxyyz;\
		ilcn.xxzz = pkdn->mom.Hxxzz;\
		ilcn.xyzz = pkdn->mom.Hxyzz;\
		ilcn.xzzz = pkdn->mom.Hxzzz;\
		ilcn.yyzz = pkdn->mom.Hyyzz;\
		ilcn.yzzz = pkdn->mom.Hyzzz;\
		ilcn.zzzz = pkdn->mom.Hzzzz;\
	case 3:\
		ilcn.xxx = pkdn->mom.Oxxx;\
		ilcn.xyy = pkdn->mom.Oxyy;\
		ilcn.xxy = pkdn->mom.Oxxy;\
		ilcn.yyy = pkdn->mom.Oyyy;\
		ilcn.xxz = pkdn->mom.Oxxz;\
		ilcn.yyz = pkdn->mom.Oyyz;\
		ilcn.xyz = pkdn->mom.Oxyz;\
		ilcn.xzz = pkdn->mom.Oxzz;\
		ilcn.yzz = pkdn->mom.Oyzz;\
		ilcn.zzz = pkdn->mom.Ozzz;\
	case 2:\
		ilcn.xx = pkdn->mom.Qxx;\
		ilcn.yy = pkdn->mom.Qyy;\
                ilcn.zz = pkdn->mom.Qzz;\
		ilcn.xy = pkdn->mom.Qxy;\
		ilcn.xz = pkdn->mom.Qxz;\
		ilcn.yz = pkdn->mom.Qyz;\
	case 1:\
	default:\
		ilcn.m = pkdn->fMass;\
		ilcn.x = x;\
		ilcn.y = y;\
		ilcn.z = z;\
		}\
	}
#else
#define SETILIST(iOrder,ilcn,pkdn,x,y,z)\
{\
	switch (iOrder) {\
        double tr;\
	case 4:\
		ilcn.xxxx = pkdn->mom.Hxxxx;\
		ilcn.xyyy = pkdn->mom.Hxyyy;\
		ilcn.xxxy = pkdn->mom.Hxxxy;\
		ilcn.yyyy = pkdn->mom.Hyyyy;\
		ilcn.xxxz = pkdn->mom.Hxxxz;\
		ilcn.yyyz = pkdn->mom.Hyyyz;\
		ilcn.xxyy = pkdn->mom.Hxxyy;\
		ilcn.xxyz = pkdn->mom.Hxxyz;\
		ilcn.xyyz = pkdn->mom.Hxyyz;\
		ilcn.xxzz = pkdn->mom.Hxxzz;\
		ilcn.xyzz = pkdn->mom.Hxyzz;\
		ilcn.xzzz = pkdn->mom.Hxzzz;\
		ilcn.yyzz = pkdn->mom.Hyyzz;\
		ilcn.yzzz = pkdn->mom.Hyzzz;\
		ilcn.zzzz = pkdn->mom.Hzzzz;\
	case 3:\
		ilcn.xxx = pkdn->mom.Oxxx;\
		ilcn.xyy = pkdn->mom.Oxyy;\
		ilcn.xxy = pkdn->mom.Oxxy;\
		ilcn.yyy = pkdn->mom.Oyyy;\
		ilcn.xxz = pkdn->mom.Oxxz;\
		ilcn.yyz = pkdn->mom.Oyyz;\
		ilcn.xyz = pkdn->mom.Oxyz;\
		ilcn.xzz = pkdn->mom.Oxzz;\
		ilcn.yzz = pkdn->mom.Oyzz;\
		ilcn.zzz = pkdn->mom.Ozzz;\
	case 2:\
		tr = pkdn->mom.Qxx + pkdn->mom.Qyy + pkdn->mom.Qzz;\
		ilcn.xx = pkdn->mom.Qxx - tr/3.0;\
		ilcn.yy = pkdn->mom.Qyy - tr/3.0;\
                ilcn.zz = pkdn->mom.Qzz - tr/3.0;\
		ilcn.xy = pkdn->mom.Qxy;\
		ilcn.xz = pkdn->mom.Qxz;\
		ilcn.yz = pkdn->mom.Qyz;\
	case 1:\
	default:\
		ilcn.m = pkdn->fMass;\
		ilcn.x = x;\
		ilcn.y = y;\
		ilcn.z = z;\
		}\
	}
#endif


void pkdBucketWalk(PKD,int,int,int);

#endif






