#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "pkd.h"
#include "dumpframe.h"


#define	DFGETCOEFF4_SCALAR( dt, a,b,c,d, t1,t2,tL,tR, v1,v2,vL,vR ) { \
double s1, s2; \
s1 = (v2-vL)/(t2-tL)*dt; \
s2 = (vR-v1)/(tR-t1)*dt; \
d = v1; \
c = s1; \
a = (s1+s2)-2*(v2-v1); \
b = v2 - a - c - d; \
}      									 
 																  
#define	DFGETCOEFF4_VEC( dt, a,b,c,d, t1,t2,tL,tR, v1,v2,vL,vR ) { \
DFGETCOEFF4_SCALAR( dt, a[0],b[0],c[0],d[0], t1,t2,tL,tR, v1[0],v2[0],vL[0],vR[0] ); \
DFGETCOEFF4_SCALAR( dt, a[1],b[1],c[1],d[1], t1,t2,tL,tR, v1[1],v2[1],vL[1],vR[1] ); \
DFGETCOEFF4_SCALAR( dt, a[2],b[2],c[2],d[2], t1,t2,tL,tR, v1[2],v2[2],vL[2],vR[2] ); \
}

void dfGetCoeff4( struct DumpFrameContext *df, int ifs ) {
	double dt = df->fs[ifs+1].dTime-df->fs[ifs].dTime; 

	df->rdt = 1./dt;

	DFGETCOEFF4_VEC( dt, df->a.eye, df->b.eye, df->c.eye, df->d.eye, 
			   df->fs[ifs].dTime, df->fs[ifs+1].dTime, df->fs[ifs-1].dTime, df->fs[ifs+2].dTime,
			   df->fs[ifs].eye, df->fs[ifs+1].eye, df->fs[ifs-1].eye, df->fs[ifs+2].eye );
	DFGETCOEFF4_VEC( dt, df->a.target, df->b.target, df->c.target, df->d.target, 
			   df->fs[ifs].dTime, df->fs[ifs+1].dTime, df->fs[ifs-1].dTime, df->fs[ifs+2].dTime,
			   df->fs[ifs].target, df->fs[ifs+1].target, df->fs[ifs-1].target, df->fs[ifs+2].target );
	DFGETCOEFF4_VEC( dt, df->a.up, df->b.up, df->c.up, df->d.up, 
			   df->fs[ifs].dTime, df->fs[ifs+1].dTime, df->fs[ifs-1].dTime, df->fs[ifs+2].dTime, 
			   df->fs[ifs].up, df->fs[ifs+1].up, df->fs[ifs-1].up, df->fs[ifs+2].up );
	DFGETCOEFF4_SCALAR( dt, df->a.FOV, df->b.FOV, df->c.FOV, df->d.FOV, 
			   df->fs[ifs].dTime, df->fs[ifs+1].dTime, df->fs[ifs-1].dTime, df->fs[ifs+2].dTime, 
			   df->fs[ifs].FOV, df->fs[ifs+1].FOV, df->fs[ifs-1].FOV, df->fs[ifs+2].FOV );
	DFGETCOEFF4_SCALAR( dt, df->a.zClipNear, df->b.zClipNear, df->c.zClipNear, df->d.zClipNear, 
			   df->fs[ifs].dTime, df->fs[ifs+1].dTime, df->fs[ifs-1].dTime, df->fs[ifs+2].dTime, 
			   df->fs[ifs].zClipNear, df->fs[ifs+1].zClipNear, df->fs[ifs-1].zClipNear, df->fs[ifs+2].zClipNear );
	DFGETCOEFF4_SCALAR( dt, df->a.zClipFar, df->b.zClipFar, df->c.zClipFar, df->d.zClipFar, 
			   df->fs[ifs].dTime, df->fs[ifs+1].dTime, df->fs[ifs-1].dTime, df->fs[ifs+2].dTime, 
			   df->fs[ifs].zClipFar, df->fs[ifs+1].zClipFar, df->fs[ifs-1].zClipFar, df->fs[ifs+2].zClipFar );
	}

#define	DFGETCOEFF3L_SCALAR( dt, a,b,c,d, t1,t2,tL, v1,v2,vL ) { \
double s1; \
s1 = (v2-vL)/(t2-tL)*dt; \
d = v1; \
c = s1; \
b = (v2-v1)-s1; \
a = 0; \
}      									 
 																  
#define	DFGETCOEFF3L_VEC( dt, a,b,c,d, t1,t2,tL, v1,v2,vL ) { \
DFGETCOEFF3L_SCALAR( dt, a[0],b[0],c[0],d[0], t1,t2,tL, v1[0],v2[0],vL[0] ); \
DFGETCOEFF3L_SCALAR( dt, a[1],b[1],c[1],d[1], t1,t2,tL, v1[1],v2[1],vL[1] ); \
DFGETCOEFF3L_SCALAR( dt, a[2],b[2],c[2],d[2], t1,t2,tL, v1[2],v2[2],vL[2] ); \
}

void dfGetCoeff3L( struct DumpFrameContext *df, int ifs ) {
	double dt = df->fs[ifs+1].dTime-df->fs[ifs].dTime; 

	df->rdt = 1./dt;

	DFGETCOEFF3L_VEC( dt, df->a.eye, df->b.eye, df->c.eye, df->d.eye, 
			   df->fs[ifs].dTime, df->fs[ifs+1].dTime, df->fs[ifs-1].dTime,
			   df->fs[ifs].eye, df->fs[ifs+1].eye, df->fs[ifs-1].eye );
	DFGETCOEFF3L_VEC( dt, df->a.target, df->b.target, df->c.target, df->d.target, 
			   df->fs[ifs].dTime, df->fs[ifs+1].dTime, df->fs[ifs-1].dTime,
			   df->fs[ifs].target, df->fs[ifs+1].target, df->fs[ifs-1].target );
	DFGETCOEFF3L_VEC( dt, df->a.up, df->b.up, df->c.up, df->d.up, 
			   df->fs[ifs].dTime, df->fs[ifs+1].dTime, df->fs[ifs-1].dTime, 
			   df->fs[ifs].up, df->fs[ifs+1].up, df->fs[ifs-1].up );
	DFGETCOEFF3L_SCALAR( dt, df->a.FOV, df->b.FOV, df->c.FOV, df->d.FOV, 
			   df->fs[ifs].dTime, df->fs[ifs+1].dTime, df->fs[ifs-1].dTime, 
			   df->fs[ifs].FOV, df->fs[ifs+1].FOV, df->fs[ifs-1].FOV );
	DFGETCOEFF3L_SCALAR( dt, df->a.zClipNear, df->b.zClipNear, df->c.zClipNear, df->d.zClipNear, 
			   df->fs[ifs].dTime, df->fs[ifs+1].dTime, df->fs[ifs-1].dTime,  
			   df->fs[ifs].zClipNear, df->fs[ifs+1].zClipNear, df->fs[ifs-1].zClipNear );
	DFGETCOEFF3L_SCALAR( dt, df->a.zClipFar, df->b.zClipFar, df->c.zClipFar, df->d.zClipFar, 
			   df->fs[ifs].dTime, df->fs[ifs+1].dTime, df->fs[ifs-1].dTime, 
			   df->fs[ifs].zClipFar, df->fs[ifs+1].zClipFar, df->fs[ifs-1].zClipFar );
	}

#define	DFGETCOEFF3R_SCALAR( dt, a,b,c,d, t1,t2,tR, v1,v2,vR ) { \
double s2; \
s2 = (vR-v1)/(tR-t1)*dt; \
d = v1; \
b = s2-(v2-v1); \
c = v2 - d - b; \
a = 0; \
}      									 
 																  
#define	DFGETCOEFF3R_VEC( dt, a,b,c,d, t1,t2,tR, v1,v2,vR ) { \
DFGETCOEFF3R_SCALAR( dt, a[0],b[0],c[0],d[0], t1,t2,tR, v1[0],v2[0],vR[0] ); \
DFGETCOEFF3R_SCALAR( dt, a[1],b[1],c[1],d[1], t1,t2,tR, v1[1],v2[1],vR[1] ); \
DFGETCOEFF3R_SCALAR( dt, a[2],b[2],c[2],d[2], t1,t2,tR, v1[2],v2[2],vR[2] ); \
}

void dfGetCoeff3R( struct DumpFrameContext *df, int ifs ) {
	double dt = df->fs[ifs+1].dTime-df->fs[ifs].dTime; 

	df->rdt = 1./dt;

	DFGETCOEFF3R_VEC( dt, df->a.eye, df->b.eye, df->c.eye, df->d.eye, 
			   df->fs[ifs].dTime, df->fs[ifs+1].dTime, df->fs[ifs+2].dTime,
			   df->fs[ifs].eye, df->fs[ifs+1].eye, df->fs[ifs+2].eye );
	DFGETCOEFF3R_VEC( dt, df->a.target, df->b.target, df->c.target, df->d.target, 
			   df->fs[ifs].dTime, df->fs[ifs+1].dTime, df->fs[ifs+2].dTime,
			   df->fs[ifs].target, df->fs[ifs+1].target, df->fs[ifs+2].target );
	DFGETCOEFF3R_VEC( dt, df->a.up, df->b.up, df->c.up, df->d.up, 
			   df->fs[ifs].dTime, df->fs[ifs+1].dTime, df->fs[ifs+2].dTime, 
			   df->fs[ifs].up, df->fs[ifs+1].up, df->fs[ifs+2].up );
	DFGETCOEFF3R_SCALAR( dt, df->a.FOV, df->b.FOV, df->c.FOV, df->d.FOV, 
			   df->fs[ifs].dTime, df->fs[ifs+1].dTime, df->fs[ifs+2].dTime, 
			   df->fs[ifs].FOV, df->fs[ifs+1].FOV, df->fs[ifs+2].FOV );
	DFGETCOEFF3R_SCALAR( dt, df->a.zClipNear, df->b.zClipNear, df->c.zClipNear, df->d.zClipNear, 
			   df->fs[ifs].dTime, df->fs[ifs+1].dTime, df->fs[ifs+2].dTime,  
			   df->fs[ifs].zClipNear, df->fs[ifs+1].zClipNear, df->fs[ifs+2].zClipNear );
	DFGETCOEFF3R_SCALAR( dt, df->a.zClipFar, df->b.zClipFar, df->c.zClipFar, df->d.zClipFar, 
			   df->fs[ifs].dTime, df->fs[ifs+1].dTime, df->fs[ifs+2].dTime, 
			   df->fs[ifs].zClipFar, df->fs[ifs+1].zClipFar, df->fs[ifs+2].zClipFar );
	}

#define	DFGETCOEFF2_SCALAR( dt, a,b,c,d, t1,t2, v1,v2 ) { \
d = v1; \
c = v2-v1; \
b = 0; \
a = 0; \
}      									 
 																  
#define	DFGETCOEFF2_VEC( dt, a,b,c,d, t1,t2, v1,v2 ) { \
DFGETCOEFF2_SCALAR( dt, a[0],b[0],c[0],d[0], t1,t2, v1[0],v2[0] ); \
DFGETCOEFF2_SCALAR( dt, a[1],b[1],c[1],d[1], t1,t2, v1[1],v2[1] ); \
DFGETCOEFF2_SCALAR( dt, a[2],b[2],c[2],d[2], t1,t2, v1[2],v2[2] ); \
}

void dfGetCoeff2( struct DumpFrameContext *df, int ifs ) {
	double dt = df->fs[ifs+1].dTime-df->fs[ifs].dTime; 

	df->rdt = 1./dt;

	DFGETCOEFF2_VEC( dt, df->a.eye, df->b.eye, df->c.eye, df->d.eye, 
			   df->fs[ifs].dTime, df->fs[ifs+1].dTime, 
			   df->fs[ifs].eye, df->fs[ifs+1].eye );
	DFGETCOEFF2_VEC( dt, df->a.target, df->b.target, df->c.target, df->d.target, 
			   df->fs[ifs].dTime, df->fs[ifs+1].dTime, 
			   df->fs[ifs].target, df->fs[ifs+1].target );
	DFGETCOEFF2_VEC( dt, df->a.up, df->b.up, df->c.up, df->d.up, 
			   df->fs[ifs].dTime, df->fs[ifs+1].dTime,
			   df->fs[ifs].up, df->fs[ifs+1].up );
	DFGETCOEFF2_SCALAR( dt, df->a.FOV, df->b.FOV, df->c.FOV, df->d.FOV, 
			   df->fs[ifs].dTime, df->fs[ifs+1].dTime,
			   df->fs[ifs].FOV, df->fs[ifs+1].FOV );
	DFGETCOEFF2_SCALAR( dt, df->a.zClipNear, df->b.zClipNear, df->c.zClipNear, df->d.zClipNear, 
			   df->fs[ifs].dTime, df->fs[ifs+1].dTime, 
			   df->fs[ifs].zClipNear, df->fs[ifs+1].zClipNear );
	DFGETCOEFF2_SCALAR( dt, df->a.zClipFar, df->b.zClipFar, df->c.zClipFar, df->d.zClipFar, 
			   df->fs[ifs].dTime, df->fs[ifs+1].dTime,  
			   df->fs[ifs].zClipFar, df->fs[ifs+1].zClipFar );
	}

void dfGetCoeff( struct DumpFrameContext *df, int ifs ) {
	if (ifs > 1 && ifs+2 < df->nFrameSetup) 
		dfGetCoeff4( df, ifs );
	else if (ifs > 1) 
		dfGetCoeff3L( df, ifs );
	else if (ifs+2 < df->nFrameSetup) 
		dfGetCoeff3R( df, ifs );
	else 
		dfGetCoeff2( df, ifs );
	}
