
#ifndef COLLISION_HINCLUDED
#define COLLISION_HINCLUDED


#include "pkd.h" /* for PARTICLE struct */
#include "ssio.h" /* for SSDATA struct */
#ifdef COLLISIONS


#include "special.h"

#ifdef AGGS
#include "aggs.h" /* for Aggregate struct */
#endif

#ifdef RUBBLE_ZML
#include "rubble.h"
#endif

#define NORMAL 0	/* drift types */
#define KEPLER 1


FLOAT RADIUS(PARTICLE *p);
#define RADIUS(p) (2.0*(p)->fSoft) /* currently PARTICLE has no radius field;
									  to ensure all interactions are strictly
									  Newtonian, set the softening length to
									  be half the physical radius. */

int REJECT(PARTICLE *p);
#define REJECT(p) ((p)->dtCol < 0)

int COLLISION(double t);
#define COLLISION(t) ((t) < DBL_MAX)

unsigned int BIT(unsigned int n);
#define BIT(n) (1 << (n))

#define COLL_LOG_NONE		0
#define COLL_LOG_VERBOSE	1
#define COLL_LOG_TERSE		2

#define MISS	      0
#define MERGE	      BIT(0)
#define BOUNCE	      BIT(1)
#define FRAG	      BIT(2)
#define BINARY_MERGE  BIT(3)

#define MAX_NUM_FRAG 4

#ifdef SAND_PILE
#define MAX_NUM_WALLS 10
#ifdef TUMBLER
typedef struct {
	double n[3];		/*unit normal to flat wall or axis of cylinder */
	double ndotp;		/*dot product of unit normal n and point p on wall */
	double radius;		/* radius of cylinder */
	double omega;		/*rotation rate of cylinder */
	double dEpsN,dEpsT;
	double hotParam;	/*if nonzero, defines a permanent velocity of the wall
						  that sets a minimum energy for colliding particles */
	int type;			/*type = 0 for flat wall, = 1 for cylinder */
	} WALL;
#else
typedef struct {
	double x1,z1;
	double x2,z2;
	double dEpsN,dEpsT;
	} WALL;
#endif
typedef struct {
	int nWalls;
	WALL wall[MAX_NUM_WALLS];
	} WALLS;
#endif

enum {ConstEps,Frosty200,Frosty120,Compacted,Glancing}; /* bounce options */

enum {EscVel,MaxTrv}; /* slide options */

typedef struct {
	int iOutcomes;
	double dDensity;
	double dBounceLimit;
	int iBounceOption;
	double dEpsN;
	double dEpsT;
	int iSlideOption;
	double dSlideLimit;
	double dSlideLimit2;
	double dSlideEpsN;
	double dSlideEpsT;
	double dCollapseLimit;
	double dCollapseEpsN;
	double dCollapseEpsT;
	double dCrushLimit;
	double dCrushEpsN;
	double dCrushEpsT;
	int bFixCollapse;
#ifdef SAND_PILE
	WALLS walls;
#endif
#ifdef RUBBLE_ZML
	int iRubForcedOutcome;
	int iRubColor;
	double dRubbleMinFracMass;
	int bDoRubbleKDKRestart;
	DUST_BINS_PARAMS DB;
#endif
	} COLLISION_PARAMS;

typedef struct {
	int iPid;
	int iOrder;
	int iIndex;
	int iOrgIdx;
	} PARTICLE_ID;

typedef struct {
	PARTICLE_ID id;
	FLOAT fMass;
	FLOAT fRadius;
	FLOAT r[3];
	FLOAT v[3];
	FLOAT w[3];
        FLOAT a[3];
	FLOAT dt;
        FLOAT fDummy;
	int iRung;
	int iColor;
        int bTinyStep;
#ifdef AGGS
	Aggregate agg;
#endif
        int iPad;  /* This is a cheat to make the code work on intel
		      64 bit machines */
	} COLLIDER;

FLOAT SOFT(COLLIDER *c);
#define SOFT(c) (0.5*(c)->fRadius)

FLOAT SOFT_FROM_SSDATA(SSDATA *d);
#define SOFT_FROM_SSDATA(d) (0.5*(d)->radius) /* used in pkd.c */

#ifdef AGGS
int COLLIDER_IS_AGG(COLLIDER *c);
int COLLIDER_AGG_IDX(COLLIDER *c);

#define COLLIDER_IS_AGG(c) ((c)->id.iOrgIdx < 0)
#define COLLIDER_AGG_IDX(c) (-1 - (c)->id.iOrgIdx)

void pkdAggsDoCollision(PKD pkd,double dt,const COLLIDER *c1,const COLLIDER *c2,
						int bPeriodic,const COLLISION_PARAMS *CP,int iAggNewIdx,
						int *piOutcome,double *dT,COLLIDER *cOut,int *pnOut);
#endif /* AGGS */

void pkdNextCollision(PKD pkd, double *dtCol, int *iOrder1, int *iOrder2);
void pkdGetColliderInfo(PKD pkd, int iOrder, COLLIDER *c);
void PutColliderInfo(PKD pkd, const COLLIDER *c,int iOrder2,PARTICLE *p,
		     double dt);
void pkdDoCollision(PKD pkd, double dt, const COLLIDER *c1, const COLLIDER *c2,
 int bPeriodic, int iTime0, double dBaseStep, double dTimeNow,
const COLLISION_PARAMS *CP, int *piOutcome,double *dT,	COLLIDER *cOut, 
int *pnOut);
void pkdResetColliders(PKD pkd, int iOrder1, int iOrder2);
double LastKickTime(int iRung, double dBaseStep, double dTimeNow);
void pkdSetBall(PKD pkd,double dDelta,double fac);
void pkdFindTightestBinary(PKD pkd,double *dBindEn,int *iOrder1,int *iOrder2,
 int *n);
void SetMergerRung(const COLLIDER *c1,const COLLIDER *c2,COLLIDER *c,
				   double dBaseStep,double dTimeNow,int iTime0);
void MergerReorder(PKD pkd,const COLLIDER *pc1,const COLLIDER *pc2,const COLLIDER *c,
				   COLLIDER *cOut,double dt,const COLLISION_PARAMS *CP,
				   int bDiagInfo,const FLOAT fOffset[]
#ifdef RUBBLE_ZML
				   ,double dMassInDust
#endif
#ifdef SLIDING_PATCH
				   ,FLOAT fShear
#endif
				   ,int bPeriodic);
void pkdMergeBinary(PKD pkd,const COLLIDER *c1,const COLLIDER *c2,COLLIDER *c,
					int bPeriodic,double dBaseStep,double dTimeNow,int iTime0,
					double dDensity,int *bool);
#ifdef SLIDING_PATCH
#define MAXLARGEMASS 25	 /* Maximum number of particles to randomize */
#define MAXNEIGHBORS 250 /* Maximum neighbors surrounding a large
			    particle */
void pkdFindLargeMasses(PKD,double,double,double,double,PARTICLE *p,double *,int *);
void pkdGetNeighborParticles(PKD,double *,double,int,double,PARTICLE *p,double *,int *);
#endif

#endif /* COLLISIONS */

#endif
