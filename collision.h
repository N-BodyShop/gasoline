#ifndef COLLISION_HINCLUDED
#define COLLISION_HINCLUDED

#ifdef COLLISIONS

#include "pkd.h"

/*#define FIX_COLLAPSE*/

#define NORMAL 0	/* drift types */
#define KEPLER 1

#ifdef OLD_KEPLER
#define HILL_SCALE 1.26 /* multiplies reduced Hill radius */
#endif

#define REJECT(p) ((p)->dtCol < 0)

#define COLLISION(t) ((t) < DBL_MAX)

#define BIT(n) (1 << (n))

#define MISS	0
#define MERGE	BIT(0)
#define BOUNCE	BIT(1)
#define FRAG	BIT(2)

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
#ifdef SAND_PILE
	WALLS walls;
#endif
	} COLLISION_PARAMS;

typedef struct {
	int iPid;
	int iIndex;
	int iOrder;
	} PARTICLE_ID;

typedef struct {
	PARTICLE_ID id;
	FLOAT fMass;
	FLOAT fRadius; /*DEBUG = 2*fSoft so spline not used in force calcs*/
	FLOAT r[3];
	FLOAT v[3];
	FLOAT w[3];
	FLOAT dt;
	int iRung;
	int iColor;
	int bTinyStep;
	} COLLIDER;

void pkdNextCollision(PKD pkd, double *dtCol, int *iOrder1, int *iOrder2);
void pkdGetColliderInfo(PKD pkd, int iOrder, COLLIDER *c);
void pkdDoCollision(PKD pkd, double dt, const COLLIDER *c1, const COLLIDER *c2,
					int bPeriodic, const COLLISION_PARAMS *CP, int *piOutcome,
					double *dT,	COLLIDER *cOut, int *pnOut);
void pkdResetColliders(PKD pkd, int iOrder1, int iOrder2);

#endif /* COLLISIONS */

#endif
