#ifndef COLLISION_HINCLUDED
#define COLLISION_HINCLUDED

#ifdef COLLISIONS

#include <limits.h> /* for DBL_MAX */

#include "pkd.h"

/*#define SMOOTH_STEP*/ /*DEBUG out of date*/

#define HILL_STEP /*DEBUG*/

#ifdef IN_A_BOX
#define BOX_HALF_SIZE 109.0 /*DEBUG hardwired*/
#endif /* IN_A_BOX */

#define NORMAL 0	/* drift types */
#define KEPLER 1

/*DEBUG HILL_SCALE should probably be a parameter*/
#define HILL_SCALE 1.26 /* multiplies reduced Hill radius */

#define REJECT(p) ((p)->dTEnc < 0)

#define COLLISION(t) ((t) < DBL_MAX)

#define BIT(n) (1 << (n))

#define NONE	0
#define MERGE	BIT(0)
#define BOUNCE	BIT(1)
#define FRAG	BIT(2)

#define MAX_NUM_FRAG 4

void partToCollider(PARTICLE *p, int iPid, int iIndex, COLLIDER *c);
void colliderToPart(COLLIDER *c, PARTICLE *p);
void pkdDoCollision(PKD pkd, int iOutcomes, double dEpsN, double dEpsT,
					COLLIDER *Collider1, COLLIDER *Collider2,
					double * pdImpactEnergy, int *piOutcome,
					COLLIDER *pOut, int *pnOut);
void pkdMerge(PKD pkd, COLLIDER **pOut, int *pnOut);
void pkdBounce(PKD pkd, double dEpsN, double dEpsT,
			   COLLIDER **pOut, int *pnOut);
void pkdFrag(PKD pkd, COLLIDER **pOut, int *pnOut);

#endif /* COLLISIONS */

#endif /* !COLLISION_HINCLUDED */
