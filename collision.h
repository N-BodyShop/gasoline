#ifndef COLLISION_HINCLUDED
#define COLLISION_HINCLUDED

#ifdef PLANETS

#include <limits.h> /* for DBL_MAX */

#include "pkd.h"

/*#define OUTPUT_DT*/ /*DEBUG*/

/*#define SMOOTH_STEP*/ /*DEBUG*/

#define COLLISION(t) ((t) < DBL_MAX)

#define BIT(n) (1 << (n))

#define MERGE	BIT(0)
#define BOUNCE	BIT(1)
#define FRAG	BIT(2)

#define MAX_NUM_FRAG 4

void partToCollider(PARTICLE *,int,int,COLLIDER *);
void colliderToPart(COLLIDER *,PARTICLE *);
void pkdDoCollision(PKD,int,double,double,COLLIDER *,COLLIDER *,double *,
					int *,COLLIDER *,int *);
void pkdMerge(PKD,COLLIDER **,int *);
void pkdBounce(PKD,double,double,COLLIDER **,int *);
void pkdFrag(PKD,COLLIDER **,int *);

#endif /* PLANETS */

#endif /* !COLLISION_HINCLUDED */
