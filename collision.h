#ifndef COLLISION_HINCLUDED
#define COLLISION_HINCLUDED

#include <limits.h> /* for DBL_MAX */

#include "pkd.h"

#define VERBOSE_COLLISION /*DEBUG*/

#define DT_TSC 0.0628319 /*DEBUG -- same as dEta*/

#define E_BOUNCE	DBL_MAX /*DEBUG -- should be parameter*/
#define E_FRAG		DBL_MAX /*DEBUG -- should be parameter*/

#define MAX_NUM_FRAG 4  /*DEBUG -- should be parameter?*/

#define EN 0.8 /*DEBUG -- should be parameter*/
#define ET 1.0 /*DEBUG -- should be parameter*/

#define COLLISION(t) ((t) < DBL_MAX)

void partToCollider(PARTICLE *,int,int,COLLIDER *);
void colliderToPart(COLLIDER *,PARTICLE *);
void pkdDoCollision(PKD,COLLIDER *,COLLIDER *,double *,COLLIDER *,int *);
void pkdMerge(PKD,COLLIDER **,int *);
void pkdBounce(PKD,COLLIDER **,int *);
void pkdFrag(PKD,COLLIDER **,int *);

#endif /* !COLLISION_HINCLUDED */
